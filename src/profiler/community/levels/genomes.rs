//! Genomes within each kept species.
//!
//! The fragments that the species-level EM gives to a species are shared
//! among the species' genomes by the same match-quality EM. The genome with the
//! largest share is the species' best-fit genome, and the mean identity of the
//! species' fragments to it estimates how close the sample's strain is to it.
//! At the genome level (`--strain`) a genome is also reported on its own when
//! the fragments that clearly belong to it and match it at strain-level
//! identity cover it evenly. A genome represented by marker genes never is:
//! its genes are conserved and too short to tell strains apart. A genome
//! reported on its own counts each of its fragments by the species' share of
//! the fragment in the species-level EM times the genome EM's posterior, so
//! that a fragment counts once within its species' rows: the species row
//! keeps the rest of the species' reads.

use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use super::{
    add_placement, candidate_errors, estimated_ani, fold_candidate_of, fold_candidates,
    longest_alignment, placement_identity, resolved_unit, tie_margin, tied_errors, with_errors,
    CoverageTests, EmSettings, Fragment, GenomeEvidence, LONG_TIE_SHARE, MIN_JUDGED_SPECIES,
};
use crate::profiler::community::MappingData;
use crate::profiler::evidence::{
    EvidenceCoverage, EvidenceMode, EvidencePolicy, TerminalCandidate,
};
use crate::profiler::likelihood_em::{self, ReadGroups};
use crate::profiler::ProfileRank;

/// Genome candidates of a long fragment whose errors exceed the best by less
/// than this share of the fragment's longest placement count as tied; short
/// fragments tie as species candidates do ([`tie_margin`]).
const GENOME_LONG_TIE_SHARE: f64 = 0.001;
/// A fragment matches a genome at strain-level identity when it aligns to the
/// genome at an identity of at least the sample's identity to its own genomes
/// less this margin: a strain within 99.5 % ANI of the genome, for every read
/// technology.
const STRAIN_MARGIN: f64 = 0.005;
/// A genome is reported on its own only when the fragments that the genome EM
/// gives to it match it within sequencing error: an estimated ANI of at least
/// this. Those fragments are the ones it fits best among the species' genomes,
/// so a relative's estimate reads high; for a genome in the sample it stays
/// within this of 1.
const GENOME_ROW_ANI: f64 = 0.999;

/// One fragment of a species: its read, mapped bases for depth, the
/// species' share of it in the species-level EM, and per genome of the
/// species `(assembly, errors, identity)`.
struct GenomeFragment {
    read: usize,
    bases: f64,
    share: f64,
    genomes: Vec<(u32, f64, f64)>,
}

/// A genome reported on its own, with its reads and their bases: each
/// fragment counted by the species' share of it times the genome EM's
/// posterior.
pub(super) struct StrainGenome {
    pub(super) assembly_id: u32,
    pub(super) reads: f64,
    pub(super) bases: f64,
    pub(super) mean_identity: f64,
    pub(super) coverage: EvidenceCoverage,
    pub(super) asm_len: f64,
}

/// What the genome level finds for one species.
pub(super) struct SpeciesGenomes {
    pub(super) best_fit: u32,
    /// Mean identity to the best-fit genome of the fragments whose best species
    /// match is this species, at any identity the mapper reports.
    pub(super) mean_identity: f64,
    pub(super) strains: Vec<StrainGenome>,
    /// With [`GenomeLevel::keep_shares`], each fragment's shares of the
    /// genomes reported on their own.
    pub(super) shares: StrainShares,
}

/// The fragments that the genome EM over the genomes reported on their own
/// shares among them: each fragment's read and read group, and each group's
/// posteriors over those genomes.
#[derive(Default)]
pub(super) struct StrainShares {
    pub(super) members: Vec<(usize, usize)>,
    pub(super) posteriors: Vec<Vec<(u32, f64)>>,
}

/// One genome's genome-level tests, for the decisions table.
pub(super) struct StrainDiagnostic {
    pub(super) taxid: u32,
    pub(super) assembly_id: u32,
    pub(super) best_fit: bool,
    pub(super) em_share: f64,
    pub(super) resolved: u64,
    pub(super) gated: u64,
    pub(super) coverage: Option<EvidenceCoverage>,
    pub(super) kept: bool,
}

/// The genome level of a sample.
pub(super) struct GenomeResult {
    pub(super) species: HashMap<u32, SpeciesGenomes>,
    /// Identity of reads to the genome they come from, measured on the
    /// species: the fragment-weighted median of their mean identity to their
    /// best-fit genome (`SpeciesGenomes::mean_identity`).
    pub(super) baseline_identity: f64,
    /// The same median over the fragments that match the best-fit genome at
    /// strain-level identity: the baseline of genome rows, whose fragments are
    /// all such matches.
    pub(super) strain_baseline_identity: f64,
    pub(super) diagnostics: Vec<StrainDiagnostic>,
}

/// Strain-level identity and the baseline that genome rows are measured
/// against.
#[derive(Clone, Copy)]
struct StrainIdentity {
    identity: f64,
    baseline: f64,
}

/// One species' genome EM.
struct SpeciesShare {
    taxid: u32,
    members: Vec<GenomeFragment>,
    group_of: Vec<Option<usize>>,
    em: likelihood_em::EmResult,
    best_fit: u32,
}

/// Settings shared by every species.
pub(super) struct GenomeLevel {
    pub(super) policy: EvidencePolicy,
    /// The sample's error rate, which sets the EM's error penalty.
    pub(super) error_rate: f64,
    pub(super) em: EmSettings,
    /// Report genomes on their own (`--strain`).
    pub(super) report_strains: bool,
    /// Keep each fragment's shares of the genomes reported on their own, for
    /// the per-read classifications.
    pub(super) keep_shares: bool,
    pub(super) chunk_cutoff: f64,
    pub(super) ratio_bounds: (f64, f64),
}

impl GenomeLevel {
    /// The coverage tests of a genome reported on its own, corrected for
    /// `multiple_tests` genomes.
    pub(super) fn tests(&self, multiple_tests: usize) -> CoverageTests {
        CoverageTests {
            multiple_tests,
            chunk_cutoff: self.chunk_cutoff,
            ratio_bounds: self.ratio_bounds,
            // Marker genomes are never reported on their own.
            markers: None,
            min_marker_genes: 0,
        }
    }
}

/// For each species of `species_of` (each fragment's species and its share of
/// the fragment), share its fragments among its genomes; then measure the
/// sample's baseline identity on the species (not on the units labelled above
/// species, `labels`) and, with `report_strains`, find the genomes reported
/// on their own.
pub(super) fn resolve(
    data: &MappingData,
    fragments: &[Fragment],
    species_of: &[Option<(u32, f64)>],
    labels: &HashMap<u32, ProfileRank>,
    level: &GenomeLevel,
) -> GenomeResult {
    let shares = genome_fragments(data, fragments, species_of, level.policy)
        .into_iter()
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|(taxid, members)| {
            let mut groups = ReadGroups::default();
            let mut candidates = Vec::new();
            let group_of = members
                .iter()
                .map(|m| {
                    candidates.clear();
                    candidates.extend(m.genomes.iter().map(|g| (g.0, g.1)));
                    groups.add(&mut candidates)
                })
                .collect::<Vec<_>>();
            let em = level.em.run(&groups);
            let best_fit = em
                .weight
                .iter()
                .max_by(|a, b| a.1.total_cmp(b.1).then(b.0.cmp(a.0)))
                .map(|(&assembly, _)| assembly)
                .expect("a species has at least one genome");
            SpeciesShare {
                taxid,
                members,
                group_of,
                em,
                best_fit,
            }
        })
        .collect::<Vec<_>>();
    let best_fit = shares.iter().map(|s| (s.taxid, s.best_fit)).collect();
    let identity = identity_to_best_fit(data, &best_fit, level.policy);
    let baseline_identity = weighted_median(
        identity
            .iter()
            .filter(|(taxid, _)| !labels.contains_key(taxid))
            .map(|(_, &value)| value),
    )
    .unwrap_or(1.0 - level.error_rate);
    let strain_identity = baseline_identity - STRAIN_MARGIN;
    let strain = StrainIdentity {
        identity: strain_identity,
        baseline: weighted_median(
            shares
                .iter()
                .filter(|s| !labels.contains_key(&s.taxid))
                .map(|s| strain_level_identity(&s.members, s.best_fit, strain_identity)),
        )
        .unwrap_or(baseline_identity),
    };
    // Species are independent; their results are taken in their order.
    let found = shares
        .into_par_iter()
        .map(|share| {
            let mut diagnostics = Vec::new();
            let (strains, shares) = if level.report_strains {
                strain_genomes(data, &share, strain, level, &mut diagnostics)
            } else {
                (Vec::new(), StrainShares::default())
            };
            let genomes = SpeciesGenomes {
                best_fit: share.best_fit,
                mean_identity: identity.get(&share.taxid).map_or(f64::NAN, |m| m.0),
                strains,
                shares,
            };
            (share.taxid, genomes, diagnostics)
        })
        .collect::<Vec<_>>();
    let mut species = HashMap::new();
    let mut diagnostics = Vec::new();
    for (taxid, genomes, rows) in found {
        species.insert(taxid, genomes);
        diagnostics.extend(rows);
    }
    GenomeResult {
        species,
        baseline_identity,
        strain_baseline_identity: strain.baseline,
        diagnostics,
    }
}

/// Mean identity to `assembly` of the fragments that match it at
/// `strain_identity` or above, and their number.
fn strain_level_identity(
    members: &[GenomeFragment],
    assembly: u32,
    strain_identity: f64,
) -> (f64, f64) {
    let (sum, count) = members
        .iter()
        .filter_map(|m| m.genomes.iter().find(|g| g.0 == assembly).map(|g| g.2))
        .filter(|&identity| identity >= strain_identity)
        .fold((0.0, 0.0), |(sum, count), identity| {
            (sum + identity, count + 1.0)
        });
    (if count > 0.0 { sum / count } else { f64::NAN }, count)
}

/// Median of `(value, weight)` pairs (one per species) with finite values
/// and positive weight; `None` for fewer than [`MIN_JUDGED_SPECIES`].
fn weighted_median(values: impl Iterator<Item = (f64, f64)>) -> Option<f64> {
    let mut values = values
        .filter(|v| v.0.is_finite() && v.1 > 0.0)
        .collect::<Vec<_>>();
    if values.len() < MIN_JUDGED_SPECIES {
        return None;
    }
    values.sort_by(|a, b| a.0.total_cmp(&b.0));
    let half = values.iter().map(|v| v.1).sum::<f64>() / 2.0;
    let mut seen = 0.0;
    values.iter().find_map(|&(value, weight)| {
        seen += weight;
        (seen >= half).then_some(value)
    })
}

/// Per species, its fragments with their errors and identity on each of its
/// genomes; genome candidates within the tie margin get the best errors.
fn genome_fragments(
    data: &MappingData,
    fragments: &[Fragment],
    species_of: &[Option<(u32, f64)>],
    policy: EvidencePolicy,
) -> HashMap<u32, Vec<GenomeFragment>> {
    let resolved = fragments
        .iter()
        .zip(species_of)
        .filter_map(|(fragment, species)| species.map(|(taxid, share)| (fragment, taxid, share)))
        .collect::<Vec<_>>();
    let mut out: HashMap<u32, Vec<GenomeFragment>> = HashMap::new();
    fold_candidate_of(
        data,
        policy,
        &resolved,
        (|r| r.0.read, |r| r.1),
        |_, candidate| Some(candidate_genomes(candidate?, policy)),
        |&(fragment, taxid, share), genomes| {
            if let Some(genomes) = genomes {
                out.entry(taxid).or_default().push(GenomeFragment {
                    read: fragment.read,
                    bases: fragment.bases,
                    share,
                    genomes,
                });
            }
        },
    );
    out
}

/// A species candidate's genomes with errors and identity; genomes within
/// the tie margin of the best get the best errors.
fn candidate_genomes(
    candidate: &TerminalCandidate,
    policy: EvidencePolicy,
) -> Vec<(u32, f64, f64)> {
    let longest = longest_alignment(candidate.placements.iter());
    let mut genomes = candidate
        .placements
        .iter()
        .map(|p| {
            (
                p.assembly_id,
                candidate_errors(p, longest),
                placement_identity(p),
            )
        })
        .collect::<Vec<_>>();
    let fewest = genomes.iter().map(|g| g.1).fold(f64::INFINITY, f64::min);
    let margin = tie_margin(policy.mode, longest, GENOME_LONG_TIE_SHARE);
    for genome in &mut genomes {
        genome.1 = tied_errors(genome.1, fewest, margin);
    }
    genomes
}

/// Per species with a best-fit genome in `best_fit`, the mean identity to
/// that genome and the number of fragments behind it. A fragment counts when
/// the species is its best species match at any identity the mapper reports
/// (fewest errors, within the species-level tie margin) and it aligns to the
/// best-fit genome. With the species identity cutoff off, the mean of a
/// distant strain is not cut off from below; fragments that match another
/// species better stay out.
fn identity_to_best_fit(
    data: &MappingData,
    best_fit: &HashMap<u32, u32>,
    policy: EvidencePolicy,
) -> HashMap<u32, (f64, f64)> {
    let ungated = EvidencePolicy {
        identity: 0.0,
        ..policy
    };
    // Only reads that align to a best-fit genome can add an identity.
    let genomes: HashSet<u32> = best_fit.values().copied().collect();
    let reads = (0..data.read_alignments.len())
        .into_par_iter()
        .filter(|&read| {
            data.read_alignments[read]
                .alignments
                .iter()
                .any(|a| genomes.contains(&a.assembly_id))
        })
        .collect::<Vec<_>>();
    let mut sums: HashMap<u32, (f64, f64)> = HashMap::new();
    fold_candidates(
        data,
        ungated,
        &reads,
        |&read| read,
        |_, candidates| best_fit_identities(candidates, best_fit, policy.mode),
        |_, identities| {
            for (taxid, identity) in identities {
                let entry = sums.entry(taxid).or_insert((0.0, 0.0));
                entry.0 += identity;
                entry.1 += 1.0;
            }
        },
    );
    sums.into_iter()
        .map(|(taxid, (sum, count))| (taxid, (sum / count, count)))
        .collect()
}

/// A fragment's identity to the best-fit genome of each unit in `best_fit`
/// that is its best match among its candidate units (fewest errors, within
/// the species-level tie margin), in candidate order.
fn best_fit_identities(
    candidates: &[TerminalCandidate],
    best_fit: &HashMap<u32, u32>,
    mode: EvidenceMode,
) -> Vec<(u32, f64)> {
    let species = with_errors(candidates, |_| true);
    // The best is the first candidate with the fewest errors, in candidate
    // order; its aligned bases set the tie margin of long reads.
    let Some(&(fewest, best)) = species.iter().min_by(|a, b| a.0.total_cmp(&b.0)) else {
        return Vec::new();
    };
    let margin = tie_margin(mode, f64::from(best.best.both_mapped), LONG_TIE_SHARE);
    species
        .iter()
        .filter(|&&(errors, _)| errors - fewest < margin)
        .filter_map(|&(_, candidate)| {
            let placement = candidate.placement_on(*best_fit.get(&candidate.taxid)?)?;
            Some((candidate.taxid, placement_identity(placement)))
        })
        .collect()
}

/// Genomes of one species reported on their own: strain-level coverage from
/// fragments that resolve to a genome (posterior >= RESOLVE_POSTERIOR) at
/// strain-level identity, the chunk-breadth cutoff and breadth ratios, and a
/// second genome EM over the kept genomes for their reads and depth, with each
/// fragment's shares of them when [`GenomeLevel::keep_shares`].
fn strain_genomes(
    data: &MappingData,
    share: &SpeciesShare,
    strain: StrainIdentity,
    level: &GenomeLevel,
    diagnostics: &mut Vec<StrainDiagnostic>,
) -> (Vec<StrainGenome>, StrainShares) {
    let (taxid, members, group_of, posteriors) = (
        share.taxid,
        &share.members,
        &share.group_of,
        &share.em.posteriors,
    );
    let strain_identity = strain.identity;
    let mut resolved: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut strain_level = Vec::new();
    for (member, group) in members.iter().zip(group_of) {
        let Some(assembly) = group.and_then(|g| resolved_unit(&posteriors[g])) else {
            continue;
        };
        if data
            .marker_model
            .as_deref()
            .is_some_and(|model| model.is_marker_genome(assembly))
        {
            continue;
        }
        let identity = member
            .genomes
            .iter()
            .find(|g| g.0 == assembly)
            .map_or(0.0, |g| g.2);
        let counts = resolved.entry(assembly).or_insert((0, 0));
        counts.0 += 1;
        if identity < strain_identity {
            continue;
        }
        counts.1 += 1;
        strain_level.push((member.read, assembly));
    }
    let mut evidence: HashMap<u32, GenomeEvidence> = HashMap::new();
    fold_candidate_of(
        data,
        level.policy,
        &strain_level,
        (|r| r.0, |_| taxid),
        |&(_, assembly), candidate| candidate.and_then(|c| c.placement_on(assembly)).copied(),
        |_, placement| {
            if let Some(placement) = placement {
                add_placement(&mut evidence, taxid, &placement, true);
            }
        },
    );
    let tests = level.tests(evidence.len().max(1));
    for (&assembly, genome) in evidence.iter_mut() {
        let coverage = tests.measure(assembly, &mut genome.assembly);
        genome.kept = tests.covers(&coverage);
        genome.coverage = Some(coverage);
    }
    let total_weight = share.em.total_weight().max(1e-12);
    let mut listed = resolved.keys().copied().collect::<Vec<_>>();
    if !resolved.contains_key(&share.best_fit) {
        listed.push(share.best_fit);
    }
    listed.sort_unstable();
    for assembly in listed {
        let (resolved, gated) = resolved.get(&assembly).copied().unwrap_or((0, 0));
        let (coverage, kept) = evidence
            .get(&assembly)
            .map_or((None, false), |genome| (genome.coverage, genome.kept));
        diagnostics.push(StrainDiagnostic {
            taxid,
            assembly_id: assembly,
            best_fit: assembly == share.best_fit,
            em_share: share.em.weight.get(&assembly).copied().unwrap_or(0.0) / total_weight,
            resolved,
            gated,
            coverage,
            kept,
        });
    }
    // Second genome EM over the kept genomes, with the fragments that match at
    // least one of them at strain-level identity; the rest stay with the
    // species. A genome whose fragments do not match it within sequencing
    // error is dropped and the EM rerun over the rest.
    let mut kept: HashMap<u32, (EvidenceCoverage, f64)> = evidence
        .into_iter()
        .filter(|(_, genome)| genome.kept)
        .map(|(assembly, genome)| {
            let coverage = genome.coverage.expect("measured above");
            (assembly, (coverage, genome.assembly.asm_len))
        })
        .collect();
    loop {
        if kept.is_empty() {
            return (Vec::new(), StrainShares::default());
        }
        let KeptShares {
            em,
            slots,
            reads,
            bases,
            identity,
        } = share_among_kept(members, &kept, strain_identity, level);
        let mean_identity = |assembly: &u32| {
            identity
                .get(assembly)
                .filter(|(_, weight)| *weight > 0.0)
                .map_or(f64::NAN, |(sum, weight)| sum / weight)
        };
        let failing = kept
            .keys()
            .filter(|assembly| {
                let ani = estimated_ani(strain.baseline, mean_identity(assembly));
                ani.is_nan() || ani < GENOME_ROW_ANI
            })
            .copied()
            .collect::<Vec<_>>();
        if !failing.is_empty() {
            for assembly in failing {
                kept.remove(&assembly);
            }
            continue;
        }
        let mut strains = kept
            .iter()
            .map(|(&assembly_id, &(coverage, asm_len))| StrainGenome {
                assembly_id,
                reads: reads.get(&assembly_id).copied().unwrap_or(0.0),
                bases: bases.get(&assembly_id).copied().unwrap_or(0.0),
                mean_identity: mean_identity(&assembly_id),
                coverage,
                asm_len,
            })
            .collect::<Vec<_>>();
        strains.sort_by(|a, b| {
            b.reads
                .total_cmp(&a.reads)
                .then(a.assembly_id.cmp(&b.assembly_id))
        });
        let shares = if level.keep_shares {
            StrainShares {
                members: members
                    .iter()
                    .zip(slots)
                    .filter_map(|(member, slot)| slot.map(|group| (member.read, group)))
                    .collect(),
                posteriors: em.posteriors,
            }
        } else {
            StrainShares::default()
        };
        return (strains, shares);
    }
}

/// The genome EM over the kept genomes of a species, sharing each fragment
/// among the kept genomes it matches at strain-level identity (each member's
/// read group, if any), with each genome's reads and their bases (each
/// fragment counted by the species' share of it times its posterior) and the
/// posterior-weighted sum and weight of its fragments' identity to it.
struct KeptShares {
    em: likelihood_em::EmResult,
    slots: Vec<Option<usize>>,
    reads: HashMap<u32, f64>,
    bases: HashMap<u32, f64>,
    identity: HashMap<u32, (f64, f64)>,
}

fn share_among_kept(
    members: &[GenomeFragment],
    kept: &HashMap<u32, (EvidenceCoverage, f64)>,
    strain_identity: f64,
    level: &GenomeLevel,
) -> KeptShares {
    let mut groups = ReadGroups::default();
    let mut candidates = Vec::new();
    let slots = members
        .iter()
        .map(|member| {
            candidates.clear();
            candidates.extend(
                member
                    .genomes
                    .iter()
                    .filter(|g| kept.contains_key(&g.0) && g.2 >= strain_identity)
                    .map(|g| (g.0, g.1)),
            );
            groups.add(&mut candidates)
        })
        .collect::<Vec<_>>();
    let em = level.em.run(&groups);
    let mut reads: HashMap<u32, f64> = HashMap::new();
    let mut bases: HashMap<u32, f64> = HashMap::new();
    let mut identity: HashMap<u32, (f64, f64)> = HashMap::new();
    for (member, slot) in members.iter().zip(&slots) {
        let Some(slot) = *slot else {
            continue;
        };
        for &(assembly, posterior) in &em.posteriors[slot] {
            let count = member.share * posterior;
            *reads.entry(assembly).or_insert(0.0) += count;
            *bases.entry(assembly).or_insert(0.0) += member.bases * count;
            if let Some(g) = member.genomes.iter().find(|g| g.0 == assembly) {
                let entry = identity.entry(assembly).or_insert((0.0, 0.0));
                entry.0 += posterior * g.2;
                entry.1 += posterior;
            }
        }
    }
    KeptShares {
        em,
        slots,
        reads,
        bases,
        identity,
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_fixtures::{data, row};
    use super::*;
    use crate::profiler::evidence::EvidenceAlignment;

    fn level() -> GenomeLevel {
        GenomeLevel {
            policy: EvidencePolicy {
                mode: EvidenceMode::Single,
                identity: 0.97,
                mapped_len: 50,
                fraction: 0.5,
            },
            error_rate: 0.003,
            em: EmSettings {
                penalty: likelihood_em::error_penalty(0.003),
                max_iterations: 200,
            },
            report_strains: false,
            keep_shares: false,
            chunk_cutoff: 0.5,
            ratio_bounds: (0.75, 0.75),
        }
    }

    /// Every read gets the species-level candidate of its first row's species,
    /// which the species-level EM gives it whole.
    fn fragments(data: &MappingData) -> (Vec<Fragment>, Vec<Option<(u32, f64)>>) {
        data.read_alignments
            .iter()
            .enumerate()
            .map(|(read, bundle)| {
                let taxid = bundle.alignments[0].profile_taxid;
                (
                    Fragment {
                        read,
                        candidates: vec![(taxid, 0.0)],
                        bases: 100.0,
                    },
                    Some((taxid, 1.0)),
                )
            })
            .unzip()
    }

    #[test]
    fn weighted_median_needs_three_values_and_follows_the_weights() {
        assert_eq!(
            weighted_median([(0.99, 1.0), (0.98, 1.0)].into_iter()),
            None
        );
        assert_eq!(
            weighted_median([(0.99, 1.0), (0.98, 1.0), (0.97, 1.0)].into_iter()),
            Some(0.98)
        );
        assert_eq!(
            weighted_median([(0.99, 10.0), (0.98, 1.0), (0.97, 1.0)].into_iter()),
            Some(0.99)
        );
        assert_eq!(
            weighted_median([(0.99, 1.0), (f64::NAN, 5.0), (0.97, 0.0), (0.95, 1.0)].into_iter()),
            None
        );
    }

    #[test]
    fn best_fit_is_the_genome_the_fragments_match_best() {
        // Species 1 has genomes 10, 11 and 12; its reads match genome 10 exactly and genome 11 with two errors;
        // genome 12 is matched by one read only.
        let mut reads = (0..6)
            .map(|_| vec![row(1, 10, 1.0), row(1, 11, 0.98)])
            .collect::<Vec<_>>();
        reads.push(vec![row(1, 12, 0.99)]);
        let data = data(reads);
        let (fragments, species_of) = fragments(&data);
        let result = resolve(&data, &fragments, &species_of, &HashMap::new(), &level());
        let species = &result.species[&1];
        assert_eq!(species.best_fit, 10);
        // Fragments whose best species is species 1 and that align to the best-fit genome: the six exact ones.
        assert!((species.mean_identity - 1.0).abs() < 1e-9);
        assert!(species.strains.is_empty());
    }

    #[test]
    fn identity_counts_fragments_below_the_species_cutoff_whose_best_species_it_is() {
        let data = data(vec![
            // Above the species cutoff of 0.97.
            vec![row(1, 10, 0.99)],
            // Below it, with species 1 still the best match: counted.
            vec![row(1, 10, 0.95)],
            // Species 2 matches better: left out of species 1.
            vec![row(1, 10, 0.93), row(2, 20, 0.99)],
            // Species 1 only on another of its genomes: no identity to the best-fit genome.
            vec![row(1, 11, 0.99)],
        ]);
        let best_fit = HashMap::from([(1, 10), (2, 20)]);
        let identity = identity_to_best_fit(&data, &best_fit, level().policy);
        let (mean, count) = identity[&1];
        assert_eq!(count, 2.0);
        assert!((mean - 0.97).abs() < 1e-9, "{mean}");
        assert_eq!(identity[&2], (0.99, 1.0));
    }

    /// A read of 100 bases at the `i`-th of 300 positions tiling a genome of 10 kb.
    fn tiled(taxid: u32, assembly: u32, identity: f32, i: u32) -> EvidenceAlignment {
        let start = i * 33;
        EvidenceAlignment {
            asm_len: 10_000.0,
            span: (start, start + 100),
            fixed_chunk: (start, start + 100),
            flex_chunk: (start, start + 100),
            ..row(taxid, assembly, identity)
        }
    }

    #[test]
    fn a_genome_whose_fragments_match_it_only_as_a_relative_is_not_reported() {
        // Strain A is genome 10. Strain B is a relative of genome 11: its reads match genome 11 better than genome 10
        // but still with errors beyond sequencing error.
        let mut reads = (0..300)
            .map(|i| vec![tiled(1, 10, 0.998, i), tiled(1, 11, 0.990, i)])
            .collect::<Vec<_>>();
        reads.extend((0..300).map(|i| vec![tiled(1, 11, 0.995, i), tiled(1, 10, 0.985, i)]));
        let data = data(reads);
        let (fragments, species_of) = fragments(&data);
        let level = GenomeLevel {
            report_strains: true,
            chunk_cutoff: 0.1,
            ratio_bounds: (0.5, 0.5),
            ..level()
        };
        let result = resolve(&data, &fragments, &species_of, &HashMap::new(), &level);
        let strains = &result.species[&1].strains;
        // Genome 11 passes the coverage tests but its fragments read about 99.8 % to it: dropped. Strain B's reads
        // do not match genome 10 at strain-level identity, so they stay with the species instead of going to 10.
        assert_eq!(strains.len(), 1);
        assert_eq!(strains[0].assembly_id, 10);
        assert!(
            (strains[0].reads - 300.0).abs() < 1e-6,
            "{}",
            strains[0].reads
        );
        assert!((strains[0].mean_identity - 0.998).abs() < 1e-9);
        let kept_by_coverage = result
            .diagnostics
            .iter()
            .filter(|d| d.kept)
            .map(|d| d.assembly_id)
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(kept_by_coverage, [10, 11].into());
    }

    #[test]
    fn genome_rows_are_measured_against_a_baseline_cut_off_like_them() {
        // Species 2, 3 and 4 are exact strains of genomes 20, 30 and 40, with a tail of reads carrying more
        // sequencing errors. Species 1 is a strain related to genome 10: half its reads match it closely, half
        // poorly. Baseline 0.9981 (all reads), strain-level identity 0.9931, strain baseline 0.999.
        let mut reads = Vec::new();
        for (taxid, assembly) in [(2, 20), (3, 30), (4, 40)] {
            reads.extend((0..270).map(|i| vec![tiled(taxid, assembly, 0.999, i)]));
            reads.extend((270..300).map(|i| vec![tiled(taxid, assembly, 0.990, i)]));
        }
        reads.extend((0..300).map(|i| vec![tiled(1, 10, 0.9975, i)]));
        reads.extend((0..300).map(|i| vec![tiled(1, 10, 0.985, i)]));
        let data = data(reads);
        let (fragments, species_of) = fragments(&data);
        let level = GenomeLevel {
            report_strains: true,
            chunk_cutoff: 0.1,
            ratio_bounds: (0.5, 0.5),
            ..level()
        };
        let result = resolve(&data, &fragments, &species_of, &HashMap::new(), &level);
        assert!((result.baseline_identity - 0.9981).abs() < 1e-9);
        assert!((result.strain_baseline_identity - 0.999).abs() < 1e-9);
        // Genome 10's strain-level fragments read 0.9975: within 0.001 of the all-read baseline, but 0.0015 below
        // the strain baseline they are comparable with.
        assert!(result.species[&1].strains.is_empty());
        for taxid in [2, 3, 4] {
            assert_eq!(result.species[&taxid].strains.len(), 1, "species {taxid}");
        }
    }

    #[test]
    fn long_reads_use_the_same_strain_level_margin() {
        // An exact strain of genome 10 and, as many reads, a strain at 99.0 % to it. Baseline 0.997, so strain-level
        // identity is 0.992: the relative's reads at 0.989 must not count as matches of genome 10 at strain level.
        let mut reads = (0..300)
            .map(|i| vec![tiled(1, 10, 0.999, i)])
            .collect::<Vec<_>>();
        reads.extend((0..300).map(|i| vec![tiled(1, 10, 0.989, i)]));
        let data = data(reads);
        let (fragments, species_of) = fragments(&data);
        let mut level = GenomeLevel {
            report_strains: true,
            chunk_cutoff: 0.1,
            ratio_bounds: (0.5, 0.5),
            ..level()
        };
        level.policy.mode = EvidenceMode::Long;
        let result = resolve(&data, &fragments, &species_of, &HashMap::new(), &level);
        let strains = &result.species[&1].strains;
        assert_eq!(strains.len(), 1);
        assert!((strains[0].mean_identity - 0.999).abs() < 1e-9);
        assert!((strains[0].reads - 300.0).abs() < 1e-6);
    }

    #[test]
    fn a_genome_row_counts_its_fragments_by_their_species_share() {
        // Species 1 is an exact strain of genome 10. The species-level EM gives it 0.8 of the first 150 of its 300
        // fragments (another species holds the rest of them) and the other 150 whole.
        let reads = (0..300)
            .map(|i| vec![tiled(1, 10, 0.999, i)])
            .collect::<Vec<_>>();
        let data = data(reads);
        let (fragments, mut species_of) = fragments(&data);
        for species in &mut species_of[..150] {
            *species = Some((1, 0.8));
        }
        let level = GenomeLevel {
            report_strains: true,
            chunk_cutoff: 0.1,
            ratio_bounds: (0.5, 0.5),
            ..level()
        };
        let result = resolve(&data, &fragments, &species_of, &HashMap::new(), &level);
        let strains = &result.species[&1].strains;
        assert_eq!(strains.len(), 1);
        // The genome row counts each fragment by the species' share of it.
        assert!((strains[0].reads - 270.0).abs() < 1e-9);
        assert!((strains[0].bases - 27_000.0).abs() < 1e-6);
        // The species also holds half of ten fragments that the genome level does not resolve: its species row keeps
        // those 5 reads, and its rows together hold its 275 reads. The other species keep their shares of the
        // fragments, which the genome level does not touch.
        let species_reads = 270.0 + 10.0 * 0.5;
        let (species_row, _) = super::super::species_remainder(
            species_reads,
            27_500.0,
            strains.iter().map(|s| (s.reads, s.bases)),
        )
        .unwrap();
        assert!((species_row - 5.0).abs() < 1e-9);
        assert!((species_row + strains[0].reads - species_reads).abs() < 1e-9);
    }
}

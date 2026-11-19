//! Genomes within each kept species.
//!
//! The fragments that the species-level EM gives to a species are shared
//! among the species' genomes by the same match-quality EM. The genome with the
//! largest share is the species' best-fit genome, and the mean identity of the
//! species' fragments to it estimates how close the sample's strain is to it.
//! At the genome level (`--strain`) a genome is also reported on its own when
//! the fragments that clearly belong to it and match it at strain-level
//! identity cover it evenly. A genome represented by marker genes never is:
//! its genes are conserved and too short to tell strains apart.

use std::collections::HashMap;

use rayon::prelude::*;

use super::{add_placement, fold_candidates, within, Fragment, GenomeEvidence, RESOLVE_POSTERIOR};
use crate::profiler::community::MappingData;
use crate::profiler::evidence::{
    calculate_coverage, EvidenceCoverage, EvidenceMode, EvidencePolicy, TerminalCandidate,
};
use crate::profiler::likelihood_em::{self, ReadGroups};
use crate::profiler::ProfileRank;

/// Genome candidates of a short fragment whose errors exceed the best by
/// less than this count as tied.
const SHORT_TIE_ERRORS: f64 = 0.5;
/// Long reads: the same as a share of the best genome's aligned bases.
const LONG_TIE_SHARE: f64 = 0.001;
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

/// One fragment of a species: its read, mapped bases for depth, and per
/// genome of the species `(assembly, errors, identity)`.
struct GenomeFragment {
    read: usize,
    bases: f64,
    genomes: Vec<(u32, f64, f64)>,
}

/// A genome reported on its own.
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
}

/// One genome's genome-level tests (development diagnostics).
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
struct StrainLevel {
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
    mean_identity: f64,
}

/// Settings shared by every species.
pub(super) struct GenomeLevel {
    pub(super) policy: EvidencePolicy,
    pub(super) error_rate: f64,
    pub(super) penalty: f64,
    /// Report genomes on their own (`--strain`).
    pub(super) report_strains: bool,
    pub(super) chunk_cutoff: f64,
    pub(super) ratio_bounds: (f64, f64),
    pub(super) max_iterations: usize,
    pub(super) tolerance: f64,
}

/// For each species of `species_of`, share its fragments among its genomes;
/// then measure the sample's baseline identity and, with `report_strains`,
/// find the genomes reported on their own.
pub(super) fn resolve(
    data: &MappingData,
    fragments: &[Fragment],
    species_of: &[Option<u32>],
    level: &GenomeLevel,
) -> GenomeResult {
    let mut shares = genome_fragments(data, fragments, species_of, level.policy)
        .into_iter()
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|(taxid, members)| {
            let mut groups = ReadGroups::default();
            let group_of = members
                .iter()
                .map(|m| groups.add(m.genomes.iter().map(|g| (g.0, g.1)).collect(), 1.0))
                .collect::<Vec<_>>();
            let em = likelihood_em::run(
                &groups,
                &HashMap::new(),
                level.penalty,
                level.max_iterations,
                level.tolerance,
            );
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
                mean_identity: f64::NAN,
            }
        })
        .collect::<Vec<_>>();
    let best_fit = shares.iter().map(|s| (s.taxid, s.best_fit)).collect();
    let identity = identity_to_best_fit(data, &best_fit, level.policy);
    for share in &mut shares {
        share.mean_identity = identity.get(&share.taxid).map_or(f64::NAN, |m| m.0);
    }
    let baseline_identity =
        weighted_median(identity.values().copied()).unwrap_or(1.0 - level.error_rate);
    let strain_identity = baseline_identity - STRAIN_MARGIN;
    let strain = StrainLevel {
        identity: strain_identity,
        baseline: weighted_median(
            shares
                .iter()
                .map(|s| strain_level_identity(&s.members, s.best_fit, strain_identity)),
        )
        .unwrap_or(baseline_identity),
    };
    let mut diagnostics = Vec::new();
    let species = shares
        .into_iter()
        .map(|share| {
            let strains = if level.report_strains {
                strain_genomes(data, &share, strain, level, &mut diagnostics)
            } else {
                Vec::new()
            };
            (
                share.taxid,
                SpeciesGenomes {
                    best_fit: share.best_fit,
                    mean_identity: share.mean_identity,
                    strains,
                },
            )
        })
        .collect();
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

/// Median of `(value, weight)` pairs with finite values and positive weight;
/// `None` for fewer than three.
fn weighted_median(values: impl Iterator<Item = (f64, f64)>) -> Option<f64> {
    let mut values = values
        .filter(|v| v.0.is_finite() && v.1 > 0.0)
        .collect::<Vec<_>>();
    if values.len() < 3 {
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
    species_of: &[Option<u32>],
    policy: EvidencePolicy,
) -> HashMap<u32, Vec<GenomeFragment>> {
    let resolved = fragments
        .iter()
        .zip(species_of)
        .filter_map(|(fragment, species)| species.map(|taxid| (fragment, taxid)))
        .collect::<Vec<_>>();
    let mut out: HashMap<u32, Vec<GenomeFragment>> = HashMap::new();
    fold_candidates(
        data,
        policy,
        &resolved,
        |r| r.0.read,
        |r, candidates| {
            let candidate = candidates.iter().find(|c| c.taxid == r.1)?;
            Some(candidate_genomes(candidate, policy))
        },
        |&(fragment, taxid), genomes| {
            if let Some(genomes) = genomes {
                out.entry(taxid).or_default().push(GenomeFragment {
                    read: fragment.read,
                    bases: fragment.bases,
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
    let longest = candidate
        .placements
        .iter()
        .map(|p| f64::from(p.both_mapped))
        .fold(0.0, f64::max);
    let mut genomes = candidate
        .placements
        .iter()
        .map(|p| {
            let mapped = f64::from(p.both_mapped).max(1.0);
            let errors = p.pair_errors
                + super::UNALIGNED_ERROR_WEIGHT * (longest - f64::from(p.both_mapped)).max(0.0);
            (p.assembly_id, errors, 1.0 - p.pair_errors / mapped)
        })
        .collect::<Vec<_>>();
    let fewest = genomes.iter().map(|g| g.1).fold(f64::INFINITY, f64::min);
    let margin = match policy.mode {
        EvidenceMode::Long => LONG_TIE_SHARE * longest,
        EvidenceMode::Paired | EvidenceMode::Single => SHORT_TIE_ERRORS,
    };
    for genome in &mut genomes {
        if genome.1 - fewest < margin {
            genome.1 = fewest;
        }
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
        collect_exact_ties: false,
        collect_identities: false,
        collect_best_placements: false,
        ..policy
    };
    let reads = (0..data.read_alignments.len()).collect::<Vec<_>>();
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

/// A fragment's identity to the best-fit genome of each species in
/// `best_fit` that is its best species match (fewest errors, within the
/// species-level tie margin), in candidate order.
fn best_fit_identities(
    candidates: &[TerminalCandidate],
    best_fit: &HashMap<u32, u32>,
    mode: EvidenceMode,
) -> Vec<(u32, f64)> {
    let longest = candidates
        .iter()
        .map(|c| f64::from(c.best.both_mapped))
        .fold(0.0, f64::max);
    let species = candidates
        .iter()
        .filter(|c| c.profile_rank == ProfileRank::Species)
        .map(|c| {
            let errors = c.best.pair_errors
                + super::UNALIGNED_ERROR_WEIGHT
                    * (longest - f64::from(c.best.both_mapped)).max(0.0);
            (c, errors)
        })
        .collect::<Vec<_>>();
    let Some(&(best, fewest)) = species.iter().min_by(|a, b| a.1.total_cmp(&b.1)) else {
        return Vec::new();
    };
    let margin = match mode {
        EvidenceMode::Long => super::LONG_TIE_SHARE * f64::from(best.best.both_mapped),
        EvidenceMode::Paired | EvidenceMode::Single => super::SHORT_TIE_ERRORS,
    };
    species
        .iter()
        .filter(|&&(_, errors)| errors - fewest < margin)
        .filter_map(|&(candidate, _)| {
            let assembly = *best_fit.get(&candidate.taxid)?;
            let placement = candidate
                .placements
                .iter()
                .find(|p| p.assembly_id == assembly)?;
            Some((
                candidate.taxid,
                1.0 - placement.pair_errors / f64::from(placement.both_mapped).max(1.0),
            ))
        })
        .collect()
}

/// Genomes of one species reported on their own: strain-level coverage from
/// fragments that resolve to a genome (posterior >= RESOLVE_POSTERIOR) at
/// strain-level identity, the chunk-breadth cutoff and breadth ratios, and a
/// second genome EM over the kept genomes for their reads and depth.
fn strain_genomes(
    data: &MappingData,
    share: &SpeciesShare,
    strain: StrainLevel,
    level: &GenomeLevel,
    diagnostics: &mut Vec<StrainDiagnostic>,
) -> Vec<StrainGenome> {
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
        let Some(assembly) = group.and_then(|g| {
            posteriors[g]
                .iter()
                .find(|p| p.1 >= RESOLVE_POSTERIOR)
                .map(|p| p.0)
        }) else {
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
    fold_candidates(
        data,
        level.policy,
        &strain_level,
        |r| r.0,
        |&(_, assembly), candidates| {
            candidates
                .iter()
                .find(|c| c.taxid == taxid)
                .and_then(|c| c.placements.iter().find(|p| p.assembly_id == assembly))
                .copied()
        },
        |_, placement| {
            if let Some(placement) = placement {
                add_placement(&mut evidence, taxid, &placement, true);
            }
        },
    );
    let tests = evidence.len().max(1);
    let total_weight: f64 = share.em.weight.values().sum::<f64>().max(1e-12);
    let mut kept: HashMap<u32, (EvidenceCoverage, f64)> = HashMap::new();
    let mut measured: HashMap<u32, (EvidenceCoverage, bool)> = HashMap::new();
    for (&assembly, genome) in evidence.iter_mut() {
        let coverage = calculate_coverage(&mut genome.assembly, tests);
        let passes = coverage.fixed_chunk > level.chunk_cutoff
            && within(coverage.oebr, level.ratio_bounds.0)
            && within(coverage.coebr, level.ratio_bounds.1);
        if passes {
            kept.insert(assembly, (coverage, genome.assembly.asm_len));
        }
        measured.insert(assembly, (coverage, passes));
    }
    let mut listed = resolved.keys().copied().collect::<Vec<_>>();
    if !resolved.contains_key(&share.best_fit) {
        listed.push(share.best_fit);
    }
    listed.sort_unstable();
    for assembly in listed {
        let (resolved, gated) = resolved.get(&assembly).copied().unwrap_or((0, 0));
        let (coverage, kept) = measured
            .remove(&assembly)
            .map_or((None, false), |(c, k)| (Some(c), k));
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
    loop {
        if kept.is_empty() {
            return Vec::new();
        }
        let KeptShares {
            em,
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
                let ani = 1.0 - (strain.baseline - mean_identity(assembly));
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
                reads: em.weight.get(&assembly_id).copied().unwrap_or(0.0),
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
        return strains;
    }
}

/// The genome EM over the kept genomes of a species, sharing each fragment
/// among the kept genomes it matches at strain-level identity, with each
/// genome's fragment bases and the posterior-weighted sum and weight of its
/// fragments' identity to it.
struct KeptShares {
    em: likelihood_em::EmResult,
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
    let slots = members
        .iter()
        .map(|member| {
            let candidates = member
                .genomes
                .iter()
                .filter(|g| kept.contains_key(&g.0) && g.2 >= strain_identity)
                .map(|g| (g.0, g.1))
                .collect::<Vec<_>>();
            groups.add(candidates, 1.0)
        })
        .collect::<Vec<_>>();
    let em = likelihood_em::run(
        &groups,
        &HashMap::new(),
        level.penalty,
        level.max_iterations,
        level.tolerance,
    );
    let mut bases: HashMap<u32, f64> = HashMap::new();
    let mut identity: HashMap<u32, (f64, f64)> = HashMap::new();
    for (member, slot) in members.iter().zip(&slots) {
        let Some(slot) = *slot else {
            continue;
        };
        for &(assembly, posterior) in &em.posteriors[slot] {
            *bases.entry(assembly).or_insert(0.0) += member.bases * posterior;
            if let Some(g) = member.genomes.iter().find(|g| g.0 == assembly) {
                let entry = identity.entry(assembly).or_insert((0.0, 0.0));
                entry.0 += posterior * g.2;
                entry.1 += posterior;
            }
        }
    }
    KeptShares {
        em,
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
                collect_exact_ties: false,
                collect_identities: false,
                collect_best_placements: false,
            },
            error_rate: 0.003,
            penalty: likelihood_em::error_penalty(0.003),
            report_strains: false,
            chunk_cutoff: 0.5,
            ratio_bounds: (0.75, 0.75),
            max_iterations: 200,
            tolerance: 1e-6,
        }
    }

    /// Every read gets the species-level candidate of its first row's species.
    fn fragments(data: &MappingData) -> (Vec<Fragment>, Vec<Option<u32>>) {
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
                    Some(taxid),
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
        let result = resolve(&data, &fragments, &species_of, &level());
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
        let result = resolve(&data, &fragments, &species_of, &level);
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
        let result = resolve(&data, &fragments, &species_of, &level);
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
        let result = resolve(&data, &fragments, &species_of, &level);
        let strains = &result.species[&1].strains;
        assert_eq!(strains.len(), 1);
        assert!((strains[0].mean_identity - 0.999).abs() < 1e-9);
        assert!((strains[0].reads - 300.0).abs() < 1e-6);
    }
}

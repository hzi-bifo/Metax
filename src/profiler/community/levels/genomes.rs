//! Genomes within each kept species.
//!
//! The fragments that the species-level EM gives to a species are shared
//! among the species' genomes by the same match-quality EM. The genome with the
//! largest share is the species' best-fit genome, and the mean identity of the
//! species' fragments to it estimates how close the sample's strain is to it.
//! At the genome level (`--strain`) a genome is also reported on its own when
//! the fragments that clearly belong to it and match it at strain-level
//! identity cover it evenly.

use std::collections::HashMap;

use super::{add_placement, within, Fragment, GenomeEvidence, RESOLVE_POSTERIOR};
use crate::profiler::community::MappingData;
use crate::profiler::evidence::{
    calculate_coverage, CandidateScratch, EvidenceCoverage, EvidenceMode, EvidencePolicy,
};
use crate::profiler::likelihood_em::{self, ReadGroups};

/// Genome candidates of a short fragment whose errors exceed the best by
/// less than this count as tied.
const SHORT_TIE_ERRORS: f64 = 0.5;
/// Long reads: the same as a share of the best genome's aligned bases.
const LONG_TIE_SHARE: f64 = 0.001;
/// A fragment counts for a genome's strain-level coverage when it aligns to
/// the genome at an identity of at least the sample's identity to its own
/// genomes less this margin.
const SHORT_STRAIN_MARGIN: f64 = 0.005;
const LONG_STRAIN_MARGIN: f64 = 0.01;
/// A genome is the sample's strain itself at an estimated ANI of at least this.
pub(super) const SAME_STRAIN_ANI: f64 = 0.995;

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
    /// Mean identity of the species' fragments to the best-fit genome.
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
    /// best-fit genome.
    pub(super) baseline_identity: f64,
    pub(super) diagnostics: Vec<StrainDiagnostic>,
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
    let shares = genome_fragments(data, fragments, species_of, level.policy)
        .into_iter()
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
            let mean_identity = mean_identity_on(&members, best_fit);
            SpeciesShare {
                taxid,
                members,
                group_of,
                em,
                best_fit,
                mean_identity,
            }
        })
        .collect::<Vec<_>>();
    let baseline_identity = weighted_median(
        shares
            .iter()
            .map(|s| (s.mean_identity, s.members.len() as f64)),
    )
    .unwrap_or(1.0 - level.error_rate);
    let mut diagnostics = Vec::new();
    let mut scratch = CandidateScratch::default();
    let species = shares
        .into_iter()
        .map(|share| {
            let strains = if level.report_strains {
                strain_genomes(
                    data,
                    &share,
                    baseline_identity,
                    level,
                    &mut scratch,
                    &mut diagnostics,
                )
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
        diagnostics,
    }
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
    let mut scratch = CandidateScratch::default();
    let mut out: HashMap<u32, Vec<GenomeFragment>> = HashMap::new();
    for (fragment, species) in fragments.iter().zip(species_of) {
        let Some(taxid) = *species else {
            continue;
        };
        let candidates = scratch.build(&data.read_alignments[fragment.read].alignments, policy);
        let Some(candidate) = candidates.iter().find(|c| c.taxid == taxid) else {
            continue;
        };
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
        out.entry(taxid).or_default().push(GenomeFragment {
            read: fragment.read,
            bases: fragment.bases,
            genomes,
        });
    }
    out
}

/// Mean identity of the fragments that align to `assembly`.
fn mean_identity_on(members: &[GenomeFragment], assembly: u32) -> f64 {
    let (sum, count) = members
        .iter()
        .filter_map(|m| m.genomes.iter().find(|g| g.0 == assembly).map(|g| g.2))
        .fold((0.0, 0usize), |(s, n), x| (s + x, n + 1));
    if count > 0 {
        sum / count as f64
    } else {
        f64::NAN
    }
}

/// Genomes of one species reported on their own: strain-level coverage from
/// fragments that resolve to a genome (posterior >= RESOLVE_POSTERIOR) at
/// strain-level identity, the chunk-breadth cutoff and breadth ratios, and a
/// second genome EM over the kept genomes for their reads and depth.
fn strain_genomes(
    data: &MappingData,
    share: &SpeciesShare,
    baseline_identity: f64,
    level: &GenomeLevel,
    scratch: &mut CandidateScratch,
    diagnostics: &mut Vec<StrainDiagnostic>,
) -> Vec<StrainGenome> {
    let (taxid, members, group_of, posteriors) = (
        share.taxid,
        &share.members,
        &share.group_of,
        &share.em.posteriors,
    );
    let margin = match level.policy.mode {
        EvidenceMode::Long => LONG_STRAIN_MARGIN,
        EvidenceMode::Paired | EvidenceMode::Single => SHORT_STRAIN_MARGIN,
    };
    let strain_identity = baseline_identity - margin;
    let mut evidence: HashMap<u32, GenomeEvidence> = HashMap::new();
    let mut resolved: HashMap<u32, (u64, u64)> = HashMap::new();
    for (member, group) in members.iter().zip(group_of) {
        let Some(assembly) = group.and_then(|g| {
            posteriors[g]
                .iter()
                .find(|p| p.1 >= RESOLVE_POSTERIOR)
                .map(|p| p.0)
        }) else {
            continue;
        };
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
        let candidates = scratch.build(&data.read_alignments[member.read].alignments, level.policy);
        if let Some(placement) = candidates
            .iter()
            .find(|c| c.taxid == taxid)
            .and_then(|c| c.placements.iter().find(|p| p.assembly_id == assembly))
        {
            add_placement(&mut evidence, taxid, placement, true);
        }
    }
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
    if kept.is_empty() {
        return Vec::new();
    }
    // Second genome EM over the kept genomes.
    let mut groups = ReadGroups::default();
    let slots = members
        .iter()
        .map(|member| {
            let candidates = member
                .genomes
                .iter()
                .filter(|g| kept.contains_key(&g.0))
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
    let mut strains = kept
        .into_iter()
        .map(|(assembly_id, (coverage, asm_len))| {
            let (sum, weight) = identity.get(&assembly_id).copied().unwrap_or((0.0, 0.0));
            StrainGenome {
                assembly_id,
                reads: em.weight.get(&assembly_id).copied().unwrap_or(0.0),
                bases: bases.get(&assembly_id).copied().unwrap_or(0.0),
                mean_identity: if weight > 0.0 { sum / weight } else { f64::NAN },
                coverage,
                asm_len,
            }
        })
        .collect::<Vec<_>>();
    strains.sort_by(|a, b| {
        b.reads
            .total_cmp(&a.reads)
            .then(a.assembly_id.cmp(&b.assembly_id))
    });
    strains
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use metam::mapper::CandidateBoundary;

    use super::*;
    use crate::profiler::community::ReadBundle;
    use crate::profiler::evidence::{EvidenceAlignment, MappingRole};
    use crate::profiler::ProfileRank;

    /// A single-end read of 100 bases aligned over its whole length.
    fn row(taxid: u32, assembly: u32, identity: f32) -> EvidenceAlignment {
        EvidenceAlignment {
            asm: Arc::from(format!("asm{assembly}")),
            acc: Arc::from(format!("ref{assembly}")),
            profile_taxid: taxid,
            profile_rank: ProfileRank::Species,
            reference_id: assembly,
            assembly_id: assembly,
            role: MappingRole::Unpaired,
            pair_id: 0,
            candidate_boundary: CandidateBoundary::COMPLETE,
            sampling_scale: 1.0,
            asm_len: 1_000_000.0,
            identity,
            fraction: 1.0,
            mapped_len: 100,
            read_len: 100,
            span: (100, 200),
            fixed_chunk: (100, 200),
            flex_chunk: (100, 200),
            marker_gene_len: 0,
        }
    }

    fn data(reads: Vec<Vec<EvidenceAlignment>>) -> MappingData {
        MappingData {
            read_alignments: reads
                .into_iter()
                .enumerate()
                .map(|(i, alignments)| ReadBundle {
                    read_name: Arc::from(format!("read{i}")),
                    alignments,
                })
                .collect(),
            asm_lengths: Arc::new(HashMap::new()),
            min_breadth: 0.0,
            min_chunk_breadth: 0.0,
            safeguard_min_chunk_breadth: 0.0,
            has_scaled_genomes: false,
            total_reads: 0,
            aligned_reads: 0,
            mapping_rate: 1.0,
            fallback_mode: false,
            species_to_genus: Arc::new(HashMap::new()),
            terminal_ranks: Arc::new(HashMap::new()),
            skipped_reference_at_load: 0,
            marker_model: None,
            genus_rollup: None,
            unclassified_species: None,
            assembly_names: Vec::new(),
        }
    }

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
        // Fragments that align to the best-fit genome: the six exact ones.
        assert!((species.mean_identity - 1.0).abs() < 1e-9);
        assert!(species.strains.is_empty());
    }
}

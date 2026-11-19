//! Level-wise profiling engine (development; selected with `--level-engine`).
//!
//! Each fragment's candidates are the species whose genomes it aligns to
//! within the read-level gates; candidates whose estimated errors exceed the
//! fragment's best by less than a margin count as tied. A first EM, weighted
//! by how well each fragment matches each species ([`likelihood_em`]), gives a
//! fragment to a species when its posterior there is at least 0.9; such a
//! fragment supports every genome of that species it aligns to. A species is
//! kept when one of its genomes shows enough evenly spread coverage and the
//! species wins enough of the fragments it could explain, and is reported on
//! the kept genome where most of its fragments fit best. Removed species are
//! then judged again, strongest first, against the species kept so far, so
//! that a removed sister no longer splits the fragments it ties on. Fragments
//! with no candidate at the species identity are offered to their families.
//! A second EM shares
//! the fragments among the kept species, and abundance is the resulting depth
//! on each species' best-fit genome; a fragment on a marker genome counts as
//! one on the whole genome ([`fragment_depth`]). Read counts are given as full-
//! index equivalents, scaled up for genomes indexed in sampled segments. A genome represented by marker genes is
//! measured gene by gene ([`MarkerModel`]); the tests are the same. With
//! `--save-raw-profile` the engine also writes `<prefix>.levels_genomes.tsv`,
//! the tests of every measured genome.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use anyhow::Result;
use rayon::prelude::*;

use super::{add_marker_genome_reads, fragment_depth, CommunityProfiler, MappingData};
use crate::profiler::evidence::{
    calculate_coverage, pvalue_passes, standard_ratio_lower_bounds, witness_is_better,
    CandidatePlacement, CandidateScratch, EvidenceCoverage, EvidenceMode, EvidencePolicy, GeneHit,
    PremapChunkGuard, PresenceAssembly, PresenceInterval, RescueKinds, TerminalCandidate,
    WitnessDecision, STANDARD_RATIO_UPPER,
};
use crate::profiler::likelihood_em::{self, ReadGroups};
use crate::profiler::markers::MarkerModel;
use crate::profiler::output::{BestFitReport, ProfileRow};
use crate::profiler::rollup::{
    GenomeLayout, GenusRollup, GenusRow, OrphanedFragments, MIN_GENUS_FRAGMENTS, MIN_WINDOW_SPREAD,
};
use crate::profiler::unclassified::GenomeReads;
use crate::profiler::{parse_reference_taxonomy, ProfileRank};
use crate::taxonomy::Taxonomy;

mod genomes;

/// Each base that a candidate aligns less than the read's longest alignment
/// counts as this many errors in the EM.
const UNALIGNED_ERROR_WEIGHT: f64 = 0.5;
/// Candidates whose likelihood is below this fraction of the read's best
/// candidate are dropped before the EM.
const MIN_RELATIVE_LIKELIHOOD: f64 = 1e-6;
/// Candidates whose errors exceed the read's best by less than this margin
/// count as tied: short reads by 0.5 estimated errors, long reads by this
/// share of the best candidate's aligned bases. Smaller differences are as
/// likely to come from a strain the database does not hold as from the
/// candidate itself.
const SHORT_TIE_ERRORS: f64 = 0.5;
const LONG_TIE_SHARE: f64 = 0.005;
/// A fragment resolves to a species, and adds to the coverage of that
/// species' genomes, when the first EM gives the species at least this
/// posterior.
const RESOLVE_POSTERIOR: f64 = 0.9;
/// A species is kept only when it wins (posterior >= RESOLVE_POSTERIOR) at
/// least this share of the fragments that list it among their candidates: a
/// species that only collects reads which a present relative explains wins
/// few of them, while a present species wins most of its reads even when its
/// strain differs from every genome of it in the database.
const MIN_WON_SHARE: f64 = 0.3;
/// Share of the read groups' weight that may still change when the EM stops.
const EM_TOLERANCE: f64 = 1e-7;
/// The sample's error rate is measured on this share of the fragments that
/// match their best candidate best.
const ERROR_RATE_SHARE: f64 = 0.75;
/// The candidate passes build the candidates of this many fragments at once on
/// the profiler's threads, then fold the results in fragment order, so that
/// interval lists and sums come out as in a single pass.
const PARALLEL_BLOCK: usize = 1 << 16;

/// A fragment's candidate species with their errors (estimated errors of the
/// species' best placement plus the unaligned-base penalty), the read it comes
/// from, and its mapped bases (first mate of a pair) for depth.
struct Fragment {
    read: usize,
    candidates: Vec<(u32, f64)>,
    bases: f64,
}

/// Fragments shared among a set of species: the EM's read groups, each
/// fragment's group, and each group's bases.
struct SharedReads {
    groups: ReadGroups,
    group_of: Vec<Option<usize>>,
    group_bases: Vec<f64>,
}

/// Coverage of one genome from the fragments that resolve to its species.
#[derive(Clone)]
struct GenomeEvidence {
    taxid: u32,
    assembly: PresenceAssembly,
    coverage: Option<EvidenceCoverage>,
    kept: bool,
}

pub(super) fn profile(profiler: &CommunityProfiler, data: MappingData) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(profiler.cfg.threads.max(1))
        .build()?
        .install(|| profile_levels(profiler, data))
}

fn profile_levels(profiler: &CommunityProfiler, data: MappingData) -> Result<()> {
    let cfg = &profiler.cfg;
    let policy = EvidencePolicy::from_config(cfg).ok_or_else(|| {
        anyhow::anyhow!(
            "the level-wise engine takes Illumina reads or long reads without --strain for now"
        )
    })?;
    let mode = policy.mode;
    let (fragments, unmatched, error_rate) = collect_fragments(&data, policy);
    let penalty = likelihood_em::error_penalty(error_rate);
    log::info!(
        target: "PROFILER",
        "Level-wise engine: {} fragments with candidates; error rate {:.5}, error penalty {:.3}.",
        fragments.len(),
        error_rate,
        penalty
    );

    // First EM over every candidate species; a fragment that it gives to one
    // species with posterior >= RESOLVE_POSTERIOR adds to the coverage of
    // that species' genomes.
    let all = share_reads(&fragments, |_| true, penalty);
    let first = likelihood_em::run(
        &all.groups,
        &HashMap::new(),
        penalty,
        cfg.max_em_iterations,
        EM_TOLERANCE,
    );
    let resolved_taxid = first
        .posteriors
        .iter()
        .map(|posteriors| {
            posteriors
                .iter()
                .find(|p| p.1 >= RESOLVE_POSTERIOR)
                .map(|p| p.0)
        })
        .collect::<Vec<_>>();
    // Fragments listing each species among their candidates, and those resolved
    // to it.
    let mut counts = FragmentCounts::default();
    for &group in all.group_of.iter().flatten() {
        let resolved = resolved_taxid[group];
        for &(taxid, _) in &first.posteriors[group] {
            *counts.candidates.entry(taxid).or_default() += 1;
            if resolved.is_some() {
                *counts.resolved_candidates.entry(taxid).or_default() += 1;
            }
        }
        if let Some(taxid) = resolved {
            *counts.won.entry(taxid).or_default() += 1;
        }
    }
    let resolved_fragments = fragments
        .iter()
        .zip(&all.group_of)
        .filter_map(|(fragment, group)| {
            group
                .and_then(|g| resolved_taxid[g])
                .map(|taxid| (fragment.read, taxid))
        })
        .collect::<Vec<_>>();
    let mut genomes: HashMap<u32, GenomeEvidence> = HashMap::new();
    fold_candidates(
        &data,
        policy,
        &resolved_fragments,
        |r| r.0,
        |r, candidates| species_placements(candidates, r.1),
        |r, hits| add_hits(&mut genomes, r.1, hits),
    );

    let (oebr_default, coebr_default) = standard_ratio_lower_bounds(ProfileRank::Species, mode);
    let min_oebr = cfg.min_oebr.unwrap_or(oebr_default);
    let min_coebr = cfg.min_coebr.unwrap_or(coebr_default);
    let tests = SpeciesTests {
        multiple_tests: genomes.len().max(1),
        chunk_cutoff: data.min_chunk_breadth,
        ratio_bounds: (min_oebr, min_coebr),
        markers: data.marker_model.clone(),
        min_marker_genes: cfg.min_marker_genes,
    };
    genomes.par_iter_mut().for_each(|(&id, genome)| {
        let won = counts.won.get(&genome.taxid).copied().unwrap_or(0);
        let candidates = counts.candidates.get(&genome.taxid).copied().unwrap_or(0);
        tests.judge(id, genome, won, candidates);
    });
    let (second_kept, second_taken) = second_look(
        &data,
        policy,
        &FirstEm {
            fragments: &fragments,
            group_of: &all.group_of,
            em: &first,
            resolved: &resolved_taxid,
        },
        &tests,
        &mut counts,
        &mut genomes,
    );
    log::info!(
        target: "PROFILER",
        "Level-wise engine: second look kept {} removed species on {} fragments their removed sisters tied on.",
        second_kept,
        second_taken
    );
    let rescued = own_evidence_look(
        &FirstEm {
            fragments: &fragments,
            group_of: &all.group_of,
            em: &first,
            resolved: &resolved_taxid,
        },
        &tests,
        &mut counts,
        &mut genomes,
        |taxid, list| coverage_of(&data, policy, &fragments, &tests, taxid, list),
    );
    let at_genus: HashSet<u32> = rescued
        .iter()
        .filter(|&(_, &genus)| genus)
        .map(|(&taxid, _)| taxid)
        .collect();
    log::info!(
        target: "PROFILER",
        "Level-wise engine: {} removed species kept on fragments that no kept species explains; {} of them \
         reported at their genus, as their reads tie with other species along their genome.",
        rescued.len(),
        at_genus.len()
    );
    if cfg.keep_raw {
        write_genome_diagnostics(profiler, &data, &genomes, &counts, (min_oebr, min_coebr))?;
    }
    // Best-fit genome of each kept species: its kept genome where most of its
    // fragments fit best, then the one with most fragments.
    let mut best_fit: HashMap<u32, u32> = HashMap::new();
    for (&assembly_id, genome) in &genomes {
        if !genome.kept {
            continue;
        }
        let key = |id: u32, g: &GenomeEvidence| {
            (g.assembly.fits, g.assembly.count, std::cmp::Reverse(id))
        };
        let better = best_fit
            .get(&genome.taxid)
            .is_none_or(|&current| key(assembly_id, genome) > key(current, &genomes[&current]));
        if better {
            best_fit.insert(genome.taxid, assembly_id);
        }
    }

    // Genus level, where the data carries the roll-up: each measured species'
    // strongest genome is its decision, and a removed species whose reads
    // cover its genome evenly is kept on that genome.
    let rollup = data.genus_rollup.as_deref();
    let decisions = rollup.map_or_else(HashMap::new, |rollup| {
        species_decisions(&mut genomes, rollup)
    });
    let decision_metrics: HashMap<u32, WitnessDecision> = decisions
        .iter()
        .map(|(&taxid, &(decision, _))| (taxid, decision))
        .collect();
    let mut even_species = 0;
    if let Some(rollup) = rollup {
        for taxid in rollup.even_coverage_species(
            &decision_metrics,
            |taxid| data.terminal_ranks.get(&taxid).copied(),
            data.min_chunk_breadth,
        ) {
            best_fit.insert(taxid, decisions[&taxid].1);
            even_species += 1;
        }
    }

    // Second EM over the kept species.
    let kept = share_reads(&fragments, |taxid| best_fit.contains_key(&taxid), penalty);
    let em = likelihood_em::run(
        &kept.groups,
        &HashMap::new(),
        penalty,
        cfg.max_em_iterations,
        EM_TOLERANCE,
    );
    let without_kept = kept.group_of.iter().filter(|g| g.is_none()).count();
    log::info!(
        target: "PROFILER",
        "Level-wise engine: {} species kept, {} read groups, EM {} + {} iterations; {} fragments without a kept species.",
        best_fit.len(),
        kept.groups.len(),
        first.iterations,
        em.iterations,
        without_kept
    );
    let genus_rows = rollup.map_or_else(Vec::new, |rollup| {
        let accepted: HashSet<u32> = best_fit.keys().copied().collect();
        orphan_genus_rows(
            &data,
            &fragments,
            &kept.group_of,
            policy,
            rollup,
            &decision_metrics,
            &accepted,
        )
    });
    if rollup.is_some() {
        log::info!(
            target: "PROFILER",
            "Level-wise engine: {} species kept for even coverage; {} genus rows from {:.0} fragments of removed species.",
            even_species,
            genus_rows.len(),
            // A float sum of no rows is -0.0, printed "-0".
            genus_rows.iter().fold(0.0, |total, row| total + row.fragments)
        );
    }
    // Family rows from fragments with no candidate at the species identity,
    // built with the default engine's genus fallback gates: they come from
    // genomes missing from the database, which their reads cannot place in a
    // genus.
    let family_level = rollup.map_or_else(Vec::new, |rollup| {
        let family_of = family_map(&profiler.taxonomy, &data.genus_of_species);
        // On a fractional database, reads from the unindexed parts of a kept
        // genome land on the sampled segments of its relatives, while the kept
        // species' reads already stand for its whole genome: its family gets
        // no family row (as its genus gets no roll-up row).
        let excluded: HashSet<u32> = if rollup.fractional {
            best_fit
                .keys()
                .filter_map(|taxid| family_of.get(taxid).copied())
                .collect()
        } else {
            HashSet::new()
        };
        unmatched_rows(
            profiler,
            &data,
            &unmatched,
            &UnmatchedLevel {
                policy: EvidencePolicy {
                    identity: cfg.genus_identity,
                    mapped_len: super::genus_fallback_min_mapped_len(cfg),
                    fraction: super::GENUS_FALLBACK_MIN_FRACTION,
                    collect_exact_ties: false,
                    collect_identities: false,
                    collect_best_placements: false,
                    ..policy
                },
                penalty,
                tests: &tests,
                layout: &rollup.layout,
                unit_of: &family_of,
                excluded: &excluded,
                rank: "family",
            },
        )
    });
    if rollup.is_some() {
        log::info!(
            target: "PROFILER",
            "Level-wise engine: {} fragments without a species-level candidate; {} family rows from {:.0} of them.",
            unmatched.len(),
            family_level.len(),
            family_level
                .iter()
                .fold(0.0, |total, (_, row)| total + row.equivalent_mapcount)
        );
    }
    let group_bases = &kept.group_bases;
    // Each species' share of the fragments' bases.
    let mut bases: HashMap<u32, f64> = HashMap::new();
    for (posteriors, &group_total) in em.posteriors.iter().zip(group_bases) {
        for &(taxid, posterior) in posteriors {
            *bases.entry(taxid).or_insert(0.0) += group_total * posterior;
        }
    }

    // Genomes within each kept species, for the fragments the species-level EM
    // gives to one species.
    let species_of = kept
        .group_of
        .iter()
        .map(|group| {
            group.and_then(|g| {
                em.posteriors[g]
                    .iter()
                    .find(|p| p.1 >= RESOLVE_POSTERIOR)
                    .map(|p| p.0)
            })
        })
        .collect::<Vec<_>>();
    let genome_result = genomes::resolve(
        &data,
        &fragments,
        &species_of,
        &genomes::GenomeLevel {
            policy,
            error_rate,
            penalty,
            report_strains: cfg.genome_level,
            chunk_cutoff: data.min_chunk_breadth,
            ratio_bounds: (min_oebr, min_coebr),
            max_iterations: cfg.max_em_iterations,
            tolerance: EM_TOLERANCE,
        },
    );
    let strain_taxids = if cfg.genome_level {
        strain_taxids(cfg)
    } else {
        HashMap::new()
    };
    let species_genomes = &genome_result.species;
    let baseline_identity = genome_result.baseline_identity;
    let strain_baseline = genome_result.strain_baseline_identity;
    log::info!(
        target: "PROFILER",
        "Level-wise engine: baseline identity {:.5} (median identity of the species to their best-fit genome); \
         {:.5} over strain-level matches.",
        baseline_identity,
        strain_baseline
    );
    if cfg.keep_raw && cfg.genome_level {
        write_strain_diagnostics(profiler, &data, &genome_result.diagnostics)?;
    }
    // Estimated ANI of the sample's strain to a genome: its fragments' mean
    // identity to the genome relative to the sample's baseline identity.
    let ani = |identity: f64| (1.0 - (baseline_identity - identity)).min(1.0);
    // The genome each kept species is reported on: its best-fit genome when
    // the genome level found one among its measured genomes, otherwise the
    // kept genome chosen above.
    let reported: HashMap<u32, u32> = best_fit
        .iter()
        .map(|(&taxid, &fallback)| {
            let id = species_genomes
                .get(&taxid)
                .map(|g| g.best_fit)
                .filter(|id| genomes.contains_key(id))
                .unwrap_or(fallback);
            (taxid, id)
        })
        .collect();
    let markers = data
        .marker_model
        .as_deref()
        .map_or_else(MarkerEvidence::default, |model| {
            marker_evidence(
                &data,
                policy,
                &fragments,
                (&kept, &em),
                &reported,
                (
                    model,
                    cfg.estimate_unclassified.then(|| cfg.coverage_gate_is_or()),
                ),
            )
        });
    let marker_depth = &markers.depth;

    // Species rows, each followed by the genomes reported on their own. The
    // rows of a species are exclusive: the species row keeps what no reported
    // genome holds and is left out when less than one read is left.
    let mut blocks: Vec<TaxonBlock> = Vec::new();
    // Species reported at their genus: (genus, species, row).
    let mut genus_parts: Vec<(u32, u32, ProfileRow)> = Vec::new();
    for (&taxid, &assembly_id) in &reported {
        let found = species_genomes.get(&taxid);
        let genome = &genomes[&assembly_id];
        // Full-index-equivalent fragments, and the depth of the species' bases
        // on its genome; a species on a marker genome has no genome rows and
        // takes the depth its fragments add there.
        let scale = genome.assembly.sampling_scale.max(1.0);
        let depth = |bases: f64| {
            marker_depth
                .get(&taxid)
                .copied()
                .unwrap_or(bases / genome.assembly.asm_len.max(1.0))
        };
        let coverage = genome.coverage.as_ref().expect("measured above");
        let reported = found.map_or(&[][..], |g| &g.strains[..]);
        let species_name = profiler.taxonomy.get_name(taxid).unwrap_or("unknown");
        // Every row of the level-wise engine has the identity columns, NA when
        // the species has no fragments of its own.
        let identity = Some(found.map_or((f64::NAN, f64::NAN), |g| {
            (g.mean_identity, ani(g.mean_identity))
        }));
        if let Some(&genus) = at_genus
            .contains(&taxid)
            .then(|| data.genus_of_species.get(&taxid))
            .flatten()
        {
            let reads = em.weight.get(&taxid).copied().unwrap_or(0.0);
            let species_bases = bases.get(&taxid).copied().unwrap_or(0.0);
            let genus_name = profiler.taxonomy.get_name(genus).unwrap_or("unknown");
            genus_parts.push((
                genus,
                taxid,
                profile_row(
                    genus_name.to_string(),
                    genus,
                    "genus",
                    (reads * scale, depth(species_bases)),
                    coverage,
                    data.assembly_names[assembly_id as usize].clone(),
                    Some((f64::NAN, f64::NAN)),
                ),
            ));
            continue;
        }
        let species = species_remainder(
            em.weight.get(&taxid).copied().unwrap_or(0.0),
            bases.get(&taxid).copied().unwrap_or(0.0),
            reported.iter().map(|s| (s.reads, s.bases)),
        )
        .map(|(reads, species_bases)| {
            profile_row(
                species_name.to_string(),
                taxid,
                "species",
                (reads * scale, depth(species_bases)),
                coverage,
                data.assembly_names[assembly_id as usize].clone(),
                identity,
            )
        });
        let mut strains = Vec::new();
        for strain in reported {
            let strain_ani = (1.0 - (strain_baseline - strain.mean_identity)).min(1.0);
            let accession = data.assembly_names[strain.assembly_id as usize].clone();
            let (name, strain_taxid) = match strain_taxids.get(&*accession) {
                Some(&t) => (
                    profiler
                        .taxonomy
                        .get_name(t)
                        .unwrap_or(species_name)
                        .to_string(),
                    t,
                ),
                None => (format!("{species_name} {accession}"), taxid),
            };
            strains.push(profile_row(
                name,
                strain_taxid,
                "strain",
                (
                    strain.reads
                        * genomes
                            .get(&strain.assembly_id)
                            .map_or(1.0, |g| g.assembly.sampling_scale.max(1.0)),
                    strain.bases / strain.asm_len.max(1.0),
                ),
                &strain.coverage,
                accession,
                Some((strain.mean_identity, strain_ani)),
            ));
        }
        blocks.push(TaxonBlock {
            taxid,
            taxon: species,
            genomes: strains,
        });
    }
    // One genus row per genus: the roll-up's, then the species reported at
    // their genus, added in taxid order.
    genus_parts.sort_by_key(|part| (part.0, part.1));
    let mut genus_blocks: Vec<(u32, ProfileRow)> = genus_rows
        .iter()
        .map(|row| {
            (
                row.genus,
                genus_profile_row(
                    profiler
                        .taxonomy
                        .get_name(row.genus)
                        .unwrap_or("unknown")
                        .to_string(),
                    row,
                    data.assembly_names[decisions[&row.donors[0]].1 as usize].clone(),
                ),
            )
        })
        .collect();
    for (genus, _, part) in genus_parts {
        match genus_blocks.iter_mut().find(|block| block.0 == genus) {
            Some((_, row)) => {
                row.equivalent_mapcount += part.equivalent_mapcount;
                row.depth += part.depth;
            }
            None => genus_blocks.push((genus, part)),
        }
    }
    for (taxid, row) in genus_blocks.into_iter().chain(family_level) {
        blocks.push(TaxonBlock {
            taxid,
            taxon: Some(row),
            genomes: Vec::new(),
        });
    }
    let total_depth: f64 = blocks.iter().map(TaxonBlock::depth).sum();
    let abundance = |depth: f64| {
        if total_depth > 0.0 {
            100.0 * depth / total_depth
        } else {
            0.0
        }
    };
    for block in &mut blocks {
        for row in block.taxon.iter_mut().chain(block.genomes.iter_mut()) {
            row.abundance = abundance(row.depth);
        }
    }
    // By the taxon's total depth, then taxid.
    blocks.sort_by(|a, b| b.depth().total_cmp(&a.depth()).then(a.taxid.cmp(&b.taxid)));
    let rows = blocks
        .into_iter()
        .flat_map(|block| block.taxon.into_iter().chain(block.genomes))
        .collect::<Vec<_>>();

    let mut raw_rows = Vec::new();
    for (&assembly_id, genome) in &genomes {
        let coverage = genome.coverage.as_ref().expect("measured above");
        raw_rows.push(profile_row(
            profiler
                .taxonomy
                .get_name(genome.taxid)
                .unwrap_or("unknown")
                .to_string(),
            genome.taxid,
            "species",
            (
                genome.assembly.count as f64 * genome.assembly.sampling_scale.max(1.0),
                coverage.depth,
            ),
            coverage,
            data.assembly_names[assembly_id as usize].clone(),
            None,
        ));
    }
    // One row per genome: genomes of a species with the same reads are ordered
    // by accession, so the order does not depend on hash-map iteration.
    raw_rows.sort_by(|a, b| {
        b.equivalent_mapcount
            .total_cmp(&a.equivalent_mapcount)
            .then(a.taxid.cmp(&b.taxid))
            .then_with(|| row_genome(a).cmp(&row_genome(b)))
    });
    let unclassified = profiler.estimate_unclassified(&data, &rows, &markers.genome_reads);
    profiler.write_profile_files(&rows, unclassified, raw_rows)
}

/// Every fragment with at least one candidate species, and the sample's error
/// rate.
fn collect_fragments(
    data: &MappingData,
    policy: EvidencePolicy,
) -> (Vec<Fragment>, Vec<usize>, f64) {
    let reads = (0..data.read_alignments.len()).collect::<Vec<_>>();
    let mut fragments = Vec::with_capacity(reads.len());
    let mut unmatched = Vec::new();
    let mut best_matches: Vec<(f32, f32)> = Vec::with_capacity(reads.len());
    fold_candidates(
        data,
        policy,
        &reads,
        |&read| read,
        |&read, candidates| fragment_of(read, candidates, policy, Some),
        |&read, found| match found {
            Some((fragment, best)) => {
                fragments.push(fragment);
                best_matches.push(best);
            }
            None => unmatched.push(read),
        },
    );
    (fragments, unmatched, error_rate(&mut best_matches))
}

/// Estimated errors of a candidate's best placement, counting the bases it
/// aligns less than the fragment's longest alignment (`longest`).
fn candidate_errors(candidate: &TerminalCandidate, longest: f64) -> f64 {
    candidate.best.pair_errors
        + UNALIGNED_ERROR_WEIGHT * (longest - f64::from(candidate.best.both_mapped)).max(0.0)
}

/// Bases of a fragment's longest candidate alignment.
fn longest_alignment(candidates: &[TerminalCandidate]) -> f64 {
    candidates
        .iter()
        .map(|c| f64::from(c.best.both_mapped))
        .fold(0.0, f64::max)
}

/// Read `read` as a fragment with its candidate units and their errors, with
/// the errors and aligned bases of its best candidate; `None` without one.
/// `unit` gives the unit a candidate species counts for (the species itself,
/// or its genus); a unit keeps the fewest errors of its species.
fn fragment_of(
    read: usize,
    candidates: &[TerminalCandidate],
    policy: EvidencePolicy,
    unit: impl Fn(u32) -> Option<u32>,
) -> Option<(Fragment, (f32, f32))> {
    let longest = longest_alignment(candidates);
    let mut units = candidates
        .iter()
        .filter(|c| c.profile_rank == ProfileRank::Species)
        .filter_map(|c| {
            let errors = candidate_errors(c, longest);
            Some((unit(c.taxid)?, errors, c.best.both_mapped, c.best.r1_mapped))
        })
        .collect::<Vec<_>>();
    if units.is_empty() {
        return None;
    }
    units.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
    units.dedup_by_key(|u| u.0);
    units.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let fewest = units[0].1;
    let margin = match policy.mode {
        EvidenceMode::Long => LONG_TIE_SHARE * f64::from(units[0].2),
        EvidenceMode::Paired | EvidenceMode::Single => SHORT_TIE_ERRORS,
    };
    let fragment = Fragment {
        read,
        candidates: units
            .iter()
            .map(|c| (c.0, if c.1 - fewest < margin { fewest } else { c.1 }))
            .collect(),
        bases: f64::from(units[0].3),
    };
    Some((fragment, (fewest as f32, units[0].2 as f32)))
}

/// Group the fragments for the EM over the species that `allowed` admits,
/// dropping candidates whose likelihood is below [`MIN_RELATIVE_LIKELIHOOD`]
/// of the fragment's best allowed candidate.
fn share_reads(fragments: &[Fragment], allowed: impl Fn(u32) -> bool, penalty: f64) -> SharedReads {
    let max_excess_errors = -MIN_RELATIVE_LIKELIHOOD.ln() / penalty;
    let mut shared = SharedReads {
        groups: ReadGroups::default(),
        group_of: Vec::with_capacity(fragments.len()),
        group_bases: Vec::new(),
    };
    for fragment in fragments {
        let mut candidates = fragment
            .candidates
            .iter()
            .filter(|c| allowed(c.0))
            .copied()
            .collect::<Vec<_>>();
        let fewest = candidates.iter().map(|c| c.1).fold(f64::INFINITY, f64::min);
        candidates.retain(|c| c.1 - fewest <= max_excess_errors);
        let slot = shared.groups.add(candidates, 1.0);
        if let Some(slot) = slot {
            if slot == shared.group_bases.len() {
                shared.group_bases.push(0.0);
            }
            shared.group_bases[slot] += fragment.bases;
        }
        shared.group_of.push(slot);
    }
    shared
}

/// Errors per aligned base of the fragments' best candidates, over the share
/// [`ERROR_RATE_SHARE`] that match best: the sample's sequencing error rate,
/// little affected by reads from relatives of the database genomes.
fn error_rate(best_matches: &mut [(f32, f32)]) -> f64 {
    best_matches.sort_by(|a, b| (a.0 / a.1.max(1.0)).total_cmp(&(b.0 / b.1.max(1.0))));
    let kept = &best_matches[..(best_matches.len() as f64 * ERROR_RATE_SHARE).ceil() as usize];
    let errors: f64 = kept.iter().map(|m| f64::from(m.0)).sum();
    let bases: f64 = kept.iter().map(|m| f64::from(m.1)).sum();
    if bases > 0.0 {
        errors / bases
    } else {
        0.01
    }
}

fn add_placement(
    genomes: &mut HashMap<u32, GenomeEvidence>,
    taxid: u32,
    placement: &CandidatePlacement,
    fits: bool,
) {
    let genome = genomes
        .entry(placement.assembly_id)
        .or_insert_with(|| GenomeEvidence {
            taxid,
            assembly: PresenceAssembly {
                asm_len: placement.asm_len,
                sampling_scale: placement.sampling_scale,
                ..PresenceAssembly::default()
            },
            coverage: None,
            kept: false,
        });
    let assembly = &mut genome.assembly;
    assembly.count += 1;
    assembly.fits += u64::from(fits);
    assembly.other_gene_hits += u64::from(placement.other_genes);
    assembly.reference_span_bases += u64::from(placement.span.1.saturating_sub(placement.span.0));
    assembly.intervals.push(PresenceInterval {
        reference_id: placement.reference_id,
        span: placement.span,
        mate_span: placement.mate_span,
        fixed_chunk: placement.fixed_chunk,
        flex_chunk: placement.flex_chunk,
    });
}

/// What the fragments of species reported on marker genomes add: each such
/// species' depth, and for the unclassified estimate the marker genome reads
/// of the fragments resolved to it.
#[derive(Default)]
struct MarkerEvidence {
    depth: HashMap<u32, f64>,
    genome_reads: HashMap<u32, GenomeReads>,
}

/// [`MarkerEvidence`] of the kept species reported on a marker genome
/// (`reported`). A species' depth is the sum over its fragments of their share
/// in the final EM (`shared`, `em`) times the depth one
/// fragment adds on that genome, [`fragment_depth`]: its read length times the
/// genes it hits over the placements that the genome's genes offer a read of
/// its length. With `unclassified` (`--estimate-unclassified`; `Some(true)` for
/// long reads), each fragment resolved to such a species
/// ([`RESOLVE_POSTERIOR`]) also records its placements on marker genomes, as in
/// the default engine. Only the fragments whose read group lists such a
/// species are read again.
fn marker_evidence(
    data: &MappingData,
    policy: EvidencePolicy,
    fragments: &[Fragment],
    (shared, em): (&SharedReads, &likelihood_em::EmResult),
    reported: &HashMap<u32, u32>,
    (model, unclassified): (&MarkerModel, Option<bool>),
) -> MarkerEvidence {
    let on_markers = |taxid: u32| {
        reported
            .get(&taxid)
            .is_some_and(|&id| model.is_marker_genome(id))
    };
    let items = fragments
        .iter()
        .zip(&shared.group_of)
        .filter_map(|(fragment, group)| {
            let group = (*group)?;
            em.posteriors[group]
                .iter()
                .any(|p| on_markers(p.0))
                .then_some((fragment.read, group))
        })
        .collect::<Vec<_>>();
    let mut evidence = MarkerEvidence::default();
    fold_candidates(
        data,
        policy,
        &items,
        |item| item.0,
        |&(_, group), candidates| {
            em.posteriors[group]
                .iter()
                .filter(|p| on_markers(p.0))
                .filter_map(|&(taxid, share)| {
                    let candidate = candidates.iter().find(|c| c.taxid == taxid)?;
                    let placement = candidate
                        .placements
                        .iter()
                        .find(|p| p.assembly_id == reported[&taxid])
                        .copied()
                        .unwrap_or(candidate.best);
                    let mapped = match policy.mode {
                        EvidenceMode::Paired | EvidenceMode::Single => placement.r1_mapped,
                        EvidenceMode::Long => placement.both_mapped,
                    };
                    let resolved =
                        (unclassified.is_some() && share >= RESOLVE_POSTERIOR).then(|| {
                            (
                                candidate.placements.clone(),
                                candidate.other_gene_hits.clone(),
                            )
                        });
                    Some((
                        taxid,
                        share * fragment_depth(placement, mapped, Some(model)),
                        resolved,
                    ))
                })
                .collect::<Vec<_>>()
        },
        |_, found| {
            for (taxid, added, resolved) in found {
                *evidence.depth.entry(taxid).or_insert(0.0) += added;
                for (placement, hits) in resolved
                    .iter()
                    .flat_map(|(p, h)| p.iter().map(move |p| (p, h)))
                {
                    add_marker_genome_reads(
                        model,
                        unclassified == Some(true),
                        placement,
                        hits,
                        &mut evidence.genome_reads,
                    );
                }
            }
        },
    );
    evidence
}

/// Build the candidates of each item's read on the profiler's threads, take
/// what a pass needs with `extract`, and hand the results to `fold` in item
/// order, [`PARALLEL_BLOCK`] items at a time.
pub(super) fn fold_candidates<I: Sync, T: Send>(
    data: &MappingData,
    policy: EvidencePolicy,
    items: &[I],
    read_of: impl Fn(&I) -> usize + Sync,
    extract: impl Fn(&I, &[TerminalCandidate]) -> T + Sync,
    mut fold: impl FnMut(&I, T),
) {
    for block in items.chunks(PARALLEL_BLOCK) {
        let results = block
            .par_iter()
            .map_init(CandidateScratch::default, |scratch, item| {
                let candidates =
                    scratch.build(&data.read_alignments[read_of(item)].alignments, policy);
                extract(item, candidates)
            })
            .collect::<Vec<_>>();
        for (item, result) in block.iter().zip(results) {
            fold(item, result);
        }
    }
}

/// A fragment's placements for one species, each with whether it is on the
/// species' best genome, and its hits on further marker genes of those genomes.
#[derive(Default)]
struct SpeciesHits {
    placements: Vec<(CandidatePlacement, bool)>,
    gene_hits: Vec<GeneHit>,
}

/// The placements and gene hits of species `taxid` among a fragment's
/// `candidates`.
fn species_placements(candidates: &[TerminalCandidate], taxid: u32) -> SpeciesHits {
    candidates
        .iter()
        .find(|c| c.taxid == taxid)
        .map_or_else(SpeciesHits::default, |candidate| SpeciesHits {
            placements: candidate
                .placements
                .iter()
                .map(|p| (*p, p.assembly_id == candidate.best.assembly_id))
                .collect(),
            gene_hits: candidate.other_gene_hits.clone(),
        })
}

/// Add a fragment's placements for species `taxid` to their genomes, and the
/// bases and intervals of its other marker genes, which a gene-by-gene
/// measure counts as covered (genomes of other kinds have none).
fn add_hits(genomes: &mut HashMap<u32, GenomeEvidence>, taxid: u32, hits: SpeciesHits) {
    for (placement, fits) in &hits.placements {
        add_placement(genomes, taxid, placement, *fits);
    }
    for hit in hits.gene_hits {
        if let Some(genome) = genomes.get_mut(&hit.assembly_id) {
            let assembly = &mut genome.assembly;
            assembly.reference_span_bases +=
                u64::from(hit.interval.span.1.saturating_sub(hit.interval.span.0));
            assembly.intervals.push(hit.interval);
        }
    }
}

/// Fragments listing each species among their candidates, those resolved to
/// any species, those resolved to it, those it took on the second look, and
/// those it won with no kept species among their candidates.
#[derive(Default)]
struct FragmentCounts {
    candidates: HashMap<u32, u64>,
    resolved_candidates: HashMap<u32, u64>,
    won: HashMap<u32, u64>,
    second_look: HashMap<u32, u64>,
    own: HashMap<u32, u64>,
}

/// The species tests: enough coverage, evenly spread, and a large enough share
/// of the fragments listing the species won. A genome represented by marker
/// genes is measured by the `markers` model, its genes standing for the chunks
/// of a whole genome, and needs at least `min_marker_genes` genes hit.
struct SpeciesTests {
    multiple_tests: usize,
    chunk_cutoff: f64,
    ratio_bounds: (f64, f64),
    markers: Option<Arc<MarkerModel>>,
    min_marker_genes: usize,
}

impl SpeciesTests {
    /// Coverage of genome `assembly_id`: gene by gene for a marker genome,
    /// otherwise over its sequence.
    fn measure(&self, assembly_id: u32, assembly: &mut PresenceAssembly) -> EvidenceCoverage {
        self.markers
            .as_deref()
            .and_then(|model| model.coverage(assembly_id, assembly, self.multiple_tests))
            .unwrap_or_else(|| calculate_coverage(assembly, self.multiple_tests))
    }

    /// Measure genome `assembly_id` and decide whether it keeps its species,
    /// which won `won` of the `candidates` fragments listing it.
    fn judge(&self, assembly_id: u32, genome: &mut GenomeEvidence, won: u64, candidates: u64) {
        let coverage = self.measure(assembly_id, &mut genome.assembly);
        genome.kept =
            self.covers(&coverage) && won as f64 >= MIN_WON_SHARE * candidates.max(1) as f64;
        genome.coverage = Some(coverage);
    }

    /// Whether a genome's coverage is enough and evenly spread.
    fn covers(&self, coverage: &EvidenceCoverage) -> bool {
        coverage.fixed_chunk > self.chunk_cutoff
            && coverage
                .marker_genes_hit
                .is_none_or(|hit| hit >= self.min_marker_genes)
            && within(coverage.oebr, self.ratio_bounds.0)
            && within(coverage.coebr, self.ratio_bounds.1)
    }
}

/// The first EM over every candidate species and the species each read group
/// resolves to.
struct FirstEm<'a> {
    fragments: &'a [Fragment],
    group_of: &'a [Option<usize>],
    em: &'a likelihood_em::EmResult,
    resolved: &'a [Option<u32>],
}

/// Second look at the removed species, largest first-EM weight first. Sister
/// species that tie on most fragments split them in the first EM by their
/// distinguishing fragments, so that neither reaches the resolving posterior
/// and both are removed, the present one with them. Judged again, a removed
/// species takes the read groups on which no other removed species is more
/// likely and that it holds with at least the resolving posterior against the
/// kept species, and is tested again on a copy of its genomes with them added.
/// Only if it passes does the copy replace its evidence; it is then kept and
/// competes with the species judged after it. A species that still fails is
/// left as the first judgment measured it, so the second look only adds
/// species. Every group it won before is still won, as its posterior against
/// fewer rivals only rises. Returns the number of species kept and of the
/// fragments they took.
fn second_look(
    data: &MappingData,
    policy: EvidencePolicy,
    first: &FirstEm,
    tests: &SpeciesTests,
    counts: &mut FragmentCounts,
    genomes: &mut HashMap<u32, GenomeEvidence>,
) -> (usize, u64) {
    let mut kept: HashSet<u32> = genomes
        .values()
        .filter(|g| g.kept)
        .map(|g| g.taxid)
        .collect();
    let mut order = first
        .em
        .weight
        .iter()
        .filter(|&(taxid, &weight)| weight > 0.0 && !kept.contains(taxid))
        .map(|(&taxid, &weight)| (taxid, weight))
        .collect::<Vec<_>>();
    order.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    let removed: HashSet<u32> = order.iter().map(|o| o.0).collect();
    let mut open_groups: HashMap<u32, Vec<usize>> = HashMap::new();
    for (group, posteriors) in first.em.posteriors.iter().enumerate() {
        for &(taxid, _) in posteriors {
            if removed.contains(&taxid) && first.resolved[group] != Some(taxid) {
                open_groups.entry(taxid).or_default().push(group);
            }
        }
    }
    let mut genomes_of: HashMap<u32, Vec<u32>> = HashMap::new();
    for (&id, genome) in genomes.iter() {
        if removed.contains(&genome.taxid) {
            genomes_of.entry(genome.taxid).or_default().push(id);
        }
    }
    let members = GroupMembers::new(first.group_of, first.em.posteriors.len());
    let (mut kept_now, mut taken_total) = (0, 0);
    for &(taxid, _) in &order {
        let Some(groups) = open_groups.get(&taxid) else {
            continue;
        };
        let taken_fragments = groups
            .iter()
            .filter(|&&group| takes_group(&first.em.posteriors[group], taxid, &kept))
            .flat_map(|&group| members.of(group))
            .copied()
            .collect::<Vec<_>>();
        if taken_fragments.is_empty() {
            continue;
        }
        let mut trial: HashMap<u32, GenomeEvidence> = genomes_of
            .get(&taxid)
            .into_iter()
            .flatten()
            .map(|id| (*id, genomes[id].clone()))
            .collect();
        let reads = taken_fragments
            .iter()
            .map(|&fragment| first.fragments[fragment].read)
            .collect::<Vec<_>>();
        fold_candidates(
            data,
            policy,
            &reads,
            |&read| read,
            |_, candidates| species_placements(candidates, taxid),
            |_, hits| add_hits(&mut trial, taxid, hits),
        );
        let taken = taken_fragments.len() as u64;
        counts.second_look.insert(taxid, taken);
        let won = counts.won.get(&taxid).copied().unwrap_or(0) + taken;
        let candidates = counts.candidates.get(&taxid).copied().unwrap_or(0);
        let mut passes = false;
        for (&id, genome) in &mut trial {
            tests.judge(id, genome, won, candidates);
            passes |= genome.kept;
        }
        if passes {
            counts.won.insert(taxid, won);
            genomes.extend(trial);
            kept.insert(taxid);
            kept_now += 1;
            taken_total += taken;
        }
    }
    (kept_now, taken_total)
}

/// Last look at the removed species that cover one of their genomes well
/// enough but win too small a share of their candidate fragments, largest
/// first-EM weight first; `coverage` measures a species' genomes on a list of
/// fragments alone.
///
/// - A species whose fragments shared with a kept species (both among the
///   candidates) alone span one of its genomes, above the chunk-breadth
///   cutoff, cannot be told apart from that relative and stays removed, its
///   reads with the relative.
/// - Otherwise it is kept when its own fragments, those it won in the first
///   judgment with no kept species among their candidates, alone pass the
///   coverage tests on one of its genomes; it then counts as kept for the
///   species judged after it. A species that collects reads of a present
///   relative's strain owns almost none of them; one that collects reads of a
///   genome missing from the database, or of an element shared across
///   genera, owns only the parts where it matches best, which leaves its own
///   coverage patchy.
/// - A kept species is reported at its genus when the fragments where it
///   shares the best match with another species alone span one of its
///   genomes: a near-identical relative ties along the whole genome, an
///   element shared with other genera in one region only.
///
/// Returns the species kept, each with whether it is reported at its genus.
fn own_evidence_look(
    first: &FirstEm,
    tests: &SpeciesTests,
    counts: &mut FragmentCounts,
    genomes: &mut HashMap<u32, GenomeEvidence>,
    mut coverage: impl FnMut(u32, &[usize]) -> Vec<(u32, EvidenceCoverage)>,
) -> HashMap<u32, bool> {
    let mut kept: HashSet<u32> = genomes
        .values()
        .filter(|g| g.kept)
        .map(|g| g.taxid)
        .collect();
    let judged: HashSet<u32> = genomes
        .values()
        .filter(|g| !kept.contains(&g.taxid))
        .filter(|g| g.coverage.as_ref().is_some_and(|c| tests.covers(c)))
        .map(|g| g.taxid)
        .collect();
    // The fragments each such species won, and those listing it.
    let mut won: HashMap<u32, Vec<usize>> = HashMap::new();
    let mut listing: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, (fragment, group)) in first.fragments.iter().zip(first.group_of).enumerate() {
        if let Some(taxid) = group.and_then(|g| first.resolved[g]) {
            if judged.contains(&taxid) {
                won.entry(taxid).or_default().push(index);
            }
        }
        for &(taxid, _) in &fragment.candidates {
            if judged.contains(&taxid) {
                listing.entry(taxid).or_default().push(index);
            }
        }
    }
    let mut order = won
        .keys()
        .map(|&taxid| (taxid, first.em.weight.get(&taxid).copied().unwrap_or(0.0)))
        .collect::<Vec<_>>();
    order.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    let spans_genome = |coverages: &[(u32, EvidenceCoverage)]| {
        coverages
            .iter()
            .any(|(_, c)| c.fixed_chunk > tests.chunk_cutoff)
    };
    let mut rescued = HashMap::new();
    for (taxid, _) in order {
        let with_kept = |f: &usize| {
            first.fragments[*f]
                .candidates
                .iter()
                .any(|c| c.0 != taxid && kept.contains(&c.0))
        };
        let listed = &listing[&taxid];
        let shared = listed.iter().copied().filter(with_kept).collect::<Vec<_>>();
        if !shared.is_empty() && spans_genome(&coverage(taxid, &shared)) {
            continue;
        }
        let own = won[&taxid]
            .iter()
            .copied()
            .filter(|f| !with_kept(f))
            .collect::<Vec<_>>();
        counts.own.insert(taxid, own.len() as u64);
        if own.is_empty() {
            continue;
        }
        let own_genomes: HashSet<u32> = coverage(taxid, &own)
            .into_iter()
            .filter(|(_, c)| tests.covers(c))
            .map(|(id, _)| id)
            .collect();
        if own_genomes.is_empty() {
            continue;
        }
        for (id, genome) in genomes.iter_mut() {
            if own_genomes.contains(id) {
                genome.kept = true;
            }
        }
        let tied = listed
            .iter()
            .copied()
            .filter(|&f| {
                let candidates = &first.fragments[f].candidates;
                let best = candidates
                    .iter()
                    .take_while(|c| c.1 == candidates[0].1)
                    .count();
                best > 1 && candidates[..best].iter().any(|c| c.0 == taxid)
            })
            .collect::<Vec<_>>();
        let at_genus = !tied.is_empty() && spans_genome(&coverage(taxid, &tied));
        kept.insert(taxid);
        rescued.insert(taxid, at_genus);
    }
    rescued
}

/// The coverage of each genome of species `taxid` from the fragments `list`
/// alone.
fn coverage_of(
    data: &MappingData,
    policy: EvidencePolicy,
    fragments: &[Fragment],
    tests: &SpeciesTests,
    taxid: u32,
    list: &[usize],
) -> Vec<(u32, EvidenceCoverage)> {
    let mut genomes: HashMap<u32, GenomeEvidence> = HashMap::new();
    fold_candidates(
        data,
        policy,
        list,
        |&f| fragments[f].read,
        |_, candidates| species_placements(candidates, taxid),
        |_, hits| add_hits(&mut genomes, taxid, hits),
    );
    genomes
        .into_iter()
        .map(|(id, mut genome)| (id, tests.measure(id, &mut genome.assembly)))
        .collect()
}

/// Whether removed species `taxid` takes a read group with first-EM
/// `posteriors` on the second look: no other removed species is more likely
/// there, and it holds at least the resolving posterior against the `kept`
/// species.
fn takes_group(posteriors: &[(u32, f64)], taxid: u32, kept: &HashSet<u32>) -> bool {
    let own = posteriors
        .iter()
        .find(|p| p.0 == taxid)
        .map_or(0.0, |p| p.1);
    let mut with_kept = own;
    for &(other, posterior) in posteriors {
        if other == taxid {
            continue;
        }
        if kept.contains(&other) {
            with_kept += posterior;
        } else if posterior > own {
            return false;
        }
    }
    own > 0.0 && own >= RESOLVE_POSTERIOR * with_kept
}

/// The fragments of each read group.
struct GroupMembers {
    start: Vec<usize>,
    fragments: Vec<usize>,
}

impl GroupMembers {
    fn new(group_of: &[Option<usize>], groups: usize) -> Self {
        let mut start = vec![0; groups + 1];
        for &group in group_of.iter().flatten() {
            start[group + 1] += 1;
        }
        for group in 0..groups {
            start[group + 1] += start[group];
        }
        let mut next = start.clone();
        let mut fragments = vec![0; start[groups]];
        for (fragment, group) in group_of.iter().enumerate() {
            if let &Some(group) = group {
                fragments[next[group]] = fragment;
                next[group] += 1;
            }
        }
        Self { start, fragments }
    }

    fn of(&self, group: usize) -> &[usize] {
        &self.fragments[self.start[group]..self.start[group + 1]]
    }
}

fn within(ratio: f64, lower: f64) -> bool {
    ratio.is_finite() && ratio >= lower && ratio <= STANDARD_RATIO_UPPER
}

/// The reads and bases of a species that its reported genomes, given as
/// `(reads, bases)`, do not hold; `None` when genomes are reported and less
/// than one read is left, so the species gets no row of its own.
fn species_remainder(
    reads: f64,
    bases: f64,
    reported: impl Iterator<Item = (f64, f64)>,
) -> Option<(f64, f64)> {
    let (mut held_reads, mut held_bases, mut any) = (0.0, 0.0, false);
    for (r, b) in reported {
        held_reads += r;
        held_bases += b;
        any = true;
    }
    let reads = (reads - held_reads).max(0.0);
    (!any || reads >= 1.0).then(|| (reads, (bases - held_bases).max(0.0)))
}

/// The rows of one taxon: a species row (unless every read is in a genome
/// row) with its genome rows, or a genus row.
struct TaxonBlock {
    taxid: u32,
    taxon: Option<ProfileRow>,
    genomes: Vec<ProfileRow>,
}

impl TaxonBlock {
    fn depth(&self) -> f64 {
        self.taxon.as_ref().map_or(0.0, |s| s.depth)
            + self.genomes.iter().map(|s| s.depth).sum::<f64>()
    }
}

/// Each measured species' decision for the genus level: its strongest genome
/// in the default engine's witness order (lowest assembly id on a tie), with
/// the window spread of its reads, and that genome.
fn species_decisions(
    genomes: &mut HashMap<u32, GenomeEvidence>,
    rollup: &GenusRollup,
) -> HashMap<u32, (WitnessDecision, u32)> {
    let mut decisions: HashMap<u32, (WitnessDecision, u32)> = HashMap::new();
    for (&assembly_id, genome) in genomes.iter_mut() {
        let Some(coverage) = genome.coverage.as_mut() else {
            continue;
        };
        coverage.window_spread = rollup
            .layout
            .window_spread(assembly_id, &genome.assembly.intervals);
        let decision = WitnessDecision {
            metrics: *coverage,
            standard: genome.kept,
            rescues: RescueKinds::default(),
            required_premap_guard: PremapChunkGuard::default(),
        };
        let better = decisions
            .get(&genome.taxid)
            .is_none_or(|&(current, current_id)| {
                witness_is_better(decision, current)
                    || (!witness_is_better(current, decision) && assembly_id < current_id)
            });
        if better {
            decisions.insert(genome.taxid, (decision, assembly_id));
        }
    }
    decisions
}

/// Where fragments without a species-level candidate go: the build `policy`
/// (the genus identity, mapped length and aligned fraction), the unit each
/// candidate species counts for (`unit_of`), the units that get no row
/// (`excluded`) and the `rank` of its rows, with the EM penalty, the species
/// tests and the genome layout.
struct UnmatchedLevel<'a> {
    policy: EvidencePolicy,
    penalty: f64,
    tests: &'a SpeciesTests,
    layout: &'a GenomeLayout,
    unit_of: &'a HashMap<u32, u32>,
    excluded: &'a HashSet<u32>,
    rank: &'static str,
}

/// Rows from the fragments `reads` that have no candidate at the species
/// identity. They are built again with the level's policy and their candidate
/// species counted for their units; an EM shares them, and a fragment resolves
/// to a unit at [`RESOLVE_POSTERIOR`]. A unit gets a row when its resolved
/// fragments number at least [`MIN_GENUS_FRAGMENTS`], reach
/// [`MIN_WINDOW_SPREAD`] of the windows of the genome that holds most of them
/// (the reads of a genome missing from the database spread over the genome of
/// a relative; reads of conserved genes from elsewhere stay in a few places),
/// and are at least [`MIN_WON_SHARE`] of the fragments listing the unit.
/// Returns the rows in unit order.
fn unmatched_rows(
    profiler: &CommunityProfiler,
    data: &MappingData,
    reads: &[usize],
    level: &UnmatchedLevel,
) -> Vec<(u32, ProfileRow)> {
    let unit_of = |taxid: u32| level.unit_of.get(&taxid).copied();
    let policy = level.policy;
    let mut fragments = Vec::new();
    fold_candidates(
        data,
        policy,
        reads,
        |&read| read,
        |&read, candidates| fragment_of(read, candidates, policy, unit_of),
        |_, found| fragments.extend(found.map(|(fragment, _)| fragment)),
    );
    let shared = share_reads(&fragments, |_| true, level.penalty);
    let em = likelihood_em::run(
        &shared.groups,
        &HashMap::new(),
        level.penalty,
        profiler.cfg.max_em_iterations,
        EM_TOLERANCE,
    );
    // Fragments listing each unit, those resolved to it, and its share of the
    // fragments' bases.
    let mut listing: HashMap<u32, f64> = HashMap::new();
    let mut won: HashMap<u32, f64> = HashMap::new();
    let mut resolved: Vec<(usize, u32)> = Vec::new();
    let mut bases: HashMap<u32, f64> = HashMap::new();
    for (fragment, group) in fragments.iter().zip(&shared.group_of) {
        let Some(group) = *group else {
            continue;
        };
        let posteriors = &em.posteriors[group];
        for &(unit, _) in posteriors {
            *listing.entry(unit).or_default() += 1.0;
        }
        if let Some(&(unit, _)) = posteriors.iter().find(|p| p.1 >= RESOLVE_POSTERIOR) {
            resolved.push((fragment.read, unit));
            *won.entry(unit).or_default() += 1.0;
        }
    }
    for (posteriors, &group_total) in em.posteriors.iter().zip(&shared.group_bases) {
        for &(unit, posterior) in posteriors {
            *bases.entry(unit).or_insert(0.0) += group_total * posterior;
        }
    }
    // Each resolved fragment on the best genome of its unit's best species.
    let mut genomes: HashMap<u32, GenomeEvidence> = HashMap::new();
    fold_candidates(
        data,
        policy,
        &resolved,
        |r| r.0,
        |&(_, unit), candidates| {
            let longest = longest_alignment(candidates);
            candidates
                .iter()
                .filter(|c| {
                    c.profile_rank == ProfileRank::Species && unit_of(c.taxid) == Some(unit)
                })
                .min_by(|a, b| {
                    candidate_errors(a, longest)
                        .total_cmp(&candidate_errors(b, longest))
                        .then(a.taxid.cmp(&b.taxid))
                })
                .map(|c| c.best)
        },
        |&(_, unit), placement| {
            if let Some(placement) = placement {
                add_placement(&mut genomes, unit, &placement, true);
            }
        },
    );
    // Per unit, the genome that holds most of its resolved fragments.
    let mut best: HashMap<u32, u32> = HashMap::new();
    for (&id, genome) in &genomes {
        let key = |id: u32| (genomes[&id].assembly.count, std::cmp::Reverse(id));
        if best
            .get(&genome.taxid)
            .is_none_or(|&current| key(id) > key(current))
        {
            best.insert(genome.taxid, id);
        }
    }
    let mut rows = Vec::new();
    for (unit, id) in best {
        if level.excluded.contains(&unit) {
            continue;
        }
        let genome = genomes.get_mut(&id).expect("measured above");
        let won = won.get(&unit).copied().unwrap_or(0.0);
        let spread = level.layout.window_spread(id, &genome.assembly.intervals);
        if won < MIN_GENUS_FRAGMENTS
            || won < MIN_WON_SHARE * listing.get(&unit).copied().unwrap_or(0.0)
            || spread.is_none_or(|s| s < MIN_WINDOW_SPREAD)
        {
            continue;
        }
        let coverage = level.tests.measure(id, &mut genome.assembly);
        let depth = bases.get(&unit).copied().unwrap_or(0.0) / genome.assembly.asm_len.max(1.0);
        rows.push((
            unit,
            profile_row(
                profiler
                    .taxonomy
                    .get_name(unit)
                    .unwrap_or("unknown")
                    .to_string(),
                unit,
                level.rank,
                (
                    em.weight.get(&unit).copied().unwrap_or(0.0)
                        * genome.assembly.sampling_scale.max(1.0),
                    depth,
                ),
                &coverage,
                data.assembly_names[id as usize].clone(),
                Some((f64::NAN, f64::NAN)),
            ),
        ));
    }
    rows.sort_by_key(|row| row.0);
    rows
}

/// The family of every species in `species_to_genus` whose genus has one.
fn family_map(taxonomy: &Taxonomy, species_to_genus: &HashMap<u32, u32>) -> HashMap<u32, u32> {
    let mut of_genus: HashMap<u32, Option<u32>> = HashMap::new();
    let mut families = HashMap::new();
    for (&species, &genus) in species_to_genus {
        let family = *of_genus.entry(genus).or_insert_with(|| {
            taxonomy
                .get_parents(genus)
                .into_iter()
                .find(|(_, _, rank)| rank == "family")
                .map(|(taxid, _, _)| taxid)
        });
        if let Some(family) = family {
            families.insert(species, family);
        }
    }
    families
}

/// Genus rows from the fragments with no kept candidate species, offered to the
/// roll-up with all their candidate species and their best one's placement.
fn orphan_genus_rows(
    data: &MappingData,
    fragments: &[Fragment],
    group_of: &[Option<usize>],
    policy: EvidencePolicy,
    rollup: &GenusRollup,
    decisions: &HashMap<u32, WitnessDecision>,
    accepted: &HashSet<u32>,
) -> Vec<GenusRow> {
    let mut orphans = OrphanedFragments::new(rollup, false);
    let mut scratch = CandidateScratch::default();
    for (fragment, group) in fragments.iter().zip(group_of) {
        if group.is_some() {
            continue;
        }
        let bundle = &data.read_alignments[fragment.read];
        // Candidates are ordered by errors; the first is the best.
        let taxids = fragment.candidates.iter().map(|c| c.0).collect::<Vec<_>>();
        let candidates = scratch.build(&bundle.alignments, policy);
        let Some(best) = candidates
            .iter()
            .find(|c| c.taxid == taxids[0])
            .map(|c| c.best)
        else {
            continue;
        };
        orphans.offer(
            &taxids,
            0,
            f64::from(best.r1_mapped) / best.asm_len.max(1.0),
            f64::from(best.r1_mapped),
            best.sampling_scale,
            &bundle.read_name,
        );
    }
    orphans.into_genus_rows(decisions, accepted)
}

/// A genus row of the roll-up: its fragments and depth, and the coverage of
/// the main donor species on `genome`, its strongest genome; identity NA.
fn genus_profile_row(name: String, row: &GenusRow, genome: std::sync::Arc<str>) -> ProfileRow {
    let coverage = &row.witness;
    ProfileRow {
        taxname: name,
        taxid: row.genus,
        rank: "genus".to_string(),
        equivalent_mapcount: row.fragments,
        depth: row.depth,
        abundance: 0.0,
        breadth: coverage.breadth,
        expected_breadth: coverage.expected_breadth,
        breadth_pvalue: coverage.breadth_pvalue,
        fixed_chunk: coverage.fixed_chunk,
        flex_chunk: coverage.flex_chunk,
        expected_flex_chunk: coverage.expected_flex_chunk,
        chunk_pvalue: coverage.chunk_pvalue,
        best_fit: Some(BestFitReport {
            genome,
            evidence_depth: coverage.depth,
            identity: Some((f64::NAN, f64::NAN)),
        }),
        annotations: None,
    }
}

/// The genome a row is measured on, if any.
fn row_genome(row: &ProfileRow) -> Option<&str> {
    row.best_fit.as_ref().map(|g| &*g.genome)
}

/// One profile row; `amount` is `(reads, depth)`.
fn profile_row(
    name: String,
    taxid: u32,
    rank: &str,
    (reads, depth): (f64, f64),
    coverage: &EvidenceCoverage,
    genome: std::sync::Arc<str>,
    identity: Option<(f64, f64)>,
) -> ProfileRow {
    ProfileRow {
        taxname: name,
        taxid,
        rank: rank.to_string(),
        equivalent_mapcount: reads,
        depth,
        abundance: 0.0,
        breadth: coverage.breadth,
        expected_breadth: coverage.expected_breadth,
        breadth_pvalue: coverage.breadth_pvalue,
        fixed_chunk: coverage.fixed_chunk,
        flex_chunk: coverage.flex_chunk,
        expected_flex_chunk: coverage.expected_flex_chunk,
        chunk_pvalue: coverage.chunk_pvalue,
        best_fit: Some(BestFitReport {
            genome,
            evidence_depth: coverage.depth,
            identity,
        }),
        annotations: None,
    }
}

/// Strain taxid of every genome that declares one.
fn strain_taxids(cfg: &crate::profiler::ProfilerConfig) -> HashMap<std::sync::Arc<str>, u32> {
    let mut out = HashMap::new();
    let Some(database) = cfg.db_metadata.as_deref() else {
        return out;
    };
    let mut previous: Option<&str> = None;
    for reference in &database.references {
        if previous == Some(&*reference.genome) {
            continue;
        }
        previous = Some(&*reference.genome);
        if let Some(strain) = parse_reference_taxonomy(&reference.taxonomy_field)
            .ok()
            .and_then(|t| t.primary_taxid)
        {
            out.insert(reference.genome.clone(), strain);
        }
    }
    out
}

/// Development diagnostics: one row per genome with resolved fragments, with
/// its coverage, the species' candidate and won fragments, and each test.
fn write_genome_diagnostics(
    profiler: &CommunityProfiler,
    data: &MappingData,
    genomes: &HashMap<u32, GenomeEvidence>,
    counts: &FragmentCounts,
    (min_oebr, min_coebr): (f64, f64),
) -> Result<()> {
    use std::io::Write;
    let path = format!("{}.levels_genomes.tsv", profiler.cfg.outprefix);
    let mut out = std::io::BufWriter::new(std::fs::File::create(&path)?);
    writeln!(
        out,
        "taxid\tname\tgenome\tfragments\tfits\tspecies_candidate_fragments\tspecies_resolved_candidate_fragments\tspecies_won_fragments\toebr\tcoebr\tfixed_chunk\tchunk_cutoff\tbreadth_pvalue\tchunk_pvalue\tamount_ok\tevenness_ok\tpvalue_ok\tsecond_look_fragments\town_won_fragments"
    )?;
    let mut ids = genomes.keys().copied().collect::<Vec<_>>();
    ids.sort_unstable();
    for id in ids {
        let genome = &genomes[&id];
        let coverage = genome.coverage.as_ref().expect("measured before");
        let p = |v: Option<f64>| v.map_or("NA".to_string(), |v| format!("{v:.3e}"));
        let count = |map: &HashMap<u32, u64>| map.get(&genome.taxid).copied().unwrap_or(0);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            genome.taxid,
            profiler
                .taxonomy
                .get_name(genome.taxid)
                .unwrap_or("unknown"),
            data.assembly_names[id as usize],
            genome.assembly.count,
            genome.assembly.fits,
            count(&counts.candidates),
            count(&counts.resolved_candidates),
            count(&counts.won),
            coverage.oebr,
            coverage.coebr,
            coverage.fixed_chunk,
            data.min_chunk_breadth,
            p(coverage.breadth_pvalue),
            p(coverage.chunk_pvalue),
            coverage.fixed_chunk > data.min_chunk_breadth,
            within(coverage.oebr, min_oebr) && within(coverage.coebr, min_coebr),
            pvalue_passes(
                coverage.breadth_pvalue,
                coverage.chunk_pvalue,
                profiler.cfg.strict
            ),
            count(&counts.second_look),
            count(&counts.own),
        )?;
    }
    Ok(())
}

/// Development diagnostics of the genome level: one row per genome that
/// fragments resolve to, and per species its best-fit genome.
fn write_strain_diagnostics(
    profiler: &CommunityProfiler,
    data: &MappingData,
    rows: &[genomes::StrainDiagnostic],
) -> Result<()> {
    use std::io::Write;
    let path = format!("{}.levels_strains.tsv", profiler.cfg.outprefix);
    let mut out = std::io::BufWriter::new(std::fs::File::create(&path)?);
    writeln!(
        out,
        "taxid\tname\tgenome\tbest_fit\tem_share\tresolved_fragments\tstrain_identity_fragments\toebr\tcoebr\tfixed_chunk\tchunk_cutoff\tkept"
    )?;
    for row in rows {
        let (oebr, coebr, chunk) = row
            .coverage
            .as_ref()
            .map_or((f64::NAN, f64::NAN, f64::NAN), |c| {
                (c.oebr, c.coebr, c.fixed_chunk)
            });
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}",
            row.taxid,
            profiler.taxonomy.get_name(row.taxid).unwrap_or("unknown"),
            data.assembly_names[row.assembly_id as usize],
            row.best_fit,
            row.em_share,
            row.resolved,
            row.gated,
            oebr,
            coebr,
            chunk,
            data.min_chunk_breadth,
            row.kept,
        )?;
    }
    Ok(())
}

/// Mapping data for the tests of the level-wise engine and its genome level.
#[cfg(test)]
mod test_fixtures {
    use std::collections::HashMap;
    use std::sync::Arc;

    use metam::mapper::CandidateBoundary;

    use crate::profiler::community::{MappingData, ReadBundle};
    use crate::profiler::evidence::{EvidenceAlignment, MappingRole};
    use crate::profiler::ProfileRank;

    /// A single-end read of 100 bases aligned over its whole length.
    pub(super) fn row(taxid: u32, assembly: u32, identity: f32) -> EvidenceAlignment {
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

    pub(super) fn data(reads: Vec<Vec<EvidenceAlignment>>) -> MappingData {
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
            genus_of_species: Arc::new(HashMap::new()),
            terminal_ranks: Arc::new(HashMap::new()),
            skipped_reference_at_load: 0,
            marker_model: None,
            genus_rollup: None,
            unclassified_species: None,
            assembly_names: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::test_fixtures::{data, row};
    use super::*;

    #[test]
    fn species_rows_keep_what_their_genome_rows_do_not_hold() {
        // No genome reported: the species row keeps everything.
        assert_eq!(
            species_remainder(100.0, 15_000.0, std::iter::empty()),
            Some((100.0, 15_000.0))
        );
        // Two genomes hold most reads: the species row keeps the rest.
        assert_eq!(
            species_remainder(
                100.0,
                15_000.0,
                [(60.0, 9_000.0), (30.0, 4_500.0)].into_iter()
            ),
            Some((10.0, 1_500.0))
        );
        // Less than one read left: no species row.
        assert_eq!(
            species_remainder(100.0, 15_000.0, [(99.6, 14_940.0)].into_iter()),
            None
        );
        // Genome rows never leave a negative remainder.
        assert_eq!(
            species_remainder(100.0, 15_000.0, [(100.0, 15_100.0)].into_iter()),
            None
        );
    }

    /// Single-reference genomes 0..=12 of 1 Mb each; species 1 and 2 in genus 100, species 3 in genus 200.
    fn rollup() -> GenusRollup {
        GenusRollup {
            genus_of: HashMap::from([(1, 100), (2, 100), (3, 200)]),
            marker_only: HashMap::new(),
            layout: crate::profiler::rollup::GenomeLayout::new((0..=12).map(|id| (id, 1_000_000))),
            fractional: false,
        }
    }

    /// A genome's evidence from `n` placements spread along its 1 Mb.
    fn evidence(taxid: u32, kept: bool, n: u32) -> GenomeEvidence {
        let mut assembly = PresenceAssembly {
            asm_len: 1_000_000.0,
            sampling_scale: 1.0,
            ..PresenceAssembly::default()
        };
        for i in 0..n {
            let start = i * (1_000_000 / n);
            assembly.count += 1;
            assembly.reference_span_bases += 100;
            assembly.intervals.push(PresenceInterval {
                reference_id: 0,
                span: (start, start + 100),
                mate_span: (0, 0),
                fixed_chunk: (start, start + 100),
                flex_chunk: (start, start + 100),
            });
        }
        let coverage = calculate_coverage(&mut assembly, 1);
        GenomeEvidence {
            taxid,
            assembly,
            coverage: Some(coverage),
            kept,
        }
    }

    #[test]
    fn a_species_decision_is_its_strongest_genome() {
        // Species 1: genome 10 kept, genome 11 not kept but with more reads; species 2: genome 12, not kept.
        let mut genomes = HashMap::from([
            (10, evidence(1, true, 40)),
            (11, evidence(1, false, 80)),
            (12, evidence(2, false, 40)),
        ]);
        let decisions = species_decisions(&mut genomes, &rollup());
        let (one, genome) = decisions[&1];
        assert_eq!(genome, 10);
        assert!(one.standard);
        let (two, genome) = decisions[&2];
        assert_eq!(genome, 12);
        assert!(!two.standard);
        assert!(two.metrics.window_spread.is_some());
    }

    #[test]
    fn fragments_of_removed_species_of_one_genus_go_to_the_genus() {
        // 250 fragments only on species 1 (genus 100) and 30 on species 1 and 3 (two genera); nothing kept.
        let mut reads = (0..250).map(|_| vec![row(1, 1, 0.98)]).collect::<Vec<_>>();
        reads.extend((0..30).map(|_| vec![row(1, 1, 0.98), row(3, 3, 0.98)]));
        let data = data(reads);
        let fragments = (0..280)
            .map(|read| Fragment {
                read,
                candidates: if read < 250 {
                    vec![(1, 2.0)]
                } else {
                    vec![(1, 2.0), (3, 2.0)]
                },
                bases: 100.0,
            })
            .collect::<Vec<_>>();
        let group_of = vec![None; fragments.len()];
        let policy = EvidencePolicy {
            mode: EvidenceMode::Single,
            identity: 0.97,
            mapped_len: 50,
            fraction: 0.5,
            collect_exact_ties: false,
            collect_identities: false,
            collect_best_placements: false,
        };
        let decision = |coebr: f64| WitnessDecision {
            metrics: EvidenceCoverage {
                coebr,
                window_spread: Some(0.9),
                ..EvidenceCoverage::default()
            },
            standard: false,
            rescues: RescueKinds::default(),
            required_premap_guard: PremapChunkGuard::default(),
        };
        let rollup = rollup();
        let rows = orphan_genus_rows(
            &data,
            &fragments,
            &group_of,
            policy,
            &rollup,
            &HashMap::from([(1, decision(0.8)), (3, decision(0.8))]),
            &HashSet::new(),
        );
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].genus, 100);
        assert_eq!(rows[0].fragments, 250.0);
        // Read-1 bases over the genome's length, per fragment.
        assert!((rows[0].depth - 250.0 * 100.0 / 1_000_000.0).abs() < 1e-12);
        // A donor whose reads do not spread over its genome gives no row.
        let rows = orphan_genus_rows(
            &data,
            &fragments,
            &group_of,
            policy,
            &rollup,
            &HashMap::from([(1, decision(0.2))]),
            &HashSet::new(),
        );
        assert!(rows.is_empty());
    }

    #[test]
    fn species_rows_write_the_identity_columns_even_when_unknown() {
        let coverage = EvidenceCoverage::default();
        let row = |name: &str, identity| {
            profile_row(
                name.to_string(),
                1,
                "species",
                (10.0, 1.0),
                &coverage,
                std::sync::Arc::from("GCF_1"),
                identity,
            )
        };
        let mut writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(Vec::new());
        row("Known", Some((0.99, 0.995)))
            .write(&mut writer)
            .unwrap();
        row("Unknown", Some((f64::NAN, f64::NAN)))
            .write(&mut writer)
            .unwrap();
        let text = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        let lines = text.lines().collect::<Vec<_>>();
        assert!(lines.iter().all(|line| line.split('\t').count() == 17));
        assert!(lines[1].ends_with("GCF_1\tNA\tNA"));
    }

    #[test]
    fn genus_rows_write_the_same_columns_as_species_rows() {
        let coverage = EvidenceCoverage::default();
        let species = profile_row(
            "Species one".to_string(),
            1,
            "species",
            (10.0, 1.0),
            &coverage,
            std::sync::Arc::from("GCF_1"),
            Some((0.99, 0.995)),
        );
        let genus = genus_profile_row(
            "Genus".to_string(),
            &GenusRow {
                genus: 100,
                fragments: 250.0,
                equivalent: 250.0,
                depth: 0.5,
                mapped_len: 25_000.0,
                witness: coverage,
                donors: vec![3],
                reads: None,
            },
            std::sync::Arc::from("GCF_3"),
        );
        let mut writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(Vec::new());
        species.write(&mut writer).unwrap();
        genus.write(&mut writer).unwrap();
        let text = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        let lines = text.lines().collect::<Vec<_>>();
        assert!(lines.iter().all(|line| line.split('\t').count() == 17));
        assert!(lines[1].ends_with("GCF_3\tNA\tNA"));
    }

    #[test]
    fn a_removed_species_takes_the_fragments_it_ties_on_with_a_weaker_removed_sister() {
        let none = HashSet::new();
        // Tied with a removed sister that the first EM gives less weight.
        assert!(takes_group(&[(1, 0.82), (2, 0.18)], 1, &none));
        // The weaker sister does not take them.
        assert!(!takes_group(&[(1, 0.82), (2, 0.18)], 2, &none));
        // Nor fragments that another removed species fits better.
        assert!(!takes_group(&[(1, 0.05), (2, 0.95)], 1, &none));
        // Against a kept species it needs the resolving posterior.
        let kept = HashSet::from([3]);
        assert!(!takes_group(&[(1, 0.82), (3, 0.18)], 1, &kept));
        assert!(takes_group(&[(1, 0.82), (2, 0.1), (3, 0.08)], 1, &kept));
        // The sister judged after a kept species keeps losing the ties.
        let after = HashSet::from([1]);
        assert!(!takes_group(&[(1, 0.82), (2, 0.18)], 2, &after));
    }

    #[test]
    fn group_members_list_each_group_s_fragments_in_order() {
        let members = GroupMembers::new(&[Some(1), None, Some(0), Some(1), Some(2)], 3);
        assert_eq!(members.of(0), &[2]);
        assert_eq!(members.of(1), &[0, 3]);
        assert_eq!(members.of(2), &[4]);
    }

    #[test]
    fn the_last_look_keeps_species_with_reads_of_their_own() {
        // Species 1 is kept. Species 2 won 10 fragments: one shared with species
        // 1, one tied with removed species 5, eight its own. Species 3 won 10
        // fragments that all list species 1.
        let fragment = |candidates: &[(u32, f64)]| Fragment {
            read: 0,
            candidates: candidates.to_vec(),
            bases: 100.0,
        };
        let mut fragments = Vec::new();
        let mut resolved = Vec::new();
        for i in 0..10 {
            fragments.push(fragment(match i {
                0 => &[(2, 0.0), (1, 1.0)],
                1 => &[(2, 0.0), (5, 0.0)],
                _ => &[(2, 0.0)],
            }));
            resolved.push(Some(2));
            fragments.push(fragment(&[(3, 0.0), (1, 0.0)]));
            resolved.push(Some(3));
        }
        let group_of = (0..fragments.len()).map(Some).collect::<Vec<_>>();
        let em = likelihood_em::EmResult {
            weight: HashMap::from([(1, 100.0), (2, 10.0), (3, 10.0)]),
            posteriors: Vec::new(),
            iterations: 0,
        };
        let tests = SpeciesTests {
            multiple_tests: 1,
            chunk_cutoff: 0.5,
            ratio_bounds: (0.5, 0.5),
            markers: None,
            min_marker_genes: 3,
        };
        let measured = |fixed_chunk| EvidenceCoverage {
            fixed_chunk,
            oebr: 1.0,
            coebr: 1.0,
            ..EvidenceCoverage::default()
        };
        let genome = |taxid, kept| GenomeEvidence {
            taxid,
            assembly: PresenceAssembly::default(),
            coverage: Some(measured(1.0)),
            kept,
        };
        let start = HashMap::from([
            (10, genome(1, true)),
            (20, genome(2, false)),
            (30, genome(3, false)),
        ]);
        let first = FirstEm {
            fragments: &fragments,
            group_of: &group_of,
            em: &em,
            resolved: &resolved,
        };
        // Coverage of genome 20 or 30 from a list: whether the shared list
        // (one fragment for species 2), the own list (nine fragments) and the
        // tied list (one fragment) reach the cutoff.
        let look = |shared: f64, own: f64, tied: f64| {
            let mut genomes = start.clone();
            let mut counts = FragmentCounts::default();
            let rescued =
                own_evidence_look(&first, &tests, &mut counts, &mut genomes, |taxid, list| {
                    let id = if taxid == 2 { 20 } else { 30 };
                    let chunk = match (taxid, list.len()) {
                        (2, 9) => own,
                        (2, _) if list[0] == 0 => shared,
                        (2, _) => tied,
                        _ => 0.0,
                    };
                    vec![(id, measured(chunk))]
                });
            (
                rescued,
                genomes[&20].kept,
                counts.own.get(&2).copied(),
                counts.own.get(&3).copied(),
            )
        };
        // Own reads covering the genome keep species 2; species 3 owns none.
        assert_eq!(
            look(0.0, 1.0, 0.0),
            (HashMap::from([(2, false)]), true, Some(9), Some(0))
        );
        // Reads shared with the kept species along the genome: stays with it.
        assert_eq!(look(1.0, 1.0, 0.0), (HashMap::new(), false, None, Some(0)));
        // Own reads in one region only: stays removed.
        assert_eq!(look(0.0, 0.0, 0.0).0, HashMap::new());
        // Ties along the genome: kept, reported at its genus.
        assert_eq!(look(0.0, 1.0, 1.0).0, HashMap::from([(2, true)]));
    }

    #[test]
    fn a_fragment_without_a_species_match_counts_for_the_genera_of_its_species() {
        // One read below the species identity on two species of genus 100 and
        // one of genus 200.
        let data = data(vec![vec![
            row(2, 20, 0.88),
            row(3, 30, 0.90),
            row(4, 40, 0.85),
        ]]);
        let policy = |identity| EvidencePolicy {
            mode: EvidenceMode::Single,
            identity,
            mapped_len: 50,
            fraction: 0.6,
            collect_exact_ties: false,
            collect_identities: false,
            collect_best_placements: false,
        };
        let genus_of = HashMap::from([(2, 100), (3, 100), (4, 200)]);
        let mut scratch = CandidateScratch::default();
        let alignments = &data.read_alignments[0].alignments;
        // No candidate at the species identity.
        let species = scratch.build(alignments, policy(0.95));
        assert!(fragment_of(0, species, policy(0.95), Some).is_none());
        // At the genus identity, one candidate per genus with its species'
        // fewest errors, fewest first.
        let genus = scratch.build(alignments, policy(0.80));
        let (fragment, _) = fragment_of(0, genus, policy(0.80), |t| genus_of.get(&t).copied())
            .expect("candidates at the genus identity");
        let units = fragment.candidates.iter().map(|c| c.0).collect::<Vec<_>>();
        assert_eq!(units, vec![100, 200]);
        assert!((fragment.candidates[0].1 - 10.0).abs() < 0.5);
        assert!((fragment.candidates[1].1 - 15.0).abs() < 0.5);
    }

    #[test]
    fn species_map_to_the_family_of_their_genus() {
        // nodes: taxid -> (parent, rank). Genus 10 lies in family 100, genus 20
        // directly under an order without a family.
        let node = |parent: u32, rank: &str| (parent, rank.to_string());
        let taxonomy = Taxonomy::from_parts(
            HashMap::from([
                (1, node(1, "no rank")),
                (1000, node(1, "order")),
                (100, node(1000, "family")),
                (10, node(100, "genus")),
                (20, node(1000, "genus")),
                (11, node(10, "species")),
                (12, node(10, "species")),
                (21, node(20, "species")),
            ]),
            HashMap::new(),
            HashMap::new(),
            HashSet::new(),
        );
        let families = family_map(&taxonomy, &HashMap::from([(11, 10), (12, 10), (21, 20)]));
        assert_eq!(families, HashMap::from([(11, 100), (12, 100)]));
    }

    #[test]
    fn a_marker_genome_is_measured_gene_by_gene_and_needs_enough_genes() {
        use crate::profiler::markers::{FragmentGeometry, FragmentLayout};
        // Genome 7: five marker genes (references 1-5) of 1,000 bases; genome 8
        // is a whole genome of 5,000 bases.
        let geometry = FragmentGeometry::new(FragmentLayout::Reads, 100, 100, 50, 0.5);
        let genes = [1_000u32; 5];
        let model = MarkerModel::from_gene_lengths(&[(7, &genes[..])], geometry);
        let tests = SpeciesTests {
            multiple_tests: 1,
            chunk_cutoff: 0.1,
            ratio_bounds: (0.0, 0.0),
            markers: Some(Arc::new(model)),
            min_marker_genes: 3,
        };
        let assembly = |genes: &[u32]| {
            let mut assembly = PresenceAssembly {
                asm_len: 5_000.0,
                ..PresenceAssembly::default()
            };
            for &gene in genes {
                for start in [100, 500] {
                    assembly.count += 1;
                    assembly.reference_span_bases += 100;
                    assembly.intervals.push(PresenceInterval {
                        reference_id: gene,
                        span: (start, start + 100),
                        mate_span: (0, 0),
                        fixed_chunk: (gene, gene),
                        flex_chunk: (gene, gene),
                    });
                }
            }
            assembly
        };
        // Two of five genes hit: chunk breadth 0.4, but fewer than three genes.
        let two = tests.measure(7, &mut assembly(&[1, 2]));
        assert_eq!(two.marker_genes_hit, Some(2));
        assert!((two.fixed_chunk - 0.4).abs() < 1e-9);
        assert!(!tests.covers(&two));
        let three = tests.measure(7, &mut assembly(&[1, 2, 3]));
        assert_eq!(three.marker_genes_hit, Some(3));
        assert!(tests.covers(&three));
        // A genome without marker genes is measured over its sequence.
        assert_eq!(
            tests.measure(8, &mut assembly(&[1, 2])).marker_genes_hit,
            None
        );
    }

    #[test]
    fn a_species_on_a_marker_genome_takes_the_depth_of_its_whole_genome() {
        use crate::profiler::evidence::EvidenceAlignment;
        use crate::profiler::markers::{FragmentGeometry, FragmentLayout};
        // Species 1 on marker genome 0 (one gene of 1,000 bases), species 2 on
        // whole genome 1; three reads on the gene, one tied with genome 1.
        let geometry = FragmentGeometry::new(FragmentLayout::Reads, 100, 100, 50, 0.5);
        let model = MarkerModel::from_gene_lengths(&[(0, &[1_000][..])], geometry);
        let gene = |taxid| EvidenceAlignment {
            marker_gene_len: 1_000,
            asm_len: 1_000.0,
            ..row(taxid, 0, 1.0)
        };
        let mut data = data(vec![
            vec![gene(1)],
            vec![gene(1)],
            vec![gene(1), row(2, 1, 1.0)],
        ]);
        data.marker_model = Some(Arc::new(model));
        let policy = EvidencePolicy {
            mode: EvidenceMode::Single,
            identity: 0.8,
            mapped_len: 50,
            fraction: 0.5,
            collect_exact_ties: false,
            collect_identities: false,
            collect_best_placements: false,
        };
        let fragments = (0..3)
            .map(|read| Fragment {
                read,
                candidates: Vec::new(),
                bases: 100.0,
            })
            .collect::<Vec<_>>();
        let shared = SharedReads {
            groups: ReadGroups::default(),
            group_of: vec![Some(0), Some(0), Some(1)],
            group_bases: Vec::new(),
        };
        let em = likelihood_em::EmResult {
            weight: HashMap::new(),
            posteriors: vec![vec![(1, 1.0)], vec![(1, 0.5), (2, 0.5)]],
            iterations: 0,
        };
        let reported = HashMap::from([(1, 0), (2, 1)]);
        let model = data.marker_model.clone().unwrap();
        let evidence = marker_evidence(
            &data,
            policy,
            &fragments,
            (&shared, &em),
            &reported,
            (&model, Some(false)),
        );
        let depth = &evidence.depth;
        // 2.5 fragments of read length 100 over the placements the gene offers;
        // species 2 is on a whole genome and keeps bases over length.
        let placements = model.fragment_placements(0).unwrap();
        assert!((depth[&1] - 2.5 * 100.0 / placements).abs() < 1e-12);
        assert!(!depth.contains_key(&2));
        // The two fragments resolved to species 1 count for the estimate; the
        // shared one does not.
        assert_eq!(evidence.genome_reads.len(), 1);
        assert_eq!(evidence.genome_reads[&0].reads(), 2);
    }
}

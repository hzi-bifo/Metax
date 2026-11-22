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
//! that a removed sister no longer splits the fragments it ties on
//! ([`second_look`]); a removed species is also kept when the fragments it won
//! that list no kept species cover one of its genomes on their own
//! ([`own_evidence_look`]), and a kept virus is removed when its own fragments
//! do not cover the part of its genome that no kept bacterium or archaeon
//! explains ([`prophage_look`]). Fragments with no candidate at the species
//! identity are offered to their families; on noisy long reads, those whose
//! best match is a kept species count for it instead ([`credited_reads`]). A
//! second EM shares the fragments among the kept species, and abundance is the
//! resulting depth on each species' best-fit genome; a fragment on a marker
//! genome counts as one on the whole genome ([`fragment_depth`]). Read counts
//! are given as full-index equivalents, scaled up for genomes indexed in
//! sampled segments. On a database with such genomes, a species whose reads
//! are far less identical than the sample's typical species is left out beside
//! a kept congener and otherwise reported at its genus; with
//! `--relative-identity-gap` (on by default for single-end reads on such a
//! database) a species is left out beside a deeper congener whose reads are
//! that much more identical ([`identity_relatives`]). A genome represented by
//! marker genes is measured gene by gene ([`MarkerModel`]); the tests are the
//! same. The genomes whose references are labelled only above species (a
//! genus, a family or higher) form one unit per labelled taxon, judged and
//! reported like a species at that rank, without the identity rules. Named
//! units come first: on the second look they are judged before the labelled
//! unit of their taxon and do not lose to it the reads they share with it
//! ([`named_units_first`]); after the looks, a kept labelled unit must cover
//! its genome beyond the part that the kept named units of its taxon explain,
//! or it is removed and its reads go to them ([`named_first_look`]). Its row
//! counts the taxon's genus or family row as well. `--host` drops the rows of taxa that the pathogen table lists for
//! other hosts only and annotates the others ([`host`]), and `--min-reads` and
//! `--min-breadth` drop a taxon whose rows hold fewer reads (a species with its
//! genome rows) or whose breadth is lower ([`drop_below_minimums`]), both
//! before abundance is normalised. The profile opens with two `#` lines: the
//! run's metadata and the column names. With `--save-raw-profile` the engine
//! also writes `<prefix>.decisions.tsv`: how every taxon it judged, at every
//! rank, was decided, with the tests behind it ([`tables`]), and with
//! `--save-classification` `<prefix>.classify.txt.gz`: each fragment's shares
//! of the profile's rows ([`classify`]).

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use anyhow::Result;
use rayon::prelude::*;

use super::{add_marker_genome_reads, fragment_depth, CommunityProfiler, MappingData};
use crate::profiler::evidence::{
    calculate_coverage, standard_ratio_lower_bounds, CandidatePlacement, CandidateScratch,
    EvidenceCoverage, EvidenceMode, EvidencePolicy, GeneHit, PremapChunkGuard, PresenceAssembly,
    PresenceInterval, RescueKinds, TerminalCandidate, WitnessDecision, STANDARD_RATIO_UPPER,
};
use crate::profiler::identity::MIN_JUDGED_SPECIES;
use crate::profiler::likelihood_em::{self, ReadGroups};
use crate::profiler::markers::MarkerModel;
use crate::profiler::output::{BestFitReport, ProfileRow};
use crate::profiler::pathogen::load_pathogen_table;
use crate::profiler::rollup::{
    GenomeLayout, GenusRollup, GenusRow, OrphanedFragments, MIN_GENUS_FRAGMENTS, MIN_WINDOW_SPREAD,
};
use crate::profiler::unclassified::GenomeReads;
use crate::profiler::{
    equivalent_read_count_passes_minimum, parse_reference_taxonomy, ProfileRank,
};
use crate::taxonomy::Taxonomy;

mod classify;
mod genomes;
mod host;
mod tables;

use tables::{
    Decisions, FinalRows, Names, SpeciesDecision, StrainLevel, TableHeader, HOST_COLUMNS,
    PROFILE_COLUMNS,
};

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
/// Lower bound of a species' base-level breadth ratio (OEBR), below the default
/// engine's: amplified libraries cover a genome unevenly base by base while its
/// chunks stay evenly covered, which the chunk-level ratio (COEBR) still tests
/// at its own bound.
const SPECIES_MIN_OEBR: f64 = 0.6;
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

impl Fragment {
    /// The candidates tied with the best: candidates are sorted by errors, and
    /// those within the tie margin of the best carry exactly its errors.
    fn best_tied(&self) -> &[(u32, f64)] {
        let Some(&(_, fewest)) = self.candidates.first() else {
            return &[];
        };
        let tied = self.candidates.iter().take_while(|c| c.1 == fewest).count();
        &self.candidates[..tied]
    }
}

/// Fragments shared among a set of species: the EM's read groups, each
/// fragment's group, and each group's bases.
struct SharedReads {
    groups: ReadGroups,
    group_of: Vec<Option<usize>>,
    group_bases: Vec<f64>,
}

impl SharedReads {
    /// Each unit's share of the fragments' bases in `em`, the EM over these
    /// read groups.
    fn unit_bases(&self, em: &likelihood_em::EmResult) -> HashMap<u32, f64> {
        let mut bases: HashMap<u32, f64> = HashMap::new();
        for (posteriors, &group_total) in em.posteriors.iter().zip(&self.group_bases) {
            for &(unit, posterior) in posteriors {
                *bases.entry(unit).or_insert(0.0) += group_total * posterior;
            }
        }
        bases
    }
}

/// The unit a read group resolves to: the first of its `posteriors` at
/// [`RESOLVE_POSTERIOR`] or above.
fn resolved_unit(posteriors: &[(u32, f64)]) -> Option<u32> {
    resolved_share(posteriors).map(|p| p.0)
}

/// The unit a read group resolves to, with its posterior.
fn resolved_share(posteriors: &[(u32, f64)]) -> Option<(u32, f64)> {
    posteriors
        .iter()
        .find(|p| p.1 >= RESOLVE_POSTERIOR)
        .copied()
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
    // With `--save-classification` (not on a database of sampled genomes),
    // each level keeps how it shares the fragments.
    let classify_path = profiler.classification_path(&data);
    let keep_shares = classify_path.is_some();
    // Units labelled above species: the taxa of references whose most
    // specific label is a genus, family or higher taxon. Their genomes form
    // one unit per taxon, judged like a species.
    let labels: HashMap<u32, ProfileRank> = data
        .terminal_ranks
        .iter()
        .filter(|(_, &rank)| rank != ProfileRank::Species)
        .map(|(&taxid, &rank)| (taxid, rank))
        .collect();
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
    let first_em = likelihood_em::run(&all.groups, penalty, cfg.max_em_iterations, EM_TOLERANCE);
    let resolved_taxid = first_em
        .posteriors
        .iter()
        .map(|posteriors| resolved_unit(posteriors))
        .collect::<Vec<_>>();
    let first = FirstEm {
        fragments: &fragments,
        group_of: &all.group_of,
        em: &first_em,
        resolved: &resolved_taxid,
    };
    let mut counts = first.counts();
    let resolved_fragments = first
        .won()
        .map(|(fragment, taxid)| (fragments[fragment].read, taxid))
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
    drop(resolved_fragments);

    let (oebr_default, coebr_default) =
        standard_ratio_lower_bounds(ProfileRank::Species, policy.mode);
    let min_coebr = cfg.min_coebr.unwrap_or(coebr_default);
    // The genome level keeps the default engine's bounds.
    let genome_bounds = (cfg.min_oebr.unwrap_or(oebr_default), min_coebr);
    let tests = SpeciesTests {
        multiple_tests: genomes.len().max(1),
        chunk_cutoff: data.min_chunk_breadth,
        ratio_bounds: (cfg.min_oebr.unwrap_or(SPECIES_MIN_OEBR), min_coebr),
        markers: data.marker_model.clone(),
        min_marker_genes: cfg.min_marker_genes,
    };
    genomes.par_iter_mut().for_each(|(&id, genome)| {
        let won = counts.won.get(&genome.taxid).copied().unwrap_or(0);
        let candidates = counts.candidates.get(&genome.taxid).copied().unwrap_or(0);
        tests.judge(id, genome, won, candidates);
    });
    let ancestors = |taxid: u32| {
        profiler
            .taxonomy
            .get_parents(taxid)
            .into_iter()
            .map(|(ancestor, _, _)| ancestor)
            .collect::<Vec<_>>()
    };
    let taken = second_look(
        &data,
        policy,
        &first,
        &tests,
        &mut counts,
        &mut genomes,
        (&labels, &ancestors),
    );
    log::info!(
        target: "PROFILER",
        "Level-wise engine: second look kept {} removed species on {} fragments their removed sisters tied on.",
        taken.len(),
        taken.values().map(Vec::len).sum::<usize>()
    );
    let rescued = own_evidence_look(&first, &tests, &mut counts, &mut genomes, |taxid, list| {
        coverage_of(&data, policy, &fragments, &tests, taxid, list)
    });
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
    let removed_viruses = prophage_look(
        &first,
        &taken,
        &tests,
        &mut genomes,
        (
            |taxid| domain(&profiler.taxonomy, taxid),
            (&labels, &ancestors),
        ),
        |taxid, list| evidence_of(&data, policy, &fragments, taxid, list),
    );
    log::info!(
        target: "PROFILER",
        "Level-wise engine: {} kept virus species removed: their own fragments do not cover the rest of their \
         genome, beyond the part that kept cellular species explain.",
        removed_viruses.len()
    );
    let named_first = if labels.is_empty() {
        HashMap::new()
    } else {
        let named_first = named_first_look(
            &first,
            &taken,
            &tests,
            &mut genomes,
            (&labels, &ancestors),
            |taxid, list| evidence_of(&data, policy, &fragments, taxid, list),
        );
        log::info!(
            target: "PROFILER",
            "Level-wise engine: {} kept units labelled above species removed: their own fragments do not cover the \
             rest of their genome, beyond the part that kept named units of their taxon explain.",
            named_first.len()
        );
        named_first
    };
    // The first EM is not read again.
    let first_iterations = first_em.iterations;
    drop((all, first_em, resolved_taxid));
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
    let decisions = rollup.map_or_else(HashMap::new, |rollup| species_decisions(&genomes, rollup));
    // The placements on the measured genomes are not read again.
    for genome in genomes.values_mut() {
        genome.assembly.intervals = Vec::new();
    }
    let decision_metrics: HashMap<u32, WitnessDecision> = decisions
        .iter()
        .map(|(&taxid, &(decision, _))| (taxid, decision))
        .collect();
    let even = rollup.map_or_else(Vec::new, |rollup| {
        rollup.even_coverage_species(
            &decision_metrics,
            |taxid| data.terminal_ranks.get(&taxid).copied(),
            data.min_chunk_breadth,
        )
    });
    for &taxid in &even {
        best_fit.insert(taxid, decisions[&taxid].1);
    }
    // What became of each measured species, for the decisions table (empty
    // when it is not written); the identity rules and the genus reports below
    // refine it.
    let mut outcomes = if cfg.keep_raw {
        tables::judged_species(
            &genomes,
            &best_fit,
            &tables::Looks {
                second_look: &taken,
                own_evidence: &rescued,
                prophage: &removed_viruses,
                named_first: &named_first,
                even_coverage: &even,
            },
        )
    } else {
        HashMap::new()
    };

    // Second EM over the kept species.
    let kept = share_reads(&fragments, |taxid| best_fit.contains_key(&taxid), penalty);
    let em = likelihood_em::run(&kept.groups, penalty, cfg.max_em_iterations, EM_TOLERANCE);
    let without_kept = kept.group_of.iter().filter(|g| g.is_none()).count();
    log::info!(
        target: "PROFILER",
        "Level-wise engine: {} species kept, {} read groups, EM {} + {} iterations; {} fragments without a kept species.",
        best_fit.len(),
        kept.groups.len(),
        first_iterations,
        em.iterations,
        without_kept
    );
    let (genus_rows, orphans) = match rollup {
        Some(rollup) => {
            let accepted: HashSet<u32> = best_fit.keys().copied().collect();
            let (rows, orphans) = orphan_genus_rows(
                &data,
                &fragments,
                &kept.group_of,
                policy,
                (rollup, &labels),
                (&decision_metrics, &accepted),
                keep_shares,
            );
            log::info!(
                target: "PROFILER",
                "Level-wise engine: {} species kept for even coverage; {} genus rows from {:.0} fragments of removed species.",
                even.len(),
                rows.len(),
                // A float sum of no rows is -0.0, printed "-0".
                rows.iter().fold(0.0, |total, row| total + row.fragments)
            );
            (rows, orphans)
        }
        None => (Vec::new(), Vec::new()),
    };
    // Fragments with no candidate at the species identity are built again with
    // the default engine's genus fallback gates.
    let genus_policy = EvidencePolicy {
        identity: cfg.genus_identity,
        mapped_len: super::genus_fallback_min_mapped_len(cfg),
        fraction: super::GENUS_FALLBACK_MIN_FRACTION,
        collect_exact_ties: false,
        collect_identities: false,
        collect_best_placements: false,
        ..policy
    };
    // On noisy long reads many reads of a present species align below the
    // species gates; those whose best match is a kept species are its reads,
    // not the family level's. Other reads keep the family level as it is.
    let noisy_kept: HashSet<u32> = if cfg.noisy_long_reads {
        best_fit.keys().copied().collect()
    } else {
        HashSet::new()
    };
    // Family rows from those fragments: they come from genomes missing from
    // the database, which their reads cannot place in a genus.
    let (family_level, family_candidates, family_shares) = match rollup {
        Some(rollup) => {
            let mut family_of = family_map(&profiler.taxonomy, &data.genus_of_species);
            // A unit labelled at a family counts for itself.
            family_of.extend(
                labels
                    .iter()
                    .filter(|(_, &rank)| rank == ProfileRank::Family)
                    .map(|(&taxid, _)| (taxid, taxid)),
            );
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
            let UnmatchedResult {
                rows,
                candidates,
                shares,
            } = unmatched_rows(
                profiler,
                &data,
                &unmatched,
                &UnmatchedLevel {
                    policy: genus_policy,
                    penalty,
                    tests: &tests,
                    layout: &rollup.layout,
                    unit_of: &family_of,
                    excluded: &excluded,
                    kept: &noisy_kept,
                    rank: "family",
                    keep_shares,
                },
            );
            log::info!(
                target: "PROFILER",
                "Level-wise engine: {} fragments without a species-level candidate; {} family rows from {:.0} of them.",
                unmatched.len(),
                rows.len(),
                rows.iter()
                    .fold(0.0, |total, (_, row)| total + row.equivalent_mapcount)
            );
            (rows, candidates, shares)
        }
        None => (Vec::new(), Vec::new(), None),
    };
    // Each species' share of the fragments' bases.
    let mut bases = kept.unit_bases(&em);

    // Genomes within each kept species, for the fragments the species-level EM
    // gives to one species, with that species' share of each.
    let species_of = kept
        .group_of
        .iter()
        .map(|group| group.and_then(|g| resolved_share(&em.posteriors[g])))
        .collect::<Vec<_>>();
    let genome_level = genomes::GenomeLevel {
        policy,
        error_rate,
        report_strains: cfg.genome_level,
        keep_shares,
        chunk_cutoff: data.min_chunk_breadth,
        ratio_bounds: genome_bounds,
        max_iterations: cfg.max_em_iterations,
        tolerance: EM_TOLERANCE,
    };
    let genome_result = genomes::resolve(&data, &fragments, &species_of, &labels, &genome_level);
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
    let MarkerEvidence {
        depth: mut marker_depth,
        genome_reads,
    } = data
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
    // Reads from the unindexed parts of a sampled genome land on its relatives:
    // the margin on a database with sampled genomes, and the relative identity
    // guard where it is on.
    let margin = cfg.identity_margin.filter(|_| data.has_scaled_genomes);
    let gap = cfg.relative_identity_gap;
    let (left_out, low_at_genus) = if margin.is_some() || gap.is_some() {
        let species = reported
            .iter()
            .filter(|(taxid, _)| !labels.contains_key(taxid))
            .filter_map(|(&taxid, &id)| {
                let identity = species_genomes.get(&taxid)?.mean_identity;
                let depth = marker_depth.get(&taxid).copied().unwrap_or_else(|| {
                    bases.get(&taxid).copied().unwrap_or(0.0)
                        / genomes[&id].assembly.asm_len.max(1.0)
                });
                identity.is_finite().then_some((taxid, (identity, depth)))
            })
            .collect::<HashMap<_, _>>();
        let found = identity_relatives(
            &data.genus_of_species,
            &species,
            baseline_identity,
            (margin, gap),
        );
        log::info!(
            target: "PROFILER",
            "Level-wise engine: {} species left out and {} reported at their genus by the identity rules for \
             relatives (margin {}, relative identity gap {}).",
            found.0.len(),
            found.1.len(),
            margin.map_or_else(|| "off".to_string(), |m| m.to_string()),
            gap.map_or_else(|| "off".to_string(), |g| g.to_string())
        );
        found
    } else {
        (HashSet::new(), HashSet::new())
    };
    // The fragments that the family level leaves to a kept species, too noisy
    // for the species gates, count for it.
    let credited = if rollup.is_some() && cfg.noisy_long_reads {
        let credited = credited_reads(
            &data,
            genus_policy,
            &unmatched,
            (&em.weight, &reported),
            keep_shares,
        );
        log::info!(
            target: "PROFILER",
            "Level-wise engine: {:.0} fragments without a species-level candidate credited to the kept species that match them best.",
            credited.reads.values().sum::<f64>()
        );
        credited
    } else {
        Credited::default()
    };
    for (taxid, added) in &credited.bases {
        *bases.entry(*taxid).or_insert(0.0) += added;
    }
    for (taxid, added) in &credited.marker_depth {
        *marker_depth.entry(*taxid).or_insert(0.0) += added;
    }
    let reads_of = |taxid: u32| {
        em.weight.get(&taxid).copied().unwrap_or(0.0)
            + credited.reads.get(&taxid).copied().unwrap_or(0.0)
    };

    // Species rows, each followed by the genomes reported on their own. The
    // rows of a species are exclusive: a genome row counts each fragment by the
    // species' share of it, the species row keeps what no reported genome
    // holds and is left out when less than one read is left.
    let mut blocks: Vec<TaxonBlock> = Vec::new();
    // Species reported at their genus: (genus, species, row); for the
    // decisions table and the per-read classifications, (genus, species,
    // reads).
    let mut genus_parts: Vec<(u32, u32, ProfileRow)> = Vec::new();
    let mut genus_members: Vec<(u32, u32, f64)> = Vec::new();
    for (&taxid, &assembly_id) in &reported {
        if left_out.contains(&taxid) {
            outcomes
                .entry(taxid)
                .and_modify(|d| *d = SpeciesDecision::LeftOut);
            continue;
        }
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
            (
                g.mean_identity,
                estimated_ani(baseline_identity, g.mean_identity),
            )
        }));
        // A unit labelled above species stays at its own rank.
        let rank = labels.get(&taxid).map_or("species", |rank| rank.as_str());
        let reported_at_genus = !labels.contains_key(&taxid)
            && (at_genus.contains(&taxid) || low_at_genus.contains(&taxid));
        if let Some(&genus) = data
            .genus_of_species
            .get(&taxid)
            .filter(|_| reported_at_genus)
        {
            outcomes.entry(taxid).and_modify(|d| {
                *d = if at_genus.contains(&taxid) {
                    SpeciesDecision::AtGenusTie
                } else {
                    SpeciesDecision::AtGenusIdentity
                }
            });
            let reads = reads_of(taxid);
            if cfg.keep_raw || keep_shares {
                genus_members.push((genus, taxid, reads));
            }
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
            reads_of(taxid),
            bases.get(&taxid).copied().unwrap_or(0.0),
            reported.iter().map(|s| (s.reads, s.bases)),
        )
        .map(|(reads, species_bases)| {
            profile_row(
                species_name.to_string(),
                taxid,
                rank,
                (reads * scale, depth(species_bases)),
                coverage,
                data.assembly_names[assembly_id as usize].clone(),
                identity,
            )
        });
        let mut strains = Vec::new();
        for strain in reported {
            let strain_ani = estimated_ani(strain_baseline, strain.mean_identity);
            let accession = data.assembly_names[strain.assembly_id as usize].clone();
            let (name, strain_taxid) =
                strain_label(&profiler.taxonomy, &strain_taxids, taxid, &accession);
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
            breadth: coverage.breadth,
            unit: true,
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
                profile_row(
                    profiler
                        .taxonomy
                        .get_name(row.genus)
                        .unwrap_or("unknown")
                        .to_string(),
                    row.genus,
                    "genus",
                    (row.fragments, row.depth),
                    &row.witness,
                    data.assembly_names[decisions[&row.donors[0]].1 as usize].clone(),
                    Some((f64::NAN, f64::NAN)),
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
    // A genus or family row of a taxon that a kept unit is labelled with is
    // counted in that unit's row.
    for (taxid, row) in genus_blocks.into_iter().chain(family_level) {
        match blocks
            .iter_mut()
            .find(|block| block.unit && block.taxid == taxid)
        {
            Some(block) => block.count(row),
            None => blocks.push(TaxonBlock {
                taxid,
                breadth: row.breadth,
                taxon: Some(row),
                genomes: Vec::new(),
                unit: false,
            }),
        }
    }
    // `--host` drops rows first, then `--min-reads` and `--min-breadth`.
    let host = profiler.host_taxid();
    let mut dropped = Dropped::default();
    if let Some(host) = host {
        let table = cfg
            .pathogen_host
            .as_deref()
            .map(load_pathogen_table)
            .transpose()?;
        let rows = |blocks: &[TaxonBlock]| {
            blocks
                .iter()
                .map(|b| usize::from(b.taxon.is_some()) + b.genomes.len())
                .sum::<usize>()
        };
        let before = rows(&blocks);
        host::HostFilter {
            host,
            table: table.as_ref(),
            taxonomy: &profiler.taxonomy,
        }
        .apply(&mut blocks, &mut dropped);
        log::info!(
            target: "PROFILER",
            "Level-wise engine: {} rows dropped by --host {host} ({} taxa entirely): not pathogens of the host.",
            before - rows(&blocks),
            dropped.taxa.len()
        );
    }
    let below = drop_below_minimums(
        &mut blocks,
        &Minimums {
            reads: cfg.min_reads,
            breadth: cfg.breadth,
        },
    );
    if !below.is_empty() {
        log::info!(
            target: "PROFILER",
            "Level-wise engine: {} taxa dropped by --min-reads or --min-breadth; their reads go to no other taxon.",
            below.len()
        );
    }
    dropped.taxa.extend(
        below
            .into_iter()
            .map(|(taxid, below)| (taxid, Removal::BelowMinimum(below))),
    );
    set_abundance(&mut blocks);
    // By the taxon's total depth, then taxid.
    blocks.sort_by(|a, b| b.depth().total_cmp(&a.depth()).then(a.taxid.cmp(&b.taxid)));
    // With `--save-raw-profile`, what each taxon holds in the profile.
    let finals = cfg.keep_raw.then(|| FinalRows::of(&blocks));
    if let Some(path) = &classify_path {
        let written = classify::Classifications {
            fragments: &fragments,
            group_of: &kept.group_of,
            posteriors: &em.posteriors,
            genome_level: species_genomes,
            genus_members: &genus_members,
            rollup_rows: &genus_rows,
            orphans: &orphans,
            credited: &credited.fragments,
            families: family_shares.as_ref(),
        }
        .write(path, &blocks, &data)?;
        log::info!(
            target: "PROFILER",
            "Level-wise engine: {written} fragments written to the per-read classifications."
        );
    }
    let rows = blocks
        .into_iter()
        .flat_map(|block| block.taxon.into_iter().chain(block.genomes))
        .collect::<Vec<_>>();

    let unclassified = profiler.estimate_unclassified(&data, &rows, &genome_reads);
    let header = TableHeader::of(cfg, host);
    let columns = if host.is_some() {
        [&PROFILE_COLUMNS[..], &HOST_COLUMNS[..]].concat()
    } else {
        PROFILE_COLUMNS.to_vec()
    };
    profiler.write_profile(
        &header.lines("profile", &columns),
        &rows,
        unclassified,
        true,
    )?;
    if let Some(finals) = finals {
        let decisions_table = Decisions {
            genomes: &genomes,
            counts: &counts,
            tests: &tests,
            outcomes: &outcomes,
            reported: &reported,
            strains: cfg.genome_level.then(|| StrainLevel {
                rows: &genome_result.diagnostics,
                species: species_genomes,
                tests: SpeciesTests {
                    multiple_tests: 1,
                    chunk_cutoff: genome_level.chunk_cutoff,
                    ratio_bounds: genome_level.ratio_bounds,
                    markers: None,
                    min_marker_genes: 0,
                },
                strain_taxids: &strain_taxids,
            }),
            rollup_rows: &genus_rows,
            rollup_genomes: &decisions,
            genus_members: &genus_members,
            families: &family_candidates,
            dropped: &dropped,
            labels: &labels,
            named_first: &named_first,
        };
        let names = Names {
            taxonomy: &profiler.taxonomy,
            genomes: &data.assembly_names,
        };
        let path = format!("{}.decisions.tsv", profiler.profile_prefix());
        tables::write_file(&path, |out| {
            decisions_table.write(out, &header, &names, &finals)
        })?;
    }
    Ok(())
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

/// Estimated errors of a candidate's placement, counting the bases it aligns
/// less than the fragment's longest alignment (`longest`).
fn candidate_errors(placement: &CandidatePlacement, longest: f64) -> f64 {
    placement.pair_errors
        + UNALIGNED_ERROR_WEIGHT * (longest - f64::from(placement.both_mapped)).max(0.0)
}

/// The margin within which a candidate's errors tie with the best: short
/// reads by [`SHORT_TIE_ERRORS`], long reads by `long_share` of `bases`.
fn tie_margin(mode: EvidenceMode, bases: f64, long_share: f64) -> f64 {
    match mode {
        EvidenceMode::Long => long_share * bases,
        EvidenceMode::Paired | EvidenceMode::Single => SHORT_TIE_ERRORS,
    }
}

/// Bases of a fragment's longest candidate alignment.
fn longest_alignment(candidates: &[TerminalCandidate]) -> f64 {
    candidates
        .iter()
        .map(|c| f64::from(c.best.both_mapped))
        .fold(0.0, f64::max)
}

/// Read `read` as a fragment with its candidate units and their errors, with
/// the errors and aligned bases of its best candidate; `None` without one. A
/// candidate is a species, or the taxon above species that its references are
/// labelled with. `unit` gives the unit a candidate counts for (the candidate
/// itself, or its family); a unit keeps the fewest errors of its candidates.
fn fragment_of(
    read: usize,
    candidates: &[TerminalCandidate],
    policy: EvidencePolicy,
    unit: impl Fn(u32) -> Option<u32>,
) -> Option<(Fragment, (f32, f32))> {
    let longest = longest_alignment(candidates);
    let mut units = candidates
        .iter()
        .filter_map(|c| {
            let errors = candidate_errors(&c.best, longest);
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
    let margin = tie_margin(policy.mode, f64::from(units[0].2), LONG_TIE_SHARE);
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
        let slot = shared.groups.add(candidates);
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
    best_matches.par_sort_by(|a, b| (a.0 / a.1.max(1.0)).total_cmp(&(b.0 / b.1.max(1.0))));
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

/// Estimated ANI of the sample's strain to a genome: its fragments' mean
/// `identity` to the genome relative to the sample's `baseline` identity.
fn estimated_ani(baseline: f64, identity: f64) -> f64 {
    (1.0 - (baseline - identity)).min(1.0)
}

/// The identity rules for relatives, over the kept `species` with their mean
/// identity to their best-fit genome and their depth. With a `margin` (on a
/// database whose genomes are indexed in sampled segments, where the reads of a
/// genome's unindexed parts land on the segments of its relatives), a species
/// whose reads are more than the margin less identical than the sample's
/// `baseline` is left out (first set) when a congener within the margin is
/// kept, whose read count already stands for its whole genome, and is
/// otherwise reported at its genus (second set); one without a genus stays,
/// and the margin judges nothing with fewer than [`MIN_JUDGED_SPECIES`]
/// species. With a `gap` (`--relative-identity-gap`), a species is left out
/// when a congener with more depth has a mean identity at least the gap
/// higher, as its reads are better explained as that relative's.
fn identity_relatives(
    genus_of: &HashMap<u32, u32>,
    species: &HashMap<u32, (f64, f64)>,
    baseline: f64,
    (margin, gap): (Option<f64>, Option<f64>),
) -> (HashSet<u32>, HashSet<u32>) {
    let (mut left_out, mut at_genus) = (HashSet::new(), HashSet::new());
    let mut congeners: HashMap<u32, Vec<u32>> = HashMap::new();
    for taxid in species.keys() {
        if let Some(&genus) = genus_of.get(taxid) {
            congeners.entry(genus).or_default().push(*taxid);
        }
    }
    if let Some(margin) = margin.filter(|_| species.len() >= MIN_JUDGED_SPECIES) {
        let low = |taxid: u32| species[&taxid].0 < baseline - margin;
        for &taxid in species.keys().filter(|&&taxid| low(taxid)) {
            match genus_of.get(&taxid) {
                Some(genus) if congeners[genus].iter().any(|&c| !low(c)) => left_out.insert(taxid),
                Some(_) => at_genus.insert(taxid),
                None => false,
            };
        }
    }
    if let Some(gap) = gap {
        for members in congeners.values() {
            for &taxid in members {
                let (identity, depth) = species[&taxid];
                if members
                    .iter()
                    .any(|c| species[c].1 > depth && species[c].0 >= identity + gap)
                {
                    at_genus.remove(&taxid);
                    left_out.insert(taxid);
                }
            }
        }
    }
    (left_out, at_genus)
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
                    let resolved =
                        (unclassified.is_some() && share >= RESOLVE_POSTERIOR).then(|| {
                            (
                                candidate.placements.clone(),
                                candidate.other_gene_hits.clone(),
                            )
                        });
                    Some((
                        taxid,
                        share
                            * marker_fragment_depth(
                                candidate,
                                reported[&taxid],
                                policy.mode,
                                model,
                            ),
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

/// The depth that a fragment adds on marker genome `assembly_id` of its
/// `candidate` species ([`fragment_depth`]), from its placement there or else
/// its best one, with the first mate's mapped bases for short reads.
fn marker_fragment_depth(
    candidate: &TerminalCandidate,
    assembly_id: u32,
    mode: EvidenceMode,
    model: &MarkerModel,
) -> f64 {
    let placement = candidate
        .placements
        .iter()
        .find(|p| p.assembly_id == assembly_id)
        .copied()
        .unwrap_or(candidate.best);
    let mapped = match mode {
        EvidenceMode::Paired | EvidenceMode::Single => placement.r1_mapped,
        EvidenceMode::Long => placement.both_mapped,
    };
    fragment_depth(placement, mapped, Some(model))
}

/// Build the candidates of each item's read on the profiler's threads, take
/// what a pass needs with `extract`, and hand the results to `fold` in item
/// order, [`PARALLEL_BLOCK`] items at a time.
fn fold_candidates<I: Sync, T: Send>(
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

/// Add the placements of the fragments `list` for species `taxid` to
/// `genomes`, in list order ([`add_hits`]).
fn add_fragment_hits(
    data: &MappingData,
    policy: EvidencePolicy,
    fragments: &[Fragment],
    taxid: u32,
    list: &[usize],
    genomes: &mut HashMap<u32, GenomeEvidence>,
) {
    fold_candidates(
        data,
        policy,
        list,
        |&f| fragments[f].read,
        |_, candidates| species_placements(candidates, taxid),
        |_, hits| add_hits(genomes, taxid, hits),
    );
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
        genome.kept = self.covers(&coverage) && wins_enough(won, candidates);
        genome.coverage = Some(coverage);
    }

    /// Whether a genome's coverage is enough and evenly spread.
    fn covers(&self, coverage: &EvidenceCoverage) -> bool {
        self.amount_ok(coverage) && self.marker_genes_ok(coverage) && self.evenness_ok(coverage)
    }

    /// Whether a marker genome has enough genes hit (any other genome has).
    fn marker_genes_ok(&self, coverage: &EvidenceCoverage) -> bool {
        coverage
            .marker_genes_hit
            .is_none_or(|hit| hit >= self.min_marker_genes)
    }

    /// Whether a genome's chunk breadth is above the cutoff.
    fn amount_ok(&self, coverage: &EvidenceCoverage) -> bool {
        coverage.fixed_chunk > self.chunk_cutoff
    }

    /// Whether both breadth ratios are within their bounds.
    fn evenness_ok(&self, coverage: &EvidenceCoverage) -> bool {
        within(coverage.oebr, self.ratio_bounds.0) && within(coverage.coebr, self.ratio_bounds.1)
    }
}

/// Whether a species won enough (`won`) of the `candidates` fragments listing
/// it: at least [`MIN_WON_SHARE`].
fn wins_enough(won: u64, candidates: u64) -> bool {
    won as f64 >= MIN_WON_SHARE * candidates.max(1) as f64
}

/// The first EM over every candidate species and the species each read group
/// resolves to.
struct FirstEm<'a> {
    fragments: &'a [Fragment],
    group_of: &'a [Option<usize>],
    em: &'a likelihood_em::EmResult,
    resolved: &'a [Option<u32>],
}

impl FirstEm<'_> {
    /// Fragments listing each species among their candidates, those resolved
    /// to any species, and those resolved to it, counted per read group.
    fn counts(&self) -> FragmentCounts {
        let mut counts = FragmentCounts::default();
        let sizes = group_sizes(self.group_of, self.em.posteriors.len());
        for ((posteriors, &resolved), size) in
            self.em.posteriors.iter().zip(self.resolved).zip(sizes)
        {
            for &(taxid, _) in posteriors {
                *counts.candidates.entry(taxid).or_default() += size;
                if resolved.is_some() {
                    *counts.resolved_candidates.entry(taxid).or_default() += size;
                }
            }
            if let Some(taxid) = resolved {
                *counts.won.entry(taxid).or_default() += size;
            }
        }
        counts
    }

    /// Each fragment resolved to a species, with that species, in fragment
    /// order.
    fn won(&self) -> impl Iterator<Item = (usize, u32)> + '_ {
        self.group_of
            .iter()
            .enumerate()
            .filter_map(|(fragment, group)| {
                group
                    .and_then(|g| self.resolved[g])
                    .map(|taxid| (fragment, taxid))
            })
    }

    /// `taxids` ordered by their weight in the first EM, largest first, then
    /// by taxid.
    fn by_weight(&self, taxids: impl Iterator<Item = u32>) -> Vec<u32> {
        let weight = |taxid: &u32| self.em.weight.get(taxid).copied().unwrap_or(0.0);
        let mut order = taxids.collect::<Vec<_>>();
        order.sort_by(|a, b| weight(b).total_cmp(&weight(a)).then(a.cmp(b)));
        order
    }
}

/// The species with a kept genome.
fn kept_species(genomes: &HashMap<u32, GenomeEvidence>) -> HashSet<u32> {
    genomes
        .values()
        .filter(|g| g.kept)
        .map(|g| g.taxid)
        .collect()
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
/// fewer rivals only rises. Named units first: a removed unit labelled above
/// species (`labels`) is judged after the removed units of its taxon
/// (`ancestors` gives a unit's lineage), which do not lose to it the groups
/// they share with it ([`named_units_first`]). Returns the fragments each
/// species it kept took.
fn second_look(
    data: &MappingData,
    policy: EvidencePolicy,
    first: &FirstEm,
    tests: &SpeciesTests,
    counts: &mut FragmentCounts,
    genomes: &mut HashMap<u32, GenomeEvidence>,
    (labels, ancestors): (&HashMap<u32, ProfileRank>, impl Fn(u32) -> Vec<u32>),
) -> HashMap<u32, Vec<usize>> {
    let mut kept = kept_species(genomes);
    let order = first.by_weight(
        first
            .em
            .weight
            .iter()
            .filter(|&(taxid, &weight)| weight > 0.0 && !kept.contains(taxid))
            .map(|(&taxid, _)| taxid),
    );
    let removed: HashSet<u32> = order.iter().copied().collect();
    let labelled_above = if labels.is_empty() {
        HashMap::new()
    } else {
        order
            .iter()
            .filter_map(|&unit| {
                let above = ancestors(unit)
                    .into_iter()
                    .filter(|&taxon| {
                        taxon != unit && labels.contains_key(&taxon) && removed.contains(&taxon)
                    })
                    .collect::<HashSet<_>>();
                (!above.is_empty()).then_some((unit, above))
            })
            .collect()
    };
    let order = named_units_first(order, &labelled_above);
    let none = HashSet::new();
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
    let mut taken_by = HashMap::new();
    for taxid in order {
        let Some(groups) = open_groups.get(&taxid) else {
            continue;
        };
        let taken_fragments = groups
            .iter()
            .filter(|&&group| {
                takes_group(
                    &first.em.posteriors[group],
                    taxid,
                    &kept,
                    labelled_above.get(&taxid).unwrap_or(&none),
                )
            })
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
        add_fragment_hits(
            data,
            policy,
            first.fragments,
            taxid,
            &taken_fragments,
            &mut trial,
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
            taken_by.insert(taxid, taken_fragments);
        }
    }
    taken_by
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
    let mut kept = kept_species(genomes);
    let judged: HashSet<u32> = genomes
        .values()
        .filter(|g| !kept.contains(&g.taxid))
        .filter(|g| g.coverage.as_ref().is_some_and(|c| tests.covers(c)))
        .map(|g| g.taxid)
        .collect();
    // The fragments each such species won, and those listing it.
    let mut won: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, taxid) in first.won() {
        if judged.contains(&taxid) {
            won.entry(taxid).or_default().push(index);
        }
    }
    let mut listing: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, fragment) in first.fragments.iter().enumerate() {
        for &(taxid, _) in &fragment.candidates {
            if judged.contains(&taxid) {
                listing.entry(taxid).or_default().push(index);
            }
        }
    }
    let order = first.by_weight(won.keys().copied());
    let spans_genome = |coverages: &[(u32, EvidenceCoverage)]| {
        coverages
            .iter()
            .any(|(_, c)| c.fixed_chunk > tests.chunk_cutoff)
    };
    let mut rescued = HashMap::new();
    for taxid in order {
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
                let best = first.fragments[f].best_tied();
                best.len() > 1 && best.iter().any(|c| c.0 == taxid)
            })
            .collect::<Vec<_>>();
        let at_genus = !tied.is_empty() && spans_genome(&coverage(taxid, &tied));
        kept.insert(taxid);
        rescued.insert(taxid, at_genus);
    }
    rescued
}

/// Kept viruses judged on their own reads. A virus integrated in a host genome
/// (a prophage) collects the reads of that region of a host strain in the
/// sample whose copy matches the virus better than the copies in the indexed
/// genomes of its species, which align the same reads a little worse. The
/// hosts are the kept bacteria and archaea (`domain` gives a unit's domain),
/// with the taxa labelled above species of their lineage ([`with_labelled`]);
/// each kept virus is judged on its own reads beside them
/// ([`own_reads_look`]). A free virus whose host carries a related prophage
/// covers the rest of its genome evenly; reads of a prophage leave it bare or
/// patchy. Returns the species removed, in taxid order.
fn prophage_look(
    first: &FirstEm,
    taken: &HashMap<u32, Vec<usize>>,
    tests: &SpeciesTests,
    genomes: &mut HashMap<u32, GenomeEvidence>,
    (domain, labelled): (
        impl Fn(u32) -> Option<u32>,
        (&HashMap<u32, ProfileRank>, impl Fn(u32) -> Vec<u32>),
    ),
    evidence: impl FnMut(u32, &[usize]) -> HashMap<u32, GenomeEvidence>,
) -> Vec<u32> {
    let (mut viruses, mut hosts) = (HashSet::new(), HashSet::new());
    for genome in genomes.values().filter(|g| g.kept) {
        match domain(genome.taxid) {
            Some(VIRUSES) => viruses.insert(genome.taxid),
            Some(BACTERIA | ARCHAEA) => hosts.insert(genome.taxid),
            _ => false,
        };
    }
    if viruses.is_empty() || hosts.is_empty() {
        return Vec::new();
    }
    let hosts = with_labelled(hosts, labelled);
    let explainers = viruses
        .into_iter()
        .map(|virus| (virus, hosts.clone()))
        .collect::<HashMap<_, _>>();
    own_reads_look(first, taken, tests, genomes, &explainers, evidence)
        .into_iter()
        .map(|(taxid, _)| taxid)
        .collect()
}

/// `units` with the taxa of their lineage (`ancestors`) that references are
/// labelled with above species (`labels`), kept or not: the unnamed genomes
/// of a taxon stand for its kept named units, and carry their prophages.
fn with_labelled(
    mut units: HashSet<u32>,
    (labels, ancestors): (&HashMap<u32, ProfileRank>, impl Fn(u32) -> Vec<u32>),
) -> HashSet<u32> {
    if labels.is_empty() {
        return units;
    }
    let taxa = units
        .iter()
        .flat_map(|&unit| ancestors(unit))
        .filter(|taxon| labels.contains_key(taxon))
        .collect::<Vec<_>>();
    units.extend(taxa);
    units
}

/// Named units first: a kept unit labelled above species (`labels`) is
/// judged on its own reads beside the kept, more specific units of its own
/// lineage (`ancestors` gives a unit's lineage): the named species inside a
/// labelled genus, or the species and genera inside a labelled family
/// ([`own_reads_look`]). A new species of the taxon, with only unnamed
/// genomes, covers the rest of its genome and stays; an unnamed copy of a
/// named species is removed, and its reads go to that species in the second
/// EM. Returns the units removed, each with the units listed by the fragments
/// it shares with them.
fn named_first_look(
    first: &FirstEm,
    taken: &HashMap<u32, Vec<usize>>,
    tests: &SpeciesTests,
    genomes: &mut HashMap<u32, GenomeEvidence>,
    (labels, ancestors): (&HashMap<u32, ProfileRank>, impl Fn(u32) -> Vec<u32>),
    evidence: impl FnMut(u32, &[usize]) -> HashMap<u32, GenomeEvidence>,
) -> HashMap<u32, Vec<u32>> {
    let kept = kept_species(genomes);
    if !kept.iter().any(|unit| labels.contains_key(unit)) {
        return HashMap::new();
    }
    let lineage: HashMap<u32, Vec<u32>> =
        kept.iter().map(|&unit| (unit, ancestors(unit))).collect();
    let explainers = kept
        .iter()
        .filter(|unit| labels.contains_key(unit))
        .filter_map(|&label| {
            let inside = kept
                .iter()
                .copied()
                .filter(|&unit| unit != label && lineage[&unit].contains(&label))
                .collect::<HashSet<_>>();
            (!inside.is_empty()).then_some((label, inside))
        })
        .collect::<HashMap<_, _>>();
    if explainers.is_empty() {
        return HashMap::new();
    }
    own_reads_look(first, taken, tests, genomes, &explainers, evidence)
        .into_iter()
        .collect()
}

/// Kept units judged on their own reads beside the kept units that may
/// explain them (`explainers`, per unit). The fragments that list a unit
/// together with one of its explainers, whichever of them won each, mark the
/// part of the unit's genomes that the explainers account for. The unit's own
/// fragments, those it holds (won in the first EM or took on the second look)
/// that list none of its explainers, must pass the coverage tests over the
/// rest of one of its genomes alone ([`without_bases`]). A unit that fails is
/// removed, its reads going to the kept units in the second EM; one that no
/// fragment lists together with an explainer is left as it is. Returns the
/// units removed, in taxid order, each with the explainers its shared
/// fragments list.
fn own_reads_look(
    first: &FirstEm,
    taken: &HashMap<u32, Vec<usize>>,
    tests: &SpeciesTests,
    genomes: &mut HashMap<u32, GenomeEvidence>,
    explainers: &HashMap<u32, HashSet<u32>>,
    mut evidence: impl FnMut(u32, &[usize]) -> HashMap<u32, GenomeEvidence>,
) -> Vec<(u32, Vec<u32>)> {
    let mut held: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, taxid) in first.won() {
        if explainers.contains_key(&taxid) {
            held.entry(taxid).or_default().push(index);
        }
    }
    for (taxid, list) in taken {
        if explainers.contains_key(taxid) {
            held.entry(*taxid).or_default().extend(list);
        }
    }
    let lists_any = |f: usize, units: &HashSet<u32>| {
        first.fragments[f]
            .candidates
            .iter()
            .any(|c| units.contains(&c.0))
    };
    let mut shared: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, fragment) in first.fragments.iter().enumerate() {
        for c in &fragment.candidates {
            if explainers
                .get(&c.0)
                .is_some_and(|units| lists_any(index, units))
            {
                shared.entry(c.0).or_default().push(index);
            }
        }
    }
    let mut shared = shared.into_iter().collect::<Vec<_>>();
    shared.sort_unstable_by_key(|entry| entry.0);
    let mut removed = Vec::new();
    for (taxid, shared) in shared {
        let units = &explainers[&taxid];
        let list = held.remove(&taxid).unwrap_or_default();
        let own = list
            .iter()
            .copied()
            .filter(|&f| !lists_any(f, units))
            .collect::<Vec<_>>();
        let covered = !own.is_empty() && {
            let explained = evidence(taxid, &shared);
            evidence(taxid, &own).into_iter().any(|(id, genome)| {
                let mut rest = match explained.get(&id) {
                    Some(other) => without_bases(&genome.assembly, &other.assembly),
                    None => genome.assembly,
                };
                tests.covers(&tests.measure(id, &mut rest))
            })
        };
        if covered {
            continue;
        }
        for genome in genomes.values_mut() {
            if genome.taxid == taxid {
                genome.kept = false;
            }
        }
        let mut listed = shared
            .iter()
            .flat_map(|&f| first.fragments[f].candidates.iter().map(|c| c.0))
            .filter(|unit| units.contains(unit))
            .collect::<Vec<_>>();
        listed.sort_unstable();
        listed.dedup();
        removed.push((taxid, listed));
    }
    removed
}

/// `assembly` without the bases that the reads of `other` cover: each position
/// moves back by the covered bases before it on its contig and the genome
/// shrinks by all of them, so that the coverage tests measure `assembly` over
/// the remaining bases alone. A window within covered bases shrinks to
/// nothing, and a fragment whose reads all lie within them no longer counts.
fn without_bases(assembly: &PresenceAssembly, other: &PresenceAssembly) -> PresenceAssembly {
    let mut reads = other
        .intervals
        .iter()
        .flat_map(|r| [(r.reference_id, r.span), (r.reference_id, r.mate_span)])
        .filter(|(_, span)| span.1 > span.0)
        .collect::<Vec<_>>();
    reads.sort_unstable();
    // Covered segments, merged per contig, each with the covered bases before
    // it on its contig.
    let mut segments: Vec<(u32, (u32, u32), u32)> = Vec::new();
    for (contig, span) in reads {
        match segments.last_mut() {
            Some(last) if last.0 == contig && span.0 <= last.1 .1 => {
                last.1 .1 = last.1 .1.max(span.1);
            }
            last => {
                let before = match last {
                    Some(&mut (c, (start, end), before)) if c == contig => before + (end - start),
                    _ => 0,
                };
                segments.push((contig, span, before));
            }
        }
    }
    let covered: u32 = segments.iter().map(|s| s.1 .1 - s.1 .0).sum();
    let shift = |contig: u32, position: u32| {
        let after = segments.partition_point(|s| (s.0, s.1 .0) <= (contig, position));
        match after.checked_sub(1).map(|i| segments[i]) {
            Some((c, (start, end), before)) if c == contig => {
                position - before - (position.min(end) - start)
            }
            _ => position,
        }
    };
    let window = |contig: u32, span: (u32, u32)| (shift(contig, span.0), shift(contig, span.1));
    let mut rest = PresenceAssembly {
        asm_len: (assembly.asm_len - f64::from(covered)).max(0.0),
        sampling_scale: assembly.sampling_scale,
        ..PresenceAssembly::default()
    };
    for record in &assembly.intervals {
        let contig = record.reference_id;
        let span = window(contig, record.span);
        let mate_span = if record.mate_span.1 > record.mate_span.0 {
            window(contig, record.mate_span)
        } else {
            record.mate_span
        };
        if span.1 <= span.0 && mate_span.1 <= mate_span.0 {
            continue;
        }
        rest.count += 1;
        rest.reference_span_bases += u64::from(span.1 - span.0);
        rest.intervals.push(PresenceInterval {
            reference_id: contig,
            span,
            mate_span,
            fixed_chunk: window(contig, record.fixed_chunk),
            flex_chunk: window(contig, record.flex_chunk),
        });
    }
    rest
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
    evidence_of(data, policy, fragments, taxid, list)
        .into_iter()
        .map(|(id, mut genome)| (id, tests.measure(id, &mut genome.assembly)))
        .collect()
}

/// The genomes of species `taxid` with the placements of the fragments `list`
/// alone.
fn evidence_of(
    data: &MappingData,
    policy: EvidencePolicy,
    fragments: &[Fragment],
    taxid: u32,
    list: &[usize],
) -> HashMap<u32, GenomeEvidence> {
    let mut genomes = HashMap::new();
    add_fragment_hits(data, policy, fragments, taxid, list, &mut genomes);
    genomes
}

/// `order` with each removed unit labelled above species moved after the
/// last of the removed units whose lineage lists it (`labelled_above`, per
/// unit): named units are judged first.
fn named_units_first(order: Vec<u32>, labelled_above: &HashMap<u32, HashSet<u32>>) -> Vec<u32> {
    if labelled_above.is_empty() {
        return order;
    }
    let mut after: HashMap<u32, usize> = HashMap::new();
    for (position, unit) in order.iter().enumerate() {
        for &label in labelled_above.get(unit).into_iter().flatten() {
            after.insert(label, position);
        }
    }
    let mut keyed = order
        .into_iter()
        .enumerate()
        .map(|(position, unit)| {
            let target = after
                .get(&unit)
                .map_or(position, |&last| last.max(position));
            ((target, target > position, position), unit)
        })
        .collect::<Vec<_>>();
    keyed.sort_unstable_by_key(|&(key, _)| key);
    keyed.into_iter().map(|(_, unit)| unit).collect()
}

/// Whether removed species `taxid` takes a read group with first-EM
/// `posteriors` on the second look: no other removed species is more likely
/// there, and it holds at least the resolving posterior against the `kept`
/// species; the removed units labelled above species of its taxon (`above`)
/// are not its rivals.
fn takes_group(
    posteriors: &[(u32, f64)],
    taxid: u32,
    kept: &HashSet<u32>,
    above: &HashSet<u32>,
) -> bool {
    let own = posteriors
        .iter()
        .find(|p| p.0 == taxid)
        .map_or(0.0, |p| p.1);
    let mut with_kept = own;
    for &(other, posterior) in posteriors {
        if other == taxid || above.contains(&other) {
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

/// The number of fragments in each of `groups` read groups.
fn group_sizes(group_of: &[Option<usize>], groups: usize) -> Vec<u64> {
    let mut sizes = vec![0; groups];
    for &group in group_of.iter().flatten() {
        sizes[group] += 1;
    }
    sizes
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

/// The rows of one taxon: a unit's row (a species row, or the row of a taxon
/// its references are labelled with; none when every read is in a genome
/// row) with its genome rows, or a genus or family row; with the taxon's
/// breadth (a unit's on the genome it is reported on, which its row shows).
struct TaxonBlock {
    taxid: u32,
    taxon: Option<ProfileRow>,
    genomes: Vec<ProfileRow>,
    breadth: f64,
    /// The block of a kept unit.
    unit: bool,
}

impl TaxonBlock {
    /// Count a genus or family `row` of the unit's own taxon in the unit's row:
    /// reads and depth add, while coverage, best-fit genome and identity stay
    /// the unit's. A unit without a row of its own takes `row`.
    fn count(&mut self, row: ProfileRow) {
        match &mut self.taxon {
            Some(own) => {
                own.equivalent_mapcount += row.equivalent_mapcount;
                own.depth += row.depth;
            }
            None => self.taxon = Some(row),
        }
    }

    fn depth(&self) -> f64 {
        self.taxon.as_ref().map_or(0.0, |s| s.depth)
            + self.genomes.iter().map(|s| s.depth).sum::<f64>()
    }

    /// Full-index-equivalent reads of the taxon's rows.
    fn reads(&self) -> f64 {
        self.taxon
            .iter()
            .chain(&self.genomes)
            .fold(0.0, |total, row| total + row.equivalent_mapcount)
    }
}

/// `--min-reads` and `--min-breadth`, as the default engine applies them to
/// its rows: a taxon is dropped from the profile when its rows hold fewer
/// full-index-equivalent reads (a species counting with its genome rows), or
/// when its breadth is lower.
struct Minimums {
    reads: Option<usize>,
    breadth: Option<f64>,
}

/// What a dropped taxon missed: its reads and the minimum, its breadth and the
/// minimum.
#[derive(Clone, Copy, Debug, PartialEq)]
struct BelowMinimum {
    reads: Option<(f64, usize)>,
    breadth: Option<(f64, f64)>,
}

impl Minimums {
    fn check(&self, reads: f64, breadth: f64) -> Option<BelowMinimum> {
        let below = BelowMinimum {
            reads: self
                .reads
                .filter(|_| !equivalent_read_count_passes_minimum(reads, self.reads))
                .map(|minimum| (reads, minimum)),
            breadth: self
                .breadth
                .filter(|&minimum| breadth.is_nan() || breadth < minimum)
                .map(|minimum| (breadth, minimum)),
        };
        (below.reads.is_some() || below.breadth.is_some()).then_some(below)
    }
}

impl BelowMinimum {
    /// The minimums missed, in words.
    fn reason(&self) -> String {
        let mut missed = Vec::new();
        if let Some((reads, minimum)) = self.reads {
            missed.push(format!(
                "removed by --min-reads: {reads:.2} reads < {minimum}"
            ));
        }
        if let Some((breadth, minimum)) = self.breadth {
            missed.push(format!(
                "removed by --min-breadth: breadth {breadth:.4} < {minimum:.4}"
            ));
        }
        missed.join("; ")
    }
}

/// Set every row's abundance: its depth over the summed depth of the rows, in
/// percent.
fn set_abundance(blocks: &mut [TaxonBlock]) {
    let total_depth: f64 = blocks.iter().map(TaxonBlock::depth).sum();
    let abundance = |depth: f64| {
        if total_depth > 0.0 {
            100.0 * depth / total_depth
        } else {
            0.0
        }
    };
    for block in blocks {
        for row in block.taxon.iter_mut().chain(block.genomes.iter_mut()) {
            row.abundance = abundance(row.depth);
        }
    }
}

/// Why rows that passed their tests were dropped from the profile.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Removal {
    /// `--min-reads` or `--min-breadth`.
    BelowMinimum(BelowMinimum),
    /// `--host`: the pathogen table lists the taxon for other hosts only.
    NotPathogen(u32),
}

impl Removal {
    /// Why, in words.
    fn reason(&self) -> String {
        match self {
            Self::BelowMinimum(below) => below.reason(),
            Self::NotPathogen(host) => format!("removed by --host: not a pathogen of host {host}"),
        }
    }
}

/// The taxa whose rows were all dropped, and the genome (strain) rows dropped,
/// by genome.
#[derive(Default)]
struct Dropped {
    taxa: HashMap<u32, Removal>,
    genomes: HashMap<Arc<str>, Removal>,
}

/// Drop the blocks of the taxa below `minimums`, before abundance is
/// normalised; their reads go to no other taxon. Returns what each dropped
/// taxon missed.
fn drop_below_minimums(
    blocks: &mut Vec<TaxonBlock>,
    minimums: &Minimums,
) -> HashMap<u32, BelowMinimum> {
    let mut dropped = HashMap::new();
    if minimums.reads.is_none() && minimums.breadth.is_none() {
        return dropped;
    }
    blocks.retain(|block| match minimums.check(block.reads(), block.breadth) {
        Some(below) => {
            dropped.insert(block.taxid, below);
            false
        }
        None => true,
    });
    dropped
}

/// Each measured species' decision for the genus level: its strongest genome
/// in the default engine's witness order (lowest assembly id on a tie), with
/// the window spread of its reads, and that genome.
fn species_decisions(
    genomes: &HashMap<u32, GenomeEvidence>,
    rollup: &GenusRollup,
) -> HashMap<u32, (WitnessDecision, u32)> {
    let mut decisions: HashMap<u32, (WitnessDecision, u32)> = HashMap::new();
    for (&assembly_id, genome) in genomes {
        let Some(coverage) = genome.coverage else {
            continue;
        };
        let decision = WitnessDecision {
            metrics: EvidenceCoverage {
                window_spread: rollup
                    .layout
                    .window_spread(assembly_id, &genome.assembly.intervals),
                ..coverage
            },
            standard: genome.kept,
            rescues: RescueKinds::default(),
            required_premap_guard: PremapChunkGuard::default(),
        };
        let (id, decision) = match decisions.get(&genome.taxid) {
            Some(&(current, current_id)) => {
                super::coherent_witness_max((current_id, current), (assembly_id, decision))
            }
            None => (assembly_id, decision),
        };
        decisions.insert(genome.taxid, (decision, id));
    }
    decisions
}

/// Where fragments without a species-level candidate go: the build `policy`
/// (the genus fallback identity, mapped length and aligned fraction), the
/// unit each candidate species counts for (`unit_of`), the units that get no
/// row (`excluded`), the kept species whose best-matched reads stay with them
/// (`kept`; noisy long reads only, else empty) and the `rank` of its rows,
/// with the EM penalty, the species tests and the genome layout; with
/// `keep_shares`, the EM's shares of each fragment are kept.
struct UnmatchedLevel<'a> {
    policy: EvidencePolicy,
    penalty: f64,
    tests: &'a SpeciesTests,
    layout: &'a GenomeLayout,
    unit_of: &'a HashMap<u32, u32>,
    excluded: &'a HashSet<u32>,
    kept: &'a HashSet<u32>,
    rank: &'static str,
    keep_shares: bool,
}

/// Fragments that an EM shares among units: each fragment's read and read
/// group, and each group's posteriors.
struct UnitShares {
    reads: Vec<usize>,
    group_of: Vec<Option<usize>>,
    posteriors: Vec<Vec<(u32, f64)>>,
}

/// What the level for fragments without a species-level candidate finds: its
/// rows in unit order, every unit with resolved fragments as a candidate,
/// kept or not, and with `keep_shares` how its EM shares the fragments.
struct UnmatchedResult {
    rows: Vec<(u32, ProfileRow)>,
    candidates: Vec<UnitCandidate>,
    shares: Option<UnitShares>,
}

/// A unit of the level for fragments without a species-level candidate, with
/// what its tests saw: the genome holding most of its resolved fragments and
/// their number there, the fragments listing the unit and those resolved to it,
/// their window spread on that genome (not measured for an excluded unit), and
/// its coverage when it got a row.
struct UnitCandidate {
    unit: u32,
    genome: u32,
    reads: u64,
    listing: f64,
    won: f64,
    spread: Option<f64>,
    excluded: bool,
    coverage: Option<EvidenceCoverage>,
}

/// Rows from the fragments `reads` that have no candidate at the species
/// identity. They are built again with the level's policy; a fragment whose
/// best match, or one tied with it, is a kept species is that species' read,
/// too noisy for the species gate, and is left out ([`best_match_kept`]). The
/// others come from genomes missing from the database: their candidate
/// species are counted for their units; an EM shares them, and a fragment
/// resolves to a unit at [`RESOLVE_POSTERIOR`]. A unit gets a row when its
/// resolved fragments number at least [`MIN_GENUS_FRAGMENTS`] (as many as a
/// genus row of the roll-up needs), reach [`MIN_WINDOW_SPREAD`] of the windows
/// of the genome that holds most of them (the reads of a genome missing from
/// the database spread over the genome of a relative; reads of conserved genes
/// from elsewhere stay in a few places), and are at least [`MIN_WON_SHARE`] of
/// the fragments listing the unit.
fn unmatched_rows(
    profiler: &CommunityProfiler,
    data: &MappingData,
    reads: &[usize],
    level: &UnmatchedLevel,
) -> UnmatchedResult {
    let unit_of = |taxid: u32| level.unit_of.get(&taxid).copied();
    let policy = level.policy;
    let mut fragments = Vec::new();
    fold_candidates(
        data,
        policy,
        reads,
        |&read| read,
        |&read, candidates| {
            if !level.kept.is_empty() && best_match_kept(read, candidates, policy, level.kept) {
                return None;
            }
            fragment_of(read, candidates, policy, unit_of)
        },
        |_, found| fragments.extend(found.map(|(fragment, _)| fragment)),
    );
    let shared = share_reads(&fragments, |_| true, level.penalty);
    let em = likelihood_em::run(
        &shared.groups,
        level.penalty,
        profiler.cfg.max_em_iterations,
        EM_TOLERANCE,
    );
    // Fragments listing each unit and those resolved to it, counted per read
    // group, the resolved fragments, and each unit's share of the fragments'
    // bases.
    let resolved_units = em
        .posteriors
        .iter()
        .map(|posteriors| resolved_unit(posteriors))
        .collect::<Vec<_>>();
    let mut listing: HashMap<u32, f64> = HashMap::new();
    let mut won: HashMap<u32, f64> = HashMap::new();
    let sizes = group_sizes(&shared.group_of, em.posteriors.len());
    for ((posteriors, &resolved), size) in em.posteriors.iter().zip(&resolved_units).zip(sizes) {
        let size = size as f64;
        for &(unit, _) in posteriors {
            *listing.entry(unit).or_default() += size;
        }
        if let Some(unit) = resolved {
            *won.entry(unit).or_default() += size;
        }
    }
    let resolved = fragments
        .iter()
        .zip(&shared.group_of)
        .filter_map(|(fragment, group)| {
            group
                .and_then(|g| resolved_units[g])
                .map(|unit| (fragment.read, unit))
        })
        .collect::<Vec<_>>();
    let bases = shared.unit_bases(&em);
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
                .filter(|c| unit_of(c.taxid) == Some(unit))
                .min_by(|a, b| {
                    candidate_errors(&a.best, longest)
                        .total_cmp(&candidate_errors(&b.best, longest))
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
    let mut candidates = Vec::new();
    for (unit, id) in best {
        let genome = genomes.get_mut(&id).expect("measured above");
        let listed = listing.get(&unit).copied().unwrap_or(0.0);
        let won = won.get(&unit).copied().unwrap_or(0.0);
        let mut candidate = UnitCandidate {
            unit,
            genome: id,
            reads: genome.assembly.count,
            listing: listed,
            won,
            spread: None,
            excluded: level.excluded.contains(&unit),
            coverage: None,
        };
        if candidate.excluded {
            candidates.push(candidate);
            continue;
        }
        let spread = level.layout.window_spread(id, &genome.assembly.intervals);
        candidate.spread = spread;
        if won < MIN_GENUS_FRAGMENTS
            || won < MIN_WON_SHARE * listed
            || spread.is_none_or(|s| s < MIN_WINDOW_SPREAD)
        {
            candidates.push(candidate);
            continue;
        }
        let coverage = level.tests.measure(id, &mut genome.assembly);
        candidate.coverage = Some(coverage);
        candidates.push(candidate);
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
    candidates.sort_by_key(|candidate| candidate.unit);
    UnmatchedResult {
        rows,
        candidates,
        shares: level.keep_shares.then(|| UnitShares {
            reads: fragments.iter().map(|fragment| fragment.read).collect(),
            group_of: shared.group_of,
            posteriors: em.posteriors,
        }),
    }
}

/// Whether the best match of a fragment's `candidates` under `policy` (fewest
/// errors, with the species-level tie margins), or one tied with it, is a
/// species in `kept`.
fn best_match_kept(
    read: usize,
    candidates: &[TerminalCandidate],
    policy: EvidencePolicy,
    kept: &HashSet<u32>,
) -> bool {
    fragment_of(read, candidates, policy, Some).is_some_and(|(fragment, _)| {
        fragment
            .best_tied()
            .iter()
            .any(|(taxid, _)| kept.contains(taxid))
    })
}

/// What the fragments credited to kept species add to each: reads, bases,
/// and depth for a species reported on a marker genome; and when kept, each
/// credited fragment's read with its shares.
#[derive(Default)]
struct Credited {
    reads: HashMap<u32, f64>,
    bases: HashMap<u32, f64>,
    marker_depth: HashMap<u32, f64>,
    fragments: Vec<(usize, Vec<(u32, f64)>)>,
}

/// The fragments `reads` without a species-level candidate whose best match
/// under `policy` (the genus gates), or one tied with it, is a kept species
/// (with an EM `weight` and a `reported` genome): each counts for the tied kept
/// species, shared by their weights, as the final EM would share it. Its bases
/// add to a species' depth, or on a marker genome the depth it adds there
/// ([`fragment_depth`]). With `keep_shares`, each fragment's shares are kept.
fn credited_reads(
    data: &MappingData,
    policy: EvidencePolicy,
    reads: &[usize],
    (weight, reported): (&HashMap<u32, f64>, &HashMap<u32, u32>),
    keep_shares: bool,
) -> Credited {
    let markers = data.marker_model.as_deref();
    let mut credited = Credited::default();
    fold_candidates(
        data,
        policy,
        reads,
        |&read| read,
        |&read, candidates| {
            let (fragment, _) = fragment_of(read, candidates, policy, Some)?;
            let tied = fragment
                .best_tied()
                .iter()
                .map(|&(taxid, _)| taxid)
                .filter(|taxid| {
                    reported.contains_key(taxid) && weight.get(taxid).is_some_and(|&w| w > 0.0)
                })
                .collect::<Vec<_>>();
            let total: f64 = tied.iter().map(|taxid| weight[taxid]).sum();
            (total > 0.0).then(|| {
                tied.iter()
                    .map(|&taxid| {
                        let share = weight[&taxid] / total;
                        let on_marker = markers
                            .filter(|model| model.is_marker_genome(reported[&taxid]))
                            .and_then(|model| {
                                let candidate = candidates.iter().find(|c| c.taxid == taxid)?;
                                Some(
                                    share
                                        * marker_fragment_depth(
                                            candidate,
                                            reported[&taxid],
                                            policy.mode,
                                            model,
                                        ),
                                )
                            });
                        (taxid, share, share * fragment.bases, on_marker)
                    })
                    .collect::<Vec<_>>()
            })
        },
        |&read, found| {
            let Some(found) = found else {
                return;
            };
            if keep_shares {
                credited
                    .fragments
                    .push((read, found.iter().map(|f| (f.0, f.1)).collect()));
            }
            for (taxid, share, bases, on_marker) in found {
                *credited.reads.entry(taxid).or_insert(0.0) += share;
                match on_marker {
                    Some(depth) => *credited.marker_depth.entry(taxid).or_insert(0.0) += depth,
                    None => *credited.bases.entry(taxid).or_insert(0.0) += bases,
                }
            }
        },
    );
    credited
}

/// NCBI taxids of the viruses, the bacteria and the archaea.
const VIRUSES: u32 = 10239;
const BACTERIA: u32 = 2;
const ARCHAEA: u32 = 2157;

/// Which of [`VIRUSES`], [`BACTERIA`] and [`ARCHAEA`] the lineage of `taxid`
/// reaches, if any.
fn domain(taxonomy: &Taxonomy, taxid: u32) -> Option<u32> {
    taxonomy
        .get_parents(taxid)
        .iter()
        .map(|&(id, _, _)| id)
        .find(|&id| matches!(id, VIRUSES | BACTERIA | ARCHAEA))
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
/// roll-up with all their candidate species and their best one's placement. A
/// fragment whose best candidate is a unit labelled above species (`labels`)
/// is not offered: a removed genus unit's reads do not make a row for that
/// genus.
fn orphan_genus_rows(
    data: &MappingData,
    fragments: &[Fragment],
    group_of: &[Option<usize>],
    policy: EvidencePolicy,
    (rollup, labels): (&GenusRollup, &HashMap<u32, ProfileRank>),
    (decisions, accepted): (&HashMap<u32, WitnessDecision>, &HashSet<u32>),
    keep_offered: bool,
) -> (Vec<GenusRow>, Vec<(usize, u32)>) {
    // The roll-up takes only fragments whose candidate species all lie in one
    // genus; the others need no candidates built.
    let orphaned = fragments
        .iter()
        .zip(group_of)
        .filter(|(fragment, group)| {
            group.is_none()
                && !labels.contains_key(&fragment.candidates[0].0)
                && rollup
                    .common_genus(fragment.candidates.iter().map(|c| c.0))
                    .is_some()
        })
        .map(|(fragment, _)| fragment)
        .collect::<Vec<_>>();
    let mut orphans = OrphanedFragments::new(rollup, false);
    let mut offered = Vec::new();
    fold_candidates(
        data,
        policy,
        &orphaned,
        |fragment| fragment.read,
        // Candidates are ordered by errors; the first is the best.
        |fragment, candidates| {
            candidates
                .iter()
                .find(|c| c.taxid == fragment.candidates[0].0)
                .map(|c| c.best)
        },
        |fragment, best| {
            let Some(best) = best else {
                return;
            };
            let taxids = fragment.candidates.iter().map(|c| c.0).collect::<Vec<_>>();
            if keep_offered {
                offered.push((fragment.read, taxids[0]));
            }
            orphans.offer(
                &taxids,
                0,
                f64::from(best.r1_mapped) / best.asm_len.max(1.0),
                f64::from(best.r1_mapped),
                best.sampling_scale,
                &data.read_alignments[fragment.read].read_name,
            );
        },
    );
    (orphans.into_genus_rows(decisions, accepted), offered)
}

/// One profile row; `amount` is `(reads, depth)`.
fn profile_row(
    name: String,
    taxid: u32,
    rank: &str,
    (reads, depth): (f64, f64),
    coverage: &EvidenceCoverage,
    genome: Arc<str>,
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
fn strain_taxids(cfg: &crate::profiler::ProfilerConfig) -> HashMap<Arc<str>, u32> {
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

/// Name and taxid of the genome row of genome `accession` of species `taxid`:
/// those of its strain when the genome declares one (named as the species when
/// the taxonomy has no name for it), otherwise the species' with the accession.
fn strain_label(
    taxonomy: &Taxonomy,
    strain_taxids: &HashMap<Arc<str>, u32>,
    taxid: u32,
    accession: &str,
) -> (String, u32) {
    let species_name = taxonomy.get_name(taxid).unwrap_or("unknown");
    match strain_taxids.get(accession) {
        Some(&strain) => (
            taxonomy
                .get_name(strain)
                .unwrap_or(species_name)
                .to_string(),
            strain,
        ),
        None => (format!("{species_name} {accession}"), taxid),
    }
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
        let genomes = HashMap::from([
            (10, evidence(1, true, 40)),
            (11, evidence(1, false, 80)),
            (12, evidence(2, false, 40)),
        ]);
        let decisions = species_decisions(&genomes, &rollup());
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
        let (rows, offered) = orphan_genus_rows(
            &data,
            &fragments,
            &group_of,
            policy,
            (&rollup, &HashMap::new()),
            (
                &HashMap::from([(1, decision(0.8)), (3, decision(0.8))]),
                &HashSet::new(),
            ),
            true,
        );
        // Every fragment of one genus is offered, with its best candidate as donor.
        assert_eq!(offered.len(), 250);
        assert!(offered
            .iter()
            .all(|&(read, donor)| read < 250 && donor == 1));
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].genus, 100);
        assert_eq!(rows[0].fragments, 250.0);
        // Read-1 bases over the genome's length, per fragment.
        assert!((rows[0].depth - 250.0 * 100.0 / 1_000_000.0).abs() < 1e-12);
        // A donor whose reads do not spread over its genome gives no row.
        let (rows, offered) = orphan_genus_rows(
            &data,
            &fragments,
            &group_of,
            policy,
            (&rollup, &HashMap::new()),
            (&HashMap::from([(1, decision(0.2))]), &HashSet::new()),
            false,
        );
        assert!(rows.is_empty() && offered.is_empty());
    }

    #[test]
    fn the_unclassified_estimate_counts_a_species_and_its_genome_rows_once() {
        // A species with 100 reads, 70 of them in two genome rows (one under a strain taxid, one under the
        // species'): its species row keeps the other 30, so its rows hold its 100 reads once.
        let (remainder, _) = species_remainder(
            100.0,
            15_000.0,
            [(50.0, 7_500.0), (20.0, 3_000.0)].into_iter(),
        )
        .unwrap();
        let row = |taxid: u32, rank: &str, reads: f64| {
            profile_row(
                rank.to_string(),
                taxid,
                rank,
                (reads, 1.0),
                &EvidenceCoverage::default(),
                Arc::from("GCF_1"),
                Some((f64::NAN, f64::NAN)),
            )
        };
        let rows = [
            row(1, "species", remainder),
            row(11, "strain", 50.0),
            row(1, "strain", 20.0),
            row(100, "genus", 40.0),
            row(1000, "family", 10.0),
        ];
        let estimate =
            crate::profiler::unclassified::estimate(&rows, &HashMap::new(), 200).unwrap();
        // 200 input fragments, 150 explained.
        assert!((estimate.fragments - 50.0).abs() < 1e-9);
        assert!((estimate.share - 0.25).abs() < 1e-12);
    }

    #[test]
    fn taxa_below_a_minimum_are_dropped_with_their_genome_rows() {
        let row = |rank: &str, reads: f64| {
            profile_row(
                rank.to_string(),
                1,
                rank,
                (reads, 1.0),
                &EvidenceCoverage::default(),
                Arc::from("GCF_1"),
                Some((f64::NAN, f64::NAN)),
            )
        };
        let blocks = || {
            vec![
                // A species with 120 reads in its species row and 30 in a genome row.
                TaxonBlock {
                    taxid: 1,
                    taxon: Some(row("species", 120.0)),
                    genomes: vec![row("strain", 30.0)],
                    breadth: 0.9,
                    unit: true,
                },
                TaxonBlock {
                    taxid: 100,
                    taxon: Some(row("genus", 10.0)),
                    genomes: Vec::new(),
                    breadth: 0.03,
                    unit: false,
                },
                TaxonBlock {
                    taxid: 1000,
                    taxon: Some(row("family", 300.0)),
                    genomes: Vec::new(),
                    breadth: f64::NAN,
                    unit: false,
                },
            ]
        };
        let taxids = |blocks: &[TaxonBlock]| blocks.iter().map(|b| b.taxid).collect::<Vec<_>>();
        // Without minimums nothing is dropped.
        let mut kept = blocks();
        let none = Minimums {
            reads: None,
            breadth: None,
        };
        assert!(drop_below_minimums(&mut kept, &none).is_empty());
        assert_eq!(taxids(&kept), [1, 100, 1000]);
        // A species counts with its genome rows: 150 reads pass 130.
        let mut kept = blocks();
        let reads = Minimums {
            reads: Some(130),
            breadth: None,
        };
        let dropped = drop_below_minimums(&mut kept, &reads);
        assert_eq!(taxids(&kept), [1, 1000]);
        assert_eq!(
            dropped[&100].reason(),
            "removed by --min-reads: 10.00 reads < 130"
        );
        // Breadth: below the minimum, or unknown.
        let mut kept = blocks();
        let breadth = Minimums {
            reads: Some(5),
            breadth: Some(0.05),
        };
        let dropped = drop_below_minimums(&mut kept, &breadth);
        assert_eq!(taxids(&kept), [1]);
        assert_eq!(
            dropped[&100].reason(),
            "removed by --min-breadth: breadth 0.0300 < 0.0500"
        );
        assert_eq!(
            dropped[&1000].reason(),
            "removed by --min-breadth: breadth NaN < 0.0500"
        );
        // Both minimums missed.
        let both = Minimums {
            reads: Some(20),
            breadth: Some(0.05),
        };
        assert_eq!(
            both.check(12.0, 0.031).unwrap().reason(),
            "removed by --min-reads: 12.00 reads < 20; removed by --min-breadth: breadth 0.0310 < 0.0500"
        );
        assert_eq!(both.check(20.0, 0.05), None);
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
                Arc::from("GCF_1"),
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
    fn a_removed_species_takes_the_fragments_it_ties_on_with_a_weaker_removed_sister() {
        let none = HashSet::new();
        // Tied with a removed sister that the first EM gives less weight.
        assert!(takes_group(&[(1, 0.82), (2, 0.18)], 1, &none, &none));
        // The weaker sister does not take them.
        assert!(!takes_group(&[(1, 0.82), (2, 0.18)], 2, &none, &none));
        // Nor fragments that another removed species fits better.
        assert!(!takes_group(&[(1, 0.05), (2, 0.95)], 1, &none, &none));
        // Against a kept species it needs the resolving posterior.
        let kept = HashSet::from([3]);
        assert!(!takes_group(&[(1, 0.82), (3, 0.18)], 1, &kept, &none));
        assert!(takes_group(
            &[(1, 0.82), (2, 0.1), (3, 0.08)],
            1,
            &kept,
            &none
        ));
        // The sister judged after a kept species keeps losing the ties.
        let after = HashSet::from([1]);
        assert!(!takes_group(&[(1, 0.82), (2, 0.18)], 2, &after, &none));
        // A removed unit labelled above species (100) of its taxon is no
        // rival; kept units still are.
        let above = HashSet::from([100]);
        assert!(!takes_group(&[(100, 0.7), (1, 0.3)], 1, &none, &none));
        assert!(takes_group(&[(100, 0.7), (1, 0.3)], 1, &none, &above));
        assert!(!takes_group(
            &[(100, 0.6), (1, 0.3), (3, 0.1)],
            1,
            &kept,
            &above
        ));
    }

    #[test]
    fn the_taxa_labelled_above_species_of_a_host_are_hosts_too() {
        let labels = HashMap::from([(100, ProfileRank::Genus), (200, ProfileRank::Family)]);
        let ancestors = |taxid: u32| match taxid {
            1 => vec![1, 100, 200, 2],
            3 => vec![3, 300],
            _ => Vec::new(),
        };
        let hosts = with_labelled(HashSet::from([1, 3]), (&labels, ancestors));
        assert_eq!(hosts, HashSet::from([1, 3, 100, 200]));
        // Without labelled references the hosts stay.
        let hosts = with_labelled(HashSet::from([1, 3]), (&HashMap::new(), ancestors));
        assert_eq!(hosts, HashSet::from([1, 3]));
    }

    #[test]
    fn named_units_are_judged_before_the_labelled_units_of_their_taxon() {
        // By weight: genus unit 100, unrelated 7, species 1 of genus 100,
        // family unit 200 above both, species 2 of genus 100.
        let order = vec![100, 7, 1, 200, 2];
        let above = HashMap::from([
            (1, HashSet::from([100, 200])),
            (2, HashSet::from([100, 200])),
            (100, HashSet::from([200])),
        ]);
        assert_eq!(
            named_units_first(order.clone(), &above),
            vec![7, 1, 2, 100, 200]
        );
        // Without labelled units the order stays.
        assert_eq!(named_units_first(order.clone(), &HashMap::new()), order);
    }

    #[test]
    fn group_members_list_each_group_s_fragments_in_order() {
        let members = GroupMembers::new(&[Some(1), None, Some(0), Some(1), Some(2)], 3);
        assert_eq!(members.of(0), &[2]);
        assert_eq!(members.of(1), &[0, 3]);
        assert_eq!(members.of(2), &[4]);
    }

    /// A fragment's placement on contig `contig` of a genome, its windows the
    /// read itself.
    fn read_at(contig: u32, span: (u32, u32)) -> PresenceInterval {
        PresenceInterval {
            reference_id: contig,
            span,
            mate_span: (0, 0),
            fixed_chunk: span,
            flex_chunk: span,
        }
    }

    #[test]
    fn bases_covered_by_other_reads_are_cut_out_of_the_genome() {
        // Other reads cover 100..300 of contig 0 and 0..50 of contig 1.
        let other = PresenceAssembly {
            intervals: vec![
                read_at(0, (100, 200)),
                read_at(0, (150, 300)),
                read_at(1, (0, 50)),
            ],
            ..PresenceAssembly::default()
        };
        let assembly = PresenceAssembly {
            asm_len: 1000.0,
            count: 4,
            reference_span_bases: 310,
            intervals: vec![
                read_at(0, (0, 50)),
                read_at(0, (250, 400)),
                read_at(0, (120, 180)),
                read_at(1, (40, 100)),
            ],
            ..PresenceAssembly::default()
        };
        let rest = without_bases(&assembly, &other);
        assert_eq!(rest.asm_len, 750.0);
        assert_eq!(rest.count, 3);
        assert_eq!(rest.reference_span_bases, 50 + 100 + 50);
        let spans = rest
            .intervals
            .iter()
            .map(|r| (r.reference_id, r.span))
            .collect::<Vec<_>>();
        assert_eq!(spans, [(0, (0, 50)), (0, (100, 200)), (1, (0, 50))]);
    }

    #[test]
    fn kept_viruses_must_cover_what_no_kept_cellular_species_explains() {
        // Bacterium 1 and eukaryote 2 are kept, as are viruses 7, 8 and 9.
        // Virus 7 won 10 fragments, 4 of them listing species 1, and took 2
        // more on the second look; virus 8 won 5 fragments that list species 2
        // only; every fragment of virus 9 lists species 1.
        let fragment = |candidates: &[(u32, f64)]| Fragment {
            read: 0,
            candidates: candidates.to_vec(),
            bases: 100.0,
        };
        let mut fragments = Vec::new();
        let mut resolved = Vec::new();
        for i in 0..10 {
            fragments.push(fragment(if i < 4 {
                &[(7, 0.0), (1, 1.0)]
            } else {
                &[(7, 0.0)]
            }));
            resolved.push(Some(7));
        }
        for _ in 0..5 {
            fragments.push(fragment(&[(8, 0.0), (2, 0.0)]));
            resolved.push(Some(8));
            fragments.push(fragment(&[(9, 0.0), (1, 0.0)]));
            resolved.push(Some(9));
        }
        let taken = HashMap::from([(7, vec![fragments.len(), fragments.len() + 1])]);
        for _ in 0..2 {
            fragments.push(fragment(&[(7, 0.0), (8, 1.0)]));
            resolved.push(None);
        }
        let group_of = (0..fragments.len()).map(Some).collect::<Vec<_>>();
        let em = likelihood_em::EmResult {
            weight: HashMap::new(),
            posteriors: Vec::new(),
            iterations: 0,
        };
        let first = FirstEm {
            fragments: &fragments,
            group_of: &group_of,
            em: &em,
            resolved: &resolved,
        };
        let tests = SpeciesTests {
            multiple_tests: 1,
            chunk_cutoff: 0.3,
            ratio_bounds: (0.5, 0.5),
            markers: None,
            min_marker_genes: 3,
        };
        let genome = |taxid, intervals: Vec<PresenceInterval>| GenomeEvidence {
            taxid,
            assembly: PresenceAssembly {
                asm_len: 1000.0,
                count: intervals.len() as u64,
                reference_span_bases: intervals
                    .iter()
                    .map(|r| u64::from(r.span.1 - r.span.0))
                    .sum(),
                intervals,
                ..PresenceAssembly::default()
            },
            coverage: None,
            kept: true,
        };
        let start = HashMap::from([
            (10, genome(1, Vec::new())),
            (20, genome(2, Vec::new())),
            (70, genome(7, Vec::new())),
            (80, genome(8, Vec::new())),
            (90, genome(9, Vec::new())),
        ]);
        // Virus 7's 4 shared fragments cover the first 600 bases of its
        // genome; its 8 own ones lie every 50 bases over the other 400, or all
        // in one place.
        let look = |spread: bool| {
            let mut genomes = start.clone();
            let removed = prophage_look(
                &first,
                &taken,
                &tests,
                &mut genomes,
                (
                    |taxid| match taxid {
                        1 => Some(BACTERIA),
                        2 => None,
                        _ => Some(VIRUSES),
                    },
                    (&HashMap::new(), |_| Vec::new()),
                ),
                |taxid, list| {
                    assert_eq!(taxid, 7);
                    let reads = match list.len() {
                        4 => (0..4)
                            .map(|i| read_at(0, (150 * i, 150 * i + 150)))
                            .collect(),
                        8 => (0..8)
                            .map(|i| {
                                let start = if spread { 600 + 50 * i } else { 600 };
                                read_at(0, (start, start + 25))
                            })
                            .collect(),
                        n => panic!("unexpected list of {n}"),
                    };
                    HashMap::from([(70, genome(7, reads))])
                },
            );
            let kept = [10, 20, 70, 80, 90].map(|id| genomes[&id].kept);
            (removed, kept)
        };
        // Own reads spread over the unshared 400 bases: kept, though over the
        // whole genome they would cover too little of it.
        assert_eq!(look(true), (vec![9], [true, true, true, true, false]));
        // Own reads in one place: removed with virus 9.
        assert_eq!(look(false), (vec![7, 9], [true, true, false, true, false]));
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
    fn references_labelled_above_species_are_candidate_units() {
        // A read on a genome of species 2 and, better, on a genome labelled only at genus 100.
        let data = data(vec![vec![
            row(2, 20, 0.97),
            crate::profiler::evidence::EvidenceAlignment {
                profile_rank: ProfileRank::Genus,
                ..row(100, 30, 0.99)
            },
        ]]);
        let policy = EvidencePolicy {
            mode: EvidenceMode::Single,
            identity: 0.95,
            mapped_len: 50,
            fraction: 0.6,
            collect_exact_ties: false,
            collect_identities: false,
            collect_best_placements: false,
        };
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&data.read_alignments[0].alignments, policy);
        let (fragment, _) = fragment_of(0, candidates, policy, Some).unwrap();
        let units = fragment.candidates.iter().map(|c| c.0).collect::<Vec<_>>();
        assert_eq!(units, vec![100, 2]);
    }

    #[test]
    fn a_labelled_unit_must_cover_what_the_named_units_of_its_taxon_do_not_explain() {
        // Genus unit 100 and its named species 1 are kept, as is species 2 of another genus. Unit 100 won 12
        // fragments, 4 of them listing species 1; unit 200 (labelled at another genus) won 4 listing species 2 only.
        let fragment = |candidates: &[(u32, f64)]| Fragment {
            read: 0,
            candidates: candidates.to_vec(),
            bases: 100.0,
        };
        let mut fragments = Vec::new();
        let mut resolved = Vec::new();
        for i in 0..12 {
            fragments.push(fragment(if i < 4 {
                &[(100, 0.0), (1, 1.0)]
            } else {
                &[(100, 0.0)]
            }));
            resolved.push(Some(100));
        }
        for _ in 0..4 {
            fragments.push(fragment(&[(200, 0.0), (2, 0.0)]));
            resolved.push(Some(200));
        }
        let group_of = (0..fragments.len()).map(Some).collect::<Vec<_>>();
        let em = likelihood_em::EmResult {
            weight: HashMap::new(),
            posteriors: Vec::new(),
            iterations: 0,
        };
        let first = FirstEm {
            fragments: &fragments,
            group_of: &group_of,
            em: &em,
            resolved: &resolved,
        };
        let tests = SpeciesTests {
            multiple_tests: 1,
            chunk_cutoff: 0.3,
            ratio_bounds: (0.5, 0.5),
            markers: None,
            min_marker_genes: 0,
        };
        let genome = |taxid, intervals: Vec<PresenceInterval>| GenomeEvidence {
            taxid,
            assembly: PresenceAssembly {
                asm_len: 1000.0,
                count: intervals.len() as u64,
                reference_span_bases: intervals
                    .iter()
                    .map(|r| u64::from(r.span.1 - r.span.0))
                    .sum(),
                intervals,
                ..PresenceAssembly::default()
            },
            coverage: None,
            kept: true,
        };
        let start = HashMap::from([
            (10, genome(1, Vec::new())),
            (20, genome(2, Vec::new())),
            (1000, genome(100, Vec::new())),
            (2000, genome(200, Vec::new())),
        ]);
        let labels = HashMap::from([(100, ProfileRank::Genus), (200, ProfileRank::Genus)]);
        // Species 1 lies in genus 100, species 2 in genus 300.
        let ancestors = |taxid: u32| match taxid {
            1 => vec![1, 100, 1000],
            2 => vec![2, 300, 1000],
            100 => vec![100, 1000],
            _ => vec![taxid],
        };
        // The 4 fragments shared with species 1 cover the first 600 bases of unit 100's genome; its 8 own ones lie
        // every 50 bases over the other 400 (a new species of the genus), or all in one place (an unnamed copy of
        // species 1 that collects a few reads of its own).
        let look = |spread: bool| {
            let mut genomes = start.clone();
            let removed = named_first_look(
                &first,
                &HashMap::new(),
                &tests,
                &mut genomes,
                (&labels, ancestors),
                |taxid, list| {
                    assert_eq!(taxid, 100);
                    let reads = match list.len() {
                        4 => (0..4)
                            .map(|i| read_at(0, (150 * i, 150 * i + 150)))
                            .collect(),
                        8 => (0..8)
                            .map(|i| {
                                let start = if spread { 600 + 50 * i } else { 600 };
                                read_at(0, (start, start + 25))
                            })
                            .collect(),
                        n => panic!("unexpected list of {n}"),
                    };
                    HashMap::from([(1000, genome(100, reads))])
                },
            );
            let kept = [10, 20, 1000, 2000].map(|id| genomes[&id].kept);
            (removed, kept)
        };
        // A new species of the genus stays; unit 200 has no named unit of its genus kept and is not judged.
        assert_eq!(look(true), (HashMap::new(), [true, true, true, true]));
        assert_eq!(
            look(false),
            (HashMap::from([(100, vec![1])]), [true, true, false, true])
        );
        // Without labelled units nothing is judged.
        let mut genomes = start.clone();
        let removed = named_first_look(
            &first,
            &HashMap::new(),
            &tests,
            &mut genomes,
            (&HashMap::new(), ancestors),
            |_, _| panic!("no unit to judge"),
        );
        assert!(removed.is_empty());
    }

    #[test]
    fn a_genus_row_of_a_labelled_unit_s_taxon_is_counted_in_the_unit_s_row() {
        let coverage = EvidenceCoverage {
            breadth: 0.8,
            ..EvidenceCoverage::default()
        };
        let row = |rank: &str, reads: f64, genome: &str, coverage: &EvidenceCoverage| {
            profile_row(
                "Genus".to_string(),
                100,
                rank,
                (reads, reads / 100.0),
                coverage,
                Arc::from(genome),
                Some((0.99, 0.995)),
            )
        };
        let rollup = row("genus", 40.0, "GCF_2", &EvidenceCoverage::default());
        let mut block = TaxonBlock {
            taxid: 100,
            taxon: Some(row("genus", 60.0, "GCF_1", &coverage)),
            genomes: Vec::new(),
            breadth: 0.8,
            unit: true,
        };
        block.count(rollup.clone());
        let merged = block.taxon.as_ref().unwrap();
        assert_eq!(merged.equivalent_mapcount, 100.0);
        assert!((merged.depth - 1.0).abs() < 1e-12);
        // The unit's coverage, genome and identity stay.
        assert_eq!(merged.breadth, 0.8);
        assert_eq!(&*merged.best_fit.as_ref().unwrap().genome, "GCF_1");
        // A unit whose reads all sit in genome rows takes the row.
        let mut block = TaxonBlock {
            taxid: 100,
            taxon: None,
            genomes: vec![row("strain", 60.0, "GCF_1", &coverage)],
            breadth: 0.8,
            unit: true,
        };
        block.count(rollup);
        assert_eq!(block.taxon.as_ref().unwrap().equivalent_mapcount, 40.0);
    }

    #[test]
    fn a_genome_of_a_labelled_unit_is_named_after_the_label() {
        let taxonomy = Taxonomy::from_parts(
            HashMap::from([
                (1, (1, "no rank".to_string())),
                (561, (1, "genus".to_string())),
                (562, (561, "species".to_string())),
                (83333, (562, "strain".to_string())),
            ]),
            HashMap::from([
                (561, "Escherichia".to_string()),
                (562, "Escherichia coli".to_string()),
                (83333, "Escherichia coli K-12".to_string()),
            ]),
            HashMap::new(),
            HashSet::new(),
        );
        let strain_taxids = HashMap::from([(Arc::from("GCF_K12"), 83333)]);
        assert_eq!(
            strain_label(&taxonomy, &strain_taxids, 561, "GCF_1"),
            ("Escherichia GCF_1".to_string(), 561)
        );
        assert_eq!(
            strain_label(&taxonomy, &strain_taxids, 562, "GCF_K12"),
            ("Escherichia coli K-12".to_string(), 83333)
        );
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

    #[test]
    fn far_less_identical_species_beside_a_close_congener_are_left_out() {
        // Genus 10: species 1 close, species 2 far below; genus 20: species 3
        // far below alone; genus 30: species 4 and 5 close; species 6 far below
        // without a genus. Depths are equal.
        let genus_of = HashMap::from([(1, 10), (2, 10), (3, 20), (4, 30), (5, 30)]);
        let species = HashMap::from([
            (1, (0.998, 1.0)),
            (2, (0.96, 1.0)),
            (3, (0.95, 1.0)),
            (4, (0.997, 1.0)),
            (5, (0.996, 1.0)),
            (6, (0.94, 1.0)),
        ]);
        let (left_out, at_genus) =
            identity_relatives(&genus_of, &species, 0.997, (Some(0.025), None));
        assert_eq!(left_out, HashSet::from([2]));
        assert_eq!(at_genus, HashSet::from([3]));
        // Too few species for a baseline.
        let few = HashMap::from([(1, (0.998, 1.0)), (2, (0.96, 1.0))]);
        let (left_out, at_genus) = identity_relatives(&genus_of, &few, 0.997, (Some(0.025), None));
        assert!(left_out.is_empty() && at_genus.is_empty());
    }

    #[test]
    fn a_species_is_left_out_beside_a_deeper_and_more_identical_congener() {
        // Genus 40: species 7 deep and identical; species 8 shallower and 0.015
        // less identical; species 9 deeper than 7 but less identical.
        let genus_of = HashMap::from([(7, 40), (8, 40), (9, 40)]);
        let species = HashMap::from([(7, (0.99, 1.0)), (8, (0.975, 0.5)), (9, (0.95, 2.0))]);
        let (left_out, at_genus) =
            identity_relatives(&genus_of, &species, 0.99, (None, Some(0.01)));
        assert_eq!(left_out, HashSet::from([8]));
        assert!(at_genus.is_empty());
        // A wider gap keeps it.
        let (left_out, _) = identity_relatives(&genus_of, &species, 0.99, (None, Some(0.02)));
        assert!(left_out.is_empty());
    }

    #[test]
    fn a_read_whose_best_match_is_a_kept_species_stays_with_it() {
        // One read below the species identity: best on species 3 (identity
        // 0.90), then species 2 (0.88) and 4 (0.85).
        let data = data(vec![vec![
            row(2, 20, 0.88),
            row(3, 30, 0.90),
            row(4, 40, 0.85),
        ]]);
        let policy = EvidencePolicy {
            mode: EvidenceMode::Single,
            identity: 0.80,
            mapped_len: 50,
            fraction: 0.6,
            collect_exact_ties: false,
            collect_identities: false,
            collect_best_placements: false,
        };
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&data.read_alignments[0].alignments, policy);
        assert!(best_match_kept(0, candidates, policy, &HashSet::from([3])));
        // Kept species that match worse than the best do not hold the read.
        assert!(!best_match_kept(
            0,
            candidates,
            policy,
            &HashSet::from([2, 4])
        ));
        assert!(!best_match_kept(0, candidates, policy, &HashSet::new()));
    }

    #[test]
    fn a_read_below_the_species_gates_counts_for_the_kept_species_it_matches_best() {
        // Read 0 ties on kept species 2 and 3 (weights 30 and 10); read 1
        // matches species 4 (not kept) best and species 2 worse.
        let data = data(vec![
            vec![row(2, 20, 0.88), row(3, 30, 0.88)],
            vec![row(4, 40, 0.90), row(2, 20, 0.85)],
        ]);
        let policy = EvidencePolicy {
            mode: EvidenceMode::Single,
            identity: 0.80,
            mapped_len: 50,
            fraction: 0.6,
            collect_exact_ties: false,
            collect_identities: false,
            collect_best_placements: false,
        };
        let weight = HashMap::from([(2, 30.0), (3, 10.0)]);
        let reported = HashMap::from([(2, 20), (3, 30)]);
        let credited = credited_reads(&data, policy, &[0, 1], (&weight, &reported), true);
        assert!((credited.reads[&2] - 0.75).abs() < 1e-12);
        assert!((credited.reads[&3] - 0.25).abs() < 1e-12);
        assert!((credited.bases[&2] - 75.0).abs() < 1e-9);
        assert!(!credited.reads.contains_key(&4));
    }
}

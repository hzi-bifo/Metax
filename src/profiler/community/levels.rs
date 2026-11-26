//! Level-wise profiling: species, then the levels below and above.
//!
//! A *fragment* is a read, or a read pair counted once. Its candidates are the
//! taxa whose genomes it aligns to within the read-level gates; candidates
//! whose estimated errors exceed the fragment's best by less than a margin
//! count as tied. A *unit* is a species, or one taxon above species holding
//! the genomes whose references are labelled only at that taxon; it is judged
//! and reported like a species at its rank, without the identity rules.
//!
//! 1. **First EM.** An EM weighted by how well each fragment matches each unit
//!    ([`likelihood_em`]) gives a fragment to a unit when its posterior there
//!    is at least 0.9; such a fragment supports every genome of that unit it
//!    aligns to.
//! 2. **Species tests.** A unit is kept when one of its genomes shows enough
//!    evenly spread coverage and the unit wins enough of the fragments whose
//!    read group in the first EM lists it.
//! 3. **Looks.** Removed units are judged again, strongest first, against the
//!    units kept so far, so that a removed sister no longer splits the
//!    fragments it ties on ([`second_look`]); named units are judged before
//!    the labelled unit of their taxon and do not lose to it the fragments
//!    they share with it ([`named_units_first`]). A removed unit is also kept
//!    when the fragments it won that list no kept unit cover one of its
//!    genomes on their own ([`own_evidence_look`]). A kept virus is removed
//!    when its own fragments do not cover the part of its genome that no kept
//!    bacterium or archaeon explains ([`prophage_look`]), and a kept labelled
//!    unit when they do not cover its genome beyond the part that the kept
//!    named units of its taxon explain ([`named_first_look`]); the reads of a
//!    removed unit go to the kept units in the second EM.
//! 4. **Even coverage.** A removed species whose fragments cover its strongest
//!    genome evenly, but too thinly for the amount test, is kept when no taxon
//!    of its genus is kept; not on a fractional database.
//! 5. **Second EM.** It shares the fragments among the kept units. A unit's
//!    abundance is the resulting depth on its best-fit genome. A fragment on a
//!    marker genome counts as one on the whole genome ([`fragment_depth`]),
//!    and such a genome is measured gene by gene with the same tests
//!    ([`MarkerModel`]).
//! 6. **Rows above species.** The fragments that no kept unit takes go to a
//!    genus when their matches lie in one and their evidence supports it,
//!    otherwise to the lowest common ancestor of their matches when enough of
//!    them spread over a genome ([`higher`]); not on a fractional database,
//!    where they are not reported above species. On noisy long reads, a
//!    fragment without a species-level candidate whose best match is a kept
//!    species counts for it instead ([`credited_reads`]).
//! 7. **Genomes.** The fragments of each kept species are shared among its
//!    genomes; the genome with the largest share is its best-fit genome, and
//!    with `--strain` a genome is also reported on its own ([`genomes`]).
//! 8. **Identity rules.** On a fractional database, a species whose reads are
//!    far less identical than the sample's typical species is left out beside
//!    a kept congener and otherwise reported at its genus; with
//!    `--relative-identity-gap` (on by default for single-end reads on such a
//!    database) a species is left out beside a deeper congener whose reads are
//!    that much more identical ([`identity_relatives`]).
//! 9. **Filters and abundance.** `--host` drops the rows of taxa that the
//!    pathogen table lists for other hosts only and annotates the others
//!    ([`host`]); `--min-reads` and `--min-breadth` drop a taxon whose rows
//!    hold fewer reads (a species with its genome rows) or whose breadth is
//!    lower ([`drop_below_minimums`]). Abundance is normalised over the rows
//!    left. Read counts are full-index equivalents: on a fractional database
//!    they are scaled up by each genome's indexed share.
//! 10. **Outputs.** The profile opens with two `#` lines: the run's metadata
//!     and the column names. `--save-decisions` writes
//!     `<prefix>.decisions.tsv`, how every taxon judged at every rank was
//!     decided, with the tests behind it ([`tables`]), and
//!     `--save-classification` writes `<prefix>.classify.txt.gz`, each
//!     fragment's shares of the profile's rows ([`classify`]).

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::sync::Arc;

use anyhow::Result;
use csv::WriterBuilder;
use rayon::prelude::*;

use super::MappingData;
use crate::profiler::evidence::{
    calculate_coverage, merge_spans, CandidatePlacement, CandidateScratch, EvidenceCoverage,
    EvidenceMode, EvidencePolicy, GeneHit, PresenceAssembly, PresenceInterval, TerminalCandidate,
};
use crate::profiler::likelihood_em::{self, ReadGroups};
use crate::profiler::markers::MarkerModel;
use crate::profiler::output::{BestFitReport, ProfileRow, HOST_COLUMNS, PROFILE_COLUMNS};
use crate::profiler::pathogen::load_pathogen_table;
use crate::profiler::rollup::GenomeLayout;
use crate::profiler::unclassified::{self, GenomeReads};
use crate::profiler::{
    equivalent_read_count_passes_minimum, ProfileRank, Profiler, ProfilerConfig,
};

mod classify;
mod genomes;
mod higher;
mod host;
mod tables;

use higher::Lineages;
use tables::{Decisions, FinalRows, Names, SpeciesDecision, StrainLevel, TableHeader};

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
/// least this share of the fragments whose read group in the first EM lists
/// it (candidates within [`MIN_RELATIVE_LIKELIHOOD`] of the best): a species that only collects reads which a present relative explains wins
/// few of them, while a present species wins most of its reads even when its
/// strain differs from every genome of it in the database.
const MIN_WON_SHARE: f64 = 0.3;
/// Lower bound of a species' base-level breadth ratio (OEBR), below its
/// chunk-level one: amplified libraries cover a genome unevenly base by base
/// while its chunks stay evenly covered, which the chunk-level ratio (COEBR)
/// still tests at its own bound.
const SPECIES_MIN_OEBR: f64 = 0.6;
/// Upper bound of both breadth ratios (OEBR and COEBR), at every level.
const MAX_BREADTH_RATIO: f64 = 1.50;
/// Lower bound of a species' chunk-level breadth ratio (COEBR) for short
/// reads.
const SPECIES_MIN_COEBR_SHORT_READS: f64 = 0.75;
/// Lower bound of a species' chunk-level breadth ratio (COEBR) for long reads.
const SPECIES_MIN_COEBR_LONG_READS: f64 = 0.70;
/// Lower bound of a genome's base-level breadth ratio (OEBR) at the genome
/// level (`--strain`).
const GENOME_MIN_OEBR: f64 = 0.70;
/// Lower bound of a genome's chunk-level breadth ratio (COEBR) at the genome
/// level, for short reads.
const GENOME_MIN_COEBR_SHORT_READS: f64 = 0.75;
/// Lower bound of a genome's chunk-level breadth ratio (COEBR) at the genome
/// level, for long reads.
const GENOME_MIN_COEBR_LONG_READS: f64 = 0.70;
/// Fewest species with best placements for an identity baseline to be taken.
const MIN_JUDGED_SPECIES: usize = 3;
/// Share of the read groups' weight that may still change when the EM stops.
const EM_TOLERANCE: f64 = 1e-7;

/// The sample's error rate is measured on this share of the fragments that
/// match their best candidate best.
const ERROR_RATE_SHARE: f64 = 0.75;
/// The candidate passes build the candidates of this many fragments at once on
/// the profiler's threads, then fold the results in fragment order, so that
/// interval lists and sums come out as in a single pass.
const PARALLEL_BLOCK: usize = 1 << 16;

/// The likelihood EM of a run: its penalty per error and its most iterations;
/// it stops at [`EM_TOLERANCE`].
#[derive(Clone, Copy)]
struct EmSettings {
    penalty: f64,
    max_iterations: usize,
}

impl EmSettings {
    fn run(&self, groups: &ReadGroups) -> likelihood_em::EmResult {
        likelihood_em::run(groups, self.penalty, self.max_iterations, EM_TOLERANCE)
    }
}

/// A fragment's candidate species with their errors (estimated errors of the
/// species' best placement plus the unaligned-base penalty), the read it comes
/// from, and its mapped bases (first mate of a pair) for depth.
struct Fragment {
    read: usize,
    candidates: Vec<(u32, f64)>,
    bases: f64,
}

impl Fragment {
    /// The candidates tied with the best ([`best_tied`]).
    fn best_tied(&self) -> &[(u32, f64)] {
        best_tied(&self.candidates)
    }
}

/// The candidates tied with the best among `ranked` candidates ([`ranked`]):
/// they are sorted by errors, and those within the tie margin of the best
/// carry exactly its errors.
fn best_tied(ranked: &[(u32, f64)]) -> &[(u32, f64)] {
    let Some(&(_, fewest)) = ranked.first() else {
        return &[];
    };
    let tied = ranked.iter().take_while(|c| c.1 == fewest).count();
    &ranked[..tied]
}

/// Fragments shared among a set of species: the EM's read groups, each
/// fragment's group, and each group's bases.
struct SharedReads {
    groups: ReadGroups,
    group_of: Vec<Option<usize>>,
    group_bases: Vec<f64>,
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

impl Profiler {
    /// Profile the loaded mappings `data` on the run's threads.
    pub(in crate::profiler) fn profile(&self, data: MappingData) -> Result<()> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.cfg.threads.max(1))
            .build()?
            .install(|| profile_levels(self, data))
    }
}

fn profile_levels(profiler: &Profiler, data: MappingData) -> Result<()> {
    let cfg = &profiler.cfg;
    let lineages = Lineages::new(&profiler.taxonomy);
    let names = Names {
        taxonomy: &profiler.taxonomy,
        genomes: &data.assembly_names,
        lineages: &lineages,
    };
    let policy = EvidencePolicy::from_config(cfg).ok_or_else(|| {
        anyhow::anyhow!(
            "the profiler takes Illumina reads, single-end or paired, or single-end long reads"
        )
    })?;
    // With `--save-classification` (not on a fractional database), each level
    // keeps how it shares the fragments.
    let classify_path = classification_path(cfg, &data);
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
    let em = EmSettings {
        penalty: likelihood_em::error_penalty(error_rate),
        max_iterations: cfg.max_em_iterations,
    };
    log::info!(
        target: "PROFILER",
        "{} fragments with candidates; error rate {:.5}, error penalty {:.3}.",
        fragments.len(),
        error_rate,
        em.penalty
    );
    let run = Run {
        cfg,
        data: &data,
        names: &names,
        policy,
        long_reads: policy.mode == EvidenceMode::Long,
        em,
        fragments: &fragments,
        labels: &labels,
        keep_shares: classify_path.is_some(),
    };

    let mut species = run.judge_species();
    let even = run.keep_even_coverage(&mut species);
    // What became of each measured species, for the decisions table (empty
    // when it is not written); the identity rules and the genus reports below
    // refine it.
    let mut outcomes = if cfg.save_decisions {
        tables::judged_species(&species.genomes, &species.best_fit, &species.looks(&even))
    } else {
        HashMap::new()
    };
    let kept = run.share_kept(&species);
    let (above, to_credit) = run.place_above_species(&species, &kept, unmatched, even.len());
    let genome_level = run.genome_level(error_rate);
    let genome_result = genomes::resolve(
        &data,
        &fragments,
        &kept.species_of(),
        &labels,
        &genome_level,
    );
    log::info!(
        target: "PROFILER",
        "Baseline identity {:.5} (median identity of the species to their best-fit genome); \
         {:.5} over strain-level matches.",
        genome_result.baseline_identity,
        genome_result.strain_baseline_identity
    );
    let reported = species.reported(&genome_result.species);
    let amounts = run.amounts(&species, &kept, &genome_result, &reported, &to_credit);

    let rows = run.species_rows(
        &species,
        (&kept, &genome_result),
        &reported,
        &amounts,
        &mut outcomes,
    );
    let mut blocks = rows.blocks;
    let higher_rows = above
        .iter()
        .flat_map(|found| &found.rows)
        .map(|row| {
            (
                row.taxid,
                profile_row(
                    names.taxon(row.taxid),
                    row.taxid,
                    row.rank.as_str(),
                    (row.reads, row.depth),
                    &row.measured.coverage,
                    names.genome_entry(row.measured.genome()),
                    (f64::NAN, f64::NAN),
                ),
            )
        })
        .collect();
    add_rows_above_species(&mut blocks, higher_rows, rows.genus_parts);
    let dropped = run.drop_rows(&mut blocks)?;
    // By the taxon's total depth, then taxid; abundance is then summed in this
    // order.
    blocks.sort_by(|a, b| b.depth().total_cmp(&a.depth()).then(a.taxid.cmp(&b.taxid)));
    set_abundance(&mut blocks);

    // The outputs.
    let finals = cfg.save_decisions.then(|| FinalRows::of(&blocks));
    if let Some(path) = &classify_path {
        let written = classify::Classifications {
            fragments: &fragments,
            group_of: &kept.group_of,
            posteriors: &kept.em.posteriors,
            genome_level: &genome_result.species,
            genus_members: &rows.genus_members,
            higher: above.as_ref().map_or(&[][..], |found| &found.shares),
            credited: &amounts.credited.fragments,
        }
        .write(path, &blocks, &data)?;
        log::info!(
            target: "PROFILER",
            "{written} fragments written to the per-read classifications."
        );
    }
    let profile_rows = blocks
        .into_iter()
        .flat_map(|block| block.taxon.into_iter().chain(block.genomes))
        .collect::<Vec<_>>();
    let unclassified = estimate_unclassified(
        cfg,
        &data,
        &profile_rows,
        &amounts.genome_reads,
        above.as_ref(),
        run.long_reads,
    );
    let header = TableHeader::of(cfg);
    let columns = if cfg.host.is_some() {
        [&PROFILE_COLUMNS[..], &HOST_COLUMNS[..]].concat()
    } else {
        PROFILE_COLUMNS.to_vec()
    };
    write_profile(
        cfg,
        &header.lines("profile", &columns),
        &profile_rows,
        unclassified,
    )?;
    if let Some(finals) = finals {
        let decisions_table = Decisions {
            genomes: &species.genomes,
            counts: &species.counts,
            tests: &species.tests,
            outcomes: &outcomes,
            reported: &reported,
            strains: cfg.genome_level.then(|| StrainLevel {
                rows: &genome_result.diagnostics,
                species: &genome_result.species,
                tests: genome_level.tests(1),
                strain_taxids: &data.strain_taxids,
            }),
            genus_members: &rows.genus_members,
            higher: above.as_ref(),
            dropped: &dropped,
            labels: &labels,
            named_first: &species.named_first,
        };
        let path = format!("{}.decisions.tsv", profile_prefix(cfg));
        tables::write_file(&path, |out| {
            decisions_table.write(out, &header, &names, &finals)
        })?;
    }
    Ok(())
}

/// What every step of a run reads: its settings and evidence gates, its
/// mappings and fragments, the names and lineages of taxa, and the units
/// labelled above species.
struct Run<'a> {
    cfg: &'a ProfilerConfig,
    data: &'a MappingData,
    names: &'a Names<'a>,
    policy: EvidencePolicy,
    long_reads: bool,
    em: EmSettings,
    fragments: &'a [Fragment],
    labels: &'a HashMap<u32, ProfileRank>,
    /// Keep how each level shares the fragments, for the per-read
    /// classifications.
    keep_shares: bool,
}

/// The species level after its looks: the measured genomes, the fragment
/// counts behind the tests and the tests, the best-fit genome of each kept
/// species, what each look kept or removed, and the first EM's iterations.
struct SpeciesJudgment {
    genomes: HashMap<u32, GenomeEvidence>,
    counts: FragmentCounts,
    tests: CoverageTests,
    best_fit: HashMap<u32, u32>,
    /// Kept on the second look, with the fragments each took.
    taken: HashMap<u32, Vec<usize>>,
    /// Kept on their own fragments ([`own_evidence_look`]), with whether each
    /// is reported at its genus.
    rescued: HashMap<u32, bool>,
    /// Kept viruses removed by the prophage check.
    removed_viruses: Vec<u32>,
    /// Kept units labelled above species removed as the kept named units of
    /// their taxon explain their reads, with those units.
    named_first: HashMap<u32, Vec<u32>>,
    first_iterations: usize,
}

impl SpeciesJudgment {
    /// The looks for the decisions table, with the species kept for even
    /// coverage (`even`).
    fn looks<'a>(&'a self, even: &'a [u32]) -> tables::Looks<'a> {
        tables::Looks {
            second_look: &self.taken,
            own_evidence: &self.rescued,
            prophage: &self.removed_viruses,
            named_first: &self.named_first,
            even_coverage: even,
        }
    }

    /// The genome each kept species is reported on: its best-fit genome when
    /// the genome level (`found`) found one among its measured genomes,
    /// otherwise its best-fit genome of the species level.
    fn reported(&self, found: &HashMap<u32, genomes::SpeciesGenomes>) -> HashMap<u32, u32> {
        self.best_fit
            .iter()
            .map(|(&taxid, &fallback)| {
                let id = found
                    .get(&taxid)
                    .map(|g| g.best_fit)
                    .filter(|id| self.genomes.contains_key(id))
                    .unwrap_or(fallback);
                (taxid, id)
            })
            .collect()
    }
}

/// The second EM, over the kept species: each fragment's read group (`None`
/// without a kept candidate), each group's bases, and the EM's result.
struct KeptShares {
    group_of: Vec<Option<usize>>,
    group_bases: Vec<f64>,
    em: likelihood_em::EmResult,
}

impl KeptShares {
    /// Each fragment's species in the second EM, when it gives the fragment
    /// to one species ([`RESOLVE_POSTERIOR`]), with that species' share.
    fn species_of(&self) -> Vec<Option<(u32, f64)>> {
        self.group_of
            .iter()
            .map(|group| group.and_then(|g| resolved_share(&self.em.posteriors[g])))
            .collect()
    }

    /// Each species' share of the fragments' bases.
    fn unit_bases(&self) -> HashMap<u32, f64> {
        let mut bases: HashMap<u32, f64> = HashMap::new();
        for (posteriors, &group_total) in self.em.posteriors.iter().zip(&self.group_bases) {
            for &(unit, posterior) in posteriors {
                *bases.entry(unit).or_insert(0.0) += group_total * posterior;
            }
        }
        bases
    }
}

/// What the kept species' rows hold beyond the second EM's reads: each
/// species' share of the fragments' bases and the depth of a species on a
/// marker genome, with the fragments credited to kept species; the marker
/// genome reads for the unclassified estimate; and the species the identity
/// rules leave out or report at their genus.
struct SpeciesAmounts {
    bases: HashMap<u32, f64>,
    marker_depth: HashMap<u32, f64>,
    genome_reads: HashMap<u32, GenomeReads>,
    left_out: HashSet<u32>,
    low_at_genus: HashSet<u32>,
    credited: Credited,
}

impl SpeciesAmounts {
    /// The reads of species `taxid`: its weight in the second EM `em` and the
    /// fragments credited to it.
    fn reads(&self, em: &likelihood_em::EmResult, taxid: u32) -> f64 {
        em.weight.get(&taxid).copied().unwrap_or(0.0)
            + self.credited.reads.get(&taxid).copied().unwrap_or(0.0)
    }
}

/// The rows of the kept species: each species' block with its genome rows;
/// the species reported at their genus as (genus, species, row); and for the
/// decisions table and the per-read classifications, (genus, species, reads).
struct SpeciesRows {
    blocks: Vec<TaxonBlock>,
    genus_parts: Vec<(u32, u32, ProfileRow)>,
    genus_members: Vec<(u32, u32, f64)>,
}

impl Run<'_> {
    /// The species level: the first EM over every candidate species, the
    /// first judgment of each measured genome, then the looks
    /// ([module documentation](self)), and each kept species' best-fit genome.
    fn judge_species(&self) -> SpeciesJudgment {
        let (cfg, data, policy, fragments) = (self.cfg, self.data, self.policy, self.fragments);
        // First EM over every candidate species; a fragment that it gives to
        // one species with posterior >= RESOLVE_POSTERIOR adds to the coverage
        // of that species' genomes.
        let all = share_reads(fragments, |_| true, self.em.penalty);
        let first_em = self.em.run(&all.groups);
        // Only the fragments' read groups are read again.
        let SharedReads { group_of, .. } = all;
        let resolved_taxid = first_em
            .posteriors
            .iter()
            .map(|posteriors| resolved_unit(posteriors))
            .collect::<Vec<_>>();
        let first = FirstEm::new(fragments, &group_of, &first_em, &resolved_taxid);
        let mut counts = first.counts();
        let resolved_fragments = first
            .won()
            .map(|(fragment, taxid)| (fragments[fragment].read, taxid))
            .collect::<Vec<_>>();
        let mut genomes: HashMap<u32, GenomeEvidence> = HashMap::new();
        fold_candidate_of(
            data,
            policy,
            &resolved_fragments,
            (|r| r.0, |r| r.1),
            |_, candidate| species_hits(candidate),
            |r, hits| add_hits(&mut genomes, r.1, hits),
        );
        drop(resolved_fragments);

        // --min-oebr and --min-coebr set the bounds of both levels.
        let min_coebr = if self.long_reads {
            SPECIES_MIN_COEBR_LONG_READS
        } else {
            SPECIES_MIN_COEBR_SHORT_READS
        };
        let tests = CoverageTests {
            multiple_tests: genomes.len().max(1),
            chunk_cutoff: data.min_chunk_breadth,
            ratio_bounds: (
                cfg.min_oebr.unwrap_or(SPECIES_MIN_OEBR),
                cfg.min_coebr.unwrap_or(min_coebr),
            ),
            markers: data.marker_model.clone(),
            min_marker_genes: cfg.min_marker_genes,
        };
        genomes.par_iter_mut().for_each(|(&id, genome)| {
            let won = counts.won.get(&genome.taxid).copied().unwrap_or(0);
            let candidates = counts.candidates.get(&genome.taxid).copied().unwrap_or(0);
            tests.judge(id, genome, won, candidates);
        });
        let labels = self.labels;
        let ancestors = |taxid: u32| self.names.lineages.ancestors(taxid);
        let taken = second_look(
            data,
            policy,
            &first,
            &tests,
            &mut counts,
            &mut genomes,
            (labels, &ancestors),
        );
        log::info!(
            target: "PROFILER",
            "Second look kept {} removed species on {} fragments their removed sisters tied on.",
            taken.len(),
            taken.values().map(Vec::len).sum::<usize>()
        );
        let rescued =
            own_evidence_look(&first, &tests, &mut counts, &mut genomes, |taxid, list| {
                coverage_of(data, policy, fragments, &tests, taxid, list)
            });
        log::info!(
            target: "PROFILER",
            "{} removed species kept on fragments that no kept species explains; {} of them \
             reported at their genus, as their reads tie with other species along their genome.",
            rescued.len(),
            rescued.values().filter(|&&genus| genus).count()
        );
        let removed_viruses = prophage_look(
            &first,
            &taken,
            &tests,
            &mut genomes,
            (
                |taxid| domain(self.names.lineages, taxid),
                (labels, &ancestors),
            ),
            |taxid, list| evidence_of(data, policy, fragments, taxid, list),
        );
        log::info!(
            target: "PROFILER",
            "{} kept virus species removed: their own fragments do not cover the rest of their \
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
                (labels, &ancestors),
                |taxid, list| evidence_of(data, policy, fragments, taxid, list),
            );
            log::info!(
                target: "PROFILER",
                "{} kept units labelled above species removed: their own fragments do not cover \
                 the rest of their genome, beyond the part that kept named units of their taxon explain.",
                named_first.len()
            );
            named_first
        };
        let first_iterations = first_em.iterations;
        // Best-fit genome of each kept species: its kept genome where most of
        // its fragments fit best, then the one with most fragments.
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
        SpeciesJudgment {
            genomes,
            counts,
            tests,
            best_fit,
            taken,
            rescued,
            removed_viruses,
            named_first,
            first_iterations,
        }
    }

    /// Where the data carries the inputs of the levels above species (not
    /// with `--strict-coverage`) and the database is not fractional: a
    /// removed species whose reads cover its strongest genome evenly is kept
    /// on that genome. The placements on the measured genomes are not read
    /// again. Returns the species kept for even coverage.
    fn keep_even_coverage(&self, species: &mut SpeciesJudgment) -> Vec<u32> {
        let data = self.data;
        let full_above = data.above_species.as_ref().filter(|_| !data.fractional);
        let strongest = full_above.map_or_else(HashMap::new, |above| {
            strongest_genomes(&species.genomes, &above.layout)
        });
        for genome in species.genomes.values_mut() {
            genome.assembly.intervals = Vec::new();
        }
        let even = full_above.map_or_else(Vec::new, |above| {
            above.even_coverage_species(
                strongest
                    .iter()
                    .map(|(&taxid, genome)| (taxid, &genome.coverage, genome.kept)),
                &data.genus_of_species,
                |taxid| data.terminal_ranks.get(&taxid).copied(),
                data.min_chunk_breadth,
            )
        });
        for &taxid in &even {
            species.best_fit.insert(taxid, strongest[&taxid].genome);
        }
        even
    }

    /// The second EM, over the kept species.
    fn share_kept(&self, species: &SpeciesJudgment) -> KeptShares {
        let shared = share_reads(
            self.fragments,
            |taxid| species.best_fit.contains_key(&taxid),
            self.em.penalty,
        );
        let em = self.em.run(&shared.groups);
        let without_kept = shared.group_of.iter().filter(|g| g.is_none()).count();
        log::info!(
            target: "PROFILER",
            "{} species kept, {} read groups, EM {} + {} iterations; {} fragments without a kept \
             species.",
            species.best_fit.len(),
            shared.groups.len(),
            species.first_iterations,
            em.iterations,
            without_kept
        );
        KeptShares {
            group_of: shared.group_of,
            group_bases: shared.group_bases,
            em,
        }
    }

    /// The fragments that no kept species takes, placed above species
    /// ([`higher`]), when the levels above species run (not with
    /// `--strict-coverage`); not on a fractional database, where the reads of a
    /// genome missing from the database land on the indexed segments of many
    /// relatives and cannot be counted as full-index equivalents. On noisy long reads, the fragments without a species-level
    /// candidate (`unmatched`) whose best match is a kept species are that
    /// species' reads, credited to it ([`credited_reads`]), not reads of the
    /// levels above species. Returns the levels above species when they ran,
    /// and the fragments to credit.
    fn place_above_species(
        &self,
        species: &SpeciesJudgment,
        kept: &KeptShares,
        unmatched: Vec<usize>,
        even: usize,
    ) -> (Option<higher::HigherLevels>, Vec<usize>) {
        let (cfg, data) = (self.cfg, self.data);
        let floor_policy = self.policy.floor(cfg.genus_identity);
        let (unmatched, to_credit) = if self.credits_noisy_reads() {
            let kept = species.best_fit.keys().copied().collect::<HashSet<_>>();
            split_best_kept(data, floor_policy, &unmatched, &kept)
        } else {
            (unmatched, Vec::new())
        };
        let above = data.above_species.as_ref();
        if above.is_some() && data.fractional {
            log::info!(
                target: "PROFILER",
                "On a fractional database, reads not placed at a species are not reported \
                 above species."
            );
        }
        let higher = above.filter(|_| !data.fractional).map(|above| {
            let not_kept = self
                .fragments
                .iter()
                .zip(&kept.group_of)
                .filter(|(_, group)| group.is_none())
                .map(|(fragment, _)| fragment.read)
                .collect::<Vec<_>>();
            // Each kept unit with the model of its best-fit genome, its
            // fragments and their whole-genome scale.
            let markers = data.marker_model.as_deref();
            let spill_bounds = higher::spill_bounds(
                species.best_fit.iter().map(|(&taxid, &genome)| {
                    (
                        taxid,
                        species.tests.model_of(genome),
                        kept.em.weight.get(&taxid).copied().unwrap_or(0.0),
                        markers.and_then(|model| model.declared_genome_scale(genome)),
                    )
                }),
                self.names.lineages,
            );
            let level = higher::Level {
                data,
                tests: &species.tests,
                layout: &above.layout,
                below_species: (&unmatched, floor_policy),
                not_kept: (&not_kept, self.policy),
                em: self.em,
                keep_shares: self.keep_shares,
                spill_bounds: &spill_bounds,
            };
            let found = higher::place(&level, self.names.lineages);
            let genus_rows = found
                .rows
                .iter()
                .filter(|row| row.rank == higher::Rank::Genus)
                .count();
            log::info!(
                target: "PROFILER",
                "{} species kept for even coverage; above species, {} fragments without a \
                 species-level match and {} of species not kept: {} genus rows ({} candidates), {} rows above \
                 genus ({} candidates), {} fragments unexplained.",
                even,
                found.fragments[0],
                found.fragments[1],
                genus_rows,
                found.genera.len(),
                found.rows.len() - genus_rows,
                found.taxa.len(),
                found.unexplained
            );
            found
        });
        (higher, to_credit)
    }

    /// Whether the fragments without a species-level candidate whose best
    /// match is a kept species count for it: on noisy long reads, where many
    /// reads of a present species align below the species gates, when the
    /// levels above species run.
    fn credits_noisy_reads(&self) -> bool {
        self.cfg.above_species && self.cfg.noisy_long_reads
    }

    /// The settings of the genome level, at the sample's `error_rate`.
    fn genome_level(&self, error_rate: f64) -> genomes::GenomeLevel {
        let cfg = self.cfg;
        let min_coebr = if self.long_reads {
            GENOME_MIN_COEBR_LONG_READS
        } else {
            GENOME_MIN_COEBR_SHORT_READS
        };
        genomes::GenomeLevel {
            policy: self.policy,
            error_rate,
            em: self.em,
            report_strains: cfg.genome_level,
            keep_shares: self.keep_shares,
            chunk_cutoff: self.data.min_chunk_breadth,
            ratio_bounds: (
                cfg.min_oebr.unwrap_or(GENOME_MIN_OEBR),
                cfg.min_coebr.unwrap_or(min_coebr),
            ),
        }
    }

    /// What the kept species' rows hold beyond the second EM's reads
    /// ([`SpeciesAmounts`]): the marker evidence of the species `reported` on
    /// marker genomes, the identity rules for relatives, and the fragments
    /// `to_credit`.
    fn amounts(
        &self,
        species: &SpeciesJudgment,
        kept: &KeptShares,
        genome_result: &genomes::GenomeResult,
        reported: &HashMap<u32, u32>,
        to_credit: &[usize],
    ) -> SpeciesAmounts {
        let (cfg, data) = (self.cfg, self.data);
        // Each species' share of the fragments' bases.
        let mut bases = kept.unit_bases();
        let MarkerEvidence {
            depth: mut marker_depth,
            genome_reads,
        } = data
            .marker_model
            .as_deref()
            .map_or_else(MarkerEvidence::default, |model| {
                marker_evidence(
                    data,
                    self.policy,
                    self.fragments,
                    (&kept.group_of, &kept.em),
                    reported,
                    (model, cfg.estimate_unclassified.then_some(self.long_reads)),
                )
            });
        let (left_out, low_at_genus) =
            self.identity_rules(species, genome_result, reported, (&bases, &marker_depth));
        // The fragments without a species-level candidate whose best match is
        // a kept species, too noisy for the species gates, count for it.
        let credited = if self.credits_noisy_reads() {
            let credited = credited_reads(
                data,
                self.policy.floor(cfg.genus_identity),
                to_credit,
                (&kept.em.weight, reported),
                self.keep_shares,
            );
            log::info!(
                target: "PROFILER",
                "{} fragments without a species-level candidate credited to the kept species \
                 that match them best.",
                credited.count
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
        SpeciesAmounts {
            bases,
            marker_depth,
            genome_reads,
            left_out,
            low_at_genus,
            credited,
        }
    }

    /// The identity rules for relatives ([`identity_relatives`]) over the
    /// kept species (not the units labelled above species) with their mean
    /// identity to their `reported` genome and their depth from their
    /// `(bases, marker_depth)`. Reads from the unindexed parts of a genome land
    /// on its relatives: the margin is set only where genomes are indexed in
    /// part, and the relative identity gap where it is on. Returns
    /// the species left out and those reported at their genus.
    fn identity_rules(
        &self,
        species: &SpeciesJudgment,
        genome_result: &genomes::GenomeResult,
        reported: &HashMap<u32, u32>,
        (bases, marker_depth): (&HashMap<u32, f64>, &HashMap<u32, f64>),
    ) -> (HashSet<u32>, HashSet<u32>) {
        let (cfg, data) = (self.cfg, self.data);
        let margin = cfg.identity_margin;
        let gap = cfg.relative_identity_gap;
        if margin.is_none() && gap.is_none() {
            return (HashSet::new(), HashSet::new());
        }
        let judged = reported
            .iter()
            .filter(|(taxid, _)| !self.labels.contains_key(taxid))
            .filter_map(|(&taxid, &id)| {
                let identity = genome_result.species.get(&taxid)?.mean_identity;
                let depth = marker_depth.get(&taxid).copied().unwrap_or_else(|| {
                    bases.get(&taxid).copied().unwrap_or(0.0)
                        / species.genomes[&id].assembly.asm_len.max(1.0)
                });
                identity.is_finite().then_some((taxid, (identity, depth)))
            })
            .collect::<HashMap<_, _>>();
        let found = identity_relatives(
            &data.genus_of_species,
            &judged,
            genome_result.baseline_identity,
            (margin, gap),
        );
        log::info!(
            target: "PROFILER",
            "{} species left out and {} reported at their genus by the identity rules for \
             relatives (margin {}, relative identity gap {}).",
            found.0.len(),
            found.1.len(),
            margin.map_or_else(|| "off".to_string(), |m| m.to_string()),
            gap.map_or_else(|| "off".to_string(), |g| g.to_string())
        );
        found
    }

    /// The rows of the kept species `reported` on their genomes, each followed
    /// by its genomes reported on their own. The rows of a species are
    /// exclusive: a genome row counts each fragment by the species' share of
    /// it, the species row keeps what no reported genome holds and is left
    /// out when less than one read is left. A species left out by the identity
    /// rules has no rows, and one reported at its genus a part of its genus
    /// row; `outcomes` records both.
    fn species_rows(
        &self,
        species: &SpeciesJudgment,
        (kept, genome_result): (&KeptShares, &genomes::GenomeResult),
        reported: &HashMap<u32, u32>,
        amounts: &SpeciesAmounts,
        outcomes: &mut HashMap<u32, SpeciesDecision>,
    ) -> SpeciesRows {
        let (data, names, labels) = (self.data, self.names, self.labels);
        let mut rows = SpeciesRows {
            blocks: Vec::new(),
            genus_parts: Vec::new(),
            genus_members: Vec::new(),
        };
        for (&taxid, &assembly_id) in reported {
            if amounts.left_out.contains(&taxid) {
                outcomes
                    .entry(taxid)
                    .and_modify(|d| *d = SpeciesDecision::LeftOut);
                continue;
            }
            let found = genome_result.species.get(&taxid);
            let genome = &species.genomes[&assembly_id];
            // Full-index-equivalent fragments, and the depth of the species'
            // bases on its genome; a species on a marker genome has no genome
            // rows and takes the depth its fragments add there.
            let scale = genome.assembly.sampling_scale.max(1.0);
            let depth = |bases: f64| {
                amounts
                    .marker_depth
                    .get(&taxid)
                    .copied()
                    .unwrap_or(bases / genome.assembly.asm_len.max(1.0))
            };
            let bases = amounts.bases.get(&taxid).copied().unwrap_or(0.0);
            let reads = amounts.reads(&kept.em, taxid);
            let coverage = genome.coverage.as_ref().expect("measured above");
            let strain_genomes = found.map_or(&[][..], |g| &g.strains[..]);
            // Every row of the profile has the identity columns, NA when the
            // species has no fragments of its own.
            let identity = found.map_or((f64::NAN, f64::NAN), |g| {
                (
                    g.mean_identity,
                    estimated_ani(genome_result.baseline_identity, g.mean_identity),
                )
            });
            let rank = unit_rank(labels.get(&taxid).copied());
            let tie_at_genus = species.rescued.get(&taxid) == Some(&true);
            let reported_at_genus = !labels.contains_key(&taxid)
                && (tie_at_genus || amounts.low_at_genus.contains(&taxid));
            if let Some(&genus) = data
                .genus_of_species
                .get(&taxid)
                .filter(|_| reported_at_genus)
            {
                outcomes.entry(taxid).and_modify(|d| {
                    *d = if tie_at_genus {
                        SpeciesDecision::AtGenusTie
                    } else {
                        SpeciesDecision::AtGenusIdentity
                    }
                });
                if self.cfg.save_decisions || self.keep_shares {
                    rows.genus_members.push((genus, taxid, reads));
                }
                rows.genus_parts.push((
                    genus,
                    taxid,
                    profile_row(
                        names.taxon(genus),
                        genus,
                        "genus",
                        (reads * scale, depth(bases)),
                        coverage,
                        names.genome_entry(assembly_id),
                        (f64::NAN, f64::NAN),
                    ),
                ));
                continue;
            }
            let species_row = species_remainder(
                reads,
                bases,
                strain_genomes.iter().map(|s| (s.reads, s.bases)),
            )
            .map(|(reads, species_bases)| {
                profile_row(
                    names.taxon(taxid),
                    taxid,
                    rank,
                    (reads * scale, depth(species_bases)),
                    coverage,
                    names.genome_entry(assembly_id),
                    identity,
                )
            });
            let mut strains = Vec::new();
            for strain in strain_genomes {
                let strain_ani =
                    estimated_ani(genome_result.strain_baseline_identity, strain.mean_identity);
                let (name, strain_taxid) = strain_label(
                    names,
                    data.strain_taxids[strain.assembly_id as usize],
                    taxid,
                    names.genome(strain.assembly_id),
                );
                strains.push(profile_row(
                    name,
                    strain_taxid,
                    "strain",
                    (
                        strain.reads
                            * species
                                .genomes
                                .get(&strain.assembly_id)
                                .map_or(1.0, |g| g.assembly.sampling_scale.max(1.0)),
                        strain.bases / strain.asm_len.max(1.0),
                    ),
                    &strain.coverage,
                    names.genome_entry(strain.assembly_id),
                    (strain.mean_identity, strain_ani),
                ));
            }
            rows.blocks.push(TaxonBlock {
                taxid,
                taxon: species_row,
                genomes: strains,
                breadth: coverage.breadth,
                unit: true,
            });
        }
        rows
    }

    /// `--host` drops rows first ([`host`]), then `--min-reads` and
    /// `--min-breadth` drop taxa ([`drop_below_minimums`]). Returns what was
    /// dropped.
    fn drop_rows(&self, blocks: &mut Vec<TaxonBlock>) -> Result<Dropped> {
        let cfg = self.cfg;
        let mut dropped = Dropped::default();
        if let Some(host) = &cfg.host {
            let table = load_pathogen_table(&host.pathogen_table)?;
            let rows = |blocks: &[TaxonBlock]| {
                blocks
                    .iter()
                    .map(|b| usize::from(b.taxon.is_some()) + b.genomes.len())
                    .sum::<usize>()
            };
            let before = rows(blocks);
            host::HostFilter {
                host: host.taxid,
                table: &table,
                lineages: self.names.lineages,
            }
            .apply(blocks, &mut dropped);
            log::info!(
                target: "PROFILER",
                "{} rows dropped by --host {} ({} taxa entirely): not pathogens of the host.",
                before - rows(blocks),
                host.taxid,
                dropped.taxa.len()
            );
        }
        let below = drop_below_minimums(
            blocks,
            &Minimums {
                reads: cfg.min_reads,
                breadth: cfg.breadth,
            },
        );
        if !below.is_empty() {
            log::info!(
                target: "PROFILER",
                "{} taxa dropped by --min-reads or --min-breadth; their reads go to no other taxon.",
                below.len()
            );
        }
        dropped.taxa.extend(
            below
                .into_iter()
                .map(|(taxid, below)| (taxid, Removal::BelowMinimum(below))),
        );
        Ok(dropped)
    }
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
        |&read, candidates| fragment_of(read, candidates, policy),
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

/// The errors of a candidate tied with the best (within `margin` of the
/// `fewest` errors) are the best's.
fn tied_errors(errors: f64, fewest: f64, margin: f64) -> f64 {
    if errors - fewest < margin {
        fewest
    } else {
        errors
    }
}

/// Aligned bases of the longest of `placements`.
fn longest_alignment<'p>(placements: impl Iterator<Item = &'p CandidatePlacement>) -> f64 {
    placements
        .map(|p| f64::from(p.both_mapped))
        .fold(0.0, f64::max)
}

/// The `candidates` of a fragment, each with its errors: its best
/// placement's estimated errors and the bases it aligns less than the longest
/// alignment of all the candidates ([`candidate_errors`]); in candidate order.
fn with_errors(candidates: &[TerminalCandidate]) -> Vec<(f64, &TerminalCandidate)> {
    let longest = longest_alignment(candidates.iter().map(|c| &c.best));
    candidates
        .iter()
        .map(|c| (candidate_errors(&c.best, longest), c))
        .collect()
}

/// A fragment's `candidates` ranked: each taxon with its errors
/// ([`with_errors`]), those within the tie margin of the best ([`tie_margin`]
/// of its aligned bases) with the best's errors, by errors, then taxid; and
/// the best candidate. `None` without a candidate.
fn ranked(
    candidates: &[TerminalCandidate],
    mode: EvidenceMode,
) -> Option<(Vec<(u32, f64)>, &TerminalCandidate)> {
    let mut ranked = with_errors(candidates);
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.taxid.cmp(&b.1.taxid)));
    let &(fewest, best) = ranked.first()?;
    let margin = tie_margin(mode, f64::from(best.best.both_mapped), LONG_TIE_SHARE);
    let taxa = ranked
        .iter()
        .map(|&(errors, candidate)| (candidate.taxid, tied_errors(errors, fewest, margin)))
        .collect();
    Some((taxa, best))
}

/// The candidate of taxon `taxid` among a fragment's `candidates`.
fn candidate_of(candidates: &[TerminalCandidate], taxid: u32) -> Option<&TerminalCandidate> {
    candidates.iter().find(|c| c.taxid == taxid)
}

/// Read `read` as a fragment with its candidate units ranked ([`ranked`]),
/// with the errors and aligned bases of its best candidate; `None` without
/// one. A candidate is a species, or the taxon above species that its
/// references are labelled with.
fn fragment_of(
    read: usize,
    candidates: &[TerminalCandidate],
    policy: EvidencePolicy,
) -> Option<(Fragment, (f32, f32))> {
    let (candidates, best) = ranked(candidates, policy.mode)?;
    let fewest = candidates[0].1;
    let fragment = Fragment {
        read,
        candidates,
        bases: f64::from(best.best.r1_mapped),
    };
    Some((fragment, (fewest as f32, best.best.both_mapped as f32)))
}

/// Group the fragments for the EM over the species that `allowed` admits,
/// dropping candidates whose likelihood is below [`MIN_RELATIVE_LIKELIHOOD`]
/// of the fragment's best allowed candidate.
fn share_reads(fragments: &[Fragment], allowed: impl Fn(u32) -> bool, penalty: f64) -> SharedReads {
    share_reads_after(ReadGroups::default(), fragments, allowed, penalty)
}

/// [`share_reads`] after the read groups `groups`, whose reads have no bases.
fn share_reads_after(
    groups: ReadGroups,
    fragments: &[Fragment],
    allowed: impl Fn(u32) -> bool,
    penalty: f64,
) -> SharedReads {
    let max_excess_errors = -MIN_RELATIVE_LIKELIHOOD.ln() / penalty;
    let mut shared = SharedReads {
        group_of: Vec::with_capacity(fragments.len()),
        group_bases: vec![0.0; groups.len()],
        groups,
    };
    let mut candidates: likelihood_em::Candidates = Vec::new();
    for fragment in fragments {
        candidates.clear();
        candidates.extend(fragment.candidates.iter().filter(|c| allowed(c.0)));
        let fewest = candidates.iter().map(|c| c.1).fold(f64::INFINITY, f64::min);
        candidates.retain(|c| c.1 - fewest <= max_excess_errors);
        let slot = shared.groups.add(&mut candidates);
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
    assembly.add_interval(PresenceInterval::from(placement));
}

/// Identity of a placement: one less its estimated errors per aligned base.
fn placement_identity(placement: &CandidatePlacement) -> f64 {
    1.0 - placement.pair_errors / f64::from(placement.both_mapped).max(1.0)
}

/// The rank of a unit's rows: species, or the rank its references are
/// labelled with (`label`).
fn unit_rank(label: Option<ProfileRank>) -> &'static str {
    label.map_or("species", ProfileRank::as_str)
}

/// Estimated ANI of the sample's strain to a genome: its fragments' mean
/// `identity` to the genome relative to the sample's `baseline` identity, at
/// most 1; unknown (NaN) when the identity is.
fn estimated_ani(baseline: f64, identity: f64) -> f64 {
    let ani = 1.0 - (baseline - identity);
    // Not `f64::min`, which would turn NaN into 1.
    if ani > 1.0 {
        1.0
    } else {
        ani
    }
}

/// The identity rules for relatives, over the kept `species` with their mean
/// identity to their best-fit genome and their depth. With a `margin` (on a
/// fractional database, where the reads of a genome's unindexed parts land on
/// the indexed segments of its relatives), a species
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
/// in the final EM (each fragment's read group `group_of`, `em`) times the depth one
/// fragment adds on that genome, [`fragment_depth`]: its read length times the
/// genes it hits over the placements that the genome's genes offer a read of
/// its length. With `unclassified` (`--estimate-unclassified`; `Some(true)` for
/// long reads), each fragment resolved to such a species
/// ([`RESOLVE_POSTERIOR`]) also records its placements on marker genomes. Only
/// the fragments whose read group lists such a species are read again.
fn marker_evidence(
    data: &MappingData,
    policy: EvidencePolicy,
    fragments: &[Fragment],
    (group_of, em): (&[Option<usize>], &likelihood_em::EmResult),
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
        .zip(group_of)
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
                    let candidate = candidate_of(candidates, taxid)?;
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
        .placement_on(assembly_id)
        .copied()
        .unwrap_or(candidate.best);
    let mapped = match mode {
        EvidenceMode::Paired | EvidenceMode::Single => placement.r1_mapped,
        EvidenceMode::Long => placement.both_mapped,
    };
    placement_depth(&placement, mapped, Some(model))
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
    fold: impl FnMut(&I, T),
) {
    fold_built(
        items,
        |scratch, item| {
            let rows = &data.read_alignments[read_of(item)].alignments;
            extract(item, scratch.build(rows, policy))
        },
        fold,
    );
}

/// [`fold_candidates`] for a pass that reads one taxon's candidate per item
/// (`taxid_of`): only that candidate is built
/// ([`CandidateScratch::build_for`]).
fn fold_candidate_of<I: Sync, T: Send>(
    data: &MappingData,
    policy: EvidencePolicy,
    items: &[I],
    (read_of, taxid_of): (impl Fn(&I) -> usize + Sync, impl Fn(&I) -> u32 + Sync),
    extract: impl Fn(&I, Option<&TerminalCandidate>) -> T + Sync,
    fold: impl FnMut(&I, T),
) {
    fold_built(
        items,
        |scratch, item| {
            let rows = &data.read_alignments[read_of(item)].alignments;
            extract(item, scratch.build_for(rows, policy, taxid_of(item)))
        },
        fold,
    );
}

/// Run `work` on each of `items` on the profiler's threads, each thread with
/// its own candidate storage, and hand the results to `fold` in item order,
/// [`PARALLEL_BLOCK`] items at a time.
fn fold_built<I: Sync, T: Send>(
    items: &[I],
    work: impl Fn(&mut CandidateScratch, &I) -> T + Sync,
    mut fold: impl FnMut(&I, T),
) {
    for block in items.chunks(PARALLEL_BLOCK) {
        let results = block
            .par_iter()
            .map_init(CandidateScratch::default, |scratch, item| {
                work(scratch, item)
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

/// The placements and gene hits of a fragment's `candidate` of a species.
fn species_hits(candidate: Option<&TerminalCandidate>) -> SpeciesHits {
    candidate.map_or_else(SpeciesHits::default, |candidate| SpeciesHits {
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
            genome.assembly.add_interval(hit.interval);
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
    fold_candidate_of(
        data,
        policy,
        list,
        (|&f| fragments[f].read, |_| taxid),
        |_, candidate| species_hits(candidate),
        |_, hits| add_hits(genomes, taxid, hits),
    );
}

/// Fragments whose read group in the first EM lists each species, those
/// resolved to it, those it took on the second look, and those it won with no
/// kept species among their candidates.
#[derive(Default)]
struct FragmentCounts {
    candidates: HashMap<u32, u64>,
    won: HashMap<u32, u64>,
    second_look: HashMap<u32, u64>,
    own: HashMap<u32, u64>,
}

/// The species tests: enough coverage, evenly spread, and a large enough share
/// of the fragments listing the species won. A genome represented by marker
/// genes is measured by the `markers` model, its genes standing for the chunks
/// of a whole genome, and needs at least `min_marker_genes` genes hit.
struct CoverageTests {
    multiple_tests: usize,
    chunk_cutoff: f64,
    ratio_bounds: (f64, f64),
    markers: Option<Arc<MarkerModel>>,
    min_marker_genes: usize,
}

impl CoverageTests {
    /// The evidence model of genome `assembly_id`.
    fn model_of(&self, assembly_id: u32) -> higher::Model {
        match &self.markers {
            Some(model) if model.is_marker_genome(assembly_id) => higher::Model::Markers,
            _ => higher::Model::Genome,
        }
    }

    /// Coverage of genome `assembly_id`: gene by gene for a marker genome,
    /// otherwise over its sequence (also on a marker genome whose genes
    /// cannot hold a read).
    fn measure(&self, assembly_id: u32, assembly: &mut PresenceAssembly) -> EvidenceCoverage {
        self.markers
            .as_deref()
            .filter(|_| self.model_of(assembly_id) == higher::Model::Markers)
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
    /// The fragments of each read group.
    members: GroupMembers,
}

impl<'a> FirstEm<'a> {
    fn new(
        fragments: &'a [Fragment],
        group_of: &'a [Option<usize>],
        em: &'a likelihood_em::EmResult,
        resolved: &'a [Option<u32>],
    ) -> Self {
        Self {
            fragments,
            group_of,
            em,
            resolved,
            members: GroupMembers::new(group_of, em.posteriors.len()),
        }
    }

    /// Fragments whose read group lists each species and those resolved to
    /// it, counted per read group.
    fn counts(&self) -> FragmentCounts {
        let mut counts = FragmentCounts::default();
        for (group, (posteriors, &resolved)) in
            self.em.posteriors.iter().zip(self.resolved).enumerate()
        {
            let size = self.members.of(group).len() as u64;
            for &(taxid, _) in posteriors {
                *counts.candidates.entry(taxid).or_default() += size;
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

    /// The fragments resolved to each species that `wanted` admits, in
    /// fragment order.
    fn won_by(&self, wanted: impl Fn(u32) -> bool) -> HashMap<u32, Vec<usize>> {
        let mut won: HashMap<u32, Vec<usize>> = HashMap::new();
        for (index, taxid) in self.won().filter(|&(_, taxid)| wanted(taxid)) {
            won.entry(taxid).or_default().push(index);
        }
        won
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
    tests: &CoverageTests,
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
            .flat_map(|&group| first.members.of(group))
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

/// Own-evidence look at the removed species that cover one of their genomes well
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
    tests: &CoverageTests,
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
    let won = first.won_by(|taxid| judged.contains(&taxid));
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
    tests: &CoverageTests,
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
        .map(|virus| (virus, &hosts))
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
    tests: &CoverageTests,
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
    let inside = kept
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
    if inside.is_empty() {
        return HashMap::new();
    }
    let explainers = inside
        .iter()
        .map(|(&label, units)| (label, units))
        .collect::<HashMap<_, _>>();
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
    tests: &CoverageTests,
    genomes: &mut HashMap<u32, GenomeEvidence>,
    explainers: &HashMap<u32, &HashSet<u32>>,
    mut evidence: impl FnMut(u32, &[usize]) -> HashMap<u32, GenomeEvidence>,
) -> Vec<(u32, Vec<u32>)> {
    let mut held = first.won_by(|taxid| explainers.contains_key(&taxid));
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
    let mut covered_spans = other
        .intervals
        .iter()
        .flat_map(|r| [(r.reference_id, r.span), (r.reference_id, r.mate_span)])
        .collect::<Vec<_>>();
    merge_spans(&mut covered_spans);
    // Covered segments, each with the covered bases before it on its contig.
    let mut before = (None, 0);
    let segments = covered_spans
        .into_iter()
        .map(|(contig, span)| {
            if before.0 != Some(contig) {
                before = (Some(contig), 0);
            }
            let segment = (contig, span, before.1);
            before.1 += span.1 - span.0;
            segment
        })
        .collect::<Vec<_>>();
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
        rest.add_interval(PresenceInterval {
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
    tests: &CoverageTests,
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
    ratio.is_finite() && ratio >= lower && ratio <= MAX_BREADTH_RATIO
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
/// row) with its genome rows, or a row above species; with the taxon's
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
    /// Count a `row` above species of the unit's own taxon in the unit's row:
    /// reads and depth add, while coverage, best-fit genome and identity stay
    /// the unit's. A unit without a row of its own takes `row`.
    fn count(&mut self, row: ProfileRow) {
        match &mut self.taxon {
            Some(own) => own.add_reads_and_depth(&row),
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

/// `--min-reads` and `--min-breadth`: a taxon is dropped from the profile when
/// its rows hold fewer full-index-equivalent reads (a species counting with its
/// genome rows), or when its breadth is lower.
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
/// percent. The depths are summed in the order of `blocks`.
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
/// by assembly id.
#[derive(Default)]
struct Dropped {
    taxa: HashMap<u32, Removal>,
    genomes: HashMap<u32, Removal>,
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

/// A measured species' strongest genome, on which the even-coverage rule
/// judges the species: its coverage, with the window spread of its reads, and
/// whether it is kept.
#[derive(Clone, Copy)]
struct StrongestGenome {
    genome: u32,
    coverage: EvidenceCoverage,
    kept: bool,
}

impl StrongestGenome {
    /// Whether this genome is stronger than `other`: a kept genome first,
    /// then the larger fixed-chunk breadth, breadth, depth and read count.
    fn stronger_than(&self, other: &Self) -> bool {
        if self.kept != other.kept {
            return self.kept;
        }
        let (own, others) = (&self.coverage, &other.coverage);
        if own.fixed_chunk != others.fixed_chunk {
            return own.fixed_chunk > others.fixed_chunk;
        }
        if own.breadth != others.breadth {
            return own.breadth > others.breadth;
        }
        if own.depth != others.depth {
            return own.depth > others.depth;
        }
        own.count > others.count
    }
}

/// Each measured species' strongest genome
/// ([`StrongestGenome::stronger_than`]; the lowest assembly id when neither
/// of two is stronger), with the window spread of its reads on `layout` when
/// it is not kept, the only case the even-coverage rule reads it.
fn strongest_genomes(
    genomes: &HashMap<u32, GenomeEvidence>,
    layout: &GenomeLayout,
) -> HashMap<u32, StrongestGenome> {
    let mut strongest: HashMap<u32, StrongestGenome> = HashMap::new();
    for (&assembly_id, genome) in genomes {
        let Some(coverage) = genome.coverage else {
            continue;
        };
        let candidate = StrongestGenome {
            genome: assembly_id,
            coverage,
            kept: genome.kept,
        };
        let stronger = strongest.get(&genome.taxid).is_none_or(|current| {
            candidate.stronger_than(current)
                || (!current.stronger_than(&candidate) && assembly_id < current.genome)
        });
        if stronger {
            strongest.insert(genome.taxid, candidate);
        }
    }
    for genome in strongest.values_mut().filter(|genome| !genome.kept) {
        let intervals = &genomes[&genome.genome].assembly.intervals;
        genome.coverage.window_spread = layout.window_spread(genome.genome, intervals);
    }
    strongest
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
    fragment_of(read, candidates, policy).is_some_and(|(fragment, _)| {
        fragment
            .best_tied()
            .iter()
            .any(|(taxid, _)| kept.contains(taxid))
    })
}

/// The `reads` without a species-level candidate whose best match under
/// `policy` (the floor gates), or one tied with it, is not a species in
/// `kept`, and those whose best match is; both in read order.
fn split_best_kept(
    data: &MappingData,
    policy: EvidencePolicy,
    reads: &[usize],
    kept: &HashSet<u32>,
) -> (Vec<usize>, Vec<usize>) {
    let (mut others, mut best_kept) = (Vec::new(), Vec::new());
    fold_candidates(
        data,
        policy,
        reads,
        |&read| read,
        |&read, candidates| best_match_kept(read, candidates, policy, kept),
        |&read, is_kept| {
            if is_kept {
                best_kept.push(read);
            } else {
                others.push(read);
            }
        },
    );
    (others, best_kept)
}

/// What the fragments credited to kept species add to each: reads, bases,
/// and depth for a species reported on a marker genome; the fragments
/// credited; and when kept, each credited fragment's read with its shares.
#[derive(Default)]
struct Credited {
    reads: HashMap<u32, f64>,
    bases: HashMap<u32, f64>,
    marker_depth: HashMap<u32, f64>,
    count: usize,
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
            let (fragment, _) = fragment_of(read, candidates, policy)?;
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
                                let candidate = candidate_of(candidates, taxid)?;
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
            credited.count += 1;
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
fn domain(lineages: &Lineages, taxid: u32) -> Option<u32> {
    lineages
        .of(taxid)
        .iter()
        .map(|node| node.0)
        .find(|&id| matches!(id, VIRUSES | BACTERIA | ARCHAEA))
}

/// Add the rows above species (`rows`, by taxon) to `blocks`: a species
/// reported at its genus (`genus_parts`: genus, species, row) joins the row of
/// its genus, in taxid order, or makes it; a row of a taxon that a kept unit
/// is labelled with is counted in that unit's row.
fn add_rows_above_species(
    blocks: &mut Vec<TaxonBlock>,
    mut rows: Vec<(u32, ProfileRow)>,
    mut genus_parts: Vec<(u32, u32, ProfileRow)>,
) {
    genus_parts.sort_by_key(|part| (part.0, part.1));
    for (genus, _, part) in genus_parts {
        match rows.iter_mut().find(|row| row.0 == genus) {
            Some((_, row)) => row.add_reads_and_depth(&part),
            None => rows.push((genus, part)),
        }
    }
    for (taxid, row) in rows {
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
}

/// One profile row; `amount` is `(reads, depth)`.
fn profile_row(
    name: String,
    taxid: u32,
    rank: &str,
    (reads, depth): (f64, f64),
    coverage: &EvidenceCoverage,
    (assembly_id, genome): (u32, Arc<str>),
    identity: (f64, f64),
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
        best_fit: BestFitReport {
            assembly_id,
            genome,
            evidence_depth: coverage.depth,
            identity,
        },
        annotations: None,
    }
}

/// Depth that one fragment adds to its candidate: its mapped bases over the
/// genome length. A fragment on a marker gene stands for one fragment of the
/// whole genome: it adds its whole read length, also when the read is clipped
/// at a gene end, over the placements that the genome's genes offer a read of
/// its length, so that marker and genome taxa have comparable depths; a long
/// read adds this for every gene of the genome it covers, as the placements
/// count every gene. Marker genes that cannot hold a read keep mapped bases
/// over their summed length.
///
/// The fragment lies on genome `assembly_id` of `asm_len` bases with a read of
/// `read_len` bases on `genes` of its genes, `mapped` of them mapped.
fn fragment_depth(
    markers: Option<&MarkerModel>,
    (assembly_id, asm_len): (u32, f64),
    read_len: u32,
    genes: u32,
    mapped: f64,
) -> f64 {
    markers
        .and_then(|model| model.read_depth(assembly_id, read_len, genes))
        .unwrap_or_else(|| mapped / asm_len.max(1.0))
}

/// [`fragment_depth`] of a fragment placed at `placement` with `mapped`
/// bases mapped.
fn placement_depth(
    placement: &CandidatePlacement,
    mapped: u32,
    markers: Option<&MarkerModel>,
) -> f64 {
    fragment_depth(
        markers,
        (placement.assembly_id, placement.asm_len),
        placement.read_len,
        1 + placement.other_genes,
        f64::from(mapped),
    )
}

/// Record a resolved fragment's `placement` on a marker genome for the
/// unclassified estimate: for long reads (`long_reads`) the whole-genome reads
/// that its hits on the genome's genes, with `other_gene_hits`, stand for; for
/// other reads one fragment with its gene hits. Placements on other genomes
/// are not recorded.
fn add_marker_genome_reads(
    model: &MarkerModel,
    long_reads: bool,
    placement: &CandidatePlacement,
    other_gene_hits: &[GeneHit],
    reads: &mut HashMap<u32, GenomeReads>,
) {
    if long_reads {
        let genes = std::iter::once(placement.reference_id).chain(
            other_gene_hits
                .iter()
                .filter(|hit| hit.assembly_id == placement.assembly_id)
                .map(|hit| hit.interval.reference_id),
        );
        if let Some((hits, per_hit)) =
            model.whole_genome_reads(placement.assembly_id, placement.read_len, genes)
        {
            reads
                .entry(placement.assembly_id)
                .or_default()
                .add(hits, per_hit);
        }
    } else if model.is_marker_genome(placement.assembly_id) {
        reads
            .entry(placement.assembly_id)
            .or_default()
            .add_fragment(1 + placement.other_genes);
    }
}

/// The unclassified share (`--estimate-unclassified`): `None` when not
/// requested, `Some(None)` when it cannot be estimated.
fn estimate_unclassified(
    cfg: &ProfilerConfig,
    data: &MappingData,
    rows: &[ProfileRow],
    marker_genome_reads: &HashMap<u32, GenomeReads>,
    above: Option<&higher::HigherLevels>,
    long_reads: bool,
) -> Option<Option<unclassified::UnclassifiedEstimate>> {
    data.unclassified_species.as_ref().map(|species| {
        // Long reads differ in length and in the genes each covers, so
        // their scales come from the reads themselves.
        let mut scales = unclassified::marker_scales(
            species,
            data.marker_model.as_deref(),
            marker_genome_reads,
            long_reads,
        );
        // Rows above species on marker genomes count like species rows.
        scales.extend(
            above
                .iter()
                .flat_map(|found| &found.rows)
                .filter_map(|row| Some((row.taxid, row.measured.unclassified_scale()?))),
        );
        let estimate = unclassified::estimate(rows, &scales, data.total_reads);
        match estimate {
            Ok(estimate) => log::info!(
                target: "PROFILER",
                "Unclassified estimate: {:.2}% of {} input {} are not explained by the profile.",
                estimate.share * 100.0,
                data.total_reads,
                if cfg.is_paired { "fragments" } else { "reads" }
            ),
            Err(taxid) => match above.and_then(|found| found.row(taxid)) {
                Some(row) => log::warn!(
                    target: "PROFILER",
                    "Unclassified estimate unavailable: the {} row {taxid} rests on a marker genome that declares no genome size, cannot hold a read of this length, or has no resolved long read.",
                    row.rank.as_str()
                ),
                None => log::warn!(
                    target: "PROFILER",
                    "Unclassified estimate unavailable: species {taxid} is represented by marker genes that declare no genome size, cannot hold a read of this length, or have no resolved long read."
                ),
            },
        }
        estimate.ok()
    })
}

/// Path prefix of the profile files: `<prefix>.pathogen` with `--host`.
fn profile_prefix(cfg: &ProfilerConfig) -> String {
    if cfg.host.is_some() {
        format!("{}.pathogen", cfg.outprefix)
    } else {
        cfg.outprefix.clone()
    }
}

/// Write `<prefix>.profile.txt`: the `header` text, then the rows and, with
/// `--estimate-unclassified`, the unclassified row.
fn write_profile(
    cfg: &ProfilerConfig,
    header: &str,
    rows: &[ProfileRow],
    unclassified: Option<Option<unclassified::UnclassifiedEstimate>>,
) -> Result<()> {
    let mut file = File::create(format!("{}.profile.txt", profile_prefix(cfg)))?;
    file.write_all(header.as_bytes())?;
    let mut profile_writer = WriterBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_writer(file);
    for row in rows {
        row.write(&mut profile_writer)?;
    }
    if let Some(estimate) = unclassified {
        unclassified::write_row(&mut profile_writer, estimate)?;
    }
    profile_writer.flush()?;
    Ok(())
}

/// Path of the per-read classifications, when they are requested and can
/// be written for this database.
fn classification_path(cfg: &ProfilerConfig, data: &MappingData) -> Option<String> {
    cfg.writes_classifications(data.has_scaled_genomes)
        .then(|| {
            let classification_path = format!("{}.classify.txt.gz", cfg.outprefix);
            log::info!(
                target: "PROFILER",
                "Writing per-read classifications to {}",
                classification_path
            );
            classification_path
        })
}

/// Name and taxid of the genome row of genome `accession` of species `taxid`:
/// those of its `strain` when the genome declares one (named as the species
/// when the taxonomy has no name for it), otherwise the species' with the
/// accession.
fn strain_label(names: &Names, strain: Option<u32>, taxid: u32, accession: &str) -> (String, u32) {
    let species_name = names.taxonomy.get_name(taxid).unwrap_or("unknown");
    match strain {
        Some(strain) => (
            names
                .taxonomy
                .get_name(strain)
                .unwrap_or(species_name)
                .to_string(),
            strain,
        ),
        None => (format!("{species_name} {accession}"), taxid),
    }
}

/// Mapping data for the tests of the species and genome levels.
#[cfg(test)]
mod test_fixtures {
    use std::collections::HashMap;

    use crate::profiler::community::{MappingData, ReadBundle};
    use crate::profiler::evidence::{EvidenceAlignment, MappingRole};

    /// A single-end read of 100 bases aligned over its whole length.
    pub(super) fn row(taxid: u32, assembly: u32, identity: f32) -> EvidenceAlignment {
        EvidenceAlignment {
            profile_taxid: taxid,
            reference_id: assembly,
            assembly_id: assembly,
            role: MappingRole::Unpaired,
            pair_id: 0,
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
                    read_name: Box::from(format!("read{i}")),
                    alignments,
                })
                .collect(),
            min_chunk_breadth: 0.0,
            has_scaled_genomes: false,
            total_reads: 0,
            genus_of_species: HashMap::new(),
            terminal_ranks: HashMap::new(),
            marker_model: None,
            fractional: false,
            above_species: None,
            unclassified_species: None,
            assembly_names: Vec::new(),
            strain_taxids: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::test_fixtures::{data, row};
    use super::*;
    use crate::profiler::markers::{FragmentGeometry, FragmentLayout};
    use crate::taxonomy::Taxonomy;

    fn approx_equal(lhs: f64, rhs: f64) {
        assert!((lhs - rhs).abs() < 1e-9, "{lhs} is not {rhs}");
    }

    fn placement(assembly_id: u32, mapped: u32, read_len: u32) -> CandidatePlacement {
        CandidatePlacement {
            taxid: 1,
            reference_id: 0,
            assembly_id,
            sampling_scale: 1.0,
            asm_len: 2_000.0,
            r1_mapped: mapped,
            both_mapped: mapped,
            read_len,
            other_genes: 0,
            pair_errors: 0.0,
            long_score: 0.0,
            span: (1, 1 + mapped),
            mate_span: (0, 0),
            fixed_chunk: (1, 1 + mapped),
            flex_chunk: (1, 1 + mapped),
            input_order: 0,
        }
    }

    #[test]
    fn a_long_read_adds_depth_for_every_marker_gene_it_covers() {
        // Long reads of 3,000 bases on genes of 1,000 bases: 500 bases inside
        // a gene, 1,000 - 1,000 + 3,001 = 3,001 placements per gene.
        let geometry = FragmentGeometry::new(FragmentLayout::LongReads, 3_000, 3_000, 500, 0.2);
        let markers = MarkerModel::from_gene_lengths(&[(7, &[1_000, 1_000])], geometry);
        let mut read = placement(7, 900, 3_000);
        approx_equal(
            placement_depth(&read, 900, Some(&markers)),
            3_000.0 / 6_002.0,
        );
        read.other_genes = 1;
        approx_equal(
            placement_depth(&read, 900, Some(&markers)),
            6_000.0 / 6_002.0,
        );
        // A read of another length is placed by its own length.
        let short = placement(7, 900, 1_000);
        approx_equal(
            placement_depth(&short, 900, Some(&markers)),
            1_000.0 / (2.0 * (1_000.0 + 1_001.0 - 400.0)),
        );
    }

    #[test]
    fn a_fragment_on_a_marker_gene_adds_its_read_over_the_genes_placements() {
        // 150 bp pairs of 400 bases with lone mates: genes of 1,000 bases
        // offer 1,191 placements each; a gene of 100 bases holds no read.
        let geometry = FragmentGeometry::new(FragmentLayout::PairsAndLoneMates, 150, 400, 50, 0.7);
        let markers =
            MarkerModel::from_gene_lengths(&[(7, &[1_000, 1_000]), (8, &[100])], geometry);
        let depth = |assembly, mapped, markers| {
            placement_depth(&placement(assembly, mapped, 150), mapped, markers)
        };
        approx_equal(depth(7, 150, Some(&markers)), 150.0 / 2_382.0);
        // A read clipped at a gene end stands for its whole read.
        approx_equal(depth(7, 110, Some(&markers)), 150.0 / 2_382.0);
        // Genome references and marker genes without placements keep mapped
        // bases over the summed length.
        approx_equal(depth(9, 140, Some(&markers)), 140.0 / 2_000.0);
        approx_equal(depth(8, 90, Some(&markers)), 90.0 / 2_000.0);
        approx_equal(depth(7, 140, None), 140.0 / 2_000.0);
    }

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

    /// Single-reference genomes 0..=12 of 1 Mb each.
    fn layout() -> GenomeLayout {
        GenomeLayout::new((0..=12).map(|id| (id, 1_000_000)))
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
        let strongest = strongest_genomes(&genomes, &layout());
        assert_eq!(strongest[&1].genome, 10);
        assert!(strongest[&1].kept);
        assert_eq!(strongest[&2].genome, 12);
        assert!(!strongest[&2].kept);
        assert!(strongest[&2].coverage.window_spread.is_some());
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
                (1, Arc::from("GCF_1")),
                (f64::NAN, f64::NAN),
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
                (1, Arc::from("GCF_1")),
                (f64::NAN, f64::NAN),
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
    fn the_estimated_ani_is_at_most_one_and_unknown_without_an_identity() {
        assert!((estimated_ani(0.99, 0.98) - 0.99).abs() < 1e-12);
        assert_eq!(estimated_ani(0.99, 0.995), 1.0);
        assert!(estimated_ani(0.99, f64::NAN).is_nan());
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
                (1, Arc::from("GCF_1")),
                identity,
            )
        };
        let mut writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(Vec::new());
        row("Known", (0.99, 0.995)).write(&mut writer).unwrap();
        row("Unknown", (f64::NAN, f64::NAN))
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
            posteriors: vec![Vec::new(); group_of.len()],
            iterations: 0,
        };
        let first = FirstEm::new(&fragments, &group_of, &em, &resolved);
        let tests = CoverageTests {
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
    fn the_own_evidence_look_keeps_species_with_reads_of_their_own() {
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
            posteriors: vec![Vec::new(); group_of.len()],
            iterations: 0,
        };
        let tests = CoverageTests {
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
        let first = FirstEm::new(&fragments, &group_of, &em, &resolved);
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
        let data = data(vec![vec![row(2, 20, 0.97), row(100, 30, 0.99)]]);
        let policy = EvidencePolicy {
            mode: EvidenceMode::Single,
            identity: 0.95,
            mapped_len: 50,
            fraction: 0.6,
        };
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&data.read_alignments[0].alignments, policy);
        let (fragment, _) = fragment_of(0, candidates, policy).unwrap();
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
            posteriors: vec![Vec::new(); group_of.len()],
            iterations: 0,
        };
        let first = FirstEm::new(&fragments, &group_of, &em, &resolved);
        let tests = CoverageTests {
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
                (1, Arc::from(genome)),
                (0.99, 0.995),
            )
        };
        let genus_row = row("genus", 40.0, "GCF_2", &EvidenceCoverage::default());
        let mut block = TaxonBlock {
            taxid: 100,
            taxon: Some(row("genus", 60.0, "GCF_1", &coverage)),
            genomes: Vec::new(),
            breadth: 0.8,
            unit: true,
        };
        block.count(genus_row.clone());
        let merged = block.taxon.as_ref().unwrap();
        assert_eq!(merged.equivalent_mapcount, 100.0);
        assert!((merged.depth - 1.0).abs() < 1e-12);
        // The unit's coverage, genome and identity stay.
        assert_eq!(merged.breadth, 0.8);
        assert_eq!(&*merged.best_fit.genome, "GCF_1");
        // A unit whose reads all sit in genome rows takes the row.
        let mut block = TaxonBlock {
            taxid: 100,
            taxon: None,
            genomes: vec![row("strain", 60.0, "GCF_1", &coverage)],
            breadth: 0.8,
            unit: true,
        };
        block.count(genus_row);
        assert_eq!(block.taxon.as_ref().unwrap().equivalent_mapcount, 40.0);
    }

    #[test]
    fn rows_above_species_merge_with_the_species_at_their_genus_and_labelled_units() {
        let row = |taxid: u32, rank: &str, reads: f64| {
            profile_row(
                format!("Taxon {taxid}"),
                taxid,
                rank,
                (reads, reads / 100.0),
                &EvidenceCoverage::default(),
                (taxid, Arc::from(format!("GCF_{taxid}"))),
                (f64::NAN, f64::NAN),
            )
        };
        // A kept unit labelled at genus 100, and species 1.
        let mut blocks = vec![
            TaxonBlock {
                taxid: 100,
                taxon: Some(row(100, "genus", 50.0)),
                genomes: Vec::new(),
                breadth: 0.5,
                unit: true,
            },
            TaxonBlock {
                taxid: 1,
                taxon: Some(row(1, "species", 10.0)),
                genomes: Vec::new(),
                breadth: 0.5,
                unit: true,
            },
        ];
        // Rows above species: genus 100, genus 200 and family 300; species 5 and 6 reported at genus 200 and
        // species 7 at genus 400.
        let rows = vec![
            (100, row(100, "genus", 20.0)),
            (200, row(200, "genus", 30.0)),
            (300, row(300, "family", 40.0)),
        ];
        let parts = vec![
            (400, 7, row(400, "genus", 3.0)),
            (200, 6, row(200, "genus", 2.0)),
            (200, 5, row(200, "genus", 1.0)),
        ];
        add_rows_above_species(&mut blocks, rows, parts);
        let found = blocks
            .iter()
            .map(|b| {
                let taxon = b.taxon.as_ref().unwrap();
                (
                    b.taxid,
                    b.unit,
                    taxon.rank.clone(),
                    taxon.equivalent_mapcount,
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            found,
            [
                (100, true, "genus".to_string(), 70.0),
                (1, true, "species".to_string(), 10.0),
                (200, false, "genus".to_string(), 33.0),
                (300, false, "family".to_string(), 40.0),
                (400, false, "genus".to_string(), 3.0),
            ]
        );
        // The unit keeps its own genome; the genus row of 200 its own.
        assert_eq!(
            &*blocks[0].taxon.as_ref().unwrap().best_fit.genome,
            "GCF_100"
        );
        assert!((blocks[2].taxon.as_ref().unwrap().depth - 0.33).abs() < 1e-12);
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
        let lineages = Lineages::new(&taxonomy);
        let names = Names {
            taxonomy: &taxonomy,
            genomes: &[],
            lineages: &lineages,
        };
        assert_eq!(
            strain_label(&names, None, 561, "GCF_1"),
            ("Escherichia GCF_1".to_string(), 561)
        );
        assert_eq!(
            strain_label(&names, Some(83333), 562, "GCF_K12"),
            ("Escherichia coli K-12".to_string(), 83333)
        );
    }

    #[test]
    fn a_marker_genome_is_measured_gene_by_gene_and_needs_enough_genes() {
        // Genome 7: five marker genes (references 1-5) of 1,000 bases; genome 8
        // is a whole genome of 5,000 bases.
        let geometry = FragmentGeometry::new(FragmentLayout::Reads, 100, 100, 50, 0.5);
        let genes = [1_000u32; 5];
        let model = MarkerModel::from_gene_lengths(&[(7, &genes[..])], geometry);
        let tests = CoverageTests {
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
        };
        let fragments = (0..3)
            .map(|read| Fragment {
                read,
                candidates: Vec::new(),
                bases: 100.0,
            })
            .collect::<Vec<_>>();
        let group_of = [Some(0), Some(0), Some(1)];
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
            (&group_of, &em),
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
    fn reads_whose_best_match_is_kept_are_split_off_in_read_order() {
        // Reads 0 and 2 tie on species 2 and 3; read 1 matches species 4 best.
        let data = data(vec![
            vec![row(2, 20, 0.88), row(3, 30, 0.88)],
            vec![row(4, 40, 0.90), row(2, 20, 0.85)],
            vec![row(3, 30, 0.88), row(2, 20, 0.88)],
        ]);
        let policy = EvidencePolicy {
            mode: EvidenceMode::Single,
            identity: 0.80,
            mapped_len: 50,
            fraction: 0.6,
        };
        let split = split_best_kept(&data, policy, &[0, 1, 2], &HashSet::from([3]));
        assert_eq!(split, (vec![1], vec![0, 2]));
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
        };
        let weight = HashMap::from([(2, 30.0), (3, 10.0)]);
        let reported = HashMap::from([(2, 20), (3, 30)]);
        let credited = credited_reads(&data, policy, &[0, 1], (&weight, &reported), true);
        assert!((credited.reads[&2] - 0.75).abs() < 1e-12);
        assert!((credited.reads[&3] - 0.25).abs() < 1e-12);
        assert!((credited.bases[&2] - 75.0).abs() < 1e-9);
        assert!(!credited.reads.contains_key(&4));
        assert_eq!(credited.count, 1);
    }
}

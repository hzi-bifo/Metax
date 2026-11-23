//! Tables of the profile: the lines that open it and, with
//! `--save-decisions`, the decisions table (`<prefix>.decisions.tsv`). Each
//! opens with two lines starting with `#`: the run's metadata, in the style of
//! the mapping file's header, then the column names.
//!
//! The decisions table has one row per taxon the engine judged, at every rank,
//! kept or not:
//!
//! - `species`: one row per measured genome of every candidate species, its
//!   coverage from the reads that resolve to the species; a unit whose
//!   references are labelled only above species has these rows at its rank
//!   (`genus`, `family`, ...);
//! - `strain` (with `--strain`): one row per genome the genome level considered;
//! - `genus`: one row per kind of fragment (without a species-level match, or
//!   of species not kept) with at least 200 fragments unique to the genus, kept
//!   or not ([`super::higher`]), and one for the species reported at the genus
//!   when no such row's evidence is the genus row's; counted in the row of a
//!   kept unit labelled at that genus when there is one;
//! - `family`, `order`, `class`, `phylum`, `domain`: one row per taxon above
//!   genus that at least 200 fragments reach, kept or not; counted in the row
//!   of a kept unit labelled at that taxon when there is one.
//!
//! Species rows come first, by the species' abundance in the profile (largest
//! first, then taxid; species not in the profile last), each species' decided
//! genome first, then by reads (most first) and accession; then strain rows,
//! by species taxid, best-fit genome first, then by EM share (largest first)
//! and accession; then genus rows, by the genus's abundance and taxid, the
//! decided row first; then the rows above genus, by rank, abundance and taxid.
//! Columns, `NA` where they do not apply:
//!
//! - `rank`, `taxid`, `name`: the taxon (a strain row is named and numbered as
//!   its row in the profile);
//! - `parent_taxid`: the species (unit) of a strain, the genus of a species,
//!   the family of a genus; `NA` for a family and higher;
//! - `genome`: the genome the row's evidence is measured on;
//! - `decided`: species rows: the genome the species' decision rests on (the
//!   genome it is reported on when kept, otherwise its genome with most
//!   reads); strain rows: each genome reported as a strain row of the
//!   profile; genus rows: the row whose evidence the genus row of the profile
//!   shows (otherwise the one with most unique reads); rows above genus:
//!   `true`;
//! - `reads`: observed reads behind the evidence, a read pair counting once,
//!   not scaled like column 4 of the profile: species rows: reads resolved to
//!   the species and placed on the genome; strain rows: reads the genome EM
//!   resolves to the genome; genus and higher rows: the fragments placed on the
//!   genome (for species reported at the genus, their reads);
//! - `depth` to `chunk_pvalue`: the coverage on `genome` and the chunk-breadth
//!   cutoff of the tests (as in the profile; ratios to four decimals);
//! - `candidate_reads`, `won_reads`: reads listing the species and those
//!   resolved to it; genus rows: reads of the kind whose tie set lists the
//!   genus and those unique to it; rows above genus: reads that reach the
//!   taxon; `second_look_reads`: reads the species took on the second look;
//!   `own_reads`: reads it won that no kept species lists;
//! - `em_share`, `strain_identity_reads`: a strain's share of its species' reads
//!   in the genome EM, and its resolved reads at strain-level identity;
//! - `window_spread`: share of the genome's windows that the reads of a genus
//!   or higher row reach;
//! - `genome_tests`: the row's own tests, `passed` or each failed test with its
//!   value and bound (on a genus row, after the kind of reads it judges);
//! - `decision`: `kept` when the taxon has rows of its own in the profile,
//!   otherwise `removed`; the same on every row of a taxon;
//! - `reason`: the taxon's outcome in words, the same on every row of a
//!   species (`removed by --min-reads: …`, `removed by --min-breadth: …` or
//!   `removed by --host: …` for a taxon those options dropped; a labelled unit
//!   says that its references are labelled only at its rank, and, when named
//!   units of its taxon explain its reads, which ones; a genus or higher taxon
//!   gives its reads, their mean identity and chunk ratio or window spread
//!   against the bounds, and where the reads of a failed taxon went);
//! - `reported_as`: the rank of the profile rows holding the taxon's reads
//!   (the unit's rank for its row and its strain rows), `none` when none;
//! - `final_reads`, `final_abundance`: columns 4 and 6 of the taxon's profile
//!   row; for a species its species row plus its strain rows; on species rows
//!   only on the decided genome; `NA` when the taxon has no row of its own.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;

use anyhow::{Context, Result};

use super::genomes::{SpeciesGenomes, StrainDiagnostic};
use super::higher::{
    self, GenusCandidate, HigherLevels, HigherRow, Rank, TaxonCandidate, GENUS_MIN_CHUNK_RATIO,
};
use super::{
    strain_label, wins_enough, within, Dropped, FragmentCounts, GenomeEvidence, SpeciesTests,
    TaxonBlock, MIN_WON_SHARE,
};
use crate::db::DatabaseMetadata;
use crate::profiler::evidence::{EvidenceCoverage, STANDARD_RATIO_UPPER};
use crate::profiler::output::format_five;
use crate::profiler::{ProfileRank, ProfilerConfig};
use crate::sequencer::Sequencer;
use crate::taxonomy::Taxonomy;

/// Columns that `--host` appends to every row of the profile.
pub(super) const HOST_COLUMNS: [&str; 3] = ["host_names", "host_taxids", "diseases"];

/// Columns of the profile. Every row has all of them, the
/// unclassified row included.
pub(super) const PROFILE_COLUMNS: [&str; 17] = [
    "name",
    "taxid",
    "rank",
    "reads",
    "depth",
    "abundance",
    "breadth",
    "expected_breadth",
    "breadth_pvalue",
    "fixed_chunk_breadth",
    "flex_chunk_breadth",
    "expected_flex_chunk_breadth",
    "chunk_pvalue",
    "evidence_depth",
    "best_fit_genome",
    "mean_identity",
    "estimated_ani",
];

/// Columns of the decisions table (see the module documentation).
const DECISION_COLUMNS: [&str; 31] = [
    "rank",
    "taxid",
    "name",
    "parent_taxid",
    "genome",
    "decided",
    "reads",
    "depth",
    "breadth",
    "expected_breadth",
    "oebr",
    "fixed_chunk_breadth",
    "flex_chunk_breadth",
    "expected_flex_chunk_breadth",
    "coebr",
    "chunk_cutoff",
    "breadth_pvalue",
    "chunk_pvalue",
    "candidate_reads",
    "won_reads",
    "second_look_reads",
    "own_reads",
    "em_share",
    "strain_identity_reads",
    "window_spread",
    "genome_tests",
    "decision",
    "reason",
    "reported_as",
    "final_reads",
    "final_abundance",
];

/// The run named on the first line of each table.
pub(super) struct TableHeader {
    /// Fingerprint of the database, as in the mapping file's header; `None`
    /// without database metadata.
    pub(super) fingerprint: Option<u128>,
    pub(super) sequencer: Sequencer,
    pub(super) paired: bool,
    /// `--host`, named only when given.
    pub(super) host: Option<u32>,
}

impl TableHeader {
    pub(super) fn of(cfg: &ProfilerConfig, host: Option<u32>) -> Self {
        Self {
            fingerprint: cfg
                .db_metadata
                .as_deref()
                .map(DatabaseMetadata::fingerprint),
            sequencer: cfg.sequencer,
            paired: cfg.is_paired,
            host,
        }
    }

    /// The two header lines of table `table`: `#metax_<table>` with the run's
    /// metadata, then `#` and the `columns`.
    pub(super) fn lines(&self, table: &str, columns: &[&str]) -> String {
        let fingerprint = self
            .fingerprint
            .map_or_else(|| "NA".to_string(), |f| format!("{f:032x}"));
        let host = self
            .host
            .map_or_else(String::new, |host| format!("\thost={host}"));
        format!(
            "#metax_{table}\tversion={}\tengine=level\tdatabase_fingerprint={fingerprint}\tsequencer={}\tpaired={}{host}\n#{}\n",
            env!("CARGO_PKG_VERSION"),
            self.sequencer,
            self.paired,
            columns.join("\t")
        )
    }
}

/// Names of taxa and genomes (by assembly id).
pub(super) struct Names<'a> {
    pub(super) taxonomy: &'a Taxonomy,
    pub(super) genomes: &'a [Arc<str>],
}

impl Names<'_> {
    fn taxon(&self, taxid: u32) -> String {
        self.taxonomy
            .get_name(taxid)
            .unwrap_or("unknown")
            .to_string()
    }

    fn genome(&self, assembly_id: u32) -> &str {
        &self.genomes[assembly_id as usize]
    }

    /// The ancestor of `taxid` at `rank`.
    fn ancestor(&self, taxid: u32, rank: &str) -> Option<u32> {
        self.taxonomy
            .get_parents(taxid)
            .into_iter()
            .find(|(_, _, r)| r == rank)
            .map(|(ancestor, _, _)| ancestor)
    }
}

/// What each taxon holds in the profile: reads and abundance (columns 4 and 6).
pub(super) struct FinalRows {
    /// Of each species with rows of its own, summed over its species row and
    /// genome rows in profile order.
    species: HashMap<u32, [f64; 2]>,
    /// Of each genome (strain) row, by genome.
    strains: HashMap<Arc<str>, [f64; 2]>,
    /// Of each row above species, by taxid.
    taxa: HashMap<u32, [f64; 2]>,
}

impl FinalRows {
    /// From the profile's taxon blocks, in profile order.
    pub(super) fn of(blocks: &[TaxonBlock]) -> Self {
        let mut finals = Self {
            species: HashMap::new(),
            strains: HashMap::new(),
            taxa: HashMap::new(),
        };
        let amount =
            |row: &crate::profiler::output::ProfileRow| [row.equivalent_mapcount, row.abundance];
        for block in blocks {
            if !block.unit {
                if let Some(row) = &block.taxon {
                    finals.taxa.insert(block.taxid, amount(row));
                }
                continue;
            }
            let mut total = [0.0; 2];
            for row in block.taxon.iter().chain(&block.genomes) {
                total[0] += row.equivalent_mapcount;
                total[1] += row.abundance;
            }
            for row in &block.genomes {
                if let Some(best_fit) = &row.best_fit {
                    finals
                        .strains
                        .insert(Arc::clone(&best_fit.genome), amount(row));
                }
            }
            finals.species.insert(block.taxid, total);
        }
        finals
    }
}

/// Write a table to `path`.
pub(super) fn write_file(
    path: &str,
    write: impl FnOnce(&mut dyn Write) -> std::io::Result<()>,
) -> Result<()> {
    let mut out =
        BufWriter::new(File::create(path).with_context(|| format!("failed to create {path}"))?);
    write(&mut out)
        .and_then(|()| out.flush())
        .with_context(|| format!("failed to write {path}"))
}

/// What became of a measured species.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum SpeciesDecision {
    /// Kept by the first judgment.
    KeptTests,
    /// Kept on the second look.
    KeptSecondLook,
    /// Kept on the last look, on its own reads.
    KeptOwnEvidence,
    /// Kept for even coverage.
    KeptEvenCoverage,
    /// Kept on the last look, reported at its genus as its reads tie with
    /// relatives along its genome.
    AtGenusTie,
    /// Reported at its genus by the identity rules for relatives.
    AtGenusIdentity,
    /// Left out by the identity rules for relatives.
    LeftOut,
    /// Kept, then removed by the prophage check.
    RemovedProphage,
    /// A unit labelled above species, kept, then removed because the kept
    /// named units of its taxon explain its reads.
    RemovedNamedFirst,
    /// Never kept.
    Removed,
}

impl SpeciesDecision {
    fn kept(self) -> bool {
        matches!(
            self,
            Self::KeptTests | Self::KeptSecondLook | Self::KeptOwnEvidence | Self::KeptEvenCoverage
        )
    }

    fn at_genus(self) -> bool {
        matches!(self, Self::AtGenusTie | Self::AtGenusIdentity)
    }

    /// The outcome in words; `tests` are those of the species' decided genome.
    fn reason(self, tests: &str) -> String {
        match self {
            Self::KeptTests => "passed amount, evenness and won share".to_string(),
            Self::KeptSecondLook => {
                "kept on the second look: a removed sister no longer splits its reads".to_string()
            }
            Self::KeptOwnEvidence => {
                "kept on its own reads, which no kept species explains".to_string()
            }
            Self::KeptEvenCoverage => {
                "kept for even coverage: its reads cover its genome evenly".to_string()
            }
            Self::AtGenusTie | Self::AtGenusIdentity => format!("reported at genus: {}", self.why_at_genus()),
            Self::LeftOut => "left out: a relative with more identical reads is kept".to_string(),
            Self::RemovedProphage => "removed by the prophage check: its own reads do not cover the genome outside the part a kept host explains".to_string(),
            Self::RemovedNamedFirst => "removed: its reads are explained by named units of its taxon".to_string(),
            Self::Removed => format!("removed: {tests}"),
        }
    }

    fn why_at_genus(self) -> &'static str {
        if self == Self::AtGenusTie {
            "its reads tie with relatives along its genome"
        } else {
            "its reads are less identical than the sample's typical species"
        }
    }
}

/// The species that the looks after the first judgment kept or removed.
pub(super) struct Looks<'a> {
    /// Kept on the second look, with the fragments each took.
    pub(super) second_look: &'a HashMap<u32, Vec<usize>>,
    /// Kept on the last look, with whether each is reported at its genus.
    pub(super) own_evidence: &'a HashMap<u32, bool>,
    /// Kept viruses removed by the prophage check.
    pub(super) prophage: &'a [u32],
    /// Kept units labelled above species removed as the kept named units of
    /// their taxon explain their reads.
    pub(super) named_first: &'a HashMap<u32, Vec<u32>>,
    /// Kept for even coverage.
    pub(super) even_coverage: &'a [u32],
}

/// The decision of each species with a measured genome once the looks are
/// done: how a species with a best-fit genome in `kept` was kept, and why one
/// without was not. The identity rules and the genus reports then turn some
/// kept species into `LeftOut` or `AtGenus*`.
pub(super) fn judged_species(
    genomes: &HashMap<u32, GenomeEvidence>,
    kept: &HashMap<u32, u32>,
    looks: &Looks,
) -> HashMap<u32, SpeciesDecision> {
    let prophage: HashSet<u32> = looks.prophage.iter().copied().collect();
    let even: HashSet<u32> = looks.even_coverage.iter().copied().collect();
    let decision = |taxid: u32| {
        if !kept.contains_key(&taxid) {
            if prophage.contains(&taxid) {
                SpeciesDecision::RemovedProphage
            } else if looks.named_first.contains_key(&taxid) {
                SpeciesDecision::RemovedNamedFirst
            } else {
                SpeciesDecision::Removed
            }
        } else if even.contains(&taxid) {
            SpeciesDecision::KeptEvenCoverage
        } else if looks.own_evidence.contains_key(&taxid) {
            SpeciesDecision::KeptOwnEvidence
        } else if looks.second_look.contains_key(&taxid) {
            SpeciesDecision::KeptSecondLook
        } else {
            SpeciesDecision::KeptTests
        }
    };
    let mut decisions = HashMap::new();
    for genome in genomes.values() {
        decisions
            .entry(genome.taxid)
            .or_insert_with(|| decision(genome.taxid));
    }
    decisions
}

/// The genome level, with `--strain`.
pub(super) struct StrainLevel<'a> {
    pub(super) rows: &'a [StrainDiagnostic],
    /// The genomes the genome level reports for each species.
    pub(super) species: &'a HashMap<u32, SpeciesGenomes>,
    /// Its coverage tests (chunk-breadth cutoff and ratio bounds).
    pub(super) tests: SpeciesTests,
    pub(super) strain_taxids: &'a HashMap<Arc<str>, u32>,
}

/// What the decisions table is written from.
pub(super) struct Decisions<'a> {
    pub(super) genomes: &'a HashMap<u32, GenomeEvidence>,
    pub(super) counts: &'a FragmentCounts,
    pub(super) tests: &'a SpeciesTests,
    pub(super) outcomes: &'a HashMap<u32, SpeciesDecision>,
    /// The genome each kept species is reported on.
    pub(super) reported: &'a HashMap<u32, u32>,
    pub(super) strains: Option<StrainLevel<'a>>,
    /// Species reported at their genus: (genus, species, reads).
    pub(super) genus_members: &'a [(u32, u32, f64)],
    /// The levels above species, when they ran.
    pub(super) higher: Option<&'a HigherLevels>,
    /// Taxa and genome rows dropped by `--host`, `--min-reads` or
    /// `--min-breadth`.
    pub(super) dropped: &'a Dropped,
    /// Units labelled above species, with their rank.
    pub(super) labels: &'a HashMap<u32, ProfileRank>,
    /// Labelled units removed as named units of their taxon explain their
    /// reads, with those units.
    pub(super) named_first: &'a HashMap<u32, Vec<u32>>,
}

/// What a genus is made of: its candidates above species, its row there, and
/// the species reported at the genus with their reads.
#[derive(Default)]
struct GenusParts<'a> {
    candidates: Vec<&'a GenusCandidate>,
    row: Option<&'a HigherRow>,
    members: Vec<(u32, f64)>,
}

/// One row of the decisions table.
struct DecisionRow<'n> {
    rank: &'static str,
    taxid: u32,
    name: String,
    parent: Option<u32>,
    genome: Option<&'n str>,
    decided: bool,
    reads: String,
    coverage: Option<EvidenceCoverage>,
    chunk_cutoff: Option<f64>,
    /// Candidate, won, second-look and own reads.
    counts: [Option<f64>; 4],
    em_share: Option<f64>,
    strain_identity_reads: Option<u64>,
    window_spread: Option<f64>,
    genome_tests: String,
    kept: bool,
    reason: String,
    reported_as: &'static str,
    finals: Option<[f64; 2]>,
}

fn na() -> String {
    "NA".to_string()
}

fn four_decimals(value: f64) -> String {
    format!("{value:.4}")
}

fn pvalue(value: Option<f64>) -> String {
    value.map_or_else(na, |v| format!("{v:.3e}"))
}

impl DecisionRow<'_> {
    fn fields(&self) -> [String; DECISION_COLUMNS.len()] {
        let coverage = self.coverage.as_ref();
        let measured = |value: fn(&EvidenceCoverage) -> String| coverage.map_or_else(na, value);
        let count = |value: Option<f64>| value.map_or_else(na, |v| format!("{v:.0}"));
        let finals = |i: usize| self.finals.map_or_else(na, |f| format_five(f[i]));
        [
            self.rank.to_string(),
            self.taxid.to_string(),
            self.name.clone(),
            self.parent.map_or_else(na, |t| t.to_string()),
            self.genome.map_or_else(na, str::to_string),
            self.decided.to_string(),
            self.reads.clone(),
            measured(|c| format_five(c.depth)),
            measured(|c| format_five(c.breadth)),
            measured(|c| format_five(c.expected_breadth)),
            measured(|c| four_decimals(c.oebr)),
            measured(|c| four_decimals(c.fixed_chunk)),
            measured(|c| format_five(c.flex_chunk)),
            measured(|c| format_five(c.expected_flex_chunk)),
            measured(|c| four_decimals(c.coebr)),
            self.chunk_cutoff.map_or_else(na, four_decimals),
            pvalue(coverage.and_then(|c| c.breadth_pvalue)),
            pvalue(coverage.and_then(|c| c.chunk_pvalue)),
            count(self.counts[0]),
            count(self.counts[1]),
            count(self.counts[2]),
            count(self.counts[3]),
            self.em_share.map_or_else(na, four_decimals),
            self.strain_identity_reads
                .map_or_else(na, |r| r.to_string()),
            self.window_spread.map_or_else(na, four_decimals),
            self.genome_tests.clone(),
            if self.kept { "kept" } else { "removed" }.to_string(),
            self.reason.clone(),
            self.reported_as.to_string(),
            finals(0),
            finals(1),
        ]
    }
}

/// `passed`, or each failed test.
fn tests_in_words(failed: &[String]) -> String {
    if failed.is_empty() {
        "passed".to_string()
    } else {
        format!("failed {}", failed.join("; failed "))
    }
}

/// The coverage tests of `tests` that `coverage` fails, in words.
fn coverage_failures(tests: &SpeciesTests, coverage: &EvidenceCoverage) -> Vec<String> {
    let mut failed = Vec::new();
    if !tests.amount_ok(coverage) {
        failed.push(format!(
            "amount: fixed chunk {:.4} <= cutoff {:.4}",
            coverage.fixed_chunk, tests.chunk_cutoff
        ));
    }
    if !tests.marker_genes_ok(coverage) {
        failed.push(format!(
            "marker genes: {} < {}",
            coverage.marker_genes_hit.unwrap_or(0),
            tests.min_marker_genes
        ));
    }
    for (name, ratio, lower) in [
        ("OEBR", coverage.oebr, tests.ratio_bounds.0),
        ("COEBR", coverage.coebr, tests.ratio_bounds.1),
    ] {
        if !within(ratio, lower) {
            failed.push(if ratio > STANDARD_RATIO_UPPER {
                format!("evenness: {name} {ratio:.4} > {STANDARD_RATIO_UPPER}")
            } else {
                format!("evenness: {name} {ratio:.4} < {lower}")
            });
        }
    }
    failed
}

/// Larger final abundance first; taxa without one last.
fn by_abundance(a: Option<[f64; 2]>, b: Option<[f64; 2]>) -> Ordering {
    match (a, b) {
        (Some(a), Some(b)) => b[1].total_cmp(&a[1]),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => Ordering::Equal,
    }
}

impl Decisions<'_> {
    /// Write the table.
    pub(super) fn write(
        &self,
        out: &mut dyn Write,
        header: &TableHeader,
        names: &Names,
        finals: &FinalRows,
    ) -> std::io::Result<()> {
        out.write_all(header.lines("decisions", &DECISION_COLUMNS).as_bytes())?;
        for row in self.rows(names, finals) {
            writeln!(out, "{}", row.fields().join("\t"))?;
        }
        Ok(())
    }

    fn rows<'n>(&self, names: &'n Names, finals: &FinalRows) -> Vec<DecisionRow<'n>> {
        let mut rows = self.species_rows(names, finals);
        rows.extend(self.strain_rows(names, finals));
        let mut genera = self.genus_rows(names, finals);
        genera.sort_by(|a, b| {
            let total = |rows: &[DecisionRow]| rows.iter().find_map(|row| row.finals);
            by_abundance(total(a), total(b)).then(a[0].taxid.cmp(&b[0].taxid))
        });
        rows.extend(genera.into_iter().flatten());
        let mut taxa = self.taxon_rows(names, finals);
        taxa.sort_by(|a, b| {
            a.0.cmp(&b.0)
                .then(by_abundance(a.1.finals, b.1.finals))
                .then(a.1.taxid.cmp(&b.1.taxid))
        });
        rows.extend(taxa.into_iter().map(|(_, row)| row));
        rows
    }

    /// The tests of genome `genome` of its species, in words.
    fn species_genome_tests(&self, genome: &GenomeEvidence) -> String {
        let coverage = genome.coverage.as_ref().expect("measured before");
        let mut failed = coverage_failures(self.tests, coverage);
        let count = |map: &HashMap<u32, u64>| map.get(&genome.taxid).copied().unwrap_or(0);
        let (won, candidates) = (count(&self.counts.won), count(&self.counts.candidates));
        if !wins_enough(won, candidates) {
            failed.push(format!(
                "won share: {:.4} < {MIN_WON_SHARE}",
                won as f64 / candidates.max(1) as f64
            ));
        }
        tests_in_words(&failed)
    }

    fn species_rows<'n>(&self, names: &'n Names, finals: &FinalRows) -> Vec<DecisionRow<'n>> {
        let genus_of: HashMap<u32, u32> = self
            .genus_members
            .iter()
            .map(|&(genus, taxid, _)| (taxid, genus))
            .collect();
        let mut by_species: HashMap<u32, Vec<u32>> = HashMap::new();
        for (&id, genome) in self.genomes {
            by_species.entry(genome.taxid).or_default().push(id);
        }
        let mut species = by_species.into_iter().collect::<Vec<_>>();
        species.sort_by(|a, b| {
            by_abundance(
                finals.species.get(&a.0).copied(),
                finals.species.get(&b.0).copied(),
            )
            .then(a.0.cmp(&b.0))
        });
        let mut rows = Vec::new();
        for (taxid, mut ids) in species {
            // Most reads first, then accession; the decided genome leads.
            ids.sort_by(|&a, &b| {
                self.genomes[&b]
                    .assembly
                    .count
                    .cmp(&self.genomes[&a].assembly.count)
                    .then_with(|| names.genome(a).cmp(names.genome(b)))
            });
            let decided = self.reported.get(&taxid).copied().unwrap_or(ids[0]);
            if let Some(position) = ids.iter().position(|&id| id == decided) {
                ids[..=position].rotate_right(1);
            }
            let decision = self.outcomes[&taxid];
            let below = self.dropped.taxa.get(&taxid);
            // A species reported at its genus whose genus row was dropped.
            let genus_dropped = genus_of
                .get(&taxid)
                .and_then(|genus| self.dropped.taxa.get(genus));
            let label = self.labels.get(&taxid).copied();
            let rank = label.map_or("species", ProfileRank::as_str);
            let mut reason = match below {
                Some(below) => below.reason(),
                None => {
                    let reason =
                        decision.reason(&self.species_genome_tests(&self.genomes[&decided]));
                    match genus_dropped {
                        Some(genus) => format!("{reason}; its genus row {}", genus.reason()),
                        None => reason,
                    }
                }
            };
            if let Some(named) = self.named_first.get(&taxid) {
                let named = named
                    .iter()
                    .map(|&unit| format!("{} ({unit})", names.taxon(unit)))
                    .collect::<Vec<_>>();
                reason = format!("{reason}: {}", named.join(", "));
            }
            if label.is_some() {
                reason = format!("{reason}; its references are labelled only at {rank}");
            }
            let total = finals.species.get(&taxid).copied();
            let reported_as = match total {
                Some(_) => rank,
                None if decision.at_genus() && genus_dropped.is_none() => "genus",
                None => "none",
            };
            let parent = match label {
                None => names.ancestor(taxid, "genus"),
                Some(ProfileRank::Genus) => names.ancestor(taxid, "family"),
                Some(_) => None,
            };
            let count =
                |map: &HashMap<u32, u64>| Some(map.get(&taxid).copied().unwrap_or(0) as f64);
            for id in ids {
                let genome = &self.genomes[&id];
                rows.push(DecisionRow {
                    rank,
                    taxid,
                    name: names.taxon(taxid),
                    parent,
                    genome: Some(names.genome(id)),
                    decided: id == decided,
                    reads: genome.assembly.count.to_string(),
                    coverage: genome.coverage,
                    chunk_cutoff: Some(self.tests.chunk_cutoff),
                    counts: [
                        count(&self.counts.candidates),
                        count(&self.counts.won),
                        count(&self.counts.second_look),
                        count(&self.counts.own),
                    ],
                    em_share: None,
                    strain_identity_reads: None,
                    window_spread: None,
                    genome_tests: self.species_genome_tests(genome),
                    kept: decision.kept() && below.is_none(),
                    reason: reason.clone(),
                    reported_as,
                    finals: total.filter(|_| id == decided),
                });
            }
        }
        rows
    }

    fn strain_rows<'n>(&self, names: &'n Names, finals: &FinalRows) -> Vec<DecisionRow<'n>> {
        let Some(level) = &self.strains else {
            return Vec::new();
        };
        let mut diagnostics = level.rows.iter().collect::<Vec<_>>();
        diagnostics.sort_by(|a, b| {
            a.taxid
                .cmp(&b.taxid)
                .then(b.best_fit.cmp(&a.best_fit))
                .then(b.em_share.total_cmp(&a.em_share))
                .then_with(|| names.genome(a.assembly_id).cmp(names.genome(b.assembly_id)))
        });
        let mut rows = Vec::new();
        for row in diagnostics {
            let accession = names.genome(row.assembly_id);
            let (name, taxid) =
                strain_label(names.taxonomy, level.strain_taxids, row.taxid, accession);
            let strain_row = finals.strains.get(accession).copied();
            let passed = level.species.get(&row.taxid).is_some_and(|species| {
                species
                    .strains
                    .iter()
                    .any(|strain| strain.assembly_id == row.assembly_id)
            });
            let (genome_tests, reason) = if strain_row.is_some() {
                (
                    "passed".to_string(),
                    "same strain: its reads match it within sequencing error and cover it evenly"
                        .to_string(),
                )
            } else if passed {
                (
                    "passed".to_string(),
                    match (
                        self.dropped.genomes.get(accession),
                        self.dropped.taxa.get(&row.taxid),
                    ) {
                        (Some(removal), _) => removal.reason(),
                        (None, Some(removal)) => {
                            format!("species not reported: {}", removal.reason())
                        }
                        (None, None) => {
                            "species not reported: its species is left out or reported at its genus"
                                .to_string()
                        }
                    },
                )
            } else if row.kept {
                (
                    "failed identity: its reads do not match it within sequencing error"
                        .to_string(),
                    "relative only: its reads match it as a relative's, not its own strain's"
                        .to_string(),
                )
            } else if let Some(coverage) = &row.coverage {
                (
                    tests_in_words(&coverage_failures(&level.tests, coverage)),
                    "coverage failed: its reads at strain-level identity do not cover it evenly enough".to_string(),
                )
            } else if row.resolved > 0 && row.gated == 0 {
                (
                    "failed identity: no read matches it at strain-level identity".to_string(),
                    "low identity: its reads match it below strain-level identity".to_string(),
                )
            } else {
                (
                    "not measured: no read resolves to it".to_string(),
                    "not measured: no read resolves to it at strain-level identity".to_string(),
                )
            };
            rows.push(DecisionRow {
                rank: "strain",
                taxid,
                name,
                parent: Some(row.taxid),
                genome: Some(accession),
                decided: strain_row.is_some(),
                reads: row.resolved.to_string(),
                coverage: row.coverage,
                chunk_cutoff: Some(level.tests.chunk_cutoff),
                counts: [None; 4],
                em_share: Some(row.em_share),
                strain_identity_reads: Some(row.gated),
                window_spread: None,
                genome_tests,
                kept: strain_row.is_some(),
                reason,
                reported_as: if strain_row.is_some() {
                    "strain"
                } else {
                    "none"
                },
                finals: strain_row,
            });
        }
        rows
    }

    /// The rows of each genus judged above species or holding species
    /// reported at it, one group per genus, the decided row first: one row per
    /// kind of fragment it is a candidate with, and one for the species
    /// reported at it when no candidate's evidence is its row's.
    fn genus_rows<'n>(&self, names: &'n Names, finals: &FinalRows) -> Vec<Vec<DecisionRow<'n>>> {
        let mut genera: BTreeMap<u32, GenusParts> = BTreeMap::new();
        if let Some(higher) = self.higher {
            for candidate in &higher.genera {
                genera
                    .entry(candidate.genus)
                    .or_default()
                    .candidates
                    .push(candidate);
            }
            for row in higher.rows.iter().filter(|row| row.rank == Rank::Genus) {
                genera.entry(row.taxid).or_default().row = Some(row);
            }
        }
        for &(genus, taxid, reads) in self.genus_members {
            genera
                .entry(genus)
                .or_default()
                .members
                .push((taxid, reads));
        }
        let mut groups = Vec::new();
        for (genus, mut parts) in genera {
            parts.members.sort_by_key(|member| member.0);
            let in_profile = parts.row.is_some() || !parts.members.is_empty();
            let below = self.dropped.taxa.get(&genus).filter(|_| in_profile);
            let kept = in_profile && below.is_none();
            let made_of = self.genus_parts_in_words(genus, &parts, names);
            let reason = match below {
                Some(below) => format!("{}; made of {made_of}", below.reason()),
                None if kept && finals.species.contains_key(&genus) => format!(
                    "kept: {made_of}; counted in the row of the unit labelled at this genus"
                ),
                None if kept => format!("kept: {made_of}"),
                None => format!("removed: {made_of}"),
            };
            let reported_as = if kept {
                "genus"
            } else if below.is_some() {
                "none"
            } else {
                self.higher
                    .and_then(|higher| higher.end_of_genus(genus))
                    .map_or("none", |end| end.rank.as_str())
            };
            // The candidate whose evidence the genus row shows; without one,
            // the species reported at the genus, or else the candidate with
            // most unique reads.
            let decided = parts
                .candidates
                .iter()
                .position(|candidate| candidate.shown)
                .or_else(|| {
                    parts.members.is_empty().then(|| {
                        let most = parts.candidates.iter().map(|c| c.unique).max();
                        parts
                            .candidates
                            .iter()
                            .position(|c| Some(c.unique) == most)
                            .expect("a genus without members has candidates")
                    })
                });
            let row = |genome: u32, decided: bool| DecisionRow {
                rank: "genus",
                taxid: genus,
                name: names.taxon(genus),
                parent: names.ancestor(genus, "family"),
                genome: Some(names.genome(genome)),
                decided,
                reads: String::new(),
                coverage: None,
                chunk_cutoff: None,
                counts: [None; 4],
                em_share: None,
                strain_identity_reads: None,
                window_spread: None,
                genome_tests: String::new(),
                kept,
                reason: reason.clone(),
                reported_as,
                finals: finals.taxa.get(&genus).copied().filter(|_| decided),
            };
            let mut rows = Vec::new();
            if decided.is_none() {
                let id = self.reported[&parts.members[0].0];
                let reads = parts.members.iter().map(|member| member.1).sum::<f64>();
                rows.push(DecisionRow {
                    reads: format_five(reads),
                    coverage: Some(self.genomes[&id].coverage.expect("measured before")),
                    genome_tests: "no genus tests: made of species reported at their genus"
                        .to_string(),
                    ..row(id, true)
                });
            }
            for (index, candidate) in parts.candidates.iter().enumerate() {
                let measured = &candidate.measured;
                rows.push(DecisionRow {
                    reads: measured.fragments.to_string(),
                    coverage: Some(measured.coverage),
                    counts: [
                        Some((candidate.unique + candidate.tied) as f64),
                        Some(candidate.unique as f64),
                        None,
                        None,
                    ],
                    window_spread: measured.spread,
                    genome_tests: format!(
                        "reads {}: {}",
                        candidate.kind.words(),
                        tests_in_words(&failures_in_words(&candidate.failures))
                    ),
                    ..row(measured.genome, decided == Some(index))
                });
            }
            rows.sort_by_key(|row| !row.decided);
            groups.push(rows);
        }
        groups
    }

    /// What a genus is made of, in words: each kind of its unique reads with
    /// their identity and chunk ratio, the reads tied with other genera that
    /// its row took, and the species reported at the genus.
    fn genus_parts_in_words(&self, genus: u32, parts: &GenusParts, names: &Names) -> String {
        let mut words = Vec::new();
        let floor = self.higher.map_or(f64::NAN, |higher| higher.identity_floor);
        // The kinds of reads whose evidence gave the genus its row.
        let passed = parts
            .candidates
            .iter()
            .filter(|candidate| candidate.failures.is_empty())
            .map(|candidate| candidate.kind.words())
            .collect::<Vec<_>>()
            .join(" and ");
        for candidate in &parts.candidates {
            let failed = |test: fn(&higher::Failure) -> bool| {
                if candidate.failures.iter().any(test) {
                    "<"
                } else {
                    ">="
                }
            };
            let mut tests = format!(
                "{} reads {} unique to it, identity {:.4} {} {floor:.4}, chunk ratio {:.4} {} \
                 {GENUS_MIN_CHUNK_RATIO}",
                candidate.unique,
                candidate.kind.words(),
                candidate.identity,
                failed(|f| matches!(f, higher::Failure::Identity { .. })),
                candidate.measured.coverage.coebr,
                failed(|f| matches!(f, higher::Failure::ChunkRatio(_))),
            );
            if !candidate.failures.is_empty() {
                if parts.row.is_some() {
                    tests.push_str(&format!("; they joined the genus row from reads {passed}"));
                } else {
                    let end = self.higher.and_then(|higher| higher.end_of_genus(genus));
                    tests.push_str(&format!(
                        "; they {}",
                        destination_in_words("go to", end, names)
                    ));
                }
            }
            words.push(tests);
        }
        if let Some(row) = parts.row {
            // Unique reads of a kind too few to be judged join the row too.
            for (kind, &unique) in [higher::Kind::BelowSpecies, higher::Kind::NotKept]
                .iter()
                .zip(&row.unique)
            {
                let judged = parts.candidates.iter().any(|c| c.kind == *kind);
                if unique > 0 && !judged {
                    words.push(format!(
                        "{unique} reads {} unique to it, too few to judge, joined the genus row",
                        kind.words()
                    ));
                }
            }
            let tied = row.reads - row.unique.iter().sum::<usize>() as f64;
            if tied > 1e-9 {
                words.push(format!("{tied:.1} reads tied with other genera"));
            }
        }
        if !parts.members.is_empty() {
            let species = parts
                .members
                .iter()
                .map(|&(taxid, reads)| {
                    format!(
                        "{taxid} ({reads:.0} reads; {})",
                        self.outcomes[&taxid].why_at_genus()
                    )
                })
                .collect::<Vec<_>>();
            words.push(format!("species reported at genus: {}", species.join(", ")));
        }
        words.join(" + ")
    }

    /// One row per taxon above genus that at least 200 fragments reach, kept
    /// or not, with its rank.
    fn taxon_rows<'n>(&self, names: &'n Names, finals: &FinalRows) -> Vec<(Rank, DecisionRow<'n>)> {
        let Some(higher) = self.higher else {
            return Vec::new();
        };
        higher
            .taxa
            .iter()
            .map(|candidate| {
                let taxid = candidate.taxid;
                let rank = candidate.rank.as_str();
                let in_profile = higher.row(taxid).is_some();
                let below = self.dropped.taxa.get(&taxid).filter(|_| in_profile);
                let kept = in_profile && below.is_none();
                let made_of = taxon_parts_in_words(candidate, names);
                let end = higher.end_of(taxid);
                let reason = match below {
                    Some(below) => below.reason(),
                    None if kept => {
                        let unit = if finals.species.contains_key(&taxid) {
                            format!("; counted in the row of the unit labelled at this {rank}")
                        } else {
                            String::new()
                        };
                        format!(
                            "kept: {made_of}; {}{unit}",
                            higher::taxon_evidence(&candidate.measured)
                        )
                    }
                    None => format!(
                        "removed: {made_of}; failed {}; they {}",
                        failures_in_words(&candidate.failures).join(" and "),
                        destination_in_words("move up to", end, names)
                    ),
                };
                let reported_as = match (kept, below) {
                    (true, _) => rank,
                    (false, Some(_)) => "none",
                    (false, None) => end.map_or("none", |end| end.rank.as_str()),
                };
                let measured = &candidate.measured;
                (
                    candidate.rank,
                    DecisionRow {
                        rank,
                        taxid,
                        name: names.taxon(taxid),
                        parent: None,
                        genome: Some(names.genome(measured.genome)),
                        decided: true,
                        reads: measured.fragments.to_string(),
                        coverage: Some(measured.coverage),
                        chunk_cutoff: None,
                        counts: [Some(candidate.fragments() as f64), None, None, None],
                        em_share: None,
                        strain_identity_reads: None,
                        window_spread: measured.spread,
                        genome_tests: tests_in_words(&failures_in_words(&candidate.failures)),
                        kept,
                        reason,
                        reported_as,
                        finals: finals.taxa.get(&taxid).copied().filter(|_| kept),
                    },
                )
            })
            .collect()
    }
}

/// Where the reads of a failed candidate above species went, in words:
/// `verb` the row that took them, or that they stay unexplained.
fn destination_in_words(verb: &str, end: Option<&HigherRow>, names: &Names) -> String {
    end.map_or_else(
        || "stay unexplained".to_string(),
        |end| format!("{verb} {} {}", end.rank.as_str(), names.taxon(end.taxid)),
    )
}

/// The failed tests of a candidate above species, in words.
fn failures_in_words(failures: &[higher::Failure]) -> Vec<String> {
    failures.iter().map(ToString::to_string).collect()
}

/// Where the reads that reach a taxon above genus come from, in words.
fn taxon_parts_in_words(candidate: &TaxonCandidate, names: &Names) -> String {
    let mut parts = Vec::new();
    if candidate.tied > 0 {
        let genera = candidate
            .tied_genera
            .iter()
            .map(|&genus| names.taxon(genus))
            .collect::<Vec<_>>();
        parts.push(format!(
            "{} tied across genera (most often {})",
            candidate.tied,
            genera.join(", ")
        ));
    }
    if candidate.from_genera > 0 {
        parts.push(format!(
            "{} unique to genera without a genus row",
            candidate.from_genera
        ));
    }
    if candidate.from_below > 0 {
        parts.push(format!(
            "{} moved up from lower ranks",
            candidate.from_below
        ));
    }
    format!("{} reads: {}", candidate.fragments(), parts.join(", "))
}

#[cfg(test)]
mod tests {
    use super::super::genomes::{StrainGenome, StrainShares};
    use super::super::higher::{Failure, Kind, Measured};
    use super::super::profile_row;
    use super::*;
    use crate::profiler::evidence::PresenceAssembly;
    use crate::profiler::output::ProfileRow;
    use crate::profiler::unclassified;

    fn header() -> TableHeader {
        TableHeader {
            fingerprint: None,
            sequencer: Sequencer::Illumina,
            paired: true,
            host: None,
        }
    }

    /// Family 1000 holds genus 100 (species 1 to 5, species 4 with strain 400) and family 2000 nothing named.
    fn taxonomy() -> Taxonomy {
        let mut nodes = HashMap::from([
            (1, (1, "no rank".to_string())),
            (1000, (1, "family".to_string())),
            (2000, (1, "family".to_string())),
            (100, (1000, "genus".to_string())),
            (400, (4, "strain".to_string())),
        ]);
        let mut names = HashMap::from([
            (1000, "Family A".to_string()),
            (2000, "Family B".to_string()),
            (100, "Genus".to_string()),
            (400, "Strain 400".to_string()),
        ]);
        for taxid in 1..=5 {
            nodes.insert(taxid, (100, "species".to_string()));
            names.insert(taxid, format!("Species {taxid}"));
        }
        Taxonomy::from_parts(nodes, names, HashMap::new(), HashSet::new())
    }

    /// Accessions `GCF_<id>` of assembly ids 0-99.
    fn accessions() -> Vec<Arc<str>> {
        (0..100).map(|id| Arc::from(format!("GCF_{id}"))).collect()
    }

    /// A profile row of `rank` for `taxid` on genome `GCF_<genome>`.
    fn final_row(taxid: u32, rank: &str, genome: u32, reads: f64, abundance: f64) -> ProfileRow {
        let mut row = profile_row(
            format!("Taxon {taxid}"),
            taxid,
            rank,
            (reads, 1.0),
            &EvidenceCoverage::default(),
            Arc::from(format!("GCF_{genome}")),
            Some((f64::NAN, f64::NAN)),
        );
        row.abundance = abundance;
        row
    }

    fn coverage(fixed_chunk: f64, coebr: f64) -> EvidenceCoverage {
        EvidenceCoverage {
            fixed_chunk,
            oebr: 1.0,
            coebr,
            window_spread: Some(0.8),
            ..EvidenceCoverage::default()
        }
    }

    fn evidence(taxid: u32, count: u64, coverage: EvidenceCoverage) -> GenomeEvidence {
        GenomeEvidence {
            taxid,
            assembly: PresenceAssembly {
                count,
                ..PresenceAssembly::default()
            },
            coverage: Some(coverage),
            kept: false,
        }
    }

    /// A small run: species 1 kept (genomes 10 and 11), with a species row and a genome row; species 2 reported at
    /// its genus; species 3 removed; species 4 kept on the second look, in a genome row under strain 400 only;
    /// species 5 removed; genus 100 with a genus row from reads of species not kept, which its reads without a
    /// species-level match join although they fail, and 5 reads tied with other genera; family 1000 kept, family
    /// 2000 not.
    struct Fixture {
        genomes: HashMap<u32, GenomeEvidence>,
        counts: FragmentCounts,
        tests: SpeciesTests,
        outcomes: HashMap<u32, SpeciesDecision>,
        reported: HashMap<u32, u32>,
        diagnostics: Vec<StrainDiagnostic>,
        species: HashMap<u32, SpeciesGenomes>,
        strain_taxids: HashMap<Arc<str>, u32>,
        higher: HigherLevels,
        genus_members: Vec<(u32, u32, f64)>,
        blocks: Vec<TaxonBlock>,
        dropped: Dropped,
        labels: HashMap<u32, ProfileRank>,
        named_first: HashMap<u32, Vec<u32>>,
    }

    fn fixture() -> Fixture {
        let genomes = HashMap::from([
            (10, evidence(1, 50, coverage(0.9, 1.0))),
            (11, evidence(1, 80, coverage(0.02, 0.61))),
            (20, evidence(2, 40, coverage(0.9, 1.0))),
            (30, evidence(3, 60, coverage(0.9, 1.0))),
            (40, evidence(4, 70, coverage(0.9, 1.0))),
            (50, evidence(5, 30, coverage(0.1, 0.9))),
        ]);
        let counts = FragmentCounts {
            candidates: HashMap::from([(1, 100), (2, 50), (3, 100), (4, 40)]),
            won: HashMap::from([(1, 90), (2, 10), (3, 5), (4, 30)]),
            second_look: HashMap::from([(4, 7)]),
            ..FragmentCounts::default()
        };
        let tests = SpeciesTests {
            multiple_tests: 1,
            chunk_cutoff: 0.5,
            ratio_bounds: (0.6, 0.75),
            markers: None,
            min_marker_genes: 0,
        };
        let reported = HashMap::from([(1, 10), (2, 20), (4, 40)]);
        let mut outcomes = judged_species(
            &genomes,
            &reported,
            &Looks {
                second_look: &HashMap::from([(4, vec![0])]),
                own_evidence: &HashMap::from([(2, true)]),
                prophage: &[],
                named_first: &HashMap::new(),
                even_coverage: &[],
            },
        );
        outcomes.insert(2, SpeciesDecision::AtGenusTie);
        let diagnostic =
            |taxid, assembly_id, best_fit, em_share, coverage, kept, gated| StrainDiagnostic {
                taxid,
                assembly_id,
                best_fit,
                em_share,
                resolved: 5,
                gated,
                coverage,
                kept,
            };
        let measured = Some(coverage(0.9, 1.0));
        let diagnostics = vec![
            diagnostic(4, 41, false, 0.1, None, false, 0),
            diagnostic(4, 40, true, 0.9, measured, true, 5),
            diagnostic(1, 10, true, 1.0, measured, true, 5),
        ];
        let strain = |assembly_id| StrainGenome {
            assembly_id,
            reads: 1.0,
            bases: 1.0,
            mean_identity: 1.0,
            coverage: EvidenceCoverage::default(),
            asm_len: 1.0,
        };
        let genomes_of = |best_fit| SpeciesGenomes {
            best_fit,
            mean_identity: 1.0,
            strains: vec![strain(best_fit)],
            shares: StrainShares::default(),
        };
        let species = HashMap::from([(1, genomes_of(10)), (4, genomes_of(40))]);
        let measured = |genome, fragments, coverage, spread| Measured {
            genome,
            fragments,
            coverage,
            spread,
            length: 1_000_000.0,
        };
        let genus = measured(50, 240, coverage(0.1, 0.9), Some(0.8));
        let family = measured(80, 300, coverage(0.2, 0.9), Some(0.7));
        let higher = HigherLevels {
            rows: vec![
                HigherRow {
                    taxid: 100,
                    rank: Rank::Genus,
                    reads: 475.0,
                    unique: [220, 250],
                    depth: 0.1,
                    measured: genus,
                },
                HigherRow {
                    taxid: 1000,
                    rank: Rank::Family,
                    reads: 320.0,
                    unique: [0, 0],
                    depth: 0.1,
                    measured: family,
                },
            ],
            genera: vec![
                GenusCandidate {
                    genus: 100,
                    kind: Kind::BelowSpecies,
                    unique: 220,
                    tied: 30,
                    identity: 0.9,
                    measured: measured(51, 200, coverage(0.05, 0.2), Some(0.6)),
                    failures: vec![Failure::ChunkRatio(0.2)],
                    shown: false,
                },
                GenusCandidate {
                    genus: 100,
                    kind: Kind::NotKept,
                    unique: 250,
                    tied: 20,
                    identity: 0.99,
                    measured: genus,
                    failures: Vec::new(),
                    shown: true,
                },
            ],
            taxa: vec![
                TaxonCandidate {
                    taxid: 1000,
                    rank: Rank::Family,
                    tied: 100,
                    from_genera: 220,
                    from_below: 0,
                    tied_genera: vec![100],
                    measured: family,
                    failures: Vec::new(),
                },
                TaxonCandidate {
                    taxid: 2000,
                    rank: Rank::Family,
                    tied: 250,
                    from_genera: 0,
                    from_below: 0,
                    tied_genera: vec![100],
                    measured: measured(90, 20, coverage(0.01, 0.5), None),
                    failures: vec![Failure::Spread(None)],
                },
            ],
            identity_floor: 0.886,
            fragments: [620, 270],
            unexplained: 250,
            shares: Vec::new(),
            moved: HashMap::from([(2000, None)]),
            genus_start: HashMap::from([(100, Some(1000))]),
        };
        let block = |taxid, taxon: Option<ProfileRow>, genomes| TaxonBlock {
            taxid,
            unit: taxon.as_ref().is_none_or(|row| row.rank == "species"),
            taxon,
            genomes,
            breadth: 0.5,
        };
        let blocks = vec![
            block(
                1,
                Some(final_row(1, "species", 10, 120.0, 60.0)),
                vec![final_row(1, "strain", 10, 30.0, 15.0)],
            ),
            block(
                100,
                Some(final_row(100, "genus", 50, 262.0, 20.0)),
                Vec::new(),
            ),
            block(4, None, vec![final_row(400, "strain", 40, 8.0, 4.0)]),
            block(
                1000,
                Some(final_row(1000, "family", 80, 10.0, 1.0)),
                Vec::new(),
            ),
        ];
        Fixture {
            genomes,
            counts,
            tests,
            outcomes,
            reported,
            diagnostics,
            species,
            strain_taxids: HashMap::from([(Arc::from("GCF_40"), 400)]),
            higher,
            genus_members: vec![(100, 2, 12.0)],
            blocks,
            dropped: Dropped::default(),
            labels: HashMap::new(),
            named_first: HashMap::new(),
        }
    }

    /// The fixture's decisions table: its two header lines and its rows by column name.
    fn decisions(fixture: &Fixture) -> (Vec<String>, Vec<HashMap<&'static str, String>>) {
        let table = Decisions {
            genomes: &fixture.genomes,
            counts: &fixture.counts,
            tests: &fixture.tests,
            outcomes: &fixture.outcomes,
            reported: &fixture.reported,
            strains: Some(StrainLevel {
                rows: &fixture.diagnostics,
                species: &fixture.species,
                tests: SpeciesTests {
                    multiple_tests: 1,
                    chunk_cutoff: 0.5,
                    ratio_bounds: (0.75, 0.75),
                    markers: None,
                    min_marker_genes: 0,
                },
                strain_taxids: &fixture.strain_taxids,
            }),
            genus_members: &fixture.genus_members,
            higher: Some(&fixture.higher),
            dropped: &fixture.dropped,
            labels: &fixture.labels,
            named_first: &fixture.named_first,
        };
        let mut out = Vec::new();
        table
            .write(
                &mut out,
                &header(),
                &Names {
                    taxonomy: &taxonomy(),
                    genomes: &accessions(),
                },
                &FinalRows::of(&fixture.blocks),
            )
            .unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines = text.lines().collect::<Vec<_>>();
        let rows = lines[2..]
            .iter()
            .map(|line| {
                let fields = line.split('\t').collect::<Vec<_>>();
                assert_eq!(fields.len(), DECISION_COLUMNS.len(), "{line}");
                DECISION_COLUMNS
                    .iter()
                    .zip(fields)
                    .map(|(&column, field)| (column, field.to_string()))
                    .collect()
            })
            .collect();
        (
            lines[..2].iter().map(|line| line.to_string()).collect(),
            rows,
        )
    }

    /// The values of `columns` in each row of `rank`.
    fn columns(
        rows: &[HashMap<&'static str, String>],
        rank: &str,
        columns: &[&str],
    ) -> Vec<Vec<String>> {
        rows.iter()
            .filter(|row| row["rank"] == rank)
            .map(|row| columns.iter().map(|&c| row[c].clone()).collect())
            .collect()
    }

    fn strings<const N: usize>(rows: &[[&str; N]]) -> Vec<Vec<String>> {
        rows.iter()
            .map(|row| row.iter().map(|v| v.to_string()).collect())
            .collect()
    }

    #[test]
    fn the_profile_opens_with_its_run_and_the_name_of_every_column() {
        let text = header().lines("profile", &PROFILE_COLUMNS);
        let lines = text.lines().collect::<Vec<_>>();
        assert_eq!(lines.len(), 2);
        assert_eq!(
            lines[0],
            format!(
                "#metax_profile\tversion={}\tengine=level\tdatabase_fingerprint=NA\tsequencer=Illumina\tpaired=true",
                env!("CARGO_PKG_VERSION")
            )
        );
        assert!(lines[1].starts_with("#name\ttaxid\trank\treads\tdepth\tabundance\t"));
        assert!(
            lines[1].ends_with("\tevidence_depth\tbest_fit_genome\tmean_identity\testimated_ani")
        );
        let with_database = TableHeader {
            fingerprint: Some(0xabc),
            ..header()
        };
        assert!(with_database
            .lines("profile", &PROFILE_COLUMNS)
            .contains("\tdatabase_fingerprint=00000000000000000000000000000abc\t"));
        // With --host: the host on the first line, its three columns last.
        let with_host = TableHeader {
            host: Some(9606),
            ..header()
        };
        let columns = [&PROFILE_COLUMNS[..], &HOST_COLUMNS[..]].concat();
        let text = with_host.lines("profile", &columns);
        let lines = text.lines().collect::<Vec<_>>();
        assert!(lines[0].ends_with("\tpaired=true\thost=9606"));
        assert!(lines[1].ends_with("\testimated_ani\thost_names\thost_taxids\tdiseases"));
        // Every row has a column per name: a species row and the unclassified row, also when no row precedes it.
        let mut writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(Vec::new());
        final_row(1, "species", 10, 5.0, 100.0)
            .write(&mut writer)
            .unwrap();
        unclassified::write_row(&mut writer, None).unwrap();
        let text = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        let rows = text.lines().collect::<Vec<_>>();
        assert_eq!(rows.len(), 2);
        assert!(rows
            .iter()
            .all(|row| row.split('\t').count() == PROFILE_COLUMNS.len()));
    }

    #[test]
    fn the_decisions_table_opens_with_its_run_and_its_columns() {
        let (lines, rows) = decisions(&fixture());
        assert_eq!(
            lines[0],
            format!(
                "#metax_decisions\tversion={}\tengine=level\tdatabase_fingerprint=NA\tsequencer=Illumina\tpaired=true",
                env!("CARGO_PKG_VERSION")
            )
        );
        assert_eq!(lines[1], format!("#{}", DECISION_COLUMNS.join("\t")));
        // Species rows, strain rows, then genus rows and the rows above genus.
        let ranks = rows
            .iter()
            .map(|row| row["rank"].as_str())
            .collect::<Vec<_>>();
        assert_eq!(
            ranks,
            [
                "species", "species", "species", "species", "species", "species", "strain",
                "strain", "strain", "genus", "genus", "family", "family"
            ]
        );
    }

    #[test]
    fn species_rows_tell_the_decided_genome_its_tests_and_the_reason() {
        let (_, rows) = decisions(&fixture());
        let found = columns(
            &rows,
            "species",
            &[
                "taxid",
                "genome",
                "parent_taxid",
                "decided",
                "genome_tests",
                "decision",
                "reported_as",
                "final_reads",
                "final_abundance",
            ],
        );
        // By abundance in the profile (species 1, then 4), then species without rows of their own by taxid; the
        // decided genome first.
        let expected = strings(&[
            ["1", "GCF_10", "100", "true", "passed", "kept", "species", "150.00000", "75.00000"],
            [
                "1",
                "GCF_11",
                "100",
                "false",
                "failed amount: fixed chunk 0.0200 <= cutoff 0.5000; failed evenness: COEBR 0.6100 < 0.75",
                "kept",
                "species",
                "NA",
                "NA",
            ],
            ["4", "GCF_40", "100", "true", "passed", "kept", "species", "8.00000", "4.00000"],
            [
                "2",
                "GCF_20",
                "100",
                "true",
                "failed won share: 0.2000 < 0.3",
                "removed",
                "genus",
                "NA",
                "NA",
            ],
            [
                "3",
                "GCF_30",
                "100",
                "true",
                "failed won share: 0.0500 < 0.3",
                "removed",
                "none",
                "NA",
                "NA",
            ],
            [
                "5",
                "GCF_50",
                "100",
                "true",
                "failed amount: fixed chunk 0.1000 <= cutoff 0.5000; failed won share: 0.0000 < 0.3",
                "removed",
                "none",
                "NA",
                "NA",
            ],
        ]);
        assert_eq!(found, expected);
        let reasons = columns(&rows, "species", &["taxid", "reason"]);
        let reason = |taxid: &str| {
            reasons
                .iter()
                .find(|row| row[0] == taxid)
                .map(|row| row[1].clone())
                .unwrap()
        };
        assert_eq!(reason("1"), "passed amount, evenness and won share");
        assert_eq!(
            reason("4"),
            "kept on the second look: a removed sister no longer splits its reads"
        );
        assert_eq!(
            reason("2"),
            "reported at genus: its reads tie with relatives along its genome"
        );
        assert_eq!(reason("3"), "removed: failed won share: 0.0500 < 0.3");
        // Counts of the species on each of its rows.
        let counts = columns(
            &rows,
            "species",
            &[
                "candidate_reads",
                "won_reads",
                "second_look_reads",
                "own_reads",
            ],
        );
        assert_eq!(counts[0], ["100", "90", "0", "0"]);
        assert_eq!(counts[2], ["40", "30", "7", "0"]);
    }

    #[test]
    fn strain_rows_come_by_species_and_name_their_profile_row() {
        let (_, rows) = decisions(&fixture());
        let found = columns(
            &rows,
            "strain",
            &[
                "taxid",
                "name",
                "parent_taxid",
                "genome",
                "decided",
                "decision",
                "reason",
                "final_reads",
            ],
        );
        let expected = strings(&[
            [
                "1",
                "Species 1 GCF_10",
                "1",
                "GCF_10",
                "true",
                "kept",
                "same strain: its reads match it within sequencing error and cover it evenly",
                "30.00000",
            ],
            [
                "400",
                "Strain 400",
                "4",
                "GCF_40",
                "true",
                "kept",
                "same strain: its reads match it within sequencing error and cover it evenly",
                "8.00000",
            ],
            [
                "4",
                "Species 4 GCF_41",
                "4",
                "GCF_41",
                "false",
                "removed",
                "low identity: its reads match it below strain-level identity",
                "NA",
            ],
        ]);
        assert_eq!(found, expected);
    }

    #[test]
    fn genus_rows_list_each_kind_of_unique_reads_and_the_species_reported_there() {
        let (_, rows) = decisions(&fixture());
        let found = columns(
            &rows,
            "genus",
            &[
                "taxid",
                "parent_taxid",
                "genome",
                "decided",
                "reads",
                "candidate_reads",
                "won_reads",
                "window_spread",
                "genome_tests",
                "decision",
                "reason",
                "reported_as",
                "final_reads",
            ],
        );
        let reason = "kept: 220 reads without a species-level match unique to it, identity 0.9000 >= 0.8860, chunk \
                      ratio 0.2000 < 0.4; they joined the genus row from reads of species not kept + 250 reads of \
                      species not kept unique to it, identity 0.9900 >= 0.8860, chunk ratio 0.9000 >= 0.4 + 5.0 \
                      reads tied with other genera + species reported at genus: 2 (12 reads; its reads tie with \
                      relatives along its genome)";
        let expected = strings(&[
            [
                "100",
                "1000",
                "GCF_50",
                "true",
                "240",
                "270",
                "250",
                "0.8000",
                "reads of species not kept: passed",
                "kept",
                reason,
                "genus",
                "262.00000",
            ],
            [
                "100",
                "1000",
                "GCF_51",
                "false",
                "200",
                "250",
                "220",
                "0.6000",
                "reads without a species-level match: failed chunk ratio: 0.2000 < 0.4",
                "kept",
                reason,
                "genus",
                "NA",
            ],
        ]);
        assert_eq!(found, expected);
        // Without a candidate whose evidence the genus row shows, the species reported at the genus are its
        // decided row.
        let mut fixture = fixture();
        fixture
            .higher
            .genera
            .retain(|c| c.kind == Kind::BelowSpecies);
        fixture.higher.rows.retain(|row| row.rank != Rank::Genus);
        let (_, rows) = decisions(&fixture);
        let found = columns(
            &rows,
            "genus",
            &["genome", "decided", "reads", "genome_tests", "final_reads"],
        );
        assert_eq!(
            found,
            strings(&[
                [
                    "GCF_20",
                    "true",
                    "12.00000",
                    "no genus tests: made of species reported at their genus",
                    "262.00000"
                ],
                [
                    "GCF_51",
                    "false",
                    "200",
                    "reads without a species-level match: failed chunk ratio: 0.2000 < 0.4",
                    "NA"
                ],
            ])
        );
    }

    #[test]
    fn genus_candidates_without_a_row_say_where_their_reads_went() {
        let mut fixture = fixture();
        fixture.genus_members.clear();
        fixture
            .higher
            .genera
            .retain(|c| c.kind == Kind::BelowSpecies);
        fixture.higher.rows.retain(|row| row.rank != Rank::Genus);
        fixture.blocks.retain(|block| block.taxid != 100);
        fixture.outcomes.insert(2, SpeciesDecision::KeptTests);
        let (_, rows) = decisions(&fixture);
        let found = columns(
            &rows,
            "genus",
            &["decided", "decision", "reason", "reported_as"],
        );
        assert_eq!(
            found,
            strings(&[[
                "true",
                "removed",
                "removed: 220 reads without a species-level match unique to it, identity 0.9000 >= 0.8860, chunk \
                 ratio 0.2000 < 0.4; they go to family Family A",
                "family"
            ]])
        );
        // Without a row above them, the reads stay unexplained.
        fixture.higher.genus_start.insert(100, None);
        let (_, rows) = decisions(&fixture);
        let found = columns(&rows, "genus", &["reason", "reported_as"]);
        assert_eq!(
            found,
            strings(&[[
                "removed: 220 reads without a species-level match unique to it, identity 0.9000 >= 0.8860, chunk \
                 ratio 0.2000 < 0.4; they stay unexplained",
                "none"
            ]])
        );
    }

    #[test]
    fn taxa_above_genus_are_listed_kept_or_not_with_their_tests() {
        let (_, rows) = decisions(&fixture());
        let found = columns(
            &rows,
            "family",
            &[
                "taxid",
                "genome",
                "reads",
                "candidate_reads",
                "window_spread",
                "genome_tests",
                "decision",
                "reason",
                "reported_as",
                "final_reads",
            ],
        );
        let expected = strings(&[
            [
                "1000",
                "GCF_80",
                "300",
                "320",
                "0.7000",
                "passed",
                "kept",
                "kept: 320 reads: 100 tied across genera (most often Genus), 220 unique to genera without a genus \
                 row; window spread 0.7000 >= 0.5",
                "family",
                "10.00000",
            ],
            [
                "2000",
                "GCF_90",
                "20",
                "250",
                "NA",
                "failed window spread: too few reads to measure",
                "removed",
                "removed: 250 reads: 250 tied across genera (most often Genus); failed window spread: too few reads \
                 to measure; they stay unexplained",
                "none",
                "NA",
            ],
        ]);
        assert_eq!(found, expected);
    }

    #[test]
    fn every_profile_row_matches_one_kept_decided_row() {
        let fixture = fixture();
        let (_, rows) = decisions(&fixture);
        let kept_decided = |rank: &str, taxid: u32| {
            rows.iter()
                .filter(|row| {
                    row["rank"] == rank
                        && row["taxid"] == taxid.to_string()
                        && row["decision"] == "kept"
                        && row["decided"] == "true"
                })
                .collect::<Vec<_>>()
        };
        let five = |value: f64| format_five(value);
        let mut matched = 0;
        for block in &fixture.blocks {
            let species_row = block.taxon.as_ref().filter(|row| row.rank == "species");
            if species_row.is_none() && block.genomes.is_empty() {
                let row = block.taxon.as_ref().unwrap();
                let found = kept_decided(&row.rank, block.taxid);
                assert_eq!(found.len(), 1, "{}", block.taxid);
                assert_eq!(found[0]["final_reads"], five(row.equivalent_mapcount));
                assert_eq!(found[0]["final_abundance"], five(row.abundance));
                matched += 1;
                continue;
            }
            let found = kept_decided("species", block.taxid);
            assert_eq!(found.len(), 1, "{}", block.taxid);
            let total = species_row
                .into_iter()
                .chain(&block.genomes)
                .fold(0.0, |total, row| total + row.equivalent_mapcount);
            assert_eq!(found[0]["final_reads"], five(total));
            matched += 1;
            for row in &block.genomes {
                let genome = &row.best_fit.as_ref().unwrap().genome;
                let found = rows
                    .iter()
                    .filter(|r| {
                        r["rank"] == "strain" && r["genome"] == **genome && r["decision"] == "kept"
                    })
                    .collect::<Vec<_>>();
                assert_eq!(found.len(), 1, "{genome}");
                assert_eq!(found[0]["taxid"], row.taxid.to_string());
                assert_eq!(found[0]["final_reads"], five(row.equivalent_mapcount));
                assert_eq!(found[0]["final_abundance"], five(row.abundance));
                matched += 1;
            }
        }
        // No kept row beyond those the profile holds.
        let kept = rows
            .iter()
            .filter(|row| row["decision"] == "kept" && row["final_reads"] != "NA")
            .count();
        assert_eq!(kept, matched);
    }

    #[test]
    fn taxa_dropped_by_a_minimum_are_removed_with_the_reason() {
        // At least 20 reads: species 4 (8 reads, in its genome row) and family 1000 (10 reads) are dropped; species
        // 1 (120 in its species row, 150 with its genome row) and genus 100 (262) stay.
        let mut fixture = fixture();
        fixture.dropped.taxa = super::super::drop_below_minimums(
            &mut fixture.blocks,
            &super::super::Minimums {
                reads: Some(20),
                breadth: None,
            },
        )
        .into_iter()
        .map(|(taxid, below)| (taxid, super::super::Removal::BelowMinimum(below)))
        .collect();
        let mut dropped = fixture.dropped.taxa.keys().copied().collect::<Vec<_>>();
        dropped.sort_unstable();
        assert_eq!(dropped, [4, 1000]);
        let (_, rows) = decisions(&fixture);
        let row = |rank: &str, genome: &str| {
            rows.iter()
                .find(|row| row["rank"] == rank && row["genome"] == genome)
                .map(|row| {
                    ["decision", "reason", "reported_as", "final_reads"]
                        .map(|column| row[column].clone())
                })
                .unwrap()
        };
        let removed = "removed by --min-reads: 8.00 reads < 20";
        assert_eq!(row("species", "GCF_40"), ["removed", removed, "none", "NA"]);
        assert_eq!(
            row("strain", "GCF_40"),
            [
                "removed",
                &format!("species not reported: {removed}"),
                "none",
                "NA"
            ]
        );
        assert_eq!(
            row("family", "GCF_80"),
            [
                "removed",
                "removed by --min-reads: 10.00 reads < 20",
                "none",
                "NA"
            ]
        );
        assert_eq!(row("species", "GCF_10")[0], "kept");
        assert_eq!(row("genus", "GCF_50")[0], "kept");
    }

    #[test]
    fn taxa_that_are_no_pathogens_of_the_host_are_removed_with_the_reason() {
        // Species 4 (in its genome row under strain 400) and family 1000 are listed for another host only; species 1,
        // genus 100 are not listed and stay.
        let table = HashMap::from([
            (
                4,
                crate::profiler::pathogen::PathogenEntry {
                    host_taxids: "10090".to_string(),
                    host_names: "Mus musculus".to_string(),
                    diseases: "unknown".to_string(),
                },
            ),
            (
                1000,
                crate::profiler::pathogen::PathogenEntry {
                    host_taxids: "10090".to_string(),
                    host_names: "Mus musculus".to_string(),
                    diseases: "unknown".to_string(),
                },
            ),
        ]);
        let mut fixture = fixture();
        let taxonomy = taxonomy();
        super::super::host::HostFilter {
            host: 9606,
            table: Some(&table),
            taxonomy: &taxonomy,
        }
        .apply(&mut fixture.blocks, &mut fixture.dropped);
        let (_, rows) = decisions(&fixture);
        let row = |rank: &str, genome: &str| {
            rows.iter()
                .find(|row| row["rank"] == rank && row["genome"] == genome)
                .map(|row| {
                    ["decision", "reason", "reported_as", "final_reads"]
                        .map(|column| row[column].clone())
                })
                .unwrap()
        };
        let removed = "removed by --host: not a pathogen of host 9606";
        assert_eq!(row("species", "GCF_40"), ["removed", removed, "none", "NA"]);
        assert_eq!(row("strain", "GCF_40"), ["removed", removed, "none", "NA"]);
        assert_eq!(row("family", "GCF_80"), ["removed", removed, "none", "NA"]);
        assert_eq!(row("species", "GCF_10")[0], "kept");
        assert_eq!(row("genus", "GCF_50")[0], "kept");
    }

    #[test]
    fn units_labelled_above_species_are_reported_at_their_rank() {
        // Unit 4 (in its genome row under strain 400) is labelled at genus; unit 3, labelled at genus too, was
        // removed as named species 1 explains its reads.
        let mut fixture = fixture();
        fixture.labels = HashMap::from([(4, ProfileRank::Genus), (3, ProfileRank::Genus)]);
        fixture.named_first = HashMap::from([(3, vec![1])]);
        fixture
            .outcomes
            .insert(3, SpeciesDecision::RemovedNamedFirst);
        let (_, rows) = decisions(&fixture);
        let row = |genome: &str| {
            rows.iter()
                .find(|row| row["genome"] == genome && row["rank"] != "strain")
                .map(|row| {
                    [
                        "rank",
                        "parent_taxid",
                        "decision",
                        "reason",
                        "reported_as",
                        "final_reads",
                    ]
                    .map(|column| row[column].clone())
                })
                .unwrap()
        };
        assert_eq!(
            row("GCF_40"),
            [
                "genus",
                "1000",
                "kept",
                "kept on the second look: a removed sister no longer splits its reads; its references are labelled \
                 only at genus",
                "genus",
                "8.00000",
            ]
        );
        assert_eq!(
            row("GCF_30"),
            [
                "genus",
                "1000",
                "removed",
                "removed: its reads are explained by named units of its taxon: Species 1 (1); its references are \
                 labelled only at genus",
                "none",
                "NA",
            ]
        );
        // A species keeps its rank and parent.
        assert_eq!(row("GCF_10")[..2], ["species", "100"]);
    }
}

//! Rows above species: the fragments that no kept species takes, placed at the
//! lowest rank that their matches and their evidence support.
//!
//! Two kinds of fragments come here:
//!
//! - fragments without a species-level candidate, placed at the floor gates
//!   (`--genus-identity`, a fixed mapped length on short reads or the species
//!   length cutoff on long reads, and a fixed aligned fraction); on noisy long
//!   reads not those whose best match, or one tied with it, is a kept
//!   species, which count for that species;
//! - fragments with species-level candidates of which none is kept, placed at
//!   the species gates.
//!
//! A fragment's tie set holds the taxa (species, or units labelled above
//! species) whose placements lie within the tie margin of its best one, on
//! genome references and marker genomes alike.
//!
//! A candidate row's fragments are measured on the genome they fit best as a
//! set, its evidence genome: the genome holding most of their best
//! placements. Its tests come from that genome's evidence model, as the
//! species tests do ([`RowRules`]): [`GENOME_ROWS`] for a genome reference,
//! [`MARKER_ROWS`] for a genome represented by marker genes, whose genes stand
//! for its chunks and on which a fragment counts on every gene it covers.
//!
//! A kept unit represented by marker genes leaves part of its own fragments
//! to its relatives: its marker set misses genes of its genome, and its genes
//! clip reads at their ends; those fragments land on relatives' genes or
//! genomes. A taxon above such kept units therefore gets a row only when its
//! fragments reach the model's spill share of theirs ([`kept_beneath`]),
//! counted in whole-genome equivalents when the row's evidence is a genome
//! reference. Kept units on genome references leave no such share.
//!
//! 1. Genus rows. A fragment whose tie set lies in one genus is unique to it.
//!    A genus gets a row when, for one kind of fragment or both, the fragments
//!    of that kind unique to it meet the amount of their evidence model, their
//!    mean identity is at least the run's species identity cutoff minus the
//!    model's margin, and their chunk ratio (COEBR) on the evidence genome
//!    reaches the model's bound; its unique fragments of both kinds must also
//!    reach the spill bound of the kept units beneath the genus. Its row then
//!    holds the fragments of both kinds unique to it; the unique fragments of
//!    a genus without a row go to the first taxon above the genus at family,
//!    order, class, phylum or domain rank (step 3).
//! 2. Ties between genera. A fragment whose tie set spans several genera, one
//!    or more of them with a genus row, goes to the genera with a genus row,
//!    shared by the likelihood EM over those genera, as the species level's EM
//!    shares fragments among kept species; to step 3 when none of its tie
//!    set's genera has one.
//! 3. Rows above genus. The other fragments go to the lowest common ancestor of
//!    their tie set's taxa, or the first taxon above it at family, order,
//!    class, phylum or domain rank; a fragment with none below the root stays
//!    unexplained. A taxon gets a row when the fragments that reach it meet
//!    the amount and the evenness of their evidence model ([`taxon_failures`])
//!    and reach the spill bound of the kept units beneath it; otherwise they
//!    move to the next taxon above it at one of those ranks, and stay
//!    unexplained above domain.
//!
//! A row's depth is that of the genome its fragments fit best as a set: the
//! genome holding most of their best placements, on which the row's tests
//! measure them (for a genus row, the unique fragments of the kind its row
//! shows). On a genome reference it is the first-read bases of all its
//! fragments over the genome's indexed length; on a marker genome, as for a
//! species reported on one, each fragment adds its read length on every gene
//! it covers over the placements that the genome's genes offer a read of that
//! length. Genomes tied exactly with it on that number share the depth: it is
//! the mean of the depths over each of them. A row on a marker genome explains
//! its fragments times that genome's whole-genome scale in the unclassified
//! estimate, as a species row does.

use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::ops::Range;
use std::sync::Arc;

use super::{
    best_tied, fold_candidates, placement_identity, ranked, share_reads, CoverageTests, EmSettings,
    Fragment,
};
use crate::profiler::community::MappingData;
use crate::profiler::evidence::{
    EvidenceCoverage, EvidencePolicy, PresenceAssembly, PresenceInterval,
};
use crate::profiler::likelihood_em;
use crate::profiler::markers::MarkerModel;
use crate::profiler::rollup::GenomeLayout;
use crate::profiler::unclassified::GenomeReads;
use crate::taxonomy::Taxonomy;

/// What a row above species must show, by the evidence model of its
/// evidence genome. A candidate's amount comes first: fragments short of it
/// are not judged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct RowRules {
    /// Fewest fragments: a genus row's unique fragments of one kind, or the
    /// fragments reaching a taxon above genus.
    pub(super) min_fragments: usize,
    /// Fewest genes hit on the evidence genome.
    pub(super) min_genes: usize,
    /// How far the mean identity of a genus row's unique fragments may lie
    /// below the run's species identity cutoff.
    pub(super) genus_identity_margin: f64,
    /// Lowest chunk ratio (COEBR) of a genus row's unique fragments.
    pub(super) genus_min_chunk_ratio: f64,
    /// How the fragments of a row above genus spread over the evidence genome.
    pub(super) above_genus: Spread,
    /// Share of a kept unit's fragments that may land on its relatives: a row
    /// of a taxon above kept units of this model needs at least this share of
    /// their fragments.
    pub(super) spill_share: f64,
}

/// The evenness test of a row above genus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Spread {
    /// Lowest share of the evidence genome's windows that the fragments reach.
    Windows(f64),
    /// Lowest chunk ratio (COEBR) of the fragments.
    ChunkRatio(f64),
}

/// Rows on a genome reference: fragments spread over its sequence.
pub(super) const GENOME_ROWS: RowRules = RowRules {
    min_fragments: 200,
    min_genes: 0,
    genus_identity_margin: 0.09,
    genus_min_chunk_ratio: 0.40,
    above_genus: Spread::Windows(0.50),
    spill_share: 0.0,
};

/// Rows on a genome represented by marker genes: fragments spread over its
/// genes, whose number stands for the amount since the genes hold a small
/// share of a genome's reads.
pub(super) const MARKER_ROWS: RowRules = RowRules {
    min_fragments: 0,
    min_genes: 10,
    genus_identity_margin: 0.09,
    genus_min_chunk_ratio: 0.40,
    above_genus: Spread::ChunkRatio(0.40),
    spill_share: 0.15,
};

/// The evidence model of a genome.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Model {
    Genome,
    Markers,
}

impl Model {
    /// The model of genome `assembly_id`.
    pub(super) fn of(data: &MappingData, assembly_id: u32) -> Self {
        if data
            .marker_model
            .as_deref()
            .is_some_and(|model| model.is_marker_genome(assembly_id))
        {
            Self::Markers
        } else {
            Self::Genome
        }
    }

    /// Its row rules.
    pub(super) fn rules(self) -> &'static RowRules {
        match self {
            Self::Genome => &GENOME_ROWS,
            Self::Markers => &MARKER_ROWS,
        }
    }
}

/// The spill bound of a taxon from the kept units beneath it: the fragments
/// they may leave to relatives, counted on marker genes and in whole-genome
/// equivalents.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct Beneath {
    pub(super) marker_fragments: f64,
    pub(super) whole_genome: f64,
}

impl Beneath {
    /// The bound of a row whose evidence has `model`.
    fn bound(&self, model: Model) -> f64 {
        match model {
            Model::Markers => self.marker_fragments,
            Model::Genome => self.whole_genome,
        }
    }
}

/// The spill bound of every taxon at genus to domain rank above the kept
/// `units`, each given as (taxid, model of its reported genome, fragments,
/// whole-genome fragments per fragment): the spill share of their model of
/// their fragments, summed; a unit without a whole-genome scale counts its
/// fragments as they are. Units are summed in taxid order.
pub(super) fn kept_beneath(
    units: impl IntoIterator<Item = (u32, Model, f64, Option<f64>)>,
    lineages: &Lineages,
) -> HashMap<u32, Beneath> {
    let mut units = units.into_iter().collect::<Vec<_>>();
    units.sort_unstable_by_key(|unit| unit.0);
    let mut beneath: HashMap<u32, Beneath> = HashMap::new();
    for (taxid, model, fragments, scale) in units {
        let spill = model.rules().spill_share * fragments;
        if spill <= 0.0 {
            continue;
        }
        for &(ancestor, rank) in lineages.of(taxid).iter() {
            if rank.is_some() {
                let bound = beneath.entry(ancestor).or_default();
                bound.marker_fragments += spill;
                bound.whole_genome += spill * scale.unwrap_or(1.0);
            }
        }
    }
    beneath
}

/// The two kinds of fragments placed above species.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum Kind {
    /// No species-level candidate; placed at the floor gates.
    BelowSpecies,
    /// Species-level candidates, none of them kept.
    NotKept,
}

impl Kind {
    pub(super) fn words(self) -> &'static str {
        match self {
            Self::BelowSpecies => "without a species-level match",
            Self::NotKept => "of species not kept",
        }
    }
}

/// The ranks of the rows above species, lowest first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum Rank {
    Genus,
    Family,
    Order,
    Class,
    Phylum,
    Domain,
}

impl Rank {
    /// The rank of a taxonomy rank name; `None` for ranks without rows.
    fn of(rank: &str) -> Option<Self> {
        match rank {
            "genus" => Some(Self::Genus),
            "family" => Some(Self::Family),
            "order" => Some(Self::Order),
            "class" => Some(Self::Class),
            "phylum" => Some(Self::Phylum),
            "domain" | "superkingdom" => Some(Self::Domain),
            _ => None,
        }
    }

    pub(super) fn as_str(self) -> &'static str {
        match self {
            Self::Genus => "genus",
            Self::Family => "family",
            Self::Order => "order",
            Self::Class => "class",
            Self::Phylum => "phylum",
            Self::Domain => "domain",
        }
    }
}

/// A taxon's lineage, the taxon first, each with its rank among [`Rank`].
pub(super) type Lineage = Arc<[(u32, Option<Rank>)]>;

/// Lineages read from the taxonomy, once per taxon.
pub(super) struct Lineages<'t> {
    taxonomy: &'t Taxonomy,
    known: RefCell<HashMap<u32, Lineage>>,
}

impl<'t> Lineages<'t> {
    pub(super) fn new(taxonomy: &'t Taxonomy) -> Self {
        Self {
            taxonomy,
            known: RefCell::new(HashMap::new()),
        }
    }

    pub(super) fn of(&self, taxid: u32) -> Lineage {
        let taxonomy = self.taxonomy;
        Arc::clone(self.known.borrow_mut().entry(taxid).or_insert_with(|| {
            taxonomy
                .get_parents(taxid)
                .into_iter()
                .map(|(ancestor, _, rank)| (ancestor, Rank::of(&rank)))
                .collect()
        }))
    }

    /// The taxa of the lineage of `taxid`, the taxon first.
    pub(super) fn ancestors(&self, taxid: u32) -> Vec<u32> {
        self.of(taxid).iter().map(|node| node.0).collect()
    }

    /// The taxon of the lineage of `taxid` at `rank`: the taxon itself, or its
    /// ancestor there.
    pub(super) fn at_rank(&self, taxid: u32, rank: Rank) -> Option<u32> {
        self.of(taxid)
            .iter()
            .find(|node| node.1 == Some(rank))
            .map(|node| node.0)
    }

    /// The first taxon of the lineage of `taxid` ranked above `rank`.
    fn above(&self, taxid: u32, rank: Rank) -> Option<(u32, Rank)> {
        first_above(&self.of(taxid), rank)
    }

    /// Where the fragments of a tie set of `taxa` that lies in no one genus
    /// start: their lowest common ancestor when it is ranked above genus,
    /// else the first taxon above it that is; `None` when there is none below
    /// the root.
    fn common_start(&self, taxa: &[u32]) -> Option<(u32, Rank)> {
        let first = self.of(taxa[0]);
        let mut at = 0;
        for &taxid in &taxa[1..] {
            let other = self.of(taxid);
            while at < first.len() && !other.iter().any(|node| node.0 == first[at].0) {
                at += 1;
            }
        }
        first_above(&first[at..], Rank::Genus)
    }
}

/// The first taxon of `lineage` ranked above `rank`.
fn first_above(lineage: &[(u32, Option<Rank>)], rank: Rank) -> Option<(u32, Rank)> {
    lineage
        .iter()
        .find_map(|&(taxid, of)| of.filter(|&of| of > rank).map(|of| (taxid, of)))
}

/// What the level places and with what.
pub(super) struct Level<'a> {
    pub(super) data: &'a MappingData,
    pub(super) tests: &'a CoverageTests,
    pub(super) layout: &'a GenomeLayout,
    /// The reads without a species-level candidate (on noisy long reads, not
    /// those whose best match is a kept species), and the floor gates.
    pub(super) below_species: (&'a [usize], EvidencePolicy),
    /// The reads of fragments whose candidate species are none kept, and the
    /// species gates.
    pub(super) not_kept: (&'a [usize], EvidencePolicy),
    pub(super) em: EmSettings,
    /// Keep each fragment's shares for the per-read classifications.
    pub(super) keep_shares: bool,
    /// The fragments of the kept units of each taxon that may land on
    /// relatives ([`kept_beneath`]).
    pub(super) kept_beneath: &'a HashMap<u32, Beneath>,
}

impl Level<'_> {
    /// The gates of the fragments of `kind`.
    fn policy(&self, kind: Kind) -> EvidencePolicy {
        match kind {
            Kind::BelowSpecies => self.below_species.1,
            Kind::NotKept => self.not_kept.1,
        }
    }
}

/// A fragment of the level: its read, its kind, its tie set (an index into
/// the level's tie sets), its best placement's genome, identity, first read's
/// mapped bases, read length and interval, and its hits on further marker
/// genes of that genome (a range of the level's gene hits).
struct Placed {
    read: usize,
    kind: Kind,
    ties: u32,
    assembly: u32,
    identity: f32,
    bases: f32,
    read_len: u32,
    interval: PresenceInterval,
    gene_hits: Range<u32>,
}

/// The level's fragments, their distinct tie sets (taxa by taxid), the
/// indexed length of each genome holding a best placement, and the fragments'
/// hits on further marker genes.
#[derive(Default)]
struct Collected {
    placed: Vec<Placed>,
    sets: Vec<Vec<u32>>,
    genomes: HashMap<u32, f64>,
    gene_hits: Vec<PresenceInterval>,
}

/// The genera of a tie set's taxa that have one, by taxid, and the genus
/// when every taxon lies in it.
struct Genera {
    genera: Vec<u32>,
    unique: Option<u32>,
}

/// A candidate row's evidence: its evidence genome (fewest assembly id among
/// genomes holding equally many of its fragments) and its model, the number
/// of fragments there, and their coverage and window spread on it; the
/// genomes holding as many of its fragments (the evidence genome among them),
/// over which its depth is averaged; the length a depth on genome references
/// is taken over: that genome's indexed length, or the harmonic mean of the
/// tied genomes' indexed lengths, so that the depth is the mean of the depths
/// on them; and on a marker genome the whole-genome fragments that each of its
/// fragments there stands for, when the genome declares its size.
#[derive(Clone, Debug)]
pub(super) struct Measured {
    pub(super) genome: u32,
    pub(super) model: Model,
    pub(super) fragments: usize,
    pub(super) coverage: EvidenceCoverage,
    pub(super) spread: Option<f64>,
    pub(super) tied: Arc<[u32]>,
    pub(super) length: f64,
    pub(super) whole_genome_scale: Option<f64>,
}

impl Measured {
    /// The row rules of its evidence model.
    pub(super) fn rules(&self) -> &'static RowRules {
        self.model.rules()
    }
}

/// A test a candidate row fails, with its bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Failure {
    Identity { identity: f64, floor: f64 },
    ChunkRatio { ratio: f64, bound: f64 },
    Windows { spread: Option<f64>, bound: f64 },
    Spill { fragments: usize, bound: f64 },
}

impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Self::Identity { identity, floor } => {
                write!(f, "identity: {identity:.4} < {floor:.4}")
            }
            Self::ChunkRatio { ratio, bound } => write!(f, "chunk ratio: {ratio:.4} < {bound}"),
            Self::Windows { spread: None, .. } => {
                write!(f, "window spread: too few reads to measure")
            }
            Self::Windows {
                spread: Some(spread),
                bound,
            } => write!(f, "window spread: {spread:.4} < {bound}"),
            Self::Spill { fragments, bound } => write!(
                f,
                "spill of kept units beneath: {fragments} fragments < {bound:.1}"
            ),
        }
    }
}

/// The spill test: a row's `fragments` at least the `bound` of the kept units
/// beneath it, when there are any.
fn spill_failure(fragments: usize, bound: Option<f64>) -> Option<Failure> {
    bound
        .filter(|&bound| (fragments as f64) < bound)
        .map(|bound| Failure::Spill { fragments, bound })
}

/// The chunk ratio test: `ratio` at least `bound`.
fn chunk_ratio_failure(ratio: f64, bound: f64) -> Option<Failure> {
    (ratio.is_nan() || ratio < bound).then_some(Failure::ChunkRatio { ratio, bound })
}

/// The tests a genus fails with its unique fragments of one kind: their mean
/// `identity` against `floor`, and their chunk ratio on the evidence genome
/// against the bound of its model.
fn genus_failures(identity: f64, floor: f64, measured: &Measured) -> Vec<Failure> {
    let mut failed = Vec::new();
    if identity < floor {
        failed.push(Failure::Identity { identity, floor });
    }
    failed.extend(chunk_ratio_failure(
        measured.coverage.coebr,
        measured.rules().genus_min_chunk_ratio,
    ));
    failed
}

/// The tests a taxon above genus fails with the fragments that reach it: how
/// they spread over the evidence genome, by its model.
pub(super) fn taxon_failures(measured: &Measured) -> Vec<Failure> {
    match measured.rules().above_genus {
        Spread::Windows(bound) => measured
            .spread
            .is_none_or(|spread| spread < bound)
            .then_some(Failure::Windows {
                spread: measured.spread,
                bound,
            })
            .into_iter()
            .collect(),
        Spread::ChunkRatio(bound) => chunk_ratio_failure(measured.coverage.coebr, bound)
            .into_iter()
            .collect(),
    }
}

/// The evidence of a taxon above genus that passes [`taxon_failures`], in
/// words.
pub(super) fn taxon_evidence(measured: &Measured) -> String {
    match measured.rules().above_genus {
        Spread::Windows(bound) => format!(
            "window spread {:.4} >= {bound}",
            measured.spread.unwrap_or(f64::NAN)
        ),
        Spread::ChunkRatio(bound) => format!(
            "{} genes hit, chunk ratio {:.4} >= {bound}",
            measured.coverage.marker_genes_hit.unwrap_or(0),
            measured.coverage.coebr
        ),
    }
}

/// A genus whose unique fragments of one kind meet the amount of their
/// evidence model: those fragments, the fragments of the kind tied between it
/// and other genera, their mean identity, its floor and their evidence, the
/// tests it fails, and whether its evidence is the one its genus row shows
/// (its passing kind with most unique fragments).
pub(super) struct GenusCandidate {
    pub(super) genus: u32,
    pub(super) kind: Kind,
    pub(super) unique: usize,
    pub(super) tied: usize,
    pub(super) identity: f64,
    pub(super) identity_floor: f64,
    pub(super) measured: Measured,
    /// The spill bound of the kept units beneath the genus, if any, which
    /// its unique fragments of both kinds must reach.
    pub(super) spill: Option<f64>,
    pub(super) failures: Vec<Failure>,
    pub(super) shown: bool,
}

/// A taxon above genus whose fragments meet the amount of their evidence
/// model: those tied across genera, those unique to genera without a genus
/// row, and those moved up from lower ranks; the genera its tied fragments
/// list most often (at most three, most first); their evidence and the tests
/// it fails.
pub(super) struct TaxonCandidate {
    pub(super) taxid: u32,
    pub(super) rank: Rank,
    pub(super) tied: usize,
    pub(super) from_genera: usize,
    pub(super) from_below: usize,
    pub(super) tied_genera: Vec<u32>,
    pub(super) measured: Measured,
    /// The spill bound of the kept units beneath the taxon, if any.
    pub(super) spill: Option<f64>,
    pub(super) failures: Vec<Failure>,
}

impl TaxonCandidate {
    pub(super) fn fragments(&self) -> usize {
        self.tied + self.from_genera + self.from_below
    }
}

/// A row above species: its fragments (a pair counting once), of a genus row
/// those unique to the genus by kind, its depth (see the [module
/// documentation](self)), and its evidence.
pub(super) struct HigherRow {
    pub(super) taxid: u32,
    pub(super) rank: Rank,
    pub(super) reads: f64,
    pub(super) unique: [usize; 2],
    pub(super) depth: f64,
    pub(super) measured: Measured,
}

/// What the levels above species find.
pub(super) struct HigherLevels {
    /// Rows by rank, then taxid.
    pub(super) rows: Vec<HigherRow>,
    /// Genus candidates by genus, then kind.
    pub(super) genera: Vec<GenusCandidate>,
    /// Candidates above genus by rank, then taxid.
    pub(super) taxa: Vec<TaxonCandidate>,
    /// Fragments of each kind, and those left unexplained.
    pub(super) fragments: [usize; 2],
    pub(super) unexplained: usize,
    /// With `keep_shares`, each placed fragment's read and its shares.
    pub(super) shares: FragmentShares,
    /// Where the fragments of each taxon without a row moved: the next taxon
    /// above it, or `None` (unexplained).
    pub(super) moved: HashMap<u32, Option<u32>>,
    /// Where the unique fragments of each genus without a genus row of their
    /// kind went: the first taxon above it, or `None` (unexplained).
    pub(super) genus_start: HashMap<u32, Option<u32>>,
}

impl HigherLevels {
    /// The row that took the fragments starting at `taxid`, if any.
    pub(super) fn end_of(&self, mut taxid: u32) -> Option<&HigherRow> {
        loop {
            if let Some(row) = self.row(taxid) {
                return Some(row);
            }
            taxid = (*self.moved.get(&taxid)?)?;
        }
    }

    /// The row that took the unique fragments of `genus` failing its tests.
    pub(super) fn end_of_genus(&self, genus: u32) -> Option<&HigherRow> {
        self.end_of((*self.genus_start.get(&genus)?)?)
    }

    pub(super) fn row(&self, taxid: u32) -> Option<&HigherRow> {
        self.rows.iter().find(|row| row.taxid == taxid)
    }
}

/// Where a fragment went.
#[derive(Clone, Copy)]
enum Destination {
    Unexplained,
    Taxon(u32),
    /// Shared by the posteriors of one of step 2's read groups.
    Shared(u32),
}

/// The fragments that reach a taxon above genus, by where they come from.
#[derive(Default)]
struct Arrivals {
    fragments: Vec<u32>,
    tied: usize,
    from_genera: usize,
    from_below: usize,
    /// Tied fragments per tie set.
    sets: HashMap<u32, usize>,
}

/// Place the level's fragments ([module documentation](self)), with the
/// taxa's `lineages`.
pub(super) fn place(level: &Level, lineages: &Lineages) -> HigherLevels {
    let Collected {
        placed,
        sets,
        genomes,
        gene_hits,
    } = collect(level);
    let genera = sets
        .iter()
        .map(|taxa| genera_of(taxa, lineages))
        .collect::<Vec<_>>();
    let mut placing = Placing {
        level,
        lineages,
        placed: &placed,
        gene_hits: &gene_hits,
        sets: &sets,
        genera: &genera,
        genomes: &genomes,
        destination: vec![Destination::Unexplained; placed.len()],
        pending: HashMap::new(),
        unexplained: 0,
    };
    let GenusRows {
        candidates,
        with_row,
        genus_start,
    } = placing.genus_rows();
    let tied_shares = placing.share_across_genera(&with_row);
    let AboveGenus {
        taxa,
        mut evidence,
        moved,
    } = placing.rows_above_genus();
    for candidate in candidates.iter().filter(|candidate| candidate.shown) {
        evidence.insert((Rank::Genus, candidate.genus), candidate.measured.clone());
    }
    let (rows, shares) = placing.rows(evidence, &tied_shares, &with_row);
    let mut fragments = [0; 2];
    for fragment in &placed {
        fragments[fragment.kind as usize] += 1;
    }
    HigherLevels {
        rows,
        genera: candidates,
        taxa,
        fragments,
        unexplained: placing.unexplained,
        shares,
        moved,
        genus_start,
    }
}

/// Each placed fragment's read and its shares of the rows' taxa.
pub(super) type FragmentShares = Vec<(usize, Vec<(u32, f64)>)>;

/// Step 1's genus rows: the genus candidates, by genus and kind; the genera
/// with a genus row, with their unique fragments of each kind; and where the
/// unique fragments of each other genus went: the first taxon above it, or
/// none (unexplained).
struct GenusRows {
    candidates: Vec<GenusCandidate>,
    with_row: BTreeMap<u32, [usize; 2]>,
    genus_start: HashMap<u32, Option<u32>>,
}

/// Step 3's rows above genus: the candidates, by rank and taxid; the
/// evidence of each taxon with a row; and where the fragments of each taxon
/// without a row moved: the next taxon above it, or none (unexplained).
struct AboveGenus {
    taxa: Vec<TaxonCandidate>,
    evidence: BTreeMap<(Rank, u32), Measured>,
    moved: HashMap<u32, Option<u32>>,
}

/// The placing of the level's fragments: where each fragment went, the
/// fragments waiting at a taxon above genus (with its rank), and the
/// fragments left unexplained.
struct Placing<'p> {
    level: &'p Level<'p>,
    lineages: &'p Lineages<'p>,
    placed: &'p [Placed],
    gene_hits: &'p [PresenceInterval],
    sets: &'p [Vec<u32>],
    genera: &'p [Genera],
    genomes: &'p HashMap<u32, f64>,
    destination: Vec<Destination>,
    pending: HashMap<u32, (Rank, Arrivals)>,
    unexplained: usize,
}

impl Placing<'_> {
    /// The evidence of the fragments `list` when they meet the amount of
    /// their evidence model ([`measure`]).
    fn measure(&self, list: &[u32]) -> Option<Measured> {
        measure(self, list)
    }

    /// The hits of `fragment` on further marker genes.
    fn hits_of(&self, fragment: &Placed) -> &[PresenceInterval] {
        &self.gene_hits[fragment.gene_hits.start as usize..fragment.gene_hits.end as usize]
    }

    /// The depth that `fragment` adds to a row on marker genomes `tied`: the
    /// mean over them of its depth on each ([`MarkerModel::read_depth`]),
    /// mapped bases over the indexed length where no gene can hold its read.
    fn marker_depth(&self, markers: &MarkerModel, fragment: &Placed, tied: &[u32]) -> f64 {
        let genes = 1 + fragment.gene_hits.len() as u32;
        tied.iter()
            .map(|&genome| {
                markers
                    .read_depth(genome, fragment.read_len, genes)
                    .unwrap_or_else(|| f64::from(fragment.bases) / self.genomes[&genome].max(1.0))
            })
            .sum::<f64>()
            / tied.len() as f64
    }

    /// The spill bound of taxon `taxid` for a row measured as `measured`.
    fn spill(&self, taxid: u32, measured: &Measured) -> Option<f64> {
        self.level
            .kept_beneath
            .get(&taxid)
            .map(|beneath| beneath.bound(measured.model))
    }

    /// The fragments waiting at `taxid` of `rank`.
    fn arrivals(&mut self, taxid: u32, rank: Rank) -> &mut Arrivals {
        &mut self
            .pending
            .entry(taxid)
            .or_insert_with(|| (rank, Arrivals::default()))
            .1
    }

    /// Step 1: genus rows from unique fragments, by kind ([`GenusRows`]).
    fn genus_rows(&mut self) -> GenusRows {
        let (placed, genera) = (self.placed, self.genera);
        let mut unique: BTreeMap<(u32, Kind), Vec<u32>> = BTreeMap::new();
        let mut tied: HashMap<(u32, Kind), usize> = HashMap::new();
        for (i, fragment) in placed.iter().enumerate() {
            let of = &genera[fragment.ties as usize];
            match of.unique {
                Some(genus) => unique
                    .entry((genus, fragment.kind))
                    .or_default()
                    .push(i as u32),
                None => {
                    for &genus in &of.genera {
                        *tied.entry((genus, fragment.kind)).or_default() += 1;
                    }
                }
            }
        }
        let mut of_genus: HashMap<u32, usize> = HashMap::new();
        for (&(genus, _), list) in &unique {
            *of_genus.entry(genus).or_default() += list.len();
        }
        let cutoff = self.level.not_kept.1.identity;
        let mut candidates = Vec::new();
        for (&(genus, kind), list) in &unique {
            let Some(measured) = self.measure(list) else {
                continue;
            };
            let identity = list
                .iter()
                .map(|&i| f64::from(placed[i as usize].identity))
                .sum::<f64>()
                / list.len() as f64;
            let identity_floor = cutoff - measured.rules().genus_identity_margin;
            let spill = self.spill(genus, &measured);
            let mut failures = genus_failures(identity, identity_floor, &measured);
            failures.extend(spill_failure(of_genus[&genus], spill));
            candidates.push(GenusCandidate {
                genus,
                kind,
                unique: list.len(),
                tied: tied.get(&(genus, kind)).copied().unwrap_or(0),
                identity,
                identity_floor,
                measured,
                spill,
                failures,
                shown: false,
            });
        }
        // A genus with a row from either kind holds the fragments of both
        // kinds unique to it; those of the others go to the first taxon above
        // their genus at family, order, class, phylum or domain rank.
        let mut with_row: BTreeMap<u32, [usize; 2]> = candidates
            .iter()
            .filter(|candidate| candidate.failures.is_empty())
            .map(|candidate| (candidate.genus, [0; 2]))
            .collect();
        let mut genus_start: HashMap<u32, Option<u32>> = HashMap::new();
        for ((genus, kind), list) in unique {
            if let Some(unique) = with_row.get_mut(&genus) {
                unique[kind as usize] = list.len();
                for &i in &list {
                    self.destination[i as usize] = Destination::Taxon(genus);
                }
                continue;
            }
            let start = self.lineages.above(genus, Rank::Genus);
            genus_start.insert(genus, start.map(|(taxid, _)| taxid));
            match start {
                Some((start, rank)) => {
                    let arrivals = self.arrivals(start, rank);
                    arrivals.from_genera += list.len();
                    arrivals.fragments.extend(list);
                }
                None => self.unexplained += list.len(),
            }
        }
        // The genus row shows the evidence of its kind with most unique
        // fragments.
        let mut shown: HashMap<u32, usize> = HashMap::new();
        for (index, candidate) in candidates.iter().enumerate() {
            if candidate.failures.is_empty() {
                let best = shown.entry(candidate.genus).or_insert(index);
                if candidates[*best].unique < candidate.unique {
                    *best = index;
                }
            }
        }
        for &index in shown.values() {
            candidates[index].shown = true;
        }
        GenusRows {
            candidates,
            with_row,
            genus_start,
        }
    }

    /// Step 2: fragments tied across genera go to the genera with a genus row
    /// (`with_row`) when their tie set holds one, shared by the EM
    /// ([`share_tied`]), and otherwise wait at the first taxon above genus of
    /// their tie set ([`Lineages::common_start`]). Returns the posteriors of
    /// the EM's read groups.
    fn share_across_genera(
        &mut self,
        with_row: &BTreeMap<u32, [usize; 2]>,
    ) -> Vec<Vec<(u32, f64)>> {
        let (placed, genera, sets, lineages) = (self.placed, self.genera, self.sets, self.lineages);
        let mut tied: Vec<u32> = Vec::new();
        let mut starts: HashMap<u32, Option<(u32, Rank)>> = HashMap::new();
        for (i, fragment) in placed.iter().enumerate() {
            let of = &genera[fragment.ties as usize];
            if of.unique.is_some() {
                continue;
            }
            if of.genera.iter().any(|genus| with_row.contains_key(genus)) {
                tied.push(i as u32);
                continue;
            }
            let start = *starts
                .entry(fragment.ties)
                .or_insert_with(|| lineages.common_start(&sets[fragment.ties as usize]));
            match start {
                Some((start, rank)) => {
                    let arrivals = self.arrivals(start, rank);
                    arrivals.tied += 1;
                    *arrivals.sets.entry(fragment.ties).or_default() += 1;
                    arrivals.fragments.push(i as u32);
                }
                None => self.unexplained += 1,
            }
        }
        if tied.is_empty() {
            return Vec::new();
        }
        let (em, groups) = share_tied(self.level, placed, &tied, with_row, lineages);
        for (i, group) in groups {
            let group = group.expect("a fragment with a genus row among its tie set's genera");
            self.destination[i as usize] = Destination::Shared(group as u32);
        }
        em.posteriors
    }

    /// Step 3: rows above genus, rank by rank ([`AboveGenus`]).
    fn rows_above_genus(&mut self) -> AboveGenus {
        let mut taxa = Vec::new();
        let mut evidence: BTreeMap<(Rank, u32), Measured> = BTreeMap::new();
        let mut moved: HashMap<u32, Option<u32>> = HashMap::new();
        for rank in [
            Rank::Family,
            Rank::Order,
            Rank::Class,
            Rank::Phylum,
            Rank::Domain,
        ] {
            let mut here = self
                .pending
                .iter()
                .filter(|(_, (of, _))| *of == rank)
                .map(|(&taxid, _)| taxid)
                .collect::<Vec<_>>();
            here.sort_unstable();
            for taxid in here {
                let (_, arrivals) = self.pending.remove(&taxid).expect("listed above");
                let passed = self.measure(&arrivals.fragments).is_some_and(|measured| {
                    let spill = self.spill(taxid, &measured);
                    let mut failures = taxon_failures(&measured);
                    failures.extend(spill_failure(arrivals.fragments.len(), spill));
                    let passed = failures.is_empty();
                    if passed {
                        evidence.insert((rank, taxid), measured.clone());
                    }
                    taxa.push(TaxonCandidate {
                        taxid,
                        rank,
                        tied: arrivals.tied,
                        from_genera: arrivals.from_genera,
                        from_below: arrivals.from_below,
                        tied_genera: most_listed(&arrivals.sets, self.genera),
                        measured,
                        spill,
                        failures,
                    });
                    passed
                });
                if passed {
                    for &i in &arrivals.fragments {
                        self.destination[i as usize] = Destination::Taxon(taxid);
                    }
                    continue;
                }
                let next = self.lineages.above(taxid, rank);
                moved.insert(taxid, next.map(|(next, _)| next));
                match next {
                    Some((next, next_rank)) => {
                        let above = self.arrivals(next, next_rank);
                        above.from_below += arrivals.fragments.len();
                        above.fragments.extend(arrivals.fragments);
                    }
                    None => self.unexplained += arrivals.fragments.len(),
                }
            }
        }
        debug_assert!(self.pending.is_empty());
        AboveGenus {
            taxa,
            evidence,
            moved,
        }
    }

    /// The rows of the taxa with `evidence`, each fragment counted where it
    /// went (`tied_shares`: the posteriors of step 2's read groups; `with_row`:
    /// each genus row's unique fragments by kind), and with
    /// [`Level::keep_shares`], each placed fragment's read and shares.
    fn rows(
        &self,
        evidence: BTreeMap<(Rank, u32), Measured>,
        tied_shares: &[Vec<(u32, f64)>],
        with_row: &BTreeMap<u32, [usize; 2]>,
    ) -> (Vec<HigherRow>, FragmentShares) {
        let on_markers = self.level.data.marker_model.as_deref().map(|markers| {
            let rows = evidence
                .iter()
                .filter(|(_, measured)| measured.model == Model::Markers)
                .map(|(&(_, taxid), measured)| (taxid, &measured.tied[..]))
                .collect::<HashMap<_, _>>();
            (markers, rows)
        });
        // Reads, first-read bases, and depth on marker genomes.
        let mut totals: HashMap<u32, [f64; 3]> = HashMap::new();
        let mut shares = Vec::new();
        for (fragment, to) in self.placed.iter().zip(&self.destination) {
            let whole;
            let parts: &[(u32, f64)] = match *to {
                Destination::Unexplained => continue,
                Destination::Taxon(taxid) => {
                    whole = [(taxid, 1.0)];
                    &whole
                }
                Destination::Shared(index) => &tied_shares[index as usize],
            };
            for &(taxid, share) in parts {
                let total = totals.entry(taxid).or_default();
                total[0] += share;
                total[1] += share * f64::from(fragment.bases);
                if let Some((markers, rows)) = &on_markers {
                    if let Some(tied) = rows.get(&taxid) {
                        total[2] += share * self.marker_depth(markers, fragment, tied);
                    }
                }
            }
            if self.level.keep_shares {
                shares.push((fragment.read, parts.to_vec()));
            }
        }
        let rows = evidence
            .into_iter()
            .map(|((rank, taxid), measured)| {
                let [reads, bases, marker_depth] = totals[&taxid];
                HigherRow {
                    taxid,
                    rank,
                    reads,
                    unique: with_row.get(&taxid).copied().unwrap_or_default(),
                    depth: match measured.model {
                        Model::Genome => bases / measured.length,
                        Model::Markers => marker_depth,
                    },
                    measured,
                }
            })
            .collect();
        (rows, shares)
    }
}

/// Step 2's shares of the fragments `list`, tied across genera: the
/// likelihood EM over the genera with a genus row (`with_row`, with their
/// unique fragments of both kinds, which count for them alone). Each fragment
/// lists those genera among the genera of its candidates under the gates of
/// its kind, each with the fewest errors of its taxa (tied candidates carrying
/// the best's errors), as the species level's EM lists kept species. The EM,
/// and each fragment with its read group.
fn share_tied(
    level: &Level,
    placed: &[Placed],
    list: &[u32],
    with_row: &BTreeMap<u32, [usize; 2]>,
    lineages: &Lineages,
) -> (likelihood_em::EmResult, Vec<(u32, Option<usize>)>) {
    let mut fragments = Vec::new();
    for (&genus, unique) in with_row {
        fragments.extend((0..unique.iter().sum::<usize>()).map(|_| Fragment {
            read: 0,
            candidates: vec![(genus, 0.0)],
            bases: 0.0,
        }));
    }
    let first_tied = fragments.len();
    let mut order = Vec::with_capacity(list.len());
    for kind in [Kind::BelowSpecies, Kind::NotKept] {
        let policy = level.policy(kind);
        let of_kind = list
            .iter()
            .copied()
            .filter(|&i| placed[i as usize].kind == kind)
            .collect::<Vec<_>>();
        fold_candidates(
            level.data,
            policy,
            &of_kind,
            |&i| placed[i as usize].read,
            |_, candidates| ranked(candidates, policy.mode).map(|found| found.0),
            |&i, found| {
                let fragment = &placed[i as usize];
                let candidates = found
                    .iter()
                    .flatten()
                    .filter_map(|&(taxid, errors)| {
                        lineages
                            .at_rank(taxid, Rank::Genus)
                            .filter(|genus| with_row.contains_key(genus))
                            .map(|genus| (genus, errors))
                    })
                    .collect();
                order.push(i);
                fragments.push(Fragment {
                    read: fragment.read,
                    candidates,
                    bases: f64::from(fragment.bases),
                });
            },
        );
    }
    let shared = share_reads(&fragments, |_| true, level.em.penalty);
    let em = level.em.run(&shared.groups);
    let groups = order
        .into_iter()
        .zip(shared.group_of[first_tied..].iter().copied())
        .collect();
    (em, groups)
}

/// Build the level's fragments: each read of both kinds with its tie set
/// ([`tie_set`]) and best placement ([`ranked`]).
fn collect(level: &Level) -> Collected {
    let mut index: HashMap<Vec<u32>, u32> = HashMap::new();
    let mut found = Collected::default();
    for (kind, (reads, policy)) in [
        (Kind::BelowSpecies, level.below_species),
        (Kind::NotKept, level.not_kept),
    ] {
        fold_candidates(
            level.data,
            policy,
            reads,
            |&read| read,
            |_, candidates| {
                ranked(candidates, policy.mode).map(|(ranked, best)| {
                    let hits = best
                        .other_gene_hits
                        .iter()
                        .filter(|hit| hit.assembly_id == best.best.assembly_id)
                        .map(|hit| hit.interval)
                        .collect::<Vec<_>>();
                    (tie_set(&ranked), best.best, hits)
                })
            },
            |&read, result| {
                let Some((taxa, best, hits)) = result else {
                    return;
                };
                let next = index.len() as u32;
                let ties = *index.entry(taxa).or_insert(next);
                found
                    .genomes
                    .entry(best.assembly_id)
                    .or_insert(best.asm_len);
                let first_hit = found.gene_hits.len() as u32;
                found.gene_hits.extend(hits);
                found.placed.push(Placed {
                    read,
                    kind,
                    ties,
                    assembly: best.assembly_id,
                    identity: placement_identity(&best) as f32,
                    bases: best.r1_mapped as f32,
                    read_len: best.read_len,
                    interval: PresenceInterval::from(&best),
                    gene_hits: first_hit..found.gene_hits.len() as u32,
                });
            },
        );
    }
    found.sets = vec![Vec::new(); index.len()];
    for (taxa, id) in index {
        found.sets[id as usize] = taxa;
    }
    found
}

/// The taxa of a fragment's tie set: those of its `ranked` candidates tied
/// with the best ([`best_tied`]), by taxid.
fn tie_set(ranked: &[(u32, f64)]) -> Vec<u32> {
    best_tied(ranked)
        .iter()
        .map(|candidate| candidate.0)
        .collect()
}

/// The genera of a tie set's `taxa`.
fn genera_of(taxa: &[u32], lineages: &Lineages) -> Genera {
    let of = taxa
        .iter()
        .map(|&taxid| lineages.at_rank(taxid, Rank::Genus))
        .collect::<Vec<_>>();
    let mut genera = of.iter().flatten().copied().collect::<Vec<_>>();
    genera.sort_unstable();
    genera.dedup();
    let unique = match genera[..] {
        [genus] if of.iter().all(Option::is_some) => Some(genus),
        _ => None,
    };
    Genera { genera, unique }
}

/// The evidence of the fragments `list` ([`Measured`]) when they meet the
/// amount of their evidence model, otherwise `None`.
fn measure(placing: &Placing, list: &[u32]) -> Option<Measured> {
    let (placed, genomes, level) = (placing.placed, placing.genomes, placing.level);
    let mut counts: HashMap<u32, usize> = HashMap::new();
    for &i in list {
        *counts.entry(placed[i as usize].assembly).or_default() += 1;
    }
    let (&genome, &fragments) = counts
        .iter()
        .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
        .expect("a candidate row has fragments");
    let model = Model::of(level.data, genome);
    if list.len() < model.rules().min_fragments {
        return None;
    }
    let mut tied = counts
        .iter()
        .filter(|&(_, &n)| n == fragments)
        .map(|(&id, _)| id)
        .collect::<Vec<_>>();
    tied.sort_unstable();
    // The indexed lengths of the tied genomes; several are summed in sorted
    // order, which the map's order cannot change.
    let mut lengths = tied
        .iter()
        .map(|id| genomes[id].max(1.0))
        .collect::<Vec<_>>();
    let length = match lengths[..] {
        [length] => length,
        _ => {
            lengths.sort_unstable_by(f64::total_cmp);
            lengths.len() as f64 / lengths.iter().map(|length| 1.0 / length).sum::<f64>()
        }
    };
    let mut assembly = PresenceAssembly {
        asm_len: genomes[&genome],
        count: fragments as u64,
        intervals: Vec::with_capacity(fragments),
        ..PresenceAssembly::default()
    };
    let on_genome = list
        .iter()
        .map(|&i| &placed[i as usize])
        .filter(|fragment| fragment.assembly == genome)
        .collect::<Vec<_>>();
    for fragment in &on_genome {
        assembly.add_interval(fragment.interval);
    }
    let spread = level.layout.window_spread(genome, &assembly.intervals);
    // A fragment counts on every further marker gene it covers.
    for fragment in &on_genome {
        let hits = placing.hits_of(fragment);
        assembly.other_gene_hits += hits.len() as u64;
        for &hit in hits {
            assembly.add_interval(hit);
        }
    }
    let coverage = level.tests.measure(genome, &mut assembly);
    if coverage.marker_genes_hit.unwrap_or(0) < model.rules().min_genes {
        return None;
    }
    let whole_genome_scale = level
        .data
        .marker_model
        .as_deref()
        .filter(|_| model == Model::Markers)
        .and_then(|markers| whole_genome_scale(placing, markers, genome, &on_genome));
    Some(Measured {
        genome,
        model,
        fragments,
        coverage,
        spread,
        tied: tied.into(),
        length,
        whole_genome_scale,
    })
}

/// The whole-genome fragments that a fragment on marker genome `genome`
/// stands for, from the fragments `on_genome` placed there, as a species
/// reported on the genome counts its fragments in the unclassified estimate:
/// for short reads the genome's scale times their gene hits per fragment, for
/// long reads the whole-genome reads of each read over their number. `None`
/// without a declared genome size or a gene that can hold a read.
fn whole_genome_scale(
    placing: &Placing,
    markers: &MarkerModel,
    genome: u32,
    on_genome: &[&Placed],
) -> Option<f64> {
    let size = markers.declared_size(genome)?;
    if !markers.long_reads() {
        let hits = on_genome
            .iter()
            .map(|fragment| 1 + fragment.gene_hits.len())
            .sum::<usize>();
        return Some(markers.genome_scale(genome, size)? * hits as f64 / on_genome.len() as f64);
    }
    let mut reads = GenomeReads::default();
    for fragment in on_genome {
        let genes = std::iter::once(fragment.interval.reference_id)
            .chain(placing.hits_of(fragment).iter().map(|hit| hit.reference_id));
        if let Some((hits, per_hit)) = markers.whole_genome_reads(genome, fragment.read_len, genes)
        {
            reads.add(hits, per_hit);
        }
    }
    reads.whole_genome_per_read(size)
}

/// The genera that the tied fragments of the tie sets `sets` (each with its
/// fragments) list most often, at most three, most first, then by taxid.
fn most_listed(sets: &HashMap<u32, usize>, genera: &[Genera]) -> Vec<u32> {
    let mut listed: HashMap<u32, usize> = HashMap::new();
    for (&set, &fragments) in sets {
        for &genus in &genera[set as usize].genera {
            *listed.entry(genus).or_default() += fragments;
        }
    }
    let mut listed = listed.into_iter().collect::<Vec<_>>();
    listed.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    listed.into_iter().take(3).map(|(genus, _)| genus).collect()
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::super::test_fixtures::{data, row};
    use super::*;
    use crate::profiler::evidence::{CandidateScratch, EvidenceAlignment, EvidenceMode};
    use crate::profiler::markers::{FragmentGeometry, FragmentLayout, MarkerModel};

    /// Domain 2 (phylum 10, class 20): order 30 holds family 40 (genera 100 and 101) and family 41 (genus 102);
    /// order 31 holds family 42 (genus 103). Species 1000 + i lie in genus 100 + i (species 1001 and 1002 in genus
    /// 100); species 5000 lies in domain 5 (no lower rank).
    fn taxonomy() -> Taxonomy {
        let node = |parent: u32, rank: &str| (parent, rank.to_string());
        Taxonomy::from_parts(
            HashMap::from([
                (1, node(1, "no rank")),
                (2, node(1, "domain")),
                (5, node(1, "domain")),
                (10, node(2, "phylum")),
                (20, node(10, "class")),
                (30, node(20, "order")),
                (31, node(20, "order")),
                (40, node(30, "family")),
                (41, node(30, "family")),
                (42, node(31, "family")),
                (100, node(40, "genus")),
                (101, node(40, "genus")),
                (102, node(41, "genus")),
                (103, node(42, "genus")),
                (1001, node(100, "species")),
                (1002, node(100, "species")),
                (1011, node(101, "species")),
                (1021, node(102, "species")),
                (1031, node(103, "species")),
                (5000, node(5, "species")),
            ]),
            HashMap::new(),
            HashMap::new(),
            HashSet::new(),
        )
    }

    /// A read of 100 bases at `start` of every genome in `on` (species, genome, identity).
    fn read(on: &[(u32, u32, f32)], start: u32) -> Vec<EvidenceAlignment> {
        on.iter()
            .map(|&(taxid, genome, identity)| EvidenceAlignment {
                span: (start, start + 100),
                fixed_chunk: (start, start + 100),
                flex_chunk: (start.saturating_sub(950), start + 1_050),
                ..row(taxid, genome, identity)
            })
            .collect()
    }

    /// A read as [`read`], with genome 5 of 100 kb.
    fn read_small(on: &[(u32, u32, f32)], start: u32) -> Vec<EvidenceAlignment> {
        read(on, start)
            .into_iter()
            .map(|a| match a.assembly_id {
                5 => EvidenceAlignment {
                    asm_len: 100_000.0,
                    ..a
                },
                _ => a,
            })
            .collect()
    }

    /// `n` reads on `on`, spread evenly along its 1 Mb genomes, or clustered in 10 kb.
    fn reads(on: &[(u32, u32, f32)], n: u32, spread: bool) -> Vec<Vec<EvidenceAlignment>> {
        (0..n)
            .map(|i| {
                read(
                    on,
                    if spread {
                        i * (990_000 / n)
                    } else {
                        500_000 + i * 25
                    },
                )
            })
            .collect()
    }

    fn policy(identity: f64) -> EvidencePolicy {
        EvidencePolicy {
            mode: EvidenceMode::Single,
            identity,
            mapped_len: 50,
            fraction: 0.6,
        }
    }

    /// The level over `below` (reads without a species-level match, at the floor identity 0.80) and `not_kept`
    /// (reads of species not kept, at the species identity 0.95).
    fn place_reads(
        below: Vec<Vec<EvidenceAlignment>>,
        not_kept: Vec<Vec<EvidenceAlignment>>,
    ) -> HigherLevels {
        let below_reads = (0..below.len()).collect::<Vec<_>>();
        let not_kept_reads = (below.len()..below.len() + not_kept.len()).collect::<Vec<_>>();
        let data = data(below.into_iter().chain(not_kept).collect());
        let tests = CoverageTests {
            multiple_tests: 1,
            chunk_cutoff: 0.0,
            ratio_bounds: (0.6, 0.75),
            markers: None,
            min_marker_genes: 0,
        };
        // Genome 5 has 100 kb, the others 1 Mb.
        let layout =
            GenomeLayout::new((0..=60).map(|id| (id, if id == 5 { 100_000 } else { 1_000_000 })));
        let taxonomy = taxonomy();
        place(
            &Level {
                data: &data,
                tests: &tests,
                layout: &layout,
                below_species: (&below_reads, policy(0.80)),
                not_kept: (&not_kept_reads, policy(0.95)),
                em: EmSettings {
                    penalty: 3.0,
                    max_iterations: 1000,
                },
                keep_shares: true,
                kept_beneath: &HashMap::new(),
            },
            &Lineages::new(&taxonomy),
        )
    }

    /// Marker genome 60 of species 1001 (genus 100): 40 genes of 1,000 bases, reference ids 1 to 40.
    fn marker_model() -> Arc<MarkerModel> {
        Arc::new(MarkerModel::from_gene_lengths(
            &[(60, &[1_000; 40])],
            FragmentGeometry::new(FragmentLayout::Reads, 100, 100, 50, 0.6),
        ))
    }

    /// A read of 100 bases of species 1001 at `start` of gene `gene` of marker genome 60.
    fn marker_read(gene: u32, start: u32, identity: f32) -> Vec<EvidenceAlignment> {
        vec![EvidenceAlignment {
            reference_id: gene,
            asm_len: 40_000.0,
            span: (start, start + 100),
            fixed_chunk: (start, start + 100),
            flex_chunk: (start, start + 100),
            marker_gene_len: 1_000,
            ..row(1001, 60, identity)
        }]
    }

    /// `per_gene` reads of species 1001 on each of the first `genes` genes of marker genome 60.
    fn marker_reads(genes: u32, per_gene: u32) -> Vec<Vec<EvidenceAlignment>> {
        (1..=genes)
            .flat_map(|gene| (0..per_gene).map(move |i| marker_read(gene, 50 + i * 150, 0.99)))
            .collect()
    }

    /// The level over reads of species not kept on marker genome 60.
    fn place_marker_reads(not_kept: Vec<Vec<EvidenceAlignment>>) -> HigherLevels {
        place_marker_reads_beside(not_kept, &HashMap::new())
    }

    /// [`place_marker_reads`] with the spill bounds `kept_beneath`.
    fn place_marker_reads_beside(
        not_kept: Vec<Vec<EvidenceAlignment>>,
        kept_beneath: &HashMap<u32, Beneath>,
    ) -> HigherLevels {
        place_marker_reads_with(not_kept, kept_beneath, marker_model())
    }

    /// [`place_marker_reads_beside`] with the marker model `model`.
    fn place_marker_reads_with(
        not_kept: Vec<Vec<EvidenceAlignment>>,
        kept_beneath: &HashMap<u32, Beneath>,
        model: Arc<MarkerModel>,
    ) -> HigherLevels {
        let reads = (0..not_kept.len()).collect::<Vec<_>>();
        let mut data = data(not_kept);
        data.marker_model = Some(Arc::clone(&model));
        let tests = CoverageTests {
            multiple_tests: 1,
            chunk_cutoff: 0.0,
            ratio_bounds: (0.6, 0.75),
            markers: Some(model),
            min_marker_genes: 3,
        };
        let layout =
            GenomeLayout::new(std::iter::once((0, 1)).chain((1..=40).map(|_| (60, 1_000))));
        let taxonomy = taxonomy();
        place(
            &Level {
                data: &data,
                tests: &tests,
                layout: &layout,
                below_species: (&[], policy(0.80)),
                not_kept: (&reads, policy(0.95)),
                em: EmSettings {
                    penalty: 3.0,
                    max_iterations: 1000,
                },
                keep_shares: false,
                kept_beneath,
            },
            &Lineages::new(&taxonomy),
        )
    }

    #[test]
    fn a_row_on_marker_genes_has_the_depth_and_whole_genome_scale_of_a_species_on_them() {
        // 72 reads of 100 bases on genome 60, whose 40 genes of 1,000 bases offer 40 x (1,000 - 2 x 60 + 101)
        // placements to a read that needs 60 of its bases inside a gene.
        let placements = 40.0 * (1_000.0 - 2.0 * 60.0 + 101.0);
        let declared = |model: Arc<MarkerModel>| {
            Arc::new(
                Arc::into_inner(model)
                    .expect("one owner")
                    .with_declared(HashMap::from([(60, 2_000_000)])),
            )
        };
        let found = place_marker_reads_with(
            marker_reads(24, 3),
            &HashMap::new(),
            declared(marker_model()),
        );
        let row = &found.rows[0];
        assert!((row.depth - 72.0 * 100.0 / placements).abs() < 1e-12);
        assert_eq!(&row.measured.tied[..], [60]);
        // In the unclassified estimate each read stands for the genome's placements over those of its genes.
        let scale = (2_000_000.0 - 120.0 + 101.0) / placements;
        assert!((row.measured.whole_genome_scale.unwrap() - scale).abs() < 1e-9);
        // Without a declared genome size there is no scale.
        let found = place_marker_reads(marker_reads(24, 3));
        assert_eq!(found.rows[0].measured.whole_genome_scale, None);
    }

    #[test]
    fn kept_units_on_marker_genes_bound_the_rows_of_the_taxa_above_them() {
        let taxonomy = taxonomy();
        let lineages = Lineages::new(&taxonomy);
        // Species 1001 (genus 100) on marker genes with 400 fragments, each standing for 20 whole-genome
        // fragments; species 1011 (genus 101) on a genome reference leaves no spill.
        let beneath = kept_beneath(
            [
                (1001, Model::Markers, 400.0, Some(20.0)),
                (1011, Model::Genome, 1_000.0, None),
            ],
            &lineages,
        );
        let share = MARKER_ROWS.spill_share;
        assert_eq!(GENOME_ROWS.spill_share, 0.0);
        let expected = Beneath {
            marker_fragments: share * 400.0,
            whole_genome: share * 8_000.0,
        };
        for taxon in [100, 40, 30, 20, 10, 2] {
            assert_eq!(beneath[&taxon], expected, "{taxon}");
        }
        assert!(!beneath.contains_key(&101) && !beneath.contains_key(&1001));
        // A row on marker genes compares with their fragments, a row on a genome reference with their whole-genome
        // equivalents.
        assert_eq!(expected.bound(Model::Markers), share * 400.0);
        assert_eq!(expected.bound(Model::Genome), share * 8_000.0);
    }

    #[test]
    fn a_genus_row_beside_kept_marker_units_needs_more_than_their_spill() {
        // 72 reads over 24 genes of genus 100: a row alone, but not beside kept units spilling 100 fragments.
        let bound = |marker_fragments| {
            HashMap::from([(
                100,
                Beneath {
                    marker_fragments,
                    whole_genome: 50.0 * marker_fragments,
                },
            )])
        };
        let found = place_marker_reads_beside(marker_reads(24, 3), &bound(100.0));
        // The genus fails; its reads move up to the family, which has no kept units beneath here.
        assert_eq!(rows(&found), [(40, Rank::Family, 72.0)]);
        let candidate = &found.genera[0];
        assert_eq!(candidate.spill, Some(100.0));
        assert_eq!(
            candidate.failures,
            [Failure::Spill {
                fragments: 72,
                bound: 100.0
            }]
        );
        assert_eq!(
            candidate.failures[0].to_string(),
            "spill of kept units beneath: 72 fragments < 100.0"
        );
        let found = place_marker_reads_beside(marker_reads(24, 3), &bound(72.0));
        assert_eq!(rows(&found), [(100, Rank::Genus, 72.0)]);
    }

    #[test]
    fn reads_on_marker_genes_make_a_genus_row_from_ten_genes_hit() {
        // 72 reads over 24 of the 40 genes: far fewer than 200 reads, but 24 genes hit, evenly.
        let found = place_marker_reads(marker_reads(24, 3));
        assert_eq!(rows(&found), [(100, Rank::Genus, 72.0)]);
        let candidate = &found.genera[0];
        assert_eq!(candidate.measured.model, Model::Markers);
        assert_eq!(candidate.measured.coverage.marker_genes_hit, Some(24));
        assert!(candidate.measured.coverage.coebr >= MARKER_ROWS.genus_min_chunk_ratio);
        // Reads on 6 genes: too few genes for a candidate at the genus or above it; they stay unexplained.
        let found = place_marker_reads(marker_reads(6, 12));
        assert!(found.rows.is_empty() && found.genera.is_empty() && found.taxa.is_empty());
        assert_eq!(found.unexplained, 72);
        // 400 reads piled on 12 genes: enough genes, but their chunk ratio fails; the family above fails it too.
        let found = place_marker_reads(marker_reads(12, 34));
        assert!(found.rows.is_empty());
        assert!(matches!(
            found.genera[0].failures[..],
            [Failure::ChunkRatio { ratio, bound }] if ratio < bound
        ));
        assert!(matches!(
            found.taxa[0].failures[..],
            [Failure::ChunkRatio { .. }]
        ));
    }

    #[test]
    fn a_fragment_counts_on_every_further_marker_gene_it_covers() {
        // Four reads on genes 1-4 of marker genome 60, each also covering two further genes: 12 genes hit.
        let interval = |gene: u32| PresenceInterval {
            reference_id: gene,
            span: (100, 200),
            mate_span: (0, 0),
            fixed_chunk: (100, 200),
            flex_chunk: (100, 200),
        };
        let placed = (0..4u32)
            .map(|i| Placed {
                read: i as usize,
                kind: Kind::NotKept,
                ties: 0,
                assembly: 60,
                identity: 0.99,
                bases: 100.0,
                read_len: 100,
                interval: interval(i + 1),
                gene_hits: 2 * i..2 * i + 2,
            })
            .collect::<Vec<_>>();
        let gene_hits = (5..=12).map(interval).collect::<Vec<_>>();
        let mut data = data(Vec::new());
        let model = marker_model();
        data.marker_model = Some(Arc::clone(&model));
        let tests = CoverageTests {
            multiple_tests: 1,
            chunk_cutoff: 0.0,
            ratio_bounds: (0.6, 0.75),
            markers: Some(model),
            min_marker_genes: 3,
        };
        let layout =
            GenomeLayout::new(std::iter::once((0, 1)).chain((1..=40).map(|_| (60, 1_000))));
        let taxonomy = taxonomy();
        let lineages = Lineages::new(&taxonomy);
        let level = Level {
            data: &data,
            tests: &tests,
            layout: &layout,
            below_species: (&[], policy(0.80)),
            not_kept: (&[], policy(0.95)),
            em: EmSettings {
                penalty: 3.0,
                max_iterations: 1000,
            },
            keep_shares: false,
            kept_beneath: &HashMap::new(),
        };
        let genomes = HashMap::from([(60, 40_000.0)]);
        let placing = Placing {
            level: &level,
            lineages: &lineages,
            placed: &placed,
            gene_hits: &gene_hits,
            sets: &[],
            genera: &[],
            genomes: &genomes,
            destination: Vec::new(),
            pending: HashMap::new(),
            unexplained: 0,
        };
        let measured = placing.measure(&[0, 1, 2, 3]).expect("12 genes hit");
        assert_eq!(measured.coverage.marker_genes_hit, Some(12));
        assert_eq!(measured.fragments, 4);
        // Without the further genes the four reads hit 4 genes, too few to be judged.
        let alone = placed
            .iter()
            .map(|fragment| Placed {
                gene_hits: 0..0,
                interval: fragment.interval,
                ..*fragment
            })
            .collect::<Vec<_>>();
        let placing = Placing {
            placed: &alone,
            ..placing
        };
        assert!(placing.measure(&[0, 1, 2, 3]).is_none());
    }

    /// Each row's taxid, rank and reads (to six decimals).
    fn rows(found: &HigherLevels) -> Vec<(u32, Rank, f64)> {
        found
            .rows
            .iter()
            .map(|row| (row.taxid, row.rank, (row.reads * 1e6).round() / 1e6))
            .collect()
    }

    #[test]
    fn a_tie_set_holds_the_taxa_within_the_margin_of_the_best() {
        let alignments = read(&[(1001, 1, 0.99), (1011, 11, 0.99), (1021, 21, 0.97)], 100);
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&alignments, policy(0.95));
        let mode = policy(0.95).mode;
        let (found, best) = ranked(candidates, mode).unwrap();
        // Species 1021 has two more errors: not tied; the best is the first tied species by taxid.
        assert_eq!(
            found.iter().map(|c| c.0).collect::<Vec<_>>(),
            [1001, 1011, 1021]
        );
        assert!((found[0].1 - 1.0).abs() < 1e-9 && found[0].1 == found[1].1);
        assert!((found[2].1 - 3.0).abs() < 1e-9);
        assert_eq!(tie_set(&found), [1001, 1011]);
        assert_eq!(best.taxid, 1001);
        assert!(ranked(&[], mode).is_none());
    }

    #[test]
    fn a_lineage_gives_the_taxa_at_each_rank_above_species() {
        let taxonomy = taxonomy();
        let lineages = Lineages::new(&taxonomy);
        assert_eq!(lineages.at_rank(1001, Rank::Genus), Some(100));
        assert_eq!(lineages.at_rank(100, Rank::Genus), Some(100));
        assert_eq!(lineages.at_rank(1001, Rank::Family), Some(40));
        assert_eq!(lineages.at_rank(5000, Rank::Genus), None);
        assert_eq!(lineages.ancestors(1001)[..3], [1001, 100, 40]);
    }

    #[test]
    fn a_tie_set_in_one_genus_is_unique_to_it() {
        let taxonomy = taxonomy();
        let lineages = Lineages::new(&taxonomy);
        let one = genera_of(&[1001, 1002], &lineages);
        assert_eq!((one.genera, one.unique), (vec![100], Some(100)));
        let two = genera_of(&[1001, 1011], &lineages);
        assert_eq!((two.genera, two.unique), (vec![100, 101], None));
        // A taxon without a genus keeps the set from lying in one genus.
        let partly = genera_of(&[1001, 5000], &lineages);
        assert_eq!((partly.genera, partly.unique), (vec![100], None));
    }

    /// Evidence measured on genome 7 of `model`, with `fragments` fragments there, `genes` marker genes hit, chunk
    /// ratio `coebr` and window spread `spread`.
    fn measured(
        model: Model,
        fragments: usize,
        genes: Option<usize>,
        coebr: f64,
        spread: Option<f64>,
    ) -> Measured {
        Measured {
            genome: 7,
            model,
            fragments,
            coverage: EvidenceCoverage {
                coebr,
                marker_genes_hit: genes,
                ..EvidenceCoverage::default()
            },
            spread,
            tied: Arc::from([7]),
            length: 1_000_000.0,
            whole_genome_scale: None,
        }
    }

    #[test]
    fn each_evidence_model_has_its_own_row_rules() {
        assert_eq!(Model::Genome.rules(), &GENOME_ROWS);
        assert_eq!(Model::Markers.rules(), &MARKER_ROWS);
        // Genomes: an amount in fragments, spread over windows above genus.
        assert_eq!((GENOME_ROWS.min_fragments, GENOME_ROWS.min_genes), (200, 0));
        assert_eq!(GENOME_ROWS.above_genus, Spread::Windows(0.5));
        // Marker genomes: an amount in genes hit, spread over genes above genus.
        assert_eq!((MARKER_ROWS.min_fragments, MARKER_ROWS.min_genes), (0, 10));
        assert_eq!(MARKER_ROWS.above_genus, Spread::ChunkRatio(0.4));
        assert_eq!(
            GENOME_ROWS.genus_min_chunk_ratio,
            MARKER_ROWS.genus_min_chunk_ratio
        );
    }

    #[test]
    fn genus_tests_take_their_bounds_from_the_evidence_model() {
        let marker = measured(Model::Markers, 40, Some(12), 0.35, None);
        assert_eq!(
            genus_failures(0.95, 0.9, &marker),
            [Failure::ChunkRatio {
                ratio: 0.35,
                bound: MARKER_ROWS.genus_min_chunk_ratio
            }]
        );
        let even = measured(Model::Markers, 40, Some(12), 0.9, None);
        assert_eq!(
            genus_failures(0.85, 0.9, &even),
            [Failure::Identity {
                identity: 0.85,
                floor: 0.9
            }]
        );
        assert!(genus_failures(0.95, 0.9, &even).is_empty());
        // A chunk ratio that cannot be measured fails.
        let unmeasured = measured(Model::Genome, 300, None, f64::NAN, Some(1.0));
        assert!(matches!(
            genus_failures(0.95, 0.9, &unmeasured)[..],
            [Failure::ChunkRatio { .. }]
        ));
    }

    #[test]
    fn rows_above_genus_spread_over_windows_on_genomes_and_over_genes_on_marker_genomes() {
        // A genome: the window spread decides, also when too few reads measure it.
        let genome = measured(Model::Genome, 300, None, 1.0, Some(0.4));
        assert_eq!(
            taxon_failures(&genome),
            [Failure::Windows {
                spread: Some(0.4),
                bound: 0.5
            }]
        );
        let unmeasured = measured(Model::Genome, 300, None, 1.0, None);
        assert_eq!(
            taxon_failures(&unmeasured)[0].to_string(),
            "window spread: too few reads to measure"
        );
        // A marker genome: the chunk ratio over its genes decides; windows play no part.
        let marker = measured(Model::Markers, 30, Some(15), 0.8, None);
        assert!(taxon_failures(&marker).is_empty());
        assert_eq!(
            taxon_evidence(&marker),
            "15 genes hit, chunk ratio 0.8000 >= 0.4"
        );
        let clumped = measured(Model::Markers, 30, Some(10), 0.3, Some(0.9));
        assert_eq!(
            taxon_failures(&clumped)[0].to_string(),
            "chunk ratio: 0.3000 < 0.4"
        );
    }

    #[test]
    fn a_genus_row_needs_enough_unique_reads_identity_and_an_even_chunk_ratio() {
        // 250 reads of a species not kept, unique to genus 100 and spread along its genome.
        let found = place_reads(Vec::new(), reads(&[(1001, 1, 0.99)], 250, true));
        assert_eq!(rows(&found), [(100, Rank::Genus, 250.0)]);
        let candidate = &found.genera[0];
        assert_eq!(
            (candidate.kind, candidate.unique, candidate.tied),
            (Kind::NotKept, 250, 0)
        );
        assert!(candidate.shown && candidate.failures.is_empty());
        assert!((candidate.identity - 0.99).abs() < 1e-6);
        assert!(
            (candidate.identity_floor - (0.95 - GENOME_ROWS.genus_identity_margin)).abs() < 1e-12
        );
        assert_eq!(candidate.measured.model, Model::Genome);
        // Depth: first-read bases over the genome's length.
        assert!((found.rows[0].depth - 250.0 * 100.0 / 1_000_000.0).abs() < 1e-12);
        // Fewer than 200 unique reads: no candidate; the reads go to the family, too few for a row there too.
        let found = place_reads(Vec::new(), reads(&[(1001, 1, 0.99)], 199, true));
        assert!(found.rows.is_empty() && found.genera.is_empty());
        assert_eq!(found.unexplained, 199);
        // Reads in one place: the chunk ratio fails, and the family's window spread too.
        let found = place_reads(Vec::new(), reads(&[(1001, 1, 0.99)], 250, false));
        assert!(found.rows.is_empty());
        assert!(
            matches!(found.genera[0].failures[..], [Failure::ChunkRatio { ratio, bound }] if ratio < bound && bound == 0.4)
        );
        assert_eq!(found.taxa[0].taxid, 40);
        assert_eq!(found.taxa[0].from_genera, 250);
        // Reads without a species-level match below the genus row's mean identity test (0.95 - 0.09): their family takes them.
        let found = place_reads(reads(&[(1001, 1, 0.85)], 250, true), Vec::new());
        assert!(matches!(
            found.genera[0].failures[..],
            [Failure::Identity { identity, .. }] if (identity - 0.85).abs() < 1e-6
        ));
        assert_eq!(rows(&found), [(40, Rank::Family, 250.0)]);
        assert_eq!(found.end_of_genus(100).map(|row| row.taxid), Some(40));
    }

    #[test]
    fn a_genus_row_holds_the_unique_reads_of_both_kinds() {
        // Genus 100: 250 even reads of species not kept, and 250 without a species match below the genus
        // identity. The first give the genus its row; the second fail on their own but join it.
        let found = place_reads(
            reads(&[(1001, 1, 0.85)], 250, true),
            reads(&[(1001, 1, 0.99)], 250, true),
        );
        assert_eq!(rows(&found), [(100, Rank::Genus, 500.0)]);
        assert_eq!(found.rows[0].unique, [250, 250]);
        let kinds = found
            .genera
            .iter()
            .map(|c| (c.kind, c.failures.is_empty(), c.shown))
            .collect::<Vec<_>>();
        assert_eq!(
            kinds,
            [
                (Kind::BelowSpecies, false, false),
                (Kind::NotKept, true, true)
            ]
        );
        // Unique reads of a kind too few to judge join the row as well.
        let found = place_reads(
            reads(&[(1001, 1, 0.85)], 50, true),
            reads(&[(1001, 1, 0.99)], 250, true),
        );
        assert_eq!(rows(&found), [(100, Rank::Genus, 300.0)]);
        assert_eq!(found.rows[0].unique, [50, 250]);
        assert_eq!(found.genera.len(), 1);
    }

    /// The shares of the placed fragments that more than one genus takes.
    fn split(found: &HigherLevels) -> Vec<&Vec<(u32, f64)>> {
        found
            .shares
            .iter()
            .map(|s| &s.1)
            .filter(|s| s.len() > 1)
            .collect()
    }

    #[test]
    fn reads_tied_across_genera_are_shared_by_the_em_among_the_genera_with_a_row() {
        let tied = reads(&[(1001, 1, 0.99), (1011, 11, 0.99)], 10, true);
        // Genus 100 with 300 unique reads, genus 101 with 200: the EM shares the tied reads 3:2.
        let mut not_kept = reads(&[(1001, 1, 0.99)], 300, true);
        not_kept.extend(reads(&[(1011, 11, 0.99)], 200, true));
        not_kept.extend(tied.clone());
        let found = place_reads(Vec::new(), not_kept);
        let reads_of = |found: &HigherLevels| {
            found
                .rows
                .iter()
                .map(|row| (row.taxid, row.reads))
                .collect::<Vec<_>>()
        };
        let rows_found = reads_of(&found);
        assert_eq!((rows_found[0].0, rows_found[1].0), (100, 101));
        assert!((rows_found[0].1 - 306.0).abs() < 1e-3 && (rows_found[1].1 - 204.0).abs() < 1e-3);
        let shared = split(&found);
        assert_eq!(shared.len(), 10);
        assert_eq!((shared[0][0].0, shared[0][1].0), (100, 101));
        assert!((shared[0][0].1 - 0.6).abs() < 1e-4);
        assert_eq!(found.genera[0].tied, 10);
        // Only genus 100 has a row: it takes the tied reads whole.
        let mut not_kept = reads(&[(1001, 1, 0.99)], 300, true);
        not_kept.extend(reads(&[(1011, 11, 0.99)], 150, true));
        not_kept.extend(tied.clone());
        let found = place_reads(Vec::new(), not_kept);
        assert_eq!(rows(&found)[0], (100, Rank::Genus, 310.0));
        // A genus with a row that matches worse than the tie set takes a share by its likelihood; a tie set
        // without a genus with a row goes above genus even when a worse match has one.
        let mut not_kept = reads(&[(1001, 1, 0.99)], 300, true);
        not_kept.extend(reads(&[(1021, 21, 0.99)], 300, true));
        not_kept.extend(reads(
            &[(1001, 1, 0.99), (1011, 11, 0.99), (1021, 21, 0.98)],
            10,
            true,
        ));
        not_kept.extend(reads(
            &[(1011, 11, 0.99), (1031, 31, 0.99), (1021, 21, 0.98)],
            5,
            true,
        ));
        let found = place_reads(Vec::new(), not_kept);
        let shared = split(&found);
        assert_eq!(shared.len(), 10);
        assert_eq!((shared[0][0].0, shared[0][1].0), (100, 102));
        assert!(shared[0][1].1 > 0.0 && shared[0][1].1 < 0.1);
        assert_eq!(found.unexplained, 5);
        // A genus row takes the tied reads of either kind: the tied reads without a species-level match go to
        // genus 100, whose row comes from reads of species not kept; genus 101's 150 unique reads are too few for
        // a row there or above.
        let mut below = reads(&[(1011, 11, 0.85)], 150, true);
        below.extend(reads(&[(1001, 1, 0.85), (1011, 11, 0.85)], 60, true));
        let found = place_reads(below, reads(&[(1001, 1, 0.99)], 300, true));
        assert_eq!(rows(&found), [(100, Rank::Genus, 360.0)]);
        assert_eq!(found.unexplained, 150);
        assert!(found.taxa.is_empty());
        // Tied reads without a genus row in their tie set go above genus, with the genera they list.
        let mut below = reads(&[(1011, 11, 0.85)], 150, true);
        below.extend(reads(&[(1001, 1, 0.85), (1011, 11, 0.85)], 60, true));
        let found = place_reads(below, Vec::new());
        assert_eq!(rows(&found), [(40, Rank::Family, 210.0)]);
        assert_eq!((found.taxa[0].tied, found.taxa[0].from_genera), (60, 150));
        assert_eq!(found.taxa[0].tied_genera, [100, 101]);
    }

    #[test]
    fn the_shares_of_the_placed_reads_sum_to_each_row() {
        let mut not_kept = reads(&[(1001, 1, 0.99)], 300, true);
        not_kept.extend(reads(&[(1011, 11, 0.99)], 200, true));
        not_kept.extend(reads(&[(1001, 1, 0.99), (1011, 11, 0.99)], 30, true));
        not_kept.extend(reads(&[(1021, 21, 0.99), (1031, 31, 0.99)], 250, true));
        let found = place_reads(reads(&[(1001, 1, 0.85)], 220, true), not_kept);
        assert_eq!(found.rows.len(), 3);
        for row in &found.rows {
            let sum = found
                .shares
                .iter()
                .flat_map(|s| &s.1)
                .filter(|share| share.0 == row.taxid)
                .map(|share| share.1)
                .sum::<f64>();
            assert!((sum - row.reads).abs() < 1e-9, "{}", row.taxid);
        }
    }

    #[test]
    fn other_reads_go_to_the_lowest_common_ancestor_of_their_matches() {
        let found_at = |on: &[(u32, u32, f32)]| {
            let found = place_reads(Vec::new(), reads(on, 250, true));
            (rows(&found), found.unexplained)
        };
        // Genera of one family, of one order, of one class.
        assert_eq!(
            found_at(&[(1001, 1, 0.99), (1011, 11, 0.99)]),
            (vec![(40, Rank::Family, 250.0)], 0)
        );
        assert_eq!(
            found_at(&[(1001, 1, 0.99), (1021, 21, 0.99)]),
            (vec![(30, Rank::Order, 250.0)], 0)
        );
        assert_eq!(
            found_at(&[(1001, 1, 0.99), (1031, 31, 0.99)]),
            (vec![(20, Rank::Class, 250.0)], 0)
        );
        // Taxa of two domains: none below the root.
        assert_eq!(
            found_at(&[(1001, 1, 0.99), (5000, 50, 0.99)]),
            (Vec::new(), 250)
        );
    }

    #[test]
    fn reads_of_a_taxon_without_a_row_move_up_rank_by_rank() {
        // Clustered reads tied across genera 100 and 101: family 40 fails its window spread, and so does every
        // rank above it; the reads stay unexplained above the domain.
        let found = place_reads(
            Vec::new(),
            reads(&[(1001, 1, 0.99), (1011, 11, 0.99)], 250, false),
        );
        assert!(found.rows.is_empty());
        let path = found
            .taxa
            .iter()
            .map(|taxon| (taxon.taxid, taxon.rank, taxon.fragments(), taxon.from_below))
            .collect::<Vec<_>>();
        assert_eq!(
            path,
            [
                (40, Rank::Family, 250, 0),
                (30, Rank::Order, 250, 250),
                (20, Rank::Class, 250, 250),
                (10, Rank::Phylum, 250, 250),
                (2, Rank::Domain, 250, 250),
            ]
        );
        assert!(
            matches!(found.taxa[0].failures[..], [Failure::Windows { spread: Some(spread), bound }] if spread < bound && bound == 0.5)
        );
        assert_eq!(found.unexplained, 250);
        assert!(found.end_of(40).is_none());
        // 150 tied reads at family 40 and 100 at family 41 are too few for either; together they make the
        // order's row, spread along the genome.
        let mut not_kept = reads(&[(1001, 1, 0.99), (1011, 11, 0.99)], 150, true);
        not_kept.extend(reads(&[(1021, 21, 0.99), (1022, 22, 0.99)], 100, true));
        let found = place_reads(Vec::new(), not_kept);
        // Species 1022 is not in the taxonomy: the second set's ancestor is the root, so its reads stay
        // unexplained; family 40's reads move up alone and are too few.
        assert!(found.rows.is_empty());
        assert_eq!(found.unexplained, 250);
        let mut not_kept = reads(&[(1001, 1, 0.99), (1011, 11, 0.99)], 150, true);
        not_kept.extend(reads(&[(1021, 21, 0.99), (1001, 1, 0.99)], 100, true));
        let found = place_reads(Vec::new(), not_kept);
        assert_eq!(rows(&found), [(30, Rank::Order, 250.0)]);
        assert_eq!(found.end_of(40).map(|row| row.taxid), Some(30));
    }

    #[test]
    fn a_row_s_depth_is_its_reads_over_the_genome_they_fit_best() {
        let depth = |found: &HigherLevels| found.rows[0].depth;
        // Genus 100: 150 reads on genome 1 (1 Mb) and 100 on genome 5 (100 kb). All 250 count over genome 1,
        // which holds most of them.
        let mut not_kept = reads(&[(1001, 1, 0.99)], 150, true);
        not_kept.extend((0..100).map(|i| read_small(&[(1002, 5, 0.99)], i * 990)));
        let found = place_reads(Vec::new(), not_kept);
        assert_eq!(rows(&found), [(100, Rank::Genus, 250.0)]);
        assert_eq!(found.rows[0].measured.genome, 1);
        assert!((depth(&found) - 250.0 * 100.0 / 1e6).abs() < 1e-12);
        // Reads tied across genera 100 and 101, 125 best on genome 5 (100 kb) and 125 on genome 6 (1 Mb): family
        // 40 is measured on genome 5 and its depth is the mean of the depths over the two genomes.
        let mut across = (0..125)
            .map(|i| read_small(&[(1001, 5, 0.99), (1011, 11, 0.99)], i * 790))
            .collect::<Vec<_>>();
        across.extend(reads(&[(1001, 6, 0.99), (1011, 11, 0.99)], 125, true));
        let found = place_reads(Vec::new(), across);
        assert_eq!(rows(&found), [(40, Rank::Family, 250.0)]);
        assert_eq!(found.rows[0].measured.genome, 5);
        assert!((depth(&found) - 250.0 * 100.0 * (1.0 / 1e5 + 1.0 / 1e6) / 2.0).abs() < 1e-12);
    }
}

//! Rows above species: the fragments that no kept species takes, placed at the
//! lowest rank that their matches and their evidence support.
//!
//! Two kinds of fragments come here:
//!
//! - fragments without a species-level candidate, placed at the floor gates
//!   (`--family-identity`, a fixed mapped length on short reads or the species
//!   length cutoff on long reads, and a fixed aligned fraction); on noisy long
//!   reads not those whose best match, or one tied with it, is a kept
//!   species, which count for that species;
//! - fragments with species-level candidates of which none is kept, placed at
//!   the species gates.
//!
//! A fragment's tie set holds the taxa (species, or units labelled above
//! species) whose placements lie within the tie margin of its best one.
//! Placements on marker genomes take no part: a marker gene holds too little
//! of a genome for its reads to stand for a taxon above species.
//!
//! 1. Genus rows. A fragment whose tie set lies in one genus is unique to it.
//!    A genus gets a row when, for one kind of fragment or both, at least
//!    [`MIN_FRAGMENTS`] fragments of that kind are unique to it, their mean
//!    identity is at least the run's species identity cutoff minus
//!    [`GENUS_IDENTITY_BELOW_SPECIES`], and their chunk ratio (COEBR) on the
//!    genome holding most of them is at least [`GENUS_MIN_CHUNK_RATIO`]. Its
//!    row then holds the fragments of both kinds unique to it; the unique
//!    fragments of a genus without a row go to its family (step 3).
//! 2. Ties between genera. A fragment whose tie set spans several genera, one
//!    or more of them with a genus row, goes to the genera with a genus row,
//!    shared by the likelihood EM over those genera, as the species level's EM
//!    shares fragments among kept species; to step 3 when none of its tie
//!    set's genera has one.
//! 3. Rows above genus. The other fragments go to the lowest common ancestor of
//!    their tie set's taxa, or the first taxon above it at family, order,
//!    class, phylum or domain rank; a fragment with none below the root stays
//!    unexplained. A taxon gets a row when at least [`MIN_FRAGMENTS`]
//!    fragments reach it and they pass [`taxon_failures`]; otherwise they move
//!    to the next taxon above it at one of those ranks, and stay unexplained
//!    above domain.
//!
//! A row's depth is the first-read bases of all its fragments over the indexed
//! length of the genome they fit best as a set: the genome holding most of
//! their best placements, on which the row's tests measure them. Genomes tied
//! exactly with it on that number share the depth: it is the mean of the
//! depths over each of them.
//!
//! On a database of sampled genomes, reads from the unindexed parts of a kept
//! genome land on the sampled segments of its relatives, while the kept
//! species already stands for its whole genome: no taxon that holds a kept
//! species gets a row there.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::sync::Arc;

use super::{
    best_match_kept, candidate_errors, fold_candidates, longest_alignment, share_reads, tie_margin,
    Fragment, SpeciesTests, EM_TOLERANCE, LONG_TIE_SHARE,
};
use crate::profiler::community::MappingData;
use crate::profiler::evidence::{
    CandidatePlacement, EvidenceCoverage, EvidencePolicy, PresenceAssembly, PresenceInterval,
    TerminalCandidate,
};
use crate::profiler::likelihood_em;
use crate::profiler::rollup::GenomeLayout;
use crate::taxonomy::Taxonomy;

/// Fewest fragments of one kind unique to a genus for a genus row, and fewest
/// fragments reaching a taxon above genus for its row.
pub(super) const MIN_FRAGMENTS: usize = 200;
/// How far the mean identity of a genus's unique fragments may lie below the
/// run's species identity cutoff.
pub(super) const GENUS_IDENTITY_BELOW_SPECIES: f64 = 0.09;
/// Lowest chunk ratio (COEBR) of a genus's unique fragments on the genome
/// holding most of them.
pub(super) const GENUS_MIN_CHUNK_RATIO: f64 = 0.40;
/// Lowest share of the windows of the genome holding most of a taxon's
/// fragments that they reach, for a row above genus.
pub(super) const MIN_WINDOW_SPREAD: f64 = 0.50;

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
type Lineage = Arc<[(u32, Option<Rank>)]>;

/// Lineages read from the taxonomy, once per taxon.
struct Lineages<'t> {
    taxonomy: &'t Taxonomy,
    known: HashMap<u32, Lineage>,
}

impl<'t> Lineages<'t> {
    fn new(taxonomy: &'t Taxonomy) -> Self {
        Self {
            taxonomy,
            known: HashMap::new(),
        }
    }

    fn of(&mut self, taxid: u32) -> Lineage {
        let taxonomy = self.taxonomy;
        Arc::clone(self.known.entry(taxid).or_insert_with(|| {
            taxonomy
                .get_parents(taxid)
                .into_iter()
                .map(|(ancestor, _, rank)| (ancestor, Rank::of(&rank)))
                .collect()
        }))
    }

    /// The genus of `taxid`: the taxon itself at genus rank, or its ancestor
    /// there.
    fn genus(&mut self, taxid: u32) -> Option<u32> {
        self.of(taxid)
            .iter()
            .find(|node| node.1 == Some(Rank::Genus))
            .map(|node| node.0)
    }

    /// The first taxon of the lineage of `taxid` ranked above `rank`.
    fn above(&mut self, taxid: u32, rank: Rank) -> Option<(u32, Rank)> {
        first_above(&self.of(taxid), rank)
    }

    /// Where the fragments of a tie set of `taxa` that lies in no one genus
    /// start: their lowest common ancestor when it is ranked above genus,
    /// else the first taxon above it that is; `None` when there is none below
    /// the root.
    fn common_start(&mut self, taxa: &[u32]) -> Option<(u32, Rank)> {
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
    pub(super) taxonomy: &'a Taxonomy,
    pub(super) tests: &'a SpeciesTests,
    pub(super) layout: &'a GenomeLayout,
    /// The reads without a species-level candidate, and the floor gates.
    pub(super) below_species: (&'a [usize], EvidencePolicy),
    /// The reads of fragments whose candidate species are none kept, and the
    /// species gates.
    pub(super) not_kept: (&'a [usize], EvidencePolicy),
    /// On noisy long reads, the kept species whose best-matched reads below
    /// the species gates count for them; otherwise empty.
    pub(super) noisy_kept: &'a HashSet<u32>,
    /// On a database of sampled genomes, every taxon that holds a kept
    /// species; otherwise empty.
    pub(super) holds_kept: &'a HashSet<u32>,
    /// The EM's error penalty and its most iterations.
    pub(super) penalty: f64,
    pub(super) max_iterations: usize,
    /// Keep each fragment's shares for the per-read classifications.
    pub(super) keep_shares: bool,
}

impl Level<'_> {
    /// The gates of the fragments of `kind`.
    fn policy(&self, kind: Kind) -> EvidencePolicy {
        match kind {
            Kind::BelowSpecies => self.below_species.1,
            Kind::NotKept => self.not_kept.1,
        }
    }

    /// Whether a candidate takes part: not on a marker genome.
    fn takes_part(&self, candidate: &TerminalCandidate) -> bool {
        self.data
            .marker_model
            .as_deref()
            .is_none_or(|model| !model.is_marker_genome(candidate.best.assembly_id))
    }
}

/// A fragment of the level: its read, its kind, its tie set (an index into
/// the level's tie sets), and its best placement's genome, identity, first
/// read's mapped bases and interval.
struct Placed {
    read: usize,
    kind: Kind,
    ties: u32,
    assembly: u32,
    identity: f32,
    bases: f32,
    interval: PresenceInterval,
}

/// The level's fragments, their distinct tie sets (taxa by taxid), and the
/// length and sampling scale of each genome holding a best placement.
#[derive(Default)]
struct Collected {
    placed: Vec<Placed>,
    sets: Vec<Vec<u32>>,
    genomes: HashMap<u32, (f64, f64)>,
}

/// The genera of a tie set's taxa that have one, by taxid, and the genus
/// when every taxon lies in it.
struct Genera {
    genera: Vec<u32>,
    unique: Option<u32>,
}

/// A candidate row's evidence: the genome holding most of its fragments
/// (fewest assembly id on equal numbers), their number there, and their
/// coverage and window spread on it; and the length its depth is taken over:
/// that genome's indexed length, or, when other genomes hold as many of its
/// fragments, the harmonic mean of their indexed lengths, so that the depth is
/// the mean of the depths on them.
#[derive(Clone, Copy, Debug)]
pub(super) struct Measured {
    pub(super) genome: u32,
    pub(super) fragments: usize,
    pub(super) coverage: EvidenceCoverage,
    pub(super) spread: Option<f64>,
    pub(super) length: f64,
}

/// A test a candidate row fails.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Failure {
    Identity { identity: f64, floor: f64 },
    ChunkRatio(f64),
    Spread(Option<f64>),
    HoldsKept,
}

impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Self::Identity { identity, floor } => {
                write!(f, "identity: {identity:.4} < {floor:.4}")
            }
            Self::ChunkRatio(ratio) => {
                write!(f, "chunk ratio: {ratio:.4} < {GENUS_MIN_CHUNK_RATIO}")
            }
            Self::Spread(None) => write!(f, "window spread: too few reads to measure"),
            Self::Spread(Some(spread)) => {
                write!(f, "window spread: {spread:.4} < {MIN_WINDOW_SPREAD}")
            }
            Self::HoldsKept => write!(f, "holds a kept species, on a database of sampled genomes"),
        }
    }
}

/// The tests a genus fails with its unique fragments of one kind: their mean
/// `identity` against `floor`, their chunk ratio in `coverage`, and whether
/// the genus `holds_kept` species on a database of sampled genomes.
fn genus_failures(
    identity: f64,
    floor: f64,
    coverage: &EvidenceCoverage,
    holds_kept: bool,
) -> Vec<Failure> {
    let mut failed = Vec::new();
    if holds_kept {
        failed.push(Failure::HoldsKept);
    }
    if identity < floor {
        failed.push(Failure::Identity { identity, floor });
    }
    if coverage.coebr.is_nan() || coverage.coebr < GENUS_MIN_CHUNK_RATIO {
        failed.push(Failure::ChunkRatio(coverage.coebr));
    }
    failed
}

/// The tests a taxon above genus fails with the fragments that reach it: the
/// share of the windows of the genome holding most of them that they reach,
/// and whether the taxon `holds_kept` species on a database of sampled
/// genomes.
pub(super) fn taxon_failures(measured: &Measured, holds_kept: bool) -> Vec<Failure> {
    let mut failed = Vec::new();
    if holds_kept {
        failed.push(Failure::HoldsKept);
    }
    if measured
        .spread
        .is_none_or(|spread| spread < MIN_WINDOW_SPREAD)
    {
        failed.push(Failure::Spread(measured.spread));
    }
    failed
}

/// The evidence of a taxon above genus that passes [`taxon_failures`], in
/// words.
pub(super) fn taxon_evidence(measured: &Measured) -> String {
    format!(
        "window spread {:.4} >= {MIN_WINDOW_SPREAD}",
        measured.spread.unwrap_or(f64::NAN)
    )
}

/// A genus with at least [`MIN_FRAGMENTS`] fragments of one kind unique to
/// it: those fragments, the fragments of the kind tied between it and other
/// genera, their mean identity and their evidence, the tests it fails, and
/// whether its evidence is the one its genus row shows (its passing kind with
/// most unique fragments).
pub(super) struct GenusCandidate {
    pub(super) genus: u32,
    pub(super) kind: Kind,
    pub(super) unique: usize,
    pub(super) tied: usize,
    pub(super) identity: f64,
    pub(super) measured: Measured,
    pub(super) failures: Vec<Failure>,
    pub(super) shown: bool,
}

/// A taxon above genus that at least [`MIN_FRAGMENTS`] fragments reach: those
/// tied across genera, those unique to genera without a genus row, and those
/// moved up from lower ranks; the genera its tied fragments list most often
/// (at most three, most first); their evidence and the tests it fails.
pub(super) struct TaxonCandidate {
    pub(super) taxid: u32,
    pub(super) rank: Rank,
    pub(super) tied: usize,
    pub(super) from_genera: usize,
    pub(super) from_below: usize,
    pub(super) tied_genera: Vec<u32>,
    pub(super) measured: Measured,
    pub(super) failures: Vec<Failure>,
}

impl TaxonCandidate {
    pub(super) fn fragments(&self) -> usize {
        self.tied + self.from_genera + self.from_below
    }
}

/// A row above species: its fragments (a pair counting once), of a genus row
/// those unique to the genus by kind, and their full-index equivalents, its
/// depth (see the [module documentation](self)), and its evidence.
pub(super) struct HigherRow {
    pub(super) taxid: u32,
    pub(super) rank: Rank,
    pub(super) reads: f64,
    pub(super) unique: [usize; 2],
    pub(super) equivalent: f64,
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
    /// The lowest mean identity of a genus's unique fragments.
    pub(super) identity_floor: f64,
    /// Fragments of each kind, and those left unexplained.
    pub(super) fragments: [usize; 2],
    pub(super) unexplained: usize,
    /// With `keep_shares`, each placed fragment's read and its shares.
    pub(super) shares: Vec<(usize, Vec<(u32, f64)>)>,
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

/// Place the level's fragments ([module documentation](self)).
pub(super) fn place(level: &Level) -> HigherLevels {
    let mut lineages = Lineages::new(level.taxonomy);
    let Collected {
        placed,
        sets,
        genomes,
    } = collect(level);
    let genera = sets
        .iter()
        .map(|taxa| genera_of(taxa, &mut lineages))
        .collect::<Vec<_>>();
    let identity_floor = level.not_kept.1.identity - GENUS_IDENTITY_BELOW_SPECIES;
    let measure = |list: &[u32]| measure(&placed, list, &genomes, level);
    let mut destination = vec![Destination::Unexplained; placed.len()];

    // Step 1: genus rows from unique fragments, by kind.
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
    let mut candidates = Vec::new();
    for (&(genus, kind), list) in &unique {
        if list.len() < MIN_FRAGMENTS {
            continue;
        }
        let identity = list
            .iter()
            .map(|&i| f64::from(placed[i as usize].identity))
            .sum::<f64>()
            / list.len() as f64;
        let measured = measure(list);
        candidates.push(GenusCandidate {
            genus,
            kind,
            unique: list.len(),
            tied: tied.get(&(genus, kind)).copied().unwrap_or(0),
            identity,
            measured,
            failures: genus_failures(
                identity,
                identity_floor,
                &measured.coverage,
                level.holds_kept.contains(&genus),
            ),
            shown: false,
        });
    }
    // A genus with a row from either kind holds the fragments of both kinds
    // unique to it; those of the others go to their family.
    let mut with_row: BTreeMap<u32, [usize; 2]> = candidates
        .iter()
        .filter(|candidate| candidate.failures.is_empty())
        .map(|candidate| (candidate.genus, [0; 2]))
        .collect();
    let mut pending: HashMap<u32, (Rank, Arrivals)> = HashMap::new();
    let mut genus_start: HashMap<u32, Option<u32>> = HashMap::new();
    let mut unexplained = 0;
    for ((genus, kind), list) in unique {
        if let Some(unique) = with_row.get_mut(&genus) {
            unique[kind as usize] = list.len();
            for &i in &list {
                destination[i as usize] = Destination::Taxon(genus);
            }
            continue;
        }
        let start = lineages.above(genus, Rank::Genus);
        genus_start.insert(genus, start.map(|(taxid, _)| taxid));
        match start {
            Some((start, rank)) => {
                let arrivals = &mut pending
                    .entry(start)
                    .or_insert_with(|| (rank, Arrivals::default()))
                    .1;
                arrivals.from_genera += list.len();
                arrivals.fragments.extend(list);
            }
            None => unexplained += list.len(),
        }
    }
    // The genus row shows the evidence of its kind with most unique fragments.
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

    // Step 2: fragments tied across genera go to the genera with a genus row
    // when their tie set holds one, shared by the EM.
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
                let arrivals = &mut pending
                    .entry(start)
                    .or_insert_with(|| (rank, Arrivals::default()))
                    .1;
                arrivals.tied += 1;
                *arrivals.sets.entry(fragment.ties).or_default() += 1;
                arrivals.fragments.push(i as u32);
            }
            None => unexplained += 1,
        }
    }
    let mut tied_shares: Vec<Vec<(u32, f64)>> = Vec::new();
    if !tied.is_empty() {
        let (em, groups) = share_tied(level, &placed, &tied, &with_row, &mut lineages);
        tied_shares = em.posteriors;
        for (i, group) in groups {
            let group = group.expect("a fragment with a genus row among its tie set's genera");
            destination[i as usize] = Destination::Shared(group as u32);
        }
    }

    // Step 3: rows above genus, rank by rank.
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
        let mut here = pending
            .iter()
            .filter(|(_, (of, _))| *of == rank)
            .map(|(&taxid, _)| taxid)
            .collect::<Vec<_>>();
        here.sort_unstable();
        for taxid in here {
            let (_, arrivals) = pending.remove(&taxid).expect("listed above");
            let passed = arrivals.fragments.len() >= MIN_FRAGMENTS && {
                let measured = measure(&arrivals.fragments);
                let failures = taxon_failures(&measured, level.holds_kept.contains(&taxid));
                let passed = failures.is_empty();
                if passed {
                    evidence.insert((rank, taxid), measured);
                }
                taxa.push(TaxonCandidate {
                    taxid,
                    rank,
                    tied: arrivals.tied,
                    from_genera: arrivals.from_genera,
                    from_below: arrivals.from_below,
                    tied_genera: most_listed(&arrivals.sets, &genera),
                    measured,
                    failures,
                });
                passed
            };
            if passed {
                for &i in &arrivals.fragments {
                    destination[i as usize] = Destination::Taxon(taxid);
                }
                continue;
            }
            let next = lineages.above(taxid, rank);
            moved.insert(taxid, next.map(|(next, _)| next));
            match next {
                Some((next, next_rank)) => {
                    let above = &mut pending
                        .entry(next)
                        .or_insert_with(|| (next_rank, Arrivals::default()))
                        .1;
                    above.from_below += arrivals.fragments.len();
                    above.fragments.extend(arrivals.fragments);
                }
                None => unexplained += arrivals.fragments.len(),
            }
        }
    }
    debug_assert!(pending.is_empty());

    // The rows: each fragment counted where it went.
    for candidate in candidates.iter().filter(|candidate| candidate.shown) {
        evidence.insert((Rank::Genus, candidate.genus), candidate.measured);
    }
    let mut totals: HashMap<u32, [f64; 3]> = HashMap::new();
    let mut shares = Vec::new();
    for (fragment, to) in placed.iter().zip(&destination) {
        let whole;
        let parts: &[(u32, f64)] = match *to {
            Destination::Unexplained => continue,
            Destination::Taxon(taxid) => {
                whole = [(taxid, 1.0)];
                &whole
            }
            Destination::Shared(index) => &tied_shares[index as usize],
        };
        let scale = genomes[&fragment.assembly].1;
        for &(taxid, share) in parts {
            let total = totals.entry(taxid).or_default();
            total[0] += share;
            total[1] += share * scale;
            total[2] += share * f64::from(fragment.bases);
        }
        if level.keep_shares {
            shares.push((fragment.read, parts.to_vec()));
        }
    }
    let rows = evidence
        .into_iter()
        .map(|((rank, taxid), measured)| {
            let [reads, equivalent, bases] = totals[&taxid];
            HigherRow {
                taxid,
                rank,
                reads,
                unique: with_row.get(&taxid).copied().unwrap_or_default(),
                equivalent,
                depth: bases / measured.length,
                measured,
            }
        })
        .collect();
    let mut fragments = [0; 2];
    for fragment in &placed {
        fragments[fragment.kind as usize] += 1;
    }
    HigherLevels {
        rows,
        genera: candidates,
        taxa,
        identity_floor,
        fragments,
        unexplained,
        shares,
        moved,
        genus_start,
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
    lineages: &mut Lineages,
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
            |_, candidates| {
                ranked(candidates, policy, |c| level.takes_part(c)).map(|found| found.0)
            },
            |&i, found| {
                let fragment = &placed[i as usize];
                let candidates = found
                    .iter()
                    .flatten()
                    .filter_map(|&(taxid, errors)| {
                        lineages
                            .genus(taxid)
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
    let shared = share_reads(&fragments, |_| true, level.penalty);
    let em = likelihood_em::run(
        &shared.groups,
        level.penalty,
        level.max_iterations,
        EM_TOLERANCE,
    );
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
        let leave_kept = kind == Kind::BelowSpecies && !level.noisy_kept.is_empty();
        fold_candidates(
            level.data,
            policy,
            reads,
            |&read| read,
            |&read, candidates| {
                if leave_kept && best_match_kept(read, candidates, policy, level.noisy_kept) {
                    return None;
                }
                ranked(candidates, policy, |c| level.takes_part(c))
                    .map(|(ranked, best)| (tie_set(&ranked), best))
            },
            |&read, result| {
                let Some((taxa, best)) = result else {
                    return;
                };
                let next = index.len() as u32;
                let ties = *index.entry(taxa).or_insert(next);
                found
                    .genomes
                    .entry(best.assembly_id)
                    .or_insert((best.asm_len, best.sampling_scale));
                found.placed.push(Placed {
                    read,
                    kind,
                    ties,
                    assembly: best.assembly_id,
                    identity: (1.0 - best.pair_errors / f64::from(best.both_mapped.max(1))) as f32,
                    bases: best.r1_mapped as f32,
                    interval: PresenceInterval {
                        reference_id: best.reference_id,
                        span: best.span,
                        mate_span: best.mate_span,
                        fixed_chunk: best.fixed_chunk,
                        flex_chunk: best.flex_chunk,
                    },
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

/// A fragment's `candidates` under `policy` among those that take part,
/// ranked as the species level ranks them ([`super::fragment_of`]): each
/// taxon with its estimated errors, those within the tie margin of the best
/// carrying the best's errors, by errors, then taxid; and the best
/// candidate's placement. `None` without a candidate.
fn ranked(
    candidates: &[TerminalCandidate],
    policy: EvidencePolicy,
    takes_part: impl Fn(&TerminalCandidate) -> bool,
) -> Option<(Vec<(u32, f64)>, CandidatePlacement)> {
    let longest = longest_alignment(candidates);
    let mut ranked = candidates
        .iter()
        .filter(|candidate| takes_part(candidate))
        .map(|candidate| (candidate_errors(&candidate.best, longest), candidate))
        .collect::<Vec<_>>();
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.taxid.cmp(&b.1.taxid)));
    let &(fewest, best) = ranked.first()?;
    let margin = tie_margin(
        policy.mode,
        f64::from(best.best.both_mapped),
        LONG_TIE_SHARE,
    );
    let taxa = ranked
        .iter()
        .map(|&(errors, candidate)| {
            let tied = errors - fewest < margin;
            (candidate.taxid, if tied { fewest } else { errors })
        })
        .collect();
    Some((taxa, best.best))
}

/// The taxa of a fragment's tie set: those of its `ranked` candidates with
/// the best's errors, by taxid.
fn tie_set(ranked: &[(u32, f64)]) -> Vec<u32> {
    let fewest = ranked[0].1;
    let mut taxa = ranked
        .iter()
        .take_while(|candidate| candidate.1 == fewest)
        .map(|candidate| candidate.0)
        .collect::<Vec<_>>();
    taxa.sort_unstable();
    taxa.dedup();
    taxa
}

/// The genera of a tie set's `taxa`.
fn genera_of(taxa: &[u32], lineages: &mut Lineages) -> Genera {
    let of = taxa
        .iter()
        .map(|&taxid| lineages.genus(taxid))
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

/// The evidence of the fragments `list` ([`Measured`]).
fn measure(
    placed: &[Placed],
    list: &[u32],
    genomes: &HashMap<u32, (f64, f64)>,
    level: &Level,
) -> Measured {
    let mut counts: HashMap<u32, usize> = HashMap::new();
    for &i in list {
        *counts.entry(placed[i as usize].assembly).or_default() += 1;
    }
    let (&genome, &fragments) = counts
        .iter()
        .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
        .expect("a candidate row has fragments");
    // The indexed lengths of the genomes that hold as many of the fragments;
    // several are summed in sorted order, which the map's order cannot change.
    let mut tied = counts
        .iter()
        .filter(|&(_, &n)| n == fragments)
        .map(|(&id, _)| genomes[&id].0.max(1.0))
        .collect::<Vec<_>>();
    let length = match tied[..] {
        [length] => length,
        _ => {
            tied.sort_unstable_by(f64::total_cmp);
            tied.len() as f64 / tied.iter().map(|length| 1.0 / length).sum::<f64>()
        }
    };
    let (asm_len, sampling_scale) = genomes[&genome];
    let mut assembly = PresenceAssembly {
        asm_len,
        sampling_scale,
        count: fragments as u64,
        intervals: Vec::with_capacity(fragments),
        ..PresenceAssembly::default()
    };
    for &i in list {
        let fragment = &placed[i as usize];
        if fragment.assembly == genome {
            let span = fragment.interval.span;
            assembly.reference_span_bases += u64::from(span.1.saturating_sub(span.0));
            assembly.intervals.push(fragment.interval);
        }
    }
    let spread = level.layout.window_spread(genome, &assembly.intervals);
    let coverage = level.tests.measure(genome, &mut assembly);
    Measured {
        genome,
        fragments,
        coverage,
        spread,
        length,
    }
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
    use super::super::test_fixtures::{data, row};
    use super::*;
    use crate::profiler::evidence::{CandidateScratch, EvidenceAlignment, EvidenceMode};

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

    /// A read as [`read`], with genome 5 of 100 kb indexed, sampled at `scale`.
    fn read_small(on: &[(u32, u32, f32)], start: u32, scale: f32) -> Vec<EvidenceAlignment> {
        read(on, start)
            .into_iter()
            .map(|a| match a.assembly_id {
                5 => EvidenceAlignment {
                    asm_len: 100_000.0,
                    sampling_scale: scale,
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
    /// (reads of species not kept, at the species identity 0.95); with `holds_kept` on a sampled database.
    fn place_reads(
        below: Vec<Vec<EvidenceAlignment>>,
        not_kept: Vec<Vec<EvidenceAlignment>>,
        holds_kept: &HashSet<u32>,
    ) -> HigherLevels {
        let below_reads = (0..below.len()).collect::<Vec<_>>();
        let not_kept_reads = (below.len()..below.len() + not_kept.len()).collect::<Vec<_>>();
        let data = data(below.into_iter().chain(not_kept).collect());
        let tests = SpeciesTests {
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
        place(&Level {
            data: &data,
            taxonomy: &taxonomy,
            tests: &tests,
            layout: &layout,
            below_species: (&below_reads, policy(0.80)),
            not_kept: (&not_kept_reads, policy(0.95)),
            noisy_kept: &HashSet::new(),
            holds_kept,
            penalty: 3.0,
            max_iterations: 1000,
            keep_shares: true,
        })
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
        let (found, best) = ranked(candidates, policy(0.95), |_| true).unwrap();
        // Species 1021 has two more errors: not tied; the best is the first tied species by taxid.
        assert_eq!(
            found.iter().map(|c| c.0).collect::<Vec<_>>(),
            [1001, 1011, 1021]
        );
        assert!((found[0].1 - 1.0).abs() < 1e-9 && found[0].1 == found[1].1);
        assert!((found[2].1 - 3.0).abs() < 1e-9);
        assert_eq!(tie_set(&found), [1001, 1011]);
        assert_eq!(best.taxid, 1001);
        // Candidates that take no part are left out before the ranking.
        let (found, best) = ranked(candidates, policy(0.95), |c| c.taxid != 1001).unwrap();
        assert_eq!((tie_set(&found), best.taxid), (vec![1011], 1011));
        assert!(ranked(candidates, policy(0.95), |_| false).is_none());
    }

    #[test]
    fn a_tie_set_in_one_genus_is_unique_to_it() {
        let taxonomy = taxonomy();
        let mut lineages = Lineages::new(&taxonomy);
        let one = genera_of(&[1001, 1002], &mut lineages);
        assert_eq!((one.genera, one.unique), (vec![100], Some(100)));
        let two = genera_of(&[1001, 1011], &mut lineages);
        assert_eq!((two.genera, two.unique), (vec![100, 101], None));
        // A taxon without a genus keeps the set from lying in one genus.
        let partly = genera_of(&[1001, 5000], &mut lineages);
        assert_eq!((partly.genera, partly.unique), (vec![100], None));
    }

    #[test]
    fn a_genus_row_needs_enough_unique_reads_identity_and_an_even_chunk_ratio() {
        let none = HashSet::new();
        // 250 reads of a species not kept, unique to genus 100 and spread along its genome.
        let found = place_reads(Vec::new(), reads(&[(1001, 1, 0.99)], 250, true), &none);
        assert_eq!(rows(&found), [(100, Rank::Genus, 250.0)]);
        let candidate = &found.genera[0];
        assert_eq!(
            (candidate.kind, candidate.unique, candidate.tied),
            (Kind::NotKept, 250, 0)
        );
        assert!(candidate.shown && candidate.failures.is_empty());
        assert!((candidate.identity - 0.99).abs() < 1e-6);
        assert!((found.identity_floor - (0.95 - GENUS_IDENTITY_BELOW_SPECIES)).abs() < 1e-12);
        // Depth: first-read bases over the genome's length.
        assert!((found.rows[0].depth - 250.0 * 100.0 / 1_000_000.0).abs() < 1e-12);
        // Fewer than 200 unique reads: no candidate; the reads go to the family, too few for a row there too.
        let found = place_reads(Vec::new(), reads(&[(1001, 1, 0.99)], 199, true), &none);
        assert!(found.rows.is_empty() && found.genera.is_empty());
        assert_eq!(found.unexplained, 199);
        // Reads in one place: the chunk ratio fails, and the family's window spread too.
        let found = place_reads(Vec::new(), reads(&[(1001, 1, 0.99)], 250, false), &none);
        assert!(found.rows.is_empty());
        assert!(
            matches!(found.genera[0].failures[..], [Failure::ChunkRatio(ratio)] if ratio < GENUS_MIN_CHUNK_RATIO)
        );
        assert_eq!(found.taxa[0].taxid, 40);
        assert_eq!(found.taxa[0].from_genera, 250);
        // Reads without a species-level match below the genus identity (0.95 - 0.09): their family takes them.
        let found = place_reads(reads(&[(1001, 1, 0.85)], 250, true), Vec::new(), &none);
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
            &HashSet::new(),
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
            &HashSet::new(),
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
        let none = HashSet::new();
        let tied = reads(&[(1001, 1, 0.99), (1011, 11, 0.99)], 10, true);
        // Genus 100 with 300 unique reads, genus 101 with 200: the EM shares the tied reads 3:2.
        let mut not_kept = reads(&[(1001, 1, 0.99)], 300, true);
        not_kept.extend(reads(&[(1011, 11, 0.99)], 200, true));
        not_kept.extend(tied.clone());
        let found = place_reads(Vec::new(), not_kept, &none);
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
        let found = place_reads(Vec::new(), not_kept, &none);
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
        let found = place_reads(Vec::new(), not_kept, &none);
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
        let found = place_reads(below, reads(&[(1001, 1, 0.99)], 300, true), &none);
        assert_eq!(rows(&found), [(100, Rank::Genus, 360.0)]);
        assert_eq!(found.unexplained, 150);
        assert!(found.taxa.is_empty());
        // Tied reads without a genus row in their tie set go above genus, with the genera they list.
        let mut below = reads(&[(1011, 11, 0.85)], 150, true);
        below.extend(reads(&[(1001, 1, 0.85), (1011, 11, 0.85)], 60, true));
        let found = place_reads(below, Vec::new(), &none);
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
        let found = place_reads(
            reads(&[(1001, 1, 0.85)], 220, true),
            not_kept,
            &HashSet::new(),
        );
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
        let none = HashSet::new();
        let found_at = |on: &[(u32, u32, f32)]| {
            let found = place_reads(Vec::new(), reads(on, 250, true), &none);
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
        let none = HashSet::new();
        // Clustered reads tied across genera 100 and 101: family 40 fails its window spread, and so does every
        // rank above it; the reads stay unexplained above the domain.
        let found = place_reads(
            Vec::new(),
            reads(&[(1001, 1, 0.99), (1011, 11, 0.99)], 250, false),
            &none,
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
            matches!(found.taxa[0].failures[..], [Failure::Spread(Some(spread))] if spread < 0.5)
        );
        assert_eq!(found.unexplained, 250);
        assert!(found.end_of(40).is_none());
        // 150 tied reads at family 40 and 100 at family 41 are too few for either; together they make the
        // order's row, spread along the genome.
        let mut not_kept = reads(&[(1001, 1, 0.99), (1011, 11, 0.99)], 150, true);
        not_kept.extend(reads(&[(1021, 21, 0.99), (1022, 22, 0.99)], 100, true));
        let found = place_reads(Vec::new(), not_kept, &none);
        // Species 1022 is not in the taxonomy: the second set's ancestor is the root, so its reads stay
        // unexplained; family 40's reads move up alone and are too few.
        assert!(found.rows.is_empty());
        assert_eq!(found.unexplained, 250);
        let mut not_kept = reads(&[(1001, 1, 0.99), (1011, 11, 0.99)], 150, true);
        not_kept.extend(reads(&[(1021, 21, 0.99), (1001, 1, 0.99)], 100, true));
        let found = place_reads(Vec::new(), not_kept, &none);
        assert_eq!(rows(&found), [(30, Rank::Order, 250.0)]);
        assert_eq!(found.end_of(40).map(|row| row.taxid), Some(30));
    }

    #[test]
    fn a_row_s_depth_is_its_reads_over_the_genome_they_fit_best() {
        let depth = |found: &HigherLevels| found.rows[0].depth;
        // Genus 100: 150 reads on genome 1 (1 Mb) and 100 on genome 5 (100 kb). All 250 count over genome 1,
        // which holds most of them.
        let mut not_kept = reads(&[(1001, 1, 0.99)], 150, true);
        not_kept.extend((0..100).map(|i| read_small(&[(1002, 5, 0.99)], i * 990, 1.0)));
        let found = place_reads(Vec::new(), not_kept, &HashSet::new());
        assert_eq!(rows(&found), [(100, Rank::Genus, 250.0)]);
        assert_eq!(found.rows[0].measured.genome, 1);
        assert!((depth(&found) - 250.0 * 100.0 / 1e6).abs() < 1e-12);
        // Reads tied across genera 100 and 101, 125 best on genome 5 (100 kb) and 125 on genome 6 (1 Mb): family
        // 40 is measured on genome 5 and its depth is the mean of the depths over the two genomes.
        let mut across = (0..125)
            .map(|i| read_small(&[(1001, 5, 0.99), (1011, 11, 0.99)], i * 790, 1.0))
            .collect::<Vec<_>>();
        across.extend(reads(&[(1001, 6, 0.99), (1011, 11, 0.99)], 125, true));
        let found = place_reads(Vec::new(), across, &HashSet::new());
        assert_eq!(rows(&found), [(40, Rank::Family, 250.0)]);
        assert_eq!(found.rows[0].measured.genome, 5);
        assert!((depth(&found) - 250.0 * 100.0 * (1.0 / 1e5 + 1.0 / 1e6) / 2.0).abs() < 1e-12);
        // On a genome indexed in sampled segments, reads count as full-index equivalents and the depth is taken
        // over the indexed length.
        let sampled = (0..250)
            .map(|i| read_small(&[(1001, 5, 0.99)], i * 390, 10.0))
            .collect();
        let found = place_reads(Vec::new(), sampled, &HashSet::new());
        assert_eq!(rows(&found), [(100, Rank::Genus, 250.0)]);
        assert!((found.rows[0].equivalent - 2_500.0).abs() < 1e-9);
        assert!((depth(&found) - 250.0 * 100.0 / 1e5).abs() < 1e-12);
    }

    #[test]
    fn on_a_sampled_database_taxa_holding_a_kept_species_get_no_row() {
        // Genus 100 and family 40 hold a kept species: the genus fails, and its reads move past the family to
        // the order.
        let found = place_reads(
            Vec::new(),
            reads(&[(1001, 1, 0.99)], 250, true),
            &HashSet::from([100, 40]),
        );
        assert_eq!(found.genera[0].failures, [Failure::HoldsKept]);
        assert_eq!(found.taxa[0].failures, [Failure::HoldsKept]);
        assert_eq!(rows(&found), [(30, Rank::Order, 250.0)]);
    }
}

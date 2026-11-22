//! Reads of species removed by the coverage evidence, reported at their genus.
//!
//! A fragment is orphaned when the coverage evidence rejects every species it
//! maps to. Often such a species is the reference of a relative that matches
//! the sample organism only in part: a divergent strain, or a new species of a
//! genus the database holds. Its reads still belong to that genus. Reads that
//! land on foreign contigs of a contaminated reference, on a misassembled
//! region, or on DNA shared across genera also leave rejected species behind,
//! but they reach only part of the reference. Two spread checks on the
//! species' own evidence keep those out: the chunk ratio (COEBR) and the share
//! of genome windows that its reads reach.
//!
//! A removed species whose reads cover its whole genome evenly, but too thinly
//! for the amount tests, is reported at species instead, when no taxon of its
//! genus is kept.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use super::evidence::{EvidenceCoverage, PresenceInterval, WitnessDecision};
use super::ProfileRank;

/// Lowest chunk ratio (COEBR) of a removed species whose reads go to its genus.
pub(crate) const MIN_CHUNK_RATIO: f64 = 0.40;
/// Lowest share of genome windows that its reads must reach.
pub(crate) const MIN_WINDOW_SPREAD: f64 = 0.50;
/// Fewest fragments a genus row needs; the level-wise engine asks as many of
/// a family row.
pub(crate) const MIN_GENUS_FRAGMENTS: f64 = 200.0;

/// Fewest fragments of a removed species reported at species for its even
/// coverage.
pub(crate) const EVEN_MIN_FRAGMENTS: usize = 200;
/// Both coverage ratios (OEBR and COEBR) of such a species lie in this range.
pub(crate) const EVEN_RATIO_MIN: f64 = 0.9;
pub(crate) const EVEN_RATIO_MAX: f64 = 1.1;
/// Lowest share of genome windows that its reads reach.
pub(crate) const EVEN_MIN_WINDOW_SPREAD: f64 = 0.9;
/// Lowest breadth and chunk p-values.
pub(crate) const EVEN_MIN_PVALUE: f64 = 1e-5;

/// At most this many windows; fewer when there are fewer reads, so that every
/// window could receive about [`FRAGMENTS_PER_WINDOW`] of them.
const MAX_WINDOWS: usize = 100;
const FRAGMENTS_PER_WINDOW: usize = 5;
/// Below this many windows the spread is not measured, and the species fails.
const MIN_WINDOWS: usize = 4;

/// Position of every reference within its genome, in database order, and the
/// length of every genome. The window spread is measured on this coordinate.
#[derive(Debug)]
pub(crate) struct GenomeLayout {
    offsets: Vec<u64>,
    lengths: Vec<u64>,
}

impl GenomeLayout {
    /// `references` lists `(assembly id, reference length)` in reference-id
    /// order; assembly ids are dense and start at zero.
    pub(crate) fn new(references: impl IntoIterator<Item = (u32, u64)>) -> Self {
        let mut offsets = Vec::new();
        let mut lengths: Vec<u64> = Vec::new();
        for (assembly, length) in references {
            let assembly = assembly as usize;
            if assembly >= lengths.len() {
                lengths.resize(assembly + 1, 0);
            }
            offsets.push(lengths[assembly]);
            lengths[assembly] += length;
        }
        Self { offsets, lengths }
    }

    /// Share of equal windows along the genome that hold the start of at least
    /// one of `intervals`. `None` when there are too few reads to measure it.
    pub(crate) fn window_spread(
        &self,
        assembly: u32,
        intervals: &[PresenceInterval],
    ) -> Option<f64> {
        let windows = (intervals.len() / FRAGMENTS_PER_WINDOW).min(MAX_WINDOWS);
        let length = *self.lengths.get(assembly as usize)?;
        if windows < MIN_WINDOWS || length == 0 {
            return None;
        }
        let mut hit = 0u128;
        for interval in intervals {
            let position =
                self.offsets[interval.reference_id as usize] + u64::from(interval.span.0);
            let window = (u128::from(position) * windows as u128 / u128::from(length)) as usize;
            hit |= 1u128 << window.min(windows - 1);
        }
        Some(f64::from(hit.count_ones()) / windows as f64)
    }
}

/// Database-wide inputs of the genus roll-up.
#[derive(Debug)]
pub(crate) struct GenusRollup {
    /// Genus of every species with genome references. Species represented only
    /// by marker genes are absent: reads failing them come mostly from the
    /// homologous genes of relatives.
    pub(crate) genus_of: HashMap<u32, u32>,
    /// Species represented only by marker genes, with their genus when known.
    /// They never take part, but a kept one keeps its genus from an
    /// even-coverage species.
    pub(crate) marker_only: HashMap<u32, Option<u32>>,
    pub(crate) layout: GenomeLayout,
    /// A fractional database: a genus with an accepted species gets no row.
    /// Reads from the unindexed parts of a genome land on the sampled segments
    /// of relatives and leave them orphaned, while the accepted species'
    /// abundance is already scaled to its full genome, so a genus row would
    /// count those reads a second time.
    pub(crate) fractional: bool,
}

impl GenusRollup {
    /// The genus of `taxids` when every one of them lies in it.
    pub(crate) fn common_genus(&self, taxids: impl IntoIterator<Item = u32>) -> Option<u32> {
        let mut genera = taxids.into_iter().map(|taxid| self.genus_of.get(&taxid));
        let genus = *genera.next()??;
        genera.all(|other| other == Some(&genus)).then_some(genus)
    }

    /// Species that the coverage evidence removed although their reads cover
    /// their genome evenly: at least [`EVEN_MIN_FRAGMENTS`] fragments, OEBR
    /// and COEBR within [`EVEN_RATIO_MIN`]-[`EVEN_RATIO_MAX`], a window spread
    /// of at least [`EVEN_MIN_WINDOW_SPREAD`], both p-values at least
    /// [`EVEN_MIN_PVALUE`], and fixed chunk breadth below the cutoff
    /// `chunk_cutoff`. A species takes part only when no taxon of its genus
    /// was accepted. Species represented only by marker genes never take
    /// part, and none does on a fractional database. Sorted by taxid.
    pub(crate) fn even_coverage_species(
        &self,
        decisions: &HashMap<u32, WitnessDecision>,
        rank_of: impl Fn(u32) -> Option<ProfileRank>,
        chunk_cutoff: f64,
    ) -> Vec<u32> {
        if self.fractional {
            return Vec::new();
        }
        let genus = |taxid: u32| match rank_of(taxid) {
            Some(ProfileRank::Genus) => Some(taxid),
            _ => self
                .genus_of
                .get(&taxid)
                .copied()
                .or_else(|| self.marker_only.get(&taxid).copied().flatten()),
        };
        let kept: HashSet<u32> = decisions
            .iter()
            .filter(|(_, decision)| decision.accepted())
            .filter_map(|(&taxid, _)| genus(taxid))
            .collect();
        let mut species = decisions
            .iter()
            .filter(|&(&taxid, decision)| {
                !decision.accepted()
                    && rank_of(taxid) == Some(ProfileRank::Species)
                    && !self.marker_only.contains_key(&taxid)
                    && has_even_coverage(&decision.metrics, chunk_cutoff)
                    && genus(taxid).is_none_or(|genus| !kept.contains(&genus))
            })
            .map(|(&taxid, _)| taxid)
            .collect::<Vec<_>>();
        species.sort_unstable();
        species
    }
}

fn has_even_coverage(metrics: &EvidenceCoverage, chunk_cutoff: f64) -> bool {
    let even = |ratio: f64| (EVEN_RATIO_MIN..=EVEN_RATIO_MAX).contains(&ratio);
    let high = |pvalue: Option<f64>| pvalue.is_some_and(|pvalue| pvalue >= EVEN_MIN_PVALUE);
    metrics.count >= EVEN_MIN_FRAGMENTS
        && even(metrics.oebr)
        && even(metrics.coebr)
        && metrics
            .window_spread
            .is_some_and(|spread| spread >= EVEN_MIN_WINDOW_SPREAD)
        && high(metrics.breadth_pvalue)
        && high(metrics.chunk_pvalue)
        && metrics.fixed_chunk < chunk_cutoff
}

/// One genus row made from orphaned fragments.
#[derive(Debug)]
pub(crate) struct GenusRow {
    pub(crate) genus: u32,
    pub(crate) fragments: f64,
    /// Fragments weighted by the sampling scale of their genome.
    pub(crate) equivalent: f64,
    pub(crate) depth: f64,
    pub(crate) mapped_len: f64,
    /// Coverage of the donor species with the most fragments; reported in the
    /// profile row, as the coherent witness is for a majority-LCA row.
    pub(crate) witness: EvidenceCoverage,
    /// Removed species whose fragments make up the row, by taxid.
    pub(crate) donors: Vec<u32>,
    /// Read names, kept only for `classify.txt.gz`.
    pub(crate) reads: Option<Vec<String>>,
}

#[derive(Debug, Default)]
struct Donor {
    genus: u32,
    fragments: f64,
    equivalent: f64,
    depth: f64,
    mapped_len: f64,
    reads: Vec<String>,
}

/// Orphaned fragments, grouped by the removed species of their best placement.
#[derive(Debug)]
pub(crate) struct OrphanedFragments<'a> {
    rollup: &'a GenusRollup,
    keep_reads: bool,
    donors: HashMap<u32, Donor>,
}

impl<'a> OrphanedFragments<'a> {
    pub(crate) fn new(rollup: &'a GenusRollup, keep_reads: bool) -> Self {
        Self {
            rollup,
            keep_reads,
            donors: HashMap::new(),
        }
    }

    /// Record a fragment none of whose candidate species was kept. It counts
    /// only when every candidate lies in one genus; `best` indexes the
    /// candidate of its best placement, whose depth, mapped length and
    /// sampling scale it carries.
    pub(crate) fn offer(
        &mut self,
        taxids: &[u32],
        best: usize,
        depth: f64,
        mapped_len: f64,
        sampling_scale: f64,
        read: &Arc<str>,
    ) {
        let Some(genus) = self.rollup.common_genus(taxids.iter().copied()) else {
            return;
        };
        let donor = self.donors.entry(taxids[best]).or_default();
        donor.genus = genus;
        donor.fragments += 1.0;
        donor.equivalent += sampling_scale;
        donor.depth += depth;
        donor.mapped_len += mapped_len;
        if self.keep_reads {
            donor.reads.push(read.to_string());
        }
    }

    /// Genus rows from the donors whose own evidence spreads over the genome,
    /// in genus order. A donor without evidence, or whose window spread could
    /// not be measured, does not pass. `accepted` holds the taxa the coverage
    /// evidence kept; on a fractional database their genera get no row.
    pub(crate) fn into_genus_rows(
        self,
        decisions: &HashMap<u32, WitnessDecision>,
        accepted: &HashSet<u32>,
    ) -> Vec<GenusRow> {
        let genus_of = &self.rollup.genus_of;
        let represented: HashSet<u32> = if self.rollup.fractional {
            accepted
                .iter()
                .filter_map(|taxid| genus_of.get(taxid).copied())
                .collect()
        } else {
            HashSet::new()
        };
        let mut by_genus: HashMap<u32, Vec<(u32, Donor)>> = HashMap::new();
        for (taxid, donor) in self.donors {
            if represented.contains(&donor.genus) {
                continue;
            }
            let Some(metrics) = decisions.get(&taxid).map(|decision| decision.metrics) else {
                continue;
            };
            let spread = metrics.window_spread.unwrap_or(0.0);
            if metrics.coebr >= MIN_CHUNK_RATIO && spread >= MIN_WINDOW_SPREAD {
                by_genus
                    .entry(donor.genus)
                    .or_default()
                    .push((taxid, donor));
            }
        }
        let keep_reads = self.keep_reads;
        let mut rows = by_genus
            .into_iter()
            .filter_map(|(genus, mut donors)| {
                let fragments = donors.iter().map(|(_, donor)| donor.fragments).sum::<f64>();
                if fragments < MIN_GENUS_FRAGMENTS {
                    return None;
                }
                // Most fragments first, then taxid: a stable witness and order.
                donors.sort_unstable_by(|(left, a), (right, b)| {
                    b.fragments.total_cmp(&a.fragments).then(left.cmp(right))
                });
                let witness = decisions[&donors[0].0].metrics;
                let mut row = GenusRow {
                    genus,
                    fragments,
                    equivalent: 0.0,
                    depth: 0.0,
                    mapped_len: 0.0,
                    witness,
                    donors: Vec::with_capacity(donors.len()),
                    reads: keep_reads.then(Vec::new),
                };
                for (taxid, donor) in donors {
                    row.equivalent += donor.equivalent;
                    row.depth += donor.depth;
                    row.mapped_len += donor.mapped_len;
                    row.donors.push(taxid);
                    if let Some(reads) = row.reads.as_mut() {
                        reads.extend(donor.reads);
                    }
                }
                Some(row)
            })
            .collect::<Vec<_>>();
        rows.sort_unstable_by_key(|row| row.genus);
        rows
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(reference_id: u32, start: u32) -> PresenceInterval {
        PresenceInterval {
            reference_id,
            span: (start, start + 150),
            mate_span: (0, 0),
            fixed_chunk: (start, start + 150),
            flex_chunk: (start, start + 150),
        }
    }

    fn decision(coebr: f64, window_spread: Option<f64>) -> WitnessDecision {
        WitnessDecision {
            metrics: EvidenceCoverage {
                coebr,
                window_spread,
                ..EvidenceCoverage::default()
            },
            standard: false,
        }
    }

    /// A removed species whose reads cover its genome evenly.
    fn even(accepted: bool) -> WitnessDecision {
        WitnessDecision {
            metrics: EvidenceCoverage {
                count: 300,
                oebr: 1.0,
                coebr: 1.0,
                window_spread: Some(0.95),
                breadth_pvalue: Some(0.5),
                chunk_pvalue: Some(0.5),
                fixed_chunk: 0.1,
                ..EvidenceCoverage::default()
            },
            standard: accepted,
        }
    }

    #[test]
    fn a_removed_species_with_even_coverage_is_reported_at_species() {
        let rollup = rollup(&[(1, 100), (2, 200), (3, 300)]);
        let ranks = |taxid: u32| {
            Some(if taxid >= 100 {
                ProfileRank::Genus
            } else {
                ProfileRank::Species
            })
        };
        let decisions = HashMap::from([(1, even(false)), (2, even(false))]);
        assert_eq!(
            rollup.even_coverage_species(&decisions, ranks, 0.3),
            vec![1, 2]
        );

        type Change = fn(&mut EvidenceCoverage);
        let failing: [(&str, Change); 10] = [
            ("fragments", |m| m.count = 199),
            ("low OEBR", |m| m.oebr = 0.89),
            ("high OEBR", |m| m.oebr = 1.11),
            ("high COEBR", |m| m.coebr = 1.11),
            ("spread", |m| m.window_spread = Some(0.89)),
            ("no spread", |m| m.window_spread = None),
            ("breadth p-value", |m| m.breadth_pvalue = Some(9e-6)),
            ("no chunk p-value", |m| m.chunk_pvalue = None),
            ("chunk amount", |m| m.fixed_chunk = 0.3),
            ("low COEBR", |m| m.coebr = 0.89),
        ];
        for (what, change) in failing {
            let mut decision = even(false);
            change(&mut decision.metrics);
            let decisions = HashMap::from([(1, decision)]);
            assert!(
                rollup
                    .even_coverage_species(&decisions, ranks, 0.3)
                    .is_empty(),
                "{what}"
            );
        }
        let mut boundary = even(false);
        boundary.metrics.count = 200;
        boundary.metrics.oebr = 0.9;
        boundary.metrics.coebr = 1.1;
        boundary.metrics.window_spread = Some(0.9);
        boundary.metrics.breadth_pvalue = Some(1e-5);
        let decisions = HashMap::from([(1, boundary)]);
        assert_eq!(
            rollup.even_coverage_species(&decisions, ranks, 0.3),
            vec![1]
        );
    }

    #[test]
    fn a_kept_taxon_of_the_genus_blocks_even_coverage_species() {
        let mut rollup = rollup(&[(1, 100), (2, 100), (3, 300)]);
        rollup.marker_only = HashMap::from([(4, Some(300)), (5, None)]);
        let ranks = |taxid: u32| {
            Some(if taxid >= 100 {
                ProfileRank::Genus
            } else {
                ProfileRank::Species
            })
        };
        // A kept species of genus 100 keeps species 1 out; species 3 of genus
        // 300 is blocked by a kept marker-only species of that genus.
        let decisions = HashMap::from([
            (1, even(false)),
            (2, even(true)),
            (3, even(false)),
            (4, even(true)),
        ]);
        assert!(rollup
            .even_coverage_species(&decisions, ranks, 0.3)
            .is_empty());
        // A kept taxon at genus rank blocks its genus too.
        let decisions = HashMap::from([(1, even(false)), (100, even(true))]);
        assert!(rollup
            .even_coverage_species(&decisions, ranks, 0.3)
            .is_empty());
        // Marker-only species and taxa above species never take part.
        let decisions = HashMap::from([(5, even(false)), (300, even(false))]);
        assert!(rollup
            .even_coverage_species(&decisions, ranks, 0.3)
            .is_empty());
        // A species without a genus is not blocked.
        let decisions = HashMap::from([(7, even(false)), (2, even(true))]);
        assert_eq!(
            rollup.even_coverage_species(&decisions, ranks, 0.3),
            vec![7]
        );
        // None on a fractional database.
        rollup.fractional = true;
        let decisions = HashMap::from([(1, even(false))]);
        assert!(rollup
            .even_coverage_species(&decisions, ranks, 0.3)
            .is_empty());
    }

    fn rollup(genus_of: &[(u32, u32)]) -> GenusRollup {
        GenusRollup {
            genus_of: genus_of.iter().copied().collect(),
            marker_only: HashMap::new(),
            layout: GenomeLayout::new([(0, 1_000)]),
            fractional: false,
        }
    }

    #[test]
    fn layout_places_references_within_their_genome() {
        // Two genomes whose references are interleaved in the database.
        let layout = GenomeLayout::new([(0, 600), (1, 100), (0, 400)]);
        assert_eq!(layout.offsets, vec![0, 0, 600]);
        assert_eq!(layout.lengths, vec![1_000, 100]);
        // 20 reads -> 4 windows of 250 bp; starts at 0, 600 and 650 (reference 2
        // at 600 + 50) fall into windows 0 and 2.
        let mut reads = vec![interval(0, 0); 18];
        reads.push(interval(2, 0));
        reads.push(interval(2, 50));
        assert_eq!(layout.window_spread(0, &reads), Some(0.5));
    }

    #[test]
    fn window_spread_needs_enough_reads_and_caps_the_windows() {
        let layout = GenomeLayout::new([(0, 1_000_000)]);
        assert_eq!(layout.window_spread(0, &vec![interval(0, 0); 19]), None);
        let even = (0..1_000)
            .map(|i| interval(0, i * 1_000))
            .collect::<Vec<_>>();
        assert_eq!(layout.window_spread(0, &even), Some(1.0));
        // A start at the very end of the genome stays in the last window.
        let end = vec![interval(0, 999_999); 1_000];
        assert_eq!(layout.window_spread(0, &end), Some(0.01));
        assert_eq!(layout.window_spread(7, &even), None);
    }

    #[test]
    fn orphans_of_one_genus_become_a_genus_row() {
        let rollup = rollup(&[(1, 10), (2, 10), (3, 20)]);
        let mut orphans = OrphanedFragments::new(&rollup, true);
        let read: Arc<str> = Arc::from("r");
        for _ in 0..150 {
            orphans.offer(&[1], 0, 0.01, 150.0, 1.0, &read);
        }
        for _ in 0..60 {
            orphans.offer(&[1, 2], 1, 0.02, 150.0, 2.0, &read);
        }
        // Candidates in two genera, or one without a genus, are left out.
        orphans.offer(&[1, 3], 0, 1.0, 1.0, 1.0, &read);
        orphans.offer(&[4], 0, 1.0, 1.0, 1.0, &read);
        let decisions =
            HashMap::from([(1, decision(0.9, Some(0.8))), (2, decision(0.6, Some(0.6)))]);
        let rows = orphans.into_genus_rows(&decisions, &HashSet::new());
        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert_eq!(
            (row.genus, row.fragments, row.donors.clone()),
            (10, 210.0, vec![1, 2])
        );
        assert!((row.depth - (1.5 + 1.2)).abs() < 1e-9);
        assert_eq!(row.equivalent, 270.0);
        assert_eq!(row.witness.coebr, 0.9);
        assert_eq!(row.reads.as_ref().map(Vec::len), Some(210));
    }

    #[test]
    fn donors_need_both_spreads_and_the_genus_enough_fragments() {
        let rollup = rollup(&[(1, 10), (2, 10), (3, 10), (5, 30)]);
        let mut orphans = OrphanedFragments::new(&rollup, false);
        let read: Arc<str> = Arc::from("r");
        for taxid in [1, 2, 3, 4, 5] {
            for _ in 0..300 {
                orphans.offer(&[taxid], 0, 0.01, 150.0, 1.0, &read);
            }
        }
        let decisions = HashMap::from([
            (1, decision(0.39, Some(1.0))), // chunk ratio too low
            (2, decision(0.9, Some(0.49))), // reads on too small a part of the genome
            (3, decision(0.9, None)),       // spread not measured
            (5, decision(0.40, Some(0.50))),
        ]);
        let rows = orphans.into_genus_rows(&decisions, &HashSet::new());
        assert_eq!(
            rows.iter().map(|row| row.genus).collect::<Vec<_>>(),
            vec![30]
        );
        assert!(rows[0].reads.is_none());

        let mut few = OrphanedFragments::new(&rollup, false);
        for _ in 0..199 {
            few.offer(&[5], 0, 0.01, 150.0, 1.0, &read);
        }
        assert!(few.into_genus_rows(&decisions, &HashSet::new()).is_empty());
    }

    #[test]
    fn a_fractional_database_rolls_up_only_genera_without_an_accepted_species() {
        // Species 1 and 2 are removed, in genera 10 and 20; species 3 of genus
        // 10 is accepted.
        let rows_of = |fractional| {
            let rollup = GenusRollup {
                fractional,
                ..rollup(&[(1, 10), (2, 20), (3, 10)])
            };
            let mut orphans = OrphanedFragments::new(&rollup, false);
            let read: Arc<str> = Arc::from("r");
            for taxid in [1, 2] {
                for _ in 0..300 {
                    orphans.offer(&[taxid], 0, 0.01, 150.0, 10.0, &read);
                }
            }
            let decisions =
                HashMap::from([(1, decision(0.9, Some(1.0))), (2, decision(0.9, Some(1.0)))]);
            orphans
                .into_genus_rows(&decisions, &HashSet::from([3]))
                .iter()
                .map(|row| row.genus)
                .collect::<Vec<_>>()
        };
        assert_eq!(rows_of(true), vec![20]);
        assert_eq!(rows_of(false), vec![10, 20]);
    }
}

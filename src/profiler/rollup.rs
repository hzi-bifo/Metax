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

use std::collections::HashMap;
use std::sync::Arc;

use super::evidence::{EvidenceCoverage, PresenceInterval, WitnessDecision};

/// Lowest chunk ratio (COEBR) of a removed species whose reads go to its genus.
pub(crate) const MIN_CHUNK_RATIO: f64 = 0.40;
/// Lowest share of genome windows that its reads must reach.
pub(crate) const MIN_WINDOW_SPREAD: f64 = 0.50;
/// Fewest fragments a genus row needs.
pub(crate) const MIN_GENUS_FRAGMENTS: f64 = 200.0;

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
    pub(crate) layout: GenomeLayout,
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
    /// Read names, kept only for `classify.txt`.
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
        let mut genera = taxids.iter().map(|taxid| self.rollup.genus_of.get(taxid));
        let Some(Some(&genus)) = genera.next() else {
            return;
        };
        if !genera.all(|other| other == Some(&genus)) {
            return;
        }
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
    /// not be measured, does not pass.
    pub(crate) fn into_genus_rows(
        self,
        decisions: &HashMap<u32, WitnessDecision>,
    ) -> Vec<GenusRow> {
        let mut by_genus: HashMap<u32, Vec<(u32, Donor)>> = HashMap::new();
        for (taxid, donor) in self.donors {
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
    use crate::profiler::evidence::{PremapChunkGuard, RescueKinds};

    fn interval(reference_id: u32, start: u32) -> PresenceInterval {
        PresenceInterval {
            reference_id,
            span: (start, start + 150),
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
            rescues: RescueKinds::default(),
            required_premap_guard: PremapChunkGuard::None,
        }
    }

    fn rollup(genus_of: &[(u32, u32)]) -> GenusRollup {
        GenusRollup {
            genus_of: genus_of.iter().copied().collect(),
            layout: GenomeLayout::new([(0, 1_000)]),
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
        let rows = orphans.into_genus_rows(&decisions);
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
        let rows = orphans.into_genus_rows(&decisions);
        assert_eq!(
            rows.iter().map(|row| row.genus).collect::<Vec<_>>(),
            vec![30]
        );
        assert!(rows[0].reads.is_none());

        let mut few = OrphanedFragments::new(&rollup, false);
        for _ in 0..199 {
            few.offer(&[5], 0, 0.01, 150.0, 1.0, &read);
        }
        assert!(few.into_genus_rows(&decisions).is_empty());
    }
}

//! Database-wide inputs of the levels above species: the genome layout on
//! which the window spread of reads is measured, the genus of every species
//! with genome references, and whether the database is fractional.
//!
//! A species that the coverage evidence removed although its reads cover its
//! whole genome evenly, but too thinly for the amount tests, is reported at
//! species instead, when no taxon of its genus is kept
//! ([`GenusRollup::even_coverage_species`]).

use std::collections::{HashMap, HashSet};

use super::evidence::{EvidenceCoverage, PresenceInterval, WitnessDecision};
use super::ProfileRank;

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

/// Database-wide inputs of the levels above species and of the even-coverage
/// species.
#[derive(Debug)]
pub(crate) struct GenusRollup {
    /// Genus of every species with genome references. Species represented only
    /// by marker genes are absent.
    pub(crate) genus_of: HashMap<u32, u32>,
    /// Species represented only by marker genes, with their genus when known.
    /// They never take part, but a kept one keeps its genus from an
    /// even-coverage species.
    pub(crate) marker_only: HashMap<u32, Option<u32>>,
    pub(crate) layout: GenomeLayout,
    /// A database of sampled genomes (fractional): reads not placed at a
    /// species are not reported above species, and no species is kept for
    /// even coverage.
    pub(crate) fractional: bool,
}

impl GenusRollup {
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
}

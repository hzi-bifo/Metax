//! Coordinate evidence of native Metam mappings: each fragment's candidate
//! taxa, with one placement per genome within the read-level gates
//! ([`CandidateScratch`]), and the coverage measurements of a genome's
//! intervals ([`calculate_coverage`]).

use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::HashMap;

use metam::format::canonical_mapping_identity;

use super::statistics::{calc_breadth_pvalue, Tail, NEUTRAL_RATIO};
use super::ProfilerConfig;

const SCORE_EPSILON: f64 = 1e-9;
const FRACTION_F32_EPSILON: f64 = f32::EPSILON as f64;

/// Role carried explicitly by one native mapping row.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MappingRole {
    /// A single-end or long-read mapping.
    Unpaired,
    /// R1 in an explicitly identified paired placement.
    Mate1,
    /// R2 in an explicitly identified paired placement.
    Mate2,
}

/// Compact representation of one mapping of a read.
#[derive(Clone, Debug)]
pub(crate) struct EvidenceAlignment {
    pub(crate) profile_taxid: u32,
    pub(crate) reference_id: u32,
    pub(crate) assembly_id: u32,
    pub(crate) role: MappingRole,
    /// Read-local mapper-selected pair identifier; zero for unpaired rows.
    pub(crate) pair_id: u32,
    pub(crate) sampling_scale: f32,
    pub(crate) asm_len: f64,
    pub(crate) identity: f32,
    pub(crate) fraction: f32,
    pub(crate) mapped_len: u32,
    /// Length of the read, aligned or not.
    pub(crate) read_len: u32,
    pub(crate) span: (u32, u32),
    pub(crate) fixed_chunk: (u32, u32),
    pub(crate) flex_chunk: (u32, u32),
    /// Length of the reference when it is a marker gene; 0 otherwise.
    pub(crate) marker_gene_len: u32,
}

impl EvidenceAlignment {
    /// The alignment reaches the first or last base of a marker gene, so the
    /// read may continue beyond the gene.
    pub(crate) fn reaches_marker_end(&self) -> bool {
        self.marker_gene_len > 0
            && (self.span.0 == 1 || self.span.1 == self.marker_gene_len.saturating_add(1))
    }

    /// The read aligns to more than half of a marker gene.
    fn aligns_over_half_of_marker_gene(&self) -> bool {
        self.marker_gene_len > 0 && 2 * u64::from(self.mapped_len) > u64::from(self.marker_gene_len)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum EvidenceMode {
    Paired,
    Single,
    Long,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct EvidencePolicy {
    pub(crate) mode: EvidenceMode,
    pub(crate) identity: f64,
    pub(crate) mapped_len: usize,
    pub(crate) fraction: f64,
}

/// Fewest mapped bases (inclusive) of a short read at the floor gates
/// ([`EvidencePolicy::floor`]).
const FLOOR_MIN_MAPPED_LEN: usize = 50;
/// Lowest aligned fraction (inclusive) of a read at the floor gates.
const FLOOR_MIN_FRACTION: f64 = 0.60;

impl EvidencePolicy {
    /// The floor gates of the fragments without a species-level candidate,
    /// for the levels above species and the noisy long reads credited to kept
    /// species: `--family-identity` (`identity`), [`FLOOR_MIN_MAPPED_LEN`]
    /// mapped bases of a short read (a long read keeps the species length
    /// cutoff), and [`FLOOR_MIN_FRACTION`] of a read aligned.
    pub(crate) fn floor(self, identity: f64) -> Self {
        Self {
            identity,
            mapped_len: match self.mode {
                EvidenceMode::Long => self.mapped_len,
                EvidenceMode::Paired | EvidenceMode::Single => FLOOR_MIN_MAPPED_LEN,
            },
            fraction: FLOOR_MIN_FRACTION,
            ..self
        }
    }

    pub(crate) fn from_config(cfg: &ProfilerConfig) -> Option<Self> {
        let mode = if cfg.is_paired && cfg.sequencer.is_illumina() {
            EvidenceMode::Paired
        } else if !cfg.is_paired && cfg.sequencer.is_illumina() {
            EvidenceMode::Single
        } else if !cfg.is_paired && cfg.sequencer.is_long_read() {
            EvidenceMode::Long
        } else {
            return None;
        };
        Some(Self {
            mode,
            identity: cfg.identity,
            mapped_len: cfg.mapped_len,
            fraction: cfg.fraction,
        })
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct CandidatePlacement {
    pub(crate) taxid: u32,
    pub(crate) reference_id: u32,
    pub(crate) assembly_id: u32,
    pub(crate) sampling_scale: f64,
    pub(crate) asm_len: f64,
    pub(crate) r1_mapped: u32,
    pub(crate) both_mapped: u32,
    /// Length of the read that gives `r1_mapped`.
    pub(crate) read_len: u32,
    /// Other marker genes of the same genome that the fragment covers
    /// ([`TerminalCandidate::other_gene_hits`]); 0 otherwise.
    pub(crate) other_genes: u32,
    pub(crate) pair_errors: f64,
    pub(crate) long_score: f64,
    pub(crate) span: (u32, u32),
    /// Span of the other mate of a pair on the same reference; empty for a
    /// single read.
    pub(crate) mate_span: (u32, u32),
    pub(crate) fixed_chunk: (u32, u32),
    pub(crate) flex_chunk: (u32, u32),
    /// Stable order of the representative mapping in the mapper output.
    /// This is biological tie-preserving metadata, unlike a numeric reference ID.
    pub(crate) input_order: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct TerminalCandidate {
    pub(crate) taxid: u32,
    pub(crate) best: CandidatePlacement,
    /// One independently evaluated placement per represented assembly.
    pub(crate) placements: Vec<CandidatePlacement>,
    /// The other marker genes that the fragment covers on the genomes of
    /// `placements`, one per gene.
    pub(crate) other_gene_hits: Vec<GeneHit>,
}

/// A marker gene that a fragment covers besides the one of its placement.
#[derive(Clone, Copy, Debug)]
pub(crate) struct GeneHit {
    pub(crate) assembly_id: u32,
    pub(crate) interval: PresenceInterval,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct PresenceInterval {
    pub(crate) reference_id: u32,
    pub(crate) span: (u32, u32),
    /// The other mate of a pair; empty for a single read.
    pub(crate) mate_span: (u32, u32),
    pub(crate) fixed_chunk: (u32, u32),
    pub(crate) flex_chunk: (u32, u32),
}

impl From<&CandidatePlacement> for PresenceInterval {
    fn from(placement: &CandidatePlacement) -> Self {
        Self {
            reference_id: placement.reference_id,
            span: placement.span,
            mate_span: placement.mate_span,
            fixed_chunk: placement.fixed_chunk,
            flex_chunk: placement.flex_chunk,
        }
    }
}

impl TerminalCandidate {
    /// The candidate's placement on genome `assembly_id`, if it has one.
    pub(crate) fn placement_on(&self, assembly_id: u32) -> Option<&CandidatePlacement> {
        self.placements
            .iter()
            .find(|p| p.assembly_id == assembly_id)
    }
}

#[derive(Clone, Debug, Default)]
pub(crate) struct PresenceAssembly {
    /// Effective indexed assembly size. This is the coordinate space that can
    /// actually receive mappings in a fractional index.
    pub(crate) asm_len: f64,
    pub(crate) count: u64,
    /// Marker genes hit beyond one per fragment: the other genes that
    /// fragments cover on this genome.
    pub(crate) other_gene_hits: u64,
    pub(crate) reference_span_bases: u64,
    pub(crate) intervals: Vec<PresenceInterval>,
    /// Resolved fragments whose best placement within their taxon is on this
    /// genome.
    pub(crate) fits: u64,
    /// Full-index equivalents of one fragment on this genome (§4.2).
    pub(crate) sampling_scale: f64,
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct EvidenceCoverage {
    pub(crate) count: usize,
    pub(crate) depth: f64,
    pub(crate) breadth: f64,
    pub(crate) expected_breadth: f64,
    pub(crate) fixed_chunk: f64,
    pub(crate) flex_chunk: f64,
    pub(crate) expected_flex_chunk: f64,
    pub(crate) oebr: f64,
    pub(crate) coebr: f64,
    pub(crate) breadth_pvalue: Option<f64>,
    pub(crate) chunk_pvalue: Option<f64>,
    /// Distinct marker genes hit, for a genome represented by marker genes.
    pub(crate) marker_genes_hit: Option<usize>,
    /// Share of genome windows that the reads reach; measured only for
    /// the even-coverage test and when there are enough reads.
    pub(crate) window_spread: Option<f64>,
}

impl PresenceAssembly {
    /// Add a fragment's `interval` and the bases of its read's span.
    pub(crate) fn add_interval(&mut self, interval: PresenceInterval) {
        self.reference_span_bases += u64::from(interval.span.1.saturating_sub(interval.span.0));
        self.intervals.push(interval);
    }
}

/// Reduce one assembly in-place. Reusing the same interval vector for the
/// three sorts avoids three full-size temporary vectors and keeps peak memory
/// proportional to the stored evidence itself.
pub(crate) fn calculate_coverage(
    assembly: &mut PresenceAssembly,
    multiple_testing_taxa: usize,
) -> EvidenceCoverage {
    let genome_size = assembly.asm_len.max(1.0);
    let (span_bases, second_mate_bases) = read_bases(&mut assembly.intervals);
    let fixed_bases = merged_interval_bases(&mut assembly.intervals, |record| record.fixed_chunk);
    let flex_bases = merged_interval_bases(&mut assembly.intervals, |record| record.flex_chunk);
    let count = usize::try_from(assembly.count).unwrap_or(usize::MAX);
    // Breadth counts the bases of both mates of a pair, and so do the depth
    // its expectation and p-value come from. The reported depth counts each
    // fragment once, by its first read, as abundance does.
    let depth = assembly.reference_span_bases as f64 / genome_size;
    let read_depth = (assembly.reference_span_bases + second_mate_bases) as f64 / genome_size;
    let breadth = (span_bases as f64 / genome_size).min(1.0);
    let fixed_chunk = (fixed_bases as f64 / genome_size).min(1.0);
    let flex_chunk = (flex_bases as f64 / genome_size).min(1.0);
    let expected_breadth = expected_breadth(read_depth, assembly.count as f64);
    // Deliberately an assembly-level expectation: it assumes a full window per
    // mapping, while observed chunks are clamped to the contig each mapping
    // landed on. The mismatch is load-bearing rather than an approximation to
    // tidy up -- coverage confined to short contigs then yields a small
    // observed chunk against a genome-wide expectation, and that collapsing
    // ratio is how contamination and cross-mapping are rejected. Deriving the
    // bin count per contig instead would raise COEBR for exactly those cases,
    // so it needs joint recalibration of the ratio floors, not a local fix.
    let flex_bins = (genome_size.sqrt() / 2.0).floor().max(1.0);
    let expected_flex_chunk = 1.0 - (1.0 - 1.0 / flex_bins).powf(assembly.count as f64);
    let oebr = ratio_or_nan(breadth, expected_breadth);
    let coebr = ratio_or_nan(flex_chunk, expected_flex_chunk);
    let tests = multiple_testing_taxa.max(1) as f64;
    let breadth_pvalue = adjusted_deviation_pvalue(
        count,
        Some(read_depth),
        breadth,
        expected_breadth,
        oebr,
        tests,
    );
    let chunk_pvalue =
        adjusted_deviation_pvalue(count, None, flex_chunk, expected_flex_chunk, coebr, tests);
    EvidenceCoverage {
        count,
        depth,
        breadth,
        expected_breadth,
        fixed_chunk,
        flex_chunk,
        expected_flex_chunk,
        oebr,
        coebr,
        breadth_pvalue,
        chunk_pvalue,
        marker_genes_hit: None,
        window_spread: None,
    }
}

/// Bases the fragments' reads cover, merged per contig, and the bases that
/// second mates add beyond their first mate, summed over fragments. Without
/// pairs this is the merged span alone.
fn read_bases(intervals: &mut [PresenceInterval]) -> (u64, u64) {
    let has_mate = |record: &PresenceInterval| record.mate_span.1 > record.mate_span.0;
    if !intervals.iter().any(has_mate) {
        return (merged_interval_bases(intervals, |record| record.span), 0);
    }
    let mut reads = Vec::with_capacity(2 * intervals.len());
    let mut second_mate_bases = 0u64;
    for record in intervals.iter() {
        reads.push((record.reference_id, record.span));
        if has_mate(record) {
            let (span, mate) = (record.span, record.mate_span);
            let overlap = mate.1.min(span.1).saturating_sub(mate.0.max(span.0));
            second_mate_bases += u64::from(mate.1 - mate.0 - overlap);
            reads.push((record.reference_id, mate));
        }
    }
    merge_spans(&mut reads);
    let covered = reads
        .iter()
        .map(|(_, span)| u64::from(span.1 - span.0))
        .sum();
    (covered, second_mate_bases)
}

/// Merge `spans`, each `(contig, (start, end))`, into the segments they cover:
/// spans of a contig that overlap or touch join, empty spans are left out, and
/// the segments come in contig and position order.
pub(crate) fn merge_spans(spans: &mut Vec<(u32, (u32, u32))>) {
    spans.retain(|(_, span)| span.1 > span.0);
    spans.sort_unstable();
    let mut segments = 0usize;
    for index in 0..spans.len() {
        let (contig, span) = spans[index];
        match segments.checked_sub(1).map(|last| &mut spans[last]) {
            Some((last_contig, last)) if *last_contig == contig && span.0 <= last.1 => {
                last.1 = last.1.max(span.1);
            }
            _ => {
                spans[segments] = (contig, span);
                segments += 1;
            }
        }
    }
    spans.truncate(segments);
}

pub(crate) fn merged_interval_bases<F>(intervals: &mut [PresenceInterval], select: F) -> u64
where
    F: Fn(&PresenceInterval) -> (u32, u32) + Copy,
{
    if intervals.is_empty() {
        return 0;
    }
    intervals.sort_unstable_by_key(|record| {
        let interval = select(record);
        (record.reference_id, interval.0, interval.1)
    });
    let mut total = 0u64;
    let mut reference = intervals[0].reference_id;
    let mut current = select(&intervals[0]);
    for record in &intervals[1..] {
        let interval = select(record);
        if record.reference_id == reference && interval.0 <= current.1 {
            current.1 = current.1.max(interval.1);
        } else {
            total = total.saturating_add(u64::from(current.1.saturating_sub(current.0)));
            reference = record.reference_id;
            current = interval;
        }
    }
    total.saturating_add(u64::from(current.1.saturating_sub(current.0)))
}

/// Expected fraction of a reference covered by `count` reads under the
/// Lander-Waterman uniform-placement occupancy model, given mean coverage
/// `depth` (mapped bases / reference length).
#[inline]
pub(crate) fn expected_breadth(depth: f64, count: f64) -> f64 {
    if count <= f64::EPSILON {
        0.0
    } else {
        let mean_span_fraction = (depth / count).clamp(0.0, 1.0);
        1.0 - (1.0 - mean_span_fraction).powf(count)
    }
}

#[inline]
fn ratio_or_nan(observed: f64, expected: f64) -> f64 {
    if expected > 0.0 {
        observed / expected
    } else {
        f64::NAN
    }
}

pub(crate) fn adjusted_deviation_pvalue(
    count: usize,
    depth: Option<f64>,
    observed: f64,
    expected: f64,
    ratio: f64,
    tests: f64,
) -> Option<f64> {
    if !ratio.is_finite() {
        return None;
    }
    if NEUTRAL_RATIO.contains(&ratio) {
        return Some(1.0);
    }
    let side = if observed < expected {
        Tail::Lower
    } else {
        Tail::Upper
    };
    calc_breadth_pvalue(count, depth, observed, Some(expected), side)
        .map(|value| (value * tests).min(1.0))
}

/// Reusable per-worker candidate storage. Hash-table capacities survive from
/// one fragment to the next, avoiding allocator traffic in the hot loop.
#[derive(Default)]
pub(crate) struct CandidateScratch {
    by_taxon_assembly: HashMap<(u32, u32), CandidatePlacement>,
    /// Best long-read placement per taxon, assembly and marker gene.
    by_taxon_gene: HashMap<(u32, u32, u32), CandidatePlacement>,
    other_gene_hits: Vec<((u32, u32), PresenceInterval)>,
    by_taxon: HashMap<u32, TerminalCandidate>,
    paired_rows: HashMap<u32, ExplicitPairRows>,
    pair_order: Vec<u32>,
    half_pair_mates: Vec<usize>,
    output: Vec<TerminalCandidate>,
}

#[derive(Clone, Copy, Debug, Default)]
struct ExplicitPairRows {
    mate1: Option<usize>,
    mate2: Option<usize>,
    malformed: bool,
}

impl CandidateScratch {
    pub(crate) fn build<'a>(
        &'a mut self,
        rows: &[EvidenceAlignment],
        policy: EvidencePolicy,
    ) -> &'a [TerminalCandidate] {
        self.by_taxon_assembly.clear();
        self.by_taxon_gene.clear();
        self.other_gene_hits.clear();
        self.by_taxon.clear();
        self.paired_rows.clear();
        self.pair_order.clear();
        self.half_pair_mates.clear();
        self.output.clear();

        match policy.mode {
            EvidenceMode::Paired => self.build_paired(rows, policy),
            EvidenceMode::Single => self.build_single(rows, policy),
            EvidenceMode::Long => self.build_long(rows, policy),
        }

        for (_, placement) in self.by_taxon_assembly.drain() {
            match self.by_taxon.entry(placement.taxid) {
                Entry::Vacant(entry) => {
                    entry.insert(TerminalCandidate {
                        taxid: placement.taxid,
                        best: placement,
                        placements: vec![placement],
                        other_gene_hits: Vec::new(),
                    });
                }
                Entry::Occupied(mut entry) => {
                    let taxon = entry.get_mut();
                    if placement_is_better(&placement, &taxon.best, policy.mode) {
                        taxon.best = placement;
                    }
                    taxon.placements.push(placement);
                }
            }
        }
        for ((taxid, assembly_id), interval) in self.other_gene_hits.drain(..) {
            if let Some(taxon) = self.by_taxon.get_mut(&taxid) {
                taxon.other_gene_hits.push(GeneHit {
                    assembly_id,
                    interval,
                });
            }
        }
        self.output
            .extend(self.by_taxon.drain().map(|(_, value)| value));
        for candidate in &mut self.output {
            candidate
                .placements
                .sort_unstable_by_key(|placement| placement.assembly_id);
        }
        self.output
            .sort_by(|left, right| candidate_order(left, right, policy.mode));
        &self.output
    }

    fn build_paired(&mut self, rows: &[EvidenceAlignment], policy: EvidencePolicy) {
        for (index, row) in rows.iter().enumerate() {
            if row.pair_id == 0 {
                continue;
            }
            let entry = match self.paired_rows.entry(row.pair_id) {
                Entry::Vacant(entry) => {
                    self.pair_order.push(row.pair_id);
                    entry.insert(ExplicitPairRows::default())
                }
                Entry::Occupied(entry) => entry.into_mut(),
            };
            let slot = match row.role {
                MappingRole::Mate1 => &mut entry.mate1,
                MappingRole::Mate2 => &mut entry.mate2,
                MappingRole::Unpaired => {
                    entry.malformed = true;
                    continue;
                }
            };
            if slot.replace(index).is_some() {
                entry.malformed = true;
            }
        }

        for order_index in 0..self.pair_order.len() {
            let pair_id = self.pair_order[order_index];
            let Some(pair) = self.paired_rows.get(&pair_id).copied() else {
                continue;
            };
            let (Some(r1_index), Some(r2_index)) = (pair.mate1, pair.mate2) else {
                continue;
            };
            if pair.malformed {
                continue;
            }
            let r1 = &rows[r1_index];
            let r2 = &rows[r2_index];
            if r1.reference_id != r2.reference_id
                || r1.profile_taxid != r2.profile_taxid
                || r1.assembly_id != r2.assembly_id
            {
                continue;
            }
            if !pair_qualifies(r1, r2, policy) {
                // A mate that runs off a marker gene is clipped there and can
                // fail its gates; its partner may still stand for the fragment.
                let (r1_alone, r2_alone) = (
                    single_mapping_qualifies(r1, policy),
                    single_mapping_qualifies(r2, policy),
                );
                if !r2_alone && r2.reaches_marker_end() && r1_alone {
                    self.half_pair_mates.push(r1_index);
                } else if !r1_alone && r1.reaches_marker_end() && r2_alone {
                    self.half_pair_mates.push(r2_index);
                }
                continue;
            }
            let placement = CandidatePlacement {
                taxid: r1.profile_taxid,
                reference_id: r1.reference_id,
                assembly_id: r1.assembly_id,
                sampling_scale: f64::from(r1.sampling_scale),
                asm_len: r1.asm_len,
                r1_mapped: r1.mapped_len,
                both_mapped: r1.mapped_len.saturating_add(r2.mapped_len),
                read_len: r1.read_len,
                other_genes: 0,
                pair_errors: estimated_errors(r1) + estimated_errors(r2),
                long_score: 0.0,
                span: r1.span,
                mate_span: r2.span,
                fixed_chunk: r1.fixed_chunk,
                flex_chunk: r1.flex_chunk,
                input_order: r1_index.min(r2_index),
            };
            self.insert_placement(placement, policy.mode);
        }

        // Without any pair, the fragment's evidence is a mate whose mate may
        // lie beyond the end of its reference, or a mate whose partner on the
        // same marker gene was clipped at the gene's end: the only single
        // mates that paired-end profiling admits.
        if self.by_taxon_assembly.is_empty() {
            let half_pair_mates = std::mem::take(&mut self.half_pair_mates);
            self.build_unpaired(rows, policy, |index, mapping| {
                mapping.pair_id == 0 || half_pair_mates.contains(&index)
            });
            self.half_pair_mates = half_pair_mates;
        }
    }

    fn build_single(&mut self, rows: &[EvidenceAlignment], policy: EvidencePolicy) {
        self.build_unpaired(rows, policy, |_, _| true);
    }

    /// One placement per qualifying single-read row accepted by `include`.
    fn build_unpaired(
        &mut self,
        rows: &[EvidenceAlignment],
        policy: EvidencePolicy,
        include: impl Fn(usize, &EvidenceAlignment) -> bool,
    ) {
        for (input_order, mapping) in rows.iter().enumerate() {
            if !include(input_order, mapping) || !single_mapping_qualifies(mapping, policy) {
                continue;
            }
            let placement = CandidatePlacement {
                taxid: mapping.profile_taxid,
                reference_id: mapping.reference_id,
                assembly_id: mapping.assembly_id,
                sampling_scale: f64::from(mapping.sampling_scale),
                asm_len: mapping.asm_len,
                r1_mapped: mapping.mapped_len,
                both_mapped: mapping.mapped_len,
                read_len: mapping.read_len,
                other_genes: 0,
                pair_errors: estimated_errors(mapping),
                long_score: 0.0,
                span: mapping.span,
                mate_span: (0, 0),
                fixed_chunk: mapping.fixed_chunk,
                flex_chunk: mapping.flex_chunk,
                input_order,
            };
            self.insert_row_placement(placement, mapping, policy.mode);
        }
        self.count_other_genes();
    }

    fn build_long(&mut self, rows: &[EvidenceAlignment], policy: EvidencePolicy) {
        for (input_order, mapping) in rows.iter().enumerate() {
            if !long_mapping_qualifies(mapping, policy) {
                continue;
            }
            let placement = CandidatePlacement {
                taxid: mapping.profile_taxid,
                reference_id: mapping.reference_id,
                assembly_id: mapping.assembly_id,
                sampling_scale: f64::from(mapping.sampling_scale),
                asm_len: mapping.asm_len,
                r1_mapped: mapping.mapped_len,
                both_mapped: mapping.mapped_len,
                read_len: mapping.read_len,
                other_genes: 0,
                pair_errors: estimated_errors(mapping),
                long_score: long_score(mapping),
                span: mapping.span,
                mate_span: (0, 0),
                fixed_chunk: mapping.fixed_chunk,
                flex_chunk: mapping.flex_chunk,
                input_order,
            };
            self.insert_row_placement(placement, mapping, policy.mode);
        }
        self.count_other_genes();
    }

    /// [`Self::insert_placement`] for the placement of one row, which is also
    /// kept as the best on its marker gene when the row lies on one.
    fn insert_row_placement(
        &mut self,
        placement: CandidatePlacement,
        mapping: &EvidenceAlignment,
        mode: EvidenceMode,
    ) {
        if mapping.marker_gene_len > 0 {
            let gene = (
                placement.taxid,
                placement.assembly_id,
                placement.reference_id,
            );
            match self.by_taxon_gene.entry(gene) {
                Entry::Vacant(entry) => {
                    entry.insert(placement);
                }
                Entry::Occupied(mut entry) => {
                    if placement_is_better(&placement, entry.get(), mode) {
                        entry.insert(placement);
                    }
                }
            }
        }
        self.insert_placement(placement, mode);
    }

    /// A fragment can lie on several marker genes of one genome: a long read
    /// that runs over them, or the two mates of a pair on two neighbouring
    /// genes. Its placement keeps the best gene and counts the others, whose
    /// coverage the fragment adds as well.
    fn count_other_genes(&mut self) {
        for ((taxid, assembly_id, reference_id), gene) in self.by_taxon_gene.drain() {
            let key = (taxid, assembly_id);
            let placement = self
                .by_taxon_assembly
                .get_mut(&key)
                .expect("a placed gene belongs to a placement");
            if reference_id != placement.reference_id {
                placement.other_genes += 1;
                self.other_gene_hits.push((
                    key,
                    PresenceInterval {
                        reference_id,
                        span: gene.span,
                        mate_span: (0, 0),
                        fixed_chunk: gene.fixed_chunk,
                        flex_chunk: gene.flex_chunk,
                    },
                ));
            }
        }
        self.other_gene_hits
            .sort_unstable_by_key(|&(key, interval)| (key, interval.reference_id));
    }

    fn insert_placement(&mut self, placement: CandidatePlacement, mode: EvidenceMode) {
        match self
            .by_taxon_assembly
            .entry((placement.taxid, placement.assembly_id))
        {
            Entry::Vacant(entry) => {
                entry.insert(placement);
            }
            Entry::Occupied(mut entry) => {
                if placement_is_better(&placement, entry.get(), mode) {
                    entry.insert(placement);
                }
            }
        }
    }
}

/// A pair is one fragment: each mate passes the length and fraction gates,
/// and the identity over the bases of both mates passes the identity gate.
#[inline]
fn pair_qualifies(r1: &EvidenceAlignment, r2: &EvidenceAlignment, policy: EvidencePolicy) -> bool {
    mate_extent_qualifies(r1, policy)
        && mate_extent_qualifies(r2, policy)
        && pair_identity(r1, r2) >= policy.identity
}

/// Identity of a pair: one minus the mismatches of both mates over their
/// mapped bases, i.e. the mapped-length-weighted mean of the two identities,
/// canonicalized like a single mapping's identity.
#[inline]
fn pair_identity(r1: &EvidenceAlignment, r2: &EvidenceAlignment) -> f64 {
    let (len1, len2) = (f64::from(r1.mapped_len), f64::from(r2.mapped_len));
    let identity =
        (evidence_identity(r1) * len1 + evidence_identity(r2) * len2) / (len1 + len2).max(1.0);
    canonical_mapping_identity(identity)
}

#[inline]
fn mate_extent_qualifies(mapping: &EvidenceAlignment, policy: EvidencePolicy) -> bool {
    mapping.mapped_len as usize >= policy.mapped_len
        && fraction_at_least(mapping.fraction, policy.fraction)
        && mapping.span.1 > mapping.span.0
}

#[inline]
fn single_mapping_qualifies(mapping: &EvidenceAlignment, policy: EvidencePolicy) -> bool {
    evidence_identity(mapping) >= policy.identity
        && mapping.mapped_len as usize >= policy.mapped_len
        && fraction_at_least(mapping.fraction, policy.fraction)
        && mapping.span.1 > mapping.span.0
}

#[inline]
/// Long reads pass with the mapped length or the read fraction, and on a
/// marker gene also with more than half of the gene, so that a read longer
/// than a short gene can count on it.
fn long_mapping_qualifies(mapping: &EvidenceAlignment, policy: EvidencePolicy) -> bool {
    evidence_identity(mapping) >= policy.identity
        && (mapping.mapped_len as usize >= policy.mapped_len
            || fraction_at_least(mapping.fraction, policy.fraction)
            || mapping.aligns_over_half_of_marker_gene())
        && mapping.span.1 > mapping.span.0
}

#[inline]
fn estimated_errors(mapping: &EvidenceAlignment) -> f64 {
    (1.0 - evidence_identity(mapping)).max(0.0) * f64::from(mapping.mapped_len)
}

#[inline]
fn long_score(mapping: &EvidenceAlignment) -> f64 {
    f64::from(mapping.mapped_len) * (2.0 * evidence_identity(mapping) - 1.0)
}

#[inline]
fn evidence_identity(mapping: &EvidenceAlignment) -> f64 {
    // Mapping identities are canonical five-decimal values at the mapper /
    // profiler boundary. The compact in-memory record stores them as f32;
    // canonicalizing after widening recovers the exact boundary value and
    // prevents f32 representation noise from splitting an exact score tie.
    canonical_mapping_identity(f64::from(mapping.identity))
}

#[inline]
fn fraction_at_least(fraction: f32, threshold: f64) -> bool {
    // An exact threshold such as 105/150 = 0.70 may round just below 0.70 in
    // f32. One f32 epsilon restores inclusive-threshold semantics without
    // materially relaxing the biological cutoff.
    f64::from(fraction) + FRACTION_F32_EPSILON >= threshold
}

#[inline]
fn placement_is_better(
    left: &CandidatePlacement,
    right: &CandidatePlacement,
    mode: EvidenceMode,
) -> bool {
    match mode {
        EvidenceMode::Paired | EvidenceMode::Single => {
            // Short-read placements are independently ranked; use the stable
            // reference ID only after all biological evidence ties.
            left.pair_errors < right.pair_errors - SCORE_EPSILON
                || ((left.pair_errors - right.pair_errors).abs() <= SCORE_EPSILON
                    && (left.both_mapped > right.both_mapped
                        || (left.both_mapped == right.both_mapped
                            && left.reference_id < right.reference_id)))
        }
        EvidenceMode::Long => {
            // Preserve the mapper's first-ranked long-read placement on an
            // exact semantic tie. A numeric reference ID is metadata, not
            // biological evidence; using it here can move an otherwise
            // identical candidate to a different contig and manufacture a
            // different coverage pattern.
            left.long_score > right.long_score + SCORE_EPSILON
                || ((left.long_score - right.long_score).abs() <= SCORE_EPSILON
                    && (left.pair_errors < right.pair_errors - SCORE_EPSILON
                        || ((left.pair_errors - right.pair_errors).abs() <= SCORE_EPSILON
                            && (left.both_mapped > right.both_mapped
                                || (left.both_mapped == right.both_mapped
                                    && left.input_order < right.input_order)))))
        }
    }
}

fn candidate_order(
    left: &TerminalCandidate,
    right: &TerminalCandidate,
    mode: EvidenceMode,
) -> Ordering {
    match mode {
        EvidenceMode::Paired | EvidenceMode::Single => left
            .best
            .pair_errors
            .partial_cmp(&right.best.pair_errors)
            .unwrap_or(Ordering::Equal)
            .then_with(|| left.taxid.cmp(&right.taxid)),
        EvidenceMode::Long => right
            .best
            .long_score
            .partial_cmp(&left.best.long_score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| left.best.input_order.cmp(&right.best.input_order))
            .then_with(|| left.taxid.cmp(&right.taxid)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mapping(
        species: u32,
        reference: u32,
        assembly: u32,
        identity: f32,
        mapped: u32,
        read_len: u32,
        role: MappingRole,
    ) -> EvidenceAlignment {
        EvidenceAlignment {
            profile_taxid: species,
            reference_id: reference,
            assembly_id: assembly,
            role,
            pair_id: if matches!(role, MappingRole::Mate1 | MappingRole::Mate2) {
                reference.max(1)
            } else {
                0
            },
            sampling_scale: 1.0,
            asm_len: 1_000_000.0,
            identity,
            fraction: mapped as f32 / read_len as f32,
            mapped_len: mapped,
            read_len,
            span: (100, 100 + mapped),
            fixed_chunk: (1, 1_100 + mapped),
            flex_chunk: (1, 2_100 + mapped),
            marker_gene_len: 0,
        }
    }

    fn policy(mode: EvidenceMode) -> EvidencePolicy {
        EvidencePolicy {
            mode,
            identity: 0.97,
            mapped_len: 50,
            fraction: 0.70,
        }
    }

    #[test]
    fn paired_candidates_use_explicit_pair_ids_not_row_adjacency() {
        let mut rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        rows[0].pair_id = 101;
        rows[2].pair_id = 101;
        rows[1].pair_id = 202;
        rows[3].pair_id = 202;

        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 2);
        assert_eq!(candidates[0].best.both_mapped, 200);
        assert_eq!(candidates[1].best.both_mapped, 200);
    }

    #[test]
    fn lone_mates_place_a_fragment_only_without_a_pair() {
        // Mates of one fragment beyond the ends of two genes of species 1, and
        // a lone mate of species 2 that fails the identity gate.
        let mut lone = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 11, 1, 0.98, 100, 100, MappingRole::Mate2),
            mapping(2, 20, 2, 0.80, 100, 100, MappingRole::Mate2),
        ];
        for row in &mut lone {
            row.pair_id = 0;
        }
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&lone, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].taxid, 1);
        assert_eq!(candidates[0].placements.len(), 1);
        assert_eq!(candidates[0].best.reference_id, 10);
        assert_eq!(candidates[0].best.both_mapped, 100);

        let mut rows = vec![
            mapping(3, 30, 3, 0.99, 100, 100, MappingRole::Mate1),
            mapping(3, 30, 3, 0.99, 100, 100, MappingRole::Mate2),
        ];
        rows[0].pair_id = 7;
        rows[1].pair_id = 7;
        rows.extend(lone);
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].taxid, 3);
    }

    #[test]
    fn lone_mates_on_two_marker_genes_count_on_both() {
        // The two mates of one fragment on two marker genes of genome 1.
        let mut rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 11, 1, 0.98, 100, 100, MappingRole::Mate2),
        ];
        for row in &mut rows {
            row.pair_id = 0;
            row.marker_gene_len = 1_000;
        }
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].best.reference_id, 10);
        assert_eq!(candidates[0].best.other_genes, 1);
        let others: Vec<_> = candidates[0]
            .other_gene_hits
            .iter()
            .map(|hit| (hit.assembly_id, hit.interval.reference_id))
            .collect();
        assert_eq!(others, [(1, 11)]);

        // On genome references the fragment keeps one placement.
        for row in &mut rows {
            row.marker_gene_len = 0;
        }
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates[0].best.other_genes, 0);
        assert!(candidates[0].other_gene_hits.is_empty());
    }

    #[test]
    fn a_long_read_counts_on_more_than_half_of_a_marker_gene() {
        let long_policy = EvidencePolicy {
            identity: 0.9,
            mapped_len: 500,
            fraction: 0.2,
            ..policy(EvidenceMode::Long)
        };
        // 250 bases of a 3,000-base read: neither 500 bases nor 20 % of it.
        let mut row = mapping(1, 10, 1, 0.95, 250, 3_000, MappingRole::Unpaired);
        assert!(!long_mapping_qualifies(&row, long_policy));
        // More than half of a 480-base gene.
        row.marker_gene_len = 480;
        assert!(long_mapping_qualifies(&row, long_policy));
        // Half of a 500-base gene is not enough.
        row.marker_gene_len = 500;
        assert!(!long_mapping_qualifies(&row, long_policy));
        // The identity cutoff still applies.
        row.marker_gene_len = 480;
        row.identity = 0.85;
        assert!(!long_mapping_qualifies(&row, long_policy));
    }

    #[test]
    fn a_long_read_counts_every_marker_gene_of_a_genome_it_covers() {
        let long_policy = EvidencePolicy {
            identity: 0.9,
            mapped_len: 500,
            fraction: 0.2,
            ..policy(EvidenceMode::Long)
        };
        // One read over genes 10 and 11 of genome 1 and gene 12 of genome 2.
        let mut rows = vec![
            mapping(1, 10, 1, 0.99, 800, 3_000, MappingRole::Unpaired),
            mapping(1, 11, 1, 0.98, 600, 3_000, MappingRole::Unpaired),
            mapping(1, 12, 2, 0.97, 700, 3_000, MappingRole::Unpaired),
        ];
        for row in &mut rows {
            row.marker_gene_len = 1_000;
        }
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, long_policy);
        assert_eq!(candidates.len(), 1);
        let placements: Vec<_> = candidates[0]
            .placements
            .iter()
            .map(|placement| {
                (
                    placement.assembly_id,
                    placement.reference_id,
                    placement.other_genes,
                )
            })
            .collect();
        assert_eq!(placements, [(1, 10, 1), (2, 12, 0)]);
        let others: Vec<_> = candidates[0]
            .other_gene_hits
            .iter()
            .map(|hit| {
                (
                    hit.assembly_id,
                    hit.interval.reference_id,
                    hit.interval.span,
                )
            })
            .collect();
        assert_eq!(others, [(1, 11, (100, 700))]);

        // Genome references keep one placement per genome.
        for row in &mut rows {
            row.marker_gene_len = 0;
        }
        let candidates = scratch.build(&rows, long_policy);
        assert!(candidates[0]
            .placements
            .iter()
            .all(|placement| placement.other_genes == 0));
        assert!(candidates[0].other_gene_hits.is_empty());
    }

    #[test]
    fn an_alignment_reaches_a_marker_gene_end_at_its_first_or_last_base() {
        let mut row = mapping(1, 10, 1, 1.0, 100, 100, MappingRole::Unpaired);
        let mut reaches = |span, gene_len| {
            row.span = span;
            row.marker_gene_len = gene_len;
            row.reaches_marker_end()
        };
        assert!(reaches((1, 101), 1_000));
        assert!(reaches((900, 1_001), 1_000));
        assert!(!reaches((2, 1_000), 1_000));
        // Not a marker gene.
        assert!(!reaches((1, 101), 0));
    }

    #[test]
    fn a_mate_clipped_at_a_marker_gene_end_leaves_its_partner_as_evidence() {
        // Species 1: mate 2 runs off a marker gene and aligns 60 of 100 bases.
        let mut rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 1.00, 60, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        assert!(scratch
            .build(&rows, policy(EvidenceMode::Paired))
            .is_empty());

        rows[1].marker_gene_len = rows[1].span.1 - 1;
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].taxid, 1);
        assert_eq!(candidates[0].best.both_mapped, 100);

        // The same with the mates swapped: mate 1 is clipped.
        let mut swapped = vec![
            mapping(1, 10, 1, 1.00, 60, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
        ];
        swapped[0].marker_gene_len = swapped[0].span.1 - 1;
        let candidates = scratch.build(&swapped, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].best.both_mapped, 100);

        // The good mate must pass the single-read gates itself.
        let mut weak = rows.clone();
        weak[0].identity = 0.90;
        assert!(scratch
            .build(&weak, policy(EvidenceMode::Paired))
            .is_empty());

        // A qualifying pair of the fragment on another taxon takes precedence.
        let mut with_pair = rows.clone();
        with_pair.push(mapping(2, 20, 2, 0.98, 100, 100, MappingRole::Mate1));
        with_pair.push(mapping(2, 20, 2, 0.98, 100, 100, MappingRole::Mate2));
        let candidates = scratch.build(&with_pair, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].taxid, 2);
    }

    #[test]
    fn paired_candidate_requires_identity_over_both_mates() {
        let pair = |mate1, mate2| {
            let rows = vec![
                mapping(1, 10, 1, mate1, 100, 100, MappingRole::Mate1),
                mapping(1, 10, 1, mate2, 100, 100, MappingRole::Mate2),
            ];
            let mut scratch = CandidateScratch::default();
            scratch
                .build(&rows, policy(EvidenceMode::Paired))
                .iter()
                .map(|candidate| candidate.best.both_mapped)
                .collect::<Vec<_>>()
        };
        // 0.9725 over both mates passes, whichever mate is the weaker one.
        assert_eq!(pair(0.99, 0.955), vec![200]);
        assert_eq!(pair(0.955, 0.99), vec![200]);
        // The cutoff is inclusive; 0.965 over both mates fails.
        assert_eq!(pair(0.97, 0.97), vec![200]);
        assert!(pair(0.98, 0.95).is_empty());
    }

    #[test]
    fn compact_identity_and_fraction_preserve_canonical_inclusive_thresholds() {
        let rows = vec![
            mapping(1, 10, 1, 0.97, 70, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.97, 70, 100, MappingRole::Mate2),
        ];
        // 0.70 is represented slightly below 0.70 as f32 on common targets;
        // the inclusive profiler cutoff must nevertheless accept it.
        assert!(f64::from(rows[0].fraction) < 0.70);
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(evidence_identity(&rows[0]), 0.97);
    }

    #[test]
    fn long_exact_placement_tie_preserves_mapper_order() {
        let mut long_policy = policy(EvidenceMode::Long);
        long_policy.identity = 0.92;
        long_policy.mapped_len = 500;
        long_policy.fraction = 0.20;
        let rows = vec![
            mapping(1, 20, 1, 0.99, 1_000, 2_000, MappingRole::Unpaired),
            mapping(1, 10, 1, 0.99, 1_000, 2_000, MappingRole::Unpaired),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, long_policy);
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].placements.len(), 1);
        assert_eq!(candidates[0].placements[0].reference_id, 20);
    }

    #[test]
    fn long_exact_tie_across_assemblies_preserves_mapper_order_deterministically() {
        let mut long_policy = policy(EvidenceMode::Long);
        long_policy.identity = 0.92;
        long_policy.mapped_len = 500;
        long_policy.fraction = 0.20;
        let rows = vec![
            mapping(1, 20, 2, 0.99, 1_000, 2_000, MappingRole::Unpaired),
            mapping(1, 10, 1, 0.99, 1_000, 2_000, MappingRole::Unpaired),
        ];
        let mut scratch = CandidateScratch::default();
        for _ in 0..100 {
            let candidates = scratch.build(&rows, long_policy);
            assert_eq!(candidates.len(), 1);
            assert_eq!(candidates[0].placements.len(), 2);
            assert_eq!(candidates[0].best.reference_id, 20);
            assert_eq!(candidates[0].best.assembly_id, 2);
        }
    }

    #[test]
    fn paired_exact_placement_tie_uses_stable_reference_id() {
        let rows = vec![
            mapping(1, 20, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 20, 1, 0.99, 100, 100, MappingRole::Mate2),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].placements.len(), 1);
        assert_eq!(candidates[0].placements[0].reference_id, 10);
    }

    #[test]
    fn placements_from_two_assemblies_remain_separate_witnesses() {
        let rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            mapping(1, 20, 2, 0.98, 100, 100, MappingRole::Mate1),
            mapping(1, 20, 2, 0.98, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].placements.len(), 2);
        assert_ne!(
            candidates[0].placements[0].assembly_id,
            candidates[0].placements[1].assembly_id
        );
    }

    #[test]
    fn paired_breadth_counts_both_mates_and_expects_their_bases() {
        let fragment = |span: (u32, u32), mate_span: (u32, u32)| PresenceInterval {
            reference_id: 1,
            span,
            mate_span,
            fixed_chunk: span,
            flex_chunk: span,
        };
        let mut assembly = PresenceAssembly {
            asm_len: 10_000.0,
            count: 2,
            reference_span_bases: 200,
            intervals: vec![
                fragment((0, 100), (300, 400)),
                // Its second mate overlaps the first mate of the other fragment.
                fragment((350, 450), (600, 700)),
            ],
            ..PresenceAssembly::default()
        };
        let coverage = calculate_coverage(&mut assembly, 1);
        // (0, 100), (300, 450) and (600, 700).
        assert!((coverage.breadth - 0.035).abs() < 1e-12);
        // Each fragment covers 200 bases: 1 - (1 - 200 / 10,000)^2.
        assert!((coverage.expected_breadth - 0.0396).abs() < 1e-12);
        // The reported depth counts each fragment once, by its first read.
        assert!((coverage.depth - 0.02).abs() < 1e-12);
    }

    #[test]
    fn coverage_merges_within_contigs_but_not_across_contigs() {
        let mut assembly = PresenceAssembly {
            asm_len: 1_000.0,
            count: 3,
            other_gene_hits: 0,
            reference_span_bases: 250,
            fits: 0,
            sampling_scale: 1.0,
            intervals: vec![
                PresenceInterval {
                    reference_id: 1,
                    span: (0, 100),
                    mate_span: (0, 0),
                    fixed_chunk: (0, 200),
                    flex_chunk: (0, 300),
                },
                PresenceInterval {
                    reference_id: 1,
                    span: (50, 150),
                    mate_span: (0, 0),
                    fixed_chunk: (100, 300),
                    flex_chunk: (200, 400),
                },
                PresenceInterval {
                    reference_id: 2,
                    span: (0, 50),
                    mate_span: (0, 0),
                    fixed_chunk: (0, 100),
                    flex_chunk: (0, 100),
                },
            ],
        };
        let coverage = calculate_coverage(&mut assembly, 1);
        assert!((coverage.breadth - 0.20).abs() < 1e-12);
        assert!((coverage.fixed_chunk - 0.40).abs() < 1e-12);
        assert!((coverage.flex_chunk - 0.50).abs() < 1e-12);
        assert!((coverage.depth - 0.25).abs() < 1e-12);
    }

    #[test]
    fn spans_join_where_they_overlap_or_touch() {
        let mut spans = vec![
            (1, (10, 20)),
            (0, (5, 8)),
            (1, (20, 25)),
            (1, (30, 30)),
            (1, (15, 18)),
            (1, (40, 50)),
        ];
        merge_spans(&mut spans);
        assert_eq!(spans, [(0, (5, 8)), (1, (10, 25)), (1, (40, 50))]);
    }

    #[test]
    fn floor_gates_keep_the_length_cutoff_of_long_reads() {
        let floor = policy(EvidenceMode::Paired).floor(0.8);
        assert_eq!(
            (floor.identity, floor.mapped_len, floor.fraction),
            (0.8, FLOOR_MIN_MAPPED_LEN, FLOOR_MIN_FRACTION)
        );
        let long = EvidencePolicy {
            mapped_len: 1_000,
            ..policy(EvidenceMode::Long)
        };
        assert_eq!(long.floor(0.8).mapped_len, 1_000);
    }
}

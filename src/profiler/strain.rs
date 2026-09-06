//! Strain Profiler - Strain-level taxonomic abundance estimation
//!
//! Experimental. This module runs an independent filter chain with its own
//! presets and does not use the community evidence policy, so its calls are
//! less well characterised than community output. See the strain note in
//! MANUAL.md.
//!
//! This module implements strain-level profiling that can distinguish between
//! closely related strains within a species. It uses a similar pipeline to
//! the community profiler but operates at the strain/assembly level rather
//! than the species level.
//!
//! # Key Differences from Community Profiler
//!
//! - **Evidence policy**: `EvidencePolicy::from_config` returns `None` here,
//!   so strain output sees no standard or rescue tiers and no guarded-MAP
//!   exact-tie resolution
//! - **Resolution**: Assigns reads to individual strain rather
//!   than aggregating to species level
//! - **Taxonomy fallback**: If strain taxonomy lookup fails, falls back to
//!   species-level assignment
//! - **Coverage tracking**: Tracks coverage per assembly/strain rather than per species
//! - **Cutoffs**: Keeps its own identity, chunk-breadth and OEBR/COEBR presets
//!   rather than the community ones

use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::{
    atomic::{AtomicUsize, Ordering as AtomicOrdering},
    Arc,
};
use std::thread;

use anyhow::{anyhow, bail, Context, Result};
use crossbeam_channel::{bounded, unbounded};
use csv::WriterBuilder;
use flate2::read::MultiGzDecoder;
use indexmap::IndexMap;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use xz2::read::XzDecoder;

use metam::mapper::CandidateBoundary;

use super::alignment::cigar_to_identity;
use super::evidence::{
    EvidenceMode, LONG_RELATIVE_SCORE_MARGIN, PAIRED_ERROR_MARGIN, SINGLE_ERROR_MARGIN,
};
use super::output::{
    ClassificationSource, ProfileAnnotations, ProfileRow, ResolvedClassifications,
};
use super::pathogen::load_pathogen_table;
use super::statistics::{coverage_deviation_pvalue, round5};
use super::{
    coverage_passes, equivalent_read_count_passes_minimum, join_alignment_worker,
    log_reference_taxonomy_adjustments, log_resolved_profiler_parameters,
    parse_reference_taxonomy, parse_sam_sequence_record, push_string_keyed_value,
    record_genome_sampling_scale, record_reference_taxonomy_adjustment,
    require_alignment_reference_dictionary, sampling_scale_from_name_parts,
    sampling_scale_from_sizes, scale_automatic_chunk_breadth, strip_pair_suffix, ProfilerConfig,
    ReadSupport, ReferenceTaxonomyAdjustments, ReferenceTaxonomyResolver,
    ResolvedProfilerParameters, StableReadMap,
};
use crate::db::DatabaseMetadata;
use crate::mapper::tsv::{
    for_each_metam_mapping_record, for_each_structured_mapping_record, load_shared_metam_profile,
    open_metam_mapping_reader, total_input_units, MetamMappingReceiver, MetamMappingRecord,
    MetamParseStats, MetamProfile,
};
use crate::progress;
use crate::taxonomy::Taxonomy;

/// Size and fractional-index metadata for a single assembly.
#[derive(Clone, Debug)]
struct AssemblyStats {
    /// Total assembly length in bases
    length: usize,

    /// Full-index-equivalent support multiplier for this assembly.
    sampling_scale: f64,
}

impl Default for AssemblyStats {
    fn default() -> Self {
        Self {
            length: 0,
            sampling_scale: 1.0,
        }
    }
}

fn intern_assembly(
    name: &str,
    assembly_ids: &mut HashMap<String, u32>,
    assembly_stats: &mut Vec<AssemblyStats>,
    source: &str,
) -> Result<u32> {
    if let Some(&assembly_id) = assembly_ids.get(name) {
        return Ok(assembly_id);
    }
    let assembly_id = u32::try_from(assembly_stats.len())
        .with_context(|| format!("{source} contains more than 2^32 assemblies"))?;
    assembly_ids.insert(name.to_owned(), assembly_id);
    assembly_stats.push(AssemblyStats::default());
    Ok(assembly_id)
}

/// Coverage metrics aggregated by taxonomy ID.
///
/// For strain-level profiling, each strain taxid has its own entry.
/// Used for filtering and statistical calculations.
#[derive(Clone, Debug, Default)]
struct TaxidCoverage {
    breadth: f64,

    chunk_breadth: f64,

    flex_chunk_breadth: f64,

    num_flex_chunks: usize,

    /// Total assembly length for this taxid
    asm_length: usize,
}

/// Metadata extracted from a reference sequence name.
///
/// Reference names follow: `asm|strain_taxid|species_taxid|acc|genome_size[|sampled]`
#[derive(Clone, Debug, Default)]
struct ReferenceMetadata {
    /// Stable numeric reference identifier used by compact coverage keys.
    reference_id: u32,

    /// Stable numeric assembly identifier used by compact coverage keys.
    assembly_id: u32,

    /// Reference sequence length from SAM header
    length: usize,

    /// Primary taxonomy ID (strain-level if available)
    taxid: Option<u32>,

    /// Fallback taxonomy ID (species-level)
    fallback_taxid: Option<u32>,
}

/// Shared context for worker threads during alignment parsing.
struct WorkerContext {
    /// Reference metadata: ref_name -> metadata
    ref_metadata: Arc<HashMap<String, ReferenceMetadata>>,

    /// Whether input is paired-end data
    is_paired: bool,
}

/// Result from a worker thread processing alignment lines.
struct WorkerResult {
    /// Parsed alignments: (read_name, record)
    alignments: Vec<(String, AlignmentRecord)>,

    /// All read names encountered
    reads_seen: Vec<String>,

    /// Count of alignments processed
    alignments_processed: usize,

    /// Mapped records whose reference was missing or unusable.
    skipped_reference: usize,
}

/// Intermediate result from parsing a single SAM line.
struct ParsedAlignment {
    /// Read name (pair suffix stripped if paired-end)
    read_name: String,

    /// Whether to count this in alignment statistics
    count_alignment: bool,

    /// Parsed alignment record (None if unmapped/filtered)
    record: Option<ParsedAlignmentRecord>,

    /// The mapped reference was absent from the sequence dictionary or had
    /// unusable metadata.
    skipped_reference: bool,
}

/// Full parsed data from a valid alignment.
struct ParsedAlignmentRecord {
    /// Alignment span coordinates (start, end)
    span: (usize, usize),

    /// Number of aligned bases
    mapped_len: usize,

    /// CIGAR string
    cigar: String,

    /// Stable reference and assembly identifiers.
    reference_id: u32,
    assembly_id: u32,

    /// Primary taxonomy ID
    taxid: Option<u32>,

    /// Fallback (species) taxonomy ID
    fallback_taxid: Option<u32>,
}

/// Absorb worker results into main aggregation structures.
///
/// Integrates results from a single worker batch into the main thread's
/// data structures. Called repeatedly as workers complete batches.
fn absorb_worker_result(
    result: WorkerResult,
    read_alignments_map: &mut StableReadMap<AlignmentRecord>,
    all_reads: &mut HashSet<String>,
    num_processed: &mut usize,
    skipped_reference: &mut usize,
    next_progress: &mut usize,
) {
    for (read, aln) in result.alignments {
        read_alignments_map.entry(read).or_default().push(aln);
    }
    for read in result.reads_seen {
        all_reads.insert(read);
    }
    *num_processed += result.alignments_processed;
    *skipped_reference += result.skipped_reference;
    while *num_processed >= *next_progress {
        log::info!(
            target: "PROFILER",
            "{} alignments processed",
            *num_processed
        );
        *next_progress += 1_000_000;
    }
}

/// Aggregated results from SAM loading phase.
///
/// Contains all data needed for the profiling algorithm.
struct SamLoadResult {
    /// Per-read alignment collections
    read_alignments: Vec<ReadAlignment>,

    /// Per-assembly coverage statistics
    /// Assembly metadata indexed by `AlignmentRecord::assembly_id`.
    assembly_stats: Vec<AssemblyStats>,

    /// Per-taxid coverage (propagated from assembly stats)
    taxid_coverage: HashMap<u32, TaxidCoverage>,

    /// Minimum breadth threshold for reporting
    min_breadth: f64,

    /// Minimum chunk breadth threshold (auto-estimated or user-set)
    min_chunk_breadth: f64,

    /// Whether at least one genome has a full/effective size multiplier > 1.
    has_scaled_genomes: bool,
}

/// Simplified alignment record for taxonomy assignment phase.
#[derive(Clone, Debug)]
struct AlignmentRecord {
    /// Number of aligned bases
    mapped_len: usize,

    /// Alignment identity. For SAM input this is derived from CIGAR; for
    /// Metam input this is read directly from the Metam TSV.
    identity: f64,

    /// Aligned fraction of the read. For SAM input this is derived from CIGAR;
    /// for Metam input this is mapped_bases / read_len.
    fraction: f64,

    /// Stable reference and assembly identifiers for coordinate evidence.
    reference_id: u32,
    assembly_id: u32,

    /// Reference coordinate span, using 1-based start and exclusive end.
    span: (usize, usize),

    /// Explicit native paired-placement metadata. Legacy SAM records use
    /// `pair_id = 0` and retain their established first-mate behavior.
    mate_number: u8,
    pair_id: u32,
    /// Taxonomic candidate levels that may be incomplete.
    candidate_boundary: CandidateBoundary,

    /// Primary taxonomy ID (strain-level if available)
    taxid: Option<u32>,

    /// Fallback taxonomy ID (species-level)
    fallback_taxid: Option<u32>,
}

/// Mapping-quality thresholds used by both strain coverage and assignment.
/// Keeping this as a small Copy value avoids sharing the full profiler config
/// with SAM workers and prevents the two phases from drifting semantically.
#[derive(Clone, Copy, Debug)]
struct AssignmentQualityThresholds {
    min_identity: f64,
    min_mapped_len: usize,
    min_fraction: f64,
    coverage_or_gate: bool,
}

impl AssignmentQualityThresholds {
    fn from_config(cfg: &ProfilerConfig) -> Self {
        Self {
            min_identity: cfg.identity,
            min_mapped_len: cfg.mapped_len,
            min_fraction: cfg.fraction,
            coverage_or_gate: cfg.coverage_gate_is_or(),
        }
    }

    #[inline]
    fn passes(self, alignment: &AlignmentRecord) -> bool {
        if alignment.identity < self.min_identity {
            return false;
        }
        self.passes_coverage(alignment)
    }

    #[inline]
    fn passes_paired_mate2(self, alignment: &AlignmentRecord) -> bool {
        // Preserve the established first-mate strain policy: R1 carries the
        // taxonomic identity gate, while R2 must still provide a complete,
        // sufficiently covered concordant placement.
        self.passes_coverage(alignment)
    }

    #[inline]
    fn passes_coverage(self, alignment: &AlignmentRecord) -> bool {
        coverage_passes(
            alignment.mapped_len >= self.min_mapped_len,
            alignment.fraction >= self.min_fraction,
            self.coverage_or_gate,
        )
    }
}

const STRAIN_SCORE_EPSILON: f64 = 1e-9;

/// One assembly-specific strain candidate retained for assignment and,
/// when clearly resolved, coordinate coverage.
#[derive(Clone, Copy, Debug)]
struct StrainCandidatePlacement {
    taxid: u32,
    fallback_taxid: Option<u32>,
    assembly_id: u32,
    reference_id: u32,
    representative_index: usize,
    errors: f64,
    long_score: f64,
    mapped_total: usize,
}

#[derive(Clone, Copy, Debug)]
struct StrainCandidate {
    taxid: u32,
    best: StrainCandidatePlacement,
}

#[derive(Clone, Copy, Debug, Default)]
struct StrainPairRows {
    mate1: Option<usize>,
    mate2: Option<usize>,
    malformed: bool,
}

/// Allocation-reusing candidate builder shared by strain assignment and
/// strain coverage. Keeping one implementation is important: abundance may
/// retain all qualified taxa, while breadth uses only a clearly resolved
/// candidate, but both views must start from exactly the same mappings.
#[derive(Default)]
struct StrainCandidateScratch {
    by_taxon_assembly: HashMap<(u32, u32), StrainCandidatePlacement>,
    by_taxon: HashMap<u32, StrainCandidatePlacement>,
    paired_rows: HashMap<u32, StrainPairRows>,
    pair_order: Vec<u32>,
    output: Vec<StrainCandidate>,
    passing_mappings: usize,
    terminal_boundary_incomplete: bool,
}

impl StrainCandidateScratch {
    fn build(
        &mut self,
        rows: &[AlignmentRecord],
        mode: EvidenceMode,
        quality: AssignmentQualityThresholds,
    ) -> &[StrainCandidate] {
        self.by_taxon_assembly.clear();
        self.by_taxon.clear();
        self.paired_rows.clear();
        self.pair_order.clear();
        self.output.clear();
        self.passing_mappings = 0;
        self.terminal_boundary_incomplete = rows
            .iter()
            .any(|mapping| mapping.candidate_boundary.terminal_incomplete());

        if mode == EvidenceMode::Paired && rows.iter().any(|row| row.pair_id > 0) {
            self.build_explicit_pairs(rows, quality);
        } else {
            self.build_individual(rows, mode, quality);
        }

        for placement in self.by_taxon_assembly.values().copied() {
            match self.by_taxon.entry(placement.taxid) {
                Entry::Vacant(entry) => {
                    entry.insert(placement);
                }
                Entry::Occupied(mut entry) => {
                    if strain_placement_is_better(&placement, entry.get(), mode) {
                        entry.insert(placement);
                    }
                }
            }
        }
        self.output
            .extend(self.by_taxon.values().copied().map(|best| StrainCandidate {
                taxid: best.taxid,
                best,
            }));
        self.output
            .sort_by(|left, right| strain_candidate_order(left, right, mode));
        &self.output
    }

    fn build_explicit_pairs(
        &mut self,
        rows: &[AlignmentRecord],
        quality: AssignmentQualityThresholds,
    ) {
        for (index, row) in rows.iter().enumerate() {
            if row.pair_id == 0 {
                continue;
            }
            let entry = match self.paired_rows.entry(row.pair_id) {
                Entry::Vacant(entry) => {
                    self.pair_order.push(row.pair_id);
                    entry.insert(StrainPairRows::default())
                }
                Entry::Occupied(entry) => entry.into_mut(),
            };
            let slot = match row.mate_number {
                1 => &mut entry.mate1,
                2 => &mut entry.mate2,
                _ => {
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
            let (Some(mate1_index), Some(mate2_index)) = (pair.mate1, pair.mate2) else {
                continue;
            };
            if pair.malformed {
                continue;
            }
            let mate1 = &rows[mate1_index];
            let mate2 = &rows[mate2_index];
            if mate1.reference_id != mate2.reference_id
                || mate1.assembly_id != mate2.assembly_id
                || mate1.taxid != mate2.taxid
                || mate1.fallback_taxid != mate2.fallback_taxid
                || !quality.passes(mate1)
                || !quality.passes_paired_mate2(mate2)
            {
                continue;
            }
            let Some(taxid) = mate1.taxid.or(mate1.fallback_taxid) else {
                continue;
            };
            self.passing_mappings = self.passing_mappings.saturating_add(2);
            self.insert_placement(
                StrainCandidatePlacement {
                    taxid,
                    fallback_taxid: mate1.fallback_taxid,
                    assembly_id: mate1.assembly_id,
                    reference_id: mate1.reference_id,
                    representative_index: mate1_index,
                    errors: estimated_strain_errors(mate1) + estimated_strain_errors(mate2),
                    long_score: 0.0,
                    mapped_total: mate1.mapped_len.saturating_add(mate2.mapped_len),
                },
                EvidenceMode::Paired,
            );
        }
    }

    fn build_individual(
        &mut self,
        rows: &[AlignmentRecord],
        mode: EvidenceMode,
        quality: AssignmentQualityThresholds,
    ) {
        for (index, row) in rows.iter().enumerate() {
            if !quality.passes(row) {
                continue;
            }
            let Some(taxid) = row.taxid.or(row.fallback_taxid) else {
                continue;
            };
            self.passing_mappings = self.passing_mappings.saturating_add(1);
            self.insert_placement(
                StrainCandidatePlacement {
                    taxid,
                    fallback_taxid: row.fallback_taxid,
                    assembly_id: row.assembly_id,
                    reference_id: row.reference_id,
                    representative_index: index,
                    errors: estimated_strain_errors(row),
                    long_score: strain_long_score(row),
                    mapped_total: row.mapped_len,
                },
                mode,
            );
        }
    }

    fn insert_placement(&mut self, placement: StrainCandidatePlacement, mode: EvidenceMode) {
        match self
            .by_taxon_assembly
            .entry((placement.taxid, placement.assembly_id))
        {
            Entry::Vacant(entry) => {
                entry.insert(placement);
            }
            Entry::Occupied(mut entry) => {
                if strain_placement_is_better(&placement, entry.get(), mode) {
                    entry.insert(placement);
                }
            }
        }
    }

    fn resolved_placements(
        &self,
        mode: EvidenceMode,
    ) -> impl Iterator<Item = StrainCandidatePlacement> + '_ {
        let resolved_taxid = if self.terminal_boundary_incomplete {
            None
        } else {
            resolved_strain_candidate_index(&self.output, mode)
                .map(|index| self.output[index].taxid)
        };
        self.by_taxon_assembly
            .values()
            .copied()
            .filter(move |placement| Some(placement.taxid) == resolved_taxid)
    }
}

#[inline]
fn estimated_strain_errors(mapping: &AlignmentRecord) -> f64 {
    (1.0 - mapping.identity).max(0.0) * mapping.mapped_len as f64
}

#[inline]
fn strain_long_score(mapping: &AlignmentRecord) -> f64 {
    mapping.mapped_len as f64 * (2.0 * mapping.identity - 1.0)
}

#[inline]
fn strain_placement_is_better(
    left: &StrainCandidatePlacement,
    right: &StrainCandidatePlacement,
    mode: EvidenceMode,
) -> bool {
    match mode {
        EvidenceMode::Paired | EvidenceMode::Single => {
            left.errors < right.errors - STRAIN_SCORE_EPSILON
                || ((left.errors - right.errors).abs() <= STRAIN_SCORE_EPSILON
                    && (left.mapped_total > right.mapped_total
                        || (left.mapped_total == right.mapped_total
                            && left.reference_id < right.reference_id)))
        }
        EvidenceMode::Long => {
            left.long_score > right.long_score + STRAIN_SCORE_EPSILON
                || ((left.long_score - right.long_score).abs() <= STRAIN_SCORE_EPSILON
                    && (left.errors < right.errors - STRAIN_SCORE_EPSILON
                        || ((left.errors - right.errors).abs() <= STRAIN_SCORE_EPSILON
                            && (left.mapped_total > right.mapped_total
                                || (left.mapped_total == right.mapped_total
                                    && left.reference_id < right.reference_id)))))
        }
    }
}

fn strain_candidate_order(
    left: &StrainCandidate,
    right: &StrainCandidate,
    mode: EvidenceMode,
) -> Ordering {
    match mode {
        EvidenceMode::Paired | EvidenceMode::Single => left
            .best
            .errors
            .partial_cmp(&right.best.errors)
            .unwrap_or(Ordering::Equal)
            .then_with(|| right.best.mapped_total.cmp(&left.best.mapped_total))
            .then_with(|| left.taxid.cmp(&right.taxid)),
        EvidenceMode::Long => right
            .best
            .long_score
            .partial_cmp(&left.best.long_score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| {
                left.best
                    .errors
                    .partial_cmp(&right.best.errors)
                    .unwrap_or(Ordering::Equal)
            })
            .then_with(|| left.taxid.cmp(&right.taxid)),
    }
}

fn resolved_strain_candidate_index(
    candidates: &[StrainCandidate],
    mode: EvidenceMode,
) -> Option<usize> {
    if candidates.is_empty() {
        return None;
    }
    if candidates.len() == 1 {
        return Some(0);
    }
    let resolved = match mode {
        EvidenceMode::Paired => {
            candidates[1].best.errors - candidates[0].best.errors >= PAIRED_ERROR_MARGIN
        }
        EvidenceMode::Single => {
            candidates[1].best.errors - candidates[0].best.errors >= SINGLE_ERROR_MARGIN
        }
        EvidenceMode::Long => {
            let best = candidates[0].best.long_score;
            best > f64::EPSILON
                && (best - candidates[1].best.long_score) / best >= LONG_RELATIVE_SCORE_MARGIN
        }
    };
    resolved.then_some(0)
}

/// Collection of alignments for a single read.
#[derive(Clone, Debug)]
struct ReadAlignment {
    /// Read identifier
    read_name: Arc<str>,

    /// All alignments for this read
    alignments: Vec<AlignmentRecord>,
}

#[derive(Clone, Copy, Debug, Hash, Ord, PartialEq, Eq, PartialOrd)]
struct StrainCoverageReferenceKey {
    taxid: u32,
    assembly_id: u32,
    reference_id: u32,
}

#[derive(Debug, Default)]
struct ResolvedStrainIntervals {
    by_reference: HashMap<StrainCoverageReferenceKey, Vec<(usize, usize)>>,
    candidate_reads: usize,
    resolved_reads: usize,
    ambiguous_reads: usize,
    boundary_ambiguous_reads: usize,
}

impl ResolvedStrainIntervals {
    fn record(&mut self, read: &ReadAlignment, placement: StrainCandidatePlacement) {
        let Some(alignment) = read.alignments.get(placement.representative_index) else {
            return;
        };
        self.record_taxid(placement.taxid, alignment);
        if let Some(fallback) = placement.fallback_taxid {
            if fallback != placement.taxid {
                self.record_taxid(fallback, alignment);
            }
        }
    }

    fn record_taxid(&mut self, taxid: u32, alignment: &AlignmentRecord) {
        self.by_reference
            .entry(StrainCoverageReferenceKey {
                taxid,
                assembly_id: alignment.assembly_id,
                reference_id: alignment.reference_id,
            })
            .or_default()
            .push(alignment.span);
    }

    fn merge(mut self, mut other: Self) -> Self {
        if self.by_reference.len() < other.by_reference.len() {
            std::mem::swap(&mut self.by_reference, &mut other.by_reference);
        }
        for (key, mut intervals) in other.by_reference {
            self.by_reference
                .entry(key)
                .or_default()
                .append(&mut intervals);
        }
        self.candidate_reads = self.candidate_reads.saturating_add(other.candidate_reads);
        self.resolved_reads = self.resolved_reads.saturating_add(other.resolved_reads);
        self.ambiguous_reads = self.ambiguous_reads.saturating_add(other.ambiguous_reads);
        self.boundary_ambiguous_reads = self
            .boundary_ambiguous_reads
            .saturating_add(other.boundary_ambiguous_reads);
        self
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct StrainAssemblyCoverage {
    breadth: f64,
    chunk_breadth: f64,
    flex_chunk_breadth: f64,
    num_flex_chunks: usize,
    asm_length: usize,
}

#[inline]
fn strain_evidence_mode(cfg: &ProfilerConfig) -> EvidenceMode {
    if cfg.is_paired {
        EvidenceMode::Paired
    } else if cfg.sequencer.is_long_read() {
        EvidenceMode::Long
    } else {
        EvidenceMode::Single
    }
}

#[inline]
fn widen_strain_interval(
    (start, end): (usize, usize),
    flank: usize,
    reference_length: usize,
) -> (usize, usize) {
    (
        start.saturating_sub(flank).max(1),
        end.saturating_add(flank)
            .min(reference_length.saturating_add(1)),
    )
}

fn finalize_resolved_strain_coverage(
    intervals: ResolvedStrainIntervals,
    references: &[ReferenceMetadata],
    assembly_stats: &[AssemblyStats],
    fixed_flank_len: usize,
) -> HashMap<u32, TaxidCoverage> {
    let mut by_taxon_assembly: HashMap<(u32, u32), StrainAssemblyCoverage> = HashMap::new();

    let mut reference_intervals = intervals.by_reference.into_iter().collect::<Vec<_>>();
    reference_intervals.sort_unstable_by_key(|(key, _)| *key);
    for (key, spans) in reference_intervals {
        let Some(reference) = references.get(key.reference_id as usize) else {
            continue;
        };
        if reference.reference_id != key.reference_id
            || reference.assembly_id != key.assembly_id
            || reference.length == 0
        {
            continue;
        }
        let Some(assembly) = assembly_stats.get(reference.assembly_id as usize) else {
            continue;
        };
        if assembly.length == 0 {
            continue;
        }
        let flex_flank_len = ((assembly.length as f64).sqrt() as usize).max(1);
        let fixed_intervals: Vec<(usize, usize)> = spans
            .iter()
            .copied()
            .map(|span| widen_strain_interval(span, fixed_flank_len, reference.length))
            .collect();
        let flex_intervals: Vec<(usize, usize)> = spans
            .iter()
            .copied()
            .map(|span| widen_strain_interval(span, flex_flank_len, reference.length))
            .collect();
        let merged = merge_intervals(spans);
        let merged_fixed = merge_intervals(fixed_intervals);
        let merged_flex = merge_intervals(flex_intervals);
        let covered: usize = merged
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum();
        let fixed_covered: usize = merged_fixed
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum();
        let flex_covered: usize = merged_flex
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum();
        let asm_len = assembly.length.max(1) as f64;
        let estimated_chunks = (asm_len.sqrt() / 2.0).max(1.0).floor() as usize;
        let entry = by_taxon_assembly
            .entry((key.taxid, key.assembly_id))
            .or_default();
        entry.breadth += covered as f64 / asm_len;
        entry.chunk_breadth += fixed_covered as f64 / asm_len;
        entry.flex_chunk_breadth += flex_covered as f64 / asm_len;
        entry.num_flex_chunks = entry
            .num_flex_chunks
            .max(estimated_chunks.max(merged_flex.len().max(1)));
        entry.asm_length = assembly.length;
    }

    let mut coverage: HashMap<u32, TaxidCoverage> = HashMap::new();
    for ((taxid, _assembly_id), assembly) in by_taxon_assembly {
        let entry = coverage.entry(taxid).or_default();
        // Preserve the established strongest-assembly aggregation while
        // preventing intervals from different candidate taxa from being
        // unioned before this point.
        entry.breadth = entry.breadth.max(assembly.breadth.min(1.0));
        entry.chunk_breadth = entry.chunk_breadth.max(assembly.chunk_breadth.min(1.0));
        entry.flex_chunk_breadth = entry
            .flex_chunk_breadth
            .max(assembly.flex_chunk_breadth.min(1.0));
        entry.num_flex_chunks = entry.num_flex_chunks.max(assembly.num_flex_chunks.max(1));
        entry.asm_length = entry.asm_length.max(assembly.asm_length);
    }
    coverage
}

/// Taxonomic information for a taxon.
///
/// Contains the essential taxonomy data for classification output.
#[derive(Clone, Debug)]
struct TaxonInfo {
    /// NCBI taxonomy ID
    taxid: u32,

    /// Scientific name
    name: String,

    /// Taxonomic rank (species, genus, etc.)
    rank: String,
}

impl TaxonInfo {
    /// Create TaxonInfo from taxonomy database lookup.
    ///
    /// Returns None if the taxid is not found in the taxonomy.
    fn from_taxonomy(taxonomy: &Taxonomy, taxid: u32) -> Option<Self> {
        let name = taxonomy.get_name(taxid)?.to_string();
        let rank = taxonomy
            .get_rank(taxid)
            .map(|r| r.to_string())
            .unwrap_or_else(|| "no rank".to_string());
        Some(Self { taxid, name, rank })
    }
}

/// Taxonomic assignment for a single read.
///
/// For strain-level profiling, may contain multiple candidate taxa
/// if the read maps equally well to multiple strains.
#[derive(Clone, Debug)]
struct ReadAssignment {
    /// Read identifier
    read_name: Arc<str>,

    /// Taxonomic information for each candidate
    taxa: Vec<TaxonInfo>,

    /// Taxonomy IDs (parallel to taxa)
    taxids: Vec<u32>,

    /// Depth contributions (mapped_len / genome_size) per candidate
    depths: Vec<f64>,

    /// Mapped lengths per candidate
    maplens: Vec<f64>,

    /// Genome-specific support multiplier for each representative candidate.
    sampling_scales: Vec<f64>,
}

#[derive(Clone, Debug)]
struct StrainTotals {
    name: String,
    rank: String,
    support: ReadSupport,
    depth: f64,
    mapped_len: f64,
}

impl StrainTotals {
    fn new(name: String, rank: String) -> Self {
        Self {
            name,
            rank,
            support: ReadSupport::default(),
            depth: 0.0,
            mapped_len: 0.0,
        }
    }

    fn prior(name: String, rank: String) -> Self {
        Self {
            name,
            rank,
            support: ReadSupport {
                observed: 0.1,
                equivalent: 0.1,
            },
            depth: 0.001,
            mapped_len: 0.001,
        }
    }
}

#[derive(Clone, Debug)]
struct StrainAmbiguousTotals {
    taxa: Vec<TaxonInfo>,
    observed_count: f64,
    depths: Vec<f64>,
    mapped_lens: Vec<f64>,
    equivalent_support_sums: Vec<f64>,
}

/// Results from the taxonomy assignment phase.
struct AssignTaxaResult {
    /// Per-read taxonomic assignments
    assignments: Vec<ReadAssignment>,

    /// Count of alignments passing quality filters
    alignments_passing: usize,
}

/// Parse a single SAM alignment line.
///
/// Similar to community profiler's parsing but preserves strain-level
/// taxonomy information from reference names.
///
/// # Returns
///
/// ParsedAlignment with strain/species taxids and alignment coordinates
fn parse_alignment_line(line: &str, ctx: &WorkerContext) -> Option<ParsedAlignment> {
    let mut iter = line.split('\t');
    let read_name_raw = iter.next()?.trim();
    let flag = iter
        .next()
        .and_then(|field| field.parse::<u16>().ok())
        .unwrap_or(0);
    if (flag & 0x200) != 0 {
        return None;
    }
    let ref_name = iter.next()?.trim();
    let start_pos = iter
        .next()
        .and_then(|field| field.parse::<isize>().ok())
        .unwrap_or(1)
        .max(1) as usize;
    let _mapq = iter.next();
    let cigar = iter.next()?.trim();
    let rnext = iter.next().unwrap_or("*").trim();
    let mate_pos = iter
        .next()
        .and_then(|field| field.parse::<isize>().ok())
        .unwrap_or(0)
        .max(1) as usize;
    let _tlen = iter.next();
    let read_seq = iter.next().unwrap_or("");
    let _qual = iter.next();

    let read_name = if ctx.is_paired {
        read_name_raw
            .rsplit_once('/')
            .map(|(base, _)| base.to_string())
            .unwrap_or_else(|| read_name_raw.to_string())
    } else {
        read_name_raw.to_string()
    };

    let discard = |skipped_reference: bool| ParsedAlignment {
        read_name: read_name.clone(),
        count_alignment: false,
        record: None,
        skipped_reference,
    };

    if ref_name == "*" || cigar == "*" {
        return Some(discard(false));
    }

    if ctx.is_paired && rnext != "=" {
        return Some(discard(false));
    }

    let metadata = match ctx.ref_metadata.get(ref_name) {
        Some(meta) => meta,
        None => return Some(discard(true)),
    };
    if metadata.length == 0 {
        return Some(discard(true));
    }

    let mapped_len = read_seq.len();
    if mapped_len == 0 {
        return Some(discard(false));
    }

    let pair_start = if ctx.is_paired {
        start_pos.min(mate_pos)
    } else {
        start_pos
    };
    let mut start = pair_start.min(metadata.length.max(1));
    if start == 0 {
        start = 1;
    }
    let ref_end = metadata.length.saturating_add(1);
    let mut end = start.saturating_add(mapped_len);
    if end > ref_end {
        end = ref_end;
    }
    let span = (start, end);

    let mut count_alignment = true;
    if ctx.is_paired {
        if (flag & 0x4) != 0 {
            count_alignment = false;
        }
        if (flag & 0x40) == 0 {
            return Some(ParsedAlignment {
                read_name,
                count_alignment,
                record: None,
                skipped_reference: false,
            });
        }
    }

    Some(ParsedAlignment {
        read_name,
        count_alignment,
        skipped_reference: false,
        record: Some(ParsedAlignmentRecord {
            span,
            mapped_len,
            cigar: cigar.to_string(),
            reference_id: metadata.reference_id,
            assembly_id: metadata.assembly_id,
            taxid: metadata.taxid,
            fallback_taxid: metadata.fallback_taxid,
        }),
    })
}

struct ParsedMetamRecord<'a> {
    read_name: &'a str,
    record: AlignmentRecord,
}

fn strain_record_from_metam<'a>(
    mapping: MetamMappingRecord<'a>,
    ref_metadata: &HashMap<String, ReferenceMetadata>,
    is_paired: bool,
) -> Option<ParsedMetamRecord<'a>> {
    if is_paired && (!mapping.mate_ref_status_same_ref || mapping.pair_id == 0) {
        return None;
    }
    let metadata = ref_metadata.get(mapping.ref_id)?;
    if metadata.length == 0 || mapping.read_len == 0 || mapping.mapped_bases == 0 {
        return None;
    }
    let start = mapping.start_1based.clamp(1, metadata.length.max(1));
    let mut end = mapping.end_1based_exclusive.max(start.saturating_add(1));
    end = end.min(metadata.length.saturating_add(1));
    let read_name = if is_paired {
        strip_pair_suffix(mapping.read_id)
    } else {
        mapping.read_id
    };
    Some(ParsedMetamRecord {
        read_name,
        record: AlignmentRecord {
            mapped_len: mapping.mapped_bases,
            identity: mapping.identity,
            fraction: (mapping.mapped_bases as f64 / mapping.read_len as f64).min(1.0),
            reference_id: metadata.reference_id,
            assembly_id: metadata.assembly_id,
            span: (start, end),
            mate_number: mapping.mate_number,
            pair_id: mapping.pair_id,
            candidate_boundary: mapping.candidate_boundary,
            taxid: metadata.taxid,
            fallback_taxid: metadata.fallback_taxid,
        },
    })
}

/// Merge overlapping genomic intervals.
///
/// Combines intervals that overlap into single intervals.
/// Takes ownership and modifies the input vector.
fn merge_intervals(mut intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if intervals.len() <= 1 {
        return intervals;
    }
    intervals.sort_unstable_by_key(|interval| interval.0);
    let mut merged: Vec<(usize, usize)> = Vec::with_capacity(intervals.len());
    merged.push(intervals[0]);
    for interval in intervals.into_iter().skip(1) {
        if let Some(last) = merged.last_mut() {
            if interval.0 < last.1 {
                last.1 = last.1.max(interval.1);
            } else {
                merged.push(interval);
            }
        }
    }
    merged
}

/// Parse reference name parts to extract taxonomy information.
///
/// Reference format: `asm|strain_taxid|species_taxid|acc|...`
///
/// # Returns
///
/// Tuple of (assembly_id, primary_taxid, fallback_taxid)
/// - primary_taxid: strain-level if available, else species-level
/// - fallback_taxid: always species-level
fn parse_reference_header(parts: &[&str]) -> (String, Option<u32>, Option<u32>) {
    let asm = parts.first().copied().unwrap_or("").to_string();
    let strain_taxid = parts
        .get(1)
        .and_then(|value| value.parse::<u32>().ok())
        .filter(|value| *value > 0);
    let species_taxid = parts
        .get(2)
        .and_then(|value| value.parse::<u32>().ok())
        .filter(|value| *value > 0);
    if let Some(strain) = strain_taxid {
        (asm, Some(strain), species_taxid)
    } else if let Some(species) = species_taxid {
        (asm, Some(species), Some(species))
    } else {
        (asm, None, None)
    }
}

const STRAIN_PROB_CUTOFF: f64 = 1e-10;

/// Whether strain profiling computes coverage probabilities for this endpoint.
///
/// Strain profiling reports a strain-level taxid, which NCBI may label `strain`,
/// `subspecies` or `no rank`, and falls back to species. Those are the rows the
/// mode exists to produce, so their coverage shape is the shape worth testing,
/// while anything above species is left untested.
///
/// The rank label alone cannot decide this: `no rank` is used throughout the
/// taxonomy, above species as well as below it, so membership is resolved
/// through the lineage instead. An endpoint qualifies when it is a species or
/// descends from one.
fn strain_probability_endpoint(taxonomy: &Taxonomy, taxid: u32) -> bool {
    matches!(taxonomy.get_rank(taxid), Some("species"))
        || taxonomy
            .get_parents(taxid)
            .iter()
            .any(|(_, _, rank)| rank == "species")
}

fn strain_pvalue_filter_passes(
    probability_expected: bool,
    cov_prob: Option<f64>,
    chunk_prob: Option<f64>,
    strict: bool,
) -> bool {
    match (cov_prob, chunk_prob) {
        (Some(cov), Some(chunk)) => {
            if strict {
                cov >= STRAIN_PROB_CUTOFF && chunk >= STRAIN_PROB_CUTOFF
            } else {
                !(cov < STRAIN_PROB_CUTOFF && chunk < STRAIN_PROB_CUTOFF)
            }
        }
        // Only endpoints above species are excused: nothing is computed for
        // them, so strict mode must not reject a test that was never run. An
        // endpoint at or below species missing both probabilities means the
        // expectation was zero, which is a real failure and still rejects.
        (None, None) if strict && !probability_expected => true,
        _ => !strict,
    }
}

fn strain_profile_row_passes(
    row: &ProfileRow,
    cfg: &ProfilerConfig,
    min_chunk_breadth: f64,
    probability_expected: bool,
) -> bool {
    let expected_ok = row.expected_breadth > 0.0 && row.expected_flex_chunk > 0.0;
    let oebr = row.breadth / row.expected_breadth;
    let coebr = row.flex_chunk / row.expected_flex_chunk;
    let min_oebr = cfg.min_oebr.unwrap_or(0.75);
    let min_coebr = cfg.min_coebr.unwrap_or(0.75);
    let lower_oebr_ok = if cfg.min_oebr.is_some() {
        oebr >= min_oebr
    } else {
        oebr > min_oebr
    };
    let lower_coebr_ok = if cfg.min_coebr.is_some() {
        coebr >= min_coebr
    } else {
        coebr > min_coebr
    };
    let min_reads_ok = equivalent_read_count_passes_minimum(row.equivalent_mapcount, cfg.min_reads);
    let fixed_chunk_ok = if cfg.min_reads.is_some() && cfg.chunk_breadth.is_none() {
        true
    } else {
        row.fixed_chunk > min_chunk_breadth
    };
    expected_ok
        && oebr.is_finite()
        && coebr.is_finite()
        && min_reads_ok
        && fixed_chunk_ok
        && lower_oebr_ok
        && oebr <= 1.5
        && lower_coebr_ok
        && coebr <= 1.5
        && strain_pvalue_filter_passes(
            probability_expected,
            row.breadth_pvalue,
            row.chunk_pvalue,
            cfg.strict,
        )
}

fn write_strain_classifications(
    path: Option<&str>,
    resolved: &ResolvedClassifications,
    ambiguous_read_counts: &HashMap<Vec<u32>, usize>,
    ambiguous_reads: &HashMap<Vec<u32>, Vec<String>>,
    ambiguous_fraction: &AmbiguousFractions,
    retained_taxa: &HashSet<u32>,
    taxonomy: &Taxonomy,
) -> Result<usize> {
    let mut writer = path
        .map(|path| {
            File::create(path)
                .with_context(|| format!("failed to create {path}"))
                .map(BufWriter::new)
        })
        .transpose()?;
    let mut classified = 0usize;

    for group in &resolved.groups {
        if retained_taxa.contains(&group.key.taxid) {
            classified = classified.saturating_add(group.count);
        }
    }

    // Preserve the allocation-free count-only path. Deterministic ordering is
    // paid for only when the optional per-read classification file is saved.
    if let Some(writer) = writer.as_mut() {
        for group in resolved.sorted_groups() {
            if !retained_taxa.contains(&group.key.taxid) {
                continue;
            }
            let final_name = taxonomy
                .get_name(group.key.taxid)
                .unwrap_or("Taxonomy deprecated");
            let final_rank = taxonomy.get_rank(group.key.taxid).unwrap_or("no rank");
            let (source_names, source_taxids) = match &group.key.source {
                ClassificationSource::Direct => {
                    (final_name.to_string(), group.key.taxid.to_string())
                }
                ClassificationSource::Candidates(taxids) => (
                    taxids
                        .iter()
                        .map(|taxid| taxonomy.get_name(*taxid).unwrap_or("Taxonomy deprecated"))
                        .collect::<Vec<_>>()
                        .join(";"),
                    taxids
                        .iter()
                        .map(u32::to_string)
                        .collect::<Vec<_>>()
                        .join(";"),
                ),
            };
            for read in group.sorted_reads() {
                writeln!(
                    writer,
                    "{read}\t{final_name}\t{}\t{final_rank}\t{source_names}\t{source_taxids}\t1",
                    group.key.taxid,
                )?;
            }
        }
    }

    let mut keys = ambiguous_read_counts.keys().collect::<Vec<_>>();
    keys.sort_unstable();
    for candidate_key in keys {
        let Some(fractions) = ambiguous_fraction.get(candidate_key) else {
            continue;
        };
        let mut retained = fractions
            .iter()
            .filter_map(|(&taxid, (name, rank, fraction))| {
                (retained_taxa.contains(&taxid) && fraction.is_finite() && *fraction > 0.0)
                    .then_some((taxid, name.as_str(), rank.as_str(), *fraction))
            })
            .collect::<Vec<_>>();
        retained.sort_unstable_by_key(|(taxid, _, _, _)| *taxid);
        let total = retained
            .iter()
            .map(|(_, _, _, fraction)| *fraction)
            .sum::<f64>();
        if total <= f64::EPSILON || !total.is_finite() {
            continue;
        }
        for (_, _, _, fraction) in &mut retained {
            *fraction /= total;
        }
        retained.sort_by(|left, right| {
            right
                .3
                .total_cmp(&left.3)
                .then_with(|| left.0.cmp(&right.0))
        });
        classified = classified.saturating_add(ambiguous_read_counts[candidate_key]);

        let Some(writer) = writer.as_mut() else {
            continue;
        };
        let final_taxid = retained[0].0;
        let final_name = retained[0].1;
        let final_rank = retained[0].2;
        let names = retained
            .iter()
            .map(|(_, name, _, _)| *name)
            .collect::<Vec<_>>()
            .join(";");
        let candidate_taxids = retained
            .iter()
            .map(|(taxid, _, _, _)| taxid.to_string())
            .collect::<Vec<_>>()
            .join(";");
        let values = retained
            .iter()
            .map(|(_, _, _, fraction)| format!("{fraction:.5}"))
            .collect::<Vec<_>>()
            .join(";");
        if let Some(reads) = ambiguous_reads.get(candidate_key) {
            let mut reads = reads.iter().map(String::as_str).collect::<Vec<_>>();
            reads.sort_unstable();
            for read in reads {
                writeln!(
                    writer,
                    "{read}\t{final_name}\t{final_taxid}\t{final_rank}\t{names}\t{candidate_taxids}\t{values}"
                )?;
            }
        }
    }

    if let Some(writer) = writer.as_mut() {
        writer.flush()?;
    }
    Ok(classified)
}

/// Automatic chunk-breadth policy for strain profiling. It is kept separate
/// from the community curve.
fn strain_chunk_breadth_cutoff(num_reads: usize, illumina: bool) -> f64 {
    if illumina {
        if num_reads <= 100_000 {
            return num_reads as f64 / 1_000_000.0;
        }
        if num_reads <= 1_000_000 {
            let reads_in_million = num_reads as f64 / 1_000_000.0;
            return 0.11 * reads_in_million + 0.09;
        }
        let reads_in_million = num_reads as f64 / 1_000_000.0;
        return (0.2 + 0.382 * reads_in_million.log10()).min(0.95);
    }
    if num_reads <= 5_000 {
        0.0
    } else if num_reads <= 1_000_000 {
        0.3
    } else {
        0.5
    }
}

#[inline]
fn scaled_strain_chunk_breadth_cutoff(num_reads: usize, illumina: bool, factor: f64) -> f64 {
    scale_automatic_chunk_breadth(strain_chunk_breadth_cutoff(num_reads, illumina), factor)
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct StrainChunkBreadthThreshold {
    value: f64,
    automatic: Option<f64>,
    source: &'static str,
}

fn resolve_strain_chunk_breadth(
    cfg: &ProfilerConfig,
    basis_reads: usize,
    use_mapped_basis: bool,
) -> StrainChunkBreadthThreshold {
    if let Some(value) = cfg.chunk_breadth {
        return StrainChunkBreadthThreshold {
            value,
            automatic: None,
            source: "user-supplied --min-chunk-breadth",
        };
    }
    if cfg.min_reads.is_some() {
        return StrainChunkBreadthThreshold {
            value: 0.0,
            automatic: None,
            source: "disabled by --min-reads",
        };
    }
    if cfg.lowbiomass {
        return StrainChunkBreadthThreshold {
            value: 0.0,
            automatic: None,
            source: "disabled by --low-biomass",
        };
    }

    let illumina = cfg.sequencer.is_illumina();
    let automatic = strain_chunk_breadth_cutoff(basis_reads, illumina);
    let value = scaled_strain_chunk_breadth_cutoff(basis_reads, illumina, cfg.chunk_breadth_scale);
    let source = match (cfg.chunk_breadth_scale != 1.0, use_mapped_basis) {
        (true, true) => "scaled auto-estimate from mapped reads",
        (true, false) => "scaled auto-estimate from total reads",
        (false, true) => "auto-estimated from mapped reads",
        (false, false) => "auto-estimated from total reads",
    };
    StrainChunkBreadthThreshold {
        value,
        automatic: Some(automatic),
        source,
    }
}

fn log_strain_chunk_breadth_scaling(cfg: &ProfilerConfig, threshold: StrainChunkBreadthThreshold) {
    let Some(automatic) = threshold.automatic else {
        return;
    };
    if cfg.chunk_breadth_scale == 1.0 || automatic == 0.0 {
        return;
    }
    let uncapped = automatic * cfg.chunk_breadth_scale;
    if uncapped > 1.0 {
        log::info!(
            target: "PROFILER",
            "User chunk-breadth scale adjusted the automatic strain cutoff: {:.6} x {} = {:.6}, capped at 1.000000.",
            automatic,
            cfg.chunk_breadth_scale,
            threshold.value,
        );
    } else {
        log::info!(
            target: "PROFILER",
            "User chunk-breadth scale adjusted the automatic strain cutoff: {:.6} x {} = {:.6}.",
            automatic,
            cfg.chunk_breadth_scale,
            threshold.value,
        );
    }
}

/// EM assignment fractions for one candidate set: taxid -> (name, rank, fraction).
///
/// Insertion-ordered so the reported candidates keep the order the EM produced.
type TaxonFractions = IndexMap<u32, (String, String, f64)>;

/// EM assignment fractions for every candidate set, keyed by its taxids.
type AmbiguousFractions = HashMap<Vec<u32>, TaxonFractions>;

/// Normalize fraction values in a map to sum to 1.0.
fn normalize_fraction_map(map: &mut TaxonFractions) {
    let sum: f64 = map.values().map(|(_, _, frac)| *frac).sum();
    if sum <= f64::EPSILON {
        return;
    }
    for (_, _, frac) in map.values_mut() {
        *frac = round5(*frac / sum);
    }
}

/// Open alignment file with automatic compression detection.
///
/// Supports plain text, gzip (.gz), and xz (.xz) compression.
fn open_alignment_reader(path: &Path) -> Result<Box<dyn Read>> {
    let mut file =
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut magic = [0u8; 6];
    let read = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    if read >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        return Ok(Box::new(MultiGzDecoder::new(file)));
    }
    if read >= 6 && magic.starts_with(&[0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00]) {
        return Ok(Box::new(XzDecoder::new(file)));
    }
    Ok(Box::new(file))
}

/// Check if a taxonomy ID is in the lineage of any given host taxids.
///
/// Used for pathogen filtering - checks if the candidate host is
/// an ancestor of any known host for the pathogen.
fn lineage_contains(taxonomy: &Taxonomy, candidate: u32, host_taxids: &[u32]) -> bool {
    if host_taxids.contains(&candidate) {
        return true;
    }
    for taxid in host_taxids {
        let parents = taxonomy.get_parents(*taxid);
        if parents
            .iter()
            .any(|(parent_id, _, _)| *parent_id == candidate)
        {
            return true;
        }
    }
    false
}

/// Resolve taxonomy information with fallback.
///
/// Tries primary taxid first, then falls back to species-level if:
/// - Primary taxid not in taxonomy database
/// - Primary taxid has insufficient coverage
///
/// # Returns
///
/// Some((resolved_taxid, taxon_info)) if a valid taxon found, None otherwise
fn resolve_taxon_info(
    taxonomy: &Taxonomy,
    coverage_map: &HashMap<u32, TaxidCoverage>,
    primary: Option<u32>,
    fallback: Option<u32>,
    min_breadth: f64,
) -> Option<(u32, TaxonInfo)> {
    let mut candidates = Vec::new();
    if let Some(taxid) = primary {
        candidates.push(taxid);
    }
    if let Some(fallback_taxid) = fallback {
        if Some(fallback_taxid) != primary {
            candidates.push(fallback_taxid);
        }
    }
    for taxid in candidates {
        let coverage = match coverage_map.get(&taxid) {
            Some(value) => value,
            None => continue,
        };
        if coverage.breadth < min_breadth {
            continue;
        }
        if let Some(taxon) = TaxonInfo::from_taxonomy(taxonomy, taxid) {
            return Some((taxid, taxon));
        }
    }
    None
}

/// Strain-level taxonomy profiler.
///
/// Provides higher resolution than community profiler by distinguishing
/// between closely related strains within a species. Uses strain-level
/// taxonomy IDs when available, with fallback to species level.
///
/// # Pipeline
///
/// 1. Load SAM and extract per-assembly coverage metrics
/// 2. Assign reads to strains (or species as fallback)
/// 3. Resolve ambiguous reads using EM algorithm
/// 4. Apply LCA for reads mapping to undetected strains
/// 5. Filter based on coverage statistics
/// 6. Output strain-level profile and classifications
pub(super) struct StrainProfiler {
    /// Profiler configuration
    cfg: ProfilerConfig,

    /// Taxonomy database (thread-safe)
    taxonomy: Arc<Taxonomy>,
}

impl StrainProfiler {
    /// Create a new strain profiler.
    pub(super) fn new(cfg: ProfilerConfig, taxonomy: Arc<Taxonomy>) -> Result<Self> {
        Ok(Self { cfg, taxonomy })
    }

    /// Run the strain-level profiling pipeline.
    pub(super) fn run(&self) -> Result<()> {
        log::info!(target: "PROFILER", "Loading alignments ...");
        let loaded = if let Some(mappings) = self.cfg.mapper_mappings.as_deref() {
            let db = self
                .cfg
                .db_metadata
                .as_deref()
                .ok_or_else(|| anyhow!("missing reference metadata for native mappings"))?;
            self.load_metam(Path::new(mappings), db)
                .context("failed to parse Metam coordinate mappings")?
        } else {
            let sam_path = Path::new(&self.cfg.sam);
            if !sam_path.exists() {
                bail!("SAM file {} does not exist", sam_path.display());
            }
            self.load_sam(sam_path)
                .context("failed to parse SAM file")?
        };
        self.profile_loaded(loaded)
    }

    pub(super) fn run_with_metam_mappings(
        &self,
        receiver: MetamMappingReceiver,
        source: &str,
    ) -> Result<()> {
        log::info!(target: "PROFILER", "Loading structured Metam coordinate mappings ...");
        let db = self
            .cfg
            .db_metadata
            .as_deref()
            .ok_or_else(|| anyhow!("missing Metax database metadata for streamed mappings"))?;
        let loaded = self
            .load_metam_mappings(receiver, source, db)
            .context("failed to consume structured Metam coordinate mappings")?;
        self.profile_loaded(loaded)
    }

    fn load_metam_profile(&self, embedded: Option<MetamProfile>) -> Result<MetamProfile> {
        if let Some(profile) = embedded {
            return Ok(profile);
        }
        if let Some(shared) = &self.cfg.mapping_profile_memory {
            return load_shared_metam_profile(shared);
        }
        bail!("mapping file is missing embedded profile counters; rerun mapping with this Metax version")
    }

    fn profile_loaded(&self, loaded: SamLoadResult) -> Result<()> {
        let verbose_logging = self.cfg.verbose || self.cfg.very_verbose;
        let SamLoadResult {
            read_alignments,
            assembly_stats,
            taxid_coverage,
            min_breadth,
            min_chunk_breadth,
            has_scaled_genomes,
        } = loaded;
        log::info!(
            target: "PROFILER",
            "Loaded {} reads with alignments across {} assemblies",
            read_alignments.len(),
            assembly_stats.len()
        );
        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "Assigning taxonomy for {} reads with {} threads ...",
                read_alignments.len(),
                self.cfg.threads
            );
        }
        let AssignTaxaResult {
            assignments,
            alignments_passing,
        } = self.assign_taxa(
            &read_alignments,
            &assembly_stats,
            &taxid_coverage,
            min_breadth,
        )?;
        // Assignment now owns every per-read value needed by EM/output.
        // Release the much larger all-mappings buffer before constructing
        // grouped abundance state, matching the community profiler's lifetime.
        drop(read_alignments);
        let assigned_read_count = assignments.len();
        if self.cfg.very_verbose {
            let avg_alignments_per_read = if assigned_read_count > 0 {
                alignments_passing as f64 / assigned_read_count as f64
            } else {
                0.0
            };
            log::info!(
                target: "PROFILER",
                "Mappings passing identity/fraction filters: {} (across {} reads, avg {:.2})",
                alignments_passing,
                assigned_read_count,
                avg_alignments_per_read
            );
        }
        log::info!(
            target: "PROFILER",
            "Assigned taxonomy to {} reads",
            assigned_read_count
        );
        let database_is_fractional = self
            .cfg
            .db_metadata
            .as_ref()
            .map(|metadata| metadata.is_fractional)
            .unwrap_or(has_scaled_genomes);
        let skip_classification = database_is_fractional || !self.cfg.save_classifications;
        if verbose_logging {
            log::info!(target: "PROFILER", "Profiling taxonomy abundance ...");
        }
        let final_classified = self.write_outputs(
            assignments,
            &taxid_coverage,
            min_chunk_breadth,
            skip_classification,
        )?;
        if verbose_logging && final_classified != assigned_read_count {
            log::info!(
                target: "PROFILER",
                "Final profile filters retained taxonomic classifications for {} of {} candidate reads.",
                final_classified,
                assigned_read_count
            );
        }
        log::info!(
            target: "PROFILER",
            "Number of taxonomically classified reads: {}",
            final_classified
        );
        log::info!(target: "PROFILER", "Taxonomy profiling finished.");
        Ok(())
    }

    /// Load and parse SAM file, computing coverage metrics.
    ///
    /// Unlike community profiler, this tracks coverage per assembly
    /// for strain-level resolution.
    fn load_sam(&self, sam_path: &Path) -> Result<SamLoadResult> {
        let reader = open_alignment_reader(sam_path)?;
        let mut reader = BufReader::new(reader);
        let mut assembly_stats: Vec<AssemblyStats> = Vec::new();
        let mut ref_metadata: HashMap<String, ReferenceMetadata> = HashMap::new();
        let mut references: Vec<ReferenceMetadata> = Vec::new();
        let mut assembly_ids: HashMap<String, u32> = HashMap::new();
        let mut genome_sampling_scales: HashMap<String, f64> = HashMap::new();
        let mut has_scaled_genomes = false;
        let mut first_alignment: Option<String> = None;

        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('@') {
                if trimmed.starts_with("@SQ") {
                    let (sn, len) = parse_sam_sequence_record(trimmed)?;
                    let parts: Vec<&str> = sn.split('|').collect();
                    let (asm_name, taxid, fallback_taxid) = parse_reference_header(&parts);
                    let assembly_id = intern_assembly(
                        &asm_name,
                        &mut assembly_ids,
                        &mut assembly_stats,
                        "SAM",
                    )?;
                    let reference_id = u32::try_from(references.len())
                        .context("SAM contains more than 2^32 reference sequences")?;
                    let sampling_scale = sampling_scale_from_name_parts(&parts);
                    record_genome_sampling_scale(
                        &mut genome_sampling_scales,
                        &asm_name,
                        sampling_scale,
                    )?;
                    has_scaled_genomes |= sampling_scale > 1.0;
                    let metadata = ReferenceMetadata {
                        reference_id,
                        assembly_id,
                        length: len,
                        taxid,
                        fallback_taxid,
                    };
                    if ref_metadata.insert(sn.clone(), metadata.clone()).is_some() {
                        bail!("SAM header contains duplicate reference name '{}'", sn);
                    }
                    references.push(metadata);
                    let stats = &mut assembly_stats[assembly_id as usize];
                    stats.length += len;
                    stats.sampling_scale = sampling_scale;
                }
                continue;
            }
            first_alignment = Some(trimmed.to_string());
            break;
        }
        require_alignment_reference_dictionary(sam_path, ref_metadata.len())?;

        let worker_context = Arc::new(WorkerContext {
            ref_metadata: Arc::new(ref_metadata.clone()),
            is_paired: self.cfg.is_paired,
        });

        let worker_count = self.cfg.threads.max(1);
        let batch_size = self.cfg.batch_size.max(1);
        let (line_tx, line_rx) = bounded::<Vec<String>>(worker_count * 2);
        let (result_tx, result_rx) = unbounded::<WorkerResult>();

        let mut handles = Vec::new();
        for _ in 0..worker_count {
            let rx = line_rx.clone();
            let tx = result_tx.clone();
            let ctx = Arc::clone(&worker_context);
            handles.push(thread::spawn(move || {
                for chunk in rx.iter() {
                    let chunk_len = chunk.len();
                    let mut reads_seen: Vec<String> = Vec::with_capacity(chunk_len);
                    let mut alignments: Vec<(String, AlignmentRecord)> =
                        Vec::with_capacity(chunk_len);
                    let mut alignments_processed = 0usize;
                    let mut skipped_reference = 0usize;
                    for line in chunk {
                        if let Some(parsed) = parse_alignment_line(&line, &ctx) {
                            let ParsedAlignment {
                                read_name,
                                count_alignment,
                                record,
                                skipped_reference: was_skipped_reference,
                            } = parsed;
                            if was_skipped_reference {
                                skipped_reference += 1;
                            }
                            reads_seen.push(read_name.clone());
                            if count_alignment {
                                alignments_processed += 1;
                            }
                            if let Some(record) = record {
                                let (identity, _, fraction) = cigar_to_identity(&record.cigar);
                                let alignment = AlignmentRecord {
                                    mapped_len: record.mapped_len,
                                    identity,
                                    fraction,
                                    reference_id: record.reference_id,
                                    assembly_id: record.assembly_id,
                                    span: record.span,
                                    mate_number: if ctx.is_paired { 1 } else { 0 },
                                    pair_id: 0,
                                    candidate_boundary: CandidateBoundary::COMPLETE,
                                    taxid: record.taxid,
                                    fallback_taxid: record.fallback_taxid,
                                };
                                alignments.push((read_name, alignment));
                            }
                        }
                    }
                    if (!reads_seen.is_empty()
                        || !alignments.is_empty()
                        || alignments_processed > 0
                        || skipped_reference > 0)
                        && tx
                            .send(WorkerResult {
                                alignments,
                                reads_seen,
                                alignments_processed,
                                skipped_reference,
                            })
                            .is_err()
                    {
                        break;
                    }
                }
            }));
        }
        drop(result_tx);
        drop(line_rx);

        let mut read_alignments_map: StableReadMap<AlignmentRecord> = StableReadMap::default();
        let mut all_reads: HashSet<String> = HashSet::new();
        let mut num_alignment = 0usize;
        let mut skipped_reference = 0usize;
        let mut next_progress: usize = 1_000_000;

        let mut pending: Vec<String> = Vec::with_capacity(batch_size);
        if let Some(first) = first_alignment.take() {
            pending.push(first);
        }
        loop {
            line.clear();
            let read = reader.read_line(&mut line)?;
            if read == 0 {
                break;
            }
            let trimmed = line.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('@') {
                continue;
            }
            pending.push(trimmed.to_string());
            if pending.len() >= batch_size {
                line_tx.send(pending)?;
                pending = Vec::with_capacity(batch_size);
                while let Ok(result) = result_rx.try_recv() {
                    absorb_worker_result(
                        result,
                        &mut read_alignments_map,
                        &mut all_reads,
                        &mut num_alignment,
                        &mut skipped_reference,
                        &mut next_progress,
                    );
                }
            }
        }
        if !pending.is_empty() {
            line_tx.send(pending)?;
        }
        drop(line_tx);

        while let Ok(result) = result_rx.try_recv() {
            absorb_worker_result(
                result,
                &mut read_alignments_map,
                &mut all_reads,
                &mut num_alignment,
                &mut skipped_reference,
                &mut next_progress,
            );
        }

        for result in result_rx.iter() {
            absorb_worker_result(
                result,
                &mut read_alignments_map,
                &mut all_reads,
                &mut num_alignment,
                &mut skipped_reference,
                &mut next_progress,
            );
        }

        for handle in handles {
            join_alignment_worker(handle, "strain SAM parsing worker")?;
        }
        if skipped_reference > 0 {
            bail!(
                "SAM contains {} mapped record(s) with a reference absent from @SQ or incompatible reference metadata",
                skipped_reference
            );
        }

        log::info!(
            target: "PROFILER",
            "{} alignments processed",
            num_alignment
        );

        let num_reads = all_reads.len();
        let aligned_reads = read_alignments_map.len();
        let mapped_unit = if self.cfg.is_paired {
            "fragments"
        } else {
            "reads"
        };
        let mapping_rate = if num_reads == 0 {
            0.0
        } else {
            aligned_reads as f64 / num_reads as f64
        };
        // Honour explicit --by-mapped; otherwise auto-force mapped-read
        // basis for low-map Illumina non-pathogen runs. Identical chunk-
        // breadth equation; only the numerator basis changes.
        let by_aligned_effective = self.cfg.resolve_by_aligned(mapping_rate);
        let by_aligned_auto = !self.cfg.by_aligned && by_aligned_effective;
        let chunk_basis = if by_aligned_effective {
            aligned_reads
        } else {
            num_reads
        };
        let min_breadth = self.cfg.breadth.unwrap_or(0.0);
        let skip_auto_chunk_filter =
            self.cfg.min_reads.is_some() && self.cfg.chunk_breadth.is_none();
        let chunk_threshold =
            resolve_strain_chunk_breadth(&self.cfg, chunk_basis, by_aligned_effective);
        let min_chunk_breadth_value = chunk_threshold.value;
        let chunk_breadth_source = chunk_threshold.source;
        log_strain_chunk_breadth_scaling(&self.cfg, chunk_threshold);
        let min_chunk_display = format!("{:.6}", min_chunk_breadth_value);
        const INDENT: &str = "                                                ";
        if self.cfg.verbose || self.cfg.very_verbose {
            log::info!(
                target: "PROFILER",
                "Profiler parameters:\n{}min_breadth        {}\n{}min_chunk_breadth  {}\n{}min_reads          {}\n{}min_oebr           {}\n{}min_coebr          {}\n{}identity           {}\n{}mapped_len   {}\n{}fraction           {}\n{}batch_size         {}\n",
                INDENT,
                min_breadth,
                INDENT,
                min_chunk_display,
                INDENT,
                self.cfg
                    .min_reads
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "None".to_string()),
                INDENT,
                self.cfg
                    .min_oebr
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "None".to_string()),
                INDENT,
                self.cfg
                    .min_coebr
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "None".to_string()),
                INDENT,
                self.cfg.identity,
                INDENT,
                self.cfg.mapped_len,
                INDENT,
                self.cfg.fraction,
                INDENT,
                self.cfg.batch_size,
            );
        }

        let read_alignments: Vec<ReadAlignment> = read_alignments_map
            .into_iter()
            .map(|(read_name, alignments)| ReadAlignment {
                read_name: Arc::<str>::from(read_name),
                alignments,
            })
            .collect();
        let total_retained_alignments: usize = read_alignments
            .iter()
            .map(|read| read.alignments.len())
            .sum();
        let taxid_coverage =
            self.resolved_strain_coverage(&read_alignments, &references, &assembly_stats)?;

        log::info!(target: "PROFILER", "Total number of reads: {}", num_reads);
        log::info!(
            target: "PROFILER",
            "Number of mapped {}: {}",
            mapped_unit,
            read_alignments.len()
        );
        log::info!(
            target: "PROFILER",
            "Number of mapping records: {}",
            num_alignment
        );
        log::info!(
            target: "PROFILER",
            "Number of mapping records retained after parsing: {}",
            total_retained_alignments
        );
        if skip_auto_chunk_filter {
            log::info!(
                target: "PROFILER",
                "Automatic minimum chunk breadth filter disabled because --min-reads is set."
            );
        } else if self.cfg.lowbiomass && self.cfg.chunk_breadth.is_none() {
            log::info!(
                target: "PROFILER",
                "Automatic minimum chunk breadth filter disabled by --low-biomass."
            );
        } else if by_aligned_effective {
            if by_aligned_auto {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from mapped {} ({}) [auto-forced: low-map Illumina run, mapping_rate {:.5} < {:.5}].",
                    mapped_unit,
                    aligned_reads,
                    mapping_rate,
                    self.cfg.low_map_rate_threshold
                );
            } else {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from mapped {} ({}).",
                    mapped_unit,
                    aligned_reads
                );
            }
        } else {
            log::info!(
                target: "PROFILER",
                "Minimum chunk breadth estimated from total reads ({}).",
                num_reads
            );
        }
        log_resolved_profiler_parameters(
            &self.cfg,
            ResolvedProfilerParameters {
                min_breadth,
                min_chunk_breadth: min_chunk_breadth_value,
                chunk_breadth_source,
                chunk_breadth_basis_reads: chunk_basis,
                total_reads: num_reads,
                mapped_reads: aligned_reads,
                mapping_rate,
            },
        );

        Ok(SamLoadResult {
            read_alignments,
            assembly_stats,
            taxid_coverage,
            min_breadth,
            min_chunk_breadth: min_chunk_breadth_value,
            has_scaled_genomes,
        })
    }

    fn load_metam(&self, mappings_path: &Path, db: &DatabaseMetadata) -> Result<SamLoadResult> {
        if !mappings_path.exists() {
            bail!(
                "Metam mapping file {} does not exist",
                mappings_path.display()
            );
        }
        let reader = open_metam_mapping_reader(mappings_path)?;
        self.load_metam_reader(reader, &mappings_path.display().to_string(), db)
    }

    fn load_metam_reader<R: BufRead>(
        &self,
        reader: R,
        source: &str,
        db: &DatabaseMetadata,
    ) -> Result<SamLoadResult> {
        self.load_metam_records(source, db, |on_mapping| {
            for_each_metam_mapping_record(
                reader,
                source,
                self.cfg.lenient_mapping_input,
                on_mapping,
            )
        })
    }

    fn load_metam_mappings(
        &self,
        receiver: MetamMappingReceiver,
        source: &str,
        db: &DatabaseMetadata,
    ) -> Result<SamLoadResult> {
        self.load_metam_records(source, db, |on_mapping| {
            for_each_structured_mapping_record(receiver, source, on_mapping)
        })
    }

    fn load_metam_records<F>(
        &self,
        source: &str,
        db: &DatabaseMetadata,
        feed: F,
    ) -> Result<SamLoadResult>
    where
        F: FnOnce(
            &mut dyn for<'a> FnMut(MetamMappingRecord<'a>) -> Result<()>,
        ) -> Result<MetamParseStats>,
    {
        let mut assembly_stats: Vec<AssemblyStats> = Vec::new();
        let mut ref_metadata: HashMap<String, ReferenceMetadata> = HashMap::new();
        let mut references: Vec<ReferenceMetadata> = Vec::with_capacity(db.references.len());
        let mut assembly_ids: HashMap<String, u32> = HashMap::new();
        let mut genome_sampling_scales: HashMap<String, f64> = HashMap::new();
        let mut taxonomy_adjustments = ReferenceTaxonomyAdjustments::new();
        let mut taxonomy_resolver = ReferenceTaxonomyResolver::new(&self.taxonomy);

        for (reference_index, reference) in db.references.iter().enumerate() {
            let declared_taxonomy = parse_reference_taxonomy(&reference.taxonomy_field)
                .with_context(|| {
                    format!(
                        "invalid taxonomy metadata for reference '{}'",
                        reference.name
                    )
                })?;
            let taxonomy = taxonomy_resolver.resolve(declared_taxonomy, &reference.name)?;
            record_reference_taxonomy_adjustment(
                &mut taxonomy_adjustments,
                declared_taxonomy,
                taxonomy,
            );
            let asm_name = reference.genome.clone();
            let assembly_id = intern_assembly(
                &asm_name,
                &mut assembly_ids,
                &mut assembly_stats,
                "database",
            )?;
            let reference_id = u32::try_from(reference_index)
                .context("database contains more than 2^32 reference sequences")?;
            let taxid = Some(taxonomy.primary_taxid.unwrap_or(taxonomy.profile_taxid));
            let fallback_taxid = Some(taxonomy.profile_taxid);
            let sampling_scale = sampling_scale_from_sizes(
                reference.full_genome_size,
                reference.effective_indexed_size,
            );
            record_genome_sampling_scale(&mut genome_sampling_scales, &asm_name, sampling_scale)?;
            let metadata = ReferenceMetadata {
                reference_id,
                assembly_id,
                length: reference.length,
                taxid,
                fallback_taxid,
            };
            if ref_metadata
                .insert(reference.name.clone(), metadata.clone())
                .is_some()
            {
                bail!(
                    "database metadata contains duplicate reference name '{}'",
                    reference.name
                );
            }
            references.push(metadata);
            let stats = &mut assembly_stats[assembly_id as usize];
            stats.length += reference.length;
            stats.sampling_scale = sampling_scale;
        }
        log_reference_taxonomy_adjustments(
            &taxonomy_adjustments,
            self.cfg.verbose || self.cfg.very_verbose,
        );

        let mut read_alignments_map: StableReadMap<AlignmentRecord> = StableReadMap::default();
        let mut num_alignment = 0usize;

        let load_progress_options = self.cfg.mapping_load_progress_options();
        let load_progress = progress::counter(load_progress_options, "Loading mapping records");
        let parse_result = feed(&mut |mapping| {
            num_alignment += 1;
            progress::tick_count(
                &load_progress,
                "Loading mapping records",
                num_alignment as u64,
                load_progress_options.interval(),
            );
            if self.cfg.is_paired && (!mapping.mate_ref_status_same_ref || mapping.pair_id == 0) {
                return Ok(());
            }
            if let Some(parsed) =
                strain_record_from_metam(mapping, &ref_metadata, self.cfg.is_paired)
            {
                push_string_keyed_value(&mut read_alignments_map, parsed.read_name, parsed.record);
            }
            Ok(())
        });
        load_progress.finish_and_clear();
        let parse_stats = parse_result?;
        if parse_stats.skipped_malformed > 0 {
            log::warn!(
                target: "PROFILER",
                "Skipped {} malformed mapping records from {}.",
                parse_stats.skipped_malformed,
                source
            );
        }

        let aligned_reads = read_alignments_map.len();
        let mapped_unit = if self.cfg.is_paired {
            "fragments"
        } else {
            "reads"
        };
        let profile = self.load_metam_profile(parse_stats.profile)?;
        let num_reads = total_input_units(&profile, self.cfg.is_paired).max(aligned_reads);
        let mapping_rate = if num_reads == 0 {
            0.0
        } else {
            aligned_reads as f64 / num_reads as f64
        };
        let by_aligned_effective = self.cfg.resolve_by_aligned(mapping_rate);
        let by_aligned_auto = !self.cfg.by_aligned && by_aligned_effective;
        let chunk_basis = if by_aligned_effective {
            aligned_reads
        } else {
            num_reads
        };
        let min_breadth = self.cfg.breadth.unwrap_or(0.0);
        let skip_auto_chunk_filter =
            self.cfg.min_reads.is_some() && self.cfg.chunk_breadth.is_none();
        let chunk_threshold =
            resolve_strain_chunk_breadth(&self.cfg, chunk_basis, by_aligned_effective);
        let min_chunk_breadth_value = chunk_threshold.value;
        let chunk_breadth_source = chunk_threshold.source;
        log_strain_chunk_breadth_scaling(&self.cfg, chunk_threshold);

        let read_alignments: Vec<ReadAlignment> = read_alignments_map
            .into_iter()
            .map(|(read_name, alignments)| ReadAlignment {
                read_name: Arc::<str>::from(read_name),
                alignments,
            })
            .collect();
        let total_retained_alignments: usize = read_alignments
            .iter()
            .map(|read| read.alignments.len())
            .sum();
        let taxid_coverage =
            self.resolved_strain_coverage(&read_alignments, &references, &assembly_stats)?;

        log::info!(target: "PROFILER", "Total number of reads: {}", num_reads);
        log::info!(
            target: "PROFILER",
            "Number of mapped {}: {}",
            mapped_unit,
            read_alignments.len()
        );
        log::info!(target: "PROFILER", "Number of mapping records: {}", num_alignment);
        log::info!(
            target: "PROFILER",
            "Number of mapping records retained after parsing: {}",
            total_retained_alignments
        );
        log::info!(target: "PROFILER", "Metam reported mapped read records: {}", profile.mapped_reads);
        log::info!(target: "PROFILER", "Metam emitted mapping records: {}", profile.emitted_records);
        if skip_auto_chunk_filter {
            log::info!(
                target: "PROFILER",
                "Automatic minimum chunk breadth filter disabled because --min-reads is set."
            );
        } else if self.cfg.lowbiomass && self.cfg.chunk_breadth.is_none() {
            log::info!(
                target: "PROFILER",
                "Automatic minimum chunk breadth filter disabled by --low-biomass."
            );
        } else if by_aligned_effective {
            if by_aligned_auto {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from mapped {} ({}) [auto-forced: low-map Illumina run, mapping_rate {:.5} < {:.5}].",
                    mapped_unit,
                    aligned_reads,
                    mapping_rate,
                    self.cfg.low_map_rate_threshold
                );
            } else {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from mapped {} ({}).",
                    mapped_unit,
                    aligned_reads
                );
            }
        } else {
            log::info!(
                target: "PROFILER",
                "Minimum chunk breadth estimated from total reads ({}).",
                num_reads
            );
        }
        log_resolved_profiler_parameters(
            &self.cfg,
            ResolvedProfilerParameters {
                min_breadth,
                min_chunk_breadth: min_chunk_breadth_value,
                chunk_breadth_source,
                chunk_breadth_basis_reads: chunk_basis,
                total_reads: num_reads,
                mapped_reads: aligned_reads,
                mapping_rate,
            },
        );

        Ok(SamLoadResult {
            read_alignments,
            assembly_stats,
            taxid_coverage,
            min_breadth,
            min_chunk_breadth: min_chunk_breadth_value,
            has_scaled_genomes: genome_sampling_scales.values().any(|scale| *scale > 1.0),
        })
    }

    /// Build strain breadth from taxonomically resolved coordinate witnesses.
    ///
    /// Every qualified candidate remains available to abundance EM later.
    /// Coordinates are stricter: an exact cross-taxon tie contributes to no
    /// strain, while a unique or platform-margin winner contributes one
    /// representative placement per assembly. This prevents conserved reads
    /// from creating breadth simultaneously on every candidate strain.
    fn resolved_strain_coverage(
        &self,
        read_alignments: &[ReadAlignment],
        references: &[ReferenceMetadata],
        assembly_stats: &[AssemblyStats],
    ) -> Result<HashMap<u32, TaxidCoverage>> {
        let mode = strain_evidence_mode(&self.cfg);
        let quality = AssignmentQualityThresholds::from_config(&self.cfg);
        let batch_size = self.cfg.batch_size.max(1);
        let pool = ThreadPoolBuilder::new()
            .num_threads(self.cfg.threads.max(1))
            .build()?;
        let pb = progress::bar(
            self.cfg.progress,
            read_alignments.len() as u64,
            "Resolving strain coverage",
        );
        let progress = pb.clone();
        let intervals = pool.install(|| {
            read_alignments
                .par_chunks(batch_size)
                .map(|chunk| {
                    let mut local = ResolvedStrainIntervals::default();
                    let mut scratch = StrainCandidateScratch::default();
                    for read in chunk {
                        let candidate_count = scratch.build(&read.alignments, mode, quality).len();
                        if candidate_count == 0 {
                            continue;
                        }
                        local.candidate_reads = local.candidate_reads.saturating_add(1);
                        if scratch.terminal_boundary_incomplete {
                            local.boundary_ambiguous_reads =
                                local.boundary_ambiguous_reads.saturating_add(1);
                        }
                        let mut resolved = false;
                        for placement in scratch.resolved_placements(mode) {
                            local.record(read, placement);
                            resolved = true;
                        }
                        if resolved {
                            local.resolved_reads = local.resolved_reads.saturating_add(1);
                        } else if candidate_count > 1 && !scratch.terminal_boundary_incomplete {
                            local.ambiguous_reads = local.ambiguous_reads.saturating_add(1);
                        }
                    }
                    progress.inc(chunk.len() as u64);
                    local
                })
                .reduce(
                    ResolvedStrainIntervals::default,
                    ResolvedStrainIntervals::merge,
                )
        });
        pb.finish_and_clear();
        log::info!(
            target: "PROFILER",
            "Strain coordinate evidence: {} candidate reads, {} resolved witnesses, {} unresolved cross-taxon ties, {} candidate-boundary reads kept ambiguous.",
            intervals.candidate_reads,
            intervals.resolved_reads,
            intervals.ambiguous_reads,
            intervals.boundary_ambiguous_reads,
        );
        let fixed_flank_len = if self.cfg.is_paired { 1000 } else { 2000 };
        Ok(finalize_resolved_strain_coverage(
            intervals,
            references,
            assembly_stats,
            fixed_flank_len,
        ))
    }

    /// Assign taxonomy to reads based on alignments.
    ///
    /// For each read:
    /// 1. Filter alignments by identity, mapped_len, and fraction
    /// 2. Try to resolve strain-level taxid, fallback to species
    /// 3. Keep best (highest identity) alignment per taxon
    /// 4. Create multi-taxon assignment if read maps to multiple strains
    fn assign_taxa(
        &self,
        read_alignments: &[ReadAlignment],
        assembly_stats: &[AssemblyStats],
        taxid_coverage: &HashMap<u32, TaxidCoverage>,
        min_breadth: f64,
    ) -> Result<AssignTaxaResult> {
        let taxonomy = self.taxonomy.clone();
        let assignment_quality = AssignmentQualityThresholds::from_config(&self.cfg);
        let evidence_mode = strain_evidence_mode(&self.cfg);
        let batch_size = std::cmp::max(1, self.cfg.batch_size);
        let pool = ThreadPoolBuilder::new()
            .num_threads(self.cfg.threads.max(1))
            .build()?;
        let pb = progress::bar(
            self.cfg.progress,
            read_alignments.len() as u64,
            "Assigning taxonomy",
        );
        let progress = pb.clone();
        let passing_counter = Arc::new(AtomicUsize::new(0));
        let assignment_chunks: Vec<Vec<ReadAssignment>> = pool.install(|| {
            read_alignments
                .par_chunks(batch_size)
                .map_with(
                    (progress.clone(), Arc::clone(&passing_counter)),
                    |(progress_bar, counter), chunk| {
                        let mut local: Vec<ReadAssignment> = Vec::new();
                        let mut scratch = StrainCandidateScratch::default();
                        for read in chunk {
                            let mut taxa_depths: IndexMap<u32, (TaxonInfo, f64, f64, f64)> =
                                IndexMap::new();
                            scratch.build(&read.alignments, evidence_mode, assignment_quality);
                            counter.fetch_add(scratch.passing_mappings, AtomicOrdering::Relaxed);
                            for candidate in &scratch.output {
                                let Some(aln) =
                                    read.alignments.get(candidate.best.representative_index)
                                else {
                                    continue;
                                };
                                let stats = match assembly_stats.get(aln.assembly_id as usize) {
                                    Some(stats) => stats,
                                    None => continue,
                                };
                                if stats.length == 0 {
                                    continue;
                                }
                                let resolved = resolve_taxon_info(
                                    &taxonomy,
                                    taxid_coverage,
                                    Some(candidate.taxid),
                                    candidate.best.fallback_taxid,
                                    min_breadth,
                                );
                                let (_taxid, taxon) = match resolved {
                                    Some(value) => value,
                                    None => continue,
                                };
                                let asm_len = stats.length.max(1) as f64;
                                let depth = aln.mapped_len as f64 / asm_len;
                                let maplen = aln.mapped_len as f64;
                                taxa_depths
                                    .entry(taxon.taxid)
                                    .and_modify(|entry| {
                                        if depth > entry.1 {
                                            entry.0 = taxon.clone();
                                            entry.1 = depth;
                                            entry.2 = maplen;
                                            entry.3 = stats.sampling_scale;
                                        }
                                    })
                                    .or_insert((taxon, depth, maplen, stats.sampling_scale));
                            }
                            if taxa_depths.is_empty() {
                                continue;
                            }
                            let mut taxids: Vec<u32> = Vec::new();
                            let mut taxa: Vec<TaxonInfo> = Vec::new();
                            let mut depths: Vec<f64> = Vec::new();
                            let mut maplens: Vec<f64> = Vec::new();
                            let mut sampling_scales: Vec<f64> = Vec::new();
                            for (taxid, (taxon, depth, maplen, sampling_scale)) in
                                taxa_depths.into_iter()
                            {
                                taxids.push(taxid);
                                taxa.push(taxon);
                                depths.push(depth);
                                maplens.push(maplen);
                                sampling_scales.push(sampling_scale);
                            }
                            local.push(ReadAssignment {
                                read_name: Arc::clone(&read.read_name),
                                taxa,
                                taxids,
                                depths,
                                maplens,
                                sampling_scales,
                            });
                        }
                        progress_bar.inc(chunk.len() as u64);
                        local
                    },
                )
                .collect()
        });
        pb.finish_and_clear();
        let assignments: Vec<ReadAssignment> = assignment_chunks.into_iter().flatten().collect();
        let alignments_passing = passing_counter.load(AtomicOrdering::Relaxed);
        Ok(AssignTaxaResult {
            assignments,
            alignments_passing,
        })
    }

    /// Write profiling results to output files.
    ///
    /// Implements the full output generation pipeline:
    /// 1. Process unambiguous (single-taxon) assignments
    /// 2. Apply LCA to ambiguous reads without detected strains
    /// 3. Run EM algorithm to redistribute remaining ambiguous reads
    /// 4. Calculate final abundances and p-values
    /// 5. Filter taxa by coverage criteria
    /// 6. Write profile.txt and, when requested, classify.txt outputs
    ///
    /// # Arguments
    /// * `assignments` - Per-read taxonomic assignments
    /// * `taxid_coverage` - Coverage metrics per taxid
    /// * `min_chunk_breadth` - Minimum chunk breadth threshold
    /// * `skip_classification` - Do not create the optional classify.txt output
    fn write_outputs(
        &self,
        assignments: Vec<ReadAssignment>,
        taxid_coverage: &HashMap<u32, TaxidCoverage>,
        min_chunk_breadth: f64,
        skip_classification: bool,
    ) -> Result<usize> {
        let classification_path = format!("{}.classify.txt", self.cfg.outprefix);
        let out_prefix = if self.cfg.host.is_some() {
            format!("{}.pathogen", self.cfg.outprefix)
        } else {
            self.cfg.outprefix.clone()
        };
        let profile_path = format!("{}.profile.txt", out_prefix);
        let write_classifications = !skip_classification;
        let mut resolved_classifications = ResolvedClassifications::new(write_classifications);
        let mut profile_writer = WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(File::create(&profile_path)?);
        let mut raw_writer = if self.cfg.keep_raw {
            Some(
                WriterBuilder::new()
                    .has_headers(false)
                    .delimiter(b'\t')
                    .from_writer(File::create(format!("{}.rprofile.txt", out_prefix))?),
            )
        } else {
            None
        };
        let pathogen_map = if let Some(path) = &self.cfg.pathogen_host {
            Some(load_pathogen_table(path)?)
        } else {
            None
        };
        let host_taxid = self.cfg.host.as_deref().and_then(|h| h.parse::<u32>().ok());
        let mut unambiguous: HashMap<u32, StrainTotals> = HashMap::new();
        let mut ambiguous_counts: HashMap<Vec<u32>, StrainAmbiguousTotals> = HashMap::new();
        let mut ambiguous_read_counts: HashMap<Vec<u32>, usize> = HashMap::new();
        let mut ambiguous_reads: HashMap<Vec<u32>, Vec<String>> = HashMap::new();
        if self.cfg.save_classifications && skip_classification {
            log::warn!(
                target: "PROFILER",
                "Cannot write classify.txt for a fractional database because per-read assignments cannot be expanded safely to full-index equivalents."
            );
        } else if self.cfg.save_classifications {
            log::info!(
                target: "PROFILER",
                "Writing per-read strain classifications to {}",
                classification_path
            );
        }
        log::info!(
            target: "PROFILER",
            "Writing abundance profile to {}",
            profile_path
        );
        for assignment in &assignments {
            if assignment.taxids.len() == 1 {
                let taxon = &assignment.taxa[0];
                let taxid = assignment.taxids[0];
                let depth = assignment.depths[0];
                let maplen = assignment.maplens[0];
                let sampling_scale = assignment.sampling_scales[0];
                resolved_classifications.record_read(taxid, &assignment.read_name);
                let entry = unambiguous
                    .entry(taxid)
                    .or_insert_with(|| StrainTotals::new(taxon.name.clone(), taxon.rank.clone()));
                entry.support.add_weighted(1.0, sampling_scale, 1.0);
                entry.depth += depth;
                entry.mapped_len += maplen;
            } else {
                let key = assignment.taxids.clone();
                *ambiguous_read_counts.entry(key.clone()).or_default() += 1;
                if write_classifications {
                    ambiguous_reads
                        .entry(key.clone())
                        .or_default()
                        .push(assignment.read_name.to_string());
                }
                ambiguous_counts
                    .entry(key)
                    .and_modify(|entry| {
                        entry.observed_count += 1.0;
                        if entry.depths.len() < assignment.depths.len() {
                            entry.depths.resize(assignment.depths.len(), 0.0);
                        }
                        if entry.mapped_lens.len() < assignment.maplens.len() {
                            entry.mapped_lens.resize(assignment.maplens.len(), 0.0);
                        }
                        if entry.equivalent_support_sums.len() < assignment.sampling_scales.len() {
                            entry
                                .equivalent_support_sums
                                .resize(assignment.sampling_scales.len(), 0.0);
                        }
                        for (idx, depth) in assignment.depths.iter().enumerate() {
                            if let Some(val) = entry.depths.get_mut(idx) {
                                *val += depth;
                            }
                        }
                        for (idx, maplen) in assignment.maplens.iter().enumerate() {
                            if let Some(val) = entry.mapped_lens.get_mut(idx) {
                                *val += *maplen;
                            }
                        }
                        for (idx, scale) in assignment.sampling_scales.iter().enumerate() {
                            if let Some(val) = entry.equivalent_support_sums.get_mut(idx) {
                                *val += *scale;
                            }
                        }
                    })
                    .or_insert_with(|| StrainAmbiguousTotals {
                        taxa: assignment.taxa.clone(),
                        observed_count: 1.0,
                        depths: assignment.depths.clone(),
                        mapped_lens: assignment.maplens.clone(),
                        equivalent_support_sums: assignment.sampling_scales.clone(),
                    });
            }
        }

        let mut taxid_count_depth: HashMap<u32, StrainTotals> = unambiguous.clone();
        for taxids in ambiguous_counts.keys() {
            for taxid in taxids {
                taxid_count_depth
                    .entry(*taxid)
                    .or_insert_with(|| StrainTotals::prior(String::new(), String::new()));
            }
        }

        let mut ambiguous_lca: Vec<Vec<u32>> = Vec::new();
        let mut lca_groups = ambiguous_counts.iter().collect::<Vec<_>>();
        lca_groups.sort_unstable_by(|left, right| left.0.cmp(right.0));
        for &(taxids, ambiguous) in &lca_groups {
            if taxids.iter().any(|t| unambiguous.contains_key(t)) {
                continue;
            }
            let taxa_set: HashSet<u32> = taxids.iter().copied().collect();
            if let Some((lca_taxid, lca_name, lca_rank)) =
                self.taxonomy.get_majority_lca(&taxa_set, 0.7)
            {
                ambiguous_lca.push(taxids.clone());
                let entry = unambiguous
                    .entry(lca_taxid)
                    .or_insert_with(|| StrainTotals::new(lca_name.clone(), lca_rank.clone()));
                let equivalent_mean = ambiguous
                    .equivalent_support_sums
                    .iter()
                    .copied()
                    .sum::<f64>()
                    / ambiguous.equivalent_support_sums.len() as f64;
                entry
                    .support
                    .add_weighted(ambiguous.observed_count, equivalent_mean, 1.0);
                entry.depth +=
                    ambiguous.depths.iter().copied().sum::<f64>() / ambiguous.depths.len() as f64;
                entry.mapped_len += ambiguous.mapped_lens.iter().copied().sum::<f64>()
                    / ambiguous.mapped_lens.len() as f64;
                resolved_classifications.record_group(
                    lca_taxid,
                    taxids.clone(),
                    ambiguous_read_counts.get(taxids).copied().unwrap_or(0),
                    ambiguous_reads.remove(taxids),
                );
            }
        }
        for taxids in ambiguous_lca {
            ambiguous_counts.remove(&taxids);
            ambiguous_read_counts.remove(&taxids);
            ambiguous_reads.remove(&taxids);
        }

        let mut ambiguous_fraction: AmbiguousFractions = HashMap::new();
        let mut em_groups = ambiguous_counts.iter().collect::<Vec<_>>();
        em_groups.sort_unstable_by(|left, right| left.0.cmp(right.0));
        let mut em_taxids = taxid_count_depth.keys().copied().collect::<Vec<_>>();
        em_taxids.extend(unambiguous.keys().copied());
        em_taxids.sort_unstable();
        em_taxids.dedup();
        let mut total_taxid_count_depth: HashMap<u32, StrainTotals>;
        let mut iterations = 0u32;
        loop {
            total_taxid_count_depth = unambiguous.clone();
            for &(taxids, ambiguous) in &em_groups {
                let sum_depth: f64 = taxids
                    .iter()
                    .map(|taxid| {
                        taxid_count_depth
                            .get(taxid)
                            .map(|totals| totals.depth)
                            .unwrap_or(0.0)
                    })
                    .sum::<f64>();
                if sum_depth <= f64::EPSILON {
                    continue;
                }
                let mut paired: Vec<(u32, &TaxonInfo, f64, f64, f64, f64)> = taxids
                    .iter()
                    .enumerate()
                    .map(|(idx, taxid)| {
                        let taxon = &ambiguous.taxa[idx];
                        let depth_val = ambiguous.depths.get(idx).copied().unwrap_or(0.0);
                        let maplen_val = ambiguous.mapped_lens.get(idx).copied().unwrap_or(0.0);
                        let equivalent_sum = ambiguous
                            .equivalent_support_sums
                            .get(idx)
                            .copied()
                            .unwrap_or(ambiguous.observed_count);
                        let fraction = taxid_count_depth
                            .get(taxid)
                            .map(|totals| totals.depth / sum_depth)
                            .unwrap_or(0.0);
                        (
                            *taxid,
                            taxon,
                            depth_val,
                            maplen_val,
                            equivalent_sum,
                            fraction,
                        )
                    })
                    .collect();
                paired.sort_by(|left, right| {
                    right
                        .5
                        .total_cmp(&left.5)
                        .then_with(|| left.0.cmp(&right.0))
                });
                let mut fractions: TaxonFractions = IndexMap::new();
                let mut consumed = false;
                for (taxid, taxon, depth, maplen, equivalent_sum, fraction) in paired {
                    if fraction >= 0.99 {
                        let entry = total_taxid_count_depth.entry(taxid).or_insert_with(|| {
                            StrainTotals::prior(taxon.name.clone(), taxon.rank.clone())
                        });
                        entry.name = taxon.name.clone();
                        entry.rank = taxon.rank.clone();
                        entry
                            .support
                            .add_weighted(ambiguous.observed_count, equivalent_sum, 1.0);
                        entry.depth += depth;
                        entry.mapped_len += maplen;
                        fractions.insert(taxid, (taxon.name.clone(), taxon.rank.clone(), 1.0));
                        consumed = true;
                        break;
                    }
                    let entry = total_taxid_count_depth.entry(taxid).or_insert_with(|| {
                        StrainTotals::prior(taxon.name.clone(), taxon.rank.clone())
                    });
                    entry.name = taxon.name.clone();
                    entry.rank = taxon.rank.clone();
                    entry
                        .support
                        .add_weighted(ambiguous.observed_count, equivalent_sum, fraction);
                    entry.depth += depth * fraction;
                    entry.mapped_len += maplen * fraction;
                    fractions.insert(
                        taxid,
                        (taxon.name.clone(), taxon.rank.clone(), round5(fraction)),
                    );
                }
                if !consumed {
                    normalize_fraction_map(&mut fractions);
                }
                ambiguous_fraction.insert(taxids.clone(), fractions);
            }
            let numerator: f64 = em_taxids
                .iter()
                .map(|taxid| {
                    let depth = total_taxid_count_depth
                        .get(taxid)
                        .map_or(0.0, |totals| totals.depth);
                    let old_depth = taxid_count_depth
                        .get(taxid)
                        .map(|previous| previous.depth)
                        .unwrap_or(0.0);
                    (depth - old_depth).abs()
                })
                .sum::<f64>();
            let denominator: f64 = em_taxids
                .iter()
                .map(|taxid| taxid_count_depth.get(taxid).map_or(0.0, |totals| totals.depth))
                .sum::<f64>()
                .max(1e-12);
            if numerator / denominator <= 1e-6 {
                taxid_count_depth = total_taxid_count_depth.clone();
                break;
            }
            iterations += 1;
            if iterations as usize >= self.cfg.max_em_iterations {
                log::warn!(
                    target: "PROFILER",
                    "Strain EM reached the maximum of {} iterations before convergence.",
                    self.cfg.max_em_iterations
                );
                taxid_count_depth = total_taxid_count_depth.clone();
                break;
            }
            taxid_count_depth = total_taxid_count_depth.clone();
        }
        log::info!(
            target: "PROFILER",
            "Finished strain EM after {} iterations",
            iterations
        );

        let norm_factor: f64 = 0.01
            * em_taxids
                .iter()
                .map(|taxid| taxid_count_depth.get(taxid).map_or(0.0, |totals| totals.depth))
                .sum::<f64>();
        let num_taxa = taxid_count_depth.len();
        let mut out_taxa_list: Vec<ProfileRow> = Vec::new();
        let mut entries: Vec<(u32, StrainTotals)> = taxid_count_depth
            .iter()
            .map(|(taxid, totals)| (*taxid, totals.clone()))
            .collect();
        entries.sort_by(|left, right| {
            right
                .1
                .depth
                .total_cmp(&left.1.depth)
                .then_with(|| left.0.cmp(&right.0))
        });

        for (taxid, totals) in entries {
            let mut name = totals.name;
            let mut rank = totals.rank;
            let mapcount = totals.support.observed;
            if mapcount < 0.5 {
                continue;
            }
            let equivalent_mapcount = totals.support.equivalent;
            let depth = totals.depth;
            let maplen = totals.mapped_len;
            if name.is_empty() {
                if let Some(taxon_name) = self.taxonomy.get_name(taxid) {
                    name = taxon_name.to_string();
                }
            }
            if rank.is_empty() {
                rank = self
                    .taxonomy
                    .get_rank(taxid)
                    .unwrap_or("no rank")
                    .to_string();
            }
            if name == "Taxonomy deprecated" {
                continue;
            }
            let abundance = if norm_factor <= f64::EPSILON {
                0.0
            } else {
                depth / norm_factor
            };
            let coverage = taxid_coverage.get(&taxid).cloned().unwrap_or_default();
            let breadth = coverage.breadth.min(1.0);
            let fixed_chunk = coverage.chunk_breadth.min(1.0);
            let flex_chunk = coverage.flex_chunk_breadth.min(1.0);
            let expected_breadth = crate::profiler::evidence::expected_breadth(depth, mapcount);
            let genome_len = if depth <= f64::EPSILON {
                coverage.asm_length as f64
            } else {
                maplen / depth
            };
            let fallback_chunks = if genome_len <= 0.0 {
                1
            } else {
                let raw = (genome_len.sqrt() / 2.0).floor();
                raw.max(1.0) as usize
            };
            let num_flex_chunks = if coverage.num_flex_chunks > 0 {
                coverage.num_flex_chunks
            } else {
                fallback_chunks
            };
            let expected_flex_chunk = if num_flex_chunks == 0 {
                0.0
            } else {
                1.0 - (1.0 - 1.0 / num_flex_chunks as f64).powf(mapcount)
            };
            let probability_expected = strain_probability_endpoint(&self.taxonomy, taxid);
            let cov_prob = if probability_expected {
                coverage_deviation_pvalue(breadth, expected_breadth, Some(depth), mapcount, num_taxa)
            } else {
                None
            };
            let chunkcov_prob = if probability_expected {
                coverage_deviation_pvalue(flex_chunk, expected_flex_chunk, None, mapcount, num_taxa)
            } else {
                None
            };
            let annotations = if let Some(host_taxid) = host_taxid {
                if let Some(map) = pathogen_map.as_ref() {
                    if let Some(entry) = map.get(&taxid) {
                        let host_taxids: Vec<&str> = entry.host_taxids.split(';').collect();
                        let host_taxid_list: Vec<u32> = host_taxids
                            .iter()
                            .filter_map(|t| t.parse::<u32>().ok())
                            .collect();
                        // An unknown host list annotates everything; otherwise the
                        // host must appear in the list directly, or through its
                        // lineage, which is the costlier check and so comes last.
                        let report_it = entry.host_taxids == "unknown"
                            || host_taxid_list.contains(&host_taxid)
                            || lineage_contains(&self.taxonomy, host_taxid, &host_taxid_list);
                        if report_it {
                            Some(ProfileAnnotations {
                                host_names: entry.host_names.clone(),
                                host_taxids: entry.host_taxids.clone(),
                                diseases: entry.diseases.clone(),
                            })
                        } else {
                            Some(ProfileAnnotations {
                                host_names: "unknown".to_string(),
                                host_taxids: "unknown".to_string(),
                                diseases: "unknown".to_string(),
                            })
                        }
                    } else {
                        Some(ProfileAnnotations {
                            host_names: "unknown".to_string(),
                            host_taxids: "unknown".to_string(),
                            diseases: "unknown".to_string(),
                        })
                    }
                } else {
                    Some(ProfileAnnotations {
                        host_names: "unknown".to_string(),
                        host_taxids: "unknown".to_string(),
                        diseases: "unknown".to_string(),
                    })
                }
            } else {
                None
            };
            out_taxa_list.push(ProfileRow {
                taxname: name,
                taxid,
                rank,
                equivalent_mapcount,
                depth,
                abundance,
                breadth,
                expected_breadth,
                breadth_pvalue: cov_prob,
                fixed_chunk,
                flex_chunk,
                expected_flex_chunk,
                chunk_pvalue: chunkcov_prob,
                annotations,
            });
        }

        let raw_taxa_list = out_taxa_list.clone();
        let mut normalized: Vec<ProfileRow> = out_taxa_list
            .into_iter()
            .filter(|row| {
                strain_profile_row_passes(
                    row,
                    &self.cfg,
                    min_chunk_breadth,
                    strain_probability_endpoint(&self.taxonomy, row.taxid),
                )
            })
            .collect();
        let total_abundance = normalized.iter().map(|row| row.abundance).sum::<f64>();
        for row in &mut normalized {
            row.normalize_abundance(total_abundance);
        }
        let retained_taxa = normalized
            .iter()
            .map(|row| row.taxid)
            .collect::<HashSet<_>>();
        let classified_reads = write_strain_classifications(
            write_classifications.then_some(classification_path.as_str()),
            &resolved_classifications,
            &ambiguous_read_counts,
            &ambiguous_reads,
            &ambiguous_fraction,
            &retained_taxa,
            &self.taxonomy,
        )?;

        for row in &normalized {
            row.write(&mut profile_writer)?;
        }
        profile_writer.flush()?;

        if let Some(writer) = raw_writer.as_mut() {
            for row in raw_taxa_list {
                row.write(writer)?;
            }
            writer.flush()?;
        }
        Ok(classified_reads)
    }
}

#[cfg(test)]
mod tests {
    use super::{
        finalize_resolved_strain_coverage, scaled_strain_chunk_breadth_cutoff,
        strain_chunk_breadth_cutoff, strain_probability_endpoint, strain_profile_row_passes,
        strain_pvalue_filter_passes, parse_alignment_line, AlignmentRecord, AssemblyStats,
        AssignmentQualityThresholds, CandidateBoundary, EvidenceMode, ProfileRow, ReadAlignment,
        ReferenceMetadata, ResolvedStrainIntervals, StrainCandidatePlacement,
        StrainCandidateScratch, Taxonomy, WorkerContext,
    };
    use crate::profiler::{ExactTiePolicy, ProfileMode, ProfilerConfig};
    use crate::sequencer::Sequencer;
    use std::collections::{HashMap, HashSet};

    fn approx_equal(lhs: f64, rhs: f64) {
        assert!(
            (lhs - rhs).abs() < 1e-9,
            "left {} not approximately equal to right {}",
            lhs,
            rhs
        );
    }

    fn quality_test_record(identity: f64, mapped_len: usize, fraction: f64) -> AlignmentRecord {
        AlignmentRecord {
            mapped_len,
            identity,
            fraction,
            reference_id: 0,
            assembly_id: 0,
            span: (100, 100 + mapped_len),
            mate_number: 0,
            pair_id: 0,
            candidate_boundary: CandidateBoundary::COMPLETE,
            taxid: Some(1),
            fallback_taxid: Some(1),
        }
    }

    fn output_test_config() -> ProfilerConfig {
        ProfilerConfig {
            sam: String::new(),
            mapper_mappings: None,
            mapping_profile_memory: None,
            db_metadata: None,
            sequencer: Sequencer::Illumina,
            batch_size: 1,
            max_em_iterations: 200,
            lenient_mapping_input: false,
            is_paired: false,
            exact_ties: ExactTiePolicy::Ambiguous,
            identity: 0.97,
            mapped_len: 50,
            breadth: Some(0.0),
            chunk_breadth: Some(0.0),
            chunk_breadth_scale: 1.0,
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            fraction: 0.8,
            lowbiomass: false,
            keep_raw: false,
            save_classifications: false,
            by_aligned: false,
            genus_identity: 0.8,
            low_map_rate_threshold: 0.3,
            genus_fallback: false,
            strict: false,
            pathogen_host: None,
            host: None,
            threads: 1,
            verbose: false,
            very_verbose: false,
            progress: crate::progress::ProgressOptions::default(),
            suppress_mapping_load_progress: false,
            outprefix: "test".to_string(),
            mode: ProfileMode::Default,
            dmp_dir: None,
            strain: true,
        }
    }

    #[test]
    fn strain_sam_parser_marks_undeclared_references() {
        let context = WorkerContext {
            ref_metadata: std::sync::Arc::new(HashMap::new()),
            is_paired: false,
        };
        let parsed = parse_alignment_line(
            "read-1\t0\tmissing-ref\t1\t60\t4=\t*\t0\t0\tACGT\tIIII",
            &context,
        )
        .expect("well-formed SAM row should produce a discard record");
        assert!(parsed.skipped_reference);
        assert!(parsed.record.is_none());
    }

    #[test]
    fn assignment_quality_gate_is_shared_by_coverage_and_assignment() {
        let short_read_gate = AssignmentQualityThresholds {
            min_identity: 0.97,
            min_mapped_len: 50,
            min_fraction: 0.80,
            coverage_or_gate: false,
        };
        assert!(short_read_gate.passes(&quality_test_record(0.97, 50, 0.80)));
        assert!(!short_read_gate.passes(&quality_test_record(0.96, 150, 1.0)));
        assert!(!short_read_gate.passes(&quality_test_record(1.0, 49, 1.0)));
        assert!(!short_read_gate.passes(&quality_test_record(1.0, 150, 0.79)));

        let long_read_gate = AssignmentQualityThresholds {
            coverage_or_gate: true,
            ..short_read_gate
        };
        assert!(long_read_gate.passes(&quality_test_record(1.0, 49, 0.80)));
        assert!(long_read_gate.passes(&quality_test_record(1.0, 50, 0.10)));
    }

    #[test]
    fn unresolved_cross_taxon_tie_contributes_no_strain_coordinates() {
        let quality = AssignmentQualityThresholds {
            min_identity: 0.90,
            min_mapped_len: 50,
            min_fraction: 0.50,
            coverage_or_gate: false,
        };
        let mut first = quality_test_record(0.99, 100, 1.0);
        first.taxid = Some(11);
        let mut second = quality_test_record(0.99, 100, 1.0);
        second.taxid = Some(12);
        second.reference_id = 1;
        second.assembly_id = 1;
        let mut scratch = StrainCandidateScratch::default();
        assert_eq!(
            scratch
                .build(&[first, second], EvidenceMode::Single, quality)
                .len(),
            2
        );
        assert_eq!(scratch.resolved_placements(EvidenceMode::Single).count(), 0);
    }

    #[test]
    fn candidate_boundary_blocks_direct_strain_coordinate_resolution() {
        let quality = AssignmentQualityThresholds {
            min_identity: 0.90,
            min_mapped_len: 50,
            min_fraction: 0.50,
            coverage_or_gate: false,
        };
        let mut record = quality_test_record(0.99, 100, 1.0);
        record.candidate_boundary = CandidateBoundary::from_bits(2).unwrap();
        let mut scratch = StrainCandidateScratch::default();
        assert_eq!(
            scratch
                .build(&[record], EvidenceMode::Single, quality)
                .len(),
            1
        );
        assert_eq!(scratch.resolved_placements(EvidenceMode::Single).count(), 0);
    }

    #[test]
    fn margin_resolved_taxon_keeps_only_its_assembly_placements() {
        let quality = AssignmentQualityThresholds {
            min_identity: 0.90,
            min_mapped_len: 50,
            min_fraction: 0.50,
            coverage_or_gate: false,
        };
        let mut best = quality_test_record(0.99, 100, 1.0);
        best.taxid = Some(11);
        let mut weaker = quality_test_record(0.98, 100, 1.0);
        weaker.taxid = Some(12);
        weaker.reference_id = 1;
        weaker.assembly_id = 1;
        let mut scratch = StrainCandidateScratch::default();
        scratch.build(&[best, weaker], EvidenceMode::Single, quality);
        let resolved: Vec<_> = scratch.resolved_placements(EvidenceMode::Single).collect();
        assert_eq!(resolved.len(), 1);
        assert_eq!(resolved[0].taxid, 11);
    }

    #[test]
    fn paired_candidate_uses_r1_identity_and_requires_covered_r2() {
        let quality = AssignmentQualityThresholds {
            min_identity: 0.97,
            min_mapped_len: 50,
            min_fraction: 0.80,
            coverage_or_gate: false,
        };
        let mut mate1 = quality_test_record(0.99, 100, 1.0);
        mate1.taxid = Some(11);
        mate1.mate_number = 1;
        mate1.pair_id = 7;
        let mut mate2 = quality_test_record(0.50, 100, 1.0);
        mate2.taxid = Some(11);
        mate2.mate_number = 2;
        mate2.pair_id = 7;
        let mut scratch = StrainCandidateScratch::default();
        assert_eq!(
            scratch
                .build(&[mate2, mate1], EvidenceMode::Paired, quality)
                .len(),
            1
        );
        let resolved: Vec<_> = scratch.resolved_placements(EvidenceMode::Paired).collect();
        assert_eq!(resolved.len(), 1);
        assert_eq!(resolved[0].representative_index, 1);
    }

    #[test]
    fn resolved_strain_interval_also_supports_its_species_fallback() {
        let alignment = AlignmentRecord {
            mapped_len: 100,
            identity: 1.0,
            fraction: 1.0,
            reference_id: 0,
            assembly_id: 0,
            span: (100, 200),
            mate_number: 0,
            pair_id: 0,
            candidate_boundary: CandidateBoundary::COMPLETE,
            taxid: Some(11),
            fallback_taxid: Some(101),
        };
        let read = ReadAlignment {
            read_name: "read".into(),
            alignments: vec![alignment],
        };
        let mut intervals = ResolvedStrainIntervals::default();
        intervals.record(
            &read,
            StrainCandidatePlacement {
                taxid: 11,
                fallback_taxid: Some(101),
                assembly_id: 0,
                reference_id: 0,
                representative_index: 0,
                errors: 0.0,
                long_score: 100.0,
                mapped_total: 100,
            },
        );
        let references = vec![ReferenceMetadata {
            reference_id: 0,
            assembly_id: 0,
            length: 1_000,
            taxid: Some(11),
            fallback_taxid: Some(101),
        }];
        let assembly_stats = vec![AssemblyStats {
            length: 1_000,
            sampling_scale: 1.0,
        }];
        let coverage =
            finalize_resolved_strain_coverage(intervals, &references, &assembly_stats, 100);
        approx_equal(coverage[&11].breadth, 0.1);
        approx_equal(coverage[&101].breadth, 0.1);
        assert_eq!(
            coverage[&11].chunk_breadth.to_bits(),
            coverage[&101].chunk_breadth.to_bits()
        );
    }

    #[test]
    fn strain_chunk_breadth_keeps_its_established_platform_policy() {
        approx_equal(strain_chunk_breadth_cutoff(0, true), 0.0);
        approx_equal(strain_chunk_breadth_cutoff(100_000, true), 0.10);
        approx_equal(strain_chunk_breadth_cutoff(1_000_000, true), 0.20);
        approx_equal(strain_chunk_breadth_cutoff(10_000_000, true), 0.582);
        approx_equal(strain_chunk_breadth_cutoff(5_000, false), 0.0);
        approx_equal(strain_chunk_breadth_cutoff(5_001, false), 0.3);
        approx_equal(strain_chunk_breadth_cutoff(1_000_001, false), 0.5);
    }

    #[test]
    fn strain_chunk_breadth_scale_preserves_default_and_caps_adjustments() {
        let raw = strain_chunk_breadth_cutoff(10_000_000, true);
        assert_eq!(
            scaled_strain_chunk_breadth_cutoff(10_000_000, true, 1.0).to_bits(),
            raw.to_bits()
        );
        approx_equal(
            scaled_strain_chunk_breadth_cutoff(10_000_000, true, 0.5),
            raw * 0.5,
        );
        approx_equal(
            scaled_strain_chunk_breadth_cutoff(10_000_000, true, 10.0),
            1.0,
        );
    }

    /// Lineage covering both placements of `no rank`: below species (1002) and
    /// above it (20), alongside the strain and subspecies endpoints.
    fn probability_endpoint_taxonomy() -> Taxonomy {
        let nodes = [
            (1u32, (1u32, "no rank")),
            (2, (1, "superkingdom")),
            (10, (2, "genus")),
            (20, (2, "no rank")),
            (100, (10, "species")),
            (1000, (100, "strain")),
            (1001, (100, "subspecies")),
            (1002, (100, "no rank")),
        ]
        .into_iter()
        .map(|(id, (parent, rank))| (id, (parent, rank.to_string())))
        .collect();
        Taxonomy::from_parts(nodes, HashMap::new(), HashMap::new(), HashSet::new())
    }

    #[test]
    fn probability_endpoints_are_resolved_through_the_lineage() {
        let taxonomy = probability_endpoint_taxonomy();
        // Species itself and everything descending from it.
        for taxid in [100, 1000, 1001, 1002] {
            assert!(strain_probability_endpoint(&taxonomy, taxid), "{taxid}");
        }
        // Above species, including the `no rank` nodes the rank label alone
        // would have mistaken for strain endpoints.
        for taxid in [1, 2, 10, 20] {
            assert!(!strain_probability_endpoint(&taxonomy, taxid), "{taxid}");
        }
        // An unknown taxid has no lineage to place it below a species.
        assert!(!strain_probability_endpoint(&taxonomy, 99_999));
    }

    #[test]
    fn strict_filter_keeps_non_species_without_pvalues() {
        // Probabilities are computed for every endpoint at or below species, so
        // a missing pair there means the expectation was zero -- a real failure
        // that still rejects.
        assert!(!strain_pvalue_filter_passes(true, None, None, true));

        // Above species nothing is computed, so the absence is expected and
        // must not be treated as evidence against the row.
        assert!(strain_pvalue_filter_passes(false, None, None, true));

        assert!(!strain_pvalue_filter_passes(true, Some(1e-3), None, true));
        assert!(!strain_pvalue_filter_passes(
            true,
            Some(1e-12),
            Some(1e-3),
            true
        ));
    }

    #[test]
    fn strain_filter_uses_unrounded_chunk_breadth() {
        let cfg = output_test_config();
        let row = |fixed_chunk| ProfileRow {
            taxname: "Species one".to_string(),
            taxid: 101,
            rank: "species".to_string(),
            equivalent_mapcount: 10.0,
            depth: 0.1,
            abundance: 1.0,
            breadth: 0.1,
            expected_breadth: 0.1,
            breadth_pvalue: Some(1.0),
            fixed_chunk,
            flex_chunk: 0.1,
            expected_flex_chunk: 0.1,
            chunk_pvalue: Some(1.0),
            annotations: None,
        };
        assert!(!strain_profile_row_passes(&row(0.499_996), &cfg, 0.5, true));
        assert!(strain_profile_row_passes(&row(0.500_004), &cfg, 0.5, true));
    }
}

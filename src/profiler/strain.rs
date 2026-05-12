//! Strain Profiler - Strain-level taxonomic abundance estimation
//!
//! This module implements strain-level profiling that can distinguish between
//! closely related strains within a species. It uses a similar pipeline to
//! the community profiler but operates at the strain/assembly level rather
//! than the species level.
//!
//! # Key Differences from Community Profiler
//!
//! - **Resolution**: Assigns reads to individual strain rather
//!   than aggregating to species level
//! - **Taxonomy fallback**: If strain taxonomy lookup fails, falls back to
//!   species-level assignment
//! - **Coverage tracking**: Tracks coverage per assembly/strain rather than per species

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::{
    atomic::{AtomicUsize, Ordering as AtomicOrdering},
    Arc,
};
use std::thread;

use anyhow::{bail, Context, Result};
use crossbeam_channel::{bounded, unbounded};
use csv::WriterBuilder;
use flate2::read::MultiGzDecoder;
use indexmap::IndexMap;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
// use once_cell::sync::Lazy;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
// use regex::Regex;
use statrs::function::beta::beta_reg;
use statrs::function::erf::erf;
use statrs::function::gamma::ln_gamma;
use xz2::read::XzDecoder;

use super::ProfilerConfig;
use crate::taxonomy::Taxonomy;

/// Coverage statistics for a single assembly.
///
/// Tracks multiple coverage metrics used for QC filtering.
#[derive(Clone, Debug, Default)]
struct AssemblyStats {
    /// Total assembly length in bases
    length: usize,

    /// Breadth of coverage: fraction of assembly with ≥1 read
    breadth: f64,

    /// Fixed-size chunk coverage
    chunk_breadth: f64,

    /// Flexible-size chunk coverage (based on √assembly_size)
    flex_chunk_breadth: f64,

    /// Number of flexible chunks used for coverage calculation
    num_flex_chunks: usize,
}

/// Coverage metrics aggregated by taxonomy ID.
///
/// For strain-level profiling, each strain taxid has its own entry.
/// Used for filtering and statistical calculations.
#[derive(Clone, Debug, Default)]
struct TaxidCoverage {
    /// Breadth of coverage
    breadth: f64,

    /// Fixed chunk coverage
    chunk_breadth: f64,

    /// Flexible chunk coverage
    flex_chunk_breadth: f64,

    /// Number of flexible chunks
    num_flex_chunks: usize,

    /// Total assembly length for this taxid
    asm_length: usize,
}

/// Metadata extracted from a reference sequence name.
///
/// Reference names follow: `asm|strain_taxid|species_taxid|acc|genome_size[|sampled]`
#[derive(Clone, Debug, Default)]
struct ReferenceMetadata {
    /// Assembly identifier (first field)
    asm: String,

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

    /// Flexible flank lengths per assembly: asm -> √(asm_length)
    asm_flex_flank_len: Arc<HashMap<String, usize>>,

    /// Fixed flank length (1000 for paired, 2000 for single)
    fixed_flank_len: usize,

    /// Whether input is paired-end data
    is_paired: bool,
}

/// Result from a worker thread processing alignment lines.
struct WorkerResult {
    /// Parsed alignments: (read_name, record)
    alignments: Vec<(String, AlignmentRecord)>,

    /// All read names encountered
    reads_seen: Vec<String>,

    /// Alignment span intervals per reference
    intervals: HashMap<String, Vec<(usize, usize)>>,

    /// Fixed chunk intervals per reference
    chunk_intervals: HashMap<String, Vec<(usize, usize)>>,

    /// Flexible chunk intervals per reference
    flex_intervals: HashMap<String, Vec<(usize, usize)>>,

    /// Count of alignments processed
    alignments_processed: usize,
}

/// Intermediate result from parsing a single SAM line.
struct ParsedAlignment {
    /// Read name (pair suffix stripped if paired-end)
    read_name: String,

    /// Whether to count this in alignment statistics
    count_alignment: bool,

    /// Parsed alignment record (None if unmapped/filtered)
    record: Option<ParsedAlignmentRecord>,
}

/// Full parsed data from a valid alignment.
struct ParsedAlignmentRecord {
    /// Full reference name from SAM
    ref_name: String,

    /// Alignment span coordinates (start, end)
    span: (usize, usize),

    /// Fixed chunk coordinates
    fixed_chunk: (usize, usize),

    /// Flexible chunk coordinates
    flex_chunk: (usize, usize),

    /// Number of aligned bases
    mapped_len: usize,

    /// CIGAR string
    cigar: String,

    /// Assembly identifier
    asm: String,

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
    read_alignments_map: &mut HashMap<String, Vec<AlignmentRecord>>,
    all_reads: &mut HashSet<String>,
    acc_intervals: &mut HashMap<String, Vec<(usize, usize)>>,
    acc_chunk_intervals: &mut HashMap<String, Vec<(usize, usize)>>,
    acc_flex_chunk_intervals: &mut HashMap<String, Vec<(usize, usize)>>,
    num_processed: &mut usize,
    next_progress: &mut usize,
) {
    for (read, aln) in result.alignments {
        read_alignments_map.entry(read).or_default().push(aln);
    }
    for read in result.reads_seen {
        all_reads.insert(read);
    }
    for (ref_name, intervals) in result.intervals {
        let entry = acc_intervals.entry(ref_name).or_default();
        entry.reserve(intervals.len());
        entry.extend(intervals);
    }
    for (ref_name, intervals) in result.chunk_intervals {
        let entry = acc_chunk_intervals.entry(ref_name).or_default();
        entry.reserve(intervals.len());
        entry.extend(intervals);
    }
    for (ref_name, intervals) in result.flex_intervals {
        let entry = acc_flex_chunk_intervals.entry(ref_name).or_default();
        entry.reserve(intervals.len());
        entry.extend(intervals);
    }
    *num_processed += result.alignments_processed;
    while *num_processed >= *next_progress {
        log::info!(
            target: "PROFILE",
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
    asm_stats: HashMap<String, AssemblyStats>,

    /// Per-taxid coverage (propagated from assembly stats)
    taxid_coverage: HashMap<u32, TaxidCoverage>,

    /// Minimum breadth threshold for reporting
    min_breadth: f64,

    /// Minimum chunk breadth threshold (auto-estimated or user-set)
    min_chunk_breadth: f64,

    /// Scaling factors for subsampled genomes
    subsample_scaling: HashMap<u32, f64>,
}

/// Simplified alignment record for taxonomy assignment phase.
#[derive(Clone, Debug)]
struct AlignmentRecord {
    /// Number of aligned bases
    mapped_len: usize,

    /// CIGAR string for identity calculation
    cigar: String,

    /// Assembly identifier
    asm: String,

    /// Primary taxonomy ID (strain-level if available)
    taxid: Option<u32>,

    /// Fallback taxonomy ID (species-level)
    fallback_taxid: Option<u32>,
}

/// Collection of alignments for a single read.
#[derive(Clone, Debug)]
struct ReadAlignment {
    /// Read identifier
    read_name: String,

    /// All alignments for this read
    alignments: Vec<AlignmentRecord>,
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
    read_name: String,

    /// Taxonomic information for each candidate
    taxa: Vec<TaxonInfo>,

    /// Taxonomy IDs (parallel to taxa)
    taxids: Vec<u32>,

    /// Depth contributions (mapped_len / genome_size) per candidate
    depths: Vec<f64>,

    /// Mapped lengths per candidate
    maplens: Vec<f64>,
}

/// Results from the taxonomy assignment phase.
struct AssignTaxaResult {
    /// Per-read taxonomic assignments
    assignments: Vec<ReadAssignment>,

    /// Count of alignments passing quality filters
    alignments_passing: usize,
}

/// Regex for parsing CIGAR strings.
/// Matches operations like 100M, 50=, 10X, 5I, 3D
// static CIGAR_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"[0-9]+[MIDNSHP=X]").unwrap());

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

    if ref_name == "*" || cigar == "*" {
        return Some(ParsedAlignment {
            read_name,
            count_alignment: false,
            record: None,
        });
    }

    if ctx.is_paired && rnext != "=" {
        return Some(ParsedAlignment {
            read_name,
            count_alignment: false,
            record: None,
        });
    }

    let metadata = match ctx.ref_metadata.get(ref_name) {
        Some(meta) => meta,
        None => {
            return Some(ParsedAlignment {
                read_name,
                count_alignment: false,
                record: None,
            });
        }
    };
    if metadata.length == 0 {
        return Some(ParsedAlignment {
            read_name,
            count_alignment: false,
            record: None,
        });
    }

    let mapped_len = read_seq.len();
    if mapped_len == 0 {
        return Some(ParsedAlignment {
            read_name,
            count_alignment: false,
            record: None,
        });
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
    let midpoint = start + ((end.saturating_sub(start)) / 2);
    let flex_len = ctx
        .asm_flex_flank_len
        .get(&metadata.asm)
        .copied()
        .unwrap_or_else(|| ((metadata.length as f64).sqrt() as usize).max(1));
    let fixed_chunk_start = start.saturating_sub(ctx.fixed_flank_len).max(1);
    let fixed_chunk_end = end.saturating_add(ctx.fixed_flank_len).min(ref_end);
    let flex_chunk_start = midpoint.saturating_sub(flex_len).max(1);
    let flex_chunk_end = midpoint.saturating_add(flex_len).min(ref_end);

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
            });
        }
    }

    Some(ParsedAlignment {
        read_name,
        count_alignment,
        record: Some(ParsedAlignmentRecord {
            ref_name: ref_name.to_string(),
            span,
            fixed_chunk: (fixed_chunk_start, fixed_chunk_end),
            flex_chunk: (flex_chunk_start, flex_chunk_end),
            mapped_len,
            cigar: cigar.to_string(),
            asm: metadata.asm.clone(),
            taxid: metadata.taxid,
            fallback_taxid: metadata.fallback_taxid,
        }),
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
    intervals.sort_unstable_by(|a, b| a.0.cmp(&b.0));
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

/// Extract subsampling scaling factor from reference name.
///
/// If reference has genome_size and sampled_size fields, computes
/// the scaling factor for abundance correction.
fn update_subsample_scaling(map: &mut HashMap<u32, f64>, parts: &[&str]) {
    if parts.len() < 6 {
        return;
    }
    let genome_size = match parts[4].parse::<f64>() {
        Ok(value) if value > 0.0 => value,
        _ => return,
    };
    let sampled_size = match parts[5].parse::<f64>() {
        Ok(value) if value > 0.0 => value,
        _ => return,
    };
    if sampled_size >= genome_size {
        return;
    }
    let ratio = genome_size / sampled_size;
    for field in parts.iter().skip(1).take(2) {
        if let Ok(taxid) = field.parse::<u32>() {
            map.entry(taxid).or_insert(ratio);
        }
    }
}

/// Write a line to an optional writer.
fn write_optional_line(writer: &mut Option<BufWriter<File>>, line: String) -> Result<()> {
    if let Some(inner) = writer.as_mut() {
        inner.write_all(line.as_bytes())?;
        inner.write_all(b"\n")?;
    }
    Ok(())
}

/// Calculate alignment identity and fraction from CIGAR string.
///
/// # Returns
///
/// Tuple of (identity, matched_len, fraction)
/// - identity: matches / gap-compressed alignment length
/// - matched_len: number of exact sequence matches (=)
/// - fraction: matched bases / total read length
// fn cigar_to_identity(cigar: &str) -> (f64, usize, f64) {
//     let mut operation_len: HashMap<char, usize> = HashMap::new();
//     let mut uncompressed_len = 0usize;
//     let mut num_gap_events = 0usize;
//     for caps in CIGAR_PATTERN.find_iter(cigar) {
//         let token = caps.as_str();
//         let (len_part, op_part) = token.split_at(token.len() - 1);
//         if let Ok(count) = len_part.parse::<usize>() {
//             let op = op_part.chars().next().unwrap();
//             *operation_len.entry(op).or_insert(0) += count;
//             uncompressed_len += count;
//             if matches!(op, 'M' | '=' | 'X' | 'H') {
//                 continue;
//             }
//             num_gap_events += 1;
//         }
//     }
//     let matched_len = *operation_len.get(&'=').unwrap_or(&0);
//     let fraction = if uncompressed_len == 0 {
//         0.0
//     } else {
//         matched_len as f64 / uncompressed_len as f64
//     };
//     let gap_compressed_alignment_len = matched_len
//         + operation_len.get(&'M').copied().unwrap_or(0)
//         + operation_len.get(&'X').copied().unwrap_or(0)
//         + num_gap_events;
//     let identity = if gap_compressed_alignment_len == 0 {
//         0.0
//     } else {
//         matched_len as f64 / gap_compressed_alignment_len as f64
//     };
//     (identity, matched_len, fraction)
// }

/// Primary implementation with an explicit switch.
///
/// - If `gap_compressed == true`:
///     * extended CIGAR (=,X present): identity = matches / (matches + mismatches + gap_events)
///     * non-extended (only M):        identity proxy = M / (M + gap_events)
///   where gap_events counts contiguous runs of I and D as 1 each.
///
/// - If `gap_compressed == false`:
///     * extended CIGAR (=,X present): identity = matches / (matches + mismatches + I_bases + D_bases)
///     * non-extended (only M):        identity proxy = M / (M + I_bases + D_bases)
///
/// `fraction` (query coverage) is always:
///     fraction = (M + = + X + I) / (M + = + X + I + S)
///
/// Returned `matched_len`:
///     * extended CIGAR: '='
///     * otherwise:      'M' (proxy)
pub fn cigar_to_identity_with_opts(cigar: &str, gap_compressed: bool) -> (f64, usize, f64) {
    const COUNT_N_AS_DELETION: bool = false;

    // Length counters
    let mut eq: usize = 0;
    let mut x: usize = 0;
    let mut i_: usize = 0;
    let mut d: usize = 0;
    let mut s: usize = 0;
    let mut m: usize = 0;

    // Gap-event counters (contiguous runs)
    let mut ins_events: usize = 0;
    let mut del_events: usize = 0;

    // Parsing state
    let mut count_buffer: usize = 0;
    let mut saw_any_op: bool = false;

    // Track whether we're inside a contiguous insertion/deletion run
    let mut in_ins_run: bool = false;
    let mut in_del_run: bool = false;

    for &byte in cigar.as_bytes() {
        if byte.is_ascii_digit() {
            count_buffer = count_buffer
                .saturating_mul(10)
                .saturating_add((byte - b'0') as usize);
            continue;
        }

        if count_buffer == 0 {
            // malformed token (e.g. "M") or "0M"; skip
            continue;
        }

        saw_any_op = true;

        match byte {
            b'=' => {
                eq += count_buffer;
                in_ins_run = false;
                in_del_run = false;
            }
            b'X' => {
                x += count_buffer;
                in_ins_run = false;
                in_del_run = false;
            }
            b'M' => {
                m += count_buffer;
                in_ins_run = false;
                in_del_run = false;
            }
            b'I' => {
                i_ += count_buffer;
                if !in_ins_run {
                    ins_events += 1;
                    in_ins_run = true;
                }
                in_del_run = false;
            }
            b'D' => {
                d += count_buffer;
                if !in_del_run {
                    del_events += 1;
                    in_del_run = true;
                }
                in_ins_run = false;
            }
            b'S' => {
                s += count_buffer;
                in_ins_run = false;
                in_del_run = false;
            }
            b'N' if COUNT_N_AS_DELETION => {
                d += count_buffer;
                if !in_del_run {
                    del_events += 1;
                    in_del_run = true;
                }
                in_ins_run = false;
            }
            _ => {
                // Ignore H, P, N (if not counted), etc.
                // Break gap runs because these terminate/interrupt the aligned segment.
                in_ins_run = false;
                in_del_run = false;
            }
        }

        count_buffer = 0;
    }

    if !saw_any_op {
        return (0.0, 0, 0.0);
    }

    // Fraction (query coverage)
    let aligned_query = m + eq + x + i_;
    let query_len = aligned_query + s;
    let fraction = if query_len == 0 {
        0.0
    } else {
        aligned_query as f64 / query_len as f64
    };

    // Identity
    let has_extended = (eq + x) > 0;
    let matched_len: usize = if has_extended { eq } else { m };

    let identity_den: usize = if gap_compressed {
        let gap_events = ins_events + del_events;
        if has_extended {
            eq + x + gap_events
        } else {
            m + gap_events
        }
    } else {
        if has_extended {
            eq + x + i_ + d
        } else {
            m + i_ + d
        }
    };

    let identity = if identity_den == 0 {
        0.0
    } else {
        matched_len as f64 / identity_den as f64
    };

    (identity, matched_len, fraction)
}

/// Backwards-compatible wrapper:
/// Call exactly like before: `cigar_to_identity(cigar)`
/// and it will compute GAP-COMPRESSED identity by default.
pub fn cigar_to_identity(cigar: &str) -> (f64, usize, f64) {
    cigar_to_identity_with_opts(cigar, true)
}

#[cfg(test)]
mod tests {
    use super::{calc_breadth_pvalue, cigar_to_identity, Tail};

    fn approx_equal(lhs: f64, rhs: f64) {
        assert!(
            (lhs - rhs).abs() < 1e-9,
            "left {} not approximately equal to right {}",
            lhs,
            rhs
        );
    }

    #[test]
    fn cigar_identity_matches() {
        // Matches the semantics of `cigar_to_identity_with_opts(_, true)`:
        //   identity = matches / (matches + mismatches + gap_events)
        //   fraction = aligned_query / (aligned_query + S)  [H is excluded]
        //   matched_len = `eq` for extended CIGAR, falls back to `m`.
        let (identity, matched_len, fraction) = cigar_to_identity("80=5I5=10X5H");
        approx_equal(identity, 85.0 / 96.0);
        assert_eq!(matched_len, 85);
        approx_equal(fraction, 1.0);

        let (identity, matched_len, fraction) = cigar_to_identity("50=50X");
        approx_equal(identity, 0.5);
        assert_eq!(matched_len, 50);
        approx_equal(fraction, 1.0);

        // Non-extended CIGAR: identity proxy = M / (M + gap_events), and
        // matched_len reports M when no '=' is present.
        let (identity, matched_len, fraction) = cigar_to_identity("100M");
        approx_equal(identity, 1.0);
        assert_eq!(matched_len, 100);
        approx_equal(fraction, 1.0);
    }

    #[test]
    fn expected_breadth_is_clamped_for_high_depths() {
        let mapcount: f64 = 5.0;
        let depth: f64 = 10.0;
        let per_read_cov = (depth / mapcount).clamp(0.0, 1.0);
        let expected_breadth = 1.0 - (1.0 - per_read_cov).powf(mapcount);

        approx_equal(expected_breadth, 1.0);

        let pvalue = calc_breadth_pvalue(
            mapcount.round() as usize,
            Some(depth),
            0.8,
            Some(expected_breadth),
            Tail::TwoSided,
        );

        let pv = pvalue.expect("expected numeric p-value even when expected breadth is ~1");
        assert!(pv.is_finite());
        assert!((0.0..=1.0).contains(&pv));
    }
}

/// Round to 5 decimal places.
fn round5(value: f64) -> f64 {
    (value * 100000.0).round() / 100000.0
}

/// Normalize fraction values in a map to sum to 1.0.
///
/// Used for EM assignment fractions.
fn normalize_fraction_map(map: &mut IndexMap<u32, (String, String, f64)>) {
    let sum: f64 = map.values().map(|(_, _, frac)| *frac).sum();
    if sum <= f64::EPSILON {
        return;
    }
    for (_, _, frac) in map.values_mut() {
        *frac = round5(*frac / sum);
    }
}

/// Format a p-value for output.
///
/// Very small values use scientific notation; None becomes "NA".
fn format_probability_field(value: Option<f64>) -> String {
    match value {
        Some(v) if v.is_finite() => {
            let adjusted = if v.abs() < 1e-12 { 0.0 } else { v };
            if adjusted.abs() > 0.0 && adjusted.abs() < 1e-4 {
                format!("{:.6e}", adjusted)
            } else {
                format!("{:.5}", adjusted)
            }
        }
        _ => "NA".to_string(),
    }
}

// ============================================================================
// Statistical Functions for Coverage P-value Calculation
// ============================================================================
//
// These functions assess whether observed coverage patterns are consistent
// with random read placement. Significant deviations may indicate false
// positives from repetitive sequences or misalignments.

/// Clamp value to [0, 1] range.
fn clamp01(x: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    }
    if x > 1.0 {
        return 1.0;
    }
    x
}

/// Check if value is a finite number.
fn is_finite_number(x: f64) -> bool {
    x.is_finite()
}

/// Infer genome size from read count and depth.
fn infer_n_from_depth(read_count: usize, depth: f64) -> Option<usize> {
    if read_count == 0 || depth <= 0.0 || !is_finite_number(depth) {
        return None;
    }
    let n = (read_count as f64 / depth).round() as isize;
    Some(n.max(1) as usize)
}

/// Infer genome size from read count and expected breadth.
fn infer_n_from_exp_breadth(read_count: usize, exp_breadth: f64) -> Option<usize> {
    if read_count == 0 {
        return None;
    }
    let p = clamp01(exp_breadth);
    if p <= 0.0 || p >= 1.0 {
        return Some(read_count.max(1));
    }
    let denom = (1.0 - p).ln();
    if !denom.is_finite() || denom == 0.0 {
        return None;
    }
    let n = (-(read_count as f64) / denom).round() as isize;
    Some(n.max(1) as usize)
}

/// Convert breadth fraction to discrete coverage count.
fn k_from_breadth(breadth: f64, n: usize) -> Option<usize> {
    if !is_finite_number(breadth) || breadth < 0.0 || breadth > 1.0 {
        return None;
    }
    let mut k = (breadth * n as f64).round() as isize;
    if k < 0 {
        k = 0;
    }
    if k > n as isize {
        k = n as isize;
    }
    Some(k as usize)
}

/// Compute Stirling numbers of the second kind for exact p-value calculation.
fn stirling2_row(m: usize, kmax: usize) -> Vec<f64> {
    let mut s_prev = vec![0.0_f64; kmax + 1];
    s_prev[0] = 1.0;

    for i in 1..=m {
        let mut s_cur = vec![0.0_f64; kmax + 1];
        let upper = i.min(kmax);
        for k in 1..=upper {
            s_cur[k] = (s_prev[k] * k as f64) + s_prev[k - 1];
        }
        s_prev = s_cur;
    }

    s_prev
}

/// Compute log(k!) using gamma function.
fn log_factorial(k: usize) -> f64 {
    ln_gamma(k as f64 + 1.0)
}

/// Compute log of binomial coefficient C(n, k).
fn log_combination(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

/// Numerically stable log-sum-exp computation.
fn logsumexp(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }
    let m = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !m.is_finite() {
        return m;
    }
    let mut acc = 0.0;
    for v in values {
        acc += (v - m).exp();
    }
    m + acc.ln()
}

/// Which tail of the distribution to test.
#[allow(dead_code)]
#[derive(Clone, Copy)]
enum Tail {
    Lower,
    Upper,
    TwoSided,
}

/// Exact p-value using Stirling numbers (for small samples).
fn pvalue_exact_stirling(
    read_count: usize,
    depth: Option<f64>,
    breadth: f64,
    exp_breadth: f64,
    side: Tail,
    max_n: usize,
    max_m: usize,
) -> Option<f64> {
    if read_count == 0 {
        return None;
    }

    let p = clamp01(exp_breadth);
    if p <= 0.0 || p >= 1.0 {
        return None;
    }

    let mut n = depth
        .filter(|d| *d > 0.0)
        .and_then(|d| infer_n_from_depth(read_count, d));
    if n.is_none() {
        n = infer_n_from_exp_breadth(read_count, p);
    }
    let n = n?;

    if n > max_n || read_count > max_m {
        return None;
    }

    let mut k_obs = k_from_breadth(breadth, n)?;
    let k_max_possible = n.min(read_count);
    if k_obs > k_max_possible {
        k_obs = k_max_possible;
    }

    if k_obs == 0 {
        return match side {
            Tail::Lower => Some(0.0),
            Tail::Upper => Some(1.0),
            Tail::TwoSided => Some(0.0),
        };
    }

    let s_values = stirling2_row(read_count, k_max_possible);
    let log_den = (read_count as f64) * (n as f64).ln();

    let pmf_log = |k: usize| -> f64 {
        if k == 0 || k > k_max_possible {
            return f64::NEG_INFINITY;
        }
        let log_comb = log_combination(n, k);
        let log_stir = s_values[k].ln();
        let log_fact = log_factorial(k);
        log_comb + log_stir + log_fact - log_den
    };

    let cdf_log = |k: usize| -> f64 {
        if k == 0 {
            return f64::NEG_INFINITY;
        }
        if k >= k_max_possible {
            return 0.0;
        }
        let mut logs = Vec::with_capacity(k);
        for t in 1..=k {
            logs.push(pmf_log(t));
        }
        logsumexp(&logs)
    };

    let value_log = match side {
        Tail::Lower => cdf_log(k_obs),
        Tail::Upper => {
            let lower_log = if k_obs > 1 {
                cdf_log(k_obs - 1)
            } else {
                f64::NEG_INFINITY
            };
            let lower_prob = if lower_log.is_finite() {
                lower_log.exp()
            } else {
                0.0
            };
            clamp01(1.0 - lower_prob).ln()
        }
        Tail::TwoSided => {
            let low = cdf_log(k_obs);
            let up = {
                let lower_log = if k_obs > 1 {
                    cdf_log(k_obs - 1)
                } else {
                    f64::NEG_INFINITY
                };
                let lower_prob = if lower_log.is_finite() {
                    lower_log.exp()
                } else {
                    0.0
                };
                clamp01(1.0 - lower_prob).ln()
            };
            (2.0_f64 * low.exp().min(up.exp())).ln()
        }
    };

    Some(clamp01(value_log.exp()))
}

/// P-value using binomial distribution approximation.
fn pvalue_binomial_approx(
    read_count: usize,
    depth: Option<f64>,
    breadth: f64,
    exp_breadth: f64,
    side: Tail,
) -> Option<f64> {
    if read_count == 0 {
        return None;
    }

    let p = clamp01(exp_breadth);
    if p <= 0.0 || p >= 1.0 {
        return None;
    }

    let mut n = depth
        .filter(|d| *d > 0.0)
        .and_then(|d| infer_n_from_depth(read_count, d));
    if n.is_none() {
        n = infer_n_from_exp_breadth(read_count, p);
    }
    let n = n?;

    let k_obs = k_from_breadth(breadth, n)?;

    let cdf = |k: isize| -> f64 {
        if k < 0 {
            return 0.0;
        }
        if k as usize >= n {
            return 1.0;
        }
        beta_reg((n - k as usize) as f64, (k as usize + 1) as f64, 1.0 - p)
    };

    let val = match side {
        Tail::Lower => cdf(k_obs as isize),
        Tail::Upper => 1.0 - cdf(k_obs as isize - 1),
        Tail::TwoSided => {
            let low = cdf(k_obs as isize);
            let up = 1.0 - cdf(k_obs as isize - 1);
            2.0 * low.min(up)
        }
    };

    Some(clamp01(val))
}

/// P-value using normal approximation (for large samples).
fn pvalue_normal_approx(
    read_count: usize,
    depth: Option<f64>,
    breadth: f64,
    exp_breadth: f64,
    side: Tail,
) -> Option<f64> {
    if read_count == 0 {
        return None;
    }

    let p = clamp01(exp_breadth);
    if p <= 0.0 || p >= 1.0 {
        return None;
    }

    let mut n = depth
        .filter(|d| *d > 0.0)
        .and_then(|d| infer_n_from_depth(read_count, d));
    if n.is_none() {
        n = infer_n_from_exp_breadth(read_count, p);
    }
    let n = n?;

    let k_obs = k_from_breadth(breadth, n)? as f64;

    let mu = n as f64 * p;
    let var = n as f64 * p * (1.0 - p);
    if var <= 0.0 {
        return None;
    }
    let sd = var.sqrt();

    let norm_cdf = |x: f64| 0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2));

    let val = match side {
        Tail::Lower => {
            let z = (k_obs + 0.5 - mu) / sd;
            norm_cdf(z)
        }
        Tail::Upper => {
            let z = (k_obs - 0.5 - mu) / sd;
            1.0 - norm_cdf(z)
        }
        Tail::TwoSided => {
            let z_low = (k_obs + 0.5 - mu) / sd;
            let z_up = (k_obs - 0.5 - mu) / sd;
            let p_low = norm_cdf(z_low);
            let p_up = 1.0 - norm_cdf(z_up);
            2.0 * p_low.min(p_up)
        }
    };

    Some(clamp01(val))
}

/// Calculate p-value for observed coverage breadth.
///
/// Selects appropriate method based on sample size:
/// exact (small), normal (moderate), or binomial (fallback).
fn calc_breadth_pvalue(
    read_count: usize,
    depth: Option<f64>,
    breadth: f64,
    exp_breadth: Option<f64>,
    side: Tail,
) -> Option<f64> {
    if read_count == 0 || !is_finite_number(breadth) || breadth < 0.0 || breadth > 1.0 {
        return None;
    }

    let mut exp_breadth = exp_breadth;
    if exp_breadth.is_none() {
        exp_breadth = depth.filter(|d| *d > 0.0).map(|d| 1.0 - (-d).exp());
    }
    let exp_breadth = exp_breadth?;
    let p = clamp01(if exp_breadth >= 1.0 {
        1.0 - 1e-12
    } else {
        exp_breadth
    });
    if p <= 0.0 || p >= 1.0 {
        return None;
    }

    let mut n = depth
        .filter(|d| *d > 0.0)
        .and_then(|d| infer_n_from_depth(read_count, d));
    if n.is_none() {
        n = infer_n_from_exp_breadth(read_count, p);
    }
    let n = n?;
    let eff = (n as f64 * p).min(n as f64 * (1.0 - p));

    if n <= 3000 && read_count <= 100 {
        if let Some(pv) = pvalue_exact_stirling(read_count, depth, breadth, p, side, 3000, 100) {
            return Some(pv);
        }
    }

    if eff >= 20.0 && (0.05..=0.95).contains(&p) {
        if let Some(pv) = pvalue_normal_approx(read_count, depth, breadth, p, side) {
            return Some(pv);
        }
    }

    pvalue_binomial_approx(read_count, depth, breadth, p, side)
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

/// Pathogen metadata from host mapping file.
struct PathogenEntry {
    /// Semicolon-separated host taxonomy IDs
    host_taxids: String,

    /// Semicolon-separated host names
    host_names: String,

    /// Semicolon-separated associated diseases
    diseases: String,
}

/// Load pathogen-host mapping from TSV file.
fn load_pathogen_table(path: &str) -> Result<HashMap<u32, PathogenEntry>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open pathogen host table: {}", path))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for (idx, line) in reader.lines().enumerate() {
        let line = line?;
        if idx == 0 {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }
        let pathogen_taxid: u32 = parts[0].parse()?;
        map.insert(
            pathogen_taxid,
            PathogenEntry {
                host_taxids: parts[1].to_string(),
                host_names: parts[2].to_string(),
                diseases: parts[3].to_string(),
            },
        );
    }
    Ok(map)
}

/// Check if a taxonomy ID is in the lineage of any given host taxids.
///
/// Used for pathogen filtering - checks if the candidate host is
/// an ancestor of any known host for the pathogen.
fn lineage_contains(taxonomy: &Taxonomy, candidate: u32, host_taxids: &[u32]) -> bool {
    if host_taxids.iter().any(|taxid| *taxid == candidate) {
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
        log::info!(target: "PROFILE", "Loading alignments ...");
        let sam_path = Path::new(&self.cfg.sam);
        if !sam_path.exists() {
            bail!("SAM file {} does not exist", sam_path.display());
        }
        let verbose_logging = self.cfg.verbose || self.cfg.very_verbose;
        let SamLoadResult {
            read_alignments,
            asm_stats,
            taxid_coverage,
            min_breadth,
            min_chunk_breadth,
            subsample_scaling,
        } = self
            .load_sam(sam_path)
            .context("failed to parse SAM file")?;
        log::info!(
            target: "PROFILE",
            "Loaded {} reads with alignments across {} reference sequences",
            read_alignments.len(),
            asm_stats.len()
        );
        if verbose_logging {
            log::info!(
                target: "PROFILE",
                "Assigning taxonomy for {} reads with {} threads ...",
                read_alignments.len(),
                self.cfg.threads
            );
        }
        let AssignTaxaResult {
            assignments,
            alignments_passing,
        } = self.assign_taxa(&read_alignments, &asm_stats, &taxid_coverage, min_breadth)?;
        let assigned_read_count = assignments.len();
        if self.cfg.very_verbose {
            let avg_alignments_per_read = if assigned_read_count > 0 {
                alignments_passing as f64 / assigned_read_count as f64
            } else {
                0.0
            };
            log::info!(
                target: "PROFILE",
                "Alignments passing identity/fraction filters: {} (across {} reads, avg {:.2})",
                alignments_passing,
                assigned_read_count,
                avg_alignments_per_read
            );
        }
        log::info!(
            target: "PROFILE",
            "Assigned taxonomy to {} reads",
            assigned_read_count
        );
        let skip_classification = !subsample_scaling.is_empty();
        if verbose_logging {
            log::info!(target: "PROFILE", "Profiling taxonomy abundance ...");
        }
        let final_classified = self.write_outputs(
            assignments,
            &taxid_coverage,
            min_chunk_breadth,
            &subsample_scaling,
            skip_classification,
        )?;
        if final_classified != assigned_read_count {
            log::warn!(
                target: "PROFILE",
                "Classified read total ({}) differed from assignment count ({})",
                final_classified,
                assigned_read_count
            );
        }
        log::info!(
            target: "PROFILE",
            "Number of taxonomically classified reads: {}",
            final_classified
        );
        log::info!(target: "PROFILE", "Taxonomy profiling finished.");
        Ok(())
    }

    /// Load and parse SAM file, computing coverage metrics.
    ///
    /// Unlike community profiler, this tracks coverage per assembly
    /// for strain-level resolution.
    fn load_sam(&self, sam_path: &Path) -> Result<SamLoadResult> {
        let reader = open_alignment_reader(sam_path)?;
        let mut reader = BufReader::new(reader);
        let mut asm_stats: HashMap<String, AssemblyStats> = HashMap::new();
        let mut ref_metadata: HashMap<String, ReferenceMetadata> = HashMap::new();
        let mut asm_to_taxids: HashMap<String, HashSet<u32>> = HashMap::new();
        let mut subsample_scaling: HashMap<u32, f64> = HashMap::new();
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
                    let mut sn: Option<String> = None;
                    let mut ln: Option<usize> = None;
                    for part in trimmed.split('\t').skip(1) {
                        if part.starts_with("SN:") {
                            sn = Some(part[3..].to_string());
                        } else if part.starts_with("LN:") {
                            ln = Some(part[3..].parse()?);
                        }
                    }
                    if let (Some(sn), Some(len)) = (sn, ln) {
                        let parts: Vec<&str> = sn.split('|').collect();
                        let (asm_name, taxid, fallback_taxid) = parse_reference_header(&parts);
                        update_subsample_scaling(&mut subsample_scaling, &parts);
                        let metadata = ReferenceMetadata {
                            asm: asm_name.clone(),
                            length: len,
                            taxid,
                            fallback_taxid,
                        };
                        ref_metadata.insert(sn.clone(), metadata);
                        let stats = asm_stats.entry(asm_name.clone()).or_default();
                        stats.length += len;
                        let taxid_set = asm_to_taxids.entry(asm_name).or_default();
                        if let Some(taxid) = taxid {
                            taxid_set.insert(taxid);
                        }
                        if let Some(fallback) = fallback_taxid {
                            taxid_set.insert(fallback);
                        }
                    }
                }
                continue;
            }
            first_alignment = Some(trimmed.to_string());
            break;
        }

        let fixed_flank_len = if self.cfg.is_paired { 1000 } else { 2000 };
        let asm_flex_flank_len: HashMap<String, usize> = asm_stats
            .iter()
            .map(|(asm, stats)| {
                let length = stats.length.max(1);
                let flank = ((length as f64).sqrt() as usize).max(1);
                (asm.clone(), flank)
            })
            .collect();

        let worker_context = Arc::new(WorkerContext {
            ref_metadata: Arc::new(ref_metadata.clone()),
            asm_flex_flank_len: Arc::new(asm_flex_flank_len),
            fixed_flank_len,
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
                    let map_capacity = (chunk_len / 4).max(1);
                    let mut intervals: HashMap<String, Vec<(usize, usize)>> =
                        HashMap::with_capacity(map_capacity);
                    let mut chunk_intervals: HashMap<String, Vec<(usize, usize)>> =
                        HashMap::with_capacity(map_capacity);
                    let mut flex_intervals: HashMap<String, Vec<(usize, usize)>> =
                        HashMap::with_capacity(map_capacity);
                    let mut alignments_processed = 0usize;
                    for line in chunk {
                        if let Some(parsed) = parse_alignment_line(&line, &ctx) {
                            let ParsedAlignment {
                                read_name,
                                count_alignment,
                                record,
                            } = parsed;
                            reads_seen.push(read_name.clone());
                            if count_alignment {
                                alignments_processed += 1;
                            }
                            if let Some(record) = record {
                                let ref_name = record.ref_name.clone();
                                intervals
                                    .entry(ref_name.clone())
                                    .or_default()
                                    .push(record.span);
                                chunk_intervals
                                    .entry(ref_name.clone())
                                    .or_default()
                                    .push(record.fixed_chunk);
                                flex_intervals
                                    .entry(ref_name)
                                    .or_default()
                                    .push(record.flex_chunk);
                                alignments.push((
                                    read_name,
                                    AlignmentRecord {
                                        mapped_len: record.mapped_len,
                                        cigar: record.cigar,
                                        asm: record.asm,
                                        taxid: record.taxid,
                                        fallback_taxid: record.fallback_taxid,
                                    },
                                ));
                            }
                        }
                    }
                    if !reads_seen.is_empty()
                        || !alignments.is_empty()
                        || alignments_processed > 0
                        || !intervals.is_empty()
                    {
                        let _ = tx.send(WorkerResult {
                            alignments,
                            reads_seen,
                            intervals,
                            chunk_intervals,
                            flex_intervals,
                            alignments_processed,
                        });
                    }
                }
            }));
        }
        drop(result_tx);
        drop(line_rx);

        let mut acc_intervals: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
        let mut acc_chunk_intervals: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
        let mut acc_flex_chunk_intervals: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
        let mut read_alignments_map: HashMap<String, Vec<AlignmentRecord>> = HashMap::new();
        let mut all_reads: HashSet<String> = HashSet::new();
        let mut num_alignment = 0usize;
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
                        &mut acc_intervals,
                        &mut acc_chunk_intervals,
                        &mut acc_flex_chunk_intervals,
                        &mut num_alignment,
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
                &mut acc_intervals,
                &mut acc_chunk_intervals,
                &mut acc_flex_chunk_intervals,
                &mut num_alignment,
                &mut next_progress,
            );
        }

        for result in result_rx.iter() {
            absorb_worker_result(
                result,
                &mut read_alignments_map,
                &mut all_reads,
                &mut acc_intervals,
                &mut acc_chunk_intervals,
                &mut acc_flex_chunk_intervals,
                &mut num_alignment,
                &mut next_progress,
            );
        }

        for handle in handles {
            let _ = handle.join();
        }

        log::info!(
            target: "PROFILE",
            "{} alignments processed",
            num_alignment
        );

        let num_reads = all_reads.len();
        let aligned_reads = read_alignments_map.len();
        let mapping_rate = if num_reads == 0 {
            0.0
        } else {
            aligned_reads as f64 / num_reads as f64
        };
        // Honour explicit --by-aligned; otherwise auto-force aligned-read
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
        let mut min_chunk_breadth = self.cfg.chunk_breadth;
        let skip_auto_chunk_filter = self.cfg.min_reads.is_some() && min_chunk_breadth.is_none();
        if min_chunk_breadth.is_none() && !skip_auto_chunk_filter {
            min_chunk_breadth = Some(self.estimate_chunk_breadth(chunk_basis));
        }
        let min_chunk_breadth_value = min_chunk_breadth.unwrap_or(0.0);
        let min_chunk_display = format!("{:.6}", min_chunk_breadth_value);
        const INDENT: &str = "                                                ";
        if self.cfg.verbose || self.cfg.very_verbose {
            log::info!(
                target: "PROFILE",
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

        for (ref_name, intervals) in acc_intervals {
            let metadata = match ref_metadata.get(&ref_name) {
                Some(meta) => meta,
                None => continue,
            };
            let merged = merge_intervals(intervals);
            let chunk_values = acc_chunk_intervals.remove(&ref_name).unwrap_or_default();
            let merged_chunk = merge_intervals(chunk_values);
            let flex_values = acc_flex_chunk_intervals
                .remove(&ref_name)
                .unwrap_or_default();
            let merged_flex = merge_intervals(flex_values);
            let asm_stats_entry = match asm_stats.get_mut(&metadata.asm) {
                Some(stats) => stats,
                None => continue,
            };
            let asm_len = asm_stats_entry.length.max(1) as f64;
            let total_span: usize = merged
                .iter()
                .map(|(start, end)| end.saturating_sub(*start))
                .sum();
            let total_chunk_span: usize = merged_chunk
                .iter()
                .map(|(start, end)| end.saturating_sub(*start))
                .sum();
            let total_flex_span: usize = merged_flex
                .iter()
                .map(|(start, end)| end.saturating_sub(*start))
                .sum();
            asm_stats_entry.breadth += total_span as f64 / asm_len;
            asm_stats_entry.chunk_breadth += total_chunk_span as f64 / asm_len;
            asm_stats_entry.flex_chunk_breadth += total_flex_span as f64 / asm_len;
            let estimated = ((asm_stats_entry.length as f64).sqrt() / 2.0)
                .max(1.0)
                .floor() as usize;
            let flex_chunks = merged_flex.len().max(1);
            asm_stats_entry.num_flex_chunks = asm_stats_entry
                .num_flex_chunks
                .max(estimated.max(flex_chunks));
        }

        for stats in asm_stats.values_mut() {
            stats.breadth = stats.breadth.min(1.0);
            stats.chunk_breadth = stats.chunk_breadth.min(1.0);
            stats.flex_chunk_breadth = stats.flex_chunk_breadth.min(1.0);
            let estimated = ((stats.length as f64).sqrt() / 2.0).max(1.0).floor() as usize;
            stats.num_flex_chunks = stats.num_flex_chunks.max(estimated.max(1));
        }

        let mut taxid_coverage: HashMap<u32, TaxidCoverage> = HashMap::new();
        for (asm, taxids) in &asm_to_taxids {
            if let Some(stats) = asm_stats.get(asm) {
                for taxid in taxids {
                    let entry = taxid_coverage.entry(*taxid).or_default();
                    entry.breadth = entry.breadth.max(stats.breadth);
                    entry.chunk_breadth = entry.chunk_breadth.max(stats.chunk_breadth);
                    entry.flex_chunk_breadth =
                        entry.flex_chunk_breadth.max(stats.flex_chunk_breadth);
                    entry.num_flex_chunks = entry.num_flex_chunks.max(stats.num_flex_chunks.max(1));
                    entry.asm_length = entry.asm_length.max(stats.length);
                }
            }
        }

        let read_alignments: Vec<ReadAlignment> = read_alignments_map
            .into_iter()
            .map(|(read_name, alignments)| ReadAlignment {
                read_name,
                alignments,
            })
            .collect();
        let total_retained_alignments: usize = read_alignments
            .iter()
            .map(|read| read.alignments.len())
            .sum();

        log::info!(target: "PROFILE", "Total number of reads: {}", num_reads);
        log::info!(
            target: "PROFILE",
            "Number of aligned reads: {}",
            read_alignments.len()
        );
        log::info!(
            target: "PROFILE",
            "Number of alignments: {}",
            num_alignment
        );
        log::info!(
            target: "PROFILE",
            "Number of alignments retained after parsing: {}",
            total_retained_alignments
        );
        if skip_auto_chunk_filter {
            log::info!(
                target: "PROFILE",
                "Automatic minimum chunk breadth filter disabled because --min-reads is set."
            );
        } else if by_aligned_effective {
            if by_aligned_auto {
                log::info!(
                    target: "PROFILE",
                    "Minimum chunk breadth estimated from aligned reads ({}) [auto-forced: low-map Illumina run, mapping_rate {:.5} < {:.5}].",
                    aligned_reads,
                    mapping_rate,
                    self.cfg.low_map_rate_threshold
                );
            } else {
                log::info!(
                    target: "PROFILE",
                    "Minimum chunk breadth estimated from aligned reads ({}).",
                    aligned_reads
                );
            }
        } else {
            log::info!(
                target: "PROFILE",
                "Minimum chunk breadth estimated from total reads ({}).",
                num_reads
            );
        }

        Ok(SamLoadResult {
            read_alignments,
            asm_stats,
            taxid_coverage,
            min_breadth,
            min_chunk_breadth: min_chunk_breadth_value,
            subsample_scaling,
        })
    }

    /// Estimate minimum chunk breadth threshold from read count.
    ///
    /// Same formula as community profiler - scales with sequencing depth.
    fn estimate_chunk_breadth(&self, num_reads: usize) -> f64 {
        if self.cfg.lowbiomass {
            return 0.0;
        }
        if self.cfg.sequencer.eq_ignore_ascii_case("illumina") {
            if num_reads <= 100_000 {
                return num_reads as f64 / 1_000_000.0;
            }
            if num_reads <= 1_000_000 {
                let reads_in_million = num_reads as f64 / 1_000_000.0;
                return 0.11 * reads_in_million + 0.09;
            }
            let reads_in_million = num_reads as f64 / 1_000_000.0;
            let value = 0.2 + 0.382 * reads_in_million.log10();
            return value.min(0.95);
        }
        if num_reads <= 5_000 {
            0.0
        } else if num_reads <= 1_000_000 {
            0.3
        } else {
            0.5
        }
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
        asm_stats: &HashMap<String, AssemblyStats>,
        taxid_coverage: &HashMap<u32, TaxidCoverage>,
        min_breadth: f64,
    ) -> Result<AssignTaxaResult> {
        let taxonomy = self.taxonomy.clone();
        let asm_stats = Arc::new(asm_stats.clone());
        let taxid_coverage = Arc::new(taxid_coverage.clone());
        let min_identity = self.cfg.identity;
        let min_mapped_len = self.cfg.mapped_len;
        let min_fraction = self.cfg.fraction;
        let batch_size = std::cmp::max(1, self.cfg.batch_size);
        let pool = ThreadPoolBuilder::new()
            .num_threads(self.cfg.threads.max(1))
            .build()?;
        let pb = ProgressBar::new(read_alignments.len() as u64);
        pb.set_draw_target(ProgressDrawTarget::stdout());
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} {msg} {bar:40.cyan/blue} {pos}/{len} [{elapsed_precise}<{eta_precise}]",
            )?
            .progress_chars("#>-")
        );
        pb.set_message("Assigning taxonomy");
        let progress = pb.clone();
        let passing_counter = Arc::new(AtomicUsize::new(0));
        let assignment_chunks: Vec<Vec<ReadAssignment>> = pool.install(|| {
            read_alignments
                .par_chunks(batch_size)
                .map_with(
                    (progress.clone(), Arc::clone(&passing_counter)),
                    |(progress_bar, counter), chunk| {
                        let mut local: Vec<ReadAssignment> = Vec::new();
                        for read in chunk {
                            let mut taxa_depths: IndexMap<u32, (TaxonInfo, f64, f64)> =
                                IndexMap::new();
                            for aln in &read.alignments {
                                if aln.mapped_len < min_mapped_len {
                                    continue;
                                }
                                let (identity, _, fraction) = cigar_to_identity(&aln.cigar);
                                if identity < min_identity || fraction < min_fraction {
                                    continue;
                                }
                                let stats = match asm_stats.get(&aln.asm) {
                                    Some(stats) => stats,
                                    None => continue,
                                };
                                if stats.length == 0 {
                                    continue;
                                }
                                counter.fetch_add(1, AtomicOrdering::Relaxed);
                                let resolved = resolve_taxon_info(
                                    &taxonomy,
                                    &taxid_coverage,
                                    aln.taxid,
                                    aln.fallback_taxid,
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
                                        }
                                    })
                                    .or_insert((taxon, depth, maplen));
                            }
                            if taxa_depths.is_empty() {
                                continue;
                            }
                            let mut taxids: Vec<u32> = Vec::new();
                            let mut taxa: Vec<TaxonInfo> = Vec::new();
                            let mut depths: Vec<f64> = Vec::new();
                            let mut maplens: Vec<f64> = Vec::new();
                            for (taxid, (taxon, depth, maplen)) in taxa_depths.into_iter() {
                                taxids.push(taxid);
                                taxa.push(taxon);
                                depths.push(depth);
                                maplens.push(maplen);
                            }
                            local.push(ReadAssignment {
                                read_name: read.read_name.clone(),
                                taxa,
                                taxids,
                                depths,
                                maplens,
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
    /// 6. Write profile.txt and classify.txt outputs
    ///
    /// # Arguments
    /// * `assignments` - Per-read taxonomic assignments
    /// * `taxid_coverage` - Coverage metrics per taxid
    /// * `min_chunk_breadth` - Minimum chunk breadth threshold
    /// * `subsample_scaling` - Scaling factors for subsampled genomes
    /// * `skip_classification` - Skip classify.txt for subsampled databases
    fn write_outputs(
        &self,
        assignments: Vec<ReadAssignment>,
        taxid_coverage: &HashMap<u32, TaxidCoverage>,
        min_chunk_breadth: f64,
        subsample_scaling: &HashMap<u32, f64>,
        skip_classification: bool,
    ) -> Result<usize> {
        let classification_path = format!("{}.classify.txt", self.cfg.outprefix);
        let out_prefix = if self.cfg.host.is_some() {
            format!("{}.pathogen", self.cfg.outprefix)
        } else {
            self.cfg.outprefix.clone()
        };
        let profile_path = format!("{}.profile.txt", out_prefix);
        let mut classification_writer = if skip_classification {
            None
        } else {
            Some(BufWriter::new(File::create(&classification_path)?))
        };
        let mut profile_writer = WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(File::create(&profile_path)?);
        let mut raw_writer = if self.cfg.keep_raw {
            Some(
                WriterBuilder::new()
                    .has_headers(false)
                    .delimiter(b'\t')
                    .from_writer(File::create(&format!("{}.rprofile.txt", out_prefix))?),
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
        let mut unambiguous: HashMap<u32, (String, String, f64, f64, f64)> = HashMap::new();
        let mut ambiguous_counts: HashMap<Vec<u32>, (Vec<TaxonInfo>, f64, Vec<f64>, Vec<f64>)> =
            HashMap::new();
        let mut ambiguous_reads: HashMap<Vec<u32>, Vec<String>> = HashMap::new();
        let mut classified_reads = 0usize;
        if skip_classification {
            log::info!(
                target: "PROFILE",
                "Skipping classify.txt generation for subsampled database."
            );
        } else {
            log::info!(
                target: "PROFILE",
                "Writing strain-level classification to {}",
                classification_path
            );
        }
        log::info!(
            target: "PROFILE",
            "Writing abundance profile to {}",
            profile_path
        );
        for assignment in &assignments {
            if assignment.taxids.len() == 1 {
                let taxon = &assignment.taxa[0];
                let taxid = assignment.taxids[0];
                let depth = assignment.depths[0];
                let maplen = assignment.maplens[0];
                classified_reads += 1;
                write_optional_line(
                    &mut classification_writer,
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        assignment.read_name, taxon.name, taxid, taxon.rank, taxon.name, taxid, 1
                    ),
                )?;
                let entry = unambiguous.entry(taxid).or_insert((
                    taxon.name.clone(),
                    taxon.rank.clone(),
                    0.0,
                    0.0,
                    0.0,
                ));
                entry.2 += 1.0;
                entry.3 += depth;
                entry.4 += maplen;
            } else {
                let key = assignment.taxids.clone();
                ambiguous_reads
                    .entry(key.clone())
                    .or_default()
                    .push(assignment.read_name.clone());
                ambiguous_counts
                    .entry(key)
                    .and_modify(|entry| {
                        entry.1 += 1.0;
                        if entry.2.len() < assignment.depths.len() {
                            entry.2.resize(assignment.depths.len(), 0.0);
                        }
                        if entry.3.len() < assignment.maplens.len() {
                            entry.3.resize(assignment.maplens.len(), 0.0);
                        }
                        for (idx, depth) in assignment.depths.iter().enumerate() {
                            if let Some(val) = entry.2.get_mut(idx) {
                                *val += depth;
                            }
                        }
                        for (idx, maplen) in assignment.maplens.iter().enumerate() {
                            if let Some(val) = entry.3.get_mut(idx) {
                                *val += *maplen;
                            }
                        }
                    })
                    .or_insert_with(|| {
                        (
                            assignment.taxa.clone(),
                            1.0,
                            assignment.depths.clone(),
                            assignment.maplens.clone(),
                        )
                    });
            }
        }

        let mut taxid_count_depth: HashMap<u32, (String, String, f64, f64, f64)> = HashMap::new();
        for (taxid, (name, rank, count, depth, maplen)) in &unambiguous {
            taxid_count_depth.insert(
                *taxid,
                (name.clone(), rank.clone(), *count, *depth, *maplen),
            );
        }
        for taxids in ambiguous_counts.keys() {
            for taxid in taxids {
                taxid_count_depth.entry(*taxid).or_insert((
                    String::new(),
                    String::new(),
                    0.1,
                    0.001,
                    0.001,
                ));
            }
        }

        let mut ambiguous_lca: Vec<Vec<u32>> = Vec::new();
        for (taxids, (taxa, mapcount, depths, maplens)) in &ambiguous_counts {
            if taxids.iter().any(|t| unambiguous.contains_key(t)) {
                continue;
            }
            let taxa_set: HashSet<u32> = taxids.iter().copied().collect();
            if let Some((lca_taxid, lca_name, lca_rank)) =
                self.taxonomy.get_majority_lca(&taxa_set, 0.7)
            {
                ambiguous_lca.push(taxids.clone());
                let entry = unambiguous.entry(lca_taxid).or_insert((
                    lca_name.clone(),
                    lca_rank.clone(),
                    0.0,
                    0.0,
                    0.0,
                ));
                entry.2 += *mapcount;
                entry.3 += depths.iter().copied().sum::<f64>() / depths.len() as f64;
                entry.4 += maplens.iter().copied().sum::<f64>() / maplens.len() as f64;
                if let Some(reads) = ambiguous_reads.get(taxids) {
                    for read in reads {
                        classified_reads += 1;
                        write_optional_line(
                            &mut classification_writer,
                            format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                read,
                                lca_name,
                                lca_taxid,
                                lca_rank,
                                taxa.iter()
                                    .map(|t| t.name.clone())
                                    .collect::<Vec<_>>()
                                    .join(";"),
                                taxids
                                    .iter()
                                    .map(|t| t.to_string())
                                    .collect::<Vec<_>>()
                                    .join(";"),
                                1
                            ),
                        )?;
                    }
                }
            }
        }
        for taxids in ambiguous_lca {
            ambiguous_counts.remove(&taxids);
            ambiguous_reads.remove(&taxids);
        }

        let mut ambiguous_fraction: HashMap<Vec<u32>, IndexMap<u32, (String, String, f64)>> =
            HashMap::new();
        let mut total_taxid_count_depth: HashMap<u32, (String, String, f64, f64, f64)>;
        let mut iterations = 0u32;
        loop {
            total_taxid_count_depth = unambiguous
                .iter()
                .map(|(taxid, (name, rank, count, depth, maplen))| {
                    (
                        *taxid,
                        (name.clone(), rank.clone(), *count, *depth, *maplen),
                    )
                })
                .collect();
            for (taxids, (taxa, mapcount, depths, maplens)) in &ambiguous_counts {
                let sum_depth: f64 = taxids
                    .iter()
                    .map(|taxid| {
                        taxid_count_depth
                            .get(taxid)
                            .map(|(_, _, _, depth, _)| *depth)
                            .unwrap_or(0.0)
                    })
                    .sum::<f64>();
                if sum_depth <= f64::EPSILON {
                    continue;
                }
                let mut paired: Vec<(u32, &TaxonInfo, f64, f64, f64)> = taxids
                    .iter()
                    .enumerate()
                    .map(|(idx, taxid)| {
                        let taxon = &taxa[idx];
                        let depth_val = depths.get(idx).copied().unwrap_or(0.0);
                        let maplen_val = maplens.get(idx).copied().unwrap_or(0.0);
                        let fraction = taxid_count_depth
                            .get(taxid)
                            .map(|(_, _, _, depth, _)| *depth / sum_depth)
                            .unwrap_or(0.0);
                        (*taxid, taxon, depth_val, maplen_val, fraction)
                    })
                    .collect();
                paired.sort_by(|a, b| b.4.partial_cmp(&a.4).unwrap_or(Ordering::Equal));
                let mut fractions: IndexMap<u32, (String, String, f64)> = IndexMap::new();
                let mut consumed = false;
                for (taxid, taxon, depth, maplen, fraction) in paired {
                    if fraction >= 0.99 {
                        let entry = total_taxid_count_depth.entry(taxid).or_insert_with(|| {
                            (taxon.name.clone(), taxon.rank.clone(), 0.1, 0.001, 0.001)
                        });
                        entry.0 = taxon.name.clone();
                        entry.1 = taxon.rank.clone();
                        entry.2 += mapcount;
                        entry.3 += depth;
                        entry.4 += maplen;
                        fractions.insert(taxid, (taxon.name.clone(), taxon.rank.clone(), 1.0));
                        consumed = true;
                        break;
                    }
                    let entry = total_taxid_count_depth.entry(taxid).or_insert_with(|| {
                        (taxon.name.clone(), taxon.rank.clone(), 0.1, 0.001, 0.001)
                    });
                    entry.0 = taxon.name.clone();
                    entry.1 = taxon.rank.clone();
                    entry.2 += mapcount * fraction;
                    entry.3 += depth * fraction;
                    entry.4 += maplen * fraction;
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
            let numerator: f64 = total_taxid_count_depth
                .iter()
                .map(|(taxid, (_, _, _, depth, _))| {
                    let old_depth = taxid_count_depth
                        .get(taxid)
                        .map(|(_, _, _, d, _)| *d)
                        .unwrap_or(0.0);
                    (depth - old_depth).abs()
                })
                .sum::<f64>();
            let denominator: f64 = taxid_count_depth
                .values()
                .map(|(_, _, _, depth, _)| *depth)
                .sum::<f64>()
                .max(1e-12);
            if numerator / denominator <= 1e-6 {
                taxid_count_depth = total_taxid_count_depth.clone();
                break;
            }
            taxid_count_depth = total_taxid_count_depth.clone();
            iterations += 1;
        }
        log::info!(
            target: "PROFILE",
            "Finished strain EM after {} iterations",
            iterations
        );

        for (taxids, reads) in &ambiguous_reads {
            if let Some(fractions) = ambiguous_fraction.get(taxids) {
                let mut sorted: Vec<(&u32, &(String, String, f64))> = fractions.iter().collect();
                sorted.sort_by(|a, b| b.1 .2.partial_cmp(&a.1 .2).unwrap_or(Ordering::Equal));
                let (final_taxid, (final_taxname, final_taxrank, _)) = sorted[0];
                let sorted_taxnames: Vec<String> = sorted
                    .iter()
                    .map(|(_, (name, _, _))| name.clone())
                    .collect();
                let sorted_taxids: Vec<String> =
                    sorted.iter().map(|(taxid, _)| taxid.to_string()).collect();
                let sorted_fractions: Vec<String> = sorted
                    .iter()
                    .map(|(_, (_, _, frac))| format!("{:.5}", frac))
                    .collect();
                for read in reads {
                    classified_reads += 1;
                    write_optional_line(
                        &mut classification_writer,
                        format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            read,
                            final_taxname,
                            final_taxid,
                            final_taxrank,
                            sorted_taxnames.join(";"),
                            sorted_taxids.join(";"),
                            sorted_fractions.join(";"),
                        ),
                    )?;
                }
            }
        }

        let norm_factor: f64 = 0.01
            * taxid_count_depth
                .values()
                .map(|(_, _, _, depth, _)| *depth)
                .sum::<f64>();
        let num_taxa = taxid_count_depth.len();
        let mut out_taxa_list: Vec<Vec<String>> = Vec::new();
        let mut entries: Vec<(u32, (String, String, f64, f64, f64))> = taxid_count_depth
            .iter()
            .map(|(taxid, (name, rank, count, depth, maplen))| {
                (
                    *taxid,
                    (name.clone(), rank.clone(), *count, *depth, *maplen),
                )
            })
            .collect();
        entries.sort_by(|a, b| b.1 .3.partial_cmp(&a.1 .3).unwrap_or(Ordering::Equal));

        for (taxid, (mut name, mut rank, mapcount, depth, maplen)) in entries {
            if mapcount < 0.5 {
                continue;
            }
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
            let scaling = subsample_scaling.get(&taxid).copied().unwrap_or(1.0);
            let scaled_mapcount = mapcount * scaling;
            let coverage = taxid_coverage.get(&taxid).cloned().unwrap_or_default();
            let breadth = coverage.breadth.min(1.0);
            let fixed_chunk = coverage.chunk_breadth.min(1.0);
            let flex_chunk = coverage.flex_chunk_breadth.min(1.0);
            let expected_breadth = if mapcount <= f64::EPSILON {
                0.0
            } else {
                let per_read_cov = (depth / mapcount).clamp(0.0, 1.0);
                1.0 - (1.0 - per_read_cov).powf(mapcount)
            };
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
            let cov_prob = if rank == "species" && expected_breadth > 0.0 {
                let ratio = breadth / expected_breadth;
                if ratio.is_finite() {
                    let read_count = mapcount.round() as usize;
                    if ratio < 0.9 || ratio > 1.2 {
                        let side = if breadth < expected_breadth {
                            Tail::Lower
                        } else if breadth > expected_breadth {
                            Tail::Upper
                        } else {
                            Tail::Upper
                        };
                        let pvalue = if (breadth - expected_breadth).abs() <= f64::EPSILON {
                            Some(1.0)
                        } else {
                            calc_breadth_pvalue(
                                read_count,
                                Some(depth),
                                breadth,
                                Some(expected_breadth),
                                side,
                            )
                        };
                        pvalue.map(|p| (p * num_taxa as f64).min(1.0))
                    } else {
                        Some(1.0)
                    }
                } else {
                    None
                }
            } else {
                None
            };
            let chunkcov_prob = if rank == "species" && expected_flex_chunk > 0.0 {
                let ratio = flex_chunk / expected_flex_chunk;
                if ratio.is_finite() {
                    let read_count = mapcount.round() as usize;
                    if ratio < 0.9 || ratio > 1.2 {
                        let side = if flex_chunk < expected_flex_chunk {
                            Tail::Lower
                        } else if flex_chunk > expected_flex_chunk {
                            Tail::Upper
                        } else {
                            Tail::Upper
                        };
                        let pvalue = if (flex_chunk - expected_flex_chunk).abs() <= f64::EPSILON {
                            Some(1.0)
                        } else {
                            calc_breadth_pvalue(
                                read_count,
                                None,
                                flex_chunk,
                                Some(expected_flex_chunk),
                                side,
                            )
                        };
                        pvalue.map(|p| (p * num_taxa as f64).min(1.0))
                    } else {
                        Some(1.0)
                    }
                } else {
                    None
                }
            } else {
                None
            };
            let cov_prob_str = format_probability_field(cov_prob);
            let chunkcov_prob_str = format_probability_field(chunkcov_prob);

            if let Some(host_taxid) = host_taxid {
                if let Some(map) = pathogen_map.as_ref() {
                    if let Some(entry) = map.get(&taxid) {
                        let host_taxids: Vec<&str> = entry.host_taxids.split(';').collect();
                        let host_taxid_list: Vec<u32> = host_taxids
                            .iter()
                            .filter_map(|t| t.parse::<u32>().ok())
                            .collect();
                        let mut report_it = false;
                        if entry.host_taxids == "unknown" {
                            report_it = true;
                        } else if host_taxid_list.contains(&host_taxid) {
                            report_it = true;
                        } else if lineage_contains(&self.taxonomy, host_taxid, &host_taxid_list) {
                            report_it = true;
                        }
                        if report_it {
                            out_taxa_list.push(vec![
                                name.clone(),
                                taxid.to_string(),
                                rank.clone(),
                                format!("{:.5}", scaled_mapcount),
                                format!("{:.5}", depth),
                                format!("{:.5}", abundance),
                                format!("{:.5}", breadth),
                                format!("{:.5}", expected_breadth),
                                cov_prob_str.clone(),
                                format!("{:.5}", fixed_chunk),
                                format!("{:.5}", flex_chunk),
                                format!("{:.5}", expected_flex_chunk),
                                chunkcov_prob_str.clone(),
                                entry.host_names.clone(),
                                entry.host_taxids.clone(),
                                entry.diseases.clone(),
                            ]);
                            continue;
                        }
                    }
                }
                out_taxa_list.push(vec![
                    name.clone(),
                    taxid.to_string(),
                    rank.clone(),
                    format!("{:.5}", scaled_mapcount),
                    format!("{:.5}", depth),
                    format!("{:.5}", abundance),
                    format!("{:.5}", breadth),
                    format!("{:.5}", expected_breadth),
                    cov_prob_str.clone(),
                    format!("{:.5}", fixed_chunk),
                    format!("{:.5}", flex_chunk),
                    format!("{:.5}", expected_flex_chunk),
                    chunkcov_prob_str.clone(),
                    "unknown".to_string(),
                    "unknown".to_string(),
                    "unknown".to_string(),
                ]);
            } else {
                out_taxa_list.push(vec![
                    name,
                    taxid.to_string(),
                    rank,
                    format!("{:.5}", scaled_mapcount),
                    format!("{:.5}", depth),
                    format!("{:.5}", abundance),
                    format!("{:.5}", breadth),
                    format!("{:.5}", expected_breadth),
                    cov_prob_str,
                    format!("{:.5}", fixed_chunk),
                    format!("{:.5}", flex_chunk),
                    format!("{:.5}", expected_flex_chunk),
                    chunkcov_prob_str,
                ]);
            }
        }

        let raw_taxa_list = out_taxa_list.clone();
        let filtered: Vec<Vec<String>> = out_taxa_list
            .iter()
            .cloned()
            .filter(|entry| {
                if entry.len() < 13 {
                    return false;
                }
                let read_count: f64 = entry[3].parse().unwrap_or(0.0);
                let breadth: f64 = entry[6].parse().unwrap_or(0.0);
                let expected_breadth: f64 = entry[7].parse().unwrap_or(0.0);
                let fixed_chunk: f64 = entry[9].parse().unwrap_or(0.0);
                let flex_chunk: f64 = entry[10].parse().unwrap_or(0.0);
                let expected_flex: f64 = entry[11].parse().unwrap_or(0.0);
                let cov_prob = entry[8].parse::<f64>().ok();
                let chunk_prob = entry[12].parse::<f64>().ok();
                // Probability cutoff for keeping a taxon.
                const PROB_CUTOFF: f64 = 1e-10;
                let pvalue_ok = match (cov_prob, chunk_prob) {
                    (Some(cov), Some(chunk)) => {
                        if self.cfg.strict {
                            cov >= PROB_CUTOFF && chunk >= PROB_CUTOFF
                        } else {
                            !(cov < PROB_CUTOFF && chunk < PROB_CUTOFF)
                        }
                    }
                    _ => !self.cfg.strict,
                };
                let min_reads_ok = match self.cfg.min_reads {
                    Some(min_reads) => read_count.is_finite() && read_count >= min_reads as f64,
                    None => true,
                };
                let fixed_chunk_ok =
                    if self.cfg.min_reads.is_some() && self.cfg.chunk_breadth.is_none() {
                        true
                    } else {
                        fixed_chunk > min_chunk_breadth
                    };
                let oebr = breadth / expected_breadth;
                let coebr = flex_chunk / expected_flex;
                let min_oebr_ok = match self.cfg.min_oebr {
                    Some(min_oebr) => oebr >= min_oebr,
                    None => oebr > 0.75,
                };
                let min_coebr_ok = match self.cfg.min_coebr {
                    Some(min_coebr) => coebr >= min_coebr,
                    None => coebr > 0.75,
                };
                expected_breadth > 0.0
                    && expected_flex > 0.0
                    && oebr.is_finite()
                    && coebr.is_finite()
                    && min_reads_ok
                    && fixed_chunk_ok
                    && min_oebr_ok
                    && oebr <= 1.5
                    && min_coebr_ok
                    && coebr <= 1.5
                    && pvalue_ok
            })
            .collect();

        let total_abundance: f64 = filtered
            .iter()
            .map(|entry| entry[5].parse::<f64>().unwrap_or(0.0))
            .sum();
        let normalized: Vec<Vec<String>> = filtered
            .into_iter()
            .map(|mut entry| {
                if let Ok(value) = entry[5].parse::<f64>() {
                    if total_abundance > 0.0 {
                        entry[5] = format!("{:.5}", value / total_abundance * 100.0);
                    } else {
                        entry[5] = "0.00000".to_string();
                    }
                }
                entry
            })
            .collect();

        for row in &normalized {
            profile_writer.write_record(row)?;
        }
        profile_writer.flush()?;

        if let Some(writer) = raw_writer.as_mut() {
            for row in raw_taxa_list {
                writer.write_record(row)?;
            }
            writer.flush()?;
        }

        if let Some(writer) = classification_writer.as_mut() {
            writer.flush()?;
        }
        Ok(classified_reads)
    }
}

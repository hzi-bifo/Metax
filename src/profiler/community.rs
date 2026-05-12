//! Community Profiler - Species-level taxonomic abundance estimation
//!
//! This module implements the community-level profiling algorithm that groups
//! reads by species taxonomy ID. It is the default profiling mode in Metax.
//!
//! # Algorithm Overview
//!
//! 1. **SAM Parsing**: Load alignments from SAM/BAM files, supporting both
//!    text SAM (optionally gzip/xz compressed) and binary BAM formats.
//!
//! 2. **Reference Metadata**: Extract assembly information from reference
//!    names in the format: `asm|taxid|species_taxid|accession|genome_size`
//!
//! 3. **Coverage Calculation**: Compute coverage metrics:
//!    - Breadth: Fraction of genome covered by at least one read
//!    - Fixed anf flex chunk: Coverage of genome chunks by at least one read
//!
//! 4. **Taxonomy Assignment**: Assign reads to species based on alignment
//!    quality, filtering by identity, mapped length, and fraction thresholds.
//!
//! 5. **Ambiguity Resolution**: Multi-mapping reads are resolved using an
//!    Expectation-Maximization (EM) algorithm that iteratively redistributes
//!    reads based on relative species abundances.
//!
//! 6. **Statistical Filtering**: Final taxa are filtered based on:
//!    - Breadth/expected breadth ratio
//!    - Chunk breadth/expected ratio
//!    - P-values from coverage distribution tests
//!
//! # Key Data Structures
//!
//! - `ReadBundle`: Groups all alignments for a single read
//! - `AlignmentEntry`: Information about a single read-to-reference alignment
//! - `SpeciesCoverage`: Coverage metrics aggregated at species level
//! - `ReadAssignment`: Final taxonomic assignment for a read

use std::cmp::Ordering;
use std::collections::{hash_map::Entry, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::str;
use std::sync::Arc;
use std::thread;

use anyhow::{anyhow, bail, Context, Result};
use crossbeam_channel::{bounded, unbounded};
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
// use once_cell::sync::Lazy;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
// use regex::Regex;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{self, Read as BamRead};
use statrs::function::beta::beta_reg;
use statrs::function::erf::erf;
use statrs::function::gamma::ln_gamma;
use xz2::read::XzDecoder;

use crate::taxonomy::Taxonomy;
use csv::WriterBuilder;

use super::ProfilerConfig;

/// Regex pattern for parsing CIGAR strings.
/// Matches operations like: 100M, 50=, 10X, 5I, 3D, etc.
/// Each match is a number followed by an operation code:
/// - M: alignment match (may include mismatches)
/// - =: sequence match
/// - X: sequence mismatch
/// - I: insertion to reference
/// - D: deletion from reference
/// - N: skipped region
/// - S: soft clipping
/// - H: hard clipping
/// - P: padding
// static CIGAR_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(r"[0-9]+[MIDNSHP=X]").unwrap());

/// Detected alignment file format.
/// Used to select the appropriate parsing strategy.
enum AlignmentFormat {
    /// Text-based SAM format (may be compressed with gzip or xz)
    Sam,
    /// Binary BAM format
    Bam,
}

/// Iterator over lines in a text-based SAM file.
///
/// Wraps a buffered reader and yields non-empty, trimmed lines.
/// Skips empty lines and properly handles line endings.
struct SamLineIter<R: BufRead> {
    reader: R,
    /// Reusable line buffer. Kept on the struct so the underlying
    /// allocation is amortised across all lines (previously this
    /// iterator allocated a fresh `String` on every `next()` call,
    /// which on ~20 M alignments dominated the main-thread cost).
    buffer: String,
}

impl<R: BufRead> SamLineIter<R> {
    /// Create a new SAM line iterator from a buffered reader.
    fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: String::with_capacity(4096),
        }
    }
}

impl<R: BufRead> Iterator for SamLineIter<R> {
    type Item = Result<String>;

    /// Read the next non-empty line from the SAM file.
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.buffer.clear();
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => return None, // EOF reached
                Ok(_) => {
                    let trimmed = self.buffer.trim_end();
                    if trimmed.is_empty() {
                        continue; // Skip empty lines
                    }
                    return Some(Ok(trimmed.to_string()));
                }
                Err(err) => return Some(Err(anyhow!(err))),
            }
        }
    }
}

/// Iterator over records in a BAM file, converting to SAM format.
///
/// Uses rust-htslib to read BAM records and converts them to SAM-format
/// strings for uniform processing with the text SAM parser.
struct BamRecordIter {
    /// BAM file reader
    reader: bam::Reader,
    /// Reusable record buffer to avoid allocations
    record: Record,
    /// BAM header containing reference sequence information
    header: Arc<bam::HeaderView>,
}

impl BamRecordIter {
    /// Create a new BAM record iterator.
    ///
    /// # Arguments
    /// * `reader` - Opened BAM file reader
    /// * `header` - BAM header view (shared for memory efficiency)
    fn new(reader: bam::Reader, header: Arc<bam::HeaderView>) -> Self {
        Self {
            reader,
            record: Record::new(),
            header,
        }
    }
}

impl Iterator for BamRecordIter {
    type Item = Result<String>;

    /// Read the next BAM record and convert to SAM format string.
    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read(&mut self.record) {
            Some(Ok(())) => Some(bam_record_to_sam_line(&self.record, &self.header)),
            Some(Err(err)) => Some(Err(anyhow!(err))),
            None => None,
        }
    }
}

/// Information about a single alignment of a read to a reference.
///
/// All "derived" fields that used to be recomputed for every alignment
/// inside `process_batch` are now parsed / computed ONCE by the SAM/BAM
/// worker (in [`parse_alignment_line`]) and stored here as ready-to-use
/// values. The downstream taxonomy-assignment loop is therefore purely
/// arithmetic and does not touch any `HashMap`, `String::split`, or
/// CIGAR scanner per alignment.
///
/// Fields marked *derived* are produced from the original SAM record
/// and are invariant for a given alignment, so pre-computing them does
/// not change any output.
#[derive(Clone, Debug)]
struct AlignmentEntry {
    /// Assembly identifier (first `|`-separated field of the reference
    /// name). *Derived.* Used as part of [`AccKey`] and for the
    /// assembly-length lookup.
    asm: String,

    /// Accession identifier (fourth `|`-separated field of the
    /// reference name, or empty if absent). *Derived.* Part of
    /// [`AccKey`].
    acc: String,

    /// Species taxid after [`build_taxonomy_lookups`] normalization
    /// (sub-species / strain-level taxids are mapped up to their
    /// species ancestor; taxids already at or above species rank pass
    /// through unchanged). *Derived.*
    species_taxid: u32,

    /// Total assembly length (sum of `LN:` values for every @SQ whose
    /// SN shares this `asm`). *Derived, always positive.*
    asm_len: f64,

    /// Gap-compressed identity from the CIGAR (matches and mismatches;
    /// see [`cigar_to_identity_with_opts`]). *Derived.*
    identity: f64,

    /// Aligned fraction from the CIGAR (fraction of the read that was
    /// actually aligned). *Derived.*
    fraction: f64,

    /// Number of read bases that aligned (excluding hard clips).
    mapped_len: usize,

    /// Genomic span of the alignment (start, end) in 1-based coordinates.
    span: (usize, usize),

    /// Fixed-size chunk around alignment for coverage QC (start, end).
    /// Extends alignment by `fixed_flank_len` (1kb for paired, 2kb for
    /// single).
    fixed_chunk: (usize, usize),

    /// Flexible-size chunk for coverage QC (start, end).
    /// Size is proportional to √(genome_size).
    flex_chunk: (usize, usize),
}

/// Bundle of all alignments for a single read.
///
/// A read may align to multiple reference sequences (multi-mapping).
/// This structure groups all alignments for downstream processing.
#[derive(Clone, Debug)]
struct ReadBundle {
    /// Read identifier (pair suffix stripped for paired-end)
    read_name: String,

    /// All alignments for this read, ordered by quality
    alignments: Vec<AlignmentEntry>,
}

/// Result from a worker thread processing a batch of alignments.
///
/// Workers parse alignment lines in parallel and return structured data
/// for the main thread to aggregate.
struct WorkerResult {
    /// Parsed alignments as (read_name, entry) pairs
    alignments: Vec<(String, AlignmentEntry)>,

    /// All read names seen (including unmapped)
    reads_seen: Vec<String>,

    /// Count of valid alignments processed (for progress tracking)
    alignments_seen: usize,

    /// Alignments rejected in the worker because the reference name
    /// didn't carry a parsable species taxid, or because the assembly
    /// total length was unknown / zero. Historically counted inside
    /// `process_batch`; we now detect these in the worker to avoid
    /// re-parsing, and propagate the count so the aggregated
    /// [`FilterStats`] keeps the same meaning.
    skipped_reference: usize,
}

/// Shared context for worker threads during alignment parsing.
///
/// Contains reference metadata needed to compute coverage coordinates.
/// Wrapped in Arc for thread-safe sharing.
struct WorkerContext {
    /// Reference sequence lengths: ref_name -> length
    ref_lengths: Arc<HashMap<String, usize>>,

    /// Per-assembly total length: asm_id -> sum(LN) over all @SQ
    /// sharing this `asm`. Used to pre-compute the per-alignment
    /// assembly length so the hot loop never has to re-lookup it.
    asm_lengths: Arc<HashMap<String, f64>>,

    /// Sparse sub-species → species normalization map. Raw taxids that
    /// are below species rank map to their species ancestor; raw
    /// taxids already at (or above) species rank are absent and the
    /// worker falls back to the raw value.
    species_normalize: Arc<HashMap<u32, u32>>,

    /// Flexible chunk flank lengths per assembly: asm_id -> √(asm_length)
    asm_flex_lengths: Arc<HashMap<String, usize>>,

    /// Fixed flank length for chunk coverage (1000 for paired, 2000 for single)
    fixed_flank_len: usize,

    /// Whether input is paired-end (affects read name parsing)
    is_paired: bool,
}

/// Intermediate result from parsing a single alignment line.
///
/// Separates the alignment entry (if valid) from metadata about whether
/// to count this alignment in statistics.
struct ParsedAlignment {
    /// Read name (pair suffix stripped if paired-end)
    read_name: String,

    /// Parsed alignment entry, or None if unmapped/filtered
    entry: Option<AlignmentEntry>,

    /// Whether to count this toward alignment statistics
    /// False for unmapped, secondary alignments, or mate-2 of pairs
    count_alignment: bool,

    /// `true` when the alignment was discarded because the reference
    /// name did not parse into the expected `asm|_|species|acc|...`
    /// layout, or because the assembly total length was unknown / zero
    /// (i.e. the failure mode that used to be counted as
    /// `FilterStats::skipped_reference` inside `process_batch`).
    skipped_reference: bool,
}

/// Coverage metrics aggregated at the species level.
#[derive(Clone, Debug, Default)]
struct SpeciesCoverage {
    /// Breadth of coverage: fraction of genome with ≥1 read
    breadth: f64,

    /// Fixed chunk coverage: uses fixed-size flanking regions
    fixed_chunk: f64,

    /// Flexible chunk coverage: uses √(genome_size)-sized regions
    flex_chunk: f64,
}

/// Unique key identifying a specific accession within an assembly.
///
/// Used as a HashMap key for aggregating coverage intervals.
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct AccKey {
    /// Species taxonomy ID
    taxid: u32,

    /// Assembly identifier (first field of reference name)
    asm: String,

    /// Accession identifier (fourth field of reference name)
    acc: String,
}

/// Statistics about alignment filtering.
///
/// Tracks how many alignments were evaluated, passed filters,
/// and why alignments were rejected.
#[derive(Clone, Copy, Debug, Default)]
struct FilterStats {
    /// Total alignments evaluated
    evaluated: usize,

    /// Alignments that passed all filters
    passed: usize,

    /// Filtered due to insufficient mapped length
    filtered_mapped_len: usize,

    /// Filtered due to low sequence identity
    filtered_identity: usize,

    /// Filtered due to low aligned fraction
    filtered_fraction: usize,

    /// Skipped due to missing reference metadata
    skipped_reference: usize,
}

impl FilterStats {
    /// Merge statistics from another FilterStats instance.
    /// Used to aggregate results from parallel workers.
    fn merge(&mut self, other: FilterStats) {
        self.evaluated += other.evaluated;
        self.passed += other.passed;
        self.filtered_mapped_len += other.filtered_mapped_len;
        self.filtered_identity += other.filtered_identity;
        self.filtered_fraction += other.filtered_fraction;
        self.skipped_reference += other.skipped_reference;
    }
}

/// Results from the taxonomy assignment phase.
///
/// Contains all data needed for abundance estimation and output generation.
#[derive(Default)]
struct AssignmentResult {
    /// Per-read taxonomic assignments
    read_assignments: Vec<ReadAssignment>,

    /// Alignment intervals per accession (for breadth calculation)
    acc_intervals: HashMap<AccKey, Vec<(usize, usize)>>,

    /// Fixed chunk intervals per accession
    acc_fixed_chunk_intervals: HashMap<AccKey, Vec<(usize, usize)>>,

    /// Flexible chunk intervals per accession
    acc_flex_chunk_intervals: HashMap<AccKey, Vec<(usize, usize)>>,

    /// Filtering statistics
    filter_stats: FilterStats,
}

/// Rank at which a read is taxonomically assigned.
///
/// Reads are normally assigned at the `Species` rank. When the genus
/// fallback mode activates (low-mapping Illumina, pathogen mode off),
/// reads that fail species thresholds but meet the fallback criteria
/// are assigned at the `Genus` rank instead.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum AssignRank {
    Species,
    Genus,
}

impl AssignRank {
    /// Return the canonical lowercase rank label used in output files.
    fn as_str(&self) -> &'static str {
        match self {
            AssignRank::Species => "species",
            AssignRank::Genus => "genus",
        }
    }
}

/// Taxonomic assignment for a single read.
///
/// A read may be assigned to one taxon (unambiguous) or multiple taxa
/// (ambiguous multi-mapping). Ambiguous assignments are resolved in
/// the EM phase.
///
/// Each `ReadAssignment` is homogeneous in rank: all candidate taxids
/// are either species (normal path) or genera (fallback path).
#[derive(Clone, Debug)]
struct ReadAssignment {
    /// Read identifier
    read: String,

    /// Assignment rank (species or genus fallback).
    rank: AssignRank,

    /// Candidate taxonomy IDs at `rank` (may be >1 for multi-mappers)
    taxids: Vec<u32>,

    /// Depth contribution for each candidate (mapped_len / genome_size)
    depths: Vec<f64>,

    /// Mapped lengths for each candidate
    mapped_lens: Vec<f64>,
}

/// Pathogen metadata from the pathogen-host mapping file.
///
/// Used to annotate results with host and disease information
/// when pathogen filtering is enabled.
#[derive(Clone, Debug)]
struct PathogenEntry {
    /// Semicolon-separated host taxonomy IDs
    host_taxids: String,

    /// Semicolon-separated host names
    host_names: String,

    /// Semicolon-separated associated diseases
    diseases: String,
}

/// Community-level taxonomy profiler.
///
/// Implements species-level profiling by grouping alignments by species
/// taxonomy ID and resolving ambiguous assignments using EM.
///
/// # Pipeline
///
/// 1. Load SAM/BAM file and extract reference metadata
/// 2. Compute coverage intervals for each reference sequence
/// 3. Assign reads to taxa based on alignment quality
/// 4. Resolve multi-mapping reads with EM algorithm
/// 5. Filter taxa based on coverage statistics
/// 6. Output abundance profile and classification files
pub struct CommunityProfiler {
    /// Profiler configuration
    cfg: ProfilerConfig,

    /// Taxonomy database (shared across threads)
    taxonomy: Arc<Taxonomy>,
}

impl CommunityProfiler {
    /// Create a new community profiler with the given configuration.
    ///
    /// # Arguments
    /// * `cfg` - Profiler configuration
    /// * `taxonomy` - NCBI taxonomy database (wrapped in Arc)
    pub fn new(cfg: ProfilerConfig, taxonomy: Arc<Taxonomy>) -> Result<Self> {
        Ok(Self { cfg, taxonomy })
    }

    /// Run the community profiling pipeline.
    ///
    /// This is the main entry point that orchestrates the full pipeline:
    /// 1. Load and parse alignments
    /// 2. Compute coverage metrics
    /// 3. Assign taxonomy
    /// 4. Write output files
    pub fn run(&self) -> Result<()> {
        let sam_path = Path::new(&self.cfg.sam);
        if !sam_path.exists() {
            bail!("SAM file {} does not exist", sam_path.display());
        }
        log::info!(target: "PROFILE", "Initializing the profiler!");

        // Load and parse alignment data
        let sam_data = self.load_sam(sam_path)?;

        // Run profiling pipeline
        self.profile(sam_data)
    }

    /// Load alignments from a SAM or BAM file.
    ///
    /// Automatically detects the file format (SAM/BAM) and handles
    /// compression (gzip/xz for SAM files).
    fn load_sam(&self, sam_path: &Path) -> Result<SamData> {
        log::info!(target: "PROFILE", "Loading alignments ...");
        match detect_alignment_format(sam_path)? {
            AlignmentFormat::Sam => self.load_text_alignments(sam_path),
            AlignmentFormat::Bam => self.load_bam_alignments(sam_path),
        }
    }

    /// Load alignments from a text-based SAM file.
    ///
    /// Handles plain text SAM and compressed formats (gzip, xz).
    /// Parses headers to extract reference sequence metadata.
    ///
    /// # SAM Header Processing
    ///
    /// @SQ lines are parsed to extract:
    /// - SN: Reference sequence name (contains taxonomy info)
    /// - LN: Reference sequence length
    ///
    /// Reference names are expected in format:
    /// `asm_id|taxid|species_taxid|accession|genome_size[|sampled_size]`
    fn load_text_alignments(&self, sam_path: &Path) -> Result<SamData> {
        let reader = open_text_alignment_reader(sam_path)?;
        let mut reader = BufReader::new(reader);

        // Maps to store reference metadata from SAM header.
        let mut ref_lengths: HashMap<String, usize> = HashMap::new();
        let mut asm_lengths: HashMap<String, f64> = HashMap::new();
        let mut subsample_scaling: HashMap<u32, f64> = HashMap::new();
        let mut species_taxids: HashSet<u32> = HashSet::new();
        // Whether the first @SQ line's SN looks like a subsampled-DB entry
        // (6 `|`-separated fields). Used later to suppress automatic
        // low-map fallback and to choose the chunk-breadth basis for
        // explicit `--genus-fallback`.
        let mut db_is_subsampled: Option<bool> = None;

        // Fixed flank length: shorter for paired-end (reads cover more area together)
        let fixed_flank_len = if self.cfg.is_paired { 1000 } else { 2000 };
        let mut first_alignment: Option<String> = None;
        let mut line = String::new();

        // Parse SAM header lines (start with @)
        loop {
            line.clear();
            let bytes = reader.read_line(&mut line)?;
            if bytes == 0 {
                break;
            }
            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('@') {
                // Parse @SQ (sequence dictionary) lines for reference info
                if trimmed.starts_with("@SQ") {
                    let mut sn: Option<String> = None;
                    let mut ln: Option<usize> = None;

                    // Extract SN (sequence name) and LN (length) tags
                    for part in trimmed.split('\t').skip(1) {
                        if part.starts_with("SN:") {
                            sn = Some(part[3..].to_string());
                        } else if part.starts_with("LN:") {
                            ln = Some(part[3..].parse()?);
                        }
                    }

                    if let (Some(sn), Some(len)) = (sn, ln) {
                        // Inspect ONLY the first @SQ line to decide whether
                        // the run is against a subsampled reference DB.
                        if db_is_subsampled.is_none() {
                            db_is_subsampled = Some(ref_name_is_subsampled(&sn));
                        }
                        ref_lengths.insert(sn.clone(), len);

                        // Parse reference name to extract assembly ID
                        let parts: Vec<&str> = sn.split('|').collect();
                        let asm = parts.first().copied().unwrap_or(&sn).to_string();

                        // Accumulate total assembly length
                        *asm_lengths.entry(asm).or_insert(0.0) += len as f64;

                        // Collect species taxid (third field) for genus lookup.
                        if let Some(species_taxid) =
                            parts.get(2).and_then(|s| s.parse::<u32>().ok())
                        {
                            species_taxids.insert(species_taxid);
                        }

                        // Extract subsample scaling factor if present
                        update_subsample_scaling(&mut subsample_scaling, &parts);
                    }
                }
                continue;
            } else {
                // First alignment line found - save and exit header parsing
                first_alignment = Some(trimmed.to_string());
                break;
            }
        }

        // Calculate flexible chunk lengths as √(assembly_length)
        // This scales coverage checking with genome size
        let asm_flex_lengths: HashMap<String, usize> = asm_lengths
            .iter()
            .map(|(asm, len)| {
                let flex = len.sqrt().round() as usize;
                (asm.clone(), flex.max(1))
            })
            .collect();

        // Precompute taxonomy lookups from the lineage. The species
        // normalization map is ALWAYS built (it applies regardless of
        // fallback mode -- sub-species still need to be merged into
        // species). The species->genus map is needed for normal databases
        // and for explicit `--genus-fallback` on subsampled databases.
        // Built BEFORE `WorkerContext` so the worker can pre-normalize
        // species taxids during parse and skip an extra HashMap lookup
        // per alignment in the hot loop.
        let db_is_subsampled = db_is_subsampled.unwrap_or(false);
        let (species_normalize, species_to_genus) = self.build_taxonomy_lookups(&species_taxids);
        let species_normalize = Arc::new(species_normalize);
        let species_to_genus = Arc::new(if db_is_subsampled && !self.cfg.genus_fallback {
            HashMap::new()
        } else {
            species_to_genus
        });
        let asm_lengths = Arc::new(asm_lengths);

        // Create shared context for worker threads (now carries the
        // pre-built taxonomy and assembly-length lookups).
        let worker_context = Arc::new(WorkerContext {
            ref_lengths: Arc::new(ref_lengths),
            asm_lengths: Arc::clone(&asm_lengths),
            species_normalize: Arc::clone(&species_normalize),
            asm_flex_lengths: Arc::new(asm_flex_lengths),
            fixed_flank_len,
            is_paired: self.cfg.is_paired,
        });

        // Process alignment records through the parallel pipeline
        let iter = SamLineIter::new(reader);
        self.consume_alignment_stream(
            first_alignment,
            iter,
            worker_context,
            asm_lengths,
            subsample_scaling,
            species_normalize,
            species_to_genus,
            db_is_subsampled,
        )
    }

    /// Load alignments from a binary BAM file.
    ///
    /// Uses rust-htslib for BAM parsing with optional multi-threaded
    /// decompression. Converts BAM records to SAM format strings for
    /// uniform processing.
    fn load_bam_alignments(&self, sam_path: &Path) -> Result<SamData> {
        let mut reader = bam::Reader::from_path(sam_path)
            .with_context(|| format!("failed to open {}", sam_path.display()))?;
        if self.cfg.threads > 1 {
            if let Err(err) = reader.set_threads(self.cfg.threads) {
                log::debug!(
                    target: "PROFILE",
                    "Unable to enable multi-threaded decoding: {}",
                    err
                );
            }
        }
        let header_view = Arc::new(reader.header().clone());

        let mut ref_lengths: HashMap<String, usize> = HashMap::new();
        let mut asm_lengths: HashMap<String, f64> = HashMap::new();
        let mut subsample_scaling: HashMap<u32, f64> = HashMap::new();
        let mut species_taxids: HashSet<u32> = HashSet::new();
        // Decide whether the DB is subsampled from the very first target
        // name (parallels the @SQ-first-line check used for text SAM).
        let mut db_is_subsampled: Option<bool> = None;
        for (idx, name_bytes) in header_view.target_names().iter().enumerate() {
            let name = str::from_utf8(name_bytes)
                .with_context(|| format!("invalid reference name at index {}", idx))?
                .to_string();
            if db_is_subsampled.is_none() {
                db_is_subsampled = Some(ref_name_is_subsampled(&name));
            }
            let len = header_view.target_len(idx as u32).unwrap_or(0) as usize;
            ref_lengths.insert(name.clone(), len);
            let parts: Vec<&str> = name.split('|').collect();
            let asm = parts.first().copied().unwrap_or(&name).to_string();
            *asm_lengths.entry(asm).or_insert(0.0) += len as f64;
            if let Some(species_taxid) = parts.get(2).and_then(|s| s.parse::<u32>().ok()) {
                species_taxids.insert(species_taxid);
            }
            update_subsample_scaling(&mut subsample_scaling, &parts);
        }

        let asm_flex_lengths: HashMap<String, usize> = asm_lengths
            .iter()
            .map(|(asm, len)| {
                let flex = len.sqrt().round() as usize;
                (asm.clone(), flex.max(1))
            })
            .collect();

        // Precompute both the species-normalization map (always needed)
        // and the species->genus lookup. For subsampled DBs the automatic
        // low-map fallback remains disabled, but explicit `--genus-fallback`
        // is allowed and therefore still needs the genus map. Built BEFORE
        // `WorkerContext` so the worker can pre-normalize.
        let db_is_subsampled = db_is_subsampled.unwrap_or(false);
        let (species_normalize, species_to_genus) = self.build_taxonomy_lookups(&species_taxids);
        let species_normalize = Arc::new(species_normalize);
        let species_to_genus = Arc::new(if db_is_subsampled && !self.cfg.genus_fallback {
            HashMap::new()
        } else {
            species_to_genus
        });
        let asm_lengths = Arc::new(asm_lengths);

        let worker_context = Arc::new(WorkerContext {
            ref_lengths: Arc::new(ref_lengths),
            asm_lengths: Arc::clone(&asm_lengths),
            species_normalize: Arc::clone(&species_normalize),
            asm_flex_lengths: Arc::new(asm_flex_lengths),
            fixed_flank_len: if self.cfg.is_paired { 1000 } else { 2000 },
            is_paired: self.cfg.is_paired,
        });

        let iter = BamRecordIter::new(reader, header_view);
        self.consume_alignment_stream(
            None,
            iter,
            worker_context,
            asm_lengths,
            subsample_scaling,
            species_normalize,
            species_to_genus,
            db_is_subsampled,
        )
    }

    /// Precompute taxonomy lookups used during per-read assignment.
    ///
    /// Given the set of raw taxids observed in the SAM/BAM @SQ headers
    /// (third `|`-separated field), this produces:
    ///
    /// - `species_normalize`: sparse `raw -> species` map. Only contains
    ///   entries where the raw taxid is **below** species rank and a
    ///   species ancestor exists. Entries where raw is already species
    ///   (or has no species ancestor) are omitted; callers treat a
    ///   missing entry as "use the raw taxid unchanged". This is the
    ///   mechanism that merges multiple sub-species / strains belonging
    ///   to the same species into a single species-level assignment.
    ///
    /// - `species_to_genus`: `species -> genus` map, keyed by the
    ///   **post-normalization** species taxid so that the genus
    ///   fallback path looks up the genus once per species regardless
    ///   of how many sub-taxon flavours share that species.
    ///
    /// Taxids without a resolvable genus are simply absent from
    /// `species_to_genus`; the genus fallback silently skips them.
    fn build_taxonomy_lookups(
        &self,
        raw_taxids: &HashSet<u32>,
    ) -> (HashMap<u32, u32>, HashMap<u32, u32>) {
        let tax = self.taxonomy.as_ref();
        let mut species_normalize: HashMap<u32, u32> = HashMap::new();
        let mut species_to_genus: HashMap<u32, u32> = HashMap::with_capacity(raw_taxids.len());
        for &raw in raw_taxids {
            // Determine the species anchor (or fall back to the raw id
            // when no species ancestor exists, to preserve behaviour for
            // references that already point above species rank).
            let species = match resolve_species_taxid(tax, raw) {
                Some(s) => {
                    if s != raw {
                        species_normalize.insert(raw, s);
                    }
                    s
                }
                None => raw,
            };
            if let Entry::Vacant(slot) = species_to_genus.entry(species) {
                if let Some(g) = resolve_genus_taxid(tax, species) {
                    slot.insert(g);
                }
            }
        }
        (species_normalize, species_to_genus)
    }

    /// Process alignment records through a parallel worker pipeline.
    ///
    /// This is the main alignment processing function that:
    /// 1. Spawns worker threads to parse alignment lines in parallel
    /// 2. Dispatches batches of lines to workers via channels
    /// 3. Collects results and aggregates read alignments
    /// 4. Computes coverage statistics and filtering thresholds
    ///
    /// # Architecture
    ///
    /// ```text
    /// Main Thread              Worker Threads
    ///     |                         |
    ///     |--[batch]--------------->| parse_alignment_line()
    ///     |<--[WorkerResult]--------|
    ///     |                         |
    ///     v                         v
    /// Aggregation              Next batch
    /// ```
    ///
    /// # Arguments
    /// * `first_record` - First alignment line (already read during header parsing)
    /// * `iter` - Iterator over remaining alignment lines
    /// * `worker_context` - Shared reference metadata for workers
    /// * `asm_lengths` - Total assembly lengths (for coverage normalization)
    /// * `subsample_scaling` - Scaling factors for subsampled genomes
    fn consume_alignment_stream<I>(
        &self,
        first_record: Option<String>,
        iter: I,
        worker_context: Arc<WorkerContext>,
        asm_lengths: Arc<HashMap<String, f64>>,
        subsample_scaling: HashMap<u32, f64>,
        species_normalize: Arc<HashMap<u32, u32>>,
        species_to_genus: Arc<HashMap<u32, u32>>,
        db_is_subsampled: bool,
    ) -> Result<SamData>
    where
        I: Iterator<Item = Result<String>>,
    {
        let worker_count = self.cfg.threads.max(1);
        let batch_size = self.cfg.batch_size.max(1);

        // Create channels for work distribution and result collection
        // bounded channel with 2x workers prevents memory bloat
        let (line_tx, line_rx) = bounded::<Vec<String>>(worker_count * 2);
        let (result_tx, result_rx) = unbounded::<WorkerResult>();

        // Spawn worker threads
        let mut handles = Vec::new();
        for _ in 0..worker_count {
            let rx = line_rx.clone();
            let tx = result_tx.clone();
            let ctx = Arc::clone(&worker_context);

            handles.push(thread::spawn(move || {
                // Process batches until channel closes
                for chunk in rx.iter() {
                    let chunk_len = chunk.len();
                    let mut reads_seen: Vec<String> = Vec::with_capacity(chunk_len);
                    let mut alignments: Vec<(String, AlignmentEntry)> =
                        Vec::with_capacity(chunk_len);
                    let mut alignments_seen = 0usize;
                    let mut skipped_reference = 0usize;

                    // Parse each alignment line in the batch
                    for line in chunk {
                        if let Some(parsed) = parse_alignment_line(&line, &ctx) {
                            let ParsedAlignment {
                                read_name,
                                entry,
                                count_alignment,
                                skipped_reference: was_skipped_ref,
                            } = parsed;

                            if count_alignment {
                                alignments_seen += 1;
                            }
                            if was_skipped_ref {
                                skipped_reference += 1;
                            }
                            // For rows that produced a kept alignment,
                            // the `read_name` is already carried as the
                            // key of the alignment tuple — it will
                            // therefore surface in `read_map.keys()` on
                            // the main thread, so we skip pushing it
                            // into `reads_seen` to avoid one `String`
                            // clone per kept alignment. Unmapped /
                            // discarded rows still push into
                            // `reads_seen` (consuming the `String`, no
                            // clone) so the main thread can count them
                            // toward the total-read denominator.
                            match entry {
                                Some(aln) => {
                                    alignments.push((read_name, aln));
                                }
                                None => {
                                    reads_seen.push(read_name);
                                }
                            }
                        }
                    }

                    // Send results back to main thread
                    if !reads_seen.is_empty() || alignments_seen > 0 || skipped_reference > 0 {
                        let _ = tx.send(WorkerResult {
                            alignments,
                            reads_seen,
                            alignments_seen,
                            skipped_reference,
                        });
                    }
                }
            }));
        }
        // Drop our copy of result_tx so channel closes when workers finish
        drop(result_tx);

        // Main-thread aggregation: one `.entry()` per alignment on
        // `read_map`, one `.insert()` per unmapped/discarded row on
        // `all_reads` (aligned rows carry `read_name` only via
        // `read_map`; see end-of-stream read-count reconciliation).
        let mut read_map: HashMap<String, Vec<AlignmentEntry>> = HashMap::new();
        let mut all_reads: HashSet<String> = HashSet::new();
        let mut num_alignments: usize = 0;
        let mut num_skipped_reference: usize = 0;
        let mut next_progress: usize = 1_000_000;

        let mut chunk: Vec<String> = Vec::with_capacity(batch_size);
        if let Some(first) = first_record {
            if !first.is_empty() {
                chunk.push(first);
            }
        }

        for line in iter {
            let record = line?;
            if record.is_empty() {
                continue;
            }
            chunk.push(record);
            if chunk.len() >= batch_size {
                line_tx
                    .send(chunk)
                    .map_err(|e| anyhow!("failed to dispatch SAM chunk: {}", e))?;
                chunk = Vec::with_capacity(batch_size);
            }
            while let Ok(result) = result_rx.try_recv() {
                absorb_worker_result(
                    result,
                    &mut read_map,
                    &mut all_reads,
                    &mut num_alignments,
                    &mut num_skipped_reference,
                    &mut next_progress,
                );
            }
        }

        if !chunk.is_empty() {
            line_tx
                .send(chunk)
                .map_err(|e| anyhow!("failed to dispatch final SAM chunk: {}", e))?;
        }
        drop(line_tx);

        for result in result_rx.iter() {
            absorb_worker_result(
                result,
                &mut read_map,
                &mut all_reads,
                &mut num_alignments,
                &mut num_skipped_reference,
                &mut next_progress,
            );
        }

        for handle in handles {
            let _ = handle.join();
        }

        // Total distinct reads = (unmapped/discarded rows tracked in
        // `all_reads`) ∪ (rows with at least one kept alignment, which
        // live in `read_map` as the sole copy of their `read_name`).
        // Workers skip pushing aligned names into `reads_seen` to save
        // one `String` clone per alignment, so we count the aligned
        // reads that are **not already** in `all_reads` here. A read
        // with both mapped and unmapped rows is correctly counted
        // once because it's present in `all_reads`.
        let aligned_unique_not_in_all_reads = read_map
            .keys()
            .filter(|k| !all_reads.contains(k.as_str()))
            .count();
        let num_reads = all_reads.len() + aligned_unique_not_in_all_reads;
        let aligned_reads = read_map.len();
        let read_alignments = read_map
            .into_iter()
            .map(|(read_name, alignments)| ReadBundle {
                read_name,
                alignments,
            })
            .collect::<Vec<_>>();
        // Mapping rate must be computed before the chunk-breadth basis can
        // be auto-selected, since the auto-rule consults it.
        let mapping_rate = if num_reads == 0 {
            0.0
        } else {
            aligned_reads as f64 / num_reads as f64
        };
        let gate_active = fallback_gate_active(&self.cfg, mapping_rate);
        let fallback_mode = genus_fallback_mode_active(&self.cfg, gate_active, db_is_subsampled);
        let chunk_plan = plan_chunk_breadth(
            &self.cfg,
            mapping_rate,
            num_reads,
            aligned_reads,
            db_is_subsampled,
            gate_active,
        );
        let min_breadth = self.cfg.breadth.unwrap_or(0.0);
        let mut min_chunk_breadth = self.cfg.chunk_breadth;
        let skip_auto_chunk_filter = self.cfg.min_reads.is_some() && min_chunk_breadth.is_none();
        if min_chunk_breadth.is_none() && !skip_auto_chunk_filter {
            min_chunk_breadth = Some(self.estimate_chunk_breadth(chunk_plan.estimate_reads));
        }

        if self.cfg.verbose || self.cfg.very_verbose {
            const INDENT: &str = "                                                ";
            log::info!(
                target: "PROFILE",
                "Profiler parameters:\n{}min_breadth        {}\n{}min_chunk_breadth  {}\n{}min_reads          {}\n{}min_oebr           {}\n{}min_coebr          {}\n{}identity           {}\n{}mapped_len   {}\n{}fraction           {}\n{}batch_size         {}\n",
                INDENT,
                min_breadth,
                INDENT,
                min_chunk_breadth.unwrap_or(0.0),
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

            log::info!(target: "PROFILE", "Total number of reads: {}", num_reads);
            log::info!(target: "PROFILE", "Number of alignments: {}", num_alignments);
            log::info!(
                target: "PROFILE",
                "Number of aligned reads: {}",
                aligned_reads
            );
            if self.cfg.chunk_breadth.is_some() {
                log::info!(
                    target: "PROFILE",
                    "Minimum chunk breadth supplied by user; read-count basis not used for estimation."
                );
            } else if skip_auto_chunk_filter {
                log::info!(
                    target: "PROFILE",
                    "Automatic minimum chunk breadth filter disabled because --min-reads is set."
                );
            } else {
                // Log the unscaled basis selected for `estimate_chunk_breadth`.
                // Any optional 1.5x boost is logged separately below.
                if !chunk_plan.use_aligned_basis {
                    log::info!(
                        target: "PROFILE",
                        "Minimum chunk breadth estimated from total reads ({}).",
                        num_reads
                    );
                    if chunk_plan.force_total_for_subsampled_fallback {
                        log::info!(
                            target: "PROFILE",
                            "--genus-fallback with a subsampled reference DB: chunk-breadth estimate uses total read count."
                        );
                    }
                } else if chunk_plan.by_aligned_auto {
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
                if chunk_plan.scaled {
                    log::info!(
                        target: "PROFILE",
                        "Low-map genus fallback active and --min-cbreadth not set: aligned-read chunk basis scaled by {:.2}x ({} -> {}) before applying the chunk-breadth equation.",
                        FALLBACK_CHUNK_BASIS_SCALE,
                        chunk_plan.basis_reads,
                        chunk_plan.estimate_reads,
                    );
                }
            }
        }

        if self.cfg.verbose || self.cfg.very_verbose {
            log::info!(
                target: "PROFILE",
                "Mapping rate: {:.5} (aligned {} / total {})",
                mapping_rate,
                aligned_reads,
                num_reads
            );
            if fallback_mode {
                if gate_active && !db_is_subsampled {
                    log::info!(
                        target: "PROFILE",
                        "Genus-level fallback mode ENABLED (low-map Illumina, host=None, mapping_rate {:.5} < {:.5}; genus_identity={:.3}).",
                        mapping_rate,
                        self.cfg.low_map_rate_threshold,
                        self.cfg.genus_identity,
                    );
                } else if self.cfg.genus_fallback && db_is_subsampled {
                    log::info!(
                        target: "PROFILE",
                        "Genus-level fallback mode ENABLED (--genus-fallback on subsampled reference DB; genus_identity={:.3}).",
                        self.cfg.genus_identity,
                    );
                } else {
                    log::info!(
                        target: "PROFILE",
                        "Genus-level fallback mode ENABLED (--genus-fallback; genus_identity={:.3}).",
                        self.cfg.genus_identity,
                    );
                }
            } else if gate_active && db_is_subsampled {
                // The automatic low-map assignment path would have been on,
                // but a subsampled reference database suppresses the
                // automatic fallback. An explicit --genus-fallback still
                // enables the assignment path.
                log::info!(
                    target: "PROFILE",
                    "Automatic low-map genus fallback mode DISABLED: reference database appears to be subsampled (first @SQ name has 6 `|`-separated fields)."
                );
            }
        }

        Ok(SamData {
            read_alignments,
            asm_lengths,
            min_breadth,
            min_chunk_breadth: min_chunk_breadth.unwrap_or(0.0),
            subsample_scaling,
            total_reads: num_reads,
            aligned_reads,
            mapping_rate,
            fallback_mode,
            species_normalize,
            species_to_genus,
            skipped_reference_at_load: num_skipped_reference,
        })
    }

    /// Estimate the minimum chunk breadth threshold from read count.
    /// The threshold calculation was determined based on the CAMI HMP toy dataset
    /// # Illumina formula
    ///
    /// - ≤100K reads: linear scaling (reads / 1M)
    /// - 100K-1M reads: 0.11 * reads_in_M + 0.09
    /// - >1M reads: 0.2 + 0.382 * log10(reads_in_M), capped at 0.95
    ///
    /// # Long reads (Nanopore/PacBio)
    ///
    /// - ≤5K reads: 0.0 (no filtering)
    /// - 5K-1M reads: 0.3
    /// - >1M reads: 0.5
    fn estimate_chunk_breadth(&self, num_reads: usize) -> f64 {
        // Low-biomass mode disables chunk breadth filtering
        if self.cfg.lowbiomass {
            return 0.0;
        }

        if self.cfg.sequencer.eq_ignore_ascii_case("illumina") {
            // Illumina: more reads -> stricter threshold
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

        // Long reads: simpler thresholds due to lower read counts
        if num_reads <= 5_000 {
            0.0
        } else if num_reads <= 1_000_000 {
            0.3
        } else {
            0.5
        }
    }

    /// Run the main profiling algorithm on loaded alignment data.
    ///
    /// Steps:
    /// 1. Assign taxonomy to reads based on alignment quality
    /// 2. Compute coverage metrics for each taxon
    /// 3. Resolve ambiguous multi-mapping reads with EM
    /// 4. Filter taxa based on coverage statistics
    /// 5. Write output files
    fn profile(&self, sam_data: SamData) -> Result<()> {
        let SamData {
            read_alignments,
            fallback_mode,
            species_to_genus,
            mapping_rate,
            total_reads,
            aligned_reads,
            skipped_reference_at_load,
            ..
        } = &sam_data;
        let verbose_logging = self.cfg.verbose || self.cfg.very_verbose;
        if verbose_logging {
            log::info!(
                target: "PROFILE",
                "Assigning taxonomy for reads with {} threads ...",
                self.cfg.threads
            );
            log::info!(
                target: "PROFILE",
                "Run mapping rate carried into assignment: {:.5} ({}/{}). Genus fallback = {}.",
                mapping_rate,
                aligned_reads,
                total_reads,
                if *fallback_mode { "enabled" } else { "disabled" },
            );
        }
        let mut assignment = self.assign_taxonomy(
            read_alignments,
            *fallback_mode,
            Arc::clone(species_to_genus),
        )?;
        // The `skipped_reference` diagnostic now accumulates two
        // streams: parser-time rejections (unparsable species taxid /
        // missing assembly length; moved here from `process_batch`)
        // plus any future process_batch-time rejections (currently
        // none).
        assignment.filter_stats.skipped_reference += *skipped_reference_at_load;
        if verbose_logging {
            log::info!(target: "PROFILE", "Assigning taxonomy finished!");
        }
        let candidate_reads = assignment.read_assignments.len();
        if verbose_logging {
            log::info!(
                target: "PROFILE",
                "Reads with candidate taxon assignments before coverage filtering: {}",
                candidate_reads
            );
            log::info!(
                target: "PROFILE",
                "Profiling taxonomy abundance ..."
            );
        }
        let final_classified = self.write_outputs(sam_data, assignment)?;
        log::info!(
            target: "PROFILE",
            "Number of taxonomically classified reads: {}",
            final_classified
        );
        log::info!(target: "PROFILE", "Taxonomy profiling finished.");
        Ok(())
    }

    /// Assign taxonomy to reads based on alignment quality.
    ///
    /// Uses rayon for parallel processing of read batches. Each read's
    /// alignments are filtered by identity, mapped length, and fraction
    /// thresholds, then the best alignment per species is selected.
    ///
    /// # Algorithm
    ///
    /// For each read:
    /// 1. Filter alignments by quality thresholds
    /// 2. Parse reference names to extract taxonomy info
    /// 3. For each species, keep only the highest-identity alignment
    /// 4. Record coverage intervals for breadth calculation
    /// 5. Store taxonomic assignment (may be multi-species)
    ///
    /// # Returns
    ///
    /// AssignmentResult containing per-read assignments and coverage intervals
    fn assign_taxonomy(
        &self,
        read_alignments: &[ReadBundle],
        fallback_mode: bool,
        species_to_genus: Arc<HashMap<u32, u32>>,
    ) -> Result<AssignmentResult> {
        let batch_size = self.cfg.batch_size.max(1);
        let cfg = Arc::new(self.cfg.clone());

        // Create rayon thread pool for parallel batch processing
        let pool = ThreadPoolBuilder::new()
            .num_threads(self.cfg.threads.max(1))
            .build()?;

        // Setup progress bar
        let pb = ProgressBar::new(read_alignments.len() as u64);
        pb.set_draw_target(ProgressDrawTarget::stdout());
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} {msg} {bar:40.cyan/blue} {pos}/{len} [{elapsed_precise}<{eta_precise}]",
            )?
            .progress_chars("#>-"),
        );
        pb.set_message("Assigning taxonomy");
        let progress = pb.clone();

        // Process read batches in parallel
        let mut results: Vec<AssignmentResult> = pool.install(|| {
            read_alignments
                .par_chunks(batch_size)
                .map_with(progress, |progress_bar, chunk| {
                    let result =
                        process_batch(chunk, &cfg, fallback_mode, species_to_genus.as_ref());
                    progress_bar.inc(chunk.len() as u64);
                    result
                })
                .collect()
        });
        pb.finish_and_clear();

        // Aggregate results from all batches
        let mut final_result = AssignmentResult::default();
        for result in results.drain(..) {
            // Merge read assignments
            final_result
                .read_assignments
                .extend(result.read_assignments.into_iter());

            // Merge coverage intervals
            for (key, intervals) in result.acc_intervals {
                final_result
                    .acc_intervals
                    .entry(key)
                    .or_default()
                    .extend(intervals);
            }
            for (key, intervals) in result.acc_fixed_chunk_intervals {
                final_result
                    .acc_fixed_chunk_intervals
                    .entry(key)
                    .or_default()
                    .extend(intervals);
            }
            for (key, intervals) in result.acc_flex_chunk_intervals {
                final_result
                    .acc_flex_chunk_intervals
                    .entry(key)
                    .or_default()
                    .extend(intervals);
            }

            // Merge filter statistics
            final_result.filter_stats.merge(result.filter_stats);
        }
        Ok(final_result)
    }

    /// Write profiling results to output files.
    ///
    /// This is a large function that:
    /// 1. Computes final coverage metrics per species
    /// 2. Performs EM algorithm for ambiguous read redistribution
    /// 3. Calculates abundance estimates and statistical p-values
    /// 4. Filters taxa based on coverage criteria
    /// 5. Writes profile.txt and classify.txt output files
    ///
    /// # Output Files
    ///
    /// - `<prefix>.profile.txt`: Filtered abundance profile with statistics
    /// - `<prefix>.rprofile.txt`: Raw (unfiltered) profile (if keep_raw=true)
    /// - `<prefix>.classify.txt`: Per-read classification assignments
    ///
    /// # EM Algorithm
    ///
    /// For multi-mapping reads, the EM algorithm iteratively:
    /// 1. E-step: Redistribute reads based on current abundance estimates
    /// 2. M-step: Recalculate abundances from redistributed reads
    /// 3. Converge when relative change < 1e-5
    ///
    /// # Returns
    ///
    /// Number of taxonomically classified reads
    fn write_outputs(&self, sam_data: SamData, assignment: AssignmentResult) -> Result<usize> {
        let SamData {
            asm_lengths,
            min_breadth,
            min_chunk_breadth,
            subsample_scaling,
            ..
        } = sam_data;
        let verbose_logging = self.cfg.verbose || self.cfg.very_verbose;

        let mut final_classified_reads = 0usize;
        let mut species_coverage: HashMap<u32, SpeciesCoverage> = HashMap::new();
        for (key, intervals) in &assignment.acc_intervals {
            let merged = merge_intervals(intervals);
            let merged_fixed = merge_intervals(
                assignment
                    .acc_fixed_chunk_intervals
                    .get(key)
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]),
            );
            let merged_flex = merge_intervals(
                assignment
                    .acc_flex_chunk_intervals
                    .get(key)
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]),
            );
            let asm_len = asm_lengths.get(&key.asm).copied().unwrap_or(1.0).max(1.0);
            let total_span: usize = merged.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let total_fixed: usize = merged_fixed.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let total_flex: usize = merged_flex.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let coverage = species_coverage.entry(key.taxid).or_default();
            coverage.breadth += total_span as f64 / asm_len;
            coverage.fixed_chunk += total_fixed as f64 / asm_len;
            coverage.flex_chunk += total_flex as f64 / asm_len;
        }

        let classification_path = format!("{}.classify.txt", self.cfg.outprefix);
        let mut classification_writer = if subsample_scaling.is_empty() {
            Some(BufWriter::new(
                File::create(&classification_path)
                    .with_context(|| format!("failed to create {}", classification_path))?,
            ))
        } else {
            None
        };

        let mut unambiguous: HashMap<u32, (f64, f64, f64)> = HashMap::new();
        let mut ambiguous_counts: HashMap<Vec<u32>, (f64, Vec<f64>, Vec<f64>)> = HashMap::new();
        let mut ambiguous_reads: HashMap<Vec<u32>, Vec<String>> = HashMap::new();

        for assignment in &assignment.read_assignments {
            let mut final_taxids = Vec::new();
            let mut final_depths = Vec::new();
            let mut final_maplens = Vec::new();
            for ((taxid, depth), maplen) in assignment
                .taxids
                .iter()
                .zip(&assignment.depths)
                .zip(&assignment.mapped_lens)
            {
                let coverage = species_coverage.get(taxid).cloned().unwrap_or_default();
                if coverage.breadth >= min_breadth {
                    final_taxids.push(*taxid);
                    final_depths.push(*depth);
                    final_maplens.push(*maplen);
                }
            }
            if !final_taxids.is_empty() {
                final_classified_reads += 1;
            }
            match final_taxids.len() {
                0 => {
                    continue;
                }
                1 => {
                    let taxid = final_taxids[0];
                    let entry = unambiguous.entry(taxid).or_insert((0.0, 0.0, 0.0));
                    entry.0 += 1.0;
                    entry.1 += final_depths[0];
                    entry.2 += final_maplens[0];
                    write_optional_line(
                        &mut classification_writer,
                        format!(
                            "{}\t{}\t{}\t{}\t1",
                            assignment.read,
                            taxid,
                            assignment.rank.as_str(),
                            taxid
                        ),
                    )?;
                }
                _ => {
                    ambiguous_reads
                        .entry(final_taxids.clone())
                        .or_default()
                        .push(assignment.read.clone());
                    ambiguous_counts
                        .entry(final_taxids.clone())
                        .and_modify(|entry| {
                            entry.0 += 1.0;
                            accumulate(&mut entry.1, &final_depths);
                            accumulate(&mut entry.2, &final_maplens);
                        })
                        .or_insert_with(|| (1.0, final_depths.clone(), final_maplens.clone()));
                }
            }
        }

        let mut taxids_with_lca: Vec<Vec<u32>> = Vec::new();
        let mut taxids_with_any_ambiguous: Vec<Vec<u32>> = Vec::new();
        let mut ambiguous_removed: HashMap<Vec<u32>, (f64, Vec<f64>, Vec<f64>)> = HashMap::new();

        for (taxids, (mapcount, depths, maplens)) in ambiguous_counts.clone() {
            let presence: Vec<bool> = taxids
                .iter()
                .map(|taxid| unambiguous.contains_key(taxid))
                .collect();
            if presence.iter().all(|&p| p) {
                continue;
            } else if presence.iter().any(|&p| p) {
                taxids_with_any_ambiguous.push(taxids.clone());
                let mut filtered_taxids = Vec::new();
                let mut filtered_depths = Vec::new();
                let mut filtered_maplens = Vec::new();
                for (idx, taxid) in taxids.iter().enumerate() {
                    if presence[idx] {
                        filtered_taxids.push(*taxid);
                        filtered_depths.push(depths[idx]);
                        filtered_maplens.push(maplens[idx]);
                    }
                }
                if filtered_taxids.len() == 1 {
                    let taxid = filtered_taxids[0];
                    let entry = unambiguous.entry(taxid).or_insert((0.0, 0.0, 0.0));
                    entry.0 += mapcount;
                    entry.1 += filtered_depths[0];
                    entry.2 += filtered_maplens[0];
                } else if !filtered_taxids.is_empty() {
                    ambiguous_removed
                        .entry(filtered_taxids.clone())
                        .and_modify(|entry| {
                            entry.0 += mapcount;
                            accumulate(&mut entry.1, &filtered_depths);
                            accumulate(&mut entry.2, &filtered_maplens);
                        })
                        .or_insert_with(|| {
                            (mapcount, filtered_depths.clone(), filtered_maplens.clone())
                        });
                }
            } else {
                taxids_with_lca.push(taxids.clone());
                let taxid_set: HashSet<u32> = taxids.iter().copied().collect();
                if let Some((lca_taxid, _lca_name, lca_rank)) =
                    self.taxonomy.get_majority_lca(&taxid_set, 0.7)
                {
                    let lca_breadth = taxids
                        .iter()
                        .filter_map(|taxid| species_coverage.get(taxid))
                        .map(|cov| cov.breadth)
                        .fold(0.0, f64::max);
                    let lca_fixed = taxids
                        .iter()
                        .filter_map(|taxid| species_coverage.get(taxid))
                        .map(|cov| cov.fixed_chunk)
                        .fold(0.0, f64::max);
                    let lca_flex = taxids
                        .iter()
                        .filter_map(|taxid| species_coverage.get(taxid))
                        .map(|cov| cov.flex_chunk)
                        .fold(0.0, f64::max);
                    species_coverage.insert(
                        lca_taxid,
                        SpeciesCoverage {
                            breadth: lca_breadth,
                            fixed_chunk: lca_fixed,
                            flex_chunk: lca_flex,
                        },
                    );
                    let depth_mean = depths.iter().sum::<f64>() / depths.len() as f64;
                    let maplen_mean = maplens.iter().sum::<f64>() / maplens.len() as f64;
                    let entry = unambiguous.entry(lca_taxid).or_insert((0.0, 0.0, 0.0));
                    entry.0 += mapcount;
                    entry.1 += depth_mean;
                    entry.2 += maplen_mean;
                    if let Some(reads) = ambiguous_reads.get(&taxids) {
                        for read in reads {
                            write_optional_line(
                                &mut classification_writer,
                                format!(
                                    "{}\t{}\t{}\t{}\t1",
                                    read,
                                    lca_taxid,
                                    lca_rank,
                                    taxids
                                        .iter()
                                        .map(|t| t.to_string())
                                        .collect::<Vec<_>>()
                                        .join(";"),
                                ),
                            )?;
                        }
                    }
                }
            }
        }

        for taxids in taxids_with_lca {
            ambiguous_counts.remove(&taxids);
            ambiguous_reads.remove(&taxids);
        }
        for taxids in taxids_with_any_ambiguous {
            ambiguous_counts.remove(&taxids);
            ambiguous_reads.remove(&taxids);
        }
        for (taxids, values) in ambiguous_removed {
            ambiguous_counts
                .entry(taxids)
                .and_modify(|entry| {
                    entry.0 += values.0;
                    accumulate(&mut entry.1, &values.1);
                    accumulate(&mut entry.2, &values.2);
                })
                .or_insert(values);
        }

        let mut species_count_depth = unambiguous.clone();
        let mut ambiguous_fraction: HashMap<Vec<u32>, HashMap<u32, f64>> = HashMap::new();
        let mut iterations = 0;
        loop {
            let mut total = unambiguous.clone();
            for (taxids, (mapcount, depths, maplens)) in &ambiguous_counts {
                let sum_depth: f64 = taxids
                    .iter()
                    .map(|taxid| species_count_depth.get(taxid).map(|v| v.1).unwrap_or(0.0))
                    .sum();
                if sum_depth <= f64::EPSILON {
                    continue;
                }
                let mut paired: Vec<(u32, f64, f64, f64)> = taxids
                    .iter()
                    .enumerate()
                    .map(|(idx, taxid)| {
                        let fraction = species_count_depth
                            .get(taxid)
                            .map(|v| v.1 / sum_depth)
                            .unwrap_or(0.0);
                        (*taxid, depths[idx], maplens[idx], fraction)
                    })
                    .collect();
                paired.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap_or(Ordering::Equal));
                let entry = ambiguous_fraction.entry(taxids.clone()).or_default();
                entry.clear();
                let mut consumed = false;
                for (taxid, depth, maplen, fraction) in paired {
                    let total_entry = total.entry(taxid).or_insert((0.0, 0.0, 0.0));
                    if fraction >= 0.99 {
                        total_entry.0 += mapcount;
                        total_entry.1 += depth;
                        total_entry.2 += maplen;
                        entry.insert(taxid, 1.0);
                        consumed = true;
                        break;
                    }
                    total_entry.0 += mapcount * fraction;
                    total_entry.1 += depth * fraction;
                    total_entry.2 += maplen * fraction;
                    entry.insert(taxid, round5(fraction));
                }
                if !consumed {
                    normalize_map(entry);
                }
            }
            if verbose_logging && iterations % 5 == 0 {
                log::info!(
                    target: "PROFILE",
                    "Profiling in EM iteration {} ...",
                    iterations
                );
            }
            let numerator: f64 = total
                .iter()
                .map(|(taxid, (_, depth, _))| {
                    let prev = species_count_depth.get(taxid).map(|v| v.1).unwrap_or(0.0);
                    (depth - prev).abs()
                })
                .sum::<f64>();
            let denominator: f64 = species_count_depth
                .values()
                .map(|(_, depth, _)| *depth)
                .sum::<f64>()
                .max(1e-12);
            if numerator / denominator <= 1e-5 {
                species_count_depth = total;
                break;
            }
            species_count_depth = total;
            iterations += 1;
        }
        if verbose_logging {
            log::info!(
                target: "PROFILE",
                "Finished profiling with {} EM iterations",
                iterations
            );
        }

        for (taxids, reads) in &ambiguous_reads {
            if let Some(fractions) = ambiguous_fraction.get(taxids) {
                let mut sorted: Vec<(&u32, &f64)> = fractions.iter().collect();
                sorted.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(Ordering::Equal));
                let final_taxid = *sorted[0].0;
                let taxid_strings: Vec<String> =
                    sorted.iter().map(|(taxid, _)| taxid.to_string()).collect();
                let fraction_strings: Vec<String> = sorted
                    .iter()
                    .map(|(_, frac)| format!("{:.5}", frac))
                    .collect();
                // All taxids within a single ambiguous group share the
                // same rank by construction: `process_batch` emits each
                // `ReadAssignment` as either all-species or all-genus,
                // and species/genus taxid namespaces are disjoint. We
                // look up the rank from the taxonomy of the chosen
                // taxid as a robust source of truth.
                let rank_label = self
                    .taxonomy
                    .get_rank(final_taxid)
                    .unwrap_or("species")
                    .to_string();
                for read in reads {
                    write_optional_line(
                        &mut classification_writer,
                        format!(
                            "{}\t{}\t{}\t{}\t{}",
                            read,
                            final_taxid,
                            rank_label,
                            taxid_strings.join(";"),
                            fraction_strings.join(";"),
                        ),
                    )?;
                }
            }
        }

        if let Some(writer) = classification_writer.as_mut() {
            writer.flush()?;
        }

        let norm_factor: f64 = 0.01
            * species_count_depth
                .values()
                .map(|(_, depth, _)| *depth)
                .sum::<f64>();
        let num_taxa = species_count_depth.len();

        let mut out_taxa_list = Vec::new();
        let host_taxid = self.cfg.host.as_deref().and_then(|h| h.parse::<u32>().ok());
        let pathogen_map = if self.cfg.host.is_some() {
            self.cfg
                .pathogen_host
                .as_deref()
                .map(|path| load_pathogen_table(path))
                .transpose()?
        } else {
            None
        };
        let mut entries: Vec<(u32, (f64, f64, f64))> = species_count_depth
            .iter()
            .map(|(taxid, (count, depth, maplen))| (*taxid, (*count, *depth, *maplen)))
            .collect();
        entries.sort_by(|a, b| b.1 .1.partial_cmp(&a.1 .1).unwrap_or(Ordering::Equal));

        for (taxid, (mapcount, depth, maplen)) in entries {
            if mapcount < 0.5 {
                continue;
            }
            let abundance = if norm_factor <= f64::EPSILON {
                0.0
            } else {
                depth / norm_factor
            };
            let scaling = subsample_scaling.get(&taxid).copied().unwrap_or(1.0);
            let scaled_mapcount = mapcount * scaling;
            let rank = self
                .taxonomy
                .get_rank(taxid)
                .unwrap_or("no rank")
                .to_string();
            let taxname = self
                .taxonomy
                .get_name(taxid)
                .unwrap_or("Taxonomy deprecated");
            if taxname == "Taxonomy deprecated" {
                continue;
            }
            let coverage = species_coverage.entry(taxid).or_default().clone();
            let breadth = coverage.breadth.min(1.0);
            let fixed_chunk = coverage.fixed_chunk.min(1.0);
            let flex_chunk = coverage.flex_chunk.min(1.0);
            let expected_breadth = if mapcount == 0.0 {
                0.0
            } else {
                1.0 - (1.0 - depth / mapcount).powf(mapcount)
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
            let genome_len = if depth <= f64::EPSILON {
                0.0
            } else {
                maplen / depth
            };
            let num_flex_chunks = if genome_len <= 0.0 {
                1
            } else {
                let raw_chunks = (genome_len.sqrt() / 2.0).floor();
                if raw_chunks < 1.0 {
                    1
                } else {
                    raw_chunks as usize
                }
            };
            let expected_flex_chunk = if num_flex_chunks == 0 {
                0.0
            } else {
                1.0 - (1.0 - 1.0 / num_flex_chunks as f64).powf(mapcount)
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
                        let (report_it, host_names, host_taxids, diseases) =
                            self.annotate_host_disease(entry, host_taxid);
                        if !report_it {
                            continue;
                        }
                        out_taxa_list.push(vec![
                            taxname.to_string(),
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
                            host_names.unwrap_or_else(|| "NA".to_string()),
                            host_taxids.unwrap_or_else(|| "NA".to_string()),
                            diseases.unwrap_or_else(|| "NA".to_string()),
                        ]);
                        continue;
                    }
                }
                out_taxa_list.push(vec![
                    taxname.to_string(),
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
                    taxname.to_string(),
                    taxid.to_string(),
                    rank.clone(),
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
                let rank = entry[2].as_str();
                let read_count: f64 = entry[3].parse().unwrap_or(0.0);
                let breadth: f64 = entry[6].parse().unwrap_or(0.0);
                let expected_breadth: f64 = entry[7].parse().unwrap_or(0.0);
                let fixed_chunk: f64 = entry[9].parse().unwrap_or(0.0);
                let flex_chunk: f64 = entry[10].parse().unwrap_or(0.0);
                let expected_flex: f64 = entry[11].parse().unwrap_or(0.0);
                let cov_prob = entry[8].parse::<f64>().ok();
                let chunk_prob = entry[12].parse::<f64>().ok();
                let pvalue_ok = pvalue_filter_passes(rank, cov_prob, chunk_prob, self.cfg.strict);
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
                min_reads_ok
                    && fixed_chunk_ok
                    && ratio_filter_passes(
                        rank,
                        breadth,
                        expected_breadth,
                        flex_chunk,
                        expected_flex,
                        self.cfg.min_oebr,
                        self.cfg.min_coebr,
                    )
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

        let out_prefix = if self.cfg.host.is_some() {
            format!("{}.pathogen", self.cfg.outprefix)
        } else {
            self.cfg.outprefix.clone()
        };
        let profile_path = format!("{}.profile.txt", out_prefix);
        let mut profile_writer = WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(File::create(&profile_path)?);
        for row in &normalized {
            profile_writer.write_record(row)?;
        }
        profile_writer.flush()?;

        if self.cfg.keep_raw {
            let raw_path = format!("{}.rprofile.txt", out_prefix);
            let mut raw_writer = WriterBuilder::new()
                .has_headers(false)
                .delimiter(b'\t')
                .from_writer(File::create(&raw_path)?);
            for row in raw_taxa_list {
                raw_writer.write_record(row)?;
            }
            raw_writer.flush()?;
        }

        Ok(final_classified_reads)
    }

    /// Annotate a pathogen with host and disease information.
    ///
    /// Checks if the pathogen can infect the specified host by examining
    /// the taxonomy lineage of the pathogen's known hosts.
    ///
    /// # Arguments
    /// * `entry` - Pathogen metadata from the pathogen-host mapping file
    /// * `host_taxid` - NCBI taxonomy ID of the host organism
    ///
    /// # Returns
    ///
    /// Tuple of (should_report, host_names, host_taxids, diseases)
    /// - should_report: true if this pathogen can infect the specified host
    /// - host_names: Semicolon-separated host names (or None)
    /// - host_taxids: Semicolon-separated host taxids (or None)
    /// - diseases: Semicolon-separated disease names (or None)
    fn annotate_host_disease(
        &self,
        entry: &PathogenEntry,
        host_taxid: u32,
    ) -> (bool, Option<String>, Option<String>, Option<String>) {
        let mut host_taxids_set: HashSet<u32> = HashSet::new();
        let mut host_names_set: HashSet<String> = HashSet::new();
        let mut diseases_set: HashSet<String> = HashSet::new();
        if !entry.host_taxids.is_empty() && entry.host_taxids != "unknown" {
            for tax in entry.host_taxids.split(';') {
                if let Ok(val) = tax.parse::<u32>() {
                    host_taxids_set.insert(val);
                }
            }
        }
        if !entry.host_names.is_empty() && entry.host_names != "unknown" {
            for name in entry.host_names.split(';') {
                host_names_set.insert(name.to_string());
            }
        }
        if !entry.diseases.is_empty() && entry.diseases != "unknown" {
            for disease in entry.diseases.split(';') {
                diseases_set.insert(disease.to_string());
            }
        }
        let host_taxids = if host_taxids_set.is_empty() {
            None
        } else {
            Some(
                host_taxids_set
                    .iter()
                    .map(|t| t.to_string())
                    .collect::<Vec<_>>()
                    .join(";"),
            )
        };
        let host_names = if host_names_set.is_empty() {
            None
        } else {
            Some(host_names_set.iter().cloned().collect::<Vec<_>>().join(";"))
        };
        let diseases = if diseases_set.is_empty() {
            None
        } else {
            Some(diseases_set.iter().cloned().collect::<Vec<_>>().join(";"))
        };
        let mut report_it = false;
        if host_taxids_set.contains(&host_taxid) {
            report_it = true;
        } else if !host_taxids_set.is_empty() {
            let mut lineage: HashSet<u32> = HashSet::new();
            for taxid_val in &host_taxids_set {
                lineage.extend(
                    self.taxonomy
                        .get_parents(*taxid_val)
                        .into_iter()
                        .map(|(t, _, _)| t),
                );
            }
            if lineage.contains(&host_taxid) {
                report_it = true;
            }
        }
        (report_it, host_names, host_taxids, diseases)
    }
}

/// Process a batch of reads for taxonomy assignment.
///
/// For each read in the batch:
/// 1. Filter alignments by mapped_len, identity, and fraction thresholds
/// 2. Parse reference names to extract taxonomy information
/// 3. Keep the best (highest identity) alignment per species
/// 4. Record coverage intervals for breadth calculation
///
/// # Arguments
/// * `chunk` - Slice of ReadBundle to process
/// * `cfg` - Profiler configuration with filtering thresholds
/// * `asm_lengths` - Assembly lengths for depth normalization
///
/// # Returns
///
/// AssignmentResult with read assignments and coverage intervals for this batch
/// Per-taxid best-alignment bookkeeping used while scanning a read's
/// alignments. `best` holds `(identity, depth, mapped_len)`; `order` keeps
/// first-seen order for deterministic output; `intervals` captures the
/// coverage interval triple from the currently-best alignment.
type BestMap = HashMap<u32, (f64, f64, f64)>;
type IntervalMap = HashMap<u32, (AccKey, (usize, usize), (usize, usize), (usize, usize))>;

/// Process a batch of reads for rank-aware taxonomy assignment.
///
/// Per-read precedence:
/// 1. Evaluate species-qualified alignments (identity/mapped_len/fraction
///    thresholds from `cfg`).
/// 2. Only when no species hit qualifies AND `fallback_mode` is `true`
///    AND a genus taxid can be resolved for the alignment's species,
///    evaluate the genus fallback thresholds (`identity >= cfg.genus_identity`,
///    mapped length >= the platform-specific genus floor (50 for Illumina,
///    `cfg.mapped_len` otherwise), and fraction >=
///    [`GENUS_FALLBACK_MIN_FRACTION`]).
///
/// Coverage intervals are attributed to the assigned taxid (species or
/// genus). Reference header genus fields are *not* consulted; the genus
/// taxid comes from `species_to_genus`, which was precomputed from the
/// NCBI taxonomy lineage.
///
/// NOTE: Each alignment in `chunk` has already been parsed and
/// normalized by [`parse_alignment_line`] — its species taxid, asm id,
/// accession, assembly length, identity, and fraction are all pre-
/// computed. This function therefore does **no** string splitting,
/// **no** `HashMap<String, _>` lookups, and **no** CIGAR scanning in
/// the per-alignment hot loop. It only does arithmetic comparisons
/// against `cfg` thresholds plus a single `species_to_genus` lookup
/// when the fallback path fires.
fn process_batch(
    chunk: &[ReadBundle],
    cfg: &ProfilerConfig,
    fallback_mode: bool,
    species_to_genus: &HashMap<u32, u32>,
) -> AssignmentResult {
    let genus_fb_mapped_min = genus_fallback_min_mapped_len(cfg);
    let mut result = AssignmentResult::default();

    // Per-read bookkeeping hoisted OUTSIDE the read loop so we reuse
    // the allocation across all reads in this batch. A fresh
    // `HashMap::new()` every iteration would allocate/drop once per
    // read (6 maps × #reads); `.clear()` keeps capacity and skips the
    // allocator traffic entirely.
    let mut species_best: BestMap = HashMap::new();
    let mut species_order: Vec<u32> = Vec::new();
    let mut species_intervals: IntervalMap = HashMap::new();
    let mut genus_best: BestMap = HashMap::new();
    let mut genus_order: Vec<u32> = Vec::new();
    let mut genus_intervals: IntervalMap = HashMap::new();

    for read in chunk {
        species_best.clear();
        species_order.clear();
        species_intervals.clear();
        genus_best.clear();
        genus_order.clear();
        genus_intervals.clear();

        for aln in &read.alignments {
            // All the pre-computed bits are already validated
            // (asm_len > 0, species_taxid parseable) by the parser, so
            // we just read them as f64 / u32 fields here.
            let depth = aln.mapped_len as f64 / aln.asm_len;
            let identity = aln.identity;
            let fraction = aln.fraction;

            // --- Step 1: species criteria (unchanged thresholds) ---
            let species_qualified = aln.mapped_len >= cfg.mapped_len
                && identity >= cfg.identity
                && fraction >= cfg.fraction;

            if species_qualified {
                update_best(
                    &mut species_best,
                    &mut species_order,
                    &mut species_intervals,
                    aln.species_taxid,
                    identity,
                    depth,
                    aln,
                );
                // Even if this alignment also satisfies the fallback
                // thresholds, precedence says the species hit wins for
                // this read; we still keep parsing other alignments to
                // find the best species candidate.
                continue;
            }

            // --- Step 2: genus fallback (only when gate is active) ---
            if !fallback_mode {
                continue;
            }
            let fallback_qualified = identity >= cfg.genus_identity
                && aln.mapped_len >= genus_fb_mapped_min
                && fraction >= GENUS_FALLBACK_MIN_FRACTION;
            if !fallback_qualified {
                continue;
            }
            let genus_taxid = match species_to_genus.get(&aln.species_taxid).copied() {
                Some(g) => g,
                // No genus ancestor resolvable => skip fallback for this alignment.
                None => continue,
            };
            update_best(
                &mut genus_best,
                &mut genus_order,
                &mut genus_intervals,
                genus_taxid,
                identity,
                depth,
                aln,
            );
        }

        // Species-first precedence: commit species if any, else genus (gated).
        if !species_order.is_empty() {
            commit_assignment(
                &mut result,
                &read.read_name,
                AssignRank::Species,
                &species_order,
                &species_best,
                &mut species_intervals,
            );
        } else if fallback_mode && !genus_order.is_empty() {
            commit_assignment(
                &mut result,
                &read.read_name,
                AssignRank::Genus,
                &genus_order,
                &genus_best,
                &mut genus_intervals,
            );
        }
    }
    result
}

/// Update per-taxid best-alignment bookkeeping, keeping the alignment with
/// the highest identity as representative for coverage intervals.
///
/// `AccKey` construction is deferred into the branches that actually
/// update the representative — the two `String` clones (`asm`, `acc`)
/// are skipped entirely when the alignment doesn't beat the current
/// best for its taxid.
#[inline]
fn update_best(
    best: &mut BestMap,
    order: &mut Vec<u32>,
    intervals: &mut IntervalMap,
    taxid: u32,
    identity: f64,
    depth: f64,
    aln: &AlignmentEntry,
) {
    match best.entry(taxid) {
        Entry::Vacant(entry) => {
            order.push(taxid);
            entry.insert((identity, depth, aln.mapped_len as f64));
            let acc_key = AccKey {
                taxid,
                asm: aln.asm.clone(),
                acc: aln.acc.clone(),
            };
            intervals.insert(taxid, (acc_key, aln.span, aln.fixed_chunk, aln.flex_chunk));
        }
        Entry::Occupied(mut entry) => {
            if identity > entry.get().0 {
                *entry.get_mut() = (identity, depth, aln.mapped_len as f64);
                let acc_key = AccKey {
                    taxid,
                    asm: aln.asm.clone(),
                    acc: aln.acc.clone(),
                };
                intervals.insert(taxid, (acc_key, aln.span, aln.fixed_chunk, aln.flex_chunk));
            }
        }
    }
}

/// Commit the best-per-taxid results for one read into the shared
/// `AssignmentResult`: appends a `ReadAssignment` and flushes coverage
/// intervals for the chosen rank.
///
/// The caller passes borrowed views of the per-read bookkeeping
/// (`order`, `best`, `intervals`) because those are hoisted outside
/// the read loop in [`process_batch`] and reused across reads. The
/// intervals map is drained (`remove(&taxid)`) so the next read can
/// `.clear()` it cheaply. To keep the three interval maps keyed by
/// the same `AccKey`, we clone the key twice when inserting into the
/// first two maps and move it into the third — same allocation count
/// as the old code but one fewer `.clone()` (the third push consumes
/// the `AccKey` instead of cloning it).
fn commit_assignment(
    result: &mut AssignmentResult,
    read_name: &str,
    rank: AssignRank,
    order: &[u32],
    best: &BestMap,
    intervals: &mut IntervalMap,
) {
    let mut taxids = Vec::with_capacity(order.len());
    let mut depths = Vec::with_capacity(order.len());
    let mut maplens = Vec::with_capacity(order.len());
    for &taxid in order {
        if let Some((_, depth, maplen)) = best.get(&taxid) {
            taxids.push(taxid);
            depths.push(*depth);
            maplens.push(*maplen);
        }
        if let Some((key, span, fixed, flex)) = intervals.remove(&taxid) {
            result
                .acc_intervals
                .entry(key.clone())
                .or_default()
                .push(span);
            result
                .acc_fixed_chunk_intervals
                .entry(key.clone())
                .or_default()
                .push(fixed);
            result
                .acc_flex_chunk_intervals
                .entry(key)
                .or_default()
                .push(flex);
        }
    }
    if taxids.is_empty() {
        return;
    }
    result.read_assignments.push(ReadAssignment {
        read: read_name.to_string(),
        rank,
        taxids,
        depths,
        mapped_lens: maplens,
    });
}

/// Absorb results from a worker thread into the main aggregation structures.
///
/// Called by the main thread to integrate worker results:
/// - Adds alignments to the read_map (keyed by read name)
/// - Records all seen reads in all_reads set
/// - Updates alignment count and logs progress
///
/// # Arguments
/// * `result` - Worker result to absorb
/// * `read_map` - Main map of read_name -> alignments
/// * `all_reads` - Set of all read names seen
/// * `num_alignments` - Running count of alignments
/// * `num_skipped_reference` - Running count of alignments discarded at
///   parse time because their reference name didn't carry a parsable
///   species taxid / assembly length (moved here from the old
///   `process_batch` counter)
/// * `next_progress` - Next milestone for progress logging
fn absorb_worker_result(
    result: WorkerResult,
    read_map: &mut HashMap<String, Vec<AlignmentEntry>>,
    all_reads: &mut HashSet<String>,
    num_alignments: &mut usize,
    num_skipped_reference: &mut usize,
    next_progress: &mut usize,
) {
    all_reads.extend(result.reads_seen);
    for (read, entry) in result.alignments {
        read_map.entry(read).or_default().push(entry);
    }
    *num_skipped_reference += result.skipped_reference;
    if result.alignments_seen > 0 {
        *num_alignments += result.alignments_seen;
        while *num_alignments >= *next_progress {
            log::info!(
                target: "PROFILE",
                "{} alignments processed",
                *num_alignments
            );
            *next_progress += 1_000_000;
        }
    }
}

/// Parse a single SAM alignment line into structured data.
///
/// Extracts key fields from the tab-delimited SAM format:
/// - Read name (with optional pair suffix stripping)
/// - Reference name (contains taxonomy info)
/// - Position and CIGAR string
/// - Alignment coordinates and chunk intervals
///
/// # SAM Fields Used
/// 1. QNAME - Read name
/// 2. FLAG - Bitwise flags (paired, unmapped, etc.)
/// 3. RNAME - Reference name
/// 4. POS - 1-based leftmost position
/// 6. CIGAR - Alignment operations
/// 7. RNEXT - Mate reference ("=" if same)
/// 8. PNEXT - Mate position
/// 10. SEQ - Read sequence
///
/// # Returns
///
/// ParsedAlignment with read name and optional AlignmentEntry
/// Returns None if the line is malformed
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
    let pnext = iter
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

    // Helper to short-circuit with "record seen but not kept".
    let discard = |reason_skipped: bool, count_alignment: bool| ParsedAlignment {
        read_name: read_name.clone(),
        entry: None,
        count_alignment,
        skipped_reference: reason_skipped,
    };

    if ref_name == "*" || cigar == "*" {
        return Some(discard(false, false));
    }
    if ctx.is_paired && rnext != "=" {
        return Some(discard(false, false));
    }

    let ref_len = match ctx.ref_lengths.get(ref_name) {
        Some(len) => *len,
        None => return Some(discard(false, false)),
    };
    if ref_len == 0 {
        return Some(discard(false, false));
    }

    let mapped_len = read_seq.len();
    if mapped_len == 0 {
        return Some(discard(false, false));
    }

    // Parse the reference-name taxonomy fields ONCE up front. Every
    // failure here used to surface as a `skipped_reference` increment
    // inside `process_batch`; we preserve that diagnostic by routing
    // the count through `ParsedAlignment::skipped_reference`.
    let mut ref_parts = ref_name.split('|');
    let asm_str = match ref_parts.next() {
        Some(v) if !v.is_empty() => v,
        _ => return Some(discard(true, false)),
    };
    let _asm_taxid = ref_parts.next();
    let raw_species_taxid = match ref_parts.next().and_then(|s| s.parse::<u32>().ok()) {
        Some(v) => v,
        None => return Some(discard(true, false)),
    };
    let species_taxid = ctx
        .species_normalize
        .get(&raw_species_taxid)
        .copied()
        .unwrap_or(raw_species_taxid);
    let acc_str = ref_parts.next().unwrap_or("");

    let asm_len = match ctx.asm_lengths.get(asm_str) {
        Some(len) if *len > 0.0 => *len,
        _ => return Some(discard(true, false)),
    };

    // Pre-compute identity / fraction once per alignment — `process_batch`
    // just reads these as f64s instead of re-scanning the CIGAR string.
    let (identity, _matched_len, fraction) = cigar_to_identity(cigar);

    let flex_len = ctx
        .asm_flex_lengths
        .get(asm_str)
        .copied()
        .unwrap_or(ctx.fixed_flank_len);

    let pair_start = if ctx.is_paired {
        start_pos.min(pnext)
    } else {
        start_pos
    };
    let mut start = pair_start.min(ref_len.max(1));
    if start == 0 {
        start = 1;
    }
    let ref_end = ref_len.saturating_add(1);
    let mut end = start.saturating_add(mapped_len);
    if end > ref_end {
        end = ref_end;
    }
    let midpoint = start + ((end.saturating_sub(start)) / 2);
    let fixed_chunk_start = start.saturating_sub(ctx.fixed_flank_len).max(1);
    let fixed_chunk_end = end.saturating_add(ctx.fixed_flank_len).min(ref_end);
    let flex_chunk_start = midpoint.saturating_sub(flex_len).max(1);
    let flex_chunk_end = midpoint.saturating_add(flex_len).min(ref_end);

    let entry = AlignmentEntry {
        asm: asm_str.to_string(),
        acc: acc_str.to_string(),
        species_taxid,
        asm_len,
        identity,
        fraction,
        mapped_len,
        span: (start, end),
        fixed_chunk: (fixed_chunk_start, fixed_chunk_end),
        flex_chunk: (flex_chunk_start, flex_chunk_end),
    };
    let mut count_alignment = true;
    if ctx.is_paired {
        if (flag & 0x4) != 0 {
            count_alignment = false;
        }
        if (flag & 0x40) == 0 {
            return Some(ParsedAlignment {
                read_name,
                entry: None,
                count_alignment,
                skipped_reference: false,
            });
        }
    }

    Some(ParsedAlignment {
        read_name,
        entry: Some(entry),
        count_alignment,
        skipped_reference: false,
    })
}

/// Merge overlapping genomic intervals.
///
/// Takes a list of (start, end) coordinate pairs and returns a list
/// where overlapping intervals have been merged into single intervals.
///
/// Used to calculate breadth of coverage - the total span of unique
/// positions covered by alignments.
///
/// # Algorithm
///
/// 1. Sort intervals by start position
/// 2. Iterate through, extending current interval or starting new one
/// 3. Intervals overlap if interval.start < current.end
///
/// # Example
///
/// Input: [(1,5), (3,7), (10,15)]
/// Output: [(1,7), (10,15)]
fn merge_intervals(intervals: &[(usize, usize)]) -> Vec<(usize, usize)> {
    if intervals.len() <= 1 {
        return intervals.to_vec();
    }
    let mut sorted = intervals.to_vec();
    sorted.sort_by_key(|i| i.0);
    let mut merged = Vec::with_capacity(sorted.len());
    merged.push(sorted[0]);
    for interval in sorted.into_iter().skip(1) {
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

/// Calculate alignment identity and fraction from a CIGAR string.
///
/// Parses CIGAR operations to compute:
/// - Identity: matches / (alignment_length + gap_events)
/// - Matched length: number of exact matches (=)
/// - Fraction: mapped_length / full_read_length
///
/// # CIGAR Operations
///
/// - M: Alignment match (may include mismatches) - NOT counted as match
/// - =: Sequence match (exact) - counted as match
/// - X: Sequence mismatch - counted in alignment
/// - I: Insertion to reference - counted in mapped length
/// - D: Deletion from reference - counts as gap event
/// - N: Skipped region - counts as gap event
/// - S: Soft clipping - not counted
/// - H: Hard clipping - counted in full length only
///
/// # Returns
///
/// Tuple of (identity, matched_len, fraction)
// fn cigar_to_identity(cigar: &str) -> (f64, usize, f64) {
//     let mut operation_len: HashMap<char, usize> = HashMap::new();
//     let mut num_gap_events = 0usize;
//     for caps in CIGAR_PATTERN.find_iter(cigar) {
//         let token = caps.as_str();
//         let (len_part, op_part) = token.split_at(token.len() - 1);
//         if let Ok(count) = len_part.parse::<usize>() {
//             let op = op_part.chars().next().unwrap();
//             *operation_len.entry(op).or_insert(0) += count;
//             if !matches!(op, 'M' | '=' | 'X' | 'H') {
//                 num_gap_events += 1;
//             }
//         }
//     }
//     let matched_len = *operation_len.get(&'=').unwrap_or(&0);
//     let aligned_len = matched_len
//         + operation_len.get(&'M').copied().unwrap_or(0)
//         + operation_len.get(&'X').copied().unwrap_or(0);
//     let insert_len = operation_len.get(&'I').copied().unwrap_or(0);
//     let mapped_len = aligned_len + insert_len;
//     let full_len = mapped_len + operation_len.get(&'H').copied().unwrap_or(0);
//     let fraction = if full_len == 0 {
//         0.0
//     } else {
//         mapped_len as f64 / full_len as f64
//     };
//     let gap_compressed = aligned_len + num_gap_events;
//     let identity = if gap_compressed == 0 {
//         0.0
//     } else {
//         matched_len as f64 / gap_compressed as f64
//     };
//     (identity, matched_len, fraction)
// }

/// Compute (identity, matched_len, fraction) from a SAM CIGAR string.
///
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

// ============================================================================
// Genus-level fallback helpers
// ============================================================================
//
// Genus fallback can assign a read at genus rank when no species hit qualifies.
// Automatic activation is limited to low-map, non-pathogen Illumina runs on
// full databases; `--genus-fallback` enables it explicitly.

/// Minimum mapped length (inclusive) for a genus-fallback candidate on
/// **Illumina** (short-read) data. Long-read platforms use `cfg.mapped_len`
/// instead (see [`genus_fallback_min_mapped_len`]).
const GENUS_FALLBACK_MIN_MAPPED_LEN: usize = 50;
/// Minimum aligned fraction (inclusive) for a genus-fallback candidate.
const GENUS_FALLBACK_MIN_FRACTION: f64 = 0.60;

/// Mapped-length floor for genus-fallback candidates: Illumina uses a fixed
/// 50 bp minimum; all other sequencers use the same cutoff as species
/// (`cfg.mapped_len`, from defaults or `-m`).
#[inline]
fn genus_fallback_min_mapped_len(cfg: &ProfilerConfig) -> usize {
    if cfg.sequencer.eq_ignore_ascii_case("illumina") {
        GENUS_FALLBACK_MIN_MAPPED_LEN
    } else {
        cfg.mapped_len
    }
}

/// Return `true` when the automatic genus-fallback gate is active.
fn fallback_gate_active(cfg: &ProfilerConfig, mapping_rate: f64) -> bool {
    cfg.is_low_map_run(mapping_rate)
}

/// Return whether genus-level assignment is active for this run.
#[inline]
fn genus_fallback_mode_active(
    cfg: &ProfilerConfig,
    low_map_gate_active: bool,
    db_is_subsampled: bool,
) -> bool {
    cfg.genus_fallback || (low_map_gate_active && !db_is_subsampled)
}

/// Run-level decision for the read count fed into `estimate_chunk_breadth`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ChunkBreadthPlan {
    /// `true` when the unscaled basis is aligned reads; `false` for total reads.
    use_aligned_basis: bool,
    /// `true` when aligned basis was selected by the automatic low-map rule.
    by_aligned_auto: bool,
    /// `true` when subsampled fallback uses total reads for threshold estimation.
    force_total_for_subsampled_fallback: bool,
    /// Unscaled read count selected as the chunk-breadth basis.
    basis_reads: usize,
    /// Read count actually passed to `estimate_chunk_breadth`.
    estimate_reads: usize,
    /// Whether `basis_reads` was multiplied by `FALLBACK_CHUNK_BASIS_SCALE`.
    scaled: bool,
}

/// Choose the read-count basis used to estimate chunk breadth.
#[inline]
fn plan_chunk_breadth(
    cfg: &ProfilerConfig,
    mapping_rate: f64,
    total_reads: usize,
    aligned_reads: usize,
    db_is_subsampled: bool,
    low_map_gate_active: bool,
) -> ChunkBreadthPlan {
    let chunk_breadth_pinned = cfg.chunk_breadth.is_some();
    let force_total_for_subsampled_fallback =
        db_is_subsampled && cfg.genus_fallback && !chunk_breadth_pinned;
    let by_aligned_effective =
        cfg.resolve_by_aligned(mapping_rate) && !force_total_for_subsampled_fallback;
    let basis_reads = if by_aligned_effective {
        aligned_reads
    } else {
        total_reads
    };

    let low_map_fallback_active = low_map_gate_active && !db_is_subsampled;
    let scaled = !chunk_breadth_pinned && by_aligned_effective && low_map_fallback_active;
    let estimate_reads = fallback_scale_chunk_basis(basis_reads, scaled);

    ChunkBreadthPlan {
        use_aligned_basis: by_aligned_effective,
        by_aligned_auto: by_aligned_effective && !cfg.by_aligned && low_map_gate_active,
        force_total_for_subsampled_fallback,
        basis_reads,
        estimate_reads,
        scaled,
    }
}

/// Multiplier used for automatic low-map genus fallback.
const FALLBACK_CHUNK_BASIS_SCALE: f64 = 1.5;

/// Scale the read count used as input to `estimate_chunk_breadth`.
#[inline]
fn fallback_scale_chunk_basis(chunk_basis_reads: usize, scale: bool) -> usize {
    if scale {
        ((chunk_basis_reads as f64) * FALLBACK_CHUNK_BASIS_SCALE).round() as usize
    } else {
        chunk_basis_reads
    }
}

/// Return `true` when a reference name has the subsampled-DB header layout.
fn ref_name_is_subsampled(sn: &str) -> bool {
    sn.split('|').count() == 6
}

/// Derive the genus taxid for a species taxid using the NCBI taxonomy
/// lineage (no header parsing).
///
/// Returns `None` if the species has no ancestor at rank "genus".
fn resolve_genus_taxid(taxonomy: &Taxonomy, species_taxid: u32) -> Option<u32> {
    taxonomy
        .get_parents(species_taxid)
        .into_iter()
        .find(|(_, _, rank)| rank == "genus")
        .map(|(taxid, _, _)| taxid)
}

/// Resolve a reference taxid to its species ancestor when one is available.
fn resolve_species_taxid(taxonomy: &Taxonomy, taxid: u32) -> Option<u32> {
    if matches!(taxonomy.get_rank(taxid), Some("species")) {
        return Some(taxid);
    }
    taxonomy
        .get_parents(taxid)
        .into_iter()
        .find(|(_, _, rank)| rank == "species")
        .map(|(id, _, _)| id)
}

/// Rank-specific default ratio bounds.
fn rank_filter_bounds(rank: &str) -> (f64, f64) {
    if rank == "genus" {
        (0.65, 1.5)
    } else {
        (0.75, 1.5)
    }
}

/// Ratio filter shared by species and genus rows.
fn ratio_filter_passes(
    rank: &str,
    breadth: f64,
    expected_breadth: f64,
    flex_chunk: f64,
    expected_flex: f64,
    min_oebr: Option<f64>,
    min_coebr: Option<f64>,
) -> bool {
    if expected_breadth <= 0.0 || expected_flex <= 0.0 {
        return false;
    }
    let (default_lo, hi) = rank_filter_bounds(rank);
    let min_oebr = min_oebr.unwrap_or(default_lo);
    let min_coebr = min_coebr.unwrap_or(default_lo);
    let br_ratio = breadth / expected_breadth;
    let fl_ratio = flex_chunk / expected_flex;
    br_ratio.is_finite()
        && fl_ratio.is_finite()
        && min_oebr.is_finite()
        && min_coebr.is_finite()
        && br_ratio >= min_oebr
        && br_ratio <= hi
        && fl_ratio >= min_coebr
        && fl_ratio <= hi
}

const TAXON_PROB_CUTOFF: f64 = 1e-5;

/// Probability filter for final taxa rows.
fn pvalue_filter_passes(
    rank: &str,
    cov_prob: Option<f64>,
    chunk_prob: Option<f64>,
    strict: bool,
) -> bool {
    match (cov_prob, chunk_prob) {
        (Some(cov), Some(chunk)) => {
            if strict {
                cov >= TAXON_PROB_CUTOFF && chunk >= TAXON_PROB_CUTOFF
            } else {
                !(cov < TAXON_PROB_CUTOFF && chunk < TAXON_PROB_CUTOFF)
            }
        }
        (None, None) if strict && rank == "genus" => true,
        _ => !strict,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profiler::ProfileMode;

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
        // "80=5I5=10X5H": extended CIGAR with 1 insertion gap-event and a
        // hard-clip. In the gap-compressed form the denominator is
        // eq + x + gap_events = 85 + 10 + 1 = 96, so identity = 85/96.
        // aligned_query = eq + x + m + i_ = 85 + 10 + 0 + 5 = 100,
        // query_len = aligned_query + s = 100, and H is excluded from the
        // query length altogether => fraction = 100/100 = 1.0.
        let (identity, matched_len, fraction) = cigar_to_identity("80=5I5=10X5H");
        approx_equal(identity, 85.0 / 96.0);
        assert_eq!(matched_len, 85);
        approx_equal(fraction, 1.0);

        // Extended CIGAR with only matches / mismatches: identity = 50/100.
        let (identity, matched_len, fraction) = cigar_to_identity("50=50X");
        approx_equal(identity, 0.5);
        assert_eq!(matched_len, 50);
        approx_equal(fraction, 1.0);

        // Non-extended CIGAR: identity proxy is M / (M + gap_events) and
        // matched_len falls back to M when no '=' operator is present.
        // "100M" => 100/(100+0) = 1.0, matched_len = 100, fraction = 1.0.
        let (identity, matched_len, fraction) = cigar_to_identity("100M");
        approx_equal(identity, 1.0);
        assert_eq!(matched_len, 100);
        approx_equal(fraction, 1.0);
    }

    // ------------------------------------------------------------------
    // Helpers for fallback / rank-aware tests.
    // ------------------------------------------------------------------

    /// Minimal `ProfilerConfig` for unit tests. Callers mutate only the
    /// fields relevant to the scenario under test.
    fn make_cfg(sequencer: &str, host: Option<&str>) -> ProfilerConfig {
        ProfilerConfig {
            sam: String::new(),
            sequencer: sequencer.to_string(),
            batch_size: 1,
            is_paired: false,
            identity: 0.95,
            mapped_len: 50,
            breadth: Some(0.0),
            chunk_breadth: Some(0.0),
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            fraction: 0.6,
            lowbiomass: false,
            keep_raw: false,
            by_aligned: false,
            genus_identity: 0.80,
            low_map_rate_threshold: 0.30,
            genus_fallback: false,
            strict: false,
            pathogen_host: None,
            host: host.map(|s| s.to_string()),
            threads: 1,
            verbose: false,
            very_verbose: false,
            outprefix: "test".to_string(),
            mode: ProfileMode::Default,
            dmp_dir: None,
            strain: false,
        }
    }

    /// Intermediate synthetic alignment used inside tests.
    ///
    /// Production code no longer stores `ref_name`/`cigar` directly on
    /// `AlignmentEntry` (the SAM worker pre-parses those into
    /// `asm`/`acc`/`species_taxid`/`identity`/`fraction`). For test
    /// readability we still describe alignments as "this ref name,
    /// this CIGAR" and let [`run_single_with_normalize`] perform the
    /// exact same pre-parsing the worker would do.
    #[derive(Clone, Debug)]
    struct RawAln {
        ref_name: String,
        cigar: String,
        mapped_len: usize,
        span: (usize, usize),
    }

    /// Construct a synthetic alignment descriptor suitable for
    /// [`run_single`] / [`run_single_with_normalize`].
    ///
    /// `ref_name` must match the expected `asm|taxid|species_taxid|acc|...`
    /// layout so the profiler can extract the species taxid.
    fn aln_entry(ref_name: &str, cigar: &str, mapped_len: usize, span: (usize, usize)) -> RawAln {
        RawAln {
            ref_name: ref_name.to_string(),
            cigar: cigar.to_string(),
            mapped_len,
            span,
        }
    }

    /// Convert a [`RawAln`] into a production [`AlignmentEntry`] by
    /// running the exact pre-parsing that
    /// [`parse_alignment_line`] does in the worker. Returns `None` when
    /// the synthetic alignment would be rejected at parse time
    /// (unparsable species taxid, unknown assembly length, etc.) so
    /// tests can still cover those error paths.
    fn build_aln(
        raw: &RawAln,
        asm_lengths: &HashMap<String, f64>,
        species_normalize: &HashMap<u32, u32>,
    ) -> Option<AlignmentEntry> {
        let mut parts = raw.ref_name.split('|');
        let asm = parts.next()?.to_string();
        if asm.is_empty() {
            return None;
        }
        let _ = parts.next();
        let raw_taxid = parts.next()?.parse::<u32>().ok()?;
        let species_taxid = species_normalize
            .get(&raw_taxid)
            .copied()
            .unwrap_or(raw_taxid);
        let acc = parts.next().unwrap_or("").to_string();
        let asm_len = asm_lengths.get(asm.as_str()).copied().unwrap_or(0.0);
        if asm_len <= 0.0 {
            return None;
        }
        let (identity, _, fraction) = cigar_to_identity(&raw.cigar);
        Some(AlignmentEntry {
            asm,
            acc,
            species_taxid,
            asm_len,
            identity,
            fraction,
            mapped_len: raw.mapped_len,
            span: raw.span,
            fixed_chunk: raw.span,
            flex_chunk: raw.span,
        })
    }

    // ------------------------------------------------------------------
    // Fallback gate tests
    // ------------------------------------------------------------------

    #[test]
    fn fallback_gate_activates_only_below_mapping_rate_threshold() {
        let cfg = make_cfg("Illumina", None);
        assert!(fallback_gate_active(&cfg, 0.0));
        assert!(fallback_gate_active(&cfg, 0.29));
        assert!(fallback_gate_active(&cfg, 0.2999999));
        // Exactly at the 0.30 threshold is NOT active (strict `<`).
        assert!(!fallback_gate_active(&cfg, 0.30));
        assert!(!fallback_gate_active(&cfg, 0.50));
        assert!(!fallback_gate_active(&cfg, 1.0));
    }

    #[test]
    fn fallback_gate_disabled_when_host_is_set() {
        // Any host string (e.g. human 9606) disables fallback regardless
        // of sequencer or mapping rate.
        let cfg_host = make_cfg("Illumina", Some("9606"));
        assert!(!fallback_gate_active(&cfg_host, 0.0));
        assert!(!fallback_gate_active(&cfg_host, 0.1));
    }

    #[test]
    fn fallback_gate_disabled_for_non_illumina() {
        for seq in ["Nanopore", "PacBio", "assembly", "OTHER"] {
            let cfg = make_cfg(seq, None);
            assert!(
                !fallback_gate_active(&cfg, 0.0),
                "expected fallback disabled for sequencer {seq}"
            );
        }
    }

    #[test]
    fn ref_name_is_subsampled_checks_six_field_layout() {
        // 6 fields (subsampled DB) => true
        assert!(ref_name_is_subsampled(
            "GCF_000172695.2|2020311|2020311|NZ_CP021128.1|6817255|600000"
        ));
        // 5 fields (full DB) => false
        assert!(!ref_name_is_subsampled(
            "GCF_000172695.2|2020311|2020311|NZ_CP021128.1|6817255"
        ));
        // Other counts => false
        assert!(!ref_name_is_subsampled("asmA"));
        assert!(!ref_name_is_subsampled("asmA|1"));
        assert!(!ref_name_is_subsampled("asmA|1|2|3"));
        assert!(!ref_name_is_subsampled("asmA|1|2|3|4|5|6|7"));
        // Empty fields still count as separators.
        assert!(ref_name_is_subsampled("|||||"));
    }

    #[test]
    fn fallback_gate_respects_custom_low_map_rate_threshold() {
        // With the default 0.30 threshold, 0.40 is NOT low-map.
        let mut cfg = make_cfg("Illumina", None);
        assert!(!fallback_gate_active(&cfg, 0.40));
        // Raise the threshold: now 0.40 IS low-map (strict `<`).
        cfg.low_map_rate_threshold = 0.50;
        assert!(fallback_gate_active(&cfg, 0.40));
        assert!(fallback_gate_active(&cfg, 0.4999));
        assert!(!fallback_gate_active(&cfg, 0.50));
        // Lower the threshold: even 0.05 no longer triggers when = threshold.
        cfg.low_map_rate_threshold = 0.05;
        assert!(!fallback_gate_active(&cfg, 0.05));
        assert!(fallback_gate_active(&cfg, 0.0499));
    }

    #[test]
    fn resolve_by_aligned_honors_explicit_flag() {
        // --by-aligned flag present: aligned basis is forced regardless of
        // sequencer, host, or mapping rate.
        let mut cfg = make_cfg("Illumina", None);
        cfg.by_aligned = true;
        assert!(cfg.resolve_by_aligned(0.99));
        cfg.sequencer = "Nanopore".to_string();
        assert!(cfg.resolve_by_aligned(0.99));
        cfg.host = Some("9606".to_string());
        assert!(cfg.resolve_by_aligned(0.99));
    }

    fn pvalue_filter(rank: &str, cov: Option<f64>, chunk: Option<f64>, strict: bool) -> bool {
        pvalue_filter_passes(rank, cov, chunk, strict)
    }

    #[test]
    fn strict_filter_requires_both_probs_above_cutoff() {
        // Both probabilities at or above the 1e-5 cutoff => keep in both modes.
        assert!(pvalue_filter("species", Some(1e-5), Some(1e-3), false));
        assert!(pvalue_filter("species", Some(1e-5), Some(1e-3), true));

        // Exactly one below cutoff:
        //   non-strict (OR): kept
        //   strict    (AND): dropped
        assert!(pvalue_filter("species", Some(1e-6), Some(1e-3), false));
        assert!(!pvalue_filter("species", Some(1e-6), Some(1e-3), true));
        assert!(pvalue_filter("species", Some(1e-3), Some(1e-6), false));
        assert!(!pvalue_filter("species", Some(1e-3), Some(1e-6), true));

        // Both below cutoff => dropped in both modes.
        assert!(!pvalue_filter("species", Some(1e-6), Some(1e-7), false));
        assert!(!pvalue_filter("species", Some(1e-6), Some(1e-7), true));
    }

    #[test]
    fn strict_filter_handles_missing_probabilities() {
        // Non-strict mode keeps rows with missing probabilities.
        assert!(pvalue_filter("species", None, Some(1e-3), false));
        assert!(pvalue_filter("species", Some(1e-3), None, false));
        assert!(pvalue_filter("species", None, None, false));

        // Species rows have p-values, so strict mode drops rows with missing
        // probability evidence.
        assert!(!pvalue_filter("species", None, Some(1e-3), true));
        assert!(!pvalue_filter("species", Some(1e-3), None, true));
        assert!(!pvalue_filter("species", None, None, true));

        // Strict mode keeps genus rows only when both p-values are absent.
        assert!(pvalue_filter("genus", None, None, true));
        assert!(!pvalue_filter("genus", None, Some(1e-3), true));
        assert!(!pvalue_filter("genus", Some(1e-3), None, true));
    }

    #[test]
    fn resolve_by_aligned_auto_forces_for_low_map_illumina_non_pathogen() {
        // Default cfg has --by-aligned unset (false); auto-rule applies.
        let cfg = make_cfg("Illumina", None);
        assert!(!cfg.by_aligned);
        // Below default 0.30 threshold => aligned basis is auto-forced.
        assert!(cfg.resolve_by_aligned(0.10));
        // At or above threshold => total-read basis.
        assert!(!cfg.resolve_by_aligned(0.30));
        assert!(!cfg.resolve_by_aligned(0.50));
        // Pathogen mode on disables auto-force.
        let cfg_host = make_cfg("Illumina", Some("9606"));
        assert!(!cfg_host.resolve_by_aligned(0.10));
        // Non-Illumina disables auto-force.
        let cfg_np = make_cfg("Nanopore", None);
        assert!(!cfg_np.resolve_by_aligned(0.10));
    }

    #[test]
    fn fallback_scale_chunk_basis_applies_only_when_scale_true() {
        // Pass-through when scaling is off.
        assert_eq!(fallback_scale_chunk_basis(0, false), 0);
        assert_eq!(fallback_scale_chunk_basis(1, false), 1);
        assert_eq!(fallback_scale_chunk_basis(100_000, false), 100_000);

        // When scaling is on, count is multiplied by FALLBACK_CHUNK_BASIS_SCALE
        // and rounded to the nearest integer.
        assert_eq!(fallback_scale_chunk_basis(0, true), 0);
        assert_eq!(fallback_scale_chunk_basis(1, true), 2); // 1.5 -> 2
        assert_eq!(fallback_scale_chunk_basis(2, true), 3);
        assert_eq!(fallback_scale_chunk_basis(10, true), 15);
        assert_eq!(fallback_scale_chunk_basis(100, true), 150);
        assert_eq!(fallback_scale_chunk_basis(100_000, true), 150_000);
        // Round-half-to-even-ish: 3 * 1.5 = 4.5 => rounds to 5.
        assert_eq!(fallback_scale_chunk_basis(3, true), 5);
    }

    #[test]
    fn fallback_chunk_basis_scale_constant_is_one_and_a_half() {
        assert!((FALLBACK_CHUNK_BASIS_SCALE - 1.5).abs() < 1e-12);
    }

    #[test]
    fn chunk_breadth_plan_scales_only_low_map_aligned_basis() {
        let total_reads = 100_000usize;
        let aligned_reads = 10_000usize;

        // --by-aligned alone: aligned basis, same equation, no 1.5x.
        let mut cfg = make_cfg("Illumina", None);
        cfg.chunk_breadth = None;
        cfg.by_aligned = true;
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // Automatic low-map fallback on a full DB: aligned basis with 1.5x.
        let mut cfg = make_cfg("Illumina", None);
        cfg.chunk_breadth = None;
        let gate = fallback_gate_active(&cfg, 0.10);
        assert!(genus_fallback_mode_active(&cfg, gate, false));
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, false, gate);
        assert!(plan.use_aligned_basis);
        assert!(plan.by_aligned_auto);
        assert!(plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, 15_000);

        // Explicit --genus-fallback on non-low-map data without --by-aligned:
        // total basis, no scale.
        let mut cfg = make_cfg("Illumina", None);
        cfg.chunk_breadth = None;
        cfg.genus_fallback = true;
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(!plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, total_reads);
        assert_eq!(plan.estimate_reads, total_reads);

        // Explicit --genus-fallback plus --by-aligned on non-low-map data:
        // aligned basis, same equation, no scale.
        cfg.by_aligned = true;
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // Explicit --genus-fallback on a subsampled DB still enables fallback
        // assignment, but chunk-breadth estimation uses total reads.
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, true, false);
        assert!(genus_fallback_mode_active(&cfg, false, true));
        assert!(!plan.use_aligned_basis);
        assert!(plan.force_total_for_subsampled_fallback);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, total_reads);
        assert_eq!(plan.estimate_reads, total_reads);

        // A low-map subsampled run without explicit --genus-fallback keeps the
        // automatic fallback assignment disabled, so it must not scale.
        let mut cfg = make_cfg("Illumina", None);
        cfg.chunk_breadth = None;
        let gate = fallback_gate_active(&cfg, 0.10);
        assert!(!genus_fallback_mode_active(&cfg, gate, true));
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, true, gate);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // User-pinned --min-cbreadth is never rescaled.
        cfg.chunk_breadth = Some(0.25);
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, false, gate);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.estimate_reads, aligned_reads);
    }

    #[test]
    fn fallback_gate_active_for_illumina_single_and_paired() {
        // Single-end
        let mut cfg = make_cfg("Illumina", None);
        cfg.is_paired = false;
        assert!(fallback_gate_active(&cfg, 0.1));
        // Paired-end
        cfg.is_paired = true;
        assert!(fallback_gate_active(&cfg, 0.1));
        // Case-insensitive sequencer match
        cfg.sequencer = "ILLUMINA".to_string();
        assert!(fallback_gate_active(&cfg, 0.1));
        cfg.sequencer = "illumina".to_string();
        assert!(fallback_gate_active(&cfg, 0.1));
    }

    // ------------------------------------------------------------------
    // Rank filter bounds
    // ------------------------------------------------------------------

    #[test]
    fn rank_filter_bounds_match_spec() {
        assert_eq!(rank_filter_bounds("species"), (0.75, 1.5));
        assert_eq!(rank_filter_bounds("genus"), (0.65, 1.5));
        // Unknown ranks default to species-level stringency.
        assert_eq!(rank_filter_bounds("family"), (0.75, 1.5));
    }

    #[test]
    fn ratio_filter_genus_passes_at_066_species_fails_at_066() {
        // Both ratios = 0.66 => genus passes (>= 0.65), species fails (< 0.75).
        let br_obs = 0.66;
        let br_exp = 1.0;
        let fx_obs = 0.66;
        let fx_exp = 1.0;
        assert!(ratio_filter_passes(
            "genus", br_obs, br_exp, fx_obs, fx_exp, None, None
        ));
        assert!(!ratio_filter_passes(
            "species", br_obs, br_exp, fx_obs, fx_exp, None, None
        ));
    }

    #[test]
    fn ratio_filter_both_fail_above_15() {
        // Ratio slightly above 1.5 => both ranks fail.
        let br_obs = 1.51;
        let br_exp = 1.0;
        let fx_obs = 1.51;
        let fx_exp = 1.0;
        assert!(!ratio_filter_passes(
            "species", br_obs, br_exp, fx_obs, fx_exp, None, None
        ));
        assert!(!ratio_filter_passes(
            "genus", br_obs, br_exp, fx_obs, fx_exp, None, None
        ));
        // At exactly 1.5 both pass (inclusive upper bound per spec).
        assert!(ratio_filter_passes(
            "species", 1.5, 1.0, 1.5, 1.0, None, None
        ));
        assert!(ratio_filter_passes("genus", 1.5, 1.0, 1.5, 1.0, None, None));
    }

    #[test]
    fn ratio_filter_requires_positive_expectations() {
        // Non-positive expected values must fail regardless of observed.
        assert!(!ratio_filter_passes(
            "species", 1.0, 0.0, 1.0, 1.0, None, None
        ));
        assert!(!ratio_filter_passes(
            "genus", 1.0, 1.0, 1.0, 0.0, None, None
        ));
    }

    #[test]
    fn ratio_filter_honours_independent_lower_bound_overrides() {
        assert!(ratio_filter_passes(
            "species",
            0.80,
            1.0,
            0.78,
            1.0,
            Some(0.80),
            Some(0.78),
        ));
        assert!(!ratio_filter_passes(
            "species",
            0.79,
            1.0,
            0.78,
            1.0,
            Some(0.80),
            Some(0.78),
        ));
        assert!(!ratio_filter_passes(
            "species",
            0.80,
            1.0,
            0.77,
            1.0,
            Some(0.80),
            Some(0.78),
        ));
    }

    // ------------------------------------------------------------------
    // process_batch precedence & fallback thresholds
    // ------------------------------------------------------------------

    /// Helper: run `process_batch` on a single read and return the
    /// resulting `ReadAssignment` (if any). Uses an empty
    /// species-normalization map (no sub-species remapping).
    fn run_single(
        cfg: &ProfilerConfig,
        asm_lengths: &HashMap<String, f64>,
        s2g: &HashMap<u32, u32>,
        fallback_mode: bool,
        read_name: &str,
        alns: Vec<RawAln>,
    ) -> Option<ReadAssignment> {
        run_single_with_normalize(
            cfg,
            asm_lengths,
            &HashMap::new(),
            s2g,
            fallback_mode,
            read_name,
            alns,
        )
    }

    /// Helper that additionally accepts a `species_normalize` map so
    /// tests can exercise the sub-species → species remapping path.
    fn run_single_with_normalize(
        cfg: &ProfilerConfig,
        asm_lengths: &HashMap<String, f64>,
        species_normalize: &HashMap<u32, u32>,
        s2g: &HashMap<u32, u32>,
        fallback_mode: bool,
        read_name: &str,
        alns: Vec<RawAln>,
    ) -> Option<ReadAssignment> {
        // Pre-parse each synthetic alignment exactly like the SAM
        // worker would — rejects are filtered out, mirroring the real
        // pipeline where unparsable records never reach `process_batch`.
        let alignments: Vec<AlignmentEntry> = alns
            .iter()
            .filter_map(|raw| build_aln(raw, asm_lengths, species_normalize))
            .collect();
        let bundle = ReadBundle {
            read_name: read_name.to_string(),
            alignments,
        };
        let mut result = process_batch(&[bundle], cfg, fallback_mode, s2g);
        result.read_assignments.pop()
    }

    fn mk_asm_map() -> HashMap<String, f64> {
        let mut m = HashMap::new();
        m.insert("asmA".to_string(), 1000.0);
        m.insert("asmB".to_string(), 1000.0);
        m
    }

    fn mk_genus_map() -> HashMap<u32, u32> {
        let mut m = HashMap::new();
        m.insert(101, 5001); // species 101 -> genus 5001
        m.insert(102, 5001); // species 102 -> genus 5001
        m.insert(103, 5002); // species 103 -> genus 5002
        m
    }

    /// CIGAR that yields identity = 1.0 (all '='): used for species-qualified hits.
    const PERFECT_CIGAR: &str = "100=";
    /// CIGAR that yields identity = 0.92, fraction = 1.0: used for genus fallback.
    const NEAR_MISS_CIGAR: &str = "92=8X";

    // ------------------------------------------------------------------
    // Sub-species -> species normalization
    // ------------------------------------------------------------------

    #[test]
    fn subspecies_taxid_is_normalized_to_species() {
        // Raw parts[2] = 9001 is a sub-species of species 101.
        // With `species_normalize = { 9001 -> 101 }` the read must be
        // assigned to 101, not 9001. This is independent of fallback.
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        let mut norm = HashMap::new();
        norm.insert(9001u32, 101u32);
        let alns = vec![aln_entry(
            "asmA|9001|9001|accA|1000",
            PERFECT_CIGAR,
            100,
            (1, 100),
        )];
        // With fallback OFF (previous-logic path): still normalizes.
        let no_fb = run_single_with_normalize(
            &cfg,
            &asm,
            &norm,
            &s2g,
            /*fallback*/ false,
            "r1",
            alns.clone(),
        )
        .expect("species assignment");
        assert_eq!(no_fb.rank, AssignRank::Species);
        assert_eq!(no_fb.taxids, vec![101]);

        // With fallback ON: species precedence still picks species 101
        // (not the raw sub-species id and not genus).
        let with_fb = run_single_with_normalize(&cfg, &asm, &norm, &s2g, true, "r2", alns)
            .expect("species assignment");
        assert_eq!(with_fb.rank, AssignRank::Species);
        assert_eq!(with_fb.taxids, vec![101]);
    }

    #[test]
    fn multiple_subspecies_of_same_species_merge_into_one() {
        // Two alignments hit different sub-species (9001 and 9002) of
        // the same species 101. They must collapse into a single
        // species-level assignment with taxid 101.
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        let mut norm = HashMap::new();
        norm.insert(9001u32, 101u32);
        norm.insert(9002u32, 101u32);
        let alns = vec![
            aln_entry("asmA|9001|9001|accA|1000", PERFECT_CIGAR, 100, (1, 100)),
            aln_entry("asmB|9002|9002|accB|1000", PERFECT_CIGAR, 100, (1, 100)),
        ];
        let assignment = run_single_with_normalize(&cfg, &asm, &norm, &s2g, false, "r1", alns)
            .expect("species assignment");
        assert_eq!(assignment.rank, AssignRank::Species);
        // One unique taxid (101), not [9001, 9002] and not [101, 101].
        assert_eq!(assignment.taxids, vec![101]);
    }

    #[test]
    fn missing_normalize_entry_keeps_raw_taxid() {
        // When the raw taxid has no entry in `species_normalize` the
        // raw id is used as-is. This preserves the original behaviour
        // for references whose third field is already a species taxid
        // (or is above species with no species ancestor).
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        let norm: HashMap<u32, u32> = HashMap::new();
        let alns = vec![aln_entry(
            "asmA|1|101|accA|1000",
            PERFECT_CIGAR,
            100,
            (1, 100),
        )];
        let assignment = run_single_with_normalize(&cfg, &asm, &norm, &s2g, false, "r1", alns)
            .expect("species assignment");
        assert_eq!(assignment.rank, AssignRank::Species);
        assert_eq!(assignment.taxids, vec![101]);
    }

    #[test]
    fn subspecies_normalization_applies_to_genus_fallback_key() {
        // Only a near-miss alignment (genus-fallback candidate) against
        // sub-species 9001, which normalizes to species 101 whose genus
        // is 5001. The genus fallback path must therefore emit genus
        // 5001 (derived via the normalized species), not whatever
        // s2g[9001] would have been (s2g is keyed by species here).
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map(); // only 101, 102, 103 are species keys
        let mut norm = HashMap::new();
        norm.insert(9001u32, 101u32);
        let alns = vec![aln_entry(
            "asmA|9001|9001|accA|1000",
            NEAR_MISS_CIGAR,
            100,
            (1, 100),
        )];
        let assignment = run_single_with_normalize(&cfg, &asm, &norm, &s2g, true, "r1", alns)
            .expect("genus assignment");
        assert_eq!(assignment.rank, AssignRank::Genus);
        assert_eq!(assignment.taxids, vec![5001]);
    }

    #[test]
    fn species_first_precedence_overrides_genus_fallback() {
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // Two alignments: first qualifies as species (perfect), second only
        // as genus fallback. Species must win even with fallback_mode=true.
        let alns = vec![
            aln_entry("asmA|1|101|accA|1000", PERFECT_CIGAR, 100, (1, 100)),
            aln_entry("asmB|1|103|accB|1000", NEAR_MISS_CIGAR, 100, (1, 100)),
        ];
        let assignment =
            run_single(&cfg, &asm, &s2g, /*fallback*/ true, "r1", alns).expect("assignment");
        assert_eq!(assignment.rank, AssignRank::Species);
        assert_eq!(assignment.taxids, vec![101]);
    }

    #[test]
    fn genus_fallback_only_when_no_species_qualifies() {
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // Only a near-miss alignment: passes genus thresholds, not species.
        let alns = vec![aln_entry(
            "asmA|1|101|accA|1000",
            NEAR_MISS_CIGAR,
            100,
            (1, 100),
        )];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns).expect("assignment");
        assert_eq!(assignment.rank, AssignRank::Genus);
        assert_eq!(assignment.taxids, vec![5001]);
        assert_eq!(assignment.rank.as_str(), "genus");
    }

    #[test]
    fn genus_fallback_skipped_when_gate_inactive() {
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // Same near-miss alignment; fallback_mode = false => no assignment.
        let alns = vec![aln_entry(
            "asmA|1|101|accA|1000",
            NEAR_MISS_CIGAR,
            100,
            (1, 100),
        )];
        let assignment = run_single(&cfg, &asm, &s2g, /*fallback*/ false, "r1", alns);
        assert!(
            assignment.is_none(),
            "expected no assignment when fallback is disabled"
        );
    }

    #[test]
    fn genus_fallback_uses_fixed_mapped_len_and_fraction_constants() {
        // Even when the user sets cfg.mapped_len and cfg.fraction higher
        // than the fallback constants, the genus path must still accept
        // mapped_len >= 50 and fraction >= 0.60 alignments.
        let mut cfg = make_cfg("Illumina", None);
        cfg.mapped_len = 200; // aggressive user threshold
        cfg.fraction = 0.95; //  aggressive user threshold
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // Identity 0.92, fraction 1.0, mapped_len 60 => species disqualified
        // (mapped_len < cfg.mapped_len), but genus-qualified (mapped_len >= 50,
        // fraction >= 0.6, identity >= cfg.genus_identity).
        let alns = vec![aln_entry(
            "asmA|1|101|accA|1000",
            NEAR_MISS_CIGAR,
            60,
            (1, 60),
        )];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns).expect("assignment");
        assert_eq!(assignment.rank, AssignRank::Genus);
        assert_eq!(assignment.taxids, vec![5001]);
    }

    #[test]
    fn genus_fallback_accepts_high_identity_alignment_that_fails_species_fraction() {
        let mut cfg = make_cfg("Illumina", None);
        cfg.identity = 0.95;
        cfg.fraction = 0.80;
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // identity = 1.0 passes species identity, but fraction = 0.70 fails
        // the species fraction cutoff while still passing the genus fallback
        // fraction cutoff. High identity must not be rejected by an upper-bound
        // species-identity check in the fallback path.
        let alns = vec![aln_entry("asmA|1|101|accA|1000", "70=30S", 70, (1, 70))];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns)
            .expect("high-identity partial alignment should fall back to genus");
        assert_eq!(assignment.rank, AssignRank::Genus);
        assert_eq!(assignment.taxids, vec![5001]);
    }

    #[test]
    fn genus_fallback_rejected_when_mapped_len_below_50() {
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // mapped_len = 49 violates the fixed 50 floor for genus fallback.
        let alns = vec![aln_entry(
            "asmA|1|101|accA|1000",
            NEAR_MISS_CIGAR,
            49,
            (1, 49),
        )];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns);
        assert!(
            assignment.is_none(),
            "expected no assignment when mapped_len < 50"
        );
    }

    #[test]
    fn genus_fallback_long_reads_use_cfg_mapped_len() {
        let mut cfg = make_cfg("Nanopore", None);
        cfg.mapped_len = 100;
        cfg.identity = 0.95;
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        let alns_short = vec![aln_entry(
            "asmA|1|101|accA|1000",
            NEAR_MISS_CIGAR,
            80,
            (1, 80),
        )];
        assert!(
            run_single(&cfg, &asm, &s2g, true, "r1", alns_short).is_none(),
            "long reads: genus fallback must require cfg.mapped_len, not the Illumina 50 bp floor"
        );
        let alns_ok = vec![aln_entry(
            "asmA|1|101|accA|1000",
            NEAR_MISS_CIGAR,
            100,
            (1, 100),
        )];
        let assignment =
            run_single(&cfg, &asm, &s2g, true, "r1", alns_ok).expect("genus assignment");
        assert_eq!(assignment.rank, AssignRank::Genus);
    }

    #[test]
    fn genus_fallback_rejected_when_fraction_below_060() {
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // CIGAR with fraction = 50/100 = 0.5 (<0.60) but identity 1.0 and
        // mapped_len above species threshold. Note this actually qualifies
        // for species (identity >= cfg.identity=0.95, mapped_len >= 50) *except*
        // fraction < cfg.fraction => species disqualified. Then fallback:
        // fraction < 0.60 => also disqualified.
        let alns = vec![aln_entry("asmA|1|101|accA|1000", "50=50S", 50, (1, 50))];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns);
        assert!(
            assignment.is_none(),
            "expected no assignment when fraction < 0.60"
        );
    }

    #[test]
    fn genus_fallback_identity_cutoff_is_cfg_genus_identity() {
        // identity of a NEAR_MISS_CIGAR = 92/100 = 0.92.
        let asm = mk_asm_map();
        let s2g = mk_genus_map();

        // With the default genus_identity=0.80, a 0.92 near-miss qualifies.
        let cfg_default = make_cfg("Illumina", None);
        assert!((cfg_default.genus_identity - 0.80).abs() < 1e-12);
        let alns = vec![aln_entry(
            "asmA|1|101|accA|1000",
            NEAR_MISS_CIGAR,
            100,
            (1, 100),
        )];
        let assignment = run_single(&cfg_default, &asm, &s2g, true, "r1", alns.clone())
            .expect("default genus_identity should accept 0.92 near-miss");
        assert_eq!(assignment.rank, AssignRank::Genus);

        // Raising genus_identity above the alignment identity rejects it.
        let mut cfg_strict = make_cfg("Illumina", None);
        cfg_strict.genus_identity = 0.93;
        let assignment = run_single(&cfg_strict, &asm, &s2g, true, "r1", alns.clone());
        assert!(
            assignment.is_none(),
            "identity 0.92 must be rejected when genus_identity=0.93"
        );

        // A 0.82 alignment: ACCEPTED at the new default 0.80 cutoff,
        // REJECTED if the cutoff is raised to 0.85. This exercises the
        // boundary around the old (0.85) and new (0.80) defaults.
        let boundary = vec![aln_entry(
            "asmA|1|101|accA|1000",
            "82=18X", // identity = 0.82, fraction = 1.0
            100,
            (1, 100),
        )];
        let assignment = run_single(&cfg_default, &asm, &s2g, true, "r1", boundary.clone())
            .expect("default genus_identity=0.80 should accept 0.82 alignment");
        assert_eq!(assignment.rank, AssignRank::Genus);

        let mut cfg_old = make_cfg("Illumina", None);
        cfg_old.genus_identity = 0.85;
        let assignment = run_single(&cfg_old, &asm, &s2g, true, "r1", boundary);
        assert!(
            assignment.is_none(),
            "0.82 identity must be rejected when genus_identity is raised to 0.85"
        );

        // Truly low-identity alignment (0.78) must still be rejected at
        // the new default 0.80 cutoff.
        let cfg_default = make_cfg("Illumina", None);
        let below_default = vec![aln_entry(
            "asmA|1|101|accA|1000",
            "78=22X", // identity = 0.78, fraction = 1.0
            100,
            (1, 100),
        )];
        let assignment = run_single(&cfg_default, &asm, &s2g, true, "r1", below_default);
        assert!(
            assignment.is_none(),
            "0.78 identity must be rejected at default genus_identity=0.80"
        );
    }

    #[test]
    fn genus_fallback_skips_alignment_without_resolvable_genus() {
        let cfg = make_cfg("Illumina", None);
        let asm = mk_asm_map();
        // species 999 is missing from the s2g map (no genus resolvable).
        let s2g: HashMap<u32, u32> = HashMap::new();
        let alns = vec![aln_entry(
            "asmA|1|999|accA|1000",
            NEAR_MISS_CIGAR,
            100,
            (1, 100),
        )];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns);
        assert!(
            assignment.is_none(),
            "expected no assignment when genus cannot be resolved"
        );
    }

    // ------------------------------------------------------------------
    // AssignRank labels drive classification output
    // ------------------------------------------------------------------

    #[test]
    fn assign_rank_str_labels() {
        assert_eq!(AssignRank::Species.as_str(), "species");
        assert_eq!(AssignRank::Genus.as_str(), "genus");
    }

    /// End-to-end-ish check: a genus-fallback `ReadAssignment` formats its
    /// rank field as `"genus"` in the same template used by
    /// `write_outputs` for unambiguous classification rows.
    #[test]
    fn genus_fallback_classification_row_uses_genus_rank() {
        let assignment = ReadAssignment {
            read: "readX".to_string(),
            rank: AssignRank::Genus,
            taxids: vec![5001],
            depths: vec![0.1],
            mapped_lens: vec![100.0],
        };
        let taxid = assignment.taxids[0];
        let line = format!(
            "{}\t{}\t{}\t{}\t1",
            assignment.read,
            taxid,
            assignment.rank.as_str(),
            taxid
        );
        assert_eq!(line, "readX\t5001\tgenus\t5001\t1");
    }
}

// ============================================================================
// Statistical Functions for Coverage P-value Calculation
// ============================================================================
//
// Three methods are used depending on sample size:
// 1. Exact Stirling numbers (small n, m)
// 2. Binomial approximation (moderate samples)
// 3. Normal approximation (large samples)

/// Clamp a value to the range [0, 1].
fn clamp01(x: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    }
    if x > 1.0 {
        return 1.0;
    }
    x
}

/// Check if a value is a finite number (not NaN or infinity).
fn is_finite_number(x: f64) -> bool {
    x.is_finite()
}

/// Infer genome size (n) from read count and depth.
///
/// n = read_count / depth
///
/// Used when genome size is not directly available but can be
/// estimated from coverage statistics.
fn infer_n_from_depth(read_count: usize, depth: f64) -> Option<usize> {
    if read_count == 0 || depth <= 0.0 || !is_finite_number(depth) {
        return None;
    }
    let n = (read_count as f64 / depth).round() as isize;
    Some(n.max(1) as usize)
}

/// Infer genome size (n) from read count and expected breadth.
///
/// Uses the relationship: exp_breadth = 1 - (1 - 1/n)^read_count
/// Solving for n: n = -read_count / ln(1 - exp_breadth)
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

/// Convert breadth (fraction) to discrete count of covered positions.
///
/// k = round(breadth * n), clamped to [0, n]
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

/// Compute Stirling numbers of the second kind S(m, k) for k in [0, kmax].
///
/// S(m, k) counts the number of ways to partition m elements into k
/// non-empty subsets. Used in the exact p-value calculation.
///
/// Uses the recurrence: S(m, k) = k * S(m-1, k) + S(m-1, k-1)
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

/// Compute log(k!) using the gamma function.
/// ln(k!) = ln(Gamma(k+1))
fn log_factorial(k: usize) -> f64 {
    ln_gamma(k as f64 + 1.0)
}

/// Compute log of binomial coefficient C(n, k) = n! / (k! * (n-k)!)
/// Uses log-space arithmetic to avoid overflow.
fn log_combination(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

/// Compute log(sum(exp(values))) in a numerically stable way.
///
/// Uses the log-sum-exp trick: log(sum(exp(x))) = m + log(sum(exp(x - m)))
/// where m = max(x). This avoids overflow/underflow.
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
    /// Test if observed is significantly lower than expected
    Lower,
    /// Test if observed is significantly higher than expected
    Upper,
    /// Test if observed differs significantly in either direction
    TwoSided,
}

/// Compute exact p-value using Stirling numbers.
///
/// For small sample sizes (n ≤ 3000, m ≤ 100), computes the exact
/// probability using Stirling numbers of the second kind.
///
/// The distribution is based on the "balls into bins" model where
/// m reads are randomly placed into n genome positions.
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
            // P(X >= k_obs) = 1 - P(X <= k_obs-1)
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

/// Compute p-value using binomial approximation.
///
/// Uses the regularized incomplete beta function to compute the
/// CDF of a binomial distribution with parameters n and p.
///
/// More efficient than exact calculation for moderate sample sizes.
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

/// Compute p-value using normal approximation to binomial.
///
/// For large samples where np(1-p) ≥ 20 and 0.05 ≤ p ≤ 0.95,
/// the binomial distribution is well-approximated by a normal.
///
/// Uses continuity correction (+/- 0.5) for better accuracy.
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
/// Selects the appropriate method based on sample size:
/// 1. Exact Stirling (n ≤ 3000, m ≤ 100) - most accurate
/// 2. Normal approximation (eff ≥ 20, 0.05 ≤ p ≤ 0.95) - fast
/// 3. Binomial approximation - fallback
///
/// # Arguments
/// * `read_count` - Number of reads mapping to this taxon
/// * `depth` - Coverage depth (optional, for inferring genome size)
/// * `breadth` - Observed breadth of coverage
/// * `exp_breadth` - Expected breadth under null hypothesis
/// * `side` - Which tail to test
///
/// # Returns
///
/// P-value in range [0, 1], or None if calculation fails
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

// ============================================================================
// Utility Functions
// ============================================================================

/// Load pathogen-host mapping from TSV file.
///
/// Expected format (tab-separated, with header):
/// pathogen_taxid  host_taxids  host_names  diseases
///
/// Multiple hosts/diseases are semicolon-separated within fields.
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

/// Round a float to 5 decimal places.
fn round5(value: f64) -> f64 {
    (value * 100000.0).round() / 100000.0
}

/// Element-wise accumulation: target[i] += values[i].
/// Extends target vector if values is longer.
fn accumulate(target: &mut Vec<f64>, values: &[f64]) {
    if target.len() < values.len() {
        target.resize(values.len(), 0.0);
    }
    for (t, v) in target.iter_mut().zip(values.iter()) {
        *t += *v;
    }
}

/// Normalize map values to sum to 1.0.
/// Used for converting counts to fractions.
fn normalize_map(map: &mut HashMap<u32, f64>) {
    let sum: f64 = map.values().copied().sum::<f64>();
    if sum <= f64::EPSILON {
        return;
    }
    for value in map.values_mut() {
        *value = round5(*value / sum);
    }
}

/// Format a p-value for output.
///
/// - Very small values (<1e-4) use scientific notation
/// - Other values use 5 decimal places
/// - None or non-finite values become "NA"
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

/// Extract subsampling scaling factor from reference name parts.
///
/// Reference names may include genome_size and sampled_size fields:
/// `asm|taxid|species_taxid|acc|genome_size|sampled_size`
///
/// If sampled_size < genome_size, stores the scaling factor
/// (genome_size / sampled_size) for abundance correction.
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
/// Does nothing if writer is None.
fn write_optional_line(writer: &mut Option<BufWriter<File>>, line: String) -> Result<()> {
    if let Some(inner) = writer.as_mut() {
        inner.write_all(line.as_bytes())?;
        inner.write_all(b"\n")?;
    }
    Ok(())
}

/// Aggregated data from SAM file loading.
///
/// Contains all information needed for the profiling phase.
struct SamData {
    /// Per-read alignment bundles
    read_alignments: Vec<ReadBundle>,

    /// Total assembly lengths: asm_id -> total_length. Shared with the
    /// worker context via `Arc` so no clone is needed when threading
    /// it through the profile / assign_taxonomy / write_outputs chain.
    asm_lengths: Arc<HashMap<String, f64>>,

    /// Minimum breadth threshold for reporting taxa
    min_breadth: f64,

    /// Minimum chunk breadth threshold (auto-estimated or user-specified)
    min_chunk_breadth: f64,

    /// Scaling factors for subsampled genomes: taxid -> scale
    subsample_scaling: HashMap<u32, f64>,

    /// Total distinct reads observed in the SAM/BAM input.
    total_reads: usize,

    /// Distinct reads that produced at least one retained alignment.
    aligned_reads: usize,

    /// Run-level mapping rate = aligned_reads / total_reads.
    /// Used to gate the genus-level fallback mode.
    mapping_rate: f64,

    /// Whether genus-level fallback assignment is active for this run.
    /// Automatic activation requires low-map Illumina non-pathogen data on a
    /// non-subsampled DB; explicit `--genus-fallback` activates it directly.
    fallback_mode: bool,

    /// Sparse map from a raw reference-level taxid (third `|` field of
    /// a reference name) to the species taxid it should be normalized
    /// to. Only populated for raw taxids that are **below** species
    /// rank (subspecies, strain, "no rank" under species, ...).
    /// Missing entries mean "keep the raw taxid as-is". This merges
    /// multiple sub-species belonging to the same species into a
    /// single species-level assignment.
    ///
    /// After the worker-side parse refactor this map is already
    /// applied INSIDE the worker (see [`parse_alignment_line`]) so
    /// every `AlignmentEntry` arrives pre-normalized. The field is
    /// kept on `SamData` for diagnostic inspection / future reuse.
    #[allow(dead_code)]
    species_normalize: Arc<HashMap<u32, u32>>,

    /// Mapping from species taxid to its enclosing genus taxid, derived from
    /// the taxonomy lineage (not parsed from reference headers). Keyed by
    /// the **post-normalization** species taxid so that sub-species
    /// aliases all map through the same entry.
    /// Species without a resolvable genus are simply absent from this map.
    species_to_genus: Arc<HashMap<u32, u32>>,

    /// Alignments that were discarded during SAM/BAM parsing because
    /// their reference name did not carry a parsable species taxid /
    /// assembly length. Seeds [`FilterStats::skipped_reference`] so
    /// the diagnostic counter preserves its historical semantics even
    /// now that the check happens in the worker instead of
    /// `process_batch`.
    skipped_reference_at_load: usize,
}

/// Convert a BAM record to SAM format string.
///
/// Reconstructs the SAM line format from BAM binary representation.
/// Preserves CIGAR operations including '=' and 'X' that are needed
/// for identity calculation.
fn bam_record_to_sam_line(record: &Record, header: &bam::HeaderView) -> Result<String> {
    let qname =
        str::from_utf8(record.qname()).with_context(|| "invalid read name in BAM record")?;
    let flag = record.flags();

    let rname = if record.tid() >= 0 {
        let name_bytes = header.tid2name(record.tid() as u32);
        str::from_utf8(name_bytes)
            .with_context(|| format!("invalid reference name for tid {}", record.tid()))?
            .to_string()
    } else {
        "*".to_string()
    };

    let pos = if record.pos() >= 0 {
        (record.pos() + 1) as i64
    } else {
        0
    };
    let mapq = record.mapq();
    let cigar = if record.cigar().len() == 0 {
        "*".to_string()
    } else {
        // rust-htslib preserves '=' and 'X' operators when stringifying BAM CIGARs,
        // ensuring we retain the match/mismatch split expected by the profiler.
        record.cigar().to_string()
    };

    let rnext = if record.mtid() < 0 {
        "*".to_string()
    } else if record.mtid() == record.tid() {
        "=".to_string()
    } else {
        let mate_bytes = header.tid2name(record.mtid() as u32);
        str::from_utf8(mate_bytes)
            .with_context(|| format!("invalid mate reference name for tid {}", record.mtid()))?
            .to_string()
    };

    let pnext = if record.mpos() >= 0 {
        (record.mpos() + 1) as i64
    } else {
        0
    };
    let tlen = record.insert_size();

    let seq = String::from_utf8(record.seq().as_bytes())
        .with_context(|| format!("invalid sequence for read {}", qname))?;
    let qual = if record.qual().is_empty() {
        "*".to_string()
    } else {
        record
            .qual()
            .iter()
            .map(|&q| (q + 33) as char)
            .collect::<String>()
    };

    Ok(format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
    ))
}

/// Report that CRAM format is unsupported and return error.
///
/// CRAM files normalize CIGAR strings, collapsing '=' and 'X' operations
/// into 'M', which prevents accurate identity calculation.
fn unsupported_cram(path: &Path) -> Result<AlignmentFormat> {
    log::warn!(
        target: "PROFILE",
        "CRAM input {} is unsupported because CRAM CIGAR records collapse '=' and 'X' operations, preventing faithful coverage accounting.",
        path.display()
    );
    bail!("CRAM inputs cannot preserve match/mismatch operations");
}

/// Detect alignment file format from filename and magic bytes.
///
/// Checks:
/// 1. File extension (.sam, .bam, .cram, .sam.gz, .sam.xz)
/// 2. Magic bytes at file start (BAM\x01, CRAM)
///
/// CRAM files are rejected as they don't preserve match/mismatch info.
fn detect_alignment_format(path: &Path) -> Result<AlignmentFormat> {
    if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".sam") || lower.ends_with(".sam.gz") || lower.ends_with(".sam.xz") {
            return Ok(AlignmentFormat::Sam);
        }
        if lower.ends_with(".bam") {
            return Ok(AlignmentFormat::Bam);
        }
        if lower.ends_with(".cram") {
            return unsupported_cram(path);
        }
    }

    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext.to_ascii_lowercase().as_str() {
            "sam" => return Ok(AlignmentFormat::Sam),
            "bam" => return Ok(AlignmentFormat::Bam),
            "cram" => {
                return unsupported_cram(path);
            }
            _ => {}
        }
    }

    let mut file = File::open(path)
        .with_context(|| format!("failed to open {} for format detection", path.display()))?;
    let mut magic = [0u8; 4];
    let read = file.read(&mut magic)?;
    if read == 4 {
        if &magic == b"BAM\x01" {
            return Ok(AlignmentFormat::Bam);
        }
        if &magic == b"CRAM" {
            return unsupported_cram(path);
        }
    }

    Ok(AlignmentFormat::Sam)
}

/// Open a text-based alignment file, handling compression.
///
/// Detects and handles:
/// - Plain text SAM
/// - Gzip-compressed SAM (.sam.gz, magic bytes 1f 8b)
/// - XZ-compressed SAM (.sam.xz, magic bytes fd 37 7a 58 5a 00)
fn open_text_alignment_reader(path: &Path) -> Result<Box<dyn Read>> {
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

//! Community Profiler - terminal-rank taxonomic abundance estimation
//!
//! This module implements the default profiling mode. It reports each
//! reference at the lowest taxonomy rank that reference's metadata
//! confidently supports, rather than rolling every assignment up to species.
//!
//! # Algorithm Overview
//!
//! 1. **Load mappings**: Read Metam coordinate mappings, either streamed from
//!    the mapper or replayed from a saved mapping file. SAM/BAM input remains
//!    supported through `--reuse-sam` for compatibility, as text (optionally
//!    gzip/xz compressed) or binary.
//!
//! 2. **Reference Metadata**: Extract assembly information from reference
//!    names, which carry
//!    `genome_or_assembly_accession|contig_or_sequence_id|taxonomy_field|genome_size`
//!    with the trailing size optional. MANUAL.md describes the format.
//!
//! 3. **Coverage Calculation**: Compute coverage metrics:
//!    - Breadth: Fraction of genome covered by at least one read
//!    - Fixed and flex chunk: Coverage of genome chunks by at least one read
//!
//! 4. **Taxonomy Assignment**: Assign reads to taxa based on alignment
//!    quality, filtering by identity, mapped length, and fraction thresholds.
//!
//! 5. **Ambiguity Resolution**: Multi-mapping reads are resolved using an
//!    Expectation-Maximization (EM) algorithm that iteratively redistributes
//!    reads based on relative taxon abundances.
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
//! - `TaxonCoverage`: Coverage metrics aggregated at the terminal taxon
//! - `ReadAssignment`: Final taxonomic assignment for a read

use std::cmp::Ordering;
use std::collections::{hash_map::Entry, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::rc::Rc;
use std::str;
use std::sync::Arc;
use std::thread;

use anyhow::{anyhow, bail, Context, Result};
use crossbeam_channel::{bounded, unbounded};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{self, Read as BamRead};
use xz2::read::XzDecoder;

use crate::taxonomy::Taxonomy;
use csv::WriterBuilder;

use crate::db::DatabaseMetadata;
use crate::mapper::tsv::{
    for_each_metam_mapping_record, for_each_structured_mapping_record, load_shared_metam_profile,
    open_metam_mapping_reader, total_input_units, MetamMappingReceiver, MetamMappingRecord,
    MetamParseStats, MetamProfile,
};
use crate::progress;

use super::alignment::cigar_to_identity;
use super::evidence::{
    calculate_coverage, evaluate_ranked_witness, guarded_map_winner, resolved_candidate_index,
    standard_ratio_lower_bounds, witness_is_better, CandidateScratch, EvidenceAlignment,
    EvidenceCoverage, EvidenceMode, EvidencePolicy, EvidenceThresholds, ExactTiePlacement,
    ExactTies, MappingRole,
    PresenceAssembly, PresenceInterval, PresenceKey, TaxonEvidence, WitnessDecision,
    STANDARD_RATIO_UPPER,
};
#[cfg(test)]
use super::evidence::{PremapChunkGuard, RescueKinds};
use super::output::{
    ClassificationSource, ProfileAnnotations, ProfileRow, ResolvedClassifications,
};
use super::pathogen::{load_pathogen_table, PathogenEntry};
use super::statistics::{coverage_deviation_pvalue, round5};
use super::{
    coverage_passes, equivalent_read_count_passes_minimum, log_reference_taxonomy_adjustments,
    join_alignment_worker, log_resolved_profiler_parameters, parse_reference_taxonomy,
    parse_sam_sequence_record, push_string_keyed_value, record_genome_sampling_scale,
    record_reference_taxonomy_adjustment, require_alignment_reference_dictionary,
    sampling_scale_from_name_parts,
    sampling_scale_from_sizes, scale_automatic_chunk_breadth, strip_pair_suffix, ProfileMode,
    ProfileRank, ProfilerConfig, ReadSupport, ReferenceTaxonomyAdjustments,
    ReferenceTaxonomyResolver, ResolvedProfilerParameters, StableReadMap,
};

/// Fixed coordinate radius used by the native single-end, paired-end, and
/// long-read evidence policies.
const NATIVE_FIXED_FLANK_LEN: usize = 1_000;

/// Fixed coordinate radius used for assembly input.
const ASSEMBLY_FIXED_FLANK_LEN: usize = 2_000;

#[inline]
fn native_fixed_flank_len(cfg: &ProfilerConfig) -> usize {
    if cfg.sequencer.is_illumina() || cfg.sequencer.is_long_read() {
        NATIVE_FIXED_FLANK_LEN
    } else {
        ASSEMBLY_FIXED_FLANK_LEN
    }
}

fn record_terminal_rank(
    ranks: &mut HashMap<u32, ProfileRank>,
    taxid: u32,
    rank: ProfileRank,
    reference: &str,
) -> Result<()> {
    match ranks.entry(taxid) {
        Entry::Vacant(entry) => {
            entry.insert(rank);
            Ok(())
        }
        Entry::Occupied(entry) if *entry.get() == rank => Ok(()),
        Entry::Occupied(entry) => bail!(
            "conflicting terminal ranks for taxid {}: '{}' and '{}' (reference '{}')",
            taxid,
            entry.get().as_str(),
            rank.as_str(),
            reference,
        ),
    }
}

/// Detected alignment file format.
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
    /// Reusable line buffer, kept on the struct so one allocation is
    /// amortised across every line rather than one per `next()`.
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
    header: Rc<bam::HeaderView>,
}

impl BamRecordIter {
    /// Create a new BAM record iterator.
    ///
    /// # Arguments
    /// * `reader` - Opened BAM file reader
    /// * `header` - BAM header view (shared for memory efficiency)
    fn new(reader: bam::Reader, header: Rc<bam::HeaderView>) -> Self {
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
/// Derived fields are computed once where the alignment enters the profiler --
/// [`CommunityProfiler::entry_from_metam_mapping`] for coordinate mappings,
/// [`parse_alignment_line`] for SAM/BAM -- and stored here ready to use. They
/// are invariant for a given alignment, so the taxonomy-assignment loop is
/// purely arithmetic and touches no `HashMap`, `String::split` or CIGAR
/// scanner per alignment.
type AlignmentEntry = EvidenceAlignment;

/// Narrow a genomic coordinate/length to `u32` for the packed
/// [`AlignmentEntry`] fields. Coordinates are bounded by `ref_len + 1`, and the
/// index rejects references `>= 2^32` bases, so the `.min` only guards the
/// unreachable 4 Gbp-contig edge — every real value is preserved exactly.
#[inline]
fn clamp_u32(value: usize) -> u32 {
    value.min(u32::MAX as usize) as u32
}

/// Widen a packed `(u32, u32)` coordinate pair back to `(usize, usize)` when
/// feeding the (much smaller) per-taxon interval maps and coverage merge.
#[inline]
fn widen_pair(pair: (u32, u32)) -> (usize, usize) {
    (pair.0 as usize, pair.1 as usize)
}

/// Bundle of all alignments for a single read.
///
/// A read may align to multiple reference sequences (multi-mapping).
/// This structure groups all alignments for downstream processing.
#[derive(Clone, Debug)]
struct ReadBundle {
    /// Read identifier (pair suffix stripped for paired-end)
    read_name: Arc<str>,

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

    /// Per-assembly full-index-equivalent support multiplier.
    asm_sampling_scales: Arc<HashMap<String, f64>>,

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

    /// `true` when the alignment was discarded because the reference name did
    /// not parse into the expected `asm|_|species|acc|...` layout, or because
    /// the assembly total length was unknown or zero.
    skipped_reference: bool,
}

/// Coverage metrics aggregated at a terminal taxonomy endpoint.
#[derive(Clone, Debug, Default)]
struct TaxonCoverage {
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
#[derive(Clone, Debug, Hash, Ord, PartialEq, Eq, PartialOrd)]
struct AccKey {
    /// Species taxonomy ID
    taxid: u32,

    /// Assembly identifier (first field of reference name)
    asm: Arc<str>,

    /// Accession identifier (fourth field of reference name)
    acc: Arc<str>,
}

/// Reference metadata pre-parsed for the Metam coordinate path.
///
/// Metam emits the reference name on every mapping row. Splitting that
/// reference name and re-looking-up assembly lengths in the hot loop creates
/// avoidable allocator and hash-table work. This cache moves that work to DB
/// load time, where each reference is parsed once.
#[derive(Clone, Debug)]
struct MetamReferenceEntry {
    ref_len: usize,
    asm: Arc<str>,
    acc: Arc<str>,
    profile_taxid: u32,
    profile_rank: ProfileRank,
    reference_id: u32,
    assembly_id: u32,
    sampling_scale: f32,
    asm_len: f64,
    flex_len: usize,
}

/// Statistics about alignment filtering.
///
/// Tracks how many alignments were evaluated, passed filters,
/// and why alignments were rejected.
#[derive(Clone, Copy, Debug, Default)]
struct FilterStats {
    evaluated: usize,

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
    /// Merge statistics from a parallel worker.
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

    /// Terminal-rank, assembly-coherent coordinate evidence for native short-
    /// and long-read profiling. Different assemblies are always independent
    /// witnesses and are never unioned.
    presence: HashMap<PresenceKey, PresenceAssembly>,

    /// Candidate, resolved-best, and direct-clear fragment counts used by the
    /// evidence-policy rescue tiers.
    taxon_evidence: HashMap<u32, TaxonEvidence>,

    /// Compact coordinate payload for exact top-score short-read ties.
    /// Populated only when guarded MAP is active.
    exact_ties: ExactTies,

    /// Terminal taxa observed as the only valid candidate for at least one
    /// fragment. This is the multiple-testing basis for the evidence policy.
    uniquely_supported_taxa: HashSet<u32>,

    /// Fragments for which a mapper boundary may hide an unrepresented
    /// profile endpoint. Emitted candidates remain in abundance EM, but
    /// cannot establish direct coordinate evidence from an incomplete set.
    boundary_ambiguous_fragments: usize,

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
    Family,
    Order,
    Class,
    Phylum,
    Kingdom,
    /// Candidate endpoints span more than one terminal rank.
    Mixed,
}

impl AssignRank {
    fn from_profile_rank(rank: ProfileRank) -> Self {
        match rank {
            ProfileRank::Species => Self::Species,
            ProfileRank::Genus => Self::Genus,
            ProfileRank::Family => Self::Family,
            ProfileRank::Order => Self::Order,
            ProfileRank::Class => Self::Class,
            ProfileRank::Phylum => Self::Phylum,
            ProfileRank::Kingdom => Self::Kingdom,
        }
    }
}

/// Taxonomic assignment for a single read.
///
/// A read may be assigned to one taxon (unambiguous) or multiple taxa
/// (ambiguous multi-mapping). Ambiguous assignments are resolved in
/// the EM phase.
///
/// Native candidates can span terminal ranks; those assignments use
/// [`AssignRank::Mixed`]. Legacy and fallback assignments remain homogeneous.
#[derive(Clone, Debug)]
struct ReadAssignment {
    read: Arc<str>,

    /// Homogeneous terminal rank, genus fallback, or `Mixed`.
    rank: AssignRank,

    /// `true` when this assignment came from the native terminal-evidence
    /// resolver and its candidates must pass the assembly-coherent presence
    /// policy. Legacy SAM/BAM and genus-fallback assignments use their
    /// existing coverage path even when their taxid is also a direct database
    /// endpoint.
    requires_terminal_evidence: bool,

    /// Candidate terminal taxids (may be >1 for multi-mappers).
    taxids: Vec<u32>,

    /// Depth contribution for each candidate (mapped_len / genome_size)
    depths: Vec<f64>,

    /// Mapped lengths for each candidate
    mapped_lens: Vec<f64>,

    /// Genome-specific multiplier for each representative candidate mapping.
    sampling_scales: Vec<f64>,
}

#[derive(Clone, Copy, Debug, Default)]
struct TaxonTotals {
    support: ReadSupport,
    depth: f64,
    mapped_len: f64,
}

#[derive(Clone, Debug, Default)]
struct AmbiguousTotals {
    observed_count: f64,
    depths: Vec<f64>,
    mapped_lens: Vec<f64>,
    /// Sum of the candidate genome scales over all reads in this group.
    equivalent_support_sums: Vec<f64>,
}

/// Optional, diagnostic abundance estimate built from every quality-qualified
/// read assignment before taxon-level coverage/evidence selection.
///
/// This is intentionally separate from the production EM state: enabling
/// `--save-raw-profile` must not add candidates to `profile.txt` or
/// `classify.txt`, or otherwise alter the established profiling policy.
struct RawProfileEstimate {
    totals: HashMap<u32, TaxonTotals>,
    coverage: HashMap<u32, TaxonCoverage>,
}

/// Community-level taxonomy profiler.
///
/// Implements terminal-rank profiling by grouping native mappings at each
/// reference-supported endpoint and resolving ambiguous assignments with EM.
///
/// # Pipeline
///
/// 1. Load coordinate mappings and extract reference metadata
/// 2. Compute coverage intervals for each reference sequence
/// 3. Assign reads to taxa based on alignment quality
/// 4. Resolve multi-mapping reads with EM algorithm
/// 5. Filter taxa based on coverage statistics
/// 6. Output abundance profile and classification files
pub struct CommunityProfiler {
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
        log::info!(target: "PROFILER", "Initializing the profiler!");

        // Load and parse alignment data. Normal Metax runs consume
        // Metam's coordinate TSV directly; SAM/BAM remains supported for
        // --reuse-sam compatibility.
        let sam_data = if let Some(mappings) = self.cfg.mapper_mappings.as_deref() {
            let db = self
                .cfg
                .db_metadata
                .as_deref()
                .ok_or_else(|| anyhow!("missing reference metadata for native mappings"))?;
            self.load_metam(Path::new(mappings), db)?
        } else {
            let sam_path = Path::new(&self.cfg.sam);
            if !sam_path.exists() {
                bail!("SAM file {} does not exist", sam_path.display());
            }
            self.load_sam(sam_path)?
        };

        // Run profiling pipeline
        self.profile(sam_data)
    }

    pub fn run_with_metam_mappings(
        &self,
        receiver: MetamMappingReceiver,
        source: &str,
    ) -> Result<()> {
        log::info!(target: "PROFILER", "Initializing the profiler!");
        let db = self
            .cfg
            .db_metadata
            .as_deref()
            .ok_or_else(|| anyhow!("missing Metax database metadata for streamed mappings"))?;
        let sam_data = self.load_metam_mappings(receiver, source, db)?;
        self.profile(sam_data)
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

    fn load_metam(&self, mappings_path: &Path, db: &DatabaseMetadata) -> Result<SamData> {
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
    ) -> Result<SamData> {
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
    ) -> Result<SamData> {
        self.load_metam_records(source, db, |on_mapping| {
            for_each_structured_mapping_record(receiver, source, on_mapping)
        })
    }

    fn load_metam_records<F>(&self, source: &str, db: &DatabaseMetadata, feed: F) -> Result<SamData>
    where
        F: FnOnce(
            &mut dyn for<'a> FnMut(MetamMappingRecord<'a>) -> Result<()>,
        ) -> Result<MetamParseStats>,
    {
        log::info!(target: "PROFILER", "Loading Metam coordinate mappings ...");
        let fixed_flank_len = native_fixed_flank_len(&self.cfg);
        let mut asm_lengths_raw: HashMap<String, f64> = HashMap::new();
        let mut genome_sampling_scales: HashMap<String, f64> = HashMap::new();
        let mut profile_taxids: HashSet<u32> = HashSet::new();
        let mut terminal_ranks: HashMap<u32, ProfileRank> = HashMap::new();
        let mut reference_taxonomies = Vec::with_capacity(db.references.len());
        let mut taxonomy_adjustments = ReferenceTaxonomyAdjustments::new();
        let mut taxonomy_resolver = ReferenceTaxonomyResolver::new(&self.taxonomy);
        for reference in &db.references {
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
            let asm = reference.genome.clone();
            *asm_lengths_raw.entry(asm).or_insert(0.0) += reference.length as f64;
            profile_taxids.insert(taxonomy.profile_taxid);
            record_terminal_rank(
                &mut terminal_ranks,
                taxonomy.profile_taxid,
                taxonomy.profile_rank,
                &reference.name,
            )?;
            let scale = sampling_scale_from_sizes(
                reference.full_genome_size,
                reference.effective_indexed_size,
            );
            record_genome_sampling_scale(&mut genome_sampling_scales, &reference.genome, scale)?;
            reference_taxonomies.push(taxonomy);
        }
        log_reference_taxonomy_adjustments(
            &taxonomy_adjustments,
            self.cfg.verbose || self.cfg.very_verbose,
        );
        let asm_flex_lengths: HashMap<String, usize> = asm_lengths_raw
            .iter()
            .map(|(asm, len)| {
                let flex = len.sqrt().round() as usize;
                (asm.clone(), flex.max(1))
            })
            .collect();
        let (species_normalize_raw, species_to_genus_raw) =
            self.build_taxonomy_lookups(&profile_taxids);
        let db_is_subsampled = db.is_fractional;
        let species_normalize = Arc::new(species_normalize_raw);
        let species_to_genus = Arc::new(if db_is_subsampled && !self.cfg.genus_fallback {
            HashMap::new()
        } else {
            species_to_genus_raw
        });
        let mut assembly_ids: HashMap<&str, u32> = HashMap::new();
        for reference in &db.references {
            if !assembly_ids.contains_key(reference.genome.as_str()) {
                let next = u32::try_from(assembly_ids.len())
                    .context("database contains more than 2^32 assemblies")?;
                assembly_ids.insert(reference.genome.as_str(), next);
            }
        }
        let mut metam_refs: HashMap<String, MetamReferenceEntry> =
            HashMap::with_capacity(db.references.len());
        for (reference_index, (reference, taxonomy)) in
            db.references.iter().zip(&reference_taxonomies).enumerate()
        {
            let asm = reference.genome.as_str();
            let raw_profile_taxid = taxonomy.profile_taxid;
            let acc = reference.contig.as_str();
            let asm_len = match asm_lengths_raw.get(asm).copied() {
                Some(value) if value > 0.0 => value,
                _ => continue,
            };
            let profile_taxid = species_normalize
                .get(&raw_profile_taxid)
                .copied()
                .unwrap_or(raw_profile_taxid);
            let flex_len = asm_flex_lengths
                .get(asm)
                .copied()
                .unwrap_or(fixed_flank_len);
            let sampling_scale = genome_sampling_scales.get(asm).copied().unwrap_or(1.0) as f32;
            let reference_id = u32::try_from(reference_index)
                .context("database contains more than 2^32 references")?;
            let assembly_id = assembly_ids[asm];
            if metam_refs.contains_key(&reference.name) {
                bail!(
                    "database metadata contains duplicate reference name '{}'",
                    reference.name
                );
            }
            metam_refs.insert(
                reference.name.clone(),
                MetamReferenceEntry {
                    ref_len: reference.length,
                    asm: Arc::<str>::from(asm),
                    acc: Arc::<str>::from(acc),
                    profile_taxid,
                    profile_rank: taxonomy.profile_rank,
                    reference_id,
                    assembly_id,
                    sampling_scale,
                    asm_len,
                    flex_len,
                },
            );
        }

        let mut read_map: StableReadMap<AlignmentEntry> = StableReadMap::default();
        let mut num_skipped_reference = 0usize;
        let mut num_alignments = 0usize;
        let load_progress_options = self.cfg.mapping_load_progress_options();
        let load_progress = progress::counter(load_progress_options, "Loading mapping records");
        let parse_result = feed(&mut |mapping| {
            num_alignments += 1;
            progress::tick_count(
                &load_progress,
                "Loading mapping records",
                num_alignments as u64,
                load_progress_options.interval(),
            );
            if self.cfg.is_paired && (!mapping.mate_ref_status_same_ref || mapping.pair_id == 0) {
                return Ok(());
            }
            match self.entry_from_metam_mapping(mapping, &metam_refs, fixed_flank_len) {
                Some((read_id, entry)) => {
                    push_string_keyed_value(&mut read_map, read_id, entry);
                }
                None => {
                    num_skipped_reference += 1;
                }
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

        let aligned_reads = read_map.len();
        let mapped_unit = if self.cfg.is_paired {
            "fragments"
        } else {
            "reads"
        };
        let profile = self.load_metam_profile(parse_stats.profile)?;
        let total_reads = total_input_units(&profile, self.cfg.is_paired).max(aligned_reads);
        let read_alignments = read_map
            .into_iter()
            .map(|(read_name, alignments)| ReadBundle {
                read_name: Arc::<str>::from(read_name),
                alignments,
            })
            .collect::<Vec<_>>();
        let evidence_policy = EvidencePolicy::from_config(&self.cfg);

        let mapping_rate = if total_reads == 0 {
            0.0
        } else {
            aligned_reads as f64 / total_reads as f64
        };
        let gate_active = fallback_gate_active(&self.cfg, mapping_rate);
        let fallback_mode = genus_fallback_mode_active(&self.cfg, gate_active, db_is_subsampled);
        let chunk_plan = plan_chunk_breadth(
            &self.cfg,
            mapping_rate,
            total_reads,
            aligned_reads,
            db_is_subsampled,
            gate_active,
        );
        let min_breadth = self.cfg.breadth.unwrap_or(0.0);
        let skip_auto_chunk_filter =
            self.cfg.min_reads.is_some() && self.cfg.chunk_breadth.is_none();
        let chunk_thresholds =
            resolve_community_chunk_breadth(&self.cfg, chunk_plan.estimate_reads);
        let min_chunk_breadth_value = chunk_thresholds.standard;
        let chunk_breadth_source =
            community_chunk_breadth_source(&self.cfg, chunk_plan, chunk_thresholds);

        log::info!(target: "PROFILER", "Total number of reads: {}", total_reads);
        log::info!(target: "PROFILER", "Number of mapping records: {}", num_alignments);
        log::info!(
            target: "PROFILER",
            "Number of mapped {}: {}",
            mapped_unit,
            aligned_reads
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
        } else if self.cfg.chunk_breadth.is_none() {
            if chunk_plan.use_aligned_basis {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from mapped {} ({}).",
                    mapped_unit,
                    aligned_reads
                );
            } else {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from total reads ({}).",
                    total_reads
                );
            }
        }
        log_chunk_breadth_scaling(chunk_thresholds);
        log_resolved_profiler_parameters(
            &self.cfg,
            ResolvedProfilerParameters {
                min_breadth,
                min_chunk_breadth: min_chunk_breadth_value,
                chunk_breadth_source,
                chunk_breadth_basis_reads: chunk_plan.estimate_reads,
                total_reads,
                mapped_reads: aligned_reads,
                mapping_rate,
            },
        );

        Ok(SamData {
            read_alignments,
            asm_lengths: Arc::new(asm_lengths_raw),
            min_breadth,
            min_chunk_breadth: min_chunk_breadth_value,
            safeguard_min_chunk_breadth: chunk_thresholds.safeguard,
            has_scaled_genomes: genome_sampling_scales.values().any(|scale| *scale > 1.0),
            total_reads,
            aligned_reads,
            mapping_rate,
            fallback_mode,
            species_to_genus,
            terminal_ranks: Arc::new(terminal_ranks),
            skipped_reference_at_load: num_skipped_reference,
            evidence_policy,
        })
    }

    fn entry_from_metam_mapping<'a>(
        &self,
        mapping: MetamMappingRecord<'a>,
        metam_refs: &HashMap<String, MetamReferenceEntry>,
        fixed_flank_len: usize,
    ) -> Option<(&'a str, AlignmentEntry)> {
        if self.cfg.is_paired && (!mapping.mate_ref_status_same_ref || mapping.pair_id == 0) {
            return None;
        }
        let reference = metam_refs.get(mapping.ref_id)?;
        if reference.ref_len == 0 || mapping.read_len == 0 || mapping.mapped_bases == 0 {
            return None;
        }
        let start = mapping.start_1based.clamp(1, reference.ref_len.max(1));
        let mut end = mapping.end_1based_exclusive.max(start.saturating_add(1));
        end = end.min(reference.ref_len.saturating_add(1));
        let midpoint = start + ((end.saturating_sub(start)) / 2);
        let fixed_chunk_start = start.saturating_sub(fixed_flank_len).max(1);
        let fixed_chunk_end = end
            .saturating_add(fixed_flank_len)
            .min(reference.ref_len.saturating_add(1));
        let flex_chunk_start = midpoint.saturating_sub(reference.flex_len).max(1);
        let flex_chunk_end = midpoint
            .saturating_add(reference.flex_len)
            .min(reference.ref_len.saturating_add(1));
        let fraction = (mapping.mapped_bases as f64 / mapping.read_len as f64).min(1.0);
        let read_name: &'a str = if self.cfg.is_paired {
            strip_pair_suffix(mapping.read_id)
        } else {
            mapping.read_id
        };
        Some((
            read_name,
            AlignmentEntry {
                asm: Arc::clone(&reference.asm),
                acc: Arc::clone(&reference.acc),
                profile_taxid: reference.profile_taxid,
                profile_rank: reference.profile_rank,
                reference_id: reference.reference_id,
                assembly_id: reference.assembly_id,
                role: if self.cfg.is_paired {
                    match mapping.mate_number {
                        1 => MappingRole::Mate1,
                        2 => MappingRole::Mate2,
                        _ => return None,
                    }
                } else {
                    MappingRole::Unpaired
                },
                pair_id: mapping.pair_id,
                candidate_boundary: mapping.candidate_boundary,
                sampling_scale: reference.sampling_scale,
                asm_len: reference.asm_len,
                identity: mapping.identity as f32,
                fraction: fraction as f32,
                mapped_len: clamp_u32(mapping.mapped_bases),
                span: (clamp_u32(start), clamp_u32(end)),
                fixed_chunk: (clamp_u32(fixed_chunk_start), clamp_u32(fixed_chunk_end)),
                flex_chunk: (clamp_u32(flex_chunk_start), clamp_u32(flex_chunk_end)),
            },
        ))
    }

    /// Load alignments from a SAM or BAM file.
    ///
    /// Automatically detects the file format (SAM/BAM) and handles
    /// compression (gzip/xz for SAM files).
    fn load_sam(&self, sam_path: &Path) -> Result<SamData> {
        log::info!(target: "PROFILER", "Loading alignments ...");
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
        let mut asm_sampling_scales: HashMap<String, f64> = HashMap::new();
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
                    let (sn, len) = parse_sam_sequence_record(trimmed)?;
                    // Inspect ONLY the first @SQ line to decide whether
                    // the run is against a subsampled reference DB.
                    if db_is_subsampled.is_none() {
                        db_is_subsampled = Some(ref_name_is_subsampled(&sn));
                    }
                    if ref_lengths.insert(sn.clone(), len).is_some() {
                        bail!("SAM header contains duplicate reference name '{}'", sn);
                    }

                    // Parse reference name to extract assembly ID
                    let parts: Vec<&str> = sn.split('|').collect();
                    let asm = parts.first().copied().unwrap_or(&sn).to_string();
                    record_genome_sampling_scale(
                        &mut asm_sampling_scales,
                        &asm,
                        sampling_scale_from_name_parts(&parts),
                    )?;

                    // Accumulate total assembly length
                    *asm_lengths.entry(asm).or_insert(0.0) += len as f64;

                    // Collect species taxid (third field) for genus lookup.
                    if let Some(species_taxid) =
                        parts.get(2).and_then(|s| s.parse::<u32>().ok())
                    {
                        species_taxids.insert(species_taxid);
                    }
                }
                continue;
            } else {
                // First alignment line found - save and exit header parsing
                first_alignment = Some(trimmed.to_string());
                break;
            }
        }
        require_alignment_reference_dictionary(sam_path, ref_lengths.len())?;

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
        let has_scaled_genomes = asm_sampling_scales.values().any(|scale| *scale > 1.0);

        // Create shared context for worker threads (now carries the
        // pre-built taxonomy and assembly-length lookups).
        let worker_context = Arc::new(WorkerContext {
            ref_lengths: Arc::new(ref_lengths),
            asm_lengths: Arc::clone(&asm_lengths),
            asm_sampling_scales: Arc::new(asm_sampling_scales),
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
            species_to_genus,
            db_is_subsampled,
            has_scaled_genomes,
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
                    target: "PROFILER",
                    "Unable to enable multi-threaded decoding: {}",
                    err
                );
            }
        }
        let header_view = Rc::new(reader.header().clone());

        let mut ref_lengths: HashMap<String, usize> = HashMap::new();
        let mut asm_lengths: HashMap<String, f64> = HashMap::new();
        let mut asm_sampling_scales: HashMap<String, f64> = HashMap::new();
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
            record_genome_sampling_scale(
                &mut asm_sampling_scales,
                &asm,
                sampling_scale_from_name_parts(&parts),
            )?;
            *asm_lengths.entry(asm).or_insert(0.0) += len as f64;
            if let Some(species_taxid) = parts.get(2).and_then(|s| s.parse::<u32>().ok()) {
                species_taxids.insert(species_taxid);
            }
        }
        require_alignment_reference_dictionary(sam_path, ref_lengths.len())?;

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
        let has_scaled_genomes = asm_sampling_scales.values().any(|scale| *scale > 1.0);

        let worker_context = Arc::new(WorkerContext {
            ref_lengths: Arc::new(ref_lengths),
            asm_lengths: Arc::clone(&asm_lengths),
            asm_sampling_scales: Arc::new(asm_sampling_scales),
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
            species_to_genus,
            db_is_subsampled,
            has_scaled_genomes,
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
    fn consume_alignment_stream<I>(
        &self,
        first_record: Option<String>,
        iter: I,
        worker_context: Arc<WorkerContext>,
        asm_lengths: Arc<HashMap<String, f64>>,
        species_to_genus: Arc<HashMap<u32, u32>>,
        db_is_subsampled: bool,
        has_scaled_genomes: bool,
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
                    if (!reads_seen.is_empty() || alignments_seen > 0 || skipped_reference > 0)
                        && tx
                            .send(WorkerResult {
                                alignments,
                                reads_seen,
                                alignments_seen,
                                skipped_reference,
                            })
                            .is_err()
                    {
                        break;
                    }
                }
            }));
        }
        // Drop our copy of result_tx so channel closes when workers finish
        drop(result_tx);
        // Only worker-held receivers may keep the work channel connected. If
        // every worker panics, the next send now fails instead of blocking on
        // this unused receiver forever.
        drop(line_rx);

        // Main-thread aggregation: one `.entry()` per alignment on
        // `read_map`, one `.insert()` per unmapped/discarded row on
        // `all_reads` (aligned rows carry `read_name` only via
        // `read_map`; see end-of-stream read-count reconciliation).
        let mut read_map: StableReadMap<AlignmentEntry> = StableReadMap::default();
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
            join_alignment_worker(handle, "community SAM parsing worker")?;
        }
        if num_skipped_reference > 0 {
            bail!(
                "SAM contains {} mapped record(s) with a reference absent from @SQ or incompatible reference metadata",
                num_skipped_reference
            );
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
        let mapped_unit = if self.cfg.is_paired {
            "fragments"
        } else {
            "reads"
        };
        let read_alignments = read_map
            .into_iter()
            .map(|(read_name, alignments)| ReadBundle {
                read_name: Arc::<str>::from(read_name),
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
        let skip_auto_chunk_filter =
            self.cfg.min_reads.is_some() && self.cfg.chunk_breadth.is_none();
        let chunk_thresholds =
            resolve_community_chunk_breadth(&self.cfg, chunk_plan.estimate_reads);
        let min_chunk_breadth_value = chunk_thresholds.standard;
        let chunk_breadth_source =
            community_chunk_breadth_source(&self.cfg, chunk_plan, chunk_thresholds);

        if self.cfg.verbose || self.cfg.very_verbose {
            const INDENT: &str = "                                                ";
            log::info!(
                target: "PROFILER",
                "Profiler parameters:\n{}min_breadth        {}\n{}min_chunk_breadth  {}\n{}min_reads          {}\n{}min_oebr           {}\n{}min_coebr          {}\n{}identity           {}\n{}mapped_len   {}\n{}fraction           {}\n{}batch_size         {}\n",
                INDENT,
                min_breadth,
                INDENT,
                min_chunk_breadth_value,
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

            log::info!(target: "PROFILER", "Total number of reads: {}", num_reads);
            log::info!(target: "PROFILER", "Number of mapping records: {}", num_alignments);
            log::info!(
                target: "PROFILER",
                "Number of mapped {}: {}",
                mapped_unit,
                aligned_reads
            );
            if self.cfg.chunk_breadth.is_some() {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth supplied by user; read-count basis not used for estimation."
                );
            } else if skip_auto_chunk_filter {
                log::info!(
                    target: "PROFILER",
                    "Automatic minimum chunk breadth filter disabled because --min-reads is set."
                );
            } else if self.cfg.lowbiomass {
                log::info!(
                    target: "PROFILER",
                    "Automatic minimum chunk breadth filter disabled by --low-biomass."
                );
            } else {
                // Log the unscaled basis selected for the chunk-breadth curve.
                // Any optional 1.5x boost is logged separately below.
                if !chunk_plan.use_aligned_basis {
                    log::info!(
                        target: "PROFILER",
                        "Minimum chunk breadth estimated from total reads ({}).",
                        num_reads
                    );
                    if chunk_plan.force_total_for_subsampled_fallback {
                        log::info!(
                            target: "PROFILER",
                            "--genus-fallback with a subsampled reference DB: chunk-breadth estimate uses total read count."
                        );
                    }
                } else if chunk_plan.by_aligned_auto {
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
                if chunk_plan.scaled {
                    log::info!(
                        target: "PROFILER",
                        "Low-map genus fallback active and --min-chunk-breadth not set: mapped-read chunk basis scaled by {:.2}x ({} -> {}) before applying the chunk-breadth equation.",
                        FALLBACK_CHUNK_BASIS_SCALE,
                        chunk_plan.basis_reads,
                        chunk_plan.estimate_reads,
                    );
                }
            }
        }
        log_chunk_breadth_scaling(chunk_thresholds);

        if self.cfg.verbose || self.cfg.very_verbose {
            log::info!(
                target: "PROFILER",
                "Mapping rate: {:.5} (mapped {} / total {})",
                mapping_rate,
                aligned_reads,
                num_reads
            );
            if fallback_mode {
                if gate_active && !db_is_subsampled {
                    log::info!(
                        target: "PROFILER",
                        "Genus-level fallback mode ENABLED (low-map Illumina, host=None, mapping_rate {:.5} < {:.5}; genus_identity={:.3}).",
                        mapping_rate,
                        self.cfg.low_map_rate_threshold,
                        self.cfg.genus_identity,
                    );
                } else if self.cfg.genus_fallback && db_is_subsampled {
                    log::info!(
                        target: "PROFILER",
                        "Genus-level fallback mode ENABLED (--genus-fallback on subsampled reference DB; genus_identity={:.3}).",
                        self.cfg.genus_identity,
                    );
                } else {
                    log::info!(
                        target: "PROFILER",
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
                    target: "PROFILER",
                    "Automatic low-map genus fallback mode DISABLED: reference database appears to be subsampled (first @SQ name has 6 `|`-separated fields)."
                );
            }
        }
        log_resolved_profiler_parameters(
            &self.cfg,
            ResolvedProfilerParameters {
                min_breadth,
                min_chunk_breadth: min_chunk_breadth_value,
                chunk_breadth_source,
                chunk_breadth_basis_reads: chunk_plan.estimate_reads,
                total_reads: num_reads,
                mapped_reads: aligned_reads,
                mapping_rate,
            },
        );

        Ok(SamData {
            read_alignments,
            asm_lengths,
            min_breadth,
            min_chunk_breadth: min_chunk_breadth_value,
            safeguard_min_chunk_breadth: chunk_thresholds.safeguard,
            has_scaled_genomes,
            total_reads: num_reads,
            aligned_reads,
            mapping_rate,
            fallback_mode,
            species_to_genus,
            terminal_ranks: Arc::new(HashMap::new()),
            skipped_reference_at_load: num_skipped_reference,
            evidence_policy: None,
        })
    }

    /// Run the main profiling algorithm on loaded alignment data.
    ///
    /// Steps:
    /// 1. Assign taxonomy to reads based on alignment quality
    /// 2. Compute coverage metrics for each taxon
    /// 3. Resolve ambiguous multi-mapping reads with EM
    /// 4. Filter taxa based on coverage statistics
    /// 5. Write output files
    fn profile(&self, mut sam_data: SamData) -> Result<()> {
        let SamData {
            read_alignments,
            fallback_mode,
            species_to_genus,
            mapping_rate,
            total_reads,
            aligned_reads,
            skipped_reference_at_load,
            evidence_policy,
            ..
        } = &sam_data;
        let verbose_logging = self.cfg.verbose || self.cfg.very_verbose;
        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "Assigning taxonomy for reads with {} threads ...",
                self.cfg.threads
            );
            log::info!(
                target: "PROFILER",
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
            *evidence_policy,
        )?;
        // The `skipped_reference` diagnostic now accumulates two
        // streams: parser-time rejections (unparsable species taxid /
        // missing assembly length; moved here from `process_batch`)
        // plus any future process_batch-time rejections (currently
        // none).
        assignment.filter_stats.skipped_reference += *skipped_reference_at_load;
        if verbose_logging {
            log::info!(target: "PROFILER", "Assigning taxonomy finished!");
        }
        let candidate_reads = assignment.read_assignments.len();
        if assignment.boundary_ambiguous_fragments > 0 {
            log::info!(
                target: "PROFILER",
                "Candidate-boundary safeguard kept {} {} in abundance ambiguity but excluded them from direct coordinate evidence.",
                assignment.boundary_ambiguous_fragments,
                if self.cfg.is_paired { "fragments" } else { "reads" },
            );
        }
        let assignment_unit = if self.cfg.is_paired {
            "fragments"
        } else {
            "reads"
        };
        // Assignment owns every value required by EM and presence filtering.
        // Release the dominant all-mappings buffer before coverage reduction
        // and output construction so mapper/profiler peak RSS does not carry
        // both representations at once.
        drop(std::mem::take(&mut sam_data.read_alignments));
        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "{} with candidate taxon assignments before coverage filtering: {}",
                if self.cfg.is_paired { "Fragments" } else { "Reads" },
                candidate_reads
            );
            log::info!(
                target: "PROFILER",
                "Profiling taxonomy abundance ..."
            );
        }
        let final_classified = self.write_outputs(sam_data, assignment)?;
        log::info!(
            target: "PROFILER",
            "Number of {} assigned to retained taxa: {}",
            assignment_unit,
            final_classified
        );
        log::info!(target: "PROFILER", "Taxonomy profiling finished.");
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
        evidence_policy: Option<EvidencePolicy>,
    ) -> Result<AssignmentResult> {
        let batch_size = self.cfg.batch_size.max(1);
        let cfg = Arc::new(self.cfg.clone());

        // Create rayon thread pool for parallel batch processing
        let pool = ThreadPoolBuilder::new()
            .num_threads(self.cfg.threads.max(1))
            .build()?;

        let pb = progress::bar(
            self.cfg.progress,
            read_alignments.len() as u64,
            "Assigning taxonomy",
        );
        let progress = pb.clone();

        // Process read batches in parallel
        let mut results: Vec<AssignmentResult> = pool.install(|| {
            read_alignments
                .par_chunks(batch_size)
                .map_with(progress, |progress_bar, chunk| {
                    let result = if let Some(policy) = evidence_policy {
                        process_evidence_batch(
                            chunk,
                            &cfg,
                            policy,
                            fallback_mode,
                            species_to_genus.as_ref(),
                        )
                    } else {
                        process_batch(chunk, &cfg, fallback_mode, species_to_genus.as_ref())
                    };
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
                .extend(result.read_assignments);

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
            for (key, mut presence) in result.presence {
                let target = final_result.presence.entry(key).or_default();
                if target.asm_len <= 0.0 {
                    target.asm_len = presence.asm_len;
                }
                target.count = target.count.saturating_add(presence.count);
                target.reference_span_bases = target
                    .reference_span_bases
                    .saturating_add(presence.reference_span_bases);
                target.intervals.append(&mut presence.intervals);
            }
            for (taxid, evidence) in result.taxon_evidence {
                let target = final_result.taxon_evidence.entry(taxid).or_default();
                target.candidate_fragments = target
                    .candidate_fragments
                    .saturating_add(evidence.candidate_fragments);
                target.resolved_candidate_fragments = target
                    .resolved_candidate_fragments
                    .saturating_add(evidence.resolved_candidate_fragments);
                target.best_fragments = target
                    .best_fragments
                    .saturating_add(evidence.best_fragments);
                target.clear_best_fragments = target
                    .clear_best_fragments
                    .saturating_add(evidence.clear_best_fragments);
            }
            final_result.exact_ties.append(result.exact_ties);
            final_result
                .uniquely_supported_taxa
                .extend(result.uniquely_supported_taxa);
            final_result.boundary_ambiguous_fragments = final_result
                .boundary_ambiguous_fragments
                .saturating_add(result.boundary_ambiguous_fragments);

            // Merge filter statistics
            final_result.filter_stats.merge(result.filter_stats);
        }
        Ok(final_result)
    }

    /// Write profiling results to output files.
    ///
    /// This is a large function that:
    /// 1. Computes final coverage metrics per terminal taxon
    /// 2. Performs EM algorithm for ambiguous read redistribution
    /// 3. Calculates abundance estimates and statistical p-values
    /// 4. Filters taxa based on coverage criteria
    /// 5. Writes profile.txt and, when requested, classify.txt output files
    ///
    /// # Output Files
    ///
    /// - `<prefix>.profile.txt`: Filtered abundance profile with statistics
    /// - `<prefix>.rprofile.txt`: Post-EM profile before taxon-level evidence
    ///   and reporting filters (if keep_raw=true)
    /// - `<prefix>.classify.txt`: Optional per-read assignments
    ///   (`--save-classification`)
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
    fn write_outputs(&self, sam_data: SamData, mut assignment: AssignmentResult) -> Result<usize> {
        let SamData {
            asm_lengths,
            min_breadth,
            min_chunk_breadth,
            safeguard_min_chunk_breadth,
            has_scaled_genomes,
            evidence_policy,
            terminal_ranks,
            ..
        } = sam_data;
        let verbose_logging = self.cfg.verbose || self.cfg.very_verbose;
        // Majority-LCA rows are evaluated later than native coordinate
        // witnesses. Resolve the same platform mode here so their default
        // ratio bounds come from the standard evidence policy rather than the
        // older output-only genus/species defaults.
        let lca_ratio_mode = EvidencePolicy::from_config(&self.cfg).map(|policy| policy.mode);

        let mut taxon_coverage: HashMap<u32, TaxonCoverage> = HashMap::new();
        let mut coverage_keys = assignment.acc_intervals.keys().collect::<Vec<_>>();
        coverage_keys.sort_unstable();
        for key in coverage_keys {
            let intervals = &assignment.acc_intervals[key];
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
            let asm_len = asm_lengths
                .get(key.asm.as_ref())
                .copied()
                .unwrap_or(1.0)
                .max(1.0);
            let total_span: usize = merged.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let total_fixed: usize = merged_fixed.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let total_flex: usize = merged_flex.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let coverage = taxon_coverage.entry(key.taxid).or_default();
            coverage.breadth += total_span as f64 / asm_len;
            coverage.fixed_chunk += total_fixed as f64 / asm_len;
            coverage.flex_chunk += total_flex as f64 / asm_len;
        }

        let mut evidence_decisions: HashMap<u32, WitnessDecision> = HashMap::new();
        let mut evidence_selected_taxa: HashSet<u32> = HashSet::new();
        if let Some(policy) = evidence_policy {
            let resolved = resolve_coordinate_evidence(
                &self.cfg,
                policy,
                &mut assignment,
                min_breadth,
                min_chunk_breadth,
                safeguard_min_chunk_breadth,
            )?;
            evidence_decisions = resolved.decisions;
            evidence_selected_taxa = resolved.selected_taxa;
        }

        // Preserve only compact coordinate measurements for the optional raw
        // profile. Its separate EM runs later, after every production value is
        // finalized, so diagnostic work cannot influence production
        // aggregation or filtering. The large coordinate stores are never
        // cloned.
        let raw_profile_coverage = self.cfg.keep_raw.then(|| taxon_coverage.clone());

        let classification_path = format!("{}.classify.txt", self.cfg.outprefix);
        let database_is_fractional = self
            .cfg
            .db_metadata
            .as_ref()
            .map(|metadata| metadata.is_fractional)
            .unwrap_or(has_scaled_genomes);
        let write_classifications = self.cfg.save_classifications && !database_is_fractional;
        let mut resolved_classifications = ResolvedClassifications::new(write_classifications);
        if self.cfg.save_classifications {
            if database_is_fractional {
                log::warn!(
                    target: "PROFILER",
                    "Cannot write classify.txt for a fractional database because per-read assignments cannot be expanded safely to full-index equivalents."
                );
            } else {
                log::info!(
                    target: "PROFILER",
                    "Writing per-read classifications to {}",
                    classification_path
                );
            }
        }

        let mut unambiguous: HashMap<u32, TaxonTotals> = HashMap::new();
        let mut ambiguous_counts: HashMap<Vec<u32>, AmbiguousTotals> = HashMap::new();
        let mut ambiguous_reads: HashMap<Vec<u32>, Vec<String>> = HashMap::new();

        for assignment in &assignment.read_assignments {
            let mut final_taxids = Vec::new();
            let mut final_depths = Vec::new();
            let mut final_maplens = Vec::new();
            let mut final_scales = Vec::new();
            for (((taxid, depth), maplen), sampling_scale) in assignment
                .taxids
                .iter()
                .zip(&assignment.depths)
                .zip(&assignment.mapped_lens)
                .zip(&assignment.sampling_scales)
            {
                let coverage = taxon_coverage.get(taxid).cloned().unwrap_or_default();
                let evidence_selected =
                    assignment.requires_terminal_evidence && evidence_selected_taxa.contains(taxid);
                let legacy_selected = !assignment.requires_terminal_evidence;
                if evidence_selected || (legacy_selected && coverage.breadth >= min_breadth) {
                    final_taxids.push(*taxid);
                    final_depths.push(*depth);
                    final_maplens.push(*maplen);
                    final_scales.push(*sampling_scale);
                }
            }
            match final_taxids.len() {
                0 => {
                    continue;
                }
                1 => {
                    let taxid = final_taxids[0];
                    let entry = unambiguous.entry(taxid).or_default();
                    entry.support.add_weighted(1.0, final_scales[0], 1.0);
                    entry.depth += final_depths[0];
                    entry.mapped_len += final_maplens[0];
                    resolved_classifications.record_read(taxid, &assignment.read);
                }
                _ => {
                    if write_classifications {
                        ambiguous_reads
                            .entry(final_taxids.clone())
                            .or_default()
                            .push(assignment.read.to_string());
                    }
                    ambiguous_counts
                        .entry(final_taxids.clone())
                        .and_modify(|entry| {
                            entry.observed_count += 1.0;
                            accumulate(&mut entry.depths, &final_depths);
                            accumulate(&mut entry.mapped_lens, &final_maplens);
                            accumulate(&mut entry.equivalent_support_sums, &final_scales);
                        })
                        .or_insert_with(|| AmbiguousTotals {
                            observed_count: 1.0,
                            depths: final_depths.clone(),
                            mapped_lens: final_maplens.clone(),
                            equivalent_support_sums: final_scales.clone(),
                        });
                }
            }
        }

        let direct_taxa_before_lca: HashSet<u32> = unambiguous.keys().copied().collect();
        let mut taxids_with_lca: Vec<Vec<u32>> = Vec::new();
        let mut lca_ratio_lower_bounds: HashMap<u32, (f64, f64)> = HashMap::new();
        // Multiple ambiguous candidate groups may fold into the same LCA. Keep
        // one globally strongest coherent descendant witness for that LCA
        // instead of letting HashMap iteration order choose the last group.
        let mut coherent_lca_witnesses: HashMap<u32, (u32, WitnessDecision)> = HashMap::new();
        // Compatibility SAM input has no native evidence decisions. Keep its
        // established per-field maxima in a separate map so an LCA group cannot
        // overwrite coverage already accumulated for the taxon itself.
        let mut legacy_lca_coverage: HashMap<u32, TaxonCoverage> = HashMap::new();
        let mut taxids_with_any_ambiguous: Vec<Vec<u32>> = Vec::new();
        let mut ambiguous_removed: HashMap<Vec<u32>, AmbiguousTotals> = HashMap::new();

        // Iterate compact candidate groups in lexical taxid order. The values
        // stay borrowed, so reproducibility costs one pointer vector rather
        // than cloning every candidate/value vector.
        let mut initial_group_keys = ambiguous_counts.keys().collect::<Vec<_>>();
        initial_group_keys.sort_unstable();
        for taxids in initial_group_keys {
            let totals = &ambiguous_counts[taxids];
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
                let mut filtered_scale_sums = Vec::new();
                for (idx, taxid) in taxids.iter().enumerate() {
                    if presence[idx] {
                        filtered_taxids.push(*taxid);
                        filtered_depths.push(totals.depths[idx]);
                        filtered_maplens.push(totals.mapped_lens[idx]);
                        filtered_scale_sums.push(totals.equivalent_support_sums[idx]);
                    }
                }
                if filtered_taxids.len() == 1 {
                    let taxid = filtered_taxids[0];
                    let entry = unambiguous.entry(taxid).or_default();
                    entry
                        .support
                        .add_weighted(totals.observed_count, filtered_scale_sums[0], 1.0);
                    entry.depth += filtered_depths[0];
                    entry.mapped_len += filtered_maplens[0];
                    resolved_classifications.record_direct_group(
                        taxid,
                        totals.observed_count.round() as usize,
                        ambiguous_reads.remove(taxids),
                    );
                } else if !filtered_taxids.is_empty() {
                    ambiguous_removed
                        .entry(filtered_taxids.clone())
                        .and_modify(|entry| {
                            entry.observed_count += totals.observed_count;
                            accumulate(&mut entry.depths, &filtered_depths);
                            accumulate(&mut entry.mapped_lens, &filtered_maplens);
                            accumulate(&mut entry.equivalent_support_sums, &filtered_scale_sums);
                        })
                        .or_insert_with(|| AmbiguousTotals {
                            observed_count: totals.observed_count,
                            depths: filtered_depths.clone(),
                            mapped_lens: filtered_maplens.clone(),
                            equivalent_support_sums: filtered_scale_sums.clone(),
                        });
                    if let Some(mut reads) = ambiguous_reads.remove(taxids) {
                        ambiguous_reads
                            .entry(filtered_taxids)
                            .or_default()
                            .append(&mut reads);
                    }
                }
            } else {
                taxids_with_lca.push(taxids.clone());
                let taxid_set: HashSet<u32> = taxids.iter().copied().collect();
                if let Some((lca_taxid, _lca_name, lca_rank)) =
                    self.taxonomy.get_majority_lca(&taxid_set, 0.7)
                {
                    if let Some(bounds) = standard_lca_ratio_lower_bounds(&lca_rank, lca_ratio_mode)
                    {
                        lca_ratio_lower_bounds.insert(lca_taxid, bounds);
                    }
                    // Carry the coverage of a single coherent descendant
                    // witness. Maximizing breadth, fixed chunk and flex chunk
                    // independently could take each from a different genome,
                    // and the reporting path then divided those by an
                    // expectation rebuilt from the group's pooled depth and
                    // count -- so OEBR/COEBR compared one genome's observation
                    // against another population's expectation. Choosing one
                    // witness keeps observed, expected and both p-values on the
                    // same coordinate system.
                    let coherent_witness = coherent_lca_witness(taxids, &evidence_decisions);
                    if let Some(witness) = coherent_witness {
                        retain_coherent_lca_witness(
                            &mut coherent_lca_witnesses,
                            lca_taxid,
                            witness,
                        );
                    } else {
                        // No member has coordinate evidence (the compatibility
                        // SAM path has no witnesses at all). Retain the
                        // established maxima over legacy coverage.
                        let pick = |select: fn(&TaxonCoverage) -> f64| {
                            taxids
                                .iter()
                                .filter_map(|taxid| taxon_coverage.get(taxid).map(select))
                                .fold(0.0, f64::max)
                        };
                        let legacy = TaxonCoverage {
                            breadth: pick(|cov| cov.breadth),
                            fixed_chunk: pick(|cov| cov.fixed_chunk),
                            flex_chunk: pick(|cov| cov.flex_chunk),
                        };
                        merge_legacy_lca_coverage(&mut legacy_lca_coverage, lca_taxid, legacy);
                    }
                    let depth_mean = totals.depths.iter().sum::<f64>() / totals.depths.len() as f64;
                    let maplen_mean =
                        totals.mapped_lens.iter().sum::<f64>() / totals.mapped_lens.len() as f64;
                    let equivalent_mean = totals.equivalent_support_sums.iter().sum::<f64>()
                        / totals.equivalent_support_sums.len() as f64;
                    let entry = unambiguous.entry(lca_taxid).or_default();
                    entry
                        .support
                        .add_weighted(totals.observed_count, equivalent_mean, 1.0);
                    entry.depth += depth_mean;
                    entry.mapped_len += maplen_mean;
                    resolved_classifications.record_group(
                        lca_taxid,
                        taxids.clone(),
                        totals.observed_count.round() as usize,
                        ambiguous_reads.remove(taxids),
                    );
                }
            }
        }

        let lca_metrics: HashMap<u32, EvidenceCoverage> = coherent_lca_witnesses
            .into_iter()
            .map(|(lca_taxid, (_, decision))| (lca_taxid, decision.metrics))
            .collect();

        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "Majority-LCA folding: {} ambiguous candidate group(s) folded; {} distinct LCA taxa carried a coherent descendant witness.",
                taxids_with_lca.len(),
                lca_metrics.len(),
            );
        }
        for taxids in taxids_with_lca {
            ambiguous_counts.remove(&taxids);
            ambiguous_reads.remove(&taxids);
        }
        for taxids in taxids_with_any_ambiguous {
            ambiguous_counts.remove(&taxids);
            ambiguous_reads.remove(&taxids);
        }
        let mut ambiguous_removed = ambiguous_removed.into_iter().collect::<Vec<_>>();
        ambiguous_removed.sort_unstable_by(|left, right| left.0.cmp(&right.0));
        for (taxids, values) in ambiguous_removed {
            ambiguous_counts
                .entry(taxids)
                .and_modify(|entry| {
                    entry.observed_count += values.observed_count;
                    accumulate(&mut entry.depths, &values.depths);
                    accumulate(&mut entry.mapped_lens, &values.mapped_lens);
                    accumulate(
                        &mut entry.equivalent_support_sums,
                        &values.equivalent_support_sums,
                    );
                })
                .or_insert(values);
        }

        let mut species_count_depth = unambiguous.clone();
        let mut ambiguous_fraction: HashMap<Vec<u32>, HashMap<u32, f64>> = HashMap::new();
        let mut em_groups = ambiguous_counts.iter().collect::<Vec<_>>();
        em_groups.sort_unstable_by(|left, right| left.0.cmp(right.0));
        let mut em_taxids = unambiguous.keys().copied().collect::<Vec<_>>();
        em_taxids.sort_unstable();
        let mut iterations = 0;
        loop {
            let mut total = unambiguous.clone();
            for &(taxids, ambiguous) in &em_groups {
                let sum_depth: f64 = taxids
                    .iter()
                    .map(|taxid| {
                        species_count_depth
                            .get(taxid)
                            .map(|totals| totals.depth)
                            .unwrap_or(0.0)
                    })
                    .sum();
                if sum_depth <= f64::EPSILON {
                    continue;
                }
                let mut paired: Vec<(u32, f64, f64, f64, f64)> = taxids
                    .iter()
                    .enumerate()
                    .map(|(idx, taxid)| {
                        let fraction = species_count_depth
                            .get(taxid)
                            .map(|totals| totals.depth / sum_depth)
                            .unwrap_or(0.0);
                        (
                            *taxid,
                            ambiguous.depths[idx],
                            ambiguous.mapped_lens[idx],
                            ambiguous.equivalent_support_sums[idx],
                            fraction,
                        )
                    })
                    .collect();
                paired.sort_by(|left, right| {
                    right
                        .4
                        .total_cmp(&left.4)
                        .then_with(|| left.0.cmp(&right.0))
                });
                // `entry` would take ownership of the key, cloning the taxid
                // vector on every EM iteration for every ambiguous group even
                // though the group is only ever inserted once. Look the key up
                // by reference and clone it only on the first iteration.
                if !ambiguous_fraction.contains_key(taxids) {
                    ambiguous_fraction.insert(taxids.clone(), HashMap::new());
                }
                let entry = ambiguous_fraction
                    .get_mut(taxids)
                    .expect("ambiguous fraction entry was just inserted");
                entry.clear();
                let mut consumed = false;
                for (taxid, depth, maplen, equivalent_sum, fraction) in paired {
                    let total_entry = total.entry(taxid).or_default();
                    if fraction >= 0.99 {
                        total_entry.support.add_weighted(
                            ambiguous.observed_count,
                            equivalent_sum,
                            1.0,
                        );
                        total_entry.depth += depth;
                        total_entry.mapped_len += maplen;
                        entry.insert(taxid, 1.0);
                        consumed = true;
                        break;
                    }
                    total_entry.support.add_weighted(
                        ambiguous.observed_count,
                        equivalent_sum,
                        fraction,
                    );
                    total_entry.depth += depth * fraction;
                    total_entry.mapped_len += maplen * fraction;
                    entry.insert(taxid, round5(fraction));
                }
                if !consumed {
                    normalize_map(entry);
                }
            }
            if verbose_logging && iterations % 5 == 0 {
                log::info!(
                    target: "PROFILER",
                    "Profiling in EM iteration {} ...",
                    iterations
                );
            }
            let numerator: f64 = em_taxids
                .iter()
                .map(|taxid| {
                    let depth = total.get(taxid).map_or(0.0, |totals| totals.depth);
                    let previous_depth = species_count_depth
                        .get(taxid)
                        .map(|previous| previous.depth)
                        .unwrap_or(0.0);
                    (depth - previous_depth).abs()
                })
                .sum::<f64>();
            let denominator: f64 = em_taxids
                .iter()
                .map(|taxid| species_count_depth.get(taxid).map_or(0.0, |totals| totals.depth))
                .sum::<f64>()
                .max(1e-12);
            if numerator / denominator <= 1e-5 {
                species_count_depth = total;
                break;
            }
            iterations += 1;
            if iterations >= self.cfg.max_em_iterations {
                log::warn!(
                    target: "PROFILER",
                    "Community EM reached the maximum of {} iterations before convergence.",
                    self.cfg.max_em_iterations
                );
                species_count_depth = total;
                break;
            }
            species_count_depth = total;
        }
        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "Finished profiling with {} EM iterations",
                iterations
            );
        }

        let norm_factor: f64 = 0.01
            * em_taxids
                .iter()
                .map(|taxid| species_count_depth.get(taxid).map_or(0.0, |totals| totals.depth))
                .sum::<f64>();
        let num_taxa = species_count_depth.len();
        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "Legacy/LCA p-value hypothesis family: {} post-EM taxa.",
                num_taxa,
            );
        }

        let mut out_taxa_list = Vec::new();
        let host_taxid = self.cfg.host.as_deref().and_then(|h| h.parse::<u32>().ok());
        let pathogen_map = if self.cfg.host.is_some() {
            self.cfg
                .pathogen_host
                .as_deref()
                .map(load_pathogen_table)
                .transpose()?
        } else {
            None
        };
        let mut entries: Vec<(u32, TaxonTotals)> = species_count_depth
            .iter()
            .map(|(taxid, totals)| (*taxid, *totals))
            .collect();
        entries.sort_by(|left, right| {
            right
                .1
                .depth
                .total_cmp(&left.1.depth)
                .then_with(|| left.0.cmp(&right.0))
        });

        for (taxid, totals) in entries {
            let mapcount = totals.support.observed;
            if mapcount < 0.5 {
                continue;
            }
            let equivalent_mapcount = totals.support.equivalent;
            let depth = totals.depth;
            let maplen = totals.mapped_len;
            let abundance = if norm_factor <= f64::EPSILON {
                0.0
            } else {
                depth / norm_factor
            };
            let rank = terminal_ranks
                .get(&taxid)
                .copied()
                .map(ProfileRank::as_str)
                .or_else(|| self.taxonomy.get_rank(taxid))
                .unwrap_or("no rank")
                .to_string();
            let taxname = self
                .taxonomy
                .get_name(taxid)
                .unwrap_or("Taxonomy deprecated");
            if taxname == "Taxonomy deprecated" {
                continue;
            }
            let evidence_metrics = if evidence_selected_taxa.contains(&taxid) {
                evidence_decisions
                    .get(&taxid)
                    .map(|decision| decision.metrics)
            } else {
                // A majority-LCA row reports the coherent coverage of its chosen
                // descendant witness instead of dividing pooled observations by
                // an expectation rebuilt from pooled depth and count. Report
                // filtering still applies to LCA rows; only the numbers change.
                lca_metrics.get(&taxid).copied()
            };
            let (
                breadth,
                fixed_chunk,
                flex_chunk,
                expected_breadth,
                cov_prob,
                expected_flex_chunk,
                chunkcov_prob,
            ) = if let Some(metrics) = evidence_metrics {
                // Presence statistics are all derived from the same resolved,
                // assembly-coherent coordinate set. EM depth remains solely
                // the abundance estimate and must not redefine its expected
                // breadth or p-values.
                (
                    metrics.breadth,
                    metrics.fixed_chunk,
                    metrics.flex_chunk,
                    metrics.expected_breadth,
                    metrics.breadth_pvalue,
                    metrics.expected_flex_chunk,
                    metrics.chunk_pvalue,
                )
            } else {
                let coverage = legacy_output_coverage(
                    taxid,
                    &direct_taxa_before_lca,
                    &taxon_coverage,
                    &legacy_lca_coverage,
                );
                let breadth = coverage.breadth.min(1.0);
                let fixed_chunk = coverage.fixed_chunk.min(1.0);
                let flex_chunk = coverage.flex_chunk.min(1.0);
                let expected_breadth = crate::profiler::evidence::expected_breadth(depth, mapcount);
                let cov_prob = if rank == "species" {
                    coverage_deviation_pvalue(breadth, expected_breadth, Some(depth), mapcount, num_taxa)
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
                    (genome_len.sqrt() / 2.0).floor().max(1.0) as usize
                };
                let expected_flex_chunk = 1.0 - (1.0 - 1.0 / num_flex_chunks as f64).powf(mapcount);
                let chunkcov_prob = if rank == "species" {
                    coverage_deviation_pvalue(flex_chunk, expected_flex_chunk, None, mapcount, num_taxa)
                } else {
                    None
                };
                (
                    breadth,
                    fixed_chunk,
                    flex_chunk,
                    expected_breadth,
                    cov_prob,
                    expected_flex_chunk,
                    chunkcov_prob,
                )
            };
            let annotations = if let Some(host_taxid) = host_taxid {
                if let Some(map) = pathogen_map.as_ref() {
                    if let Some(entry) = map.get(&taxid) {
                        let (report_it, host_names, host_taxids, diseases) =
                            self.annotate_host_disease(entry, host_taxid);
                        if !report_it {
                            continue;
                        }
                        Some(ProfileAnnotations {
                            host_names: host_names.unwrap_or_else(|| "NA".to_string()),
                            host_taxids: host_taxids.unwrap_or_else(|| "NA".to_string()),
                            diseases: diseases.unwrap_or_else(|| "NA".to_string()),
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
                None
            };
            out_taxa_list.push(ProfileRow {
                taxname: taxname.to_string(),
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

        let mut normalized: Vec<ProfileRow> = out_taxa_list
            .into_iter()
            .filter(|row| {
                community_profile_row_passes(
                    row,
                    &self.cfg,
                    min_chunk_breadth,
                    &evidence_selected_taxa,
                    &lca_ratio_lower_bounds,
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
        let final_classified_reads = write_community_classifications(
            write_classifications.then_some(classification_path.as_str()),
            &resolved_classifications,
            &ambiguous_counts,
            &ambiguous_reads,
            &ambiguous_fraction,
            &retained_taxa,
            &terminal_ranks,
            &self.taxonomy,
        )?;

        // Production assignment, EM, report filtering, and abundance
        // normalization are complete. Release their grouped working sets
        // before the opt-in shadow EM to avoid overlapping the two memory
        // peaks. Keep only per-read candidate vectors needed by raw EM, then
        // release the remaining coordinate-heavy AssignmentResult fields.
        drop(ambiguous_counts);
        drop(ambiguous_reads);
        drop(ambiguous_fraction);
        drop(unambiguous);
        drop(species_count_depth);
        let raw_assignments = raw_profile_coverage
            .as_ref()
            .map(|_| std::mem::take(&mut assignment.read_assignments));
        drop(assignment);

        let raw_taxa_list = if let (Some(mut coverage), Some(assignments)) =
            (raw_profile_coverage, raw_assignments)
        {
            let totals = estimate_unfiltered_taxon_totals(
                &assignments,
                &self.taxonomy,
                &mut coverage,
                self.cfg.max_em_iterations,
                verbose_logging,
            )?;
            self.build_raw_profile_rows(
                &RawProfileEstimate { totals, coverage },
                &evidence_decisions,
                &terminal_ranks,
                host_taxid,
                pathogen_map.as_ref(),
            )
        } else {
            Vec::new()
        };

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
            row.write(&mut profile_writer)?;
        }
        profile_writer.flush()?;

        if self.cfg.keep_raw {
            let raw_path = format!("{}.rprofile.txt", out_prefix);
            let mut raw_writer = WriterBuilder::new()
                .has_headers(false)
                .delimiter(b'\t')
                .from_writer(File::create(&raw_path)?);
            for row in raw_taxa_list {
                row.write(&mut raw_writer)?;
            }
            raw_writer.flush()?;
        }

        Ok(final_classified_reads)
    }

    /// Format the optional raw community profile from its independent EM
    /// estimate.
    ///
    /// Unlike the production table, these rows are not filtered by minimum
    /// read support, breadth/chunk amount, OEBR/COEBR, p-values, or the
    /// coordinate-evidence accept/reject decision. Mapping-quality gates,
    /// taxonomy normalization, ambiguity/LCA handling, and EM refinement have
    /// already been applied. This makes the table useful for reproducible
    /// downstream threshold exploration without changing production calls.
    fn build_raw_profile_rows(
        &self,
        estimate: &RawProfileEstimate,
        evidence_decisions: &HashMap<u32, WitnessDecision>,
        terminal_ranks: &HashMap<u32, ProfileRank>,
        host_taxid: Option<u32>,
        pathogen_map: Option<&HashMap<u32, PathogenEntry>>,
    ) -> Vec<ProfileRow> {
        let norm_factor = 0.01
            * estimate
                .totals
                .values()
                .map(|totals| totals.depth)
                .sum::<f64>();
        let num_taxa = estimate.totals.len().max(1);
        let mut entries = estimate
            .totals
            .iter()
            .map(|(&taxid, &totals)| (taxid, totals))
            .collect::<Vec<_>>();
        entries.sort_by(|left, right| {
            right
                .1
                .depth
                .partial_cmp(&left.1.depth)
                .unwrap_or(Ordering::Equal)
                .then_with(|| left.0.cmp(&right.0))
        });

        let mut rows = Vec::with_capacity(entries.len());
        for (taxid, totals) in entries {
            let mapcount = totals.support.observed;
            let equivalent_mapcount = totals.support.equivalent;
            let depth = totals.depth;
            if mapcount <= f64::EPSILON
                && equivalent_mapcount <= f64::EPSILON
                && depth <= f64::EPSILON
            {
                continue;
            }
            let maplen = totals.mapped_len;
            let abundance = if norm_factor <= f64::EPSILON {
                0.0
            } else {
                depth / norm_factor
            };
            let rank = terminal_ranks
                .get(&taxid)
                .copied()
                .map(ProfileRank::as_str)
                .or_else(|| self.taxonomy.get_rank(taxid))
                .unwrap_or("no rank")
                .to_string();
            let Some(taxname) = self.taxonomy.get_name(taxid) else {
                continue;
            };

            let evidence_metrics = evidence_decisions
                .get(&taxid)
                .map(|decision| decision.metrics);
            let (
                breadth,
                fixed_chunk,
                flex_chunk,
                expected_breadth,
                cov_prob,
                expected_flex_chunk,
                chunkcov_prob,
            ) = if let Some(metrics) = evidence_metrics {
                (
                    metrics.breadth,
                    metrics.fixed_chunk,
                    metrics.flex_chunk,
                    metrics.expected_breadth,
                    metrics.breadth_pvalue,
                    metrics.expected_flex_chunk,
                    metrics.chunk_pvalue,
                )
            } else {
                let coverage = estimate.coverage.get(&taxid).cloned().unwrap_or_default();
                let breadth = coverage.breadth.min(1.0);
                let fixed_chunk = coverage.fixed_chunk.min(1.0);
                let flex_chunk = coverage.flex_chunk.min(1.0);
                let expected_breadth = crate::profiler::evidence::expected_breadth(depth, mapcount);
                let cov_prob = if rank == "species" {
                    coverage_deviation_pvalue(breadth, expected_breadth, Some(depth), mapcount, num_taxa)
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
                    (genome_len.sqrt() / 2.0).floor().max(1.0) as usize
                };
                let expected_flex_chunk = 1.0 - (1.0 - 1.0 / num_flex_chunks as f64).powf(mapcount);
                let chunkcov_prob = if rank == "species" {
                    coverage_deviation_pvalue(flex_chunk, expected_flex_chunk, None, mapcount, num_taxa)
                } else {
                    None
                };
                (
                    breadth,
                    fixed_chunk,
                    flex_chunk,
                    expected_breadth,
                    cov_prob,
                    expected_flex_chunk,
                    chunkcov_prob,
                )
            };

            let annotations = if let Some(host_taxid) = host_taxid {
                if let Some(entry) = pathogen_map.and_then(|map| map.get(&taxid)) {
                    let (report_it, host_names, host_taxids, diseases) =
                        self.annotate_host_disease(entry, host_taxid);
                    if !report_it {
                        continue;
                    }
                    Some(ProfileAnnotations {
                        host_names: host_names.unwrap_or_else(|| "NA".to_string()),
                        host_taxids: host_taxids.unwrap_or_else(|| "NA".to_string()),
                        diseases: diseases.unwrap_or_else(|| "NA".to_string()),
                    })
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
            rows.push(ProfileRow {
                taxname: taxname.to_string(),
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
        rows
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
/// alignments. `best` holds `(identity, depth, mapped_len, sampling_scale)`;
/// `order` keeps first-seen order for deterministic output; `intervals` captures the
/// coverage interval triple from the currently-best alignment.
type BestMap = HashMap<u32, (f64, f64, f64, f64)>;
type IntervalMap = HashMap<u32, (AccKey, (usize, usize), (usize, usize), (usize, usize))>;

struct CoordinateEvidenceResult {
    decisions: HashMap<u32, WitnessDecision>,
    selected_taxa: HashSet<u32>,
}

/// Choose the one descendant witness whose coverage a majority-LCA row should
/// report.
///
/// A majority-LCA row has no coordinates of its own. Maximizing breadth, fixed
/// chunk and flex chunk independently across descendants could take each from a
/// different genome, and the report path then divides those by an expectation
/// rebuilt from the group's pooled depth and count -- comparing one genome's
/// observation against another population's expectation. Reporting a single
/// witness keeps observed, expected, OEBR, COEBR and both p-values on one
/// coordinate system.
///
/// Selection is the codebase's existing witness preference, made total with a
/// lowest-taxid tiebreak so the choice cannot depend on candidate order.
/// Returns `None` when no member has coordinate evidence at all, which is the
/// normal case for the compatibility SAM path.
fn coherent_lca_witness(
    taxids: &[u32],
    evidence_decisions: &HashMap<u32, WitnessDecision>,
) -> Option<(u32, WitnessDecision)> {
    taxids
        .iter()
        .filter_map(|taxid| {
            evidence_decisions
                .get(taxid)
                .map(|decision| (*taxid, *decision))
        })
        .reduce(coherent_witness_max)
}

/// Return the preferred coherent witness using a total, order-independent
/// comparison. The taxid tiebreak is required because `witness_is_better`
/// intentionally compares evidence strength only.
#[inline]
fn coherent_witness_max(
    left: (u32, WitnessDecision),
    right: (u32, WitnessDecision),
) -> (u32, WitnessDecision) {
    let right_is_better = witness_is_better(right.1, left.1);
    let equivalent = !right_is_better && !witness_is_better(left.1, right.1);
    if right_is_better || (equivalent && right.0 < left.0) {
        right
    } else {
        left
    }
}

/// Retain one coherent coordinate witness across every ambiguous group that
/// resolves to the same LCA. This is O(1) per folded group and cannot depend on
/// the randomized iteration order of `ambiguous_counts`.
fn retain_coherent_lca_witness(
    witnesses: &mut HashMap<u32, (u32, WitnessDecision)>,
    lca_taxid: u32,
    candidate: (u32, WitnessDecision),
) {
    match witnesses.entry(lca_taxid) {
        Entry::Vacant(entry) => {
            entry.insert(candidate);
        }
        Entry::Occupied(mut entry) => {
            let retained = coherent_witness_max(*entry.get(), candidate);
            entry.insert(retained);
        }
    }
}

/// Merge compatibility-SAM LCA coverage without mutating direct taxon
/// coverage. This preserves the established component-wise maxima while making
/// the result independent of ambiguous-group iteration order.
fn merge_legacy_lca_coverage(
    coverage: &mut HashMap<u32, TaxonCoverage>,
    lca_taxid: u32,
    candidate: TaxonCoverage,
) {
    let retained = coverage.entry(lca_taxid).or_default();
    retained.breadth = retained.breadth.max(candidate.breadth);
    retained.fixed_chunk = retained.fixed_chunk.max(candidate.fixed_chunk);
    retained.flex_chunk = retained.flex_chunk.max(candidate.flex_chunk);
}

/// Select output coverage for the legacy path. Evidence directly assigned to
/// a taxon takes precedence; otherwise an LCA row uses the deterministic
/// descendant fallback, then ordinary taxon coverage.
fn legacy_output_coverage(
    taxid: u32,
    direct_taxa_before_lca: &HashSet<u32>,
    direct: &HashMap<u32, TaxonCoverage>,
    lca: &HashMap<u32, TaxonCoverage>,
) -> TaxonCoverage {
    if direct_taxa_before_lca.contains(&taxid) {
        return direct.get(&taxid).cloned().unwrap_or_default();
    }
    lca.get(&taxid)
        .or_else(|| direct.get(&taxid))
        .cloned()
        .unwrap_or_default()
}

fn strongest_witnesses(
    evaluations: &[(PresenceKey, WitnessDecision)],
) -> HashMap<u32, WitnessDecision> {
    let mut strongest: HashMap<u32, (PresenceKey, WitnessDecision)> = HashMap::new();
    for &(key, decision) in evaluations {
        match strongest.entry(key.taxid) {
            Entry::Vacant(entry) => {
                entry.insert((key, decision));
            }
            Entry::Occupied(mut entry) => {
                if keyed_witness_is_better((key, decision), *entry.get()) {
                    entry.insert((key, decision));
                }
            }
        }
    }
    strongest
        .into_iter()
        .map(|(taxid, (_, decision))| (taxid, decision))
        .collect()
}

/// Compare coherent assembly witnesses with a stable final key. The policy
/// comparison remains authoritative; assembly identity is consulted only when
/// neither witness is stronger. This is constant-time and avoids sorting all
/// witnesses merely to make equal-evidence output reproducible.
#[inline]
fn keyed_witness_is_better(
    left: (PresenceKey, WitnessDecision),
    right: (PresenceKey, WitnessDecision),
) -> bool {
    if witness_is_better(left.1, right.1) {
        return true;
    }
    if witness_is_better(right.1, left.1) {
        return false;
    }
    left.0 < right.0
}

/// Select MAP-only species additions without allowing MAP evidence to satisfy
/// its own independent-coverage guard.
///
/// The lookup is keyed by both species and assembly. This is intentional: a
/// broad pre-MAP witness on one assembly cannot authorize MAP-created
/// coverage on another assembly of the same species.
fn guarded_map_additions(
    baseline_selected: &HashSet<u32>,
    baseline_by_key: &HashMap<PresenceKey, WitnessDecision>,
    post_evaluations: &[(PresenceKey, WitnessDecision)],
    min_chunk_breadth: f64,
) -> HashMap<u32, WitnessDecision> {
    let mut additions: HashMap<u32, (PresenceKey, WitnessDecision)> = HashMap::new();
    for &(key, post) in post_evaluations {
        if !key.profile_rank.is_species() {
            continue;
        }
        if baseline_selected.contains(&key.taxid) {
            continue;
        }
        let Some(premap_fraction) = post.required_premap_guard.fraction() else {
            continue;
        };
        let Some(pre) = baseline_by_key.get(&key) else {
            continue;
        };
        let premap_threshold = min_chunk_breadth * premap_fraction;
        if pre.metrics.fixed_chunk.partial_cmp(&premap_threshold) != Some(Ordering::Greater) {
            continue;
        }
        match additions.entry(key.taxid) {
            Entry::Vacant(entry) => {
                entry.insert((key, post));
            }
            Entry::Occupied(mut entry) => {
                if keyed_witness_is_better((key, post), *entry.get()) {
                    entry.insert((key, post));
                }
            }
        }
    }
    additions
        .into_iter()
        .map(|(taxid, (_, decision))| (taxid, decision))
        .collect()
}

/// Build the seeded prior used by guarded short-read exact-tie MAP.
///
/// Only species with at least one unique-candidate fragment are eligible and
/// only their unique-candidate read/mate-1 depth seeds the prior. Every ordinary
/// ambiguous candidate group then updates that prior by fractional EM. This
/// is deliberately separate from the final abundance EM: MAP can contribute
/// coordinate evidence, but it never hard-assigns abundance.
fn seeded_short_read_map_priors(
    assignments: &[ReadAssignment],
    eligible_species: &HashSet<u32>,
    max_iterations: usize,
) -> (HashMap<u32, f64>, usize) {
    if eligible_species.is_empty() {
        return (HashMap::new(), 0);
    }
    let mut species = eligible_species.iter().copied().collect::<Vec<_>>();
    species.sort_unstable();
    let indices = species
        .iter()
        .enumerate()
        .map(|(index, &taxid)| (taxid, index))
        .collect::<HashMap<_, _>>();
    let mut base = vec![0.0f64; species.len()];
    let mut grouped: HashMap<Vec<usize>, Vec<f64>> = HashMap::new();

    for assignment in assignments {
        if assignment.rank != AssignRank::Species || assignment.taxids.is_empty() {
            continue;
        }
        if assignment.taxids.len() == 1 {
            if let Some(&index) = indices.get(&assignment.taxids[0]) {
                let depth = assignment.depths.first().copied().unwrap_or(0.0);
                if depth.is_finite() && depth > 0.0 {
                    base[index] += depth;
                }
            }
            continue;
        }

        let mut key = Vec::with_capacity(assignment.taxids.len());
        let mut depths = Vec::with_capacity(assignment.taxids.len());
        for (&taxid, &depth) in assignment.taxids.iter().zip(&assignment.depths) {
            if let Some(&index) = indices.get(&taxid) {
                key.push(index);
                depths.push(if depth.is_finite() && depth > 0.0 {
                    depth
                } else {
                    0.0
                });
            }
        }
        if key.is_empty() {
            continue;
        }
        if let Some(total) = grouped.get_mut(key.as_slice()) {
            for (target, value) in total.iter_mut().zip(depths) {
                *target += value;
            }
        } else {
            grouped.insert(key, depths);
        }
    }

    let mut groups = grouped.into_iter().collect::<Vec<_>>();
    groups.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    let mut current = base.clone();
    let mut iterations = 0usize;
    for iteration in 0..max_iterations.max(1) {
        let mut next = base.clone();
        for (candidate_indices, depths) in &groups {
            if candidate_indices.len() == 1 {
                next[candidate_indices[0]] += depths[0];
                continue;
            }
            let denominator = candidate_indices
                .iter()
                .map(|&index| current[index])
                .sum::<f64>();
            if denominator <= f64::EPSILON || !denominator.is_finite() {
                continue;
            }
            for (&index, &depth) in candidate_indices.iter().zip(depths) {
                next[index] += depth * (current[index] / denominator);
            }
        }
        let numerator = next
            .iter()
            .zip(&current)
            .map(|(next, current)| (next - current).abs())
            .sum::<f64>();
        let denominator = current.iter().sum::<f64>().max(1e-12);
        current = next;
        iterations = iteration + 1;
        if numerator / denominator <= 1e-5 {
            break;
        }
    }

    let priors = species
        .into_iter()
        .zip(current)
        .filter(|(_, depth)| depth.is_finite() && *depth > 0.0)
        .collect();
    (priors, iterations)
}

/// Add one guarded-MAP placement and report the `presence` index it landed on,
/// so the caller can mark that witness's coverage as needing recomputation.
fn add_exact_tie_placement(
    presence: &mut Vec<(PresenceKey, PresenceAssembly)>,
    presence_indices: &mut HashMap<PresenceKey, usize>,
    taxid: u32,
    placement: ExactTiePlacement,
) -> Result<usize> {
    let key = PresenceKey {
        taxid,
        profile_rank: ProfileRank::Species,
        assembly_id: placement.assembly_id,
    };
    let index = if let Some(&index) = presence_indices.get(&key) {
        index
    } else {
        let index = presence.len();
        presence.push((key, PresenceAssembly::default()));
        presence_indices.insert(key, index);
        index
    };
    let assembly = &mut presence[index].1;
    if assembly.asm_len <= 0.0 {
        assembly.asm_len = placement.asm_len;
    } else {
        let tolerance = 1e-12 * assembly.asm_len.abs().max(placement.asm_len.abs()).max(1.0);
        if (assembly.asm_len - placement.asm_len).abs() > tolerance {
            bail!(
                "inconsistent effective assembly size for species {} assembly {}: {} and {}",
                taxid,
                placement.assembly_id,
                assembly.asm_len,
                placement.asm_len
            );
        }
    }
    assembly.count = assembly.count.saturating_add(1);
    assembly.reference_span_bases = assembly
        .reference_span_bases
        .saturating_add(u64::from(placement.span.1.saturating_sub(placement.span.0)));
    assembly.intervals.push(PresenceInterval {
        reference_id: placement.reference_id,
        span: placement.span,
        fixed_chunk: placement.fixed_chunk,
        flex_chunk: placement.flex_chunk,
    });
    Ok(index)
}

fn resolve_coordinate_evidence(
    cfg: &ProfilerConfig,
    policy: EvidencePolicy,
    assignment: &mut AssignmentResult,
    min_breadth: f64,
    standard_min_chunk_breadth: f64,
    safeguard_min_chunk_breadth: f64,
) -> Result<CoordinateEvidenceResult> {
    let thresholds = EvidenceThresholds {
        mode: policy.mode,
        min_breadth,
        standard_min_chunk_breadth,
        safeguard_min_chunk_breadth,
        chunk_filter_enabled: cfg.min_reads.is_none() || cfg.chunk_breadth.is_some(),
        min_oebr: cfg.min_oebr,
        min_coebr: cfg.min_coebr,
        strict: cfg.strict,
    };
    let mut presence = std::mem::take(&mut assignment.presence)
        .into_iter()
        .collect::<Vec<_>>();
    let baseline_witness_count = presence.len();
    // Correct for every hypothesis actually tested: one breadth test and one
    // chunk test run per assembly-coherent witness. `uniquely_supported_taxa`
    // was a proxy that matched neither the tests performed nor the taxa
    // reported, and it under-corrected -- leaving adjusted p-values smaller
    // than the testing burden justifies, so the gate rejected coherent taxa.
    // Both the baseline and post-MAP passes use this same family so the two
    // remain directly comparable.
    let multiple_testing_taxa = baseline_witness_count.max(1);
    if cfg.verbose || cfg.very_verbose {
        log::info!(
            target: "PROFILER",
            "Coverage p-value hypothesis family: {} assembly-coherent witnesses tested (uniquely supported taxa: {}).",
            multiple_testing_taxa,
            assignment.uniquely_supported_taxa.len(),
        );
    }
    let pool = ThreadPoolBuilder::new()
        .num_threads(cfg.threads.max(1))
        .build()?;
    let baseline_evaluations = pool.install(|| {
        presence
            .par_iter_mut()
            .map(|(key, assembly)| {
                let metrics = calculate_coverage(assembly, multiple_testing_taxa);
                let evidence = assignment
                    .taxon_evidence
                    .get(&key.taxid)
                    .copied()
                    .unwrap_or_default();
                (
                    *key,
                    evaluate_ranked_witness(metrics, evidence, thresholds, key.profile_rank),
                )
            })
            .collect::<Vec<_>>()
    });
    let baseline_decisions = strongest_witnesses(&baseline_evaluations);
    let baseline_selected = baseline_decisions
        .iter()
        .filter_map(|(&taxid, decision)| decision.accepted().then_some(taxid))
        .collect::<HashSet<_>>();

    if !policy.collect_exact_ties {
        if cfg.verbose || cfg.very_verbose {
            let standard = baseline_decisions
                .values()
                .filter(|decision| decision.standard)
                .count();
            let rescued = baseline_decisions
                .values()
                .filter(|decision| decision.rescued())
                .count();
            log::info!(
                target: "PROFILER",
                "Resolved coordinate evidence ({}) selected {} terminal taxa ({} standard, {} rescue) from {} assembly-coherent witnesses.",
                policy.mode.label(),
                baseline_selected.len(),
                standard,
                rescued,
                baseline_witness_count,
            );
        }
        return Ok(CoordinateEvidenceResult {
            decisions: baseline_decisions,
            selected_taxa: baseline_selected,
        });
    }

    // This is also checked by the CLI, but library callers receive a safe
    // baseline rather than an unguarded MAP pass if the amount guard is absent.
    if !policy.mode.supports_guarded_map()
        || !thresholds.chunk_filter_enabled
        || safeguard_min_chunk_breadth <= 0.0
    {
        log::warn!(
            target: "PROFILER",
            "Guarded short-read MAP requires native Illumina evidence and a positive active chunk-breadth cutoff; exact ties remain ambiguous."
        );
        return Ok(CoordinateEvidenceResult {
            decisions: baseline_decisions,
            selected_taxa: baseline_selected,
        });
    }

    let exact_ties = std::mem::take(&mut assignment.exact_ties);
    let exact_tie_count = exact_ties.fragment_count();
    if exact_tie_count == 0 {
        return Ok(CoordinateEvidenceResult {
            decisions: baseline_decisions,
            selected_taxa: baseline_selected,
        });
    }
    // Guarded MAP is a calibrated species-only policy. Higher-rank terminal
    // endpoints remain in abundance EM but must never seed or receive MAP
    // coordinate evidence.
    let map_eligible_species = presence
        .iter()
        .filter_map(|(key, _)| {
            (key.profile_rank.is_species()
                && assignment.uniquely_supported_taxa.contains(&key.taxid))
            .then_some(key.taxid)
        })
        .collect::<HashSet<_>>();
    let (priors, prior_iterations) = seeded_short_read_map_priors(
        &assignment.read_assignments,
        &map_eligible_species,
        cfg.max_em_iterations,
    );
    if priors.is_empty() {
        return Ok(CoordinateEvidenceResult {
            decisions: baseline_decisions,
            selected_taxa: baseline_selected,
        });
    }
    let baseline_by_key = baseline_evaluations
        .iter()
        .copied()
        .collect::<HashMap<_, _>>();
    let mut presence_indices = presence
        .iter()
        .enumerate()
        .map(|(index, (key, _))| (*key, index))
        .collect::<HashMap<_, _>>();
    // Update the cloned counters directly. This avoids two temporary per-taxon
    // hash maps while preserving the immutable baseline evidence.
    let mut mapped_taxon_evidence = assignment.taxon_evidence.clone();
    let mut resolved_ties = 0usize;
    let mut posterior_sum = 0.0f64;
    // Coverage is a pure function of an assembly's intervals, count and length,
    // so only witnesses that actually receive a MAP placement need their
    // intervals merged and sorted again. Track those; witnesses appended past
    // the baseline are new and always need a first computation.
    let mut coverage_dirty = vec![false; baseline_witness_count];
    for fragment in exact_ties.fragments() {
        let Some((winner, posterior)) = guarded_map_winner(fragment.tied(), &priors) else {
            continue;
        };
        resolved_ties += 1;
        posterior_sum += posterior;
        // Every candidate of this fragment was compared against the winner, so
        // the whole fragment now belongs to a resolved comparison — not only
        // the tied prefix. Crediting just the prefix left a taxon that lost the
        // tie outright with a too-small denominator, inflating the
        // resolved-best fraction that platform and low-abundance rescue read.
        for candidate in fragment.all() {
            let evidence = mapped_taxon_evidence.entry(candidate.taxid).or_default();
            evidence.resolved_candidate_fragments =
                evidence.resolved_candidate_fragments.saturating_add(1);
        }
        let winner_evidence = mapped_taxon_evidence.entry(winner.taxid).or_default();
        winner_evidence.best_fragments = winner_evidence.best_fragments.saturating_add(1);
        for &placement in exact_ties.placements(winner) {
            let index = add_exact_tie_placement(
                &mut presence,
                &mut presence_indices,
                winner.taxid,
                placement,
            )?;
            if let Some(dirty) = coverage_dirty.get_mut(index) {
                *dirty = true;
            }
        }
    }
    let post_evaluations = pool.install(|| {
        presence
            .par_iter_mut()
            .enumerate()
            .map(|(index, (key, assembly))| {
                // Reuse the baseline coverage unless this witness received new
                // MAP coordinates. `calculate_coverage` depends only on the
                // assembly contents and the (unchanged) hypothesis count, so
                // the reused value is exactly what recomputing would produce.
                let metrics = match coverage_dirty.get(index) {
                    Some(false) => baseline_evaluations[index].1.metrics,
                    _ => calculate_coverage(assembly, multiple_testing_taxa),
                };
                // The decision itself must always be re-derived: guarded MAP
                // updates evidence denominators for every candidate of a
                // resolved fragment, so a witness whose coordinates were never
                // touched can still change tier.
                let evidence = mapped_taxon_evidence
                    .get(&key.taxid)
                    .copied()
                    .unwrap_or_default();
                (
                    *key,
                    evaluate_ranked_witness(metrics, evidence, thresholds, key.profile_rank),
                )
            })
            .collect::<Vec<_>>()
    });
    let mut post_decisions = strongest_witnesses(&post_evaluations);
    let mut selected_taxa = baseline_selected.clone();
    let additions = guarded_map_additions(
        &baseline_selected,
        &baseline_by_key,
        &post_evaluations,
        safeguard_min_chunk_breadth,
    );
    let guarded_additions = additions.len();
    let low_abundance_additions = additions
        .values()
        .filter(|decision| decision.rescues.low_abundance_pe)
        .count();
    for (taxid, decision) in additions {
        selected_taxa.insert(taxid);
        // The globally strongest post-MAP witness may belong to an assembly
        // that did not satisfy the pre-MAP guard. Report the strongest
        // qualifying same-assembly witness instead.
        post_decisions.insert(taxid, decision);
    }
    // MAP is monotonic. If added conserved intervals make an upper-ratio or
    // p-value statistic fail, retain the baseline witness and its metrics.
    for &taxid in &baseline_selected {
        let post_still_passes = post_decisions
            .get(&taxid)
            .is_some_and(|decision| decision.accepted());
        if !post_still_passes {
            if let Some(&baseline) = baseline_decisions.get(&taxid) {
                post_decisions.insert(taxid, baseline);
            }
        }
    }

    if cfg.verbose || cfg.very_verbose {
        let posterior_mean = if resolved_ties == 0 {
            0.0
        } else {
            posterior_sum / resolved_ties as f64
        };
        let premap_guard = match policy.mode {
            EvidenceMode::Paired => 0.5 * safeguard_min_chunk_breadth,
            EvidenceMode::Single => safeguard_min_chunk_breadth,
            EvidenceMode::Long => unreachable!("long-read MAP is disabled"),
        };
        let post_requirement = match policy.mode {
            EvidenceMode::Paired => "standard-or-rescue",
            EvidenceMode::Single => "standard-only",
            EvidenceMode::Long => unreachable!("long-read MAP is disabled"),
        };
        log::info!(
            target: "PROFILER",
            "Guarded {} exact-tie MAP: eligible_species={} prior_iterations={} exact_ties={} resolved={} mean_posterior={:.5}; baseline_species={} guarded_additions={} (low_abundance={}) final_species={} (ordinary same-assembly pre-MAP fixed_chunk > {:.6}; PE low-abundance guard > {:.6}; post-MAP {}).",
            policy.mode.label(),
            priors.len(),
            prior_iterations,
            exact_tie_count,
            resolved_ties,
            posterior_mean,
            baseline_selected.len(),
            guarded_additions,
            low_abundance_additions,
            selected_taxa.len(),
            premap_guard,
            0.05 * safeguard_min_chunk_breadth,
            post_requirement,
        );
    }

    Ok(CoordinateEvidenceResult {
        decisions: post_decisions,
        selected_taxa,
    })
}

/// Native short-/long-read assignment policy.
///
/// All valid terminal-taxon candidates are retained through evidence
/// resolution; endpoints with independent presence support remain available
/// to abundance EM. Baseline coordinate evidence is deliberately stricter:
/// only a platform-specific, margin-resolved best endpoint contributes
/// intervals. Guarded MAP can later add an exact short-read tie only when its
/// posterior guard passes. Each represented assembly remains an independent
/// witness.
fn process_evidence_batch(
    chunk: &[ReadBundle],
    cfg: &ProfilerConfig,
    policy: EvidencePolicy,
    fallback_mode: bool,
    species_to_genus: &HashMap<u32, u32>,
) -> AssignmentResult {
    let mut result = AssignmentResult::default();
    let mut scratch = CandidateScratch::default();
    let mut indexed_candidates: Vec<usize> = Vec::new();

    // Genus fallback remains a separate, lower-rank path and runs only when
    // no terminal candidate survives the native policy.
    let genus_fb_mapped_min = genus_fallback_min_mapped_len(cfg);
    let mut genus_best: BestMap = HashMap::new();
    let mut genus_order: Vec<u32> = Vec::new();
    let mut genus_intervals: IntervalMap = HashMap::new();

    for read in chunk {
        let candidates = scratch.build(&read.alignments, policy);
        result.filter_stats.evaluated = result
            .filter_stats
            .evaluated
            .saturating_add(read.alignments.len());

        if !candidates.is_empty() {
            let resolved_index = resolved_candidate_index(candidates, policy.mode);
            let profile_boundary_incomplete = candidates
                .iter()
                .any(|candidate| candidate.candidate_boundary.profile_incomplete());
            if profile_boundary_incomplete {
                result.boundary_ambiguous_fragments =
                    result.boundary_ambiguous_fragments.saturating_add(1);
            }
            result.filter_stats.passed =
                result.filter_stats.passed.saturating_add(candidates.len());
            for candidate in candidates {
                let evidence = result.taxon_evidence.entry(candidate.taxid).or_default();
                evidence.candidate_fragments = evidence.candidate_fragments.saturating_add(1);
                if resolved_index.is_some() {
                    evidence.resolved_candidate_fragments =
                        evidence.resolved_candidate_fragments.saturating_add(1);
                }
            }
            if candidates.len() == 1 && !profile_boundary_incomplete {
                result.uniquely_supported_taxa.insert(candidates[0].taxid);
            }

            // Stable taxid order canonicalizes ambiguous EM keys independent
            // of mapping row order and candidate score ties.
            indexed_candidates.clear();
            indexed_candidates.extend(0..candidates.len());
            indexed_candidates.sort_unstable_by_key(|&index| {
                (candidates[index].taxid, candidates[index].profile_rank)
            });
            let mut taxids = Vec::with_capacity(candidates.len());
            let mut depths = Vec::with_capacity(candidates.len());
            let mut mapped_lens = Vec::with_capacity(candidates.len());
            let mut sampling_scales = Vec::with_capacity(candidates.len());
            for &index in &indexed_candidates {
                let candidate = &candidates[index];
                let placement = candidate.best;
                let mapped = match policy.mode {
                    EvidenceMode::Paired | EvidenceMode::Single => placement.r1_mapped,
                    EvidenceMode::Long => placement.both_mapped,
                };
                taxids.push(candidate.taxid);
                depths.push(f64::from(mapped) / placement.asm_len.max(1.0));
                mapped_lens.push(f64::from(mapped));
                sampling_scales.push(placement.sampling_scale);
            }
            result.read_assignments.push(ReadAssignment {
                read: Arc::clone(&read.read_name),
                rank: homogeneous_candidate_rank(candidates),
                requires_terminal_evidence: true,
                taxids,
                depths,
                mapped_lens,
                sampling_scales,
            });

            if policy.collect_exact_ties
                && candidates
                    .iter()
                    .all(|candidate| candidate.profile_rank.is_species())
            {
                result.exact_ties.record(candidates);
            }

            if let Some(index) = resolved_index {
                // Direct unique or margin-resolved winners count as both
                // resolved-best and clear-best support. An unresolved tie
                // counts for neither taxon; guarded MAP may later add only its
                // successfully resolved winner to `best_fragments`.
                let evidence = result
                    .taxon_evidence
                    .entry(candidates[index].taxid)
                    .or_default();
                evidence.best_fragments = evidence.best_fragments.saturating_add(1);
                evidence.clear_best_fragments = evidence.clear_best_fragments.saturating_add(1);
                for placement in &candidates[index].placements {
                    let key = PresenceKey {
                        taxid: placement.taxid,
                        profile_rank: placement.profile_rank,
                        assembly_id: placement.assembly_id,
                    };
                    let presence = result.presence.entry(key).or_default();
                    if presence.asm_len <= 0.0 {
                        presence.asm_len = placement.asm_len;
                    }
                    presence.count = presence.count.saturating_add(1);
                    presence.reference_span_bases = presence.reference_span_bases.saturating_add(
                        u64::from(placement.span.1.saturating_sub(placement.span.0)),
                    );
                    presence.intervals.push(PresenceInterval {
                        reference_id: placement.reference_id,
                        span: placement.span,
                        fixed_chunk: placement.fixed_chunk,
                        flex_chunk: placement.flex_chunk,
                    });
                }
            }
            continue;
        }

        if !fallback_mode {
            continue;
        }
        genus_best.clear();
        genus_order.clear();
        genus_intervals.clear();
        for aln in &read.alignments {
            // Paired fallback uses only a validated first mate from a complete
            // concordant group. Long reads use their sole unpaired row.
            if !matches!(aln.role, MappingRole::Mate1 | MappingRole::Unpaired) {
                continue;
            }
            if !aln.profile_rank.is_species() {
                continue;
            }
            let identity = f64::from(aln.identity);
            let fraction = f64::from(aln.fraction);
            let length_ok = aln.mapped_len as usize >= genus_fb_mapped_min;
            let fraction_ok = fraction >= GENUS_FALLBACK_MIN_FRACTION;
            if identity < cfg.genus_identity
                || !coverage_passes(length_ok, fraction_ok, cfg.coverage_gate_is_or())
            {
                continue;
            }
            let Some(genus_taxid) = species_to_genus.get(&aln.profile_taxid).copied() else {
                continue;
            };
            update_best(
                &mut genus_best,
                &mut genus_order,
                &mut genus_intervals,
                genus_taxid,
                identity,
                aln.mapped_len as f64 / aln.asm_len,
                aln,
            );
        }
        if !genus_order.is_empty() {
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

fn homogeneous_candidate_rank(candidates: &[super::evidence::TerminalCandidate]) -> AssignRank {
    let Some(first) = candidates.first() else {
        return AssignRank::Mixed;
    };
    if candidates
        .iter()
        .all(|candidate| candidate.profile_rank == first.profile_rank)
    {
        AssignRank::from_profile_rank(first.profile_rank)
    } else {
        AssignRank::Mixed
    }
}

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
            let identity = f64::from(aln.identity);
            let fraction = f64::from(aln.fraction);

            // --- Step 1: species criteria (unchanged thresholds) ---
            let length_ok = aln.mapped_len as usize >= cfg.mapped_len;
            let fraction_ok = fraction >= cfg.fraction;
            let species_qualified = identity >= cfg.identity
                && coverage_passes(length_ok, fraction_ok, cfg.coverage_gate_is_or());

            if species_qualified {
                update_best(
                    &mut species_best,
                    &mut species_order,
                    &mut species_intervals,
                    aln.profile_taxid,
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
            let fb_length_ok = aln.mapped_len as usize >= genus_fb_mapped_min;
            let fb_fraction_ok = fraction >= GENUS_FALLBACK_MIN_FRACTION;
            let fallback_qualified = identity >= cfg.genus_identity
                && coverage_passes(fb_length_ok, fb_fraction_ok, cfg.coverage_gate_is_or());
            if !fallback_qualified {
                continue;
            }
            let genus_taxid = match species_to_genus.get(&aln.profile_taxid).copied() {
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
            entry.insert((
                identity,
                depth,
                aln.mapped_len as f64,
                f64::from(aln.sampling_scale),
            ));
            let acc_key = AccKey {
                taxid,
                asm: aln.asm.clone(),
                acc: aln.acc.clone(),
            };
            intervals.insert(
                taxid,
                (
                    acc_key,
                    widen_pair(aln.span),
                    widen_pair(aln.fixed_chunk),
                    widen_pair(aln.flex_chunk),
                ),
            );
        }
        Entry::Occupied(mut entry) => {
            if identity > entry.get().0 {
                *entry.get_mut() = (
                    identity,
                    depth,
                    aln.mapped_len as f64,
                    f64::from(aln.sampling_scale),
                );
                let acc_key = AccKey {
                    taxid,
                    asm: aln.asm.clone(),
                    acc: aln.acc.clone(),
                };
                intervals.insert(
                    taxid,
                    (
                        acc_key,
                        widen_pair(aln.span),
                        widen_pair(aln.fixed_chunk),
                        widen_pair(aln.flex_chunk),
                    ),
                );
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
/// the same `AccKey`, the key is cloned into the first two maps and moved
/// into the third.
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
    let mut sampling_scales = Vec::with_capacity(order.len());
    for &taxid in order {
        if let Some((_, depth, maplen, sampling_scale)) = best.get(&taxid) {
            taxids.push(taxid);
            depths.push(*depth);
            maplens.push(*maplen);
            sampling_scales.push(*sampling_scale);
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
        read: Arc::<str>::from(read_name),
        rank,
        requires_terminal_evidence: false,
        taxids,
        depths,
        mapped_lens: maplens,
        sampling_scales,
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
///   parse time because their reference name did not carry a parsable
///   species taxid or assembly length
/// * `next_progress` - Next milestone for progress logging
fn absorb_worker_result(
    result: WorkerResult,
    read_map: &mut StableReadMap<AlignmentEntry>,
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
                target: "PROFILER",
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
        None => return Some(discard(true, false)),
    };
    if ref_len == 0 {
        return Some(discard(true, false));
    }

    let mapped_len = read_seq.len();
    if mapped_len == 0 {
        return Some(discard(false, false));
    }

    // Parse the reference-name taxonomy fields once up front. Parse failures
    // are reported through `ParsedAlignment::skipped_reference`.
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
    let sampling_scale = ctx.asm_sampling_scales.get(asm_str).copied().unwrap_or(1.0) as f32;

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
        asm: Arc::<str>::from(asm_str),
        acc: Arc::<str>::from(acc_str),
        profile_taxid: species_taxid,
        profile_rank: ProfileRank::Species,
        reference_id: u32::MAX,
        assembly_id: u32::MAX,
        role: MappingRole::Unpaired,
        pair_id: 0,
        candidate_boundary: metam::mapper::CandidateBoundary::COMPLETE,
        sampling_scale,
        asm_len,
        identity: identity as f32,
        fraction: fraction as f32,
        mapped_len: clamp_u32(mapped_len),
        span: (clamp_u32(start), clamp_u32(end)),
        fixed_chunk: (clamp_u32(fixed_chunk_start), clamp_u32(fixed_chunk_end)),
        flex_chunk: (clamp_u32(flex_chunk_start), clamp_u32(flex_chunk_end)),
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
    if cfg.sequencer.is_illumina() {
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

/// Run-level decision for the read count fed into the automatic chunk-breadth
/// curve.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ChunkBreadthPlan {
    /// `true` when the unscaled basis is mapped reads; `false` for total reads.
    use_aligned_basis: bool,
    /// `true` when aligned basis was selected by the automatic low-map rule.
    by_aligned_auto: bool,
    /// `true` when subsampled fallback uses total reads for threshold estimation.
    force_total_for_subsampled_fallback: bool,
    /// Unscaled read count selected as the chunk-breadth basis.
    basis_reads: usize,
    /// Read count actually passed to the automatic chunk-breadth curve.
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

/// Recall-mode scaling starts only once the automatic cutoff is based on more
/// than one million effective reads/fragments. This leaves shallow samples on
/// the default curve, where lowering an already-small amount gate
/// mostly adds noise.
const RECALL_CHUNK_SCALE_MIN_BASIS: usize = 1_000_000;
const SHORT_RECALL_CHUNK_SCALE_FLOOR: f64 = 0.70;
const LONG_RECALL_CHUNK_SCALE: f64 = 0.80;

/// Standard and safeguard cutoffs are deliberately distinct. A user scale
/// first adjusts the shared automatic safeguard; recall mode may then lower
/// only the ordinary standard-evidence gate.
#[derive(Clone, Copy, Debug, PartialEq)]
struct CommunityChunkBreadthThresholds {
    /// Calibrated read-count curve before any user or mode adjustment.
    automatic: f64,
    /// Automatic cutoff after the explicit user scale. Rescue and guarded-MAP
    /// safeguards use this value.
    standard: f64,
    safeguard: f64,
    user_scale: f64,
    recall_scale: f64,
}

impl CommunityChunkBreadthThresholds {
    #[inline]
    fn is_user_scaled(self) -> bool {
        self.user_scale != 1.0 && self.automatic > 0.0
    }

    #[inline]
    fn is_recall_scaled(self) -> bool {
        self.recall_scale < 1.0 && self.safeguard > 0.0
    }

    #[inline]
    fn is_scaled(self) -> bool {
        self.is_user_scaled() || self.is_recall_scaled()
    }
}

/// Resolve community chunk-breadth thresholds once for both native mappings
/// and the compatibility SAM path.
fn resolve_community_chunk_breadth(
    cfg: &ProfilerConfig,
    basis_reads: usize,
) -> CommunityChunkBreadthThresholds {
    // Explicit user policy is authoritative in every mode, including
    // low-biomass mode and recall mode.
    if let Some(cutoff) = cfg.chunk_breadth {
        return CommunityChunkBreadthThresholds {
            automatic: cutoff,
            standard: cutoff,
            safeguard: cutoff,
            user_scale: 1.0,
            recall_scale: 1.0,
        };
    }

    // Preserve the established --min-reads and low-biomass behavior: no
    // automatic chunk-breadth filtering, and therefore nothing to scale.
    if cfg.min_reads.is_some() || cfg.lowbiomass {
        return CommunityChunkBreadthThresholds {
            automatic: 0.0,
            standard: 0.0,
            safeguard: 0.0,
            user_scale: 1.0,
            recall_scale: 1.0,
        };
    }

    let automatic = community_chunk_breadth_cutoff(basis_reads, cfg.sequencer.is_illumina());
    let safeguard = scale_automatic_chunk_breadth(automatic, cfg.chunk_breadth_scale);
    let recall_scale = automatic_recall_chunk_breadth_scale(cfg, basis_reads);
    let standard = if recall_scale == 1.0 {
        safeguard
    } else {
        safeguard * recall_scale
    };
    CommunityChunkBreadthThresholds {
        automatic,
        standard,
        safeguard,
        user_scale: cfg.chunk_breadth_scale,
        recall_scale,
    }
}

/// Scale only an automatic community cutoff in recall mode. The short-read
/// curve is continuous at one million reads and approaches 0.70 smoothly:
/// `scale = 0.70 + 0.30 / m`, where `m = basis_reads / 1e6`. Long reads use
/// a fixed 0.80 multiplier above the same depth gate.
#[inline]
fn automatic_recall_chunk_breadth_scale(cfg: &ProfilerConfig, basis_reads: usize) -> f64 {
    if cfg.mode != ProfileMode::Recall
        || cfg.strain
        || cfg.lowbiomass
        || cfg.chunk_breadth.is_some()
        || cfg.min_reads.is_some()
        || basis_reads <= RECALL_CHUNK_SCALE_MIN_BASIS
    {
        return 1.0;
    }

    if cfg.sequencer.is_illumina() {
        let millions = basis_reads as f64 / 1_000_000.0;
        (SHORT_RECALL_CHUNK_SCALE_FLOOR + (1.0 - SHORT_RECALL_CHUNK_SCALE_FLOOR) / millions)
            .clamp(SHORT_RECALL_CHUNK_SCALE_FLOOR, 1.0)
    } else if cfg.sequencer.is_long_read() {
        LONG_RECALL_CHUNK_SCALE
    } else {
        1.0
    }
}

fn community_chunk_breadth_source(
    cfg: &ProfilerConfig,
    plan: ChunkBreadthPlan,
    thresholds: CommunityChunkBreadthThresholds,
) -> &'static str {
    if cfg.chunk_breadth.is_some() {
        "user-supplied --min-chunk-breadth"
    } else if cfg.min_reads.is_some() {
        "disabled by --min-reads"
    } else if cfg.lowbiomass {
        "disabled by --low-biomass"
    } else {
        match (thresholds.is_scaled(), plan.scaled, plan.use_aligned_basis) {
            (true, true, _) => "recall-scaled auto-estimate from scaled mapped reads",
            (true, false, true) => "recall-scaled auto-estimate from mapped reads",
            (true, false, false) => "recall-scaled auto-estimate from total reads",
            (false, true, _) => "auto-estimated from scaled mapped reads",
            (false, false, true) => "auto-estimated from mapped reads",
            (false, false, false) => "auto-estimated from total reads",
        }
    }
}

fn log_chunk_breadth_scaling(thresholds: CommunityChunkBreadthThresholds) {
    if thresholds.is_user_scaled() {
        let uncapped = thresholds.automatic * thresholds.user_scale;
        if uncapped > 1.0 {
            log::info!(
                target: "PROFILER",
                "User chunk-breadth scale adjusted the automatic safeguard: {:.6} x {} = {:.6}, capped at 1.000000. Rescue and guarded-MAP safeguards use the adjusted value.",
                thresholds.automatic,
                thresholds.user_scale,
                thresholds.safeguard,
            );
        } else {
            log::info!(
                target: "PROFILER",
                "User chunk-breadth scale adjusted the automatic safeguard: {:.6} x {} = {:.6}. Rescue and guarded-MAP safeguards use the adjusted value.",
                thresholds.automatic,
                thresholds.user_scale,
                thresholds.safeguard,
            );
        }
    }
    if thresholds.is_recall_scaled() {
        log::info!(
            target: "PROFILER",
            "Recall mode scaled only the automatic standard chunk-breadth threshold: {:.6} x {:.4} = {:.6}; rescue and guarded-MAP safeguards retain {:.6}, and OEBR/COEBR thresholds are unchanged.",
            thresholds.safeguard,
            thresholds.recall_scale,
            thresholds.standard,
            thresholds.safeguard,
        );
    }
}

/// Scale the read count used as input to the automatic chunk-breadth curve.
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

/// Depth-scaled minimum chunk-breadth cutoff calibrated for community
/// profiling.
///
/// Both platforms use a bounded, saturating Michaelis-Menten form in
/// `m = num_reads / 1e6`:
///
/// - Illumina: `0.66*m/(m+1) + 0.24*m/(m+100)`.
/// - Nanopore/PacBio: `0.85*m/(m+0.78)`.
///
/// The cutoff rises smoothly with sequencing depth and then plateaus. The
/// 0.95 clamp is defensive; both fitted plateaus are lower.
fn community_chunk_breadth_cutoff(num_reads: usize, illumina: bool) -> f64 {
    let m = num_reads as f64 / 1_000_000.0;
    let value = if illumina {
        0.66 * m / (m + 1.0) + 0.24 * m / (m + 100.0)
    } else {
        0.85 * m / (m + 0.78)
    };
    value.min(0.95)
}

/// Historical default lower ratio bound for non-LCA output rows.
fn rank_filter_lower_bound(rank: &str) -> f64 {
    if rank == "genus" {
        0.65
    } else {
        0.75
    }
}

/// Resolve the standard evidence-tier ratio defaults for a majority-LCA rank.
/// Unsupported taxonomy ranks and non-read modes retain the established
/// compatibility filter rather than guessing a new policy.
#[inline]
fn standard_lca_ratio_lower_bounds(rank: &str, mode: Option<EvidenceMode>) -> Option<(f64, f64)> {
    Some(standard_ratio_lower_bounds(
        ProfileRank::from_taxonomy_rank(rank)?,
        mode?,
    ))
}

/// Established ratio filter for non-LCA output rows.
fn ratio_filter_passes(
    rank: &str,
    breadth: f64,
    expected_breadth: f64,
    flex_chunk: f64,
    expected_flex: f64,
    min_oebr: Option<f64>,
    min_coebr: Option<f64>,
) -> bool {
    let default_lower = rank_filter_lower_bound(rank);
    ratio_filter_passes_with_defaults(
        breadth,
        expected_breadth,
        flex_chunk,
        expected_flex,
        (default_lower, default_lower),
        min_oebr,
        min_coebr,
    )
}

/// Apply OEBR/COEBR filtering with independently selected default lower
/// bounds. Explicit CLI overrides remain authoritative. The inclusive 1.50
/// upper bound is shared by the established output path and standard evidence.
fn ratio_filter_passes_with_defaults(
    breadth: f64,
    expected_breadth: f64,
    flex_chunk: f64,
    expected_flex: f64,
    default_lowers: (f64, f64),
    min_oebr: Option<f64>,
    min_coebr: Option<f64>,
) -> bool {
    if expected_breadth <= 0.0 || expected_flex <= 0.0 {
        return false;
    }
    let min_oebr = min_oebr.unwrap_or(default_lowers.0);
    let min_coebr = min_coebr.unwrap_or(default_lowers.1);
    let br_ratio = breadth / expected_breadth;
    let fl_ratio = flex_chunk / expected_flex;
    br_ratio.is_finite()
        && fl_ratio.is_finite()
        && min_oebr.is_finite()
        && min_coebr.is_finite()
        && br_ratio >= min_oebr
        && br_ratio <= STANDARD_RATIO_UPPER
        && fl_ratio >= min_coebr
        && fl_ratio <= STANDARD_RATIO_UPPER
}

/// Estimate abundance from every mapping-quality-qualified assignment without
/// applying taxon-level evidence/report filters.
///
/// This mirrors the established majority-LCA and EM mechanics used by the
/// production branch, but intentionally has no classification side effects and
/// never changes the production assignment state. It is invoked only for
/// `--save-raw-profile`.
fn estimate_unfiltered_taxon_totals(
    assignments: &[ReadAssignment],
    taxonomy: &Taxonomy,
    coverage: &mut HashMap<u32, TaxonCoverage>,
    max_iterations: usize,
    verbose_logging: bool,
) -> Result<HashMap<u32, TaxonTotals>> {
    let mut unambiguous: HashMap<u32, TaxonTotals> = HashMap::new();
    let mut ambiguous: HashMap<Vec<u32>, AmbiguousTotals> = HashMap::new();

    for assignment in assignments {
        let candidate_count = assignment.taxids.len();
        if assignment.depths.len() != candidate_count
            || assignment.mapped_lens.len() != candidate_count
            || assignment.sampling_scales.len() != candidate_count
        {
            bail!(
                "inconsistent raw-profile candidate vectors for read '{}': taxids={}, depths={}, mapped_lens={}, sampling_scales={}",
                assignment.read,
                candidate_count,
                assignment.depths.len(),
                assignment.mapped_lens.len(),
                assignment.sampling_scales.len(),
            );
        }
        if candidate_count == 0 {
            continue;
        }
        for index in 0..candidate_count {
            let depth = assignment.depths[index];
            let mapped_len = assignment.mapped_lens[index];
            let sampling_scale = assignment.sampling_scales[index];
            if !depth.is_finite()
                || depth < 0.0
                || !mapped_len.is_finite()
                || mapped_len < 0.0
                || !sampling_scale.is_finite()
                || sampling_scale <= 0.0
            {
                bail!(
                    "invalid raw-profile candidate values for read '{}' taxid {}: depth={}, mapped_len={}, sampling_scale={}",
                    assignment.read,
                    assignment.taxids[index],
                    depth,
                    mapped_len,
                    sampling_scale,
                );
            }
        }

        if candidate_count == 1 {
            let entry = unambiguous.entry(assignment.taxids[0]).or_default();
            entry
                .support
                .add_weighted(1.0, assignment.sampling_scales[0], 1.0);
            entry.depth += assignment.depths[0];
            entry.mapped_len += assignment.mapped_lens[0];
            continue;
        }

        if let Some(entry) = ambiguous.get_mut(assignment.taxids.as_slice()) {
            entry.observed_count += 1.0;
            accumulate(&mut entry.depths, &assignment.depths);
            accumulate(&mut entry.mapped_lens, &assignment.mapped_lens);
            accumulate(
                &mut entry.equivalent_support_sums,
                &assignment.sampling_scales,
            );
        } else {
            ambiguous.insert(
                assignment.taxids.clone(),
                AmbiguousTotals {
                    observed_count: 1.0,
                    depths: assignment.depths.clone(),
                    mapped_lens: assignment.mapped_lens.clone(),
                    equivalent_support_sums: assignment.sampling_scales.clone(),
                },
            );
        }
    }

    // Remove candidates with no independent support exactly as the production
    // EM does. If an entire candidate set lacks such support, fold it to the
    // established 70%-majority LCA instead of inventing a species prior.
    let mut em_groups: HashMap<Vec<u32>, AmbiguousTotals> = HashMap::new();
    let mut ambiguous = ambiguous.into_iter().collect::<Vec<_>>();
    ambiguous.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    for (taxids, totals) in ambiguous {
        let independently_supported = taxids
            .iter()
            .map(|taxid| unambiguous.contains_key(taxid))
            .collect::<Vec<_>>();
        if independently_supported.iter().all(|&present| present) {
            merge_ambiguous_totals(&mut em_groups, taxids, totals);
            continue;
        }
        if independently_supported.iter().any(|&present| present) {
            let mut filtered_taxids = Vec::with_capacity(taxids.len());
            let mut filtered_depths = Vec::with_capacity(taxids.len());
            let mut filtered_maplens = Vec::with_capacity(taxids.len());
            let mut filtered_scale_sums = Vec::with_capacity(taxids.len());
            for (index, taxid) in taxids.iter().enumerate() {
                if independently_supported[index] {
                    filtered_taxids.push(*taxid);
                    filtered_depths.push(totals.depths[index]);
                    filtered_maplens.push(totals.mapped_lens[index]);
                    filtered_scale_sums.push(totals.equivalent_support_sums[index]);
                }
            }
            if filtered_taxids.len() == 1 {
                let entry = unambiguous.entry(filtered_taxids[0]).or_default();
                entry
                    .support
                    .add_weighted(totals.observed_count, filtered_scale_sums[0], 1.0);
                entry.depth += filtered_depths[0];
                entry.mapped_len += filtered_maplens[0];
            } else if !filtered_taxids.is_empty() {
                merge_ambiguous_totals(
                    &mut em_groups,
                    filtered_taxids,
                    AmbiguousTotals {
                        observed_count: totals.observed_count,
                        depths: filtered_depths,
                        mapped_lens: filtered_maplens,
                        equivalent_support_sums: filtered_scale_sums,
                    },
                );
            }
            continue;
        }

        let taxid_set = taxids.iter().copied().collect::<HashSet<_>>();
        let Some((lca_taxid, _lca_name, _lca_rank)) = taxonomy.get_majority_lca(&taxid_set, 0.7)
        else {
            continue;
        };
        let lca_breadth = taxids
            .iter()
            .filter_map(|taxid| coverage.get(taxid))
            .map(|value| value.breadth)
            .fold(0.0, f64::max);
        let lca_fixed = taxids
            .iter()
            .filter_map(|taxid| coverage.get(taxid))
            .map(|value| value.fixed_chunk)
            .fold(0.0, f64::max);
        let lca_flex = taxids
            .iter()
            .filter_map(|taxid| coverage.get(taxid))
            .map(|value| value.flex_chunk)
            .fold(0.0, f64::max);
        coverage.insert(
            lca_taxid,
            TaxonCoverage {
                breadth: lca_breadth,
                fixed_chunk: lca_fixed,
                flex_chunk: lca_flex,
            },
        );
        let candidate_count = totals.depths.len() as f64;
        let depth_mean = totals.depths.iter().sum::<f64>() / candidate_count;
        let maplen_mean = totals.mapped_lens.iter().sum::<f64>() / candidate_count;
        let equivalent_mean = totals.equivalent_support_sums.iter().sum::<f64>() / candidate_count;
        let entry = unambiguous.entry(lca_taxid).or_default();
        entry
            .support
            .add_weighted(totals.observed_count, equivalent_mean, 1.0);
        entry.depth += depth_mean;
        entry.mapped_len += maplen_mean;
    }

    let base = unambiguous;
    let mut current = base.clone();
    let mut em_groups = em_groups.iter().collect::<Vec<_>>();
    em_groups.sort_unstable_by(|left, right| left.0.cmp(right.0));
    let mut em_taxids = base.keys().copied().collect::<Vec<_>>();
    for (taxids, _) in &em_groups {
        em_taxids.extend(taxids.iter().copied());
    }
    em_taxids.sort_unstable();
    em_taxids.dedup();
    let max_iterations = max_iterations.max(1);
    for iteration in 0..max_iterations {
        let mut next = base.clone();
        for &(taxids, totals) in &em_groups {
            let sum_depth = taxids
                .iter()
                .map(|taxid| current.get(taxid).map_or(0.0, |value| value.depth))
                .sum::<f64>();
            if sum_depth <= f64::EPSILON || !sum_depth.is_finite() {
                continue;
            }
            let mut weighted = taxids
                .iter()
                .enumerate()
                .map(|(index, &taxid)| {
                    let fraction = current.get(&taxid).map_or(0.0, |value| value.depth) / sum_depth;
                    (taxid, index, fraction)
                })
                .collect::<Vec<_>>();
            weighted.sort_by(|left, right| {
                right
                    .2
                    .partial_cmp(&left.2)
                    .unwrap_or(Ordering::Equal)
                    .then_with(|| left.0.cmp(&right.0))
            });
            for (taxid, index, fraction) in weighted {
                let entry = next.entry(taxid).or_default();
                if fraction >= 0.99 {
                    entry.support.add_weighted(
                        totals.observed_count,
                        totals.equivalent_support_sums[index],
                        1.0,
                    );
                    entry.depth += totals.depths[index];
                    entry.mapped_len += totals.mapped_lens[index];
                    break;
                }
                entry.support.add_weighted(
                    totals.observed_count,
                    totals.equivalent_support_sums[index],
                    fraction,
                );
                entry.depth += totals.depths[index] * fraction;
                entry.mapped_len += totals.mapped_lens[index] * fraction;
            }
        }

        if verbose_logging && iteration % 5 == 0 {
            log::info!(
                target: "PROFILER",
                "Raw-profile EM iteration {} ...",
                iteration
            );
        }
        let numerator = em_taxids
            .iter()
            .map(|taxid| {
                let depth = next.get(taxid).map_or(0.0, |totals| totals.depth);
                let previous = current.get(taxid).map_or(0.0, |value| value.depth);
                (depth - previous).abs()
            })
            .sum::<f64>();
        let denominator = em_taxids
            .iter()
            .map(|taxid| current.get(taxid).map_or(0.0, |totals| totals.depth))
            .sum::<f64>()
            .max(1e-12);
        current = next;
        if numerator / denominator <= 1e-5 {
            if verbose_logging {
                log::info!(
                    target: "PROFILER",
                    "Finished raw-profile EM with {} iterations",
                    iteration
                );
            }
            return Ok(current);
        }
    }

    log::warn!(
        target: "PROFILER",
        "Raw-profile EM reached the maximum of {} iterations before convergence.",
        max_iterations
    );
    Ok(current)
}

fn merge_ambiguous_totals(
    groups: &mut HashMap<Vec<u32>, AmbiguousTotals>,
    taxids: Vec<u32>,
    totals: AmbiguousTotals,
) {
    match groups.entry(taxids) {
        Entry::Vacant(entry) => {
            entry.insert(totals);
        }
        Entry::Occupied(mut entry) => {
            let current = entry.get_mut();
            current.observed_count += totals.observed_count;
            accumulate(&mut current.depths, &totals.depths);
            accumulate(&mut current.mapped_lens, &totals.mapped_lens);
            accumulate(
                &mut current.equivalent_support_sums,
                &totals.equivalent_support_sums,
            );
        }
    }
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
        // Legacy coverage probabilities are only computed for species rows, so
        // for every higher rank their absence is expected rather than evidence
        // of anything, and strict mode must not reject a test that was never
        // run. A species row missing them is a different matter and still
        // fails.
        (None, None) if strict && rank != "species" => true,
        _ => !strict,
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

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
    let mut taxids = map.keys().copied().collect::<Vec<_>>();
    taxids.sort_unstable();
    let sum: f64 = taxids.iter().map(|taxid| map[taxid]).sum::<f64>();
    if sum <= f64::EPSILON {
        return;
    }
    for value in map.values_mut() {
        *value = round5(*value / sum);
    }
}

fn community_profile_row_passes(
    row: &ProfileRow,
    cfg: &ProfilerConfig,
    min_chunk_breadth: f64,
    evidence_selected_taxa: &HashSet<u32>,
    lca_ratio_lower_bounds: &HashMap<u32, (f64, f64)>,
) -> bool {
    let min_reads_ok = equivalent_read_count_passes_minimum(row.equivalent_mapcount, cfg.min_reads);
    if evidence_selected_taxa.contains(&row.taxid) {
        return min_reads_ok;
    }
    let fixed_chunk_ok = if cfg.min_reads.is_some() && cfg.chunk_breadth.is_none() {
        true
    } else {
        row.fixed_chunk > min_chunk_breadth
    };
    let ratio_ok = if let Some(default_lowers) = lca_ratio_lower_bounds.get(&row.taxid).copied() {
        ratio_filter_passes_with_defaults(
            row.breadth,
            row.expected_breadth,
            row.flex_chunk,
            row.expected_flex_chunk,
            default_lowers,
            cfg.min_oebr,
            cfg.min_coebr,
        )
    } else {
        ratio_filter_passes(
            &row.rank,
            row.breadth,
            row.expected_breadth,
            row.flex_chunk,
            row.expected_flex_chunk,
            cfg.min_oebr,
            cfg.min_coebr,
        )
    };
    min_reads_ok
        && fixed_chunk_ok
        && ratio_ok
        && pvalue_filter_passes(&row.rank, row.breadth_pvalue, row.chunk_pvalue, cfg.strict)
}

fn write_community_classifications(
    path: Option<&str>,
    resolved: &ResolvedClassifications,
    ambiguous_counts: &HashMap<Vec<u32>, AmbiguousTotals>,
    ambiguous_reads: &HashMap<Vec<u32>, Vec<String>>,
    ambiguous_fraction: &HashMap<Vec<u32>, HashMap<u32, f64>>,
    retained_taxa: &HashSet<u32>,
    terminal_ranks: &HashMap<u32, ProfileRank>,
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

    // Classification rows are optional and potentially numerous. Sort only
    // when they are requested, and sort borrowed pointers rather than cloning
    // read names. The normal count-only path above remains allocation-free.
    if let Some(writer) = writer.as_mut() {
        for group in resolved.sorted_groups() {
            if !retained_taxa.contains(&group.key.taxid) {
                continue;
            }
            let rank = terminal_ranks
                .get(&group.key.taxid)
                .copied()
                .map(ProfileRank::as_str)
                .or_else(|| taxonomy.get_rank(group.key.taxid))
                .unwrap_or("no rank");
            let sources = match &group.key.source {
                ClassificationSource::Direct => group.key.taxid.to_string(),
                ClassificationSource::Candidates(taxids) => taxids
                    .iter()
                    .map(u32::to_string)
                    .collect::<Vec<_>>()
                    .join(";"),
            };
            for read in group.sorted_reads() {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t1",
                    read, group.key.taxid, rank, sources
                )?;
            }
        }
    }

    let mut keys = ambiguous_counts.keys().collect::<Vec<_>>();
    keys.sort_unstable();
    for taxids in keys {
        let Some(fractions) = ambiguous_fraction.get(taxids) else {
            continue;
        };
        let mut retained = fractions
            .iter()
            .filter_map(|(&taxid, &fraction)| {
                (retained_taxa.contains(&taxid) && fraction.is_finite() && fraction > 0.0)
                    .then_some((taxid, fraction))
            })
            .collect::<Vec<_>>();
        retained.sort_unstable_by_key(|(taxid, _)| *taxid);
        let total = retained.iter().map(|(_, fraction)| *fraction).sum::<f64>();
        if total <= f64::EPSILON || !total.is_finite() {
            continue;
        }
        for (_, fraction) in &mut retained {
            *fraction /= total;
        }
        retained.sort_by(|left, right| {
            right
                .1
                .total_cmp(&left.1)
                .then_with(|| left.0.cmp(&right.0))
        });
        let count = ambiguous_counts[taxids].observed_count.round() as usize;
        classified = classified.saturating_add(count);

        let Some(writer) = writer.as_mut() else {
            continue;
        };
        let final_taxid = retained[0].0;
        let rank = terminal_ranks
            .get(&final_taxid)
            .copied()
            .map(ProfileRank::as_str)
            .or_else(|| taxonomy.get_rank(final_taxid))
            .unwrap_or("no rank");
        let candidate_taxids = retained
            .iter()
            .map(|(taxid, _)| taxid.to_string())
            .collect::<Vec<_>>()
            .join(";");
        let candidate_fractions = retained
            .iter()
            .map(|(_, fraction)| format!("{fraction:.5}"))
            .collect::<Vec<_>>()
            .join(";");
        if let Some(reads) = ambiguous_reads.get(taxids) {
            let mut reads = reads.iter().map(String::as_str).collect::<Vec<_>>();
            reads.sort_unstable();
            for read in reads {
                writeln!(
                    writer,
                    "{read}\t{final_taxid}\t{rank}\t{candidate_taxids}\t{candidate_fractions}"
                )?;
            }
        }
    }

    if let Some(writer) = writer.as_mut() {
        writer.flush()?;
    }
    Ok(classified)
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

    /// Effective standard chunk-breadth threshold. In deep recall-mode
    /// community runs this may be lower than the original automatic cutoff.
    min_chunk_breadth: f64,

    /// Original unscaled cutoff retained by rescue and guarded-MAP safeguards.
    safeguard_min_chunk_breadth: f64,

    /// Whether at least one genome has a full/effective size multiplier > 1.
    has_scaled_genomes: bool,

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

    /// Mapping from species taxid to its enclosing genus taxid, derived from
    /// the taxonomy lineage (not parsed from reference headers). Keyed by
    /// the **post-normalization** species taxid so that sub-species
    /// aliases all map through the same entry.
    /// Species without a resolvable genus are simply absent from this map.
    species_to_genus: Arc<HashMap<u32, u32>>,

    /// Most-specific reference-supported rank for every native terminal taxid.
    /// Legacy SAM/BAM inputs leave this empty and use taxonomy lookup instead.
    terminal_ranks: Arc<HashMap<u32, ProfileRank>>,

    /// Alignments that were discarded during SAM/BAM parsing because
    /// their reference name did not carry a parsable species taxid /
    /// assembly length. Seeds [`FilterStats::skipped_reference`] so
    /// the diagnostic counter preserves its historical semantics even
    /// now that the check happens in the worker instead of
    /// `process_batch`.
    skipped_reference_at_load: usize,

    /// Native coordinate-evidence policy. `None` disables it and preserves
    /// the legacy behavior.
    evidence_policy: Option<EvidencePolicy>,
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
        record.pos() + 1
    } else {
        0
    };
    let mapq = record.mapq();
    let cigar = if record.cigar().is_empty() {
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
        record.mpos() + 1
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
        target: "PROFILER",
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profiler::evidence::EvidenceCoverage;
    use crate::profiler::{ExactTiePolicy, ProfileMode};
    use crate::sequencer::Sequencer;

    fn approx_equal(lhs: f64, rhs: f64) {
        assert!(
            (lhs - rhs).abs() < 1e-9,
            "left {} not approximately equal to right {}",
            lhs,
            rhs
        );
    }

    #[test]
    fn expected_breadth_is_clamped_for_high_depths() {
        let expected = crate::profiler::evidence::expected_breadth(10.0, 5.0);
        approx_equal(expected, 1.0);
        assert!(expected.is_finite());

        let zero = crate::profiler::evidence::expected_breadth(10.0, 0.0);
        approx_equal(zero, 0.0);
    }

    // ------------------------------------------------------------------
    // Helpers for fallback / rank-aware tests.
    // ------------------------------------------------------------------

    /// Minimal `ProfilerConfig` for unit tests. Callers mutate only the
    /// fields relevant to the scenario under test.
    fn make_cfg(sequencer: Sequencer, host: Option<&str>) -> ProfilerConfig {
        ProfilerConfig {
            sam: String::new(),
            mapper_mappings: None,
            mapping_profile_memory: None,
            db_metadata: None,
            sequencer,
            batch_size: 1,
            max_em_iterations: 200,
            lenient_mapping_input: false,
            is_paired: false,
            exact_ties: ExactTiePolicy::Ambiguous,
            identity: 0.95,
            mapped_len: 50,
            breadth: Some(0.0),
            chunk_breadth: Some(0.0),
            chunk_breadth_scale: 1.0,
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            fraction: 0.6,
            lowbiomass: false,
            keep_raw: false,
            save_classifications: false,
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
            progress: crate::progress::ProgressOptions::default(),
            suppress_mapping_load_progress: false,
            outprefix: "test".to_string(),
            mode: ProfileMode::Default,
            dmp_dir: None,
            strain: false,
        }
    }

    fn test_taxonomy(directory: &Path) -> Arc<Taxonomy> {
        let nodes = directory.join("nodes.dmp");
        let names = directory.join("names.dmp");
        let merged = directory.join("merged.dmp");
        std::fs::write(
            &nodes,
            "1\t|\t1\t|\tno rank\t|\n101\t|\t1\t|\tspecies\t|\n102\t|\t1\t|\tspecies\t|\n",
        )
        .unwrap();
        std::fs::write(
            &names,
            "1\t|\troot\t|\t\t|\tscientific name\t|\n101\t|\tSpecies one\t|\t\t|\tscientific name\t|\n102\t|\tSpecies two\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();
        std::fs::write(&merged, "").unwrap();
        Arc::new(Taxonomy::from_files(&nodes, &names, &merged, None).unwrap())
    }

    #[test]
    fn headerless_sam_is_rejected_before_alignment_workers_start() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let sam = directory.path().join("headerless.sam");
        std::fs::write(
            &sam,
            "read-1\t0\tref-a\t1\t60\t4=\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();
        let profiler = CommunityProfiler::new(make_cfg(Sequencer::Illumina, None), taxonomy)
            .expect("test profiler should construct");
        let error = match profiler.load_text_alignments(&sam) {
            Ok(_) => panic!("headerless SAM input must fail"),
            Err(error) => error,
        };
        assert!(format!("{error:#}").contains("no usable reference sequence dictionary"));
    }

    #[test]
    fn sam_mapping_to_undeclared_reference_is_rejected() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let sam = directory.path().join("unknown-reference.sam");
        std::fs::write(
            &sam,
            "@SQ\tSN:asm-a|0|101|contig-a|1000\tLN:1000\n\
             read-1\t0\tasm-b|0|101|contig-b|1000\t1\t60\t4=\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();
        let profiler = CommunityProfiler::new(make_cfg(Sequencer::Illumina, None), taxonomy)
            .expect("test profiler should construct");
        let error = match profiler.load_text_alignments(&sam) {
            Ok(_) => panic!("mapped RNAME absent from @SQ must fail"),
            Err(error) => error,
        };
        let message = format!("{error:#}");
        assert!(message.contains("1 mapped record(s)"));
        assert!(message.contains("absent from @SQ"));
    }

    fn raw_output_fixture() -> (SamData, AssignmentResult) {
        let good_key = AccKey {
            taxid: 101,
            asm: Arc::<str>::from("asm101"),
            acc: Arc::<str>::from("ref101"),
        };
        let rejected_key = AccKey {
            taxid: 102,
            asm: Arc::<str>::from("asm102"),
            acc: Arc::<str>::from("ref102"),
        };
        let assignment = AssignmentResult {
            read_assignments: vec![
                ReadAssignment {
                    read: "accepted-read".into(),
                    rank: AssignRank::Species,
                    requires_terminal_evidence: false,
                    taxids: vec![101],
                    depths: vec![0.1],
                    mapped_lens: vec![100.0],
                    sampling_scales: vec![1.0],
                },
                ReadAssignment {
                    read: "rejected-read".into(),
                    rank: AssignRank::Species,
                    requires_terminal_evidence: false,
                    taxids: vec![102],
                    depths: vec![0.1],
                    mapped_lens: vec![100.0],
                    sampling_scales: vec![1.0],
                },
            ],
            acc_intervals: HashMap::from([
                (good_key.clone(), vec![(0, 100)]),
                (rejected_key.clone(), vec![(0, 10)]),
            ]),
            acc_fixed_chunk_intervals: HashMap::from([
                (good_key.clone(), vec![(0, 100)]),
                (rejected_key.clone(), vec![(0, 10)]),
            ]),
            acc_flex_chunk_intervals: HashMap::from([
                (good_key, vec![(0, 67)]),
                (rejected_key, vec![(0, 10)]),
            ]),
            ..AssignmentResult::default()
        };
        let sam_data = SamData {
            read_alignments: Vec::new(),
            asm_lengths: Arc::new(HashMap::from([
                ("asm101".to_string(), 1_000.0),
                ("asm102".to_string(), 1_000.0),
            ])),
            min_breadth: 0.05,
            min_chunk_breadth: 0.05,
            safeguard_min_chunk_breadth: 0.05,
            has_scaled_genomes: false,
            total_reads: 2,
            aligned_reads: 2,
            mapping_rate: 1.0,
            fallback_mode: false,
            species_to_genus: Arc::new(HashMap::new()),
            terminal_ranks: Arc::new(HashMap::new()),
            skipped_reference_at_load: 0,
            evidence_policy: None,
        };
        (sam_data, assignment)
    }

    #[test]
    fn raw_em_refines_all_quality_qualified_candidates() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let assignments = vec![
            ReadAssignment {
                read: "unique-one".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: false,
                taxids: vec![101],
                depths: vec![1.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
            },
            ReadAssignment {
                read: "unique-two".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: false,
                taxids: vec![102],
                depths: vec![3.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
            },
            ReadAssignment {
                read: "ambiguous".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: false,
                taxids: vec![101, 102],
                depths: vec![2.0, 2.0],
                mapped_lens: vec![100.0, 100.0],
                sampling_scales: vec![1.0, 1.0],
            },
        ];
        let totals = estimate_unfiltered_taxon_totals(
            &assignments,
            &taxonomy,
            &mut HashMap::new(),
            200,
            false,
        )
        .unwrap();

        approx_equal(totals[&101].support.observed, 1.25);
        approx_equal(totals[&102].support.observed, 1.75);
        approx_equal(totals[&101].depth, 1.5);
        approx_equal(totals[&102].depth, 4.5);
    }

    #[test]
    fn raw_em_distributes_across_mixed_terminal_ranks_without_parent_rollup() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let mut assignments = vec![ReadAssignment {
            read: "unique-genus".into(),
            rank: AssignRank::Genus,
            requires_terminal_evidence: true,
            taxids: vec![111],
            depths: vec![2.0],
            mapped_lens: vec![100.0],
            sampling_scales: vec![1.0],
        }];
        for taxid in [1111, 1112, 1113] {
            assignments.push(ReadAssignment {
                read: format!("unique-species-{taxid}").into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![taxid],
                depths: vec![1.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
            });
        }
        assignments.push(ReadAssignment {
            read: "mixed-rank-ambiguous".into(),
            rank: AssignRank::Mixed,
            requires_terminal_evidence: true,
            taxids: vec![111, 1111, 1112, 1113],
            depths: vec![1.0; 4],
            mapped_lens: vec![100.0; 4],
            sampling_scales: vec![1.0; 4],
        });

        let totals = estimate_unfiltered_taxon_totals(
            &assignments,
            &taxonomy,
            &mut HashMap::new(),
            200,
            false,
        )
        .unwrap();

        approx_equal(totals[&111].depth, 2.4);
        for taxid in [1111, 1112, 1113] {
            approx_equal(totals[&taxid].depth, 1.2);
        }
    }

    #[test]
    fn raw_em_rejects_inconsistent_candidate_vectors() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let assignments = vec![ReadAssignment {
            read: "broken".into(),
            rank: AssignRank::Species,
            requires_terminal_evidence: false,
            taxids: vec![101, 102],
            depths: vec![1.0],
            mapped_lens: vec![100.0, 100.0],
            sampling_scales: vec![1.0, 1.0],
        }];
        let error = estimate_unfiltered_taxon_totals(
            &assignments,
            &taxonomy,
            &mut HashMap::new(),
            200,
            false,
        )
        .unwrap_err();
        assert!(error
            .to_string()
            .contains("inconsistent raw-profile candidate vectors"));
    }

    #[test]
    fn raw_rows_retain_rejected_native_witness_measurements() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let profiler =
            CommunityProfiler::new(make_cfg(Sequencer::Illumina, None), taxonomy).unwrap();
        let estimate = RawProfileEstimate {
            totals: HashMap::from([(
                102,
                TaxonTotals {
                    support: ReadSupport {
                        observed: 0.25,
                        equivalent: 0.25,
                    },
                    depth: 0.02,
                    mapped_len: 20.0,
                },
            )]),
            coverage: HashMap::new(),
        };
        let decisions = HashMap::from([(
            102,
            WitnessDecision {
                metrics: EvidenceCoverage {
                    breadth: 0.01234,
                    expected_breadth: 0.05678,
                    fixed_chunk: 0.02345,
                    flex_chunk: 0.03456,
                    expected_flex_chunk: 0.06789,
                    breadth_pvalue: Some(0.0002),
                    chunk_pvalue: Some(0.25),
                    ..EvidenceCoverage::default()
                },
                standard: false,
                rescues: RescueKinds::default(),
                required_premap_guard: PremapChunkGuard::None,
            },
        )]);

        let rows =
            profiler.build_raw_profile_rows(&estimate, &decisions, &HashMap::new(), None, None);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].taxid, 102);
        assert_eq!(rows[0].equivalent_mapcount, 0.25);
        assert_eq!(rows[0].abundance, 100.0);
        assert_eq!(rows[0].breadth, 0.01234);
        assert_eq!(rows[0].expected_breadth, 0.05678);
        assert_eq!(rows[0].fixed_chunk, 0.02345);
        assert_eq!(rows[0].flex_chunk, 0.03456);
        assert_eq!(rows[0].expected_flex_chunk, 0.06789);
    }

    #[test]
    fn enabling_raw_profile_does_not_change_production_outputs() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let baseline_prefix = directory.path().join("baseline");
        let raw_prefix = directory.path().join("with-raw");

        let mut baseline_cfg = make_cfg(Sequencer::Illumina, None);
        baseline_cfg.outprefix = baseline_prefix.to_string_lossy().into_owned();
        baseline_cfg.save_classifications = true;
        baseline_cfg.breadth = Some(0.05);
        baseline_cfg.chunk_breadth = Some(0.05);
        let baseline = CommunityProfiler::new(baseline_cfg, Arc::clone(&taxonomy)).unwrap();
        let (sam_data, assignment) = raw_output_fixture();
        assert_eq!(baseline.write_outputs(sam_data, assignment).unwrap(), 1);

        let mut raw_cfg = make_cfg(Sequencer::Illumina, None);
        raw_cfg.outprefix = raw_prefix.to_string_lossy().into_owned();
        raw_cfg.save_classifications = true;
        raw_cfg.keep_raw = true;
        raw_cfg.breadth = Some(0.05);
        raw_cfg.chunk_breadth = Some(0.05);
        let with_raw = CommunityProfiler::new(raw_cfg, taxonomy).unwrap();
        let (sam_data, assignment) = raw_output_fixture();
        assert_eq!(with_raw.write_outputs(sam_data, assignment).unwrap(), 1);

        let baseline_profile =
            std::fs::read(format!("{}.profile.txt", baseline_prefix.to_string_lossy())).unwrap();
        let raw_profile =
            std::fs::read(format!("{}.profile.txt", raw_prefix.to_string_lossy())).unwrap();
        assert_eq!(baseline_profile, raw_profile);

        let baseline_classify = std::fs::read(format!(
            "{}.classify.txt",
            baseline_prefix.to_string_lossy()
        ))
        .unwrap();
        let raw_classify =
            std::fs::read(format!("{}.classify.txt", raw_prefix.to_string_lossy())).unwrap();
        assert_eq!(baseline_classify, raw_classify);
        let classify_text = String::from_utf8(raw_classify).unwrap();
        assert!(classify_text.contains("accepted-read\t101"));
        assert!(!classify_text.contains("rejected-read"));

        let raw_table =
            std::fs::read_to_string(format!("{}.rprofile.txt", raw_prefix.to_string_lossy()))
                .unwrap();
        assert!(raw_table.contains("Species one\t101\t"));
        assert!(raw_table.contains("Species two\t102\t"));
    }

    #[test]
    fn final_minimum_support_controls_profile_classification_and_count_together() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let output = directory.path().join("minimum-support");
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.outprefix = output.to_string_lossy().into_owned();
        cfg.save_classifications = true;
        cfg.min_reads = Some(2);
        let profiler = CommunityProfiler::new(cfg, taxonomy).unwrap();
        let (sam_data, assignment) = raw_output_fixture();

        assert_eq!(profiler.write_outputs(sam_data, assignment).unwrap(), 0);
        assert!(
            std::fs::read(format!("{}.profile.txt", output.to_string_lossy()))
                .unwrap()
                .is_empty()
        );
        assert!(
            std::fs::read(format!("{}.classify.txt", output.to_string_lossy()))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn community_filter_uses_unrounded_chunk_breadth() {
        let cfg = make_cfg(Sequencer::Illumina, None);
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
        let selected = HashSet::new();
        let lca = HashMap::new();
        assert!(!community_profile_row_passes(
            &row(0.499_996),
            &cfg,
            0.5,
            &selected,
            &lca,
        ));
        assert!(community_profile_row_passes(
            &row(0.500_004),
            &cfg,
            0.5,
            &selected,
            &lca,
        ));
    }

    #[test]
    fn unresolved_native_terminal_does_not_enter_legacy_zero_breadth_path() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let output = directory.path().join("unresolved-native");
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.outprefix = output.to_string_lossy().into_owned();
        cfg.breadth = Some(0.0);
        cfg.chunk_breadth = Some(0.0);
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let profiler = CommunityProfiler::new(cfg, taxonomy).unwrap();
        let sam_data = SamData {
            read_alignments: Vec::new(),
            asm_lengths: Arc::new(HashMap::from([("asm101".to_string(), 1_000.0)])),
            min_breadth: 0.0,
            min_chunk_breadth: 0.0,
            safeguard_min_chunk_breadth: 0.0,
            has_scaled_genomes: false,
            total_reads: 1,
            aligned_reads: 1,
            mapping_rate: 1.0,
            fallback_mode: false,
            species_to_genus: Arc::new(HashMap::new()),
            terminal_ranks: Arc::new(HashMap::from([(101, ProfileRank::Species)])),
            skipped_reference_at_load: 0,
            evidence_policy: Some(policy),
        };
        let assignment = AssignmentResult {
            read_assignments: vec![ReadAssignment {
                read: "unresolved".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![101],
                depths: vec![0.1],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
            }],
            ..AssignmentResult::default()
        };

        assert_eq!(profiler.write_outputs(sam_data, assignment).unwrap(), 0);
        assert!(
            std::fs::read(format!("{}.profile.txt", output.to_string_lossy()))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn matching_terminal_taxid_does_not_reclassify_legacy_evidence() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let output = directory.path().join("legacy-taxid-collision");
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.outprefix = output.to_string_lossy().into_owned();
        cfg.breadth = Some(0.0);
        cfg.chunk_breadth = Some(0.0);
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let profiler = CommunityProfiler::new(cfg, taxonomy).unwrap();
        let key = AccKey {
            taxid: 101,
            asm: Arc::<str>::from("asm101"),
            acc: Arc::<str>::from("ref101"),
        };
        let sam_data = SamData {
            read_alignments: Vec::new(),
            asm_lengths: Arc::new(HashMap::from([("asm101".to_string(), 1_000.0)])),
            min_breadth: 0.0,
            min_chunk_breadth: 0.0,
            safeguard_min_chunk_breadth: 0.0,
            has_scaled_genomes: false,
            total_reads: 1,
            aligned_reads: 1,
            mapping_rate: 1.0,
            fallback_mode: true,
            species_to_genus: Arc::new(HashMap::new()),
            terminal_ranks: Arc::new(HashMap::from([(101, ProfileRank::Species)])),
            skipped_reference_at_load: 0,
            evidence_policy: Some(policy),
        };
        let assignment = AssignmentResult {
            read_assignments: vec![ReadAssignment {
                read: "legacy-fallback".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: false,
                taxids: vec![101],
                depths: vec![0.1],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
            }],
            acc_intervals: HashMap::from([(key.clone(), vec![(0, 100)])]),
            acc_fixed_chunk_intervals: HashMap::from([(key.clone(), vec![(0, 100)])]),
            acc_flex_chunk_intervals: HashMap::from([(key, vec![(0, 67)])]),
            ..AssignmentResult::default()
        };

        assert_eq!(profiler.write_outputs(sam_data, assignment).unwrap(), 1);
        let profile =
            std::fs::read_to_string(format!("{}.profile.txt", output.to_string_lossy())).unwrap();
        assert_eq!(profile.lines().count(), 1);
        assert!(profile.contains("Species one\t101\tspecies\t"));
    }

    #[test]
    fn native_fixed_flank_uses_validated_radius_for_all_read_inputs() {
        let mut illumina = make_cfg(Sequencer::Illumina, None);
        assert_eq!(
            native_fixed_flank_len(&illumina),
            NATIVE_FIXED_FLANK_LEN,
            "single-end Illumina uses the validated one-kilobase radius"
        );

        illumina.is_paired = true;
        assert_eq!(
            native_fixed_flank_len(&illumina),
            NATIVE_FIXED_FLANK_LEN,
            "paired-end evidence uses the validated one-kilobase radius"
        );

        for sequencer in [Sequencer::Nanopore, Sequencer::PacBio] {
            let long_read = make_cfg(sequencer, None);
            assert_eq!(
                native_fixed_flank_len(&long_read),
                NATIVE_FIXED_FLANK_LEN,
                "{sequencer} evidence must not inherit the two-kilobase single-end radius"
            );
        }

        let assembly = make_cfg(Sequencer::Assembly, None);
        assert_eq!(
            native_fixed_flank_len(&assembly),
            ASSEMBLY_FIXED_FLANK_LEN,
            "assembly profiling is outside the recalibrated native evidence policy"
        );
    }

    #[test]
    fn conflicting_terminal_rank_for_one_taxid_is_rejected() {
        let mut ranks = HashMap::new();
        record_terminal_rank(&mut ranks, 111, ProfileRank::Genus, "genus-ref").unwrap();
        record_terminal_rank(&mut ranks, 111, ProfileRank::Genus, "second-genus-ref").unwrap();
        let error =
            record_terminal_rank(&mut ranks, 111, ProfileRank::Species, "species-ref").unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("conflicting terminal ranks for taxid 111"));
        assert!(message.contains("genus"));
        assert!(message.contains("species"));
    }

    fn native_mapping(
        species: u32,
        reference_id: u32,
        assembly_id: u32,
        identity: f32,
        mapped_len: u32,
        read_len: u32,
        role: MappingRole,
    ) -> AlignmentEntry {
        native_mapping_at_rank(
            species,
            ProfileRank::Species,
            reference_id,
            assembly_id,
            identity,
            mapped_len,
            read_len,
            role,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn native_mapping_at_rank(
        species: u32,
        profile_rank: ProfileRank,
        reference_id: u32,
        assembly_id: u32,
        identity: f32,
        mapped_len: u32,
        read_len: u32,
        role: MappingRole,
    ) -> AlignmentEntry {
        AlignmentEntry {
            asm: Arc::<str>::from(format!("asm{assembly_id}")),
            acc: Arc::<str>::from(format!("ref{reference_id}")),
            profile_taxid: species,
            profile_rank,
            reference_id,
            assembly_id,
            role,
            pair_id: if matches!(role, MappingRole::Mate1 | MappingRole::Mate2) {
                reference_id.max(1)
            } else {
                0
            },
            candidate_boundary: metam::mapper::CandidateBoundary::COMPLETE,
            sampling_scale: 1.0,
            asm_len: 1_000.0,
            identity,
            fraction: mapped_len as f32 / read_len as f32,
            mapped_len,
            span: (100, 100 + mapped_len),
            fixed_chunk: (1, 300 + mapped_len),
            flex_chunk: (1, 500 + mapped_len),
        }
    }

    fn witness(
        fixed_chunk: f64,
        standard: bool,
        rescue: bool,
        required_premap_guard: PremapChunkGuard,
    ) -> WitnessDecision {
        WitnessDecision {
            metrics: EvidenceCoverage {
                fixed_chunk,
                ..EvidenceCoverage::default()
            },
            standard,
            rescues: RescueKinds {
                platform: rescue,
                ..RescueKinds::default()
            },
            required_premap_guard,
        }
    }

    fn low_abundance_witness(fixed_chunk: f64) -> WitnessDecision {
        WitnessDecision {
            metrics: EvidenceCoverage {
                fixed_chunk,
                ..EvidenceCoverage::default()
            },
            standard: false,
            rescues: RescueKinds {
                low_abundance_pe: true,
                ..RescueKinds::default()
            },
            required_premap_guard: PremapChunkGuard::FivePercent,
        }
    }

    #[test]
    fn guarded_map_requires_independent_coverage_on_the_same_assembly() {
        let species = 101;
        let assembly_one = PresenceKey {
            taxid: species,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        };
        let assembly_two = PresenceKey {
            taxid: species,
            profile_rank: ProfileRank::Species,
            assembly_id: 2,
        };
        let baseline_selected = HashSet::new();
        let mut baseline = HashMap::from([
            (
                assembly_one,
                witness(0.70, false, false, PremapChunkGuard::None),
            ),
            (
                assembly_two,
                witness(0.20, false, false, PremapChunkGuard::None),
            ),
        ]);
        let post = vec![(
            assembly_two,
            witness(0.70, true, false, PremapChunkGuard::Full),
        )];

        let additions = guarded_map_additions(&baseline_selected, &baseline, &post, 0.60);
        assert!(
            additions.is_empty(),
            "coverage on assembly one must not authorize MAP evidence on assembly two"
        );

        baseline.insert(
            assembly_two,
            witness(0.61, false, false, PremapChunkGuard::None),
        );
        let additions = guarded_map_additions(&baseline_selected, &baseline, &post, 0.60);
        assert!(additions.contains_key(&species));
    }

    #[test]
    fn guarded_map_uses_mode_specific_post_map_and_amount_guards() {
        let key = PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        };
        let baseline_selected = HashSet::new();
        let se_rescue_only = vec![(key, witness(0.70, false, true, PremapChunkGuard::None))];

        let se_baseline =
            HashMap::from([(key, witness(0.61, false, false, PremapChunkGuard::None))]);
        assert!(
            guarded_map_additions(&baseline_selected, &se_baseline, &se_rescue_only, 0.60,)
                .is_empty(),
            "SE MAP-only calls must pass the full post-MAP standard tier"
        );

        let pe_rescue_only = vec![(key, witness(0.70, false, true, PremapChunkGuard::Half))];
        let pe_baseline =
            HashMap::from([(key, witness(0.31, false, false, PremapChunkGuard::None))]);
        assert!(
            guarded_map_additions(&baseline_selected, &pe_baseline, &pe_rescue_only, 0.60,)
                .contains_key(&101),
            "PE retains its validated half-cutoff plus rescue-tier guard"
        );
    }

    #[test]
    fn equivalent_assembly_witnesses_use_stable_key_tiebreak() {
        let low_key = PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 3,
        };
        let high_key = PresenceKey {
            assembly_id: 9,
            ..low_key
        };
        let mut low = witness(0.50, true, false, PremapChunkGuard::Half);
        low.metrics.flex_chunk = 0.30;
        let mut high = low;
        high.metrics.flex_chunk = 0.90;

        let forward = strongest_witnesses(&[(high_key, high), (low_key, low)]);
        let reverse = strongest_witnesses(&[(low_key, low), (high_key, high)]);
        assert_eq!(forward[&101].metrics.flex_chunk, 0.30);
        assert_eq!(reverse[&101].metrics.flex_chunk, 0.30);
    }

    #[test]
    fn paired_low_abundance_map_uses_relaxed_same_assembly_guard_only_for_its_tier() {
        let species = 101;
        let assembly_one = PresenceKey {
            taxid: species,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        };
        let assembly_two = PresenceKey {
            assembly_id: 2,
            ..assembly_one
        };
        let selected = HashSet::new();
        let low_post = vec![(assembly_one, low_abundance_witness(0.29))];

        let above_guard = HashMap::from([(
            assembly_one,
            witness(0.031, false, false, PremapChunkGuard::None),
        )]);
        assert!(
            guarded_map_additions(&selected, &above_guard, &low_post, 0.60)
                .contains_key(&species),
            "the calibrated PE low-abundance tier accepts independent pre-MAP coverage above 0.05*C"
        );

        let at_guard = HashMap::from([(
            assembly_one,
            witness(0.030, false, false, PremapChunkGuard::None),
        )]);
        assert!(
            guarded_map_additions(&selected, &at_guard, &low_post, 0.60).is_empty(),
            "the independent-coverage comparison is strict"
        );

        let wrong_assembly = HashMap::from([(
            assembly_two,
            witness(0.50, false, false, PremapChunkGuard::None),
        )]);
        assert!(
            guarded_map_additions(&selected, &wrong_assembly, &low_post, 0.60).is_empty(),
            "coverage from another assembly must not authorize the rescue"
        );

        let ordinary_post = vec![(
            assembly_one,
            witness(0.70, false, true, PremapChunkGuard::Half),
        )];
        assert!(
            guarded_map_additions(&selected, &above_guard, &ordinary_post, 0.60).is_empty(),
            "ordinary PE standard/platform evidence retains its 0.50*C pre-MAP guard"
        );
    }

    #[test]
    fn evidence_policy_is_limited_to_validated_community_inputs() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        assert_eq!(
            EvidencePolicy::from_config(&cfg).map(|policy| policy.mode),
            Some(EvidenceMode::Single)
        );

        cfg.is_paired = true;
        assert_eq!(
            EvidencePolicy::from_config(&cfg).map(|policy| policy.mode),
            Some(EvidenceMode::Paired)
        );

        cfg.strain = true;
        assert!(EvidencePolicy::from_config(&cfg).is_none());

        cfg.strain = false;
        cfg.is_paired = false;
        cfg.sequencer = Sequencer::Nanopore;
        assert_eq!(
            EvidencePolicy::from_config(&cfg).map(|policy| policy.mode),
            Some(EvidenceMode::Long)
        );

        cfg.is_paired = true;
        assert!(EvidencePolicy::from_config(&cfg).is_none());
        cfg.is_paired = false;
        cfg.sequencer = Sequencer::Assembly;
        assert!(EvidencePolicy::from_config(&cfg).is_none());
    }

    #[test]
    fn fraction_normalization_is_independent_of_hash_insertion_order() {
        let mut forward = HashMap::new();
        forward.insert(1, 1.0e16);
        forward.insert(2, 1.0);
        forward.insert(3, 1.0);

        let mut reverse = HashMap::new();
        reverse.insert(3, 1.0);
        reverse.insert(2, 1.0);
        reverse.insert(1, 1.0e16);

        normalize_map(&mut forward);
        normalize_map(&mut reverse);
        for taxid in [1, 2, 3] {
            assert_eq!(forward[&taxid].to_bits(), reverse[&taxid].to_bits());
        }
    }

    #[test]
    fn exact_tie_storage_follows_the_resolved_short_read_policy() {
        let rows = vec![
            native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let bundle = || ReadBundle {
            read_name: "tie".into(),
            alignments: rows.clone(),
        };
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = true;

        let default_policy = EvidencePolicy::from_config(&cfg).unwrap();
        let default_result =
            process_evidence_batch(&[bundle()], &cfg, default_policy, false, &HashMap::new());
        assert_eq!(default_result.exact_ties.fragment_count(), 0);

        cfg.exact_ties = ExactTiePolicy::GuardedMap;
        let guarded_policy = EvidencePolicy::from_config(&cfg).unwrap();
        assert!(guarded_policy.collect_exact_ties);
        let guarded_result =
            process_evidence_batch(&[bundle()], &cfg, guarded_policy, false, &HashMap::new());
        assert_eq!(guarded_result.exact_ties.fragment_count(), 1);
        assert!(guarded_result.presence.is_empty());
        for taxid in [101, 102] {
            assert_eq!(
                guarded_result.taxon_evidence[&taxid].resolved_candidate_fragments,
                0
            );
            assert_eq!(guarded_result.taxon_evidence[&taxid].best_fragments, 0);
            assert_eq!(
                guarded_result.taxon_evidence[&taxid].clear_best_fragments,
                0
            );
        }
    }

    #[test]
    fn single_end_policy_keeps_em_candidates_and_uses_half_error_for_coverage() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let bundle = ReadBundle {
            read_name: "single".into(),
            alignments: vec![
                native_mapping(101, 10, 1, 1.000, 100, 100, MappingRole::Unpaired),
                native_mapping(102, 20, 2, 0.995, 100, 100, MappingRole::Unpaired),
            ],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        assert_eq!(policy.mode, EvidenceMode::Single);
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new());

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].taxids, vec![101, 102]);
        assert!(result.presence.contains_key(&PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        }));
        assert!(!result.presence.contains_key(&PresenceKey {
            taxid: 102,
            profile_rank: ProfileRank::Species,
            assembly_id: 2,
        }));
        assert_eq!(result.taxon_evidence[&101].best_fragments, 1);
        assert_eq!(result.taxon_evidence[&101].clear_best_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].clear_best_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].resolved_candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].resolved_candidate_fragments, 1);
    }

    #[test]
    fn candidate_boundary_singleton_stays_in_em_without_unique_breadth() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        let mut mapping = native_mapping(101, 10, 1, 1.0, 100, 100, MappingRole::Unpaired);
        mapping.candidate_boundary = metam::mapper::CandidateBoundary::from_bits(1).unwrap();
        let bundle = ReadBundle {
            read_name: "candidate-boundary".into(),
            alignments: vec![mapping],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new());

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].taxids, vec![101]);
        assert_eq!(result.boundary_ambiguous_fragments, 1);
        assert!(result.uniquely_supported_taxa.is_empty());
        assert!(result.presence.is_empty());
        assert_eq!(result.taxon_evidence[&101].candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&101].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].clear_best_fragments, 0);
    }

    #[test]
    fn mixed_rank_exact_tie_stays_in_em_without_manufacturing_breadth() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.exact_ties = ExactTiePolicy::GuardedMap;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let bundle = ReadBundle {
            read_name: "mixed-rank".into(),
            alignments: vec![
                native_mapping_at_rank(
                    111,
                    ProfileRank::Genus,
                    10,
                    1,
                    0.99,
                    100,
                    100,
                    MappingRole::Unpaired,
                ),
                native_mapping_at_rank(
                    1111,
                    ProfileRank::Species,
                    20,
                    2,
                    0.99,
                    100,
                    100,
                    MappingRole::Unpaired,
                ),
            ],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new());

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].rank, AssignRank::Mixed);
        assert_eq!(result.read_assignments[0].taxids, vec![111, 1111]);
        assert!(result.presence.is_empty());
        assert_eq!(result.exact_ties.fragment_count(), 0);
    }

    #[test]
    fn long_read_mixed_rank_tie_stays_in_em_without_manufacturing_breadth() {
        let mut cfg = make_cfg(Sequencer::Nanopore, None);
        cfg.identity = 0.92;
        cfg.mapped_len = 500;
        cfg.fraction = 0.20;
        let bundle = ReadBundle {
            read_name: "mixed-rank-long".into(),
            alignments: vec![
                native_mapping_at_rank(
                    111,
                    ProfileRank::Genus,
                    10,
                    1,
                    0.98,
                    1_000,
                    2_000,
                    MappingRole::Unpaired,
                ),
                native_mapping_at_rank(
                    1111,
                    ProfileRank::Species,
                    20,
                    2,
                    0.98,
                    1_000,
                    2_000,
                    MappingRole::Unpaired,
                ),
            ],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        assert_eq!(policy.mode, EvidenceMode::Long);
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new());

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].rank, AssignRank::Mixed);
        assert_eq!(result.read_assignments[0].taxids, vec![111, 1111]);
        assert!(result.presence.is_empty());
        assert_eq!(result.exact_ties.fragment_count(), 0);
    }

    #[test]
    fn direct_genus_reference_keeps_genus_rank_and_coordinates() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let bundle = ReadBundle {
            read_name: "genus-only".into(),
            alignments: vec![native_mapping_at_rank(
                111,
                ProfileRank::Genus,
                10,
                1,
                0.99,
                100,
                100,
                MappingRole::Unpaired,
            )],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new());

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].rank, AssignRank::Genus);
        assert!(result.presence.contains_key(&PresenceKey {
            taxid: 111,
            profile_rank: ProfileRank::Genus,
            assembly_id: 1,
        }));
    }

    #[test]
    fn guarded_map_collects_single_end_exact_ties_without_premature_coverage() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.exact_ties = ExactTiePolicy::GuardedMap;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let bundle = ReadBundle {
            read_name: "single-tie".into(),
            alignments: vec![
                native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Unpaired),
                native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Unpaired),
            ],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new());

        assert_eq!(result.exact_ties.fragment_count(), 1);
        assert!(result.presence.is_empty());
        assert_eq!(result.read_assignments[0].taxids, vec![101, 102]);
    }

    #[test]
    fn seeded_map_prior_uses_unique_depth_and_fractional_ambiguous_em() {
        let assignments = vec![
            ReadAssignment {
                read: "unique-1".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![101],
                depths: vec![95.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
            },
            ReadAssignment {
                read: "unique-2".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![102],
                depths: vec![5.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
            },
            ReadAssignment {
                read: "ambiguous".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![101, 102, 103],
                depths: vec![100.0, 100.0, 100.0],
                mapped_lens: vec![100.0; 3],
                sampling_scales: vec![1.0; 3],
            },
        ];
        let eligible = HashSet::from([101, 102]);
        let (priors, iterations) = seeded_short_read_map_priors(&assignments, &eligible, 200);
        assert!(iterations > 0);
        assert!(!priors.contains_key(&103));
        let total = priors[&101] + priors[&102];
        assert!((priors[&101] / total - 0.95).abs() < 1e-12);
        assert!((priors[&102] / total - 0.05).abs() < 1e-12);
    }

    #[test]
    fn paired_exact_tie_stays_in_em_but_creates_no_coordinate_evidence() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = true;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let rows = vec![
            native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let bundle = ReadBundle {
            read_name: "tie".into(),
            alignments: rows,
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new());

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].taxids, vec![101, 102]);
        assert_eq!(result.read_assignments[0].mapped_lens, vec![100.0, 100.0]);
        assert!(result.presence.is_empty());
        assert_eq!(result.taxon_evidence[&101].candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&101].resolved_candidate_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].resolved_candidate_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].clear_best_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].clear_best_fragments, 0);
        assert_eq!(result.exact_ties.fragment_count(), 0);
    }

    #[test]
    fn paired_margin_resolution_uses_mate1_for_abundance_and_coverage() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = true;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let rows = vec![
            native_mapping(101, 10, 1, 1.00, 90, 100, MappingRole::Mate1),
            native_mapping(101, 10, 1, 1.00, 100, 100, MappingRole::Mate2),
            native_mapping(102, 20, 2, 0.98, 100, 100, MappingRole::Mate1),
            native_mapping(102, 20, 2, 0.98, 100, 100, MappingRole::Mate2),
        ];
        let bundle = ReadBundle {
            read_name: "resolved".into(),
            alignments: rows,
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new());

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].taxids, vec![101, 102]);
        assert_eq!(result.read_assignments[0].mapped_lens, vec![90.0, 100.0]);
        let presence = &result.presence[&PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        }];
        assert_eq!(presence.count, 1);
        assert_eq!(presence.reference_span_bases, 90);
        assert_eq!(presence.intervals.len(), 1);
        assert_eq!(result.taxon_evidence[&101].best_fragments, 1);
        assert_eq!(result.taxon_evidence[&101].clear_best_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].clear_best_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].resolved_candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].resolved_candidate_fragments, 1);
        assert!(!result.presence.contains_key(&PresenceKey {
            taxid: 102,
            profile_rank: ProfileRank::Species,
            assembly_id: 2,
        }));
    }

    #[test]
    fn resolved_species_keeps_assemblies_as_independent_witnesses() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = true;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let rows = vec![
            native_mapping(101, 10, 1, 1.00, 100, 100, MappingRole::Mate1),
            native_mapping(101, 10, 1, 1.00, 100, 100, MappingRole::Mate2),
            native_mapping(101, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(101, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let bundle = ReadBundle {
            read_name: "two-assemblies".into(),
            alignments: rows,
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new());

        assert_eq!(result.presence.len(), 2);
        assert!(result.presence.contains_key(&PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        }));
        assert!(result.presence.contains_key(&PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 2,
        }));
    }

    /// Intermediate synthetic alignment used inside tests.
    ///
    /// `AlignmentEntry` stores the pre-parsed
    /// `asm`/`acc`/`species_taxid`/`identity`/`fraction` rather than the raw
    /// `ref_name`/`cigar`. Tests describe alignments as "this ref name, this
    /// CIGAR" for readability and let [`run_single_with_normalize`] perform the
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
        let asm = parts.next()?;
        if asm.is_empty() {
            return None;
        }
        let _ = parts.next();
        let raw_taxid = parts.next()?.parse::<u32>().ok()?;
        let species_taxid = species_normalize
            .get(&raw_taxid)
            .copied()
            .unwrap_or(raw_taxid);
        let acc = parts.next().unwrap_or("");
        let asm_len = asm_lengths.get(asm).copied().unwrap_or(0.0);
        if asm_len <= 0.0 {
            return None;
        }
        let (identity, _, fraction) = cigar_to_identity(&raw.cigar);
        Some(AlignmentEntry {
            asm: Arc::<str>::from(asm),
            acc: Arc::<str>::from(acc),
            profile_taxid: species_taxid,
            profile_rank: ProfileRank::Species,
            reference_id: u32::MAX,
            assembly_id: u32::MAX,
            role: MappingRole::Unpaired,
            pair_id: 0,
            candidate_boundary: metam::mapper::CandidateBoundary::COMPLETE,
            sampling_scale: 1.0,
            asm_len,
            identity: identity as f32,
            fraction: fraction as f32,
            mapped_len: clamp_u32(raw.mapped_len),
            span: (clamp_u32(raw.span.0), clamp_u32(raw.span.1)),
            fixed_chunk: (clamp_u32(raw.span.0), clamp_u32(raw.span.1)),
            flex_chunk: (clamp_u32(raw.span.0), clamp_u32(raw.span.1)),
        })
    }

    // ------------------------------------------------------------------
    // Fallback gate tests
    // ------------------------------------------------------------------

    #[test]
    fn fallback_gate_activates_only_below_mapping_rate_threshold() {
        let cfg = make_cfg(Sequencer::Illumina, None);
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
        let cfg_host = make_cfg(Sequencer::Illumina, Some("9606"));
        assert!(!fallback_gate_active(&cfg_host, 0.0));
        assert!(!fallback_gate_active(&cfg_host, 0.1));
    }

    #[test]
    fn fallback_gate_covers_long_reads_but_not_assembly() {
        // Nanopore and PacBio honour the low-map auto-trigger, so a low
        // mapping rate enables genus fallback for them too.
        for seq in [Sequencer::Nanopore, Sequencer::PacBio] {
            let cfg = make_cfg(seq, None);
            assert!(
                fallback_gate_active(&cfg, 0.0),
                "expected fallback enabled for low-map {seq}"
            );
            assert!(
                !fallback_gate_active(&cfg, 0.99),
                "expected fallback disabled for high-map {seq}"
            );
        }
        // Assembly input is not a read platform -> never auto-triggers.
        let cfg = make_cfg(Sequencer::Assembly, None);
        assert!(!fallback_gate_active(&cfg, 0.0));
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
        let mut cfg = make_cfg(Sequencer::Illumina, None);
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
        // --by-mapped flag present: mapped basis is forced regardless of
        // sequencer, host, or mapping rate.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.by_aligned = true;
        assert!(cfg.resolve_by_aligned(0.99));
        cfg.sequencer = Sequencer::Nanopore;
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

        // Legacy probabilities are computed for species only, so every higher
        // rank is kept under strict when both are absent -- not just genus.
        for rank in ["genus", "family", "order", "class", "phylum"] {
            assert!(
                pvalue_filter(rank, None, None, true),
                "{rank} rows must not be rejected for a test that was never run"
            );
            assert!(!pvalue_filter(rank, None, Some(1e-3), true));
            assert!(!pvalue_filter(rank, Some(1e-3), None, true));
        }
    }

    #[test]
    fn resolve_by_aligned_auto_forces_for_low_map_read_platforms_non_pathogen() {
        // Default cfg has --by-mapped unset (false); auto-rule applies.
        let cfg = make_cfg(Sequencer::Illumina, None);
        assert!(!cfg.by_aligned);
        // Below default 0.30 threshold => aligned basis is auto-forced.
        assert!(cfg.resolve_by_aligned(0.10));
        // At or above threshold => total-read basis.
        assert!(!cfg.resolve_by_aligned(0.30));
        assert!(!cfg.resolve_by_aligned(0.50));
        // Pathogen mode on disables auto-force.
        let cfg_host = make_cfg(Sequencer::Illumina, Some("9606"));
        assert!(!cfg_host.resolve_by_aligned(0.10));
        // Long-read platforms now also auto-force on low map; assembly does not.
        let cfg_np = make_cfg(Sequencer::Nanopore, None);
        assert!(cfg_np.resolve_by_aligned(0.10));
        let cfg_asm = make_cfg(Sequencer::Assembly, None);
        assert!(!cfg_asm.resolve_by_aligned(0.10));
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

        // --by-mapped alone: mapped basis, same equation, no 1.5x.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.by_aligned = true;
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // Automatic low-map fallback on a full DB: aligned basis with 1.5x.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        let gate = fallback_gate_active(&cfg, 0.10);
        assert!(genus_fallback_mode_active(&cfg, gate, false));
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, false, gate);
        assert!(plan.use_aligned_basis);
        assert!(plan.by_aligned_auto);
        assert!(plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, 15_000);

        // Explicit --genus-fallback on non-low-map data without --by-mapped:
        // total basis, no scale.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.genus_fallback = true;
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(!plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, total_reads);
        assert_eq!(plan.estimate_reads, total_reads);

        // Explicit --genus-fallback plus --by-mapped on non-low-map data:
        // mapped basis, same equation, no scale.
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
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        let gate = fallback_gate_active(&cfg, 0.10);
        assert!(!genus_fallback_mode_active(&cfg, gate, true));
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, true, gate);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // User-pinned --min-chunk-breadth is never rescaled.
        cfg.chunk_breadth = Some(0.25);
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, false, gate);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.estimate_reads, aligned_reads);
    }

    #[test]
    fn fallback_gate_active_for_illumina_single_and_paired() {
        // Single-end
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = false;
        assert!(fallback_gate_active(&cfg, 0.1));
        // Paired-end
        cfg.is_paired = true;
        assert!(fallback_gate_active(&cfg, 0.1));
    }

    // ------------------------------------------------------------------
    // Legacy non-LCA rank filter bounds
    // ------------------------------------------------------------------

    #[test]
    fn legacy_rank_filter_lower_bounds_remain_unchanged() {
        assert_eq!(rank_filter_lower_bound("species"), 0.75);
        assert_eq!(rank_filter_lower_bound("genus"), 0.65);
        // Unknown ranks default to species-level stringency.
        assert_eq!(rank_filter_lower_bound("family"), 0.75);
    }

    #[test]
    fn community_chunk_breadth_curve_matches_spec() {
        // Short (Illumina): 0.66*m/(m+1) + 0.24*m/(m+100), plateau 0.90.
        let short = |n| community_chunk_breadth_cutoff(n, true);
        assert!((short(0) - 0.0).abs() < 1e-12);
        assert!((short(10_000) - 0.00656).abs() < 1e-4);
        assert!((short(100_000) - 0.06024).abs() < 1e-4);
        assert!((short(1_000_000) - 0.33238).abs() < 1e-4); // C = 0.33 at 1M reads
        assert!((short(3_300_000) - 0.51418).abs() < 1e-4);
        assert!((short(16_600_000) - 0.65667).abs() < 1e-4);
        // approaches the 0.90 plateau at extreme depth, always below the clamp
        let short_deep = short(100_000_000_000);
        assert!(short_deep > 0.899 && short_deep < 0.95);

        // Long (Nanopore/PacBio): 0.85*m/(m+0.78), plateau 0.85.
        let long = |n| community_chunk_breadth_cutoff(n, false);
        assert!((long(0) - 0.0).abs() < 1e-12);
        assert!((long(10_000) - 0.01076).abs() < 1e-4);
        assert!((long(100_000) - 0.09659).abs() < 1e-4);
        assert!((long(1_000_000) - 0.47753).abs() < 1e-4);
        assert!((long(3_300_000) - 0.68750).abs() < 1e-4);
        assert!((long(16_600_000) - 0.81185).abs() < 1e-4);
        let long_deep = long(100_000_000_000);
        assert!(long_deep > 0.849 && long_deep < 0.85);

        // Monotonic non-decreasing in read count on both platforms.
        for &(a, b) in &[
            (1_000usize, 10_000usize),
            (100_000, 1_000_000),
            (1_000_000, 16_600_000),
        ] {
            assert!(
                community_chunk_breadth_cutoff(b, true) >= community_chunk_breadth_cutoff(a, true)
            );
            assert!(
                community_chunk_breadth_cutoff(b, false)
                    >= community_chunk_breadth_cutoff(a, false)
            );
        }
    }

    #[test]
    fn recall_scales_only_deep_automatic_community_chunk_cutoffs() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.mode = ProfileMode::Recall;

        let at_gate = resolve_community_chunk_breadth(&cfg, 1_000_000);
        approx_equal(at_gate.recall_scale, 1.0);
        approx_equal(at_gate.standard, at_gate.safeguard);
        assert_eq!(at_gate.standard.to_bits(), at_gate.safeguard.to_bits());

        let two_million = resolve_community_chunk_breadth(&cfg, 2_000_000);
        approx_equal(two_million.recall_scale, 0.85);
        approx_equal(
            two_million.standard,
            two_million.safeguard * two_million.recall_scale,
        );

        let five_million = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(five_million.recall_scale, 0.76);
        let deep = resolve_community_chunk_breadth(&cfg, 16_600_000);
        approx_equal(deep.recall_scale, 0.70 + 0.30 / 16.6);
        assert!(deep.recall_scale > SHORT_RECALL_CHUNK_SCALE_FLOOR);

        cfg.sequencer = Sequencer::Nanopore;
        let long_at_gate = resolve_community_chunk_breadth(&cfg, 1_000_000);
        approx_equal(long_at_gate.recall_scale, 1.0);
        let long_deep = resolve_community_chunk_breadth(&cfg, 1_000_001);
        approx_equal(long_deep.recall_scale, LONG_RECALL_CHUNK_SCALE);
        approx_equal(
            long_deep.standard,
            long_deep.safeguard * LONG_RECALL_CHUNK_SCALE,
        );

        cfg.sequencer = Sequencer::PacBio;
        approx_equal(
            resolve_community_chunk_breadth(&cfg, 5_000_000).recall_scale,
            LONG_RECALL_CHUNK_SCALE,
        );
    }

    #[test]
    fn recall_chunk_scaling_respects_user_and_nonstandard_policy_boundaries() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.mode = ProfileMode::Recall;
        cfg.chunk_breadth = Some(0.37);
        let explicit = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(explicit.standard, 0.37);
        approx_equal(explicit.safeguard, 0.37);
        approx_equal(explicit.recall_scale, 1.0);
        assert_eq!(explicit.standard.to_bits(), explicit.safeguard.to_bits());

        cfg.chunk_breadth = None;
        cfg.lowbiomass = true;
        let low_biomass = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(low_biomass.standard, 0.0);
        approx_equal(low_biomass.safeguard, 0.0);
        approx_equal(low_biomass.recall_scale, 1.0);

        cfg.lowbiomass = false;
        cfg.min_reads = Some(3);
        let min_reads = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(min_reads.standard, 0.0);
        approx_equal(min_reads.safeguard, 0.0);
        approx_equal(min_reads.recall_scale, 1.0);

        cfg.min_reads = None;
        cfg.mode = ProfileMode::Default;
        let default = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(default.recall_scale, 1.0);
        assert_eq!(default.standard.to_bits(), default.safeguard.to_bits());
        cfg.mode = ProfileMode::Precision;
        let precision = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(precision.recall_scale, 1.0);
        assert_eq!(precision.standard.to_bits(), precision.safeguard.to_bits());

        cfg.mode = ProfileMode::Recall;
        cfg.sequencer = Sequencer::Assembly;
        approx_equal(
            resolve_community_chunk_breadth(&cfg, 20_000_000).recall_scale,
            1.0,
        );
        cfg.sequencer = Sequencer::Illumina;
        cfg.strain = true;
        approx_equal(
            resolve_community_chunk_breadth(&cfg, 20_000_000).recall_scale,
            1.0,
        );
    }

    #[test]
    fn user_chunk_scale_adjusts_automatic_standard_and_safeguard_cutoffs() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;

        let unchanged = resolve_community_chunk_breadth(&cfg, 5_000_000);
        assert_eq!(unchanged.automatic.to_bits(), unchanged.safeguard.to_bits());
        assert_eq!(unchanged.automatic.to_bits(), unchanged.standard.to_bits());

        cfg.chunk_breadth_scale = 0.5;
        let relaxed = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(relaxed.safeguard, relaxed.automatic * 0.5);
        assert_eq!(relaxed.standard.to_bits(), relaxed.safeguard.to_bits());

        cfg.chunk_breadth_scale = 10.0;
        let capped = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(capped.safeguard, 1.0);
        approx_equal(capped.standard, 1.0);
    }

    #[test]
    fn recall_scaling_composes_after_user_chunk_scaling() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.chunk_breadth_scale = 0.5;
        cfg.mode = ProfileMode::Recall;

        let thresholds = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(thresholds.user_scale, 0.5);
        approx_equal(thresholds.recall_scale, 0.76);
        approx_equal(thresholds.safeguard, thresholds.automatic * 0.5);
        approx_equal(
            thresholds.standard,
            thresholds.safeguard * thresholds.recall_scale,
        );
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
    fn coherent_lca_witness_is_single_strongest_and_order_independent() {
        // Weak-but-standard beats strong-but-rescued, matching witness_is_better.
        let decisions = HashMap::from([
            (10u32, witness(0.90, false, true, PremapChunkGuard::Half)),
            (20u32, witness(0.10, true, false, PremapChunkGuard::Half)),
        ]);
        let (taxid, chosen) = coherent_lca_witness(&[10, 20], &decisions).unwrap();
        assert_eq!(taxid, 20, "a standard witness outranks a rescued one");
        assert!(chosen.standard);
        // Candidate order must not matter.
        assert_eq!(
            coherent_lca_witness(&[20, 10], &decisions).unwrap().0,
            taxid
        );

        // Fully equivalent witnesses resolve to the lowest taxid either way.
        let tied = HashMap::from([
            (77u32, witness(0.50, true, false, PremapChunkGuard::Half)),
            (33u32, witness(0.50, true, false, PremapChunkGuard::Half)),
        ]);
        assert_eq!(coherent_lca_witness(&[77, 33], &tied).unwrap().0, 33);
        assert_eq!(coherent_lca_witness(&[33, 77], &tied).unwrap().0, 33);

        // The whole coverage set comes from that one witness, never a mixture.
        let mixed = HashMap::from([
            (
                1u32,
                WitnessDecision {
                    metrics: EvidenceCoverage {
                        fixed_chunk: 0.80,
                        breadth: 0.10,
                        flex_chunk: 0.10,
                        ..EvidenceCoverage::default()
                    },
                    standard: true,
                    rescues: RescueKinds::default(),
                    required_premap_guard: PremapChunkGuard::Half,
                },
            ),
            (
                2u32,
                WitnessDecision {
                    metrics: EvidenceCoverage {
                        fixed_chunk: 0.20,
                        breadth: 0.99,
                        flex_chunk: 0.99,
                        ..EvidenceCoverage::default()
                    },
                    standard: true,
                    rescues: RescueKinds::default(),
                    required_premap_guard: PremapChunkGuard::Half,
                },
            ),
        ]);
        let (_, chosen) = coherent_lca_witness(&[1, 2], &mixed).unwrap();
        assert_eq!(chosen.metrics.fixed_chunk, 0.80);
        assert_eq!(
            chosen.metrics.breadth, 0.10,
            "breadth must come from the chosen witness, not the per-field maximum"
        );
        assert_eq!(chosen.metrics.flex_chunk, 0.10);

        // No member with coordinate evidence -> no coherent witness.
        assert!(coherent_lca_witness(&[404, 405], &decisions).is_none());
        assert!(coherent_lca_witness(&[], &decisions).is_none());
    }

    #[test]
    fn repeated_lca_groups_retain_one_global_coherent_witness() {
        let weak = (20u32, witness(0.20, true, false, PremapChunkGuard::Half));
        let strong = (10u32, witness(0.80, true, false, PremapChunkGuard::Half));

        let mut forward = HashMap::new();
        retain_coherent_lca_witness(&mut forward, 100, weak);
        retain_coherent_lca_witness(&mut forward, 100, strong);

        let mut reverse = HashMap::new();
        retain_coherent_lca_witness(&mut reverse, 100, strong);
        retain_coherent_lca_witness(&mut reverse, 100, weak);

        assert_eq!(forward[&100].0, 10);
        assert_eq!(reverse[&100].0, 10);
        assert_eq!(forward[&100].1.metrics.fixed_chunk, 0.80);
        assert_eq!(reverse[&100].1.metrics.fixed_chunk, 0.80);

        // Equivalent witnesses use the source taxid as a deterministic final
        // tiebreak, including when they arrive from different LCA groups.
        let tied_high = (77u32, witness(0.50, true, false, PremapChunkGuard::Half));
        let tied_low = (33u32, witness(0.50, true, false, PremapChunkGuard::Half));
        retain_coherent_lca_witness(&mut forward, 200, tied_high);
        retain_coherent_lca_witness(&mut forward, 200, tied_low);
        assert_eq!(forward[&200].0, 33);
    }

    #[test]
    fn legacy_lca_coverage_is_merged_without_replacing_direct_coverage() {
        let direct = HashMap::from([(
            100u32,
            TaxonCoverage {
                breadth: 0.40,
                fixed_chunk: 0.50,
                flex_chunk: 0.60,
            },
        )]);
        let direct_taxa = HashSet::from([100u32]);
        let mut lca = HashMap::new();
        merge_legacy_lca_coverage(
            &mut lca,
            100,
            TaxonCoverage {
                breadth: 0.90,
                fixed_chunk: 0.10,
                flex_chunk: 0.20,
            },
        );
        merge_legacy_lca_coverage(
            &mut lca,
            100,
            TaxonCoverage {
                breadth: 0.30,
                fixed_chunk: 0.80,
                flex_chunk: 0.70,
            },
        );

        assert_eq!(lca[&100].breadth, 0.90);
        assert_eq!(lca[&100].fixed_chunk, 0.80);
        assert_eq!(lca[&100].flex_chunk, 0.70);

        let retained_direct = legacy_output_coverage(100, &direct_taxa, &direct, &lca);
        assert_eq!(retained_direct.breadth, 0.40);
        assert_eq!(retained_direct.fixed_chunk, 0.50);
        assert_eq!(retained_direct.flex_chunk, 0.60);

        let lca_only = legacy_output_coverage(100, &HashSet::new(), &direct, &lca);
        assert_eq!(lca_only.breadth, 0.90);
        assert_eq!(lca_only.fixed_chunk, 0.80);
        assert_eq!(lca_only.flex_chunk, 0.70);
    }

    #[test]
    fn majority_lca_uses_standard_rank_and_platform_ratio_defaults() {
        let short_species =
            standard_lca_ratio_lower_bounds("species", Some(EvidenceMode::Paired)).unwrap();
        assert_eq!(short_species, (0.70, 0.75));
        assert!(ratio_filter_passes_with_defaults(
            0.70,
            1.0,
            0.75,
            1.0,
            short_species,
            None,
            None,
        ));
        assert!(!ratio_filter_passes_with_defaults(
            0.70,
            1.0,
            0.74,
            1.0,
            short_species,
            None,
            None,
        ));

        let long_species =
            standard_lca_ratio_lower_bounds("species", Some(EvidenceMode::Long)).unwrap();
        assert_eq!(long_species, (0.70, 0.70));
        assert!(ratio_filter_passes_with_defaults(
            0.70,
            1.0,
            0.70,
            1.0,
            long_species,
            None,
            None,
        ));

        for rank in [
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "superkingdom",
        ] {
            let defaults =
                standard_lca_ratio_lower_bounds(rank, Some(EvidenceMode::Single)).unwrap();
            assert_eq!(defaults, (0.65, 0.65));
            assert!(ratio_filter_passes_with_defaults(
                0.65, 1.0, 0.65, 1.0, defaults, None, None,
            ));
        }

        assert_eq!(
            standard_lca_ratio_lower_bounds("no rank", Some(EvidenceMode::Long)),
            None
        );
        assert_eq!(standard_lca_ratio_lower_bounds("genus", None), None);
    }

    #[test]
    fn majority_lca_ratio_defaults_still_honor_cli_overrides_and_upper_bound() {
        let defaults =
            standard_lca_ratio_lower_bounds("family", Some(EvidenceMode::Paired)).unwrap();
        assert!(!ratio_filter_passes_with_defaults(
            0.79,
            1.0,
            0.65,
            1.0,
            defaults,
            Some(0.80),
            None,
        ));
        assert!(!ratio_filter_passes_with_defaults(
            1.51, 1.0, 0.65, 1.0, defaults, None, None,
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
            read_name: Arc::<str>::from(read_name),
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
    // identity = 50/50 = 1.0, fraction = 50/(50+50) = 0.5 (below GENUS_FALLBACK_MIN_FRACTION).
    const LOW_FRACTION_CIGAR: &str = "50=50S";

    // ------------------------------------------------------------------
    // Sub-species -> species normalization
    // ------------------------------------------------------------------

    #[test]
    fn subspecies_taxid_is_normalized_to_species() {
        // Raw parts[2] = 9001 is a sub-species of species 101.
        // With `species_normalize = { 9001 -> 101 }` the read must be
        // assigned to 101, not 9001. This is independent of fallback.
        let cfg = make_cfg(Sequencer::Illumina, None);
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
        let cfg = make_cfg(Sequencer::Illumina, None);
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
        let cfg = make_cfg(Sequencer::Illumina, None);
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
        let cfg = make_cfg(Sequencer::Illumina, None);
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
        let cfg = make_cfg(Sequencer::Illumina, None);
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
        let cfg = make_cfg(Sequencer::Illumina, None);
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
    }

    #[test]
    fn genus_fallback_skipped_when_gate_inactive() {
        let cfg = make_cfg(Sequencer::Illumina, None);
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
        let mut cfg = make_cfg(Sequencer::Illumina, None);
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
        let mut cfg = make_cfg(Sequencer::Illumina, None);
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
        let cfg = make_cfg(Sequencer::Illumina, None);
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
    fn genus_fallback_long_reads_pass_on_length_or_fraction() {
        // For long reads the coverage gate is length OR fraction. A Nanopore
        // read whose identity misses the species cutoff falls to the genus
        // fallback, which now accepts it on EITHER the mapped-length floor
        // (cfg.mapped_len) or the fraction floor (GENUS_FALLBACK_MIN_FRACTION),
        // and rejects it only when it meets neither.
        let mut cfg = make_cfg(Sequencer::Nanopore, None);
        cfg.mapped_len = 100;
        cfg.fraction = 0.6;
        cfg.identity = 0.95; // NEAR_MISS identity 0.92 < 0.95 -> falls to genus
        let asm = mk_asm_map();
        let s2g = mk_genus_map();

        // mapped_len 80 < 100 (length fails) but NEAR_MISS fraction = 1.0:
        // accepted via fraction.
        let via_fraction = vec![aln_entry(
            "asmA|1|101|accA|1000",
            NEAR_MISS_CIGAR,
            80,
            (1, 80),
        )];
        assert_eq!(
            run_single(&cfg, &asm, &s2g, true, "r1", via_fraction)
                .expect("accepted on fraction")
                .rank,
            AssignRank::Genus
        );

        // mapped_len 100 >= 100: accepted via length.
        let via_length = vec![aln_entry(
            "asmA|1|101|accA|1000",
            NEAR_MISS_CIGAR,
            100,
            (1, 100),
        )];
        assert_eq!(
            run_single(&cfg, &asm, &s2g, true, "r1", via_length)
                .expect("accepted on length")
                .rank,
            AssignRank::Genus
        );

        // mapped_len 80 < 100 AND LOW_FRACTION fraction 0.5 < 0.60: both fail ->
        // rejected.
        let neither = vec![aln_entry(
            "asmA|1|101|accA|1000",
            LOW_FRACTION_CIGAR,
            80,
            (1, 80),
        )];
        assert!(
            run_single(&cfg, &asm, &s2g, true, "r1", neither).is_none(),
            "long reads: genus fallback rejects an alignment failing both length and fraction"
        );
    }

    #[test]
    fn genus_fallback_rejected_when_fraction_below_060() {
        let cfg = make_cfg(Sequencer::Illumina, None);
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
        let cfg_default = make_cfg(Sequencer::Illumina, None);
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
        let mut cfg_strict = make_cfg(Sequencer::Illumina, None);
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

        let mut cfg_old = make_cfg(Sequencer::Illumina, None);
        cfg_old.genus_identity = 0.85;
        let assignment = run_single(&cfg_old, &asm, &s2g, true, "r1", boundary);
        assert!(
            assignment.is_none(),
            "0.82 identity must be rejected when genus_identity is raised to 0.85"
        );

        // Truly low-identity alignment (0.78) must still be rejected at
        // the new default 0.80 cutoff.
        let cfg_default = make_cfg(Sequencer::Illumina, None);
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
        let cfg = make_cfg(Sequencer::Illumina, None);
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

    #[test]
    fn representative_mapping_carries_its_genome_scale() {
        let lower_identity = AlignmentEntry {
            asm: Arc::<str>::from("full-genome"),
            acc: Arc::<str>::from("contig-a"),
            profile_taxid: 100,
            profile_rank: ProfileRank::Species,
            reference_id: u32::MAX,
            assembly_id: u32::MAX,
            role: MappingRole::Unpaired,
            pair_id: 0,
            candidate_boundary: metam::mapper::CandidateBoundary::COMPLETE,
            sampling_scale: 1.0,
            asm_len: 1_000.0,
            identity: 0.98,
            fraction: 1.0,
            mapped_len: 100,
            span: (1, 101),
            fixed_chunk: (1, 101),
            flex_chunk: (1, 101),
        };
        let mut higher_identity = lower_identity.clone();
        higher_identity.asm = Arc::<str>::from("fractional-genome");
        higher_identity.sampling_scale = 10.0;
        higher_identity.identity = 0.99;

        let mut best = BestMap::new();
        let mut order = Vec::new();
        let mut intervals = IntervalMap::new();
        update_best(
            &mut best,
            &mut order,
            &mut intervals,
            100,
            f64::from(lower_identity.identity),
            0.1,
            &lower_identity,
        );
        update_best(
            &mut best,
            &mut order,
            &mut intervals,
            100,
            f64::from(higher_identity.identity),
            0.1,
            &higher_identity,
        );

        assert_eq!(best.get(&100).map(|entry| entry.3), Some(10.0));
        assert_eq!(
            intervals.get(&100).map(|entry| entry.0.asm.as_ref()),
            Some("fractional-genome")
        );
    }
}

//! Profiler Module - Taxonomic abundance estimation from alignment data
//!
//! This module contains the core profiling logic for Metax. It provides two
//! profiling strategies:
//!
//! - **Community Profiler** (`community.rs`): Species-level profiling that
//!   groups reads by species taxonomy ID. Best for typical metagenomic samples.
//!
//! - **Strain Profiler** (`strain.rs`): Strain-level profiling that can
//!   distinguish between closely related strains within a species. More
//!   computationally intensive but provides finer-grained resolution.
//!
//! # Profiling Pipeline
//!
//! Both profilers follow a similar pipeline:
//!
//! 1. **Load SAM/BAM**: Parse alignment file and extract reference metadata
//! 2. **Compute coverage**: Calculate breadth and depth for each reference
//! 3. **Assign taxonomy**: Map alignments to taxonomic units
//! 4. **Resolve ambiguity**: Use EM algorithm for multi-mapping reads
//! 5. **Filter and output**: Apply coverage filters and write results
//!
//! # Coverage Metrics
//!
//! The profiler computes several coverage metrics:
//!
//! - **Breadth**: Fraction of genome covered by at least one read
//! - **Depth**: Average number of reads covering each position
//! - **Chunk breadth**: Coverage of fixed-size genome chunks (for QC)
//! - **Flex chunk breadth**: Coverage of variable-size chunks (based on √genome_size)

use std::sync::Arc;

use anyhow::Result;
use clap::ValueEnum;

use crate::taxonomy::Taxonomy;

// Sub-modules for different profiling strategies
pub(crate) mod community;
pub(crate) mod strain;

/// Profiling mode that controls the sensitivity/specificity trade-off.
///
/// Each mode sets different default thresholds for alignment filtering:
///
/// - **Recall**: Lower thresholds for higher sensitivity
///   - Accepts more alignments
///
/// - **Precision**: Higher thresholds for higher specificity
///   - Stricter filtering
///
/// - **Default**: Balanced settings for general use
///   - Good trade-off between sensitivity and specificity
///   - Recommended for most applications
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ProfileMode {
    /// High sensitivity mode - more permissive filtering
    Recall,
    /// High specificity mode - stricter filtering
    Precision,
    /// Balanced mode - recommended for general use
    Default,
}

/// Configuration for the taxonomy profiler.
///
/// This structure holds all parameters that control profiling behavior,
/// from input/output paths to filtering thresholds to algorithm options.
/// Values here are the resolved defaults used by the profilers.
///
/// # Filtering Thresholds
///
/// The key thresholds that determine which alignments are retained:
///
/// - `identity`: Minimum sequence identity (matches / gap-compressed length)
/// - `mapped_len`: Minimum number of bases that must align
/// - `fraction`: Minimum fraction of read that must be aligned
/// - `breadth`: Minimum genome breadth of coverage to report a taxon
/// - `chunk_breadth`: Minimum chunked coverage
/// - `min_reads`: Minimum profile read count to report a taxon
/// - `min_oebr`: Minimum observed/expected breadth ratio
/// - `min_coebr`: Minimum observed/expected chunk breadth ratio
#[derive(Debug, Clone)]
pub struct ProfilerConfig {
    /// Path to the input SAM/BAM alignment file
    pub sam: String,

    /// Sequencing platform (Illumina, Nanopore, PacBio)
    /// Affects default threshold selection
    pub sequencer: String,

    /// Number of reads to process per batch
    /// Larger batches improve throughput but use more memory
    pub batch_size: usize,

    /// Whether input is paired-end Illumina data
    /// Affects read name parsing and coverage calculation
    pub is_paired: bool,

    /// Minimum alignment identity threshold (0.0-1.0)
    /// Calculated as: matches / (matches + mismatches + gap_events)
    pub identity: f64,

    /// Minimum number of aligned bases to keep an alignment
    pub mapped_len: usize,

    /// Minimum genome breadth of coverage to report a taxon
    /// If None, defaults to 0.0 (report all)
    pub breadth: Option<f64>,

    /// Minimum chunk breadth threshold
    /// If None, estimated automatically from read count
    pub chunk_breadth: Option<f64>,

    /// Minimum profile read count to report a taxon.
    pub min_reads: Option<usize>,

    /// Minimum observed/expected breadth ratio.
    pub min_oebr: Option<f64>,

    /// Minimum observed/expected chunk breadth ratio.
    pub min_coebr: Option<f64>,

    /// Minimum fraction of read that must align (0.0-1.0)
    pub fraction: f64,

    /// Enable low-biomass sample heuristics
    /// Disables chunk breadth threshold
    pub lowbiomass: bool,

    /// Keep raw (unfiltered) profile output
    pub keep_raw: bool,

    /// Force aligned-read basis for the chunk-breadth estimate.
    pub by_aligned: bool,

    /// Minimum identity accepted for the genus-level fallback path.
    /// Comparison is inclusive (`identity >= genus_identity`). Default 0.80.
    pub genus_identity: f64,

    /// Mapping-rate cut-off for the low-map path. A run is considered
    /// low-map (genus-fallback gate active, auto `--by-aligned` forced on)
    /// when `aligned_reads / total_reads < low_map_rate_threshold`. The
    /// comparison is strict. Default 0.30.
    pub low_map_rate_threshold: f64,

    /// Enable genus-level assignment for alignments that miss species cutoffs.
    pub genus_fallback: bool,

    /// Require both p-value columns to meet the cutoff when both are available.
    pub strict: bool,

    /// Path to pathogen-host mapping TSV (optional)
    pub pathogen_host: Option<String>,

    /// Host taxonomy ID for pathogen filtering (optional)
    pub host: Option<String>,

    /// Number of threads for parallel processing
    pub threads: usize,

    /// Enable verbose logging
    pub verbose: bool,

    /// Enable very verbose logging with full commands
    pub very_verbose: bool,

    /// Output file prefix
    pub outprefix: String,

    /// Profiling mode preset (recall, precision, default)
    pub mode: ProfileMode,

    /// Directory containing taxonomy dump files (optional)
    pub dmp_dir: Option<String>,

    /// Enable strain-level profiling
    pub strain: bool,
}

impl ProfilerConfig {
    /// Return whether the current run should be treated as "low-map" for
    /// the chunk-breadth basis / genus-fallback gate.
    ///
    /// Requires all of:
    /// 1. sequencer is Illumina (case-insensitive),
    /// 2. pathogen-detection mode is off (`host` is `None`),
    /// 3. `mapping_rate < low_map_rate_threshold` (strict `<`).
    ///
    /// Non-finite `mapping_rate` values always return `false`.
    pub fn is_low_map_run(&self, mapping_rate: f64) -> bool {
        mapping_rate.is_finite()
            && self.sequencer.eq_ignore_ascii_case("illumina")
            && self.host.is_none()
            && mapping_rate < self.low_map_rate_threshold
    }

    /// Resolve the effective aligned-vs-total basis for chunk-breadth
    /// estimation.
    ///
    /// - If the user passed `--by-aligned`, always use aligned basis.
    /// - Otherwise auto-force aligned basis for low-map Illumina
    ///   non-pathogen runs (see [`Self::is_low_map_run`]).
    /// - Fall back to total-read basis in all other cases.
    pub fn resolve_by_aligned(&self, mapping_rate: f64) -> bool {
        self.by_aligned || self.is_low_map_run(mapping_rate)
    }

    /// Read count passed into community chunk-breadth estimation when
    /// `--min-cbreadth` is unset.
    ///
    /// Uses aligned vs total reads according to [`Self::resolve_by_aligned`].
    /// Community-specific fallback/scaling rules that depend on reference DB
    /// layout are applied by the community profiler, not by this generic
    /// config helper.
    pub fn chunk_breadth_basis_reads(
        &self,
        mapping_rate: f64,
        num_reads: usize,
        aligned_reads: usize,
    ) -> usize {
        if self.resolve_by_aligned(mapping_rate) {
            aligned_reads
        } else {
            num_reads
        }
    }
}

/// Main profiler that dispatches to the appropriate profiling strategy.
///
/// This is a facade that delegates to either `CommunityProfiler` or
/// `StrainProfiler` based on the `strain` configuration option.
///
/// # Usage
///
/// ```ignore
/// let profiler = Profiler::new(config, taxonomy)?;
/// profiler.run()?;
/// ```
pub struct Profiler {
    /// Configuration for the profiling run
    cfg: ProfilerConfig,

    /// Taxonomy database (shared via Arc for thread safety)
    taxonomy: Arc<Taxonomy>,
}

impl Profiler {
    /// Create a new profiler with the given configuration and taxonomy.
    ///
    /// The taxonomy is wrapped in Arc for efficient sharing across threads.
    ///
    /// # Arguments
    ///
    /// * `cfg` - Profiler configuration
    /// * `taxonomy` - NCBI taxonomy database
    ///
    /// # Returns
    ///
    /// A new Profiler instance ready to run.
    pub fn new(cfg: ProfilerConfig, taxonomy: Taxonomy) -> Result<Self> {
        Ok(Self {
            cfg,
            taxonomy: Arc::new(taxonomy),
        })
    }

    /// Run the profiling pipeline.
    ///
    /// Dispatches to either StrainProfiler or CommunityProfiler based on
    /// the `strain` configuration option.
    ///
    /// # Returns
    ///
    /// Ok(()) on success, or an error if profiling fails.
    pub fn run(&self) -> Result<()> {
        if self.cfg.strain {
            // Strain-level profiling for higher resolution
            strain::StrainProfiler::new(self.cfg.clone(), self.taxonomy.clone())?.run()
        } else {
            // Species-level community profiling (default)
            community::CommunityProfiler::new(self.cfg.clone(), self.taxonomy.clone())?.run()
        }
    }
}

#[cfg(test)]
mod chunk_basis_tests {
    use super::ProfileMode;
    use super::ProfilerConfig;

    fn minimal_cfg() -> ProfilerConfig {
        ProfilerConfig {
            sam: String::new(),
            sequencer: "Illumina".to_string(),
            batch_size: 1,
            is_paired: false,
            identity: 0.97,
            mapped_len: 50,
            breadth: None,
            chunk_breadth: None,
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            fraction: 0.7,
            lowbiomass: false,
            keep_raw: false,
            by_aligned: false,
            genus_identity: 0.80,
            low_map_rate_threshold: 0.30,
            genus_fallback: false,
            strict: false,
            pathogen_host: None,
            host: None,
            threads: 1,
            verbose: false,
            very_verbose: false,
            outprefix: "t".to_string(),
            mode: ProfileMode::Default,
            dmp_dir: None,
            strain: false,
        }
    }

    #[test]
    fn chunk_basis_genus_fallback_non_low_map_uses_total_reads_without_by_aligned() {
        let mut cfg = minimal_cfg();
        cfg.genus_fallback = true;
        // High mapping rate => not low-map; no --by-aligned flag.
        let basis = cfg.chunk_breadth_basis_reads(0.99, 1_000_000, 500_000);
        assert_eq!(
            basis, 1_000_000,
            "explicit genus-fallback alone must not force aligned-read basis on non-low-map runs"
        );
    }

    #[test]
    fn chunk_basis_genus_fallback_non_low_map_honours_by_aligned() {
        let mut cfg = minimal_cfg();
        cfg.genus_fallback = true;
        cfg.by_aligned = true;
        let basis = cfg.chunk_breadth_basis_reads(0.99, 1_000_000, 500_000);
        assert_eq!(
            basis, 500_000,
            "--genus-fallback must not override an explicit --by-aligned basis"
        );
    }

    #[test]
    fn chunk_basis_low_map_unaffected_by_genus_fallback_flag_for_basis_rule() {
        let mut cfg = minimal_cfg();
        cfg.genus_fallback = true;
        // Low-map Illumina: auto by-aligned => aligned basis unless user overrode.
        let basis = cfg.chunk_breadth_basis_reads(0.10, 1_000_000, 100_000);
        assert_eq!(basis, 100_000);
    }
}

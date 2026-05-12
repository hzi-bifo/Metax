//! Metax - Metagenomic Taxonomy Profiler
//!
//! This is the main entry point for the Metax command-line application.
//! Metax is a taxonomy profiler for metagenomic data that uses alignment-based
//! methods with coverage-informed statistics to accurately identify and quantify
//! microbial taxa in sequencing samples.
//!
//! # Subcommands
//!
//! - `profile`: Profile metagenomic reads against a reference database to
//!   generate taxonomic abundance estimates
//! - `index`: Build a maCMD alignment index from a reference FASTA file
//!
//! # Workflow
//!
//! 1. Reads are aligned to a reference database using maCMD aligner
//! 2. Alignments are parsed and filtered based on quality thresholds
//! 3. Reads are assigned to taxa using coverage and identity metrics
//! 4. Ambiguous assignments are resolved using EM (Expectation-Maximization)
//! 5. Final abundance profiles are output with statistical confidence measures

use anyhow::Result;
use clap::{ArgAction, Args, Parser, Subcommand};

/// Parser for a probability-valued CLI option. Rejects non-finite values
/// and values outside `[0.0, 1.0]` at parse time so misuse is caught
/// before any I/O work begins.
fn parse_probability(s: &str) -> Result<f64, String> {
    let value: f64 = s
        .parse()
        .map_err(|e: std::num::ParseFloatError| format!("not a floating-point number: {e}"))?;
    if !value.is_finite() {
        return Err(format!("expected a finite value, got {s}"));
    }
    if !(0.0..=1.0).contains(&value) {
        return Err(format!("expected a value in [0.0, 1.0], got {value}"));
    }
    Ok(value)
}

/// Parser for non-negative finite floating-point thresholds.
fn parse_nonnegative_float(s: &str) -> Result<f64, String> {
    let value: f64 = s
        .parse()
        .map_err(|e: std::num::ParseFloatError| format!("not a floating-point number: {e}"))?;
    if !value.is_finite() {
        return Err(format!("expected a finite value, got {s}"));
    }
    if value < 0.0 {
        return Err(format!("expected a non-negative value, got {value}"));
    }
    Ok(value)
}

use metax::profiler::ProfileMode;
use metax::{build_index, AppConfig, MetaxApp};

/// Main command-line interface structure parsed by clap.
/// Uses derive macros for automatic argument parsing and help generation.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Metax taxonomy profiler",
    propagate_version = true
)]
struct Cli {
    // Top-level CLI dispatch target.
    /// The subcommand to execute (profile or index)
    #[command(subcommand)]
    command: Commands,
}

/// Available subcommands for Metax CLI.
///
/// Each variant corresponds to a major workflow in the application.
#[derive(Subcommand, Debug)]
enum Commands {
    /// Profile metagenomic reads against a reference database.
    ///
    /// This subcommand runs the full profiling pipeline:
    /// 1. Align reads to the reference database (unless reusing existing SAM)
    /// 2. Parse alignments and filter by quality thresholds
    /// 3. Assign taxonomic labels to reads
    /// 4. Generate abundance profile with coverage statistics
    Profile(ProfileArgs),

    /// Build a maCMD index for a reference FASTA.
    ///
    /// Creates the alignment index required for profiling. Optionally
    /// supports subsampling large genomes to reduce index size while
    /// maintaining profiling accuracy.
    Index(IndexArgs),
}

/// Command-line arguments for the `profile` subcommand.
///
/// These arguments control all aspects of the profiling pipeline, from input
/// files to alignment parameters to output filtering thresholds.
///
/// # Threshold Parameters
///
/// The profiler uses several key thresholds to filter alignments:
/// - `identity`: Minimum sequence identity (matches / alignment length)
/// - `mapped_len`: Minimum number of bases that must align
/// - `fraction`: Minimum fraction of the read that must be aligned
/// - `min_breadth`: Minimum genome coverage breadth to report a taxon
/// - `min_cbreadth`: Minimum chunked coverage
///
/// # Profiling Modes
///
/// - `default`: Balanced sensitivity and precision
/// - `recall`: Higher sensitivity
/// - `precision`: Higher specificity
#[derive(Args, Debug)]
struct ProfileArgs {
    /// Path to the maCMD reference database index (metax_db.json).
    /// This is the index created by the `index` subcommand.
    #[arg(long, help = "Path to the maCMD reference database (metax_db.json).")]
    db: Option<String>,

    /// Directory containing NCBI taxonomy dump files (nodes.dmp, names.dmp, merged.dmp).
    /// If not specified, defaults to a `data` subdirectory relative to the executable.
    #[arg(
        long,
        help = "Directory containing the NCBI-style taxonomy dump (dmp files)."
    )]
    dmp_dir: Option<String>,

    /// Input sequence file(s). For Illumina paired-end, provide two files
    /// separated by comma: "read1.fq,read2.fq"
    #[arg(
        short = 'i',
        long,
        value_name = "READS",
        help = "Comma-separated list of input read files (one or two for Illumina paired-end)."
    )]
    in_seq: Option<String>,

    /// Output file prefix. All output files will be named with this prefix
    /// (e.g., prefix.sam, prefix.profile.txt, prefix.classify.txt)
    #[arg(
        short = 'o',
        long,
        value_name = "PREFIX",
        help = "Prefix for output files."
    )]
    outprefix: String,

    /// Number of threads for parallel processing.
    /// Used for both alignment and profiling stages.
    #[arg(
        short = 't',
        long = "threads",
        value_name = "THREADS",
        default_value_t = 20,
        help = "Number of threads to use for alignment and profiling."
    )]
    threads: usize,

    /// Resume mode: skip alignment if output SAM file already exists.
    /// Useful for rerunning profiling with different parameters.
    #[arg(
        short = 'r',
        long = "resume",
        action = ArgAction::SetTrue,
        help = "Resume profiling by reusing existing alignment output if present."
    )]
    resume: bool,

    /// Path to an existing SAM/BAM file to use instead of running alignment.
    /// Skips the alignment step entirely.
    #[arg(
        long = "reuse-sam",
        value_name = "SAM",
        help = "Existing SAM (or compressed SAM) file to reuse instead of running maCMD."
    )]
    reuse_sam: Option<String>,

    /// Sequencing platform. Affects default alignment and filtering parameters.
    /// Supported: Illumina, Nanopore, PacBio
    #[arg(
        long = "sequencer",
        value_name = "TYPE",
        default_value = "Illumina",
        help = "Sequencer type (e.g. Nanopore, PacBio, Illumina)."
    )]
    sequencer: String,

    /// Enable paired-end mode for Illumina reads.
    /// Requires exactly two input files (forward and reverse reads).
    #[arg(
        short = 'p',
        long = "is-paired",
        action = ArgAction::SetTrue,
        help = "Treat Illumina inputs as paired-end reads (expects two files)."
    )]
    is_paired: bool,

    /// Enable strain-level profiling.
    /// Uses a different algorithm that can distinguish between closely related strains.
    #[arg(
        long = "strain",
        action = ArgAction::SetTrue,
        help = "Enable strain-level profiling outputs."
    )]
    strain: bool,

    /// Alignment/profiling mode preset.
    /// Controls the trade-off between sensitivity and specificity.
    #[arg(
        long = "mode",
        value_name = "MODE",
        default_value = "default",
        help = "Alignment mode preset: default, recall, or precision."
    )]
    mode: ProfileMode,

    /// Batch size for parallel read processing.
    /// Larger batches improve throughput but use more memory.
    #[arg(
        long = "batch-size",
        value_name = "N",
        help = "Maximum number of reads to process per batch."
    )]
    batch_size: Option<usize>,

    /// Minimum alignment identity threshold (0.0 to 1.0).
    /// Calculated as: matches / (matches + mismatches + gap_events)
    #[arg(
        long = "identity",
        value_name = "FLOAT",
        help = "Minimum alignment identity threshold for retaining a read."
    )]
    identity: Option<f64>,

    /// Minimum number of aligned bases required to keep an alignment.
    #[arg(
        short = 'm',
        long = "mapped-len",
        value_name = "LEN",
        help = "Minimum mapped read length threshold."
    )]
    mapped_len: Option<usize>,

    /// Minimum genome breadth of coverage (0.0 to 1.0) required to report a taxon.
    /// Breadth = fraction of genome covered by at least one read.
    #[arg(
        short = 'b',
        long = "min-breadth",
        value_name = "FRACTION",
        help = "Minimum breadth of coverage required to report a genome."
    )]
    breadth: Option<f64>,

    /// Minimum chunk breadth threshold.
    /// If not set, estimated automatically.
    #[arg(
        long = "min-cbreadth",
        value_name = "FRACTION",
        help = "Manually set the minimum chunk breadth (overrides automatic estimate)."
    )]
    chunk_breadth: Option<f64>,

    /// Minimum profile read count required to report a taxon.
    /// When set, the automatic minimum chunk-breadth filter is not applied.
    #[arg(
        long = "min-reads",
        value_name = "N",
        help = "Minimum read count required to report a taxon; skips the automatic minimum chunk-breadth filter."
    )]
    min_reads: Option<usize>,

    /// Minimum observed/expected breadth ratio required to report a taxon.
    #[arg(
        long = "min-oebr",
        value_name = "RATIO",
        value_parser = parse_nonnegative_float,
        help = "Minimum observed/expected breadth ratio required to report a taxon."
    )]
    min_oebr: Option<f64>,

    /// Minimum observed/expected chunk breadth ratio required to report a taxon.
    #[arg(
        long = "min-coebr",
        value_name = "RATIO",
        value_parser = parse_nonnegative_float,
        help = "Minimum observed/expected chunk breadth ratio required to report a taxon."
    )]
    min_coebr: Option<f64>,

    /// Minimum fraction of the read that must align (0.0 to 1.0).
    /// Helps filter partial alignments.
    #[arg(
        short = 'f',
        long = "fraction",
        value_name = "FRACTION",
        help = "Minimum aligned fraction a read must cover to be considered."
    )]
    fraction: Option<f64>,

    /// Enable heuristics for low-biomass samples.
    /// Disables minimum chunk breadth threshold (sets to 0).
    #[arg(
        short = 'l',
        long = "lowbiomass",
        action = ArgAction::SetTrue,
        help = "Apply heuristics tuned for low biomass samples."
    )]
    lowbiomass: bool,

    /// Require both p-value columns to meet the cutoff when both are available.
    #[arg(
        long = "strict",
        action = ArgAction::SetTrue,
        help = "Require both cov_prob and chunk_prob to pass the cutoff for rows with p-values."
    )]
    strict: bool,

    /// Keep the raw (unfiltered) profile output file.
    /// Produces an additional rprofile.txt with all taxa before filtering.
    #[arg(
        short = 'k',
        long = "keep-raw",
        action = ArgAction::SetTrue,
        help = "Retain the unfiltered rprofile.txt output alongside the final profile."
    )]
    keep_raw: bool,

    /// Force chunk-breadth estimation to use the aligned-read count.
    #[arg(
        long = "by-aligned",
        action = ArgAction::SetTrue,
        help = "Force aligned-read basis for the minimum chunk breadth estimate (otherwise auto-selected)."
    )]
    by_aligned: bool,

    /// Lower bound on alignment identity accepted by the genus-level
    /// fallback path. Comparison is `identity >= genus_identity`.
    #[arg(
        long = "genus-identity",
        value_name = "FLOAT",
        default_value_t = 0.80,
        value_parser = parse_probability,
        help = "Minimum identity accepted for genus-level fallback candidates (0.0-1.0)."
    )]
    genus_identity: f64,

    /// Mapping-rate threshold below which a run is treated as "low-map".
    /// Comparison is strict: `aligned_reads / total_reads < threshold`.
    #[arg(
        long = "low-map-rate-threshold",
        value_name = "FLOAT",
        default_value_t = 0.30,
        value_parser = parse_probability,
        help = "Strict upper bound on mapping rate for activating the low-map (genus fallback / auto --by-aligned) path."
    )]
    low_map_rate_threshold: f64,

    /// Enable genus-level assignment for alignments that miss species cutoffs.
    #[arg(
        long = "genus-fallback",
        action = ArgAction::SetTrue,
        help = "Force genus-level fallback for alignments that miss species cutoffs (not only low-map Illumina runs)."
    )]
    genus_fallback: bool,

    /// Compress the output SAM file with gzip after profiling.
    #[arg(
        short = 'z',
        long = "compress-sam",
        action = ArgAction::SetTrue,
        help = "Compress the generated SAM file after profiling completes."
    )]
    compress_sam: bool,

    /// Path to TSV file mapping pathogen taxids to host information.
    /// Used with --host for pathogen-focused profiling.
    #[arg(
        long = "pathogen-host",
        value_name = "TSV",
        help = "Optional TSV mapping pathogen taxids to host metadata for annotation."
    )]
    pathogen_host: Option<String>,

    /// NCBI taxonomy ID of the host organism.
    /// When set, filters output to pathogens that can infect this host.
    #[arg(
        long = "host",
        value_name = "TAXID",
        help = "NCBI taxid of the host organism (enables pathogen-specific profiling)."
    )]
    host: Option<String>,

    /// Enable verbose logging output.
    #[arg(
        long = "verbose",
        action = ArgAction::SetTrue,
        help = "Show detailed profiling progress."
    )]
    verbose: bool,

    /// Enable very verbose logging including full command lines.
    #[arg(
        long = "very-verbose",
        action = ArgAction::SetTrue,
        help = "Show extremely detailed logs including full commands."
    )]
    very_verbose: bool,

    /// Additional arguments to pass to maCMD aligner.
    /// Must be specified after `--` on the command line.
    #[arg(
        last = true,
        help = "Additional arguments passed directly to maCMD (use after `--`)."
    )]
    extra_args: Vec<String>,
}

/// Command-line arguments for the `index` subcommand.
///
/// Creates a maCMD alignment index from a reference FASTA file.
/// Optionally supports genome subsampling to reduce index size while
/// maintaining profiling accuracy for large databases.
///
/// # Subsampling Strategy
///
/// When subsampling is enabled (via `--fraction`):
/// 1. Each genome is divided into segments of `segment_length` bases
/// 2. Segments are randomly selected to achieve the target fraction
/// 3. Short segments below `min_length` are filtered for large genomes
/// 4. Original genome size is preserved in headers for count scaling
#[derive(Args, Debug)]
struct IndexArgs {
    /// Subsample fraction for each genome (0.0 to 1.0).
    /// If set, only this fraction of each genome is indexed.
    /// Useful for reducing index size while maintaining sensitivity.
    #[arg(
        short = 'f',
        long = "fraction",
        value_name = "FRACTION",
        help = "Subsample each genome to the given fraction of its size before indexing."
    )]
    fraction: Option<f64>,

    /// Random seed for reproducible subsampling.
    /// Same seed produces identical subsampled FASTA.
    #[arg(
        short = 's',
        long = "seed",
        value_name = "SEED",
        default_value_t = 42u64,
        help = "Seed to control random subsampling (for reproducible indexes)."
    )]
    seed: u64,

    /// Number of threads for parallel subsampling.
    /// Only used during the subsampling phase, not for index building.
    #[arg(
        short = 't',
        long = "threads",
        value_name = "THREADS",
        help = "Number of threads to use when subsampling genomes."
    )]
    threads: Option<usize>,

    /// Length of each subsampled segment in bases.
    /// Genomes are divided into segments of this size during subsampling.
    #[arg(
        short = 'l',
        long = "segment-length",
        value_name = "LEN",
        default_value_t = 50_000usize,
        help = "Length of each subsampled segment (in bases)."
    )]
    segment_length: usize,

    /// Minimum segment length to retain during subsampling.
    /// Segments shorter than this are discarded for genomes > 10x this size.
    #[arg(
        short = 'm',
        long = "min-length",
        value_name = "LEN",
        default_value_t = 3_000usize,
        help = "Minimum segment length to keep for large genomes (in bases)."
    )]
    min_length: usize,

    /// Compress the subsampled FASTA file with gzip after indexing.
    /// The uncompressed file is removed after compression.
    #[arg(
        short = 'z',
        long = "compress",
        action = ArgAction::SetTrue,
        help = "Compress the subsampled FASTA after indexing completes."
    )]
    compress_subsampled: bool,

    /// Input reference FASTA file.
    /// Headers must follow the format: >genome_id|txid|species_txid|seq_id|genome_size
    #[arg(help = "Reference FASTA file to index.")]
    fasta: String,

    /// Output directory for the index files.
    /// Creates a database named `metax_db` in this directory.
    #[arg(
        short = 'o',
        long = "outdir",
        help = "Output directory for database files."
    )]
    outdir: String,
}

/// Conversion from CLI ProfileArgs to the internal AppConfig structure.
///
/// This allows seamless integration between the CLI argument parsing
/// and the library's configuration types. All fields are mapped directly.
impl From<ProfileArgs> for AppConfig {
    fn from(args: ProfileArgs) -> Self {
        AppConfig {
            db: args.db,
            dmp_dir: args.dmp_dir,
            input_sequences: args.in_seq,
            outprefix: args.outprefix,
            threads: args.threads,
            resume: args.resume,
            reuse_sam: args.reuse_sam,
            sequencer: args.sequencer,
            is_paired: args.is_paired,
            strain: args.strain,
            mode: args.mode,
            batch_size: args.batch_size,
            identity: args.identity,
            mapped_len: args.mapped_len,
            breadth: args.breadth,
            chunk_breadth: args.chunk_breadth,
            min_reads: args.min_reads,
            min_oebr: args.min_oebr,
            min_coebr: args.min_coebr,
            fraction: args.fraction,
            lowbiomass: args.lowbiomass,
            keep_raw: args.keep_raw,
            by_aligned: args.by_aligned,
            genus_identity: args.genus_identity,
            low_map_rate_threshold: args.low_map_rate_threshold,
            genus_fallback: args.genus_fallback,
            strict: args.strict,
            compress_sam: args.compress_sam,
            pathogen_host: args.pathogen_host,
            host: args.host,
            verbose: args.verbose,
            very_verbose: args.very_verbose,
            extra_args: args.extra_args,
        }
    }
}

/// Main entry point for the Metax application.
///
/// Parses command-line arguments using clap and dispatches to the
/// appropriate subcommand handler:
///
/// - `profile`: Creates a MetaxApp instance and runs the full profiling pipeline
/// - `index`: Builds a maCMD alignment index from a reference FASTA
///
/// # Returns
///
/// Returns `Ok(())` on success, or an error if any step fails.
/// Errors are formatted with context using anyhow.
fn main() -> Result<()> {
    // Parse command-line arguments using clap's derive macros.
    let cli = Cli::parse();

    // Dispatch to the appropriate subcommand handler.
    match cli.command {
        Commands::Profile(args) => {
            // Convert CLI args to internal config and run profiling
            let app = MetaxApp {
                config: args.into(),
            };
            app.run()
        }
        Commands::Index(args) => {
            // Build alignment index from reference FASTA
            build_index(
                &args.fasta,
                &args.outdir,
                args.fraction,
                args.threads,
                args.seed,
                args.segment_length,
                args.min_length,
                args.compress_subsampled,
            )
        }
    }
}

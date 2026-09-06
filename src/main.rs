//! Metax - Metagenomic Taxonomy Profiler
//!
//! This is the main entry point for the Metax command-line application.
//! Metax is a taxonomy profiler for metagenomic data that uses mapping-based
//! methods with coverage-informed statistics to accurately identify and quantify
//! microbial taxa in sequencing samples.
//!
//! # Subcommands
//!
//! - `profile`: Profile metagenomic reads against a reference database to
//!   generate taxonomic abundance estimates
//! - `index`: Build an index from a reference FASTA file
//! - `info`: Summarize a database or mapping index
//!
//! # Workflow
//!
//! 1. Reads are mapped to a reference database
//! 2. Mappings are parsed and filtered based on quality thresholds
//! 3. Reads are assigned to taxa using coverage and identity metrics
//! 4. Ambiguous assignments are resolved using EM (Expectation-Maximization)
//! 5. Final abundance profiles are output with statistical confidence measures

use std::{
    io::{self, Write},
    path::PathBuf,
};

use anyhow::{bail, Context, Result};
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

/// Parser for positive finite multiplicative factors.
fn parse_positive_float(s: &str) -> Result<f64, String> {
    let value: f64 = s
        .parse()
        .map_err(|e: std::num::ParseFloatError| format!("not a floating-point number: {e}"))?;
    if !value.is_finite() {
        return Err(format!("expected a finite value, got {s}"));
    }
    if value <= 0.0 {
        return Err(format!("expected a value greater than zero, got {value}"));
    }
    Ok(value)
}

use metax::db::{MetamIndexOptions, DEFAULT_SEED_CHUNK_MIB};
use metax::mapper::mapper_options_help;
use metax::mapper::tsv::load_mapping_run_metadata;
use metax::profiler::{ExactTiePolicy, ProfileMode};
use metax::sequencer::Sequencer;
use metax::{build_index, info::run_info, AppConfig, MetaxApp};

/// Main command-line interface structure parsed by clap.
/// Uses derive macros for automatic argument parsing and help generation.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Coverage-aware taxonomic profiling of metagenomic reads",
    long_about = "Metax builds microbial reference databases and profiles short or long metagenomic reads with coverage-aware taxonomic evidence.",
    after_help = "Typical workflow:\n  1. metax index references.fna.gz -o database -t 20\n  2. metax profile -d database -i reads_R1.fq.gz -I reads_R2.fq.gz -o sample -t 20\n  3. metax info -d database",
    propagate_version = true,
    arg_required_else_help = true
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
    /// Map reads and estimate taxonomic abundance.
    #[command(
        display_order = 1,
        override_usage = "metax profile [OPTIONS] --output <PREFIX> [-- <MAPPER_ARG>...]\n       metax profile --mapper-help",
        long_about = "Map metagenomic reads, evaluate taxonomic and coverage evidence, and estimate relative abundance.\n\nPipeline stages:\n  1. Generate candidate mappings\n  2. Evaluate species/strain evidence\n  3. Consider genus fallback for otherwise unassigned reads\n  4. Apply coverage-based report filters",
        after_help = "Examples:\n  Single-end: metax profile -d DB -i reads.fq.gz -o sample -t 20\n  Paired-end: metax profile -d DB -i R1.fq.gz -I R2.fq.gz -o sample -t 20\n  Nanopore:   metax profile -d DB -i reads.fq.gz -o sample -s Nanopore -t 20\n\nSpecies at >=0.95 identity, with genus fallback at >=0.80:\n  metax profile ... --min-identity 0.95 --genus-fallback --genus-identity 0.80\n\nMapper-only candidate-generation overrides must follow '--':\n  metax profile ... -- --top-n 20 --score-ratio 0.85 --min-identity 0.80"
    )]
    Profile(ProfileArgs),

    /// Build a profiling database from reference genomes.
    #[command(
        display_order = 2,
        long_about = "Validate a reference FASTA, build the mapping index, and write the genome/taxonomy metadata required by 'metax profile'. Optional fractional indexing samples large genomes reproducibly before index construction.",
        after_help = "Examples:\n  Full database:       metax index references.fna.gz -o DB -t 20\n  Fractional database: metax index references.fna.gz -o DB -t 20 -f 0.10"
    )]
    Index(IndexArgs),

    /// Summarize a database or mapping index.
    #[command(
        display_order = 3,
        long_about = "Inspect database metadata or the mapping-index header without loading the full index into memory.",
        after_help = "Examples:\n  metax info -d DB\n  metax info -i DB/metax_db.mxi"
    )]
    Info(InfoArgs),
}

/// Command-line arguments for the `profile` subcommand.
///
/// These arguments control all aspects of the profiling pipeline, from input
/// files to mapping parameters to output filtering thresholds.
///
/// # Threshold Parameters
///
/// The profiler uses several key thresholds to evaluate taxonomic evidence:
/// - `identity`: Minimum identity for species/strain evidence
/// - `mapped_len`: Minimum mapped bases for species/strain evidence
/// - `fraction`: Minimum mapped fraction for species/strain evidence
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
    #[arg(
        short = 'd',
        long,
        value_name = "DB",
        help_heading = "Input and output",
        help = "Database directory created by 'metax index' (or its metax_db.tsv.gz metadata file)."
    )]
    db: Option<String>,

    #[arg(
        long = "dmp-dir",
        alias = "taxonomy-dir",
        value_name = "DIR",
        help_heading = "Input and output",
        help = "Directory containing NCBI taxonomy dump files. If omitted, Metax uses ~/.metax/taxdmp and downloads a validated dump there when needed."
    )]
    dmp_dir: Option<String>,

    #[arg(
        short = 'i',
        long = "reads",
        aliases = ["in-seq", "input"],
        value_name = "READS",
        help_heading = "Input and output",
        help = "Single-end/long-read file, or mate 1 for paired-end Illumina data. FASTA/FASTQ may be gzip-compressed."
    )]
    in_seq: Option<String>,

    #[arg(
        short = 'I',
        long = "reads2",
        aliases = ["in-seq2", "input2"],
        value_name = "READS2",
        help_heading = "Input and output",
        help = "Mate-2 Illumina file. Supplying -I automatically enables paired-end mode."
    )]
    in_seq2: Option<String>,

    #[arg(
        short = 'o',
        long = "output",
        alias = "outprefix",
        value_name = "PREFIX",
        required_unless_present = "mapper_help",
        help_heading = "Input and output",
        help = "Output prefix used for the profile, log, and any optional per-read files."
    )]
    outprefix: Option<String>,

    #[arg(
        short = 's',
        long = "sequencer",
        value_name = "PLATFORM",
        value_enum,
        ignore_case = true,
        help_heading = "Analysis mode and resources",
        help = "Sequencing platform; selects suitable mapper and profiler defaults. Defaults to Illumina. For --reuse-mapping it comes from the required file header."
    )]
    sequencer: Option<Sequencer>,

    #[arg(
        long = "mode",
        value_name = "MODE",
        help_heading = "Analysis mode and resources",
        help = "Evidence preset: default (validated balance), recall, or precision. For deep community samples, recall also relaxes only the automatic standard chunk-breadth gate; explicit coverage thresholds remain authoritative. Reused mappings default to their recorded mapper mode; an explicit override changes profiler thresholds only."
    )]
    mode: Option<ProfileMode>,

    #[arg(
        short = 'S',
        long = "strain",
        action = ArgAction::SetTrue,
        help_heading = "Analysis mode and resources",
        help = "[EXPERIMENTAL] Enable strain-aware profiling when strain taxids are present in the database. Strain evidence keeps its own established presets and does not use the validated community evidence policy, so its calls are less well characterised than community output."
    )]
    strain: bool,

    #[arg(
        short = 'l',
        long = "low-biomass",
        alias = "lowbiomass",
        action = ArgAction::SetTrue,
        help_heading = "Analysis mode and resources",
        help = "Use low-biomass reporting heuristics and disable the automatic chunk-breadth filter. An explicit --min-chunk-breadth remains authoritative."
    )]
    lowbiomass: bool,

    #[arg(
        short = 't',
        long = "threads",
        value_name = "N",
        default_value_t = 20,
        help_heading = "Analysis mode and resources",
        help = "Worker threads shared by mapping and profiling."
    )]
    threads: usize,

    #[arg(
        long = "min-identity",
        alias = "identity",
        value_name = "FRACTION",
        value_parser = parse_probability,
        help_heading = "Species and strain evidence",
        help = "Minimum identity for species/strain assignment (0-1). A preset default is used when omitted; the long-read default auto-calibrates from base quality. It does not discard lower-identity mapper output or change --genus-identity."
    )]
    identity: Option<f64>,

    #[arg(
        short = 'm',
        long = "mapped-len",
        value_name = "BASES",
        help_heading = "Species and strain evidence",
        help = "Minimum mapped bases for species/strain assignment. A preset default is used when omitted. On long reads this value also sets the genus-fallback length floor; short reads use a fixed genus floor."
    )]
    mapped_len: Option<usize>,

    #[arg(
        short = 'f',
        long = "min-mapped-fraction",
        aliases = ["fraction", "mapped-fraction"],
        value_name = "FRACTION",
        value_parser = parse_probability,
        help_heading = "Species and strain evidence",
        help = "Minimum mapped read fraction for species/strain assignment (0-1). A preset default is used when omitted. Genus fallback uses its own platform-aware coverage gate."
    )]
    fraction: Option<f64>,

    #[arg(
        long = "exact-ties",
        alias = "pe-exact-ties",
        value_name = "POLICY",
        value_enum,
        help_heading = "Species and strain evidence",
        help = "Treatment of exact equal-score cross-species short-read ties: guarded-map or ambiguous. Guarded MAP is the default for native Illumina community profiling when its coverage safeguard is active."
    )]
    exact_ties: Option<ExactTiePolicy>,

    #[arg(
        short = 'g',
        long = "genus-fallback",
        action = ArgAction::SetTrue,
        help_heading = "Genus fallback and low-map behavior",
        help = "Enable genus assignment in community profiles when a read has no qualifying species mapping. Eligible low-map runs on full databases enable this automatically."
    )]
    genus_fallback: bool,

    #[arg(
        long = "genus-identity",
        value_name = "FRACTION",
        default_value_t = 0.80,
        value_parser = parse_probability,
        help_heading = "Genus fallback and low-map behavior",
        help = "Minimum identity for community genus-fallback evidence (0-1), independent of --min-identity."
    )]
    genus_identity: f64,

    #[arg(
        long = "low-map-rate-threshold",
        value_name = "FRACTION",
        default_value_t = 0.30,
        value_parser = parse_probability,
        help_heading = "Genus fallback and low-map behavior",
        help = "Pre-assignment mapping rate below which low-map behavior activates (0-1). This can enable automatic genus fallback and mapped-read coverage estimation."
    )]
    low_map_rate_threshold: f64,

    #[arg(
        short = 'b',
        long = "min-breadth",
        value_name = "FRACTION",
        value_parser = parse_probability,
        help_heading = "Coverage evidence and filters",
        help = "Minimum fraction of a genome covered by taxonomically accepted mappings (0-1)."
    )]
    breadth: Option<f64>,

    #[arg(
        short = 'c',
        long = "min-chunk-breadth",
        alias = "min-cbreadth",
        value_name = "FRACTION",
        value_parser = parse_probability,
        help_heading = "Coverage evidence and filters",
        help = "Minimum chunked breadth from taxonomically accepted mappings (0-1). Automatically estimated when omitted; an explicit value is never mode-scaled."
    )]
    chunk_breadth: Option<f64>,

    #[arg(
        long = "chunk-breadth-scale",
        value_name = "FACTOR",
        value_parser = parse_positive_float,
        conflicts_with_all = ["chunk_breadth", "min_reads", "lowbiomass"],
        help_heading = "Coverage evidence and filters",
        help = "Multiply the automatically estimated minimum chunk breadth by FACTOR (default 1.0). Values below 1 relax the filter; values above 1 make it stricter. The final cutoff is capped at 1. Conflicts with --min-chunk-breadth, --min-reads, and --low-biomass."
    )]
    chunk_breadth_scale: Option<f64>,

    #[arg(
        long = "by-mapped",
        alias = "by-aligned",
        action = ArgAction::SetTrue,
        help_heading = "Coverage evidence and filters",
        help = "Estimate the automatic chunk-breadth cutoff from mapped reads/fragments instead of all input. Auto-applied on low-map runs (see --low-map-rate-threshold); on a subsampled database, community --genus-fallback estimates an unset cutoff from total reads."
    )]
    by_aligned: bool,

    #[arg(
        short = 'n',
        long = "min-reads",
        value_name = "N",
        help_heading = "Coverage evidence and filters",
        help = "Minimum full-index-equivalent assigned reads for a taxon. Disables the automatic chunk-breadth cutoff unless set explicitly."
    )]
    min_reads: Option<usize>,

    #[arg(
        long = "min-oebr",
        value_name = "RATIO",
        value_parser = parse_nonnegative_float,
        help_heading = "Coverage evidence and filters",
        help = "Minimum observed/expected genome-breadth ratio."
    )]
    min_oebr: Option<f64>,

    #[arg(
        long = "min-coebr",
        value_name = "RATIO",
        value_parser = parse_nonnegative_float,
        help_heading = "Coverage evidence and filters",
        help = "Minimum observed/expected chunk-breadth ratio."
    )]
    min_coebr: Option<f64>,

    #[arg(
        long = "strict",
        action = ArgAction::SetTrue,
        help_heading = "Coverage evidence and filters",
        help = "Require both coverage p-values to pass when both are available."
    )]
    strict: bool,

    #[arg(
        long = "host",
        value_name = "TAXID",
        help_heading = "Pathogen-focused analysis",
        help = "Host NCBI taxid; report only pathogens associated with this host."
    )]
    host: Option<String>,

    #[arg(
        long = "pathogen-table",
        alias = "pathogen-host",
        value_name = "TSV",
        help_heading = "Pathogen-focused analysis",
        help = "TSV that associates pathogen taxids with host metadata; used with --host."
    )]
    pathogen_host: Option<String>,

    #[arg(
        short = 'r',
        long = "resume",
        action = ArgAction::SetTrue,
        help_heading = "Saved and reused data",
        help = "Reuse <PREFIX>.mapping.tsv.gz when present; otherwise warn and map from the supplied reads."
    )]
    resume: bool,

    #[arg(
        long = "reuse-mapping",
        alias = "reuse-mappings",
        value_name = "MAPPING_TSV",
        requires = "db",
        conflicts_with_all = [
            "in_seq",
            "in_seq2",
            "sequencer",
            "is_paired",
            "resume",
            "reuse_sam",
            "keep_mapping_output",
            "stream_mappings",
            "temp_dir"
        ],
        help_heading = "Saved and reused data",
        help = "Profile a versioned Metax mapping.tsv or mapping.tsv.gz without remapping. Requires the matching --db. Platform, pairing, default mode, and the resolved read-level profiler gates come from the file header; explicit profiler options override the recorded gates."
    )]
    reuse_mappings: Option<String>,

    #[arg(
        long = "reuse-sam",
        value_name = "SAM_OR_BAM",
        help_heading = "Saved and reused data",
        help = "Profile an existing SAM/BAM instead of running the built-in mapper (legacy compatibility path)."
    )]
    reuse_sam: Option<String>,

    #[arg(
        long = "save-mapping",
        aliases = ["save-mappings", "keep-mapping-output", "keep-metam-output"],
        action = ArgAction::SetTrue,
        help_heading = "Saved and reused data",
        help = "Write <PREFIX>.mapping.tsv.gz for diagnostics, --resume, or later profiling with --reuse-mapping."
    )]
    keep_mapping_output: bool,

    #[arg(
        long = "stream-mappings",
        action = ArgAction::SetTrue,
        conflicts_with = "keep_mapping_output",
        help_heading = "Performance and robustness",
        help = "Stream mappings directly into the profiler without temporary disk space. This can increase peak memory and mapper contention."
    )]
    stream_mappings: bool,

    #[arg(
        long = "temp-dir",
        value_name = "DIR",
        conflicts_with_all = ["stream_mappings", "keep_mapping_output"],
        help_heading = "Performance and robustness",
        help = "Directory for temporary staged mappings. When omitted, Metax prefers suitable local temporary storage and falls back to network temporary storage or the output directory."
    )]
    temp_dir: Option<String>,

    #[arg(
        short = 'k',
        long = "save-raw-profile",
        alias = "keep-raw",
        action = ArgAction::SetTrue,
        help_heading = "Saved and reused data",
        help = "Also write a post-EM <PREFIX>.rprofile.txt table before taxon-level coverage, evidence-ratio, p-value, and minimum-support filters."
    )]
    keep_raw: bool,

    #[arg(
        long = "save-classification",
        alias = "save-classifications",
        action = ArgAction::SetTrue,
        help_heading = "Saved and reused data",
        help = "Write per-read taxonomic assignments to <PREFIX>.classify.txt. Disabled by default; fractional databases cannot produce this output safely."
    )]
    save_classifications: bool,

    #[arg(
        long = "batch-size",
        value_name = "RECORDS",
        help_heading = "Performance and robustness",
        help = "Profiler work-batch size for file-backed input. Larger batches may use more memory."
    )]
    batch_size: Option<usize>,

    #[arg(
        long = "max-em-iterations",
        value_name = "N",
        default_value_t = 200usize,
        help_heading = "Performance and robustness",
        help = "Maximum iterations for redistributing ambiguous reads; normal runs stop earlier after convergence."
    )]
    max_em_iterations: usize,

    #[arg(
        long = "lenient-mapping-input",
        action = ArgAction::SetTrue,
        help_heading = "Performance and robustness",
        help = "Skip malformed rows in reused mapping files with warnings instead of stopping."
    )]
    lenient_mapping_input: bool,

    #[arg(
        short = 'v',
        long = "verbose",
        action = ArgAction::SetTrue,
        help_heading = "Progress and logging",
        help = "Show detailed resolved settings and profiling progress."
    )]
    verbose: bool,

    #[arg(
        long = "very-verbose",
        action = ArgAction::SetTrue,
        help_heading = "Progress and logging",
        help = "Show all verbose details, including complete mapper settings."
    )]
    very_verbose: bool,

    #[arg(
        long = "progress",
        action = ArgAction::SetTrue,
        conflicts_with = "no_progress",
        help_heading = "Progress and logging",
        help = "Force live progress displays even when stderr is not an interactive terminal."
    )]
    progress: bool,

    #[arg(
        long = "no-progress",
        action = ArgAction::SetTrue,
        help_heading = "Progress and logging",
        help = "Disable live progress displays and periodic mapper progress messages."
    )]
    no_progress: bool,

    #[arg(
        long = "progress-interval",
        value_name = "N",
        default_value_t = 100_000u64,
        help_heading = "Progress and logging",
        help = "Records between count-based progress refreshes for file-backed phases."
    )]
    progress_interval: u64,

    #[arg(
        short = 'p',
        long = "is-paired",
        action = ArgAction::SetTrue,
        hide = true
    )]
    is_paired: bool,

    #[arg(
        long = "mapper-help",
        action = ArgAction::SetTrue,
        help_heading = "Advanced mapper candidate generation",
        help = "Show every supported mapper-only override, its effect, and platform-specific defaults, then exit."
    )]
    mapper_help: bool,

    #[arg(
        last = true,
        value_name = "MAPPER_ARG",
        conflicts_with = "reuse_mappings",
        help_heading = "Advanced mapper candidate generation",
        help = "Mapper-only overrides applied before profiler evidence thresholds. Place them after '--'; run 'metax profile --mapper-help' for the complete list."
    )]
    extra_args: Vec<String>,
}

/// Command-line arguments for the `index` subcommand.
///
/// Creates an index from a reference FASTA file.
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
    #[arg(
        value_name = "REFERENCE_FASTA",
        help = "Reference FASTA to index. Headers: genome|contig|taxonomy[|genome_size]. Gzip input is supported."
    )]
    fasta: String,

    #[arg(
        short = 'o',
        long = "output",
        alias = "outdir",
        value_name = "DIR",
        help_heading = "Input and output",
        help = "Database directory to create. Contains metax_db.mxi and metax_db.tsv.gz."
    )]
    outdir: String,

    #[arg(
        short = 't',
        long = "threads",
        value_name = "N",
        default_value = "1",
        help_heading = "Compute resources",
        help = "Worker threads used for fractional sampling and parallel seed generation."
    )]
    threads: Option<usize>,

    #[arg(
        long = "seed-chunk-mib",
        value_name = "MiB",
        default_value_t = DEFAULT_SEED_CHUNK_MIB,
        help_heading = "Compute resources",
        help = "Approximate reference bases held per parallel seeding batch. Larger values may improve utilization but increase peak RAM."
    )]
    seed_chunk_mib: usize,

    #[arg(
        long = "temp-dir",
        value_name = "DIR",
        help_heading = "Compute resources",
        help = "Existing directory for large temporary seed partitions and sort stores. Use fast local scratch for best indexing performance. Defaults to the hidden database staging area beside --output."
    )]
    temp_dir: Option<String>,

    #[arg(
        short = 'k',
        long = "seed-length",
        value_name = "BASES",
        default_value_t = 15,
        help_heading = "Seed index configuration",
        help = "Length of each syncmer component used by the seed scheme."
    )]
    k: usize,

    #[arg(
        short = 'w',
        long = "seed-window",
        value_name = "N",
        default_value_t = 15,
        help_heading = "Seed index configuration",
        help = "Seed density/window parameter; larger values generally produce fewer seeds and a smaller index."
    )]
    w: usize,

    #[arg(
        long = "max-occ",
        value_name = "N",
        default_value_t = 100,
        help_heading = "Seed index configuration",
        help = "Discard seed keys occurring at more than N reference positions. Lower values shrink the index but may reduce sensitivity."
    )]
    max_occ: usize,

    #[arg(
        long = "seed-mode",
        value_name = "MODE",
        default_value = "randstrobe",
        value_parser = ["randstrobe", "minimizer"],
        ignore_case = true,
        help_heading = "Seed index configuration",
        help = "Seed scheme. Randstrobes are recommended for read mapping."
    )]
    seed_mode: String,

    #[arg(
        long = "hash-bits",
        value_name = "BITS",
        default_value_t = 32,
        help_heading = "Seed index configuration",
        help = "High hash bits stored as the lookup key (24-64). More bits improve discrimination but can enlarge the index."
    )]
    hash_bits: u8,

    #[arg(
        long = "match-both-strobes",
        action = ArgAction::SetTrue,
        conflicts_with = "first_strobe_only",
        hide = true
    )]
    match_both_strobes: bool,

    #[arg(
        long = "first-strobe-only",
        action = ArgAction::SetTrue,
        conflicts_with = "match_both_strobes",
        help_heading = "Randstrobe matching",
        help = "Compatibility mode: match only the first strobe. More sensitive but less specific than the default full-pair verification."
    )]
    first_strobe_only: bool,

    #[arg(
        long = "strobe-start-dist-min",
        value_name = "BASES",
        default_value_t = 30,
        help_heading = "Randstrobe matching",
        help = "Minimum start-to-start distance from strobe 1 to strobe 2."
    )]
    strobe_start_dist_min: u16,

    #[arg(
        long = "strobe-start-dist-max",
        value_name = "BASES",
        default_value_t = 90,
        help_heading = "Randstrobe matching",
        help = "Maximum start-to-start distance from strobe 1 to strobe 2."
    )]
    strobe_start_dist_max: u16,

    #[arg(
        short = 'f',
        long = "fraction",
        value_name = "FRACTION",
        value_parser = parse_probability,
        help_heading = "Optional fractional-reference stage",
        help = "Index approximately this fraction of each eligible genome (0-1) to reduce database size."
    )]
    fraction: Option<f64>,

    #[arg(
        short = 's',
        long = "subsample-seed",
        alias = "seed",
        value_name = "SEED",
        default_value_t = 42u64,
        help_heading = "Optional fractional-reference stage",
        help = "Random seed that makes fractional genome sampling reproducible."
    )]
    seed: u64,

    #[arg(
        short = 'l',
        long = "segment-length",
        value_name = "BASES",
        default_value_t = 50_000usize,
        help_heading = "Optional fractional-reference stage",
        help = "Target length of uniformly distributed sampled genome segments."
    )]
    segment_length: usize,

    #[arg(
        short = 'm',
        long = "min-segment-length",
        alias = "min-length",
        value_name = "BASES",
        default_value_t = 3_000usize,
        help_heading = "Optional fractional-reference stage",
        help = "Discard shorter sampled segments for sufficiently large genomes."
    )]
    min_length: usize,

    #[arg(short = 'z', long = "compress", action = ArgAction::SetTrue, hide = true)]
    compress_subsampled: bool,

    #[arg(
        short = 'q',
        long = "quiet",
        action = ArgAction::SetTrue,
        help_heading = "Progress and logging",
        help = "Suppress index progress and summary messages."
    )]
    quiet: bool,

    #[arg(
        long = "progress",
        action = ArgAction::SetTrue,
        conflicts_with = "no_progress",
        help_heading = "Progress and logging",
        help = "Force live progress displays even when stderr is not an interactive terminal."
    )]
    progress: bool,

    #[arg(
        long = "no-progress",
        action = ArgAction::SetTrue,
        help_heading = "Progress and logging",
        help = "Disable live index progress displays."
    )]
    no_progress: bool,

    #[arg(
        long = "progress-interval",
        value_name = "N",
        default_value_t = 1000u64,
        help_heading = "Progress and logging",
        help = "Reference records between count-based progress refreshes."
    )]
    progress_interval: u64,
}

/// Command-line arguments for the `info` subcommand.
#[derive(Args, Debug)]
struct InfoArgs {
    #[arg(
        short = 'd',
        long,
        value_name = "DB",
        conflicts_with = "index",
        required_unless_present = "index",
        help_heading = "Inspection target",
        help = "Database directory (recommended) or metax_db.tsv.gz metadata file."
    )]
    db: Option<String>,

    #[arg(
        short = 'i',
        long,
        value_name = "INDEX",
        conflicts_with = "db",
        required_unless_present = "db",
        help_heading = "Inspection target",
        help = "Inspect a metax_db.mxi mapping index directly."
    )]
    index: Option<String>,
}

fn resolve_profile_inputs(
    in_seq: Option<String>,
    in_seq2: Option<String>,
    is_paired_flag: bool,
) -> Result<(Option<String>, bool)> {
    let mut reads = Vec::new();
    if let Some(value) = in_seq {
        let parts: Vec<String> = value.split(',').map(|s| s.trim().to_string()).collect();
        if parts.iter().any(|part| part.is_empty()) {
            bail!("input read list contains an empty path");
        }
        if parts.len() > 2 {
            bail!("profile accepts at most two input read files");
        }
        reads.extend(parts);
    }
    if let Some(value) = in_seq2 {
        let value = value.trim().to_string();
        if value.is_empty() {
            bail!("-I/--reads2 cannot be empty");
        }
        if reads.len() != 1 {
            bail!("use either '-i R1 -I R2' or legacy '-i R1,R2', not both");
        }
        reads.push(value);
    }
    if is_paired_flag && reads.len() != 2 {
        bail!("--is-paired requires exactly two input files; use '-i R1 -I R2'");
    }
    let is_paired = reads.len() == 2;
    let input_sequences = if reads.is_empty() {
        None
    } else {
        Some(reads.join(","))
    };
    Ok((input_sequences, is_paired))
}

fn resolve_index_match_both_strobes(
    seed_mode: &str,
    explicit_match_both: bool,
    first_strobe_only: bool,
) -> Result<bool> {
    let is_randstrobe = seed_mode.eq_ignore_ascii_case("randstrobe");
    if !is_randstrobe && (explicit_match_both || first_strobe_only) {
        bail!("--match-both-strobes and --first-strobe-only require --seed-mode randstrobe");
    }
    Ok(is_randstrobe && !first_strobe_only)
}

fn parse_profile_mode_name(value: &str) -> Result<ProfileMode> {
    match value {
        "default" => Ok(ProfileMode::Default),
        "recall" => Ok(ProfileMode::Recall),
        "precision" => Ok(ProfileMode::Precision),
        _ => bail!("unsupported mapper mode '{}' in --reuse-mapping", value),
    }
}

/// Conversion from CLI ProfileArgs to the internal AppConfig structure.
impl TryFrom<ProfileArgs> for AppConfig {
    type Error = anyhow::Error;

    fn try_from(args: ProfileArgs) -> Result<Self> {
        if args.mapper_help {
            bail!("--mapper-help displays mapper options and cannot start profiling");
        }
        let outprefix = args
            .outprefix
            .context("--output is required for profiling")?;
        let reuse_metadata = if let Some(path) = args.reuse_mappings.as_deref() {
            if args.db.is_none() {
                bail!("--reuse-mapping requires the matching --db directory");
            }
            Some(load_mapping_run_metadata(std::path::Path::new(path))?)
        } else {
            None
        };
        let (input_sequences, is_paired) = if args.reuse_mappings.is_some() {
            let metadata = reuse_metadata
                .as_ref()
                .expect("reuse metadata is loaded above");
            (None, metadata.is_paired)
        } else {
            resolve_profile_inputs(args.in_seq, args.in_seq2, args.is_paired)?
        };
        let sequencer = if let Some(metadata) = reuse_metadata.as_ref() {
            metadata.sequencer
        } else {
            args.sequencer.unwrap_or(Sequencer::Illumina)
        };
        let recorded_mode = reuse_metadata
            .as_ref()
            .map(|metadata| parse_profile_mode_name(&metadata.mapper_mode))
            .transpose()?;
        let mode = args
            .mode
            .unwrap_or_else(|| recorded_mode.unwrap_or(ProfileMode::Default));
        let requested_exact_ties = args.exact_ties;
        let chunk_guard_available = !args.lowbiomass
            && !args.chunk_breadth.is_some_and(|value| value <= 0.0)
            && !(args.min_reads.is_some() && args.chunk_breadth.is_none());
        if requested_exact_ties == Some(ExactTiePolicy::GuardedMap) {
            if !sequencer.is_illumina() {
                bail!("--exact-ties guarded-map requires Illumina short-read input");
            }
            if args.strain {
                bail!("--exact-ties guarded-map is supported only by community profiling, not --strain");
            }
            if args.reuse_sam.is_some() {
                bail!("--exact-ties guarded-map requires native Metax coordinate mappings and cannot be used with --reuse-sam");
            }
            if args.lowbiomass {
                bail!("--exact-ties guarded-map requires an active positive chunk-breadth safeguard and cannot be used with --low-biomass");
            }
            if args.chunk_breadth.is_some_and(|value| value <= 0.0) {
                bail!(
                    "--exact-ties guarded-map requires --min-chunk-breadth to be greater than zero"
                );
            }
            if args.min_reads.is_some() && args.chunk_breadth.is_none() {
                bail!("--exact-ties guarded-map requires an automatic or explicit chunk-breadth cutoff; when --min-reads is set, also provide --min-chunk-breadth");
            }
        }
        let exact_ties = requested_exact_ties.unwrap_or_else(|| {
            if sequencer.is_illumina()
                && !args.strain
                && args.reuse_sam.is_none()
                && chunk_guard_available
            {
                ExactTiePolicy::GuardedMap
            } else {
                ExactTiePolicy::Ambiguous
            }
        });
        Ok(AppConfig {
            db: args.db,
            dmp_dir: args.dmp_dir,
            input_sequences,
            outprefix,
            threads: args.threads,
            resume: args.resume,
            reuse_sam: args.reuse_sam,
            reuse_mappings: args.reuse_mappings,
            sequencer,
            is_paired,
            exact_ties,
            strain: args.strain,
            mode,
            batch_size: args.batch_size,
            max_em_iterations: args.max_em_iterations,
            lenient_mapping_input: args.lenient_mapping_input,
            identity: args.identity,
            mapped_len: args.mapped_len,
            breadth: args.breadth,
            chunk_breadth: args.chunk_breadth,
            chunk_breadth_scale: args.chunk_breadth_scale.unwrap_or(1.0),
            min_reads: args.min_reads,
            min_oebr: args.min_oebr,
            min_coebr: args.min_coebr,
            fraction: args.fraction,
            lowbiomass: args.lowbiomass,
            keep_raw: args.keep_raw,
            save_classifications: args.save_classifications,
            by_aligned: args.by_aligned,
            genus_identity: args.genus_identity,
            low_map_rate_threshold: args.low_map_rate_threshold,
            genus_fallback: args.genus_fallback,
            strict: args.strict,
            keep_mapping_output: args.keep_mapping_output,
            stream_mappings: args.stream_mappings,
            temp_dir: args.temp_dir,
            pathogen_host: args.pathogen_host,
            host: args.host,
            verbose: args.verbose,
            very_verbose: args.very_verbose,
            progress: args.progress,
            no_progress: args.no_progress,
            progress_interval: args.progress_interval,
            extra_args: args.extra_args,
        })
    }
}

/// Main entry point for the Metax application.
///
/// Parses command-line arguments using clap and dispatches to the
/// appropriate subcommand handler:
///
/// - `profile`: Creates a MetaxApp instance and runs the full profiling pipeline
/// - `index`: Builds an index from a reference FASTA
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
            if args.mapper_help {
                let mut stdout = io::BufWriter::new(io::stdout().lock());
                stdout.write_all(mapper_options_help().as_bytes())?;
                stdout.flush()?;
                return Ok(());
            }
            // Convert CLI args to internal config and run profiling
            let app = MetaxApp {
                config: args.try_into()?,
            };
            app.run()
        }
        Commands::Index(args) => {
            let match_both_strobes = resolve_index_match_both_strobes(
                &args.seed_mode,
                args.match_both_strobes,
                args.first_strobe_only,
            )?;
            // Build mapping index from reference FASTA
            build_index(
                &args.fasta,
                &args.outdir,
                args.fraction,
                args.threads,
                args.seed,
                args.segment_length,
                args.min_length,
                args.compress_subsampled,
                MetamIndexOptions {
                    k: args.k,
                    w: args.w,
                    max_occ: args.max_occ,
                    hash_bits: args.hash_bits,
                    seed_mode: args.seed_mode,
                    match_both_strobes,
                    strobe_start_dist_min: args.strobe_start_dist_min,
                    strobe_start_dist_max: args.strobe_start_dist_max,
                    seed_chunk_mib: args.seed_chunk_mib,
                    temp_dir: args.temp_dir.map(PathBuf::from),
                    quiet: args.quiet,
                    progress: args.progress,
                    no_progress: args.no_progress,
                    progress_interval: args.progress_interval,
                },
            )
        }
        Commands::Info(args) => run_info(args.db.as_deref(), args.index.as_deref()),
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;
    use std::path::PathBuf;

    use clap::{CommandFactory, Parser};
    use tempfile::NamedTempFile;

    use metax::db::{DatabaseMetadata, ReferenceMetadata};
    use metax::mapper::tsv::{write_mapping_run_metadata, MappingProfilerMetadata};

    use super::{
        mapper_options_help, resolve_index_match_both_strobes, resolve_profile_inputs, AppConfig,
        Cli, Commands, ExactTiePolicy, ProfileMode, Sequencer,
    };

    fn clean_mapping(sequencer: Sequencer, paired: bool, mode: ProfileMode) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        let database = DatabaseMetadata {
            metam_index: PathBuf::new(),
            references: vec![ReferenceMetadata {
                name: "G|ctg|s:8|100".to_string(),
                length: 100,
                genome: "G".to_string(),
                contig: "ctg".to_string(),
                full_genome_size: 100,
                effective_indexed_size: 100,
                taxonomy_field: "s:8".to_string(),
            }],
            is_fractional: false,
        };
        write_mapping_run_metadata(
            file.as_file_mut(),
            sequencer,
            paired,
            mode.as_str(),
            MappingProfilerMetadata {
                strain: false,
                identity: if sequencer.is_illumina() { 0.97 } else { 0.92 },
                mapped_len: if sequencer.is_illumina() { 50 } else { 500 },
                fraction: if sequencer.is_illumina() { 0.70 } else { 0.20 },
                sampled_mean_base_accuracy: None,
            },
            &database,
        )
        .unwrap();
        writeln!(
            file,
            "r1\tG|ctg|s:8|100\t0\t100\t1\t100\t100\t{}",
            if paired { "=" } else { "NA" }
        )
        .unwrap();
        file.flush().unwrap();
        file
    }

    fn long_help(subcommand_name: &str) -> String {
        let mut command = Cli::command();
        command
            .find_subcommand_mut(subcommand_name)
            .unwrap_or_else(|| panic!("missing {subcommand_name} subcommand"))
            .render_long_help()
            .to_string()
    }

    fn assert_appears_in_order(text: &str, needles: &[&str]) {
        let mut previous = 0;
        for needle in needles {
            let relative = text[previous..]
                .find(needle)
                .unwrap_or_else(|| panic!("expected help text to contain {needle:?}\n{text}"));
            previous += relative + needle.len();
        }
    }

    #[test]
    fn paired_input_is_inferred_from_mate2_flag() {
        let (input, is_paired) = resolve_profile_inputs(
            Some("r1.fq.gz".to_string()),
            Some("r2.fq.gz".to_string()),
            false,
        )
        .unwrap();
        assert_eq!(input.as_deref(), Some("r1.fq.gz,r2.fq.gz"));
        assert!(is_paired);
    }

    #[test]
    fn legacy_comma_paired_input_still_works() {
        let (input, is_paired) =
            resolve_profile_inputs(Some("r1.fq.gz,r2.fq.gz".to_string()), None, false).unwrap();
        assert_eq!(input.as_deref(), Some("r1.fq.gz,r2.fq.gz"));
        assert!(is_paired);
    }

    #[test]
    fn is_paired_without_two_inputs_is_rejected() {
        let err = resolve_profile_inputs(Some("r1.fq.gz".to_string()), None, true).unwrap_err();
        assert!(err
            .to_string()
            .contains("--is-paired requires exactly two input files"));
    }

    #[test]
    fn mixed_paired_input_forms_are_rejected() {
        let err = resolve_profile_inputs(
            Some("r1.fq.gz,r2.fq.gz".to_string()),
            Some("other.fq.gz".to_string()),
            false,
        )
        .unwrap_err();
        assert!(err.to_string().contains("not both"));
    }

    #[test]
    fn randstrobe_index_matches_both_strobes_by_default() {
        assert!(resolve_index_match_both_strobes("randstrobe", false, false).unwrap());
        assert!(resolve_index_match_both_strobes("RANDSTROBE", true, false).unwrap());
        assert!(!resolve_index_match_both_strobes("randstrobe", false, true).unwrap());
    }

    #[test]
    fn strobe_matching_flags_are_rejected_for_other_seed_modes() {
        assert!(!resolve_index_match_both_strobes("syncmer", false, false).unwrap());
        assert!(resolve_index_match_both_strobes("syncmer", true, false).is_err());
        assert!(resolve_index_match_both_strobes("minimizer", false, true).is_err());
    }

    #[test]
    fn cli_normalizes_supported_sequencers_and_rejects_unknown_values() {
        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "-o",
            "result",
            "--sequencer",
            "nanopore",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        assert_eq!(args.sequencer, Some(Sequencer::Nanopore));

        let error = Cli::try_parse_from([
            "metax",
            "profile",
            "-o",
            "result",
            "--sequencer",
            "nanoporre",
        ])
        .unwrap_err();
        assert!(error.to_string().contains("invalid value 'nanoporre'"));
    }

    #[test]
    fn reuse_mappings_requires_database_and_infers_run_shape_and_mode() {
        let file = clean_mapping(Sequencer::Illumina, true, ProfileMode::Recall);
        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "--reuse-mapping",
            file.path().to_str().unwrap(),
            "--db",
            "db",
            "-o",
            "result",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        let config = AppConfig::try_from(args).unwrap();
        assert_eq!(config.db.as_deref(), Some("db"));
        assert!(config.input_sequences.is_none());
        assert_eq!(config.sequencer, Sequencer::Illumina);
        assert!(config.is_paired);
        assert_eq!(config.mode, ProfileMode::Recall);
        assert_eq!(config.exact_ties, ExactTiePolicy::GuardedMap);
        assert_eq!(config.reuse_mappings.as_deref(), file.path().to_str());
    }

    #[test]
    fn reuse_mappings_infers_long_read_platform() {
        let file = clean_mapping(Sequencer::Nanopore, false, ProfileMode::Default);
        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "--reuse-mapping",
            file.path().to_str().unwrap(),
            "--db",
            "db",
            "-o",
            "result",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        let config = AppConfig::try_from(args).unwrap();
        assert_eq!(config.sequencer, Sequencer::Nanopore);
        assert!(!config.is_paired);
        assert_eq!(config.exact_ties, ExactTiePolicy::Ambiguous);
    }

    #[test]
    fn reuse_mappings_rejects_mapping_stage_inputs_and_overrides() {
        let file = clean_mapping(Sequencer::Illumina, false, ProfileMode::Default);
        for conflicting in [
            vec!["--reads", "reads.fq.gz"],
            vec!["--save-mapping"],
            vec!["--sequencer", "Illumina"],
            vec!["--is-paired"],
        ] {
            let mut argv = vec![
                "metax",
                "profile",
                "--reuse-mapping",
                file.path().to_str().unwrap(),
                "--db",
                "db",
                "-o",
                "result",
            ];
            argv.extend(conflicting);
            let error = Cli::try_parse_from(argv).unwrap_err();
            assert!(error.to_string().contains("cannot be used with"));
        }

        let error = Cli::try_parse_from([
            "metax",
            "profile",
            "--reuse-mapping",
            file.path().to_str().unwrap(),
            "--db",
            "db",
            "-o",
            "result",
            "--",
            "--top-n",
            "20",
        ])
        .unwrap_err();
        assert!(error.to_string().contains("cannot be used with"));
    }

    #[test]
    fn reuse_mappings_requires_database() {
        let file = clean_mapping(Sequencer::Illumina, false, ProfileMode::Default);
        let error = Cli::try_parse_from([
            "metax",
            "profile",
            "--reuse-mapping",
            file.path().to_str().unwrap(),
            "-o",
            "result",
        ])
        .unwrap_err();
        assert!(error.to_string().contains("--db <DB>"));
    }

    #[test]
    fn reuse_mappings_rejects_files_without_versioned_metadata() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "r1\tref\t0\t10\t1\t10\t10\t=").unwrap();
        file.flush().unwrap();
        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "--reuse-mapping",
            file.path().to_str().unwrap(),
            "--db",
            "db",
            "-o",
            "result",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        let error = AppConfig::try_from(args).unwrap_err();
        assert!(error.to_string().contains("missing the required"));
        assert!(error.to_string().contains("--save-mapping"));
    }

    #[test]
    fn profile_cli_accepts_canonical_workflow_names() {
        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "--db",
            "db",
            "--dmp-dir",
            "taxonomy",
            "--reads",
            "r1.fq.gz",
            "--reads2",
            "r2.fq.gz",
            "--output",
            "sample",
            "--sequencer",
            "Illumina",
            "--mode",
            "recall",
            "--min-identity",
            "0.96",
            "--mapped-len",
            "100",
            "--min-mapped-fraction",
            "0.60",
            "--genus-fallback",
            "--min-breadth",
            "0.10",
            "--min-chunk-breadth",
            "0.20",
            "--min-reads",
            "10",
            "--save-mapping",
            "--save-classification",
            "--save-raw-profile",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        assert_eq!(args.db.as_deref(), Some("db"));
        assert_eq!(args.dmp_dir.as_deref(), Some("taxonomy"));
        assert_eq!(args.in_seq.as_deref(), Some("r1.fq.gz"));
        assert_eq!(args.in_seq2.as_deref(), Some("r2.fq.gz"));
        assert_eq!(args.outprefix.as_deref(), Some("sample"));
        assert_eq!(args.identity, Some(0.96));
        assert_eq!(args.mapped_len, Some(100));
        assert_eq!(args.fraction, Some(0.60));
        assert_eq!(args.breadth, Some(0.10));
        assert_eq!(args.chunk_breadth, Some(0.20));
        assert_eq!(args.chunk_breadth_scale, None);
        assert_eq!(args.min_reads, Some(10));
        assert!(args.genus_fallback);
        assert!(args.keep_mapping_output);
        assert!(args.save_classifications);
        assert!(args.keep_raw);
    }

    #[test]
    fn chunk_breadth_scale_is_positive_optional_and_mutually_exclusive() {
        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "--reads",
            "reads.fq",
            "--output",
            "sample",
            "--chunk-breadth-scale",
            "0.5",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        assert_eq!(args.chunk_breadth_scale, Some(0.5));
        assert_eq!(AppConfig::try_from(args).unwrap().chunk_breadth_scale, 0.5);

        for invalid in ["0", "-1", "NaN", "inf"] {
            assert!(Cli::try_parse_from([
                "metax",
                "profile",
                "--reads",
                "reads.fq",
                "--output",
                "sample",
                "--chunk-breadth-scale",
                invalid,
            ])
            .is_err());
        }

        for conflicting in [
            ["--min-chunk-breadth", "0.1"],
            ["--min-reads", "3"],
            ["--low-biomass", ""],
        ] {
            let mut argv = vec![
                "metax",
                "profile",
                "--reads",
                "reads.fq",
                "--output",
                "sample",
                "--chunk-breadth-scale",
                "0.5",
                conflicting[0],
            ];
            if !conflicting[1].is_empty() {
                argv.push(conflicting[1]);
            }
            assert!(Cli::try_parse_from(argv).is_err());
        }

        let default_cli = Cli::try_parse_from([
            "metax", "profile", "--reads", "reads.fq", "--output", "sample",
        ])
        .unwrap();
        let Commands::Profile(default_args) = default_cli.command else {
            panic!("expected profile command");
        };
        assert_eq!(default_args.chunk_breadth_scale, None);
        assert_eq!(
            AppConfig::try_from(default_args)
                .unwrap()
                .chunk_breadth_scale
                .to_bits(),
            1.0f64.to_bits()
        );
    }

    #[test]
    fn per_read_classifications_are_opt_in() {
        let default_cli = Cli::try_parse_from([
            "metax", "profile", "--reads", "reads.fq", "--output", "sample",
        ])
        .unwrap();
        let Commands::Profile(default_args) = default_cli.command else {
            panic!("expected profile command");
        };
        assert!(!default_args.save_classifications);

        let enabled_cli = Cli::try_parse_from([
            "metax",
            "profile",
            "--reads",
            "reads.fq",
            "--output",
            "sample",
            "--save-classification",
        ])
        .unwrap();
        let Commands::Profile(enabled_args) = enabled_cli.command else {
            panic!("expected profile command");
        };
        assert!(enabled_args.save_classifications);
    }

    #[test]
    fn guarded_map_is_the_default_for_single_and_paired_illumina() {
        let cli =
            Cli::try_parse_from(["metax", "profile", "-i", "reads.fq.gz", "-o", "sample"]).unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        assert_eq!(args.exact_ties, None);
        let config: AppConfig = args.try_into().unwrap();
        assert_eq!(config.exact_ties, ExactTiePolicy::GuardedMap);
        assert!(!config.is_paired);

        let cli = Cli::try_parse_from([
            "metax", "profile", "-i", "r1.fq.gz", "-I", "r2.fq.gz", "-o", "sample",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        let config: AppConfig = args.try_into().unwrap();
        assert_eq!(config.exact_ties, ExactTiePolicy::GuardedMap);
        assert!(config.is_paired);

        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "-i",
            "reads.fq.gz",
            "-o",
            "sample",
            "--exact-ties",
            "ambiguous",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        let config: AppConfig = args.try_into().unwrap();
        assert_eq!(config.exact_ties, ExactTiePolicy::Ambiguous);
    }

    #[test]
    fn automatic_exact_tie_policy_falls_back_safely_when_map_is_inapplicable() {
        for extra in [
            vec!["--sequencer", "Nanopore"],
            vec!["--low-biomass"],
            vec!["--min-reads", "5"],
        ] {
            let mut argv = vec!["metax", "profile", "-i", "reads.fq.gz", "-o", "sample"];
            argv.extend(extra);
            let cli = Cli::try_parse_from(argv).unwrap();
            let Commands::Profile(args) = cli.command else {
                panic!("expected profile command");
            };
            let config: AppConfig = args.try_into().unwrap();
            assert_eq!(config.exact_ties, ExactTiePolicy::Ambiguous);
        }
    }

    #[test]
    fn guarded_map_accepts_se_but_rejects_nonillumina_and_strain_inputs() {
        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "-i",
            "reads.fq.gz",
            "-o",
            "sample",
            "--exact-ties",
            "guarded-map",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        let config = AppConfig::try_from(args).unwrap();
        assert!(!config.is_paired);
        assert_eq!(config.exact_ties, ExactTiePolicy::GuardedMap);

        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "-i",
            "reads.fq.gz",
            "-o",
            "sample",
            "--sequencer",
            "Nanopore",
            "--exact-ties",
            "guarded-map",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        assert!(AppConfig::try_from(args)
            .unwrap_err()
            .to_string()
            .contains("Illumina short-read"));

        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "-i",
            "reads.fq.gz",
            "-o",
            "sample",
            "--strain",
            "--exact-ties",
            "guarded-map",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        assert!(AppConfig::try_from(args)
            .unwrap_err()
            .to_string()
            .contains("community profiling"));
    }

    #[test]
    fn guarded_map_requires_an_independent_chunk_breadth_guard() {
        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "-i",
            "r1.fq.gz",
            "-I",
            "r2.fq.gz",
            "-o",
            "sample",
            "--min-reads",
            "5",
            "--exact-ties",
            "guarded-map",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        assert!(AppConfig::try_from(args)
            .unwrap_err()
            .to_string()
            .contains("also provide --min-chunk-breadth"));
    }

    #[test]
    fn profile_cli_preserves_legacy_long_option_aliases() {
        let cli = Cli::try_parse_from([
            "metax",
            "profile",
            "--db",
            "db",
            "--taxonomy-dir",
            "taxonomy",
            "--in-seq",
            "r1.fq.gz",
            "--in-seq2",
            "r2.fq.gz",
            "--outprefix",
            "sample",
            "--identity",
            "0.97",
            "--fraction",
            "0.80",
            "--min-cbreadth",
            "0.10",
            "--lowbiomass",
            "--pathogen-host",
            "pathogens.tsv",
            "--keep-mapping-output",
            "--keep-raw",
        ])
        .unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        assert_eq!(args.dmp_dir.as_deref(), Some("taxonomy"));
        assert_eq!(args.in_seq.as_deref(), Some("r1.fq.gz"));
        assert_eq!(args.in_seq2.as_deref(), Some("r2.fq.gz"));
        assert_eq!(args.outprefix.as_deref(), Some("sample"));
        assert_eq!(args.identity, Some(0.97));
        assert_eq!(args.fraction, Some(0.80));
        assert_eq!(args.chunk_breadth, Some(0.10));
        assert_eq!(args.pathogen_host.as_deref(), Some("pathogens.tsv"));
        assert!(args.lowbiomass);
        assert!(args.keep_mapping_output);
        assert!(args.keep_raw);
    }

    #[test]
    fn profile_cli_accepts_explicit_staging_and_streaming_modes() {
        let staged = Cli::try_parse_from([
            "metax",
            "profile",
            "--reads",
            "reads.fq",
            "--output",
            "sample",
            "--temp-dir",
            "/scratch",
        ])
        .unwrap();
        let Commands::Profile(staged) = staged.command else {
            panic!("expected profile command");
        };
        assert_eq!(staged.temp_dir.as_deref(), Some("/scratch"));
        assert!(!staged.stream_mappings);

        let streamed = Cli::try_parse_from([
            "metax",
            "profile",
            "--reads",
            "reads.fq",
            "--output",
            "sample",
            "--stream-mappings",
        ])
        .unwrap();
        let Commands::Profile(streamed) = streamed.command else {
            panic!("expected profile command");
        };
        assert!(streamed.stream_mappings);
        assert!(streamed.temp_dir.is_none());
    }

    #[test]
    fn streaming_conflicts_with_saved_or_temporary_mapping_output() {
        for conflicting in ["--save-mapping", "--temp-dir"] {
            let mut args = vec![
                "metax",
                "profile",
                "--reads",
                "reads.fq",
                "--output",
                "sample",
                "--stream-mappings",
                conflicting,
            ];
            if conflicting == "--temp-dir" {
                args.push("/scratch");
            }
            assert!(Cli::try_parse_from(args).is_err());
        }
    }

    #[test]
    fn profile_help_is_grouped_by_pipeline_stage() {
        let help = long_help("profile");
        assert_appears_in_order(
            &help,
            &[
                "Input and output:",
                "Analysis mode and resources:",
                "Species and strain evidence:",
                "Genus fallback and low-map behavior:",
                "Coverage evidence and filters:",
                "Pathogen-focused analysis:",
                "Saved and reused data:",
                "Performance and robustness:",
                "Progress and logging:",
            ],
        );
        for canonical in [
            "--dmp-dir <DIR>",
            "--reads <READS>",
            "--reads2 <READS2>",
            "--min-identity <FRACTION>",
            "--min-mapped-fraction <FRACTION>",
            "--exact-ties <POLICY>",
            "--min-chunk-breadth <FRACTION>",
            "--chunk-breadth-scale <FACTOR>",
            "--reuse-mapping <MAPPING_TSV>",
            "--save-mapping",
            "--save-classification",
            "--stream-mappings",
            "--temp-dir <DIR>",
            "--mapper-help",
        ] {
            assert!(
                help.contains(canonical),
                "missing {canonical:?} from help\n{help}"
            );
        }
        assert!(!help.contains("--taxonomy-dir"));
        assert!(!help.contains("--keep-mapping-output"));
        assert!(!help.contains("Mapping-stage acceptance:"));
        assert!(help.contains("~/.metax/taxdmp"));
        assert!(help.contains("downloads a validated dump"));
        assert!(help.contains("does not discard lower-identity mapper output"));
        assert!(help.contains("independent of --min-identity"));
        assert!(help.contains("Advanced mapper candidate generation:"));
        assert_appears_in_order(
            &help,
            &[
                "--min-chunk-breadth <FRACTION>",
                "--chunk-breadth-scale <FACTOR>",
                "--by-mapped",
                "--min-reads <N>",
            ],
        );
    }

    #[test]
    fn mapper_help_is_available_without_profile_run_arguments() {
        let cli = Cli::try_parse_from(["metax", "profile", "--mapper-help"]).unwrap();
        let Commands::Profile(args) = cli.command else {
            panic!("expected profile command");
        };
        assert!(args.mapper_help);
        assert!(args.outprefix.is_none());
        assert!(args.extra_args.is_empty());

        let error = Cli::try_parse_from(["metax", "profile"]).unwrap_err();
        assert!(error.to_string().contains("--output <PREFIX>"));
    }

    #[test]
    fn mapper_help_documents_all_user_facing_mapper_overrides() {
        let help = mapper_options_help();
        for option in [
            "--top-n",
            "--score-ratio",
            "--min-seeds",
            "--allow-single-strobe-matches",
            "--match-both-strobes",
            "--query-seed-stride",
            "--min-identity",
            "--min-mapped-frac",
            "--no-identity-filter",
            "--ignore-n",
            "--pair-aware",
            "--paired-emission",
            "--min-insert",
            "--max-insert",
            "--identity-mode",
            "--extend",
            "--fast-extend",
            "--no-long-read-rescue",
            "--max-rescue-candidates",
            "--batch-size",
            "--index-cache",
            "--preset",
        ] {
            assert!(help.contains(option), "missing {option:?} from mapper help");
        }
        assert!(help.contains("top-n=50, score-ratio=0.60, min-seeds=1"));
        assert!(help.contains("mapped-length OR gate=300 bases"));
        assert!(help.contains("options before '--' configure taxonomic"));
    }

    #[test]
    fn index_cli_defaults_to_exact_randstrobe_pairs_with_opt_out() {
        let cli = Cli::try_parse_from(["metax", "index", "reference.fa", "-o", "db"]).unwrap();
        let Commands::Index(args) = cli.command else {
            panic!("expected index command");
        };
        assert!(resolve_index_match_both_strobes(
            &args.seed_mode,
            args.match_both_strobes,
            args.first_strobe_only,
        )
        .unwrap());

        let cli = Cli::try_parse_from([
            "metax",
            "index",
            "reference.fa",
            "-o",
            "db",
            "--first-strobe-only",
        ])
        .unwrap();
        let Commands::Index(args) = cli.command else {
            panic!("expected index command");
        };
        assert!(args.first_strobe_only);
        assert!(!resolve_index_match_both_strobes(
            &args.seed_mode,
            args.match_both_strobes,
            args.first_strobe_only,
        )
        .unwrap());
    }

    #[test]
    fn index_cli_accepts_canonical_names_and_legacy_aliases() {
        let cli = Cli::try_parse_from([
            "metax",
            "index",
            "reference.fa.gz",
            "--output",
            "db",
            "--threads",
            "8",
            "--temp-dir",
            "/local/scratch",
            "--seed-length",
            "17",
            "--seed-window",
            "21",
            "--fraction",
            "0.10",
            "--subsample-seed",
            "7",
            "--segment-length",
            "1000",
            "--min-segment-length",
            "100",
        ])
        .unwrap();
        let Commands::Index(args) = cli.command else {
            panic!("expected index command");
        };
        assert_eq!(args.fasta, "reference.fa.gz");
        assert_eq!(args.outdir, "db");
        assert_eq!(args.threads, Some(8));
        assert_eq!(args.temp_dir.as_deref(), Some("/local/scratch"));
        assert_eq!(args.k, 17);
        assert_eq!(args.w, 21);
        assert_eq!(args.fraction, Some(0.10));
        assert_eq!(args.seed, 7);
        assert_eq!(args.segment_length, 1000);
        assert_eq!(args.min_length, 100);

        let cli = Cli::try_parse_from([
            "metax",
            "index",
            "reference.fa",
            "--outdir",
            "legacy-db",
            "--seed",
            "9",
            "--min-length",
            "250",
        ])
        .unwrap();
        let Commands::Index(args) = cli.command else {
            panic!("expected index command");
        };
        assert_eq!(args.outdir, "legacy-db");
        assert_eq!(args.seed, 9);
        assert_eq!(args.min_length, 250);
    }

    #[test]
    fn index_help_is_grouped_by_build_stage() {
        let help = long_help("index");
        assert_appears_in_order(
            &help,
            &[
                "Input and output:",
                "Compute resources:",
                "Seed index configuration:",
                "Randstrobe matching:",
                "Optional fractional-reference stage:",
                "Progress and logging:",
            ],
        );
        for canonical in [
            "--output <DIR>",
            "--seed-length <BASES>",
            "--seed-window <N>",
            "--temp-dir <DIR>",
            "--subsample-seed <SEED>",
            "--min-segment-length <BASES>",
        ] {
            assert!(
                help.contains(canonical),
                "missing {canonical:?} from help\n{help}"
            );
        }
        assert!(!help.contains("--outdir"));
    }
}

pub mod aligner;
pub mod logging;
pub mod profiler;
pub mod taxonomy;

use aligner::Aligner;
use anyhow::{anyhow, bail, Context, Result};
use gzp::{deflate::Gzip, ZBuilder};
use indexmap::IndexMap;
use profiler::{Profiler, ProfilerConfig};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::cmp::Reverse;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use taxonomy::Taxonomy;

use crate::logging::init_logging;
use crate::profiler::ProfileMode;

pub struct MetaxApp {
    pub config: AppConfig,
}

// CLI-facing configuration that gets resolved into profiler settings.
#[derive(Debug, Clone)]
pub struct AppConfig {
    // Input database and taxonomy locations.
    pub db: Option<String>,
    pub dmp_dir: Option<String>,
    // Input reads (comma-separated for paired-end).
    pub input_sequences: Option<String>,
    pub outprefix: String,
    pub threads: usize,
    pub resume: bool,
    pub reuse_sam: Option<String>,
    pub sequencer: String,
    pub is_paired: bool,
    pub strain: bool,
    pub mode: ProfileMode,
    // Optional overrides for profiling thresholds.
    pub batch_size: Option<usize>,
    pub identity: Option<f64>,
    pub mapped_len: Option<usize>,
    pub breadth: Option<f64>,
    pub chunk_breadth: Option<f64>,
    pub min_reads: Option<usize>,
    pub min_oebr: Option<f64>,
    pub min_coebr: Option<f64>,
    pub fraction: Option<f64>,
    pub lowbiomass: bool,
    pub keep_raw: bool,
    /// Force aligned-read basis for the chunk-breadth estimate.
    ///
    /// `false` = auto (aligned basis auto-forced for low-map Illumina
    /// non-pathogen runs, else total-read basis).
    /// `true` = use aligned-read basis regardless of mapping rate, except
    /// community `--genus-fallback` on subsampled DBs estimates an unset
    /// chunk-breadth threshold from total reads.
    pub by_aligned: bool,
    /// Lower-bound identity for genus-fallback candidates (default 0.80).
    pub genus_identity: f64,
    /// Strict upper bound on mapping rate for the low-map path (default 0.30).
    pub low_map_rate_threshold: f64,
    /// Enable genus-level fallback regardless of mapping rate (community profiler).
    /// Does not by itself apply the 1.5x low-map chunk-basis scaling.
    pub genus_fallback: bool,
    /// Require BOTH `cov_prob` and `chunk_prob` to pass the cutoff when
    /// filtering the final taxa list (AND instead of OR).
    pub strict: bool,
    pub compress_sam: bool,
    pub pathogen_host: Option<String>,
    pub host: Option<String>,
    pub verbose: bool,
    pub very_verbose: bool,
    // Extra aligner args passed after `--`.
    pub extra_args: Vec<String>,
}

impl MetaxApp {
    pub fn run(&self) -> Result<()> {
        // Orchestrate alignment, profiling, and optional compression.
        init_logging("MAIN", Some(format!("{}.log", self.config.outprefix)))?;
        if self.config.very_verbose {
            let main_command = self.build_main_command();
            log::info!(target: "MAIN", "Command:\n{}", main_command);
        }
        // Prepare SAM path (reuse/resume) or run alignment.
        let (sam_path, newly_aligned) = self.prepare_alignment()?;
        if self.config.sequencer.eq_ignore_ascii_case("illumina") {
            log::info!(
                target: "MAIN",
                "Taxonomy profiling for Illumina short reads."
            );
        } else {
            log::info!(
                target: "MAIN",
                "Taxonomy profiling for {} reads.",
                self.config.sequencer
            );
        }
        let taxonomy = Taxonomy::load_from_dmp(self.config.dmp_dir.as_deref())?;
        let profiler_config = self.build_profiler_config(&sam_path)?;
        let profiler = Profiler::new(profiler_config, taxonomy)?;
        // Run the chosen profiler (community or strain).
        profiler.run()?;
        drop(profiler);
        self.compress_alignment_if_requested(&sam_path, newly_aligned)?;
        Ok(())
    }

    fn compress_alignment_if_requested(&self, sam_path: &str, newly_aligned: bool) -> Result<()> {
        // Only compress when requested and the file is present/uncompressed.
        if !self.config.compress_sam {
            return Ok(());
        }
        let sam_path_ref = Path::new(sam_path);
        if !sam_path_ref.exists() {
            return Ok(());
        }
        let already_compressed = sam_path_ref
            .extension()
            .map(|ext| ext == "gz")
            .unwrap_or(false);
        if already_compressed {
            if self.config.verbose || self.config.very_verbose {
                log::info!(
                    target: "MAIN",
                    "Compression requested but {} already appears compressed; skipping.",
                    sam_path_ref.display()
                );
            }
            return Ok(());
        }
        let threads = self.config.threads.max(1);
        if newly_aligned {
            log::info!(
                target: "MAIN",
                "Generating compressed alignment file ..."
            );
        } else {
            log::info!(
                target: "MAIN",
                "Compressing reused alignment file ..."
            );
        }
        let gz_path = compress_to_gz(sam_path_ref, threads)?;
        log::info!(
            target: "MAIN",
            "Compressed alignment to {}",
            gz_path.display()
        );
        Ok(())
    }

    fn prepare_alignment(&self) -> Result<(String, bool)> {
        // Respect explicit reuse and resume-before-aligning behavior.
        if let Some(ref reuse) = self.config.reuse_sam {
            log::info!(
                target: "MAIN",
                "Taxonomy profiling based on specified alignment file."
            );
            log::info!(
                target: "MAIN",
                "Reusing existing SAM file at {} for profiling.",
                reuse
            );
            return Ok((reuse.clone(), false));
        }
        let sam_path = format!("{}.sam", self.config.outprefix);
        if self.config.resume {
            let sam_path_ref = std::path::Path::new(&sam_path);
            if sam_path_ref.exists() {
                log::info!(
                    target: "MAIN",
                    "Resume taxonomy profiling with your last alignment file."
                );
                return Ok((sam_path, false));
            }
            let sam_gz_path = format!("{}.sam.gz", self.config.outprefix);
            if std::path::Path::new(&sam_gz_path).exists() {
                log::info!(
                    target: "MAIN",
                    "Resume taxonomy profiling with your last compressed alignment file."
                );
                return Ok((sam_gz_path, false));
            }
        }
        let input_sequences = self
            .config
            .input_sequences
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("Please specify the input reads files."))?;
        let db = self
            .config
            .db
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("Please specify the reference database file."))?;
        let sequences: Vec<String> = input_sequences
            .split(',')
            .map(|s| s.trim().to_string())
            .collect();
        // Invoke the aligner with preset + mode parameters.
        let aligner = Aligner::new(db.clone(), sequences, self.config.outprefix.clone());
        if self.config.verbose || self.config.very_verbose {
            log::info!(
                target: "MAIN",
                "Running maCMD alignment with {} threads against {}.",
                self.config.threads,
                db
            );
        } else {
            log::info!(
                target: "MAIN",
                "Running alignment with {} threads.",
                self.config.threads
            );
        }
        aligner.run(
            self.config.threads,
            &self.config.sequencer,
            self.config.is_paired,
            self.config.mode,
            self.config.very_verbose,
            &self.config.extra_args,
        )?;
        Ok((sam_path, true))
    }

    fn build_main_command(&self) -> String {
        fn fmt_option<T: ToString>(value: &Option<T>) -> String {
            value
                .as_ref()
                .map(|v| v.to_string())
                .unwrap_or_else(|| "None".to_string())
        }

        // Render boolean flags as CLI switches for logging/debugging.
        let bool_params = {
            let mut flags = Vec::new();
            if self.config.resume {
                flags.push("-r".to_string());
            }
            if self.config.is_paired {
                flags.push("--is-paired".to_string());
            }
            if self.config.strain {
                flags.push("--strain".to_string());
            }
            if self.config.lowbiomass {
                flags.push("--lowbiomass".to_string());
            }
            if self.config.keep_raw {
                flags.push("--keep-raw".to_string());
            }
            if self.config.by_aligned {
                flags.push("--by-aligned".to_string());
            }
            if self.config.strict {
                flags.push("--strict".to_string());
            }
            if self.config.genus_fallback {
                flags.push("--genus-fallback".to_string());
            }
            if self.config.compress_sam {
                flags.push("--compress-sam".to_string());
            }
            if self.config.verbose {
                flags.push("--verbose".to_string());
            }
            if self.config.very_verbose {
                flags.push("--very-verbose".to_string());
            }
            flags.join(" ")
        };

        let extra_args = if self.config.extra_args.is_empty() {
            String::new()
        } else {
            self.config.extra_args.join(" ")
        };

        // Render mode as CLI string for logging/debugging.
        let mode = match self.config.mode {
            ProfileMode::Recall => "recall",
            ProfileMode::Precision => "precision",
            ProfileMode::Default => "default",
        };

        format!(
            concat!(
                "metax profile --db {} \\\n",
                "    --dmp-dir {} \\\n",
                "    -i {} \\\n",
                "    -o {} \\\n",
                "    -t {} \\\n",
                "    --reuse-sam {} \\\n",
                "    --sequencer {} \\\n",
                "    --mode {} \\\n",
                "    --batch-size {} \\\n",
                "    --identity {} \\\n",
                "    -m {} \\\n",
                "    --min-breadth {} \\\n",
                "    --min-cbreadth {} \\\n",
                "    --min-reads {} \\\n",
                "    --min-oebr {} \\\n",
                "    --min-coebr {} \\\n",
                "    -f {} \\\n",
                "    --pathogen-host {} \\\n",
                "    --host {} \\\n",
                "    {} \\\n",
                "    {}"
            ),
            fmt_option(&self.config.db),
            fmt_option(&self.config.dmp_dir),
            fmt_option(&self.config.input_sequences),
            self.config.outprefix,
            self.config.threads,
            fmt_option(&self.config.reuse_sam),
            self.config.sequencer,
            mode,
            fmt_option(&self.config.batch_size),
            fmt_option(&self.config.identity),
            fmt_option(&self.config.mapped_len),
            fmt_option(&self.config.breadth),
            fmt_option(&self.config.chunk_breadth),
            fmt_option(&self.config.min_reads),
            fmt_option(&self.config.min_oebr),
            fmt_option(&self.config.min_coebr),
            fmt_option(&self.config.fraction),
            fmt_option(&self.config.pathogen_host),
            fmt_option(&self.config.host),
            bool_params,
            extra_args
        )
    }

    fn resolve_thresholds(&self) -> (usize, f64, usize, f64) {
        // Fill in default thresholds based on sequencer and profile mode.
        let sequencer = self.config.sequencer.as_str();
        let mut batch_size = self.config.batch_size.unwrap_or(5000);
        let mut identity = self.config.identity;
        let mut mapped_len = self.config.mapped_len;
        let mut fraction = self.config.fraction;

        match sequencer {
            "Illumina" => {
                if batch_size == 0 {
                    batch_size = 5000;
                }
                match self.config.mode {
                    ProfileMode::Recall => {
                        if self.config.strain {
                            identity.get_or_insert(0.90);
                            mapped_len.get_or_insert(40);
                            fraction.get_or_insert(0.5);
                        } else {
                            identity.get_or_insert(0.95);
                            mapped_len.get_or_insert(50);
                            fraction.get_or_insert(0.6);
                        }
                    }
                    ProfileMode::Precision => {
                        if self.config.strain {
                            identity.get_or_insert(0.96);
                            mapped_len.get_or_insert(60);
                            fraction.get_or_insert(0.7);
                        } else {
                            identity.get_or_insert(0.98);
                            mapped_len.get_or_insert(50);
                            fraction.get_or_insert(0.8);
                        }
                    }
                    ProfileMode::Default => {
                        // Short-read (Illumina) defaults are tuned
                        // separately for strain vs community profiling.
                        //
                        // Strain mode uses the strictest short-read defaults.
                        //
                        // Community mode uses slightly more permissive defaults.
                        if self.config.strain {
                            identity.get_or_insert(0.98);
                            mapped_len.get_or_insert(50);
                            fraction.get_or_insert(0.8);
                        } else {
                            identity.get_or_insert(0.97);
                            mapped_len.get_or_insert(50);
                            fraction.get_or_insert(0.7);
                        }
                    }
                }
            }
            _ => {
                if batch_size == 0 {
                    batch_size = 5000;
                }
                match self.config.mode {
                    ProfileMode::Recall => {
                        if !self.config.strain && (sequencer == "Nanopore" || sequencer == "PacBio")
                        {
                            identity.get_or_insert(0.85);
                            mapped_len.get_or_insert(250);
                            fraction.get_or_insert(0.5);
                        } else {
                            let default_identity = match sequencer {
                                "Nanopore" => 0.84,
                                "PacBio" => 0.85,
                                _ => 0.90,
                            };
                            identity.get_or_insert(default_identity);
                            mapped_len.get_or_insert(200);
                            fraction.get_or_insert(0.55);
                        }
                    }
                    ProfileMode::Precision => {
                        if !self.config.strain && (sequencer == "Nanopore" || sequencer == "PacBio")
                        {
                            identity.get_or_insert(0.90);
                            mapped_len.get_or_insert(250);
                            fraction.get_or_insert(0.6);
                        } else {
                            let default_identity =
                                if sequencer == "assembly" { 0.95 } else { 0.87 };
                            identity.get_or_insert(default_identity);
                            mapped_len.get_or_insert(300);
                            fraction.get_or_insert(0.7);
                        }
                    }
                    ProfileMode::Default => {
                        let default_identity = if sequencer == "assembly" { 0.93 } else { 0.86 };
                        identity.get_or_insert(default_identity);
                        mapped_len.get_or_insert(250);
                        fraction.get_or_insert(0.6);
                    }
                }
            }
        }

        (
            batch_size,
            identity.unwrap_or(0.0),
            mapped_len.unwrap_or(0),
            fraction.unwrap_or(0.0),
        )
    }

    fn build_profiler_config(&self, sam: &str) -> Result<ProfilerConfig> {
        let (batch_size, identity, mapped_len, fraction) = self.resolve_thresholds();
        // Reject invalid resolved thresholds early.
        if identity <= 0.0 || mapped_len == 0 || fraction <= 0.0 {
            return Err(anyhow::anyhow!(
                "Invalid thresholds resolved for profiler: identity {}, mapped length {}, fraction {}",
                identity,
                mapped_len,
                fraction
            ));
        }
        // NOTE: `--genus-identity` and `--low-map-rate-threshold` are
        // validated eagerly at CLI parse time (see main::parse_probability),
        // so any library caller bypassing the CLI is trusted to supply
        // values in [0.0, 1.0].
        Ok(ProfilerConfig {
            sam: sam.to_string(),
            sequencer: self.config.sequencer.clone(),
            batch_size,
            is_paired: self.config.is_paired,
            identity,
            mapped_len,
            breadth: self.config.breadth,
            chunk_breadth: self.config.chunk_breadth,
            min_reads: self.config.min_reads,
            min_oebr: self.config.min_oebr,
            min_coebr: self.config.min_coebr,
            fraction,
            lowbiomass: self.config.lowbiomass,
            keep_raw: self.config.keep_raw,
            by_aligned: self.config.by_aligned,
            genus_identity: self.config.genus_identity,
            low_map_rate_threshold: self.config.low_map_rate_threshold,
            genus_fallback: self.config.genus_fallback,
            strict: self.config.strict,
            pathogen_host: self.config.pathogen_host.clone(),
            host: self.config.host.clone(),
            threads: self.config.threads,
            verbose: self.config.verbose,
            very_verbose: self.config.very_verbose,
            outprefix: self.config.outprefix.clone(),
            mode: self.config.mode,
            dmp_dir: self.config.dmp_dir.clone(),
            strain: self.config.strain,
        })
    }
}

fn compress_to_gz(path: &Path, threads: usize) -> Result<PathBuf> {
    // Write gzip alongside the input and remove the original when done.
    let file_name = path
        .file_name()
        .ok_or_else(|| anyhow!("failed to determine file name for {}", path.display()))?;
    let gz_name = format!("{}.gz", file_name.to_string_lossy());
    let gz_path = path.with_file_name(gz_name);
    let input = File::open(path)
        .with_context(|| format!("failed to open {} for compression", path.display()))?;
    let mut reader = BufReader::new(input);
    let output = File::create(&gz_path)
        .with_context(|| format!("failed to create {}", gz_path.display()))?;
    let mut encoder = ZBuilder::<Gzip, _>::new()
        .num_threads(threads)
        .from_writer(output);
    std::io::copy(&mut reader, &mut encoder)
        .with_context(|| format!("failed to compress {}", path.display()))?;
    encoder
        .finish()
        .with_context(|| format!("failed to finalize compression for {}", path.display()))?;
    fs::remove_file(path)
        .with_context(|| format!("failed to remove {} after compression", path.display()))?;
    Ok(gz_path)
}

pub fn build_index(
    fasta: &str,
    outdir: &str,
    subsample: Option<f64>,
    threads: Option<usize>,
    seed: u64,
    segment_length: usize,
    min_length: usize,
    compress_subsampled: bool,
) -> Result<()> {
    // Build a maCMD index, optionally subsampling large genomes first.
    let outdir_path = Path::new(outdir);
    fs::create_dir_all(outdir_path)
        .with_context(|| format!("failed to create output directory {outdir}"))?;
    let fasta_path = Path::new(fasta);
    if !fasta_path.exists() {
        return Err(anyhow!("reference FASTA {} does not exist", fasta));
    }

    if segment_length == 0 {
        bail!("segment length must be greater than zero");
    }
    if min_length == 0 {
        bail!("minimum length must be greater than zero");
    }

    let (fasta_for_index, subsampled_path) = if let Some(fraction) = subsample {
        let params = SubsampleParams {
            fraction,
            seed,
            segment_length,
            min_length,
        };
        let path = create_subsampled_fasta(fasta_path, outdir_path, params, threads)?;
        (path.clone(), Some(path))
    } else {
        (fasta_path.to_path_buf(), None)
    };

    let fasta_arg = format!(
        "{},{},metax_db",
        fasta_for_index.to_string_lossy(),
        outdir_path.to_string_lossy()
    );
    let status = Command::new("maCMD")
        .arg("--Create_Index")
        .arg(&fasta_arg)
        .status()
        .with_context(|| {
            format!(
                "failed to execute maCMD command: maCMD --Create_Index {}",
                fasta_arg
            )
        })?;
    if !status.success() {
        bail!("maCMD command failed with status {}", status);
    }
    if compress_subsampled {
        if let Some(path) = subsampled_path {
            if path.exists() {
                log::info!(
                    target: "MAIN",
                    "Generating compressed subsampled FASTA for {} ...",
                    path.display()
                );
                let thread_count = threads.unwrap_or(1).max(1);
                let gz_path = compress_to_gz(&path, thread_count)?;
                log::info!(
                    target: "MAIN",
                    "Compressed subsampled FASTA to {}",
                    gz_path.display()
                );
            }
        } else {
            log::warn!(
                target: "MAIN",
                "Compression requested but no subsampled FASTA was generated; skipping."
            );
        }
    }
    Ok(())
}

struct ContigRecord {
    fields: Vec<String>,
    sequence: String,
}

struct GenomeRecord {
    genome_size: usize,
    contigs: Vec<ContigRecord>,
}

#[derive(Clone, Copy)]
struct SubsampleParams {
    fraction: f64,
    seed: u64,
    segment_length: usize,
    min_length: usize,
}

#[derive(Clone)]
struct Segment {
    contig_index: usize,
    start: usize,
    end: usize,
}

impl Segment {
    fn len(&self) -> usize {
        self.end.saturating_sub(self.start).saturating_add(1)
    }
}

fn create_subsampled_fasta(
    fasta: &Path,
    outdir: &Path,
    params: SubsampleParams,
    threads: Option<usize>,
) -> Result<PathBuf> {
    // Validate subsampling parameters before doing any work.
    if !(0.0 < params.fraction && params.fraction < 1.0) {
        bail!("subsample fraction must be between 0 and 1 (exclusive)");
    }
    if let Some(t) = threads {
        if t == 0 {
            bail!("thread count must be at least 1");
        }
    }
    let genomes = parse_fasta(fasta)?;
    let output_path = outdir.join("metax_subsampled.fasta");
    let subsampled = generate_subsamples(&genomes, params, threads)?;
    let mut writer = BufWriter::new(
        File::create(&output_path)
            .with_context(|| format!("failed to create {}", output_path.display()))?,
    );
    for records in subsampled {
        for record in records {
            writer.write_all(record.header.as_bytes())?;
            writer.write_all(b"\n")?;
            for chunk in record.sequence.as_bytes().chunks(60) {
                writer.write_all(chunk)?;
                writer.write_all(b"\n")?;
            }
        }
    }
    writer.flush()?;
    Ok(output_path)
}

fn generate_subsamples(
    genomes: &[GenomeRecord],
    params: SubsampleParams,
    threads: Option<usize>,
) -> Result<Vec<Vec<SubsampledRecord>>> {
    // Generate subsampled contig segments per genome (parallel when requested).
    if genomes.is_empty() {
        return Ok(Vec::new());
    }

    let mut seed_rng = StdRng::seed_from_u64(params.seed);
    let genome_seeds: Vec<u64> = genomes.iter().map(|_| seed_rng.gen()).collect();

    let sequential = || -> Result<Vec<Vec<SubsampledRecord>>> {
        genomes
            .iter()
            .zip(genome_seeds.iter())
            .map(|(genome, &seed)| {
                let mut rng = StdRng::seed_from_u64(seed);
                subsample_genome(genome, params, &mut rng)
            })
            .collect()
    };

    if let Some(thread_count) = threads {
        if thread_count > 1 {
            let pool = ThreadPoolBuilder::new()
                .num_threads(thread_count)
                .build()
                .with_context(|| {
                    format!("failed to build thread pool with {} threads", thread_count)
                })?;
            let mut results = pool.install(|| -> Result<Vec<(usize, Vec<SubsampledRecord>)>> {
                (0..genomes.len())
                    .into_par_iter()
                    .map(|idx| {
                        let mut rng = StdRng::seed_from_u64(genome_seeds[idx]);
                        subsample_genome(&genomes[idx], params, &mut rng)
                            .map(|records| (idx, records))
                    })
                    .collect()
            })?;
            results.sort_by_key(|(idx, _)| *idx);
            Ok(results.into_iter().map(|(_, records)| records).collect())
        } else {
            sequential()
        }
    } else {
        sequential()
    }
}

fn parse_fasta(path: &Path) -> Result<Vec<GenomeRecord>> {
    // Read a multi-FASTA into per-genome records keyed by header fields.
    let file =
        File::open(path).with_context(|| format!("failed to open FASTA {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut genomes: IndexMap<String, GenomeRecord> = IndexMap::new();
    let mut current_header: Option<String> = None;
    let mut sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if let Some(header) = current_header.take() {
                store_record(&mut genomes, &header, &sequence)?;
                sequence.clear();
            }
            current_header = Some(line[1..].trim().to_string());
        } else {
            sequence.push_str(line.trim());
        }
    }
    if let Some(header) = current_header {
        store_record(&mut genomes, &header, &sequence)?;
    }
    Ok(genomes.into_values().collect())
}

fn store_record(
    genomes: &mut IndexMap<String, GenomeRecord>,
    header: &str,
    sequence: &str,
) -> Result<()> {
    // Parse header fields and append contig sequence to its genome.
    if header.is_empty() {
        return Ok(());
    }
    let mut fields: Vec<String> = header.split('|').map(|s| s.to_string()).collect();
    if fields.len() < 5 {
        return Err(anyhow!(
            "header '{}' does not contain at least five fields",
            header
        ));
    }
    let genome_id = fields[0].clone();
    let genome_size: usize = fields[4]
        .parse()
        .with_context(|| format!("invalid genome size in header '{}'", header))?;
    fields.truncate(5);
    let entry = genomes.entry(genome_id).or_insert_with(|| GenomeRecord {
        genome_size,
        contigs: Vec::new(),
    });
    if entry.genome_size != genome_size {
        entry.genome_size = genome_size;
    }
    entry.contigs.push(ContigRecord {
        fields,
        sequence: sequence.to_string(),
    });
    Ok(())
}

struct SubsampledRecord {
    header: String,
    sequence: String,
}

fn subsample_genome<R: rand::Rng + ?Sized>(
    genome: &GenomeRecord,
    params: SubsampleParams,
    rng: &mut R,
) -> Result<Vec<SubsampledRecord>> {
    // Sample segments to reach target fraction while preserving genome size metadata.
    const SHORT_GENOME_THRESHOLD: usize = 500_000;

    if genome.contigs.is_empty() {
        return Ok(Vec::new());
    }

    let mut segments: Vec<Segment> = Vec::new();
    if genome.genome_size < SHORT_GENOME_THRESHOLD {
        for (idx, contig) in genome.contigs.iter().enumerate() {
            let len = contig.sequence.len();
            if len == 0 {
                continue;
            }
            segments.push(Segment {
                contig_index: idx,
                start: 1,
                end: len,
            });
        }
    } else {
        let target_len = ((params.fraction * genome.genome_size as f64).ceil() as usize).max(1);
        let mut contig_candidates: Vec<(usize, Vec<Segment>)> = Vec::new();
        let mut total_slots = 0usize;

        for (idx, contig) in genome.contigs.iter().enumerate() {
            let len = contig.sequence.len();
            if len < params.segment_length {
                continue;
            }
            let max_segments = len / params.segment_length;
            if max_segments == 0 {
                continue;
            }
            let partition_size = len / max_segments;
            let mut contig_segments: Vec<Segment> = Vec::new();
            let mut partition_start = 0usize;

            for partition in 0..max_segments {
                let mut partition_end = if partition == max_segments - 1 {
                    len
                } else {
                    (partition + 1) * partition_size
                };
                if partition_end > len {
                    partition_end = len;
                }
                let part_start = partition_start.min(len);
                partition_start = partition_end;
                if partition_end <= part_start {
                    continue;
                }
                let available = partition_end.saturating_sub(part_start);
                if available < params.segment_length {
                    continue;
                }
                let max_start = partition_end - params.segment_length;
                if max_start < part_start {
                    continue;
                }
                let start = rng.gen_range(part_start..=max_start) + 1;
                let end = start + params.segment_length - 1;
                contig_segments.push(Segment {
                    contig_index: idx,
                    start,
                    end,
                });
            }

            if !contig_segments.is_empty() {
                contig_segments.sort_by_key(|seg| seg.start);
                total_slots += contig_segments.len();
                contig_candidates.push((idx, contig_segments));
            }
        }

        if total_slots == 0 || total_slots * params.segment_length < target_len {
            segments = select_longest_contigs(genome, target_len, params);
        } else {
            let target_segments = (target_len + params.segment_length - 1) / params.segment_length;
            let mut allocations = vec![0usize; contig_candidates.len()];
            if target_segments >= total_slots {
                for (i, (_, segs)) in contig_candidates.iter().enumerate() {
                    allocations[i] = segs.len();
                }
            } else {
                let total_slots_f = total_slots as f64;
                for (i, (_, segs)) in contig_candidates.iter().enumerate() {
                    let share = segs.len() as f64 / total_slots_f;
                    let allocated = (share * target_segments as f64).floor() as usize;
                    allocations[i] = allocated.min(segs.len());
                }
                let used = allocations.iter().sum::<usize>();
                let mut remainder = target_segments.saturating_sub(used);
                if remainder > 0 {
                    let mut fractional: Vec<(usize, f64)> = contig_candidates
                        .iter()
                        .enumerate()
                        .map(|(i, (_, segs))| {
                            let share = segs.len() as f64 / total_slots_f;
                            let exact = share * target_segments as f64;
                            let frac = exact - allocations[i] as f64;
                            (i, frac)
                        })
                        .collect();
                    fractional
                        .sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
                    for (idx, _) in fractional.iter() {
                        if remainder == 0 {
                            break;
                        }
                        if allocations[*idx] < contig_candidates[*idx].1.len() {
                            allocations[*idx] += 1;
                            remainder -= 1;
                        }
                    }
                    if remainder > 0 {
                        for (i, (_, segs)) in contig_candidates.iter().enumerate() {
                            if remainder == 0 {
                                break;
                            }
                            while allocations[i] < segs.len() && remainder > 0 {
                                allocations[i] += 1;
                                remainder -= 1;
                            }
                        }
                    }
                }
            }

            for (allocation, (_, segs)) in allocations.into_iter().zip(contig_candidates.iter()) {
                if allocation == 0 {
                    continue;
                }
                let total = segs.len();
                if allocation >= total {
                    segments.extend(segs.iter().cloned());
                } else {
                    for j in 0..allocation {
                        let idx = ((j + 1) * total) / allocation - 1;
                        segments.push(segs[idx].clone());
                    }
                }
            }
        }
    }

    if segments.is_empty() {
        return Ok(Vec::new());
    }

    segments.sort_by_key(|seg| (seg.contig_index, seg.start));
    let segments = enforce_min_segment_length(segments, genome, params.min_length);
    if segments.is_empty() {
        return Ok(Vec::new());
    }

    let total_len: usize = segments.iter().map(|s| s.len()).sum();
    if total_len == 0 {
        return Ok(Vec::new());
    }

    let mut records: Vec<SubsampledRecord> = Vec::new();
    for segment in segments {
        let contig = &genome.contigs[segment.contig_index];
        if contig.sequence.is_empty() {
            continue;
        }
        let mut fields = contig.fields.clone();
        let accession = fields
            .get(3)
            .cloned()
            .unwrap_or_else(|| String::from("unknown"));
        fields[3] = format!("{}:{}-{}", accession, segment.start, segment.end);
        fields[4] = genome.genome_size.to_string();
        if fields.len() > 5 {
            fields.truncate(5);
        }
        fields.push(total_len.to_string());
        let header = format!(">{}", fields.join("|"));
        let start = segment.start.saturating_sub(1);
        let end = segment.end.min(contig.sequence.len());
        if start >= end {
            continue;
        }
        let sequence = contig.sequence[start..end].to_string();
        records.push(SubsampledRecord { header, sequence });
    }
    Ok(records)
}

fn select_longest_contigs(
    genome: &GenomeRecord,
    target_len: usize,
    params: SubsampleParams,
) -> Vec<Segment> {
    // Fallback: take longest contigs until target length is reached.
    let mut indices: Vec<usize> = (0..genome.contigs.len()).collect();
    indices.sort_by_key(|&idx| Reverse(genome.contigs[idx].sequence.len()));

    let mut segments: Vec<Segment> = Vec::new();
    let mut fallback: Vec<Segment> = Vec::new();
    let apply_filter = genome.genome_size > params.min_length.saturating_mul(10);
    let mut accumulated = 0usize;

    for idx in indices {
        if accumulated >= target_len {
            break;
        }
        let contig = &genome.contigs[idx];
        let len = contig.sequence.len();
        if len == 0 {
            continue;
        }
        let remaining = target_len.saturating_sub(accumulated);
        let end = if len > remaining && remaining > 0 {
            remaining
        } else {
            len
        };
        if end == 0 {
            continue;
        }
        let segment = Segment {
            contig_index: idx,
            start: 1,
            end,
        };
        if apply_filter && segment.len() < params.min_length {
            fallback.push(segment);
        } else {
            accumulated += segment.len();
            segments.push(segment);
        }
    }

    if accumulated < target_len {
        for segment in fallback.iter() {
            if accumulated >= target_len {
                break;
            }
            accumulated += segment.len();
            segments.push(segment.clone());
        }
    }

    if segments.is_empty() {
        segments = fallback;
    }

    segments
}

fn enforce_min_segment_length(
    mut segments: Vec<Segment>,
    genome: &GenomeRecord,
    min_length: usize,
) -> Vec<Segment> {
    // Drop short segments for large genomes to avoid noisy fragments.
    if segments.is_empty() {
        return segments;
    }

    if min_length == 0 {
        segments.sort_by_key(|seg| (seg.contig_index, seg.start));
        return segments;
    }

    let threshold = min_length.saturating_mul(10);
    if genome.genome_size <= threshold {
        segments.sort_by_key(|seg| (seg.contig_index, seg.start));
        return segments;
    }

    let mut filtered: Vec<Segment> = segments
        .iter()
        .cloned()
        .filter(|segment| segment.len() >= min_length)
        .collect();

    if filtered.is_empty() {
        segments.sort_by_key(|seg| (seg.contig_index, seg.start));
        return segments;
    }

    filtered.sort_by_key(|seg| (seg.contig_index, seg.start));
    filtered
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Construct a minimal `AppConfig` suitable for exercising
    /// `resolve_thresholds`. Only the fields consulted by that function
    /// matter; everything else gets plausible defaults. Tests override
    /// the specific fields they care about.
    fn test_app_config() -> AppConfig {
        AppConfig {
            db: None,
            dmp_dir: None,
            input_sequences: None,
            outprefix: "test".to_string(),
            threads: 1,
            resume: false,
            reuse_sam: None,
            sequencer: "Illumina".to_string(),
            is_paired: false,
            strain: false,
            mode: ProfileMode::Default,
            batch_size: None,
            identity: None,
            mapped_len: None,
            breadth: None,
            chunk_breadth: None,
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            fraction: None,
            lowbiomass: false,
            keep_raw: false,
            by_aligned: false,
            genus_identity: 0.80,
            low_map_rate_threshold: 0.30,
            genus_fallback: false,
            strict: false,
            compress_sam: false,
            pathogen_host: None,
            host: None,
            verbose: false,
            very_verbose: false,
            extra_args: Vec::new(),
        }
    }

    fn resolve(cfg: AppConfig) -> (usize, f64, usize, f64) {
        MetaxApp { config: cfg }.resolve_thresholds()
    }

    // ------------------------------------------------------------------
    // Illumina + non-strain + Default mode => identity 0.97 / fraction 0.7
    // ------------------------------------------------------------------

    #[test]
    fn illumina_default_non_strain_uses_new_community_defaults() {
        let mut cfg = test_app_config();
        cfg.sequencer = "Illumina".to_string();
        cfg.strain = false;
        cfg.mode = ProfileMode::Default;
        let (_, identity, mapped_len, fraction) = resolve(cfg);
        assert!((identity - 0.97).abs() < 1e-12);
        assert_eq!(mapped_len, 50);
        assert!((fraction - 0.70).abs() < 1e-12);
    }

    #[test]
    fn illumina_default_paired_non_strain_also_uses_new_community_defaults() {
        // Paired-end must follow the same non-strain branch as single-end.
        let mut cfg = test_app_config();
        cfg.sequencer = "Illumina".to_string();
        cfg.is_paired = true;
        cfg.strain = false;
        cfg.mode = ProfileMode::Default;
        let (_, identity, _, fraction) = resolve(cfg);
        assert!((identity - 0.97).abs() < 1e-12);
        assert!((fraction - 0.70).abs() < 1e-12);
    }

    #[test]
    fn illumina_default_strain_uses_strain_specific_defaults() {
        // Strain mode on the Default preset uses the tightest
        // thresholds: identity 0.98, fraction 0.8 (mapped_len 50).
        let mut cfg = test_app_config();
        cfg.sequencer = "Illumina".to_string();
        cfg.strain = true;
        cfg.mode = ProfileMode::Default;
        let (_, identity, mapped_len, fraction) = resolve(cfg);
        assert!((identity - 0.98).abs() < 1e-12);
        assert_eq!(mapped_len, 50);
        assert!((fraction - 0.80).abs() < 1e-12);
    }

    #[test]
    fn illumina_default_strain_paired_also_uses_strain_specific_defaults() {
        // Paired-end short-read strain mode must match single-end.
        let mut cfg = test_app_config();
        cfg.sequencer = "Illumina".to_string();
        cfg.is_paired = true;
        cfg.strain = true;
        cfg.mode = ProfileMode::Default;
        let (_, identity, _, fraction) = resolve(cfg);
        assert!((identity - 0.98).abs() < 1e-12);
        assert!((fraction - 0.80).abs() < 1e-12);
    }

    #[test]
    fn strain_defaults_apply_only_to_illumina() {
        // Non-Illumina + strain must keep the long-read defaults; the
        // 0.98 / 0.8 tightening is Illumina-only.
        let mut cfg = test_app_config();
        cfg.sequencer = "Nanopore".to_string();
        cfg.strain = true;
        cfg.mode = ProfileMode::Default;
        let (_, identity, mapped_len, _) = resolve(cfg);
        assert_eq!(mapped_len, 250);
        assert!(
            identity < 0.90,
            "Nanopore strain identity should still be ~0.86, got {identity}"
        );
    }

    #[test]
    fn illumina_community_precision_and_recall_use_community_defaults() {
        let mut cfg = test_app_config();
        cfg.sequencer = "Illumina".to_string();
        cfg.strain = false;
        cfg.mode = ProfileMode::Precision;
        let (_, identity, mapped_len, fraction) = resolve(cfg.clone());
        assert!((identity - 0.98).abs() < 1e-12);
        assert_eq!(mapped_len, 50);
        assert!((fraction - 0.80).abs() < 1e-12);

        let mut cfg_recall = cfg;
        cfg_recall.mode = ProfileMode::Recall;
        let (_, identity, mapped_len, fraction) = resolve(cfg_recall);
        assert!((identity - 0.95).abs() < 1e-12);
        assert_eq!(mapped_len, 50);
        assert!((fraction - 0.60).abs() < 1e-12);
    }

    #[test]
    fn illumina_strain_precision_and_recall_keep_existing_defaults() {
        let mut cfg = test_app_config();
        cfg.sequencer = "Illumina".to_string();
        cfg.strain = true;
        cfg.mode = ProfileMode::Precision;
        let (_, identity, mapped_len, fraction) = resolve(cfg.clone());
        assert!((identity - 0.96).abs() < 1e-12);
        assert_eq!(mapped_len, 60);
        assert!((fraction - 0.70).abs() < 1e-12);

        let mut cfg_recall = cfg;
        cfg_recall.mode = ProfileMode::Recall;
        let (_, identity, mapped_len, fraction) = resolve(cfg_recall);
        assert!((identity - 0.90).abs() < 1e-12);
        assert_eq!(mapped_len, 40);
        assert!((fraction - 0.50).abs() < 1e-12);
    }

    #[test]
    fn user_overrides_beat_any_default() {
        // Explicit CLI values must pass through untouched regardless of
        // the Illumina-non-strain-Default tightening.
        let mut cfg = test_app_config();
        cfg.sequencer = "Illumina".to_string();
        cfg.strain = false;
        cfg.mode = ProfileMode::Default;
        cfg.identity = Some(0.80);
        cfg.fraction = Some(0.50);
        cfg.mapped_len = Some(123);
        let (_, identity, mapped_len, fraction) = resolve(cfg);
        assert!((identity - 0.80).abs() < 1e-12);
        assert_eq!(mapped_len, 123);
        assert!((fraction - 0.50).abs() < 1e-12);
    }

    #[test]
    fn non_illumina_defaults_unchanged() {
        // Nanopore / PacBio must not be affected by the Illumina-only
        // non-strain tightening.
        for seq in ["Nanopore", "PacBio"] {
            let mut cfg = test_app_config();
            cfg.sequencer = seq.to_string();
            cfg.strain = false;
            cfg.mode = ProfileMode::Default;
            let (_, identity, mapped_len, _) = resolve(cfg);
            // Long-read defaults: identity=0.86/0.85, mapped_len=250.
            // We only check mapped_len — identity differs per sequencer.
            assert_eq!(mapped_len, 250, "{}", seq);
            assert!(identity < 0.90);
        }
    }

    #[test]
    fn long_read_community_precision_and_recall_use_community_defaults() {
        for seq in ["Nanopore", "PacBio"] {
            let mut cfg = test_app_config();
            cfg.sequencer = seq.to_string();
            cfg.strain = false;
            cfg.mode = ProfileMode::Precision;
            let (_, identity, mapped_len, fraction) = resolve(cfg.clone());
            assert!((identity - 0.90).abs() < 1e-12, "{}", seq);
            assert_eq!(mapped_len, 250, "{}", seq);
            assert!((fraction - 0.60).abs() < 1e-12, "{}", seq);

            let mut cfg_recall = cfg;
            cfg_recall.mode = ProfileMode::Recall;
            let (_, identity, mapped_len, fraction) = resolve(cfg_recall);
            assert!((identity - 0.85).abs() < 1e-12, "{}", seq);
            assert_eq!(mapped_len, 250, "{}", seq);
            assert!((fraction - 0.50).abs() < 1e-12, "{}", seq);
        }
    }

    #[test]
    fn long_read_strain_precision_and_recall_keep_existing_defaults() {
        for (seq, recall_identity) in [("Nanopore", 0.84), ("PacBio", 0.85)] {
            let mut cfg = test_app_config();
            cfg.sequencer = seq.to_string();
            cfg.strain = true;
            cfg.mode = ProfileMode::Precision;
            let (_, identity, mapped_len, fraction) = resolve(cfg.clone());
            assert!((identity - 0.87).abs() < 1e-12, "{}", seq);
            assert_eq!(mapped_len, 300, "{}", seq);
            assert!((fraction - 0.70).abs() < 1e-12, "{}", seq);

            let mut cfg_recall = cfg;
            cfg_recall.mode = ProfileMode::Recall;
            let (_, identity, mapped_len, fraction) = resolve(cfg_recall);
            assert!((identity - recall_identity).abs() < 1e-12, "{}", seq);
            assert_eq!(mapped_len, 200, "{}", seq);
            assert!((fraction - 0.55).abs() < 1e-12, "{}", seq);
        }
    }
}

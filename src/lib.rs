pub mod db;
pub mod info;
pub mod logging;
pub mod mapper;
pub mod profiler;
pub mod progress;
pub mod sequencer;
pub mod taxonomy;

use std::path::Path;
use std::sync::Arc;

use anyhow::{anyhow, Result};
use mapper::{
    existing_mapping_output, is_mapping_consumer_closed, Mapper, MappingOutput, MappingProgressLog,
    StreamingStructuredMappingOutput,
};
use profiler::{ExactTiePolicy, Profiler, ProfilerConfig};
use taxonomy::Taxonomy;

use crate::db::{build_database, DatabaseMetadata, MetamIndexOptions};
use crate::logging::init_logging;
use crate::mapper::staging::{human_bytes, is_storage_exhausted, plan_staging_directory};
use crate::mapper::tsv::{
    load_mapping_run_metadata, new_shared_metam_profile, validate_mapping_database_summary,
    MappingProfilerMetadata, MappingRunMetadata, SharedMetamProfile,
};
use crate::profiler::ProfileMode;
use crate::progress::ProgressOptions;
use crate::sequencer::Sequencer;

pub struct MetaxApp {
    pub config: AppConfig,
}

enum PreparedMapping {
    Sam {
        path: String,
    },
    Metam {
        output: MappingOutput,
        db: String,
    },
    MetamStream {
        mapper: Mapper,
        profile: SharedMetamProfile,
        db: String,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct ResolvedProfilerThresholds {
    batch_size: usize,
    identity: f64,
    mapped_len: usize,
    fraction: f64,
    sampled_mean_base_accuracy: Option<f64>,
}

impl ResolvedProfilerThresholds {
    fn mapping_metadata(self, strain: bool) -> MappingProfilerMetadata {
        MappingProfilerMetadata {
            strain,
            identity: self.identity,
            mapped_len: self.mapped_len,
            fraction: self.fraction,
            sampled_mean_base_accuracy: self.sampled_mean_base_accuracy,
        }
    }
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
    pub reuse_mappings: Option<String>,
    pub sequencer: Sequencer,
    pub is_paired: bool,
    pub exact_ties: ExactTiePolicy,
    pub strain: bool,
    pub mode: ProfileMode,
    // Optional overrides for profiling thresholds.
    pub batch_size: Option<usize>,
    pub max_em_iterations: usize,
    pub lenient_mapping_input: bool,
    pub identity: Option<f64>,
    pub mapped_len: Option<usize>,
    pub breadth: Option<f64>,
    pub chunk_breadth: Option<f64>,
    /// Multiplicative adjustment for an automatically estimated chunk-breadth
    /// cutoff. A value of 1.0 preserves the calibrated policy exactly.
    pub chunk_breadth_scale: f64,
    pub min_reads: Option<usize>,
    pub min_oebr: Option<f64>,
    pub min_coebr: Option<f64>,
    pub fraction: Option<f64>,
    pub lowbiomass: bool,
    pub keep_raw: bool,
    /// Write the optional per-read classification table.
    pub save_classifications: bool,
    /// Force mapped-read basis for the chunk-breadth estimate.
    ///
    /// `false` = auto (mapped basis auto-forced for low-map Illumina
    /// non-pathogen runs, else total-read basis).
    /// `true` = use mapped-read basis regardless of mapping rate, except
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
    pub keep_mapping_output: bool,
    /// Use the direct mapper-to-profiler channel instead of disk staging.
    pub stream_mappings: bool,
    /// Explicit directory for temporary staged mappings.
    pub temp_dir: Option<String>,
    pub pathogen_host: Option<String>,
    pub host: Option<String>,
    pub verbose: bool,
    pub very_verbose: bool,
    pub progress: bool,
    pub no_progress: bool,
    pub progress_interval: u64,
    // Extra mapper args passed after `--`.
    pub extra_args: Vec<String>,
}

/// Head reads sampled to estimate long-read base quality for the adaptive
/// identity cutoff (fewer if the sample is smaller).
const QUALITY_SAMPLE_READS: usize = 10_000;
/// Safety cap on quality bytes tallied while sampling, so a few pathologically
/// long reads cannot stall startup.
const QUALITY_SAMPLE_MAX_BASES: u64 = 250_000_000;
/// Mean predicted per-base accuracy at or above which long reads keep the
/// fixed preset identity cutoff. Below it, the cutoff self-calibrates.
const ADAPTIVE_IDENTITY_CEILING: f64 = 0.95;
/// Hard lower bound for the self-calibrated long-read identity cutoff.
const ADAPTIVE_IDENTITY_FLOOR: f64 = 0.84;

/// Map a sampled mean predicted base accuracy to a self-calibrated long-read
/// identity cutoff, or `None` to keep the preset cutoff.
///
/// High-quality samples (accuracy >= [`ADAPTIVE_IDENTITY_CEILING`]) keep the
/// preset. Otherwise the cutoff is `accuracy - margin`, where the margin is
/// wider for recall and narrower for precision, and never below
/// [`ADAPTIVE_IDENTITY_FLOOR`]. Note that base-quality accuracy reflects
/// sequencing error only, not reference divergence, so a low-quality sample's
/// reads align a little below their base quality; the margin absorbs that.
fn adaptive_identity_from_accuracy(accuracy: f64, mode: ProfileMode) -> Option<f64> {
    if !accuracy.is_finite() || accuracy >= ADAPTIVE_IDENTITY_CEILING {
        return None;
    }
    let margin = match mode {
        ProfileMode::Recall => 0.04,
        ProfileMode::Default => 0.03,
        ProfileMode::Precision => 0.02,
    };
    Some((accuracy - margin).max(ADAPTIVE_IDENTITY_FLOOR))
}

impl MetaxApp {
    pub fn run(&self) -> Result<()> {
        // Orchestrate mapping, profiling, and optional compression.
        self.validate_chunk_breadth_configuration()?;
        init_logging("MAIN", Some(format!("{}.log", self.config.outprefix)))?;
        if self.config.very_verbose {
            let main_command = self.build_main_command();
            log::info!(target: "MAIN", "Command:\n{}", main_command);
        }
        // Resolve or install the taxonomy dump before starting an expensive
        // mapping run. The taxonomy itself is loaded only after mapping so it
        // does not increase mapper peak memory.
        let taxonomy_source = Taxonomy::prepare_dmp(self.config.dmp_dir.as_deref())?;
        // Resolve read-level profiler gates before mapping so saved mapping
        // metadata records the exact values used by this run. Resume/reuse
        // loads those values first, making profiling reproducible without the
        // original FASTQ quality strings.
        let mapping_metadata = self.mapping_metadata_for_thresholds()?;
        let thresholds = self.resolve_thresholds(mapping_metadata.as_ref());
        let prepared_mapping =
            self.prepare_mapping(thresholds.mapping_metadata(self.config.strain))?;
        if self.config.sequencer.is_illumina() {
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
        let taxonomy = taxonomy_source.load()?;
        let profiler_config = self.build_profiler_config(&prepared_mapping, thresholds)?;
        let profiler = Profiler::new(profiler_config, taxonomy)?;
        match prepared_mapping {
            PreparedMapping::MetamStream {
                mapper, profile, ..
            } => {
                self.run_structured_mapping(profiler, mapper, profile)?;
            }
            PreparedMapping::Sam { .. } => {
                profiler.run()?;
                drop(profiler);
            }
            PreparedMapping::Metam { output, .. } => {
                let profile_result = profiler.run();
                drop(profiler);
                if let Err(profile_error) = profile_result {
                    match output.preserve_after_failure() {
                        Ok(Some(path)) => log::warn!(
                            target: "MAPPER",
                            "Profiling failed; preserving completed temporary mapping output at {}.",
                            path.display()
                        ),
                        Ok(None) => {}
                        Err(preserve_error) => log::warn!(
                            target: "MAPPER",
                            "Profiling failed and the temporary mapping output could not be preserved: {:#}",
                            preserve_error
                        ),
                    }
                    return Err(profile_error);
                }
                output.cleanup_after_success()?;
            }
        }
        Ok(())
    }

    fn validate_chunk_breadth_configuration(&self) -> Result<()> {
        let factor = self.config.chunk_breadth_scale;
        if !factor.is_finite() || factor <= 0.0 {
            return Err(anyhow!(
                "chunk-breadth scale must be a finite value greater than zero, got {factor}"
            ));
        }
        if factor != 1.0 {
            if self.config.chunk_breadth.is_some() {
                return Err(anyhow!(
                    "chunk-breadth scale cannot be combined with an explicit minimum chunk breadth"
                ));
            }
            if self.config.min_reads.is_some() {
                return Err(anyhow!(
                    "chunk-breadth scale cannot be combined with a minimum-read cutoff, which disables automatic chunk-breadth estimation"
                ));
            }
            if self.config.lowbiomass {
                return Err(anyhow!(
                    "chunk-breadth scale cannot be combined with low-biomass mode, which disables automatic chunk-breadth estimation"
                ));
            }
        }
        Ok(())
    }

    fn run_structured_mapping(
        &self,
        profiler: Profiler,
        mapper: Mapper,
        profile: SharedMetamProfile,
    ) -> Result<()> {
        let stream = mapper.spawn_structured_stream(
            self.config.threads,
            self.config.sequencer,
            self.config.is_paired,
            self.config.mode,
            self.config.very_verbose,
            &self.config.extra_args,
            self.profile_progress_options(),
            MappingProgressLog::streamed_default(!self.config.no_progress),
            profile,
        )?;
        let StreamingStructuredMappingOutput {
            receiver,
            profile: _profile,
            worker,
        } = stream;
        let profile_result =
            profiler.run_with_metam_mappings(receiver, "in-process structured mapping stream");
        let worker_result = worker
            .join()
            .map_err(|_| anyhow!("in-process mapper worker panicked"))?;
        match (profile_result, worker_result) {
            (Ok(()), Ok(())) => {
                log::info!(target: "MAPPER", "Structured streaming mapping finished.");
                Ok(())
            }
            (Ok(()), Err(worker_err)) => Err(worker_err),
            (Err(profile_err), Ok(())) => Err(profile_err),
            (Err(profile_err), Err(worker_err)) if is_mapping_consumer_closed(&worker_err) => {
                Err(anyhow!(
                    "profiler failed while consuming mappings: {profile_err:#}\nmapper then stopped: {worker_err:#}"
                ))
            }
            (Err(profile_err), Err(worker_err)) => Err(anyhow!(
                "in-process mapper failed: {worker_err:#}\nprofiler then stopped: {profile_err:#}"
            )),
        }
    }

    fn mapping_metadata_for_thresholds(&self) -> Result<Option<MappingRunMetadata>> {
        if let Some(path) = self.config.reuse_mappings.as_deref() {
            return load_mapping_run_metadata(Path::new(path)).map(Some);
        }
        if self.config.resume {
            if let Some(output) = existing_mapping_output(&self.config.outprefix) {
                return load_mapping_run_metadata(Path::new(&output.mappings_tsv)).map(Some);
            }
        }
        Ok(None)
    }

    fn prepare_mapping(&self, profiler: MappingProfilerMetadata) -> Result<PreparedMapping> {
        // Respect explicit reuse and resume-before-aligning behavior.
        if let Some(ref reuse) = self.config.reuse_sam {
            log::info!(
                target: "MAIN",
                "Taxonomy profiling based on specified mapping file."
            );
            log::info!(
                target: "MAIN",
                "Reusing existing SAM file at {} for profiling.",
                reuse
            );
            return Ok(PreparedMapping::Sam {
                path: reuse.clone(),
            });
        }
        if let Some(ref reuse) = self.config.reuse_mappings {
            let db =
                self.config.db.as_ref().ok_or_else(|| {
                    anyhow!("--reuse-mapping requires the matching --db directory")
                })?;
            log::info!(
                target: "MAIN",
                "Profiling mapping records from {} against database metadata in {} without running the mapper.",
                reuse,
                db
            );
            return Ok(PreparedMapping::Metam {
                output: MappingOutput::reused(Path::new(reuse).to_path_buf()),
                db: db.clone(),
            });
        }
        if self.config.resume {
            if let Some(output) = existing_mapping_output(&self.config.outprefix) {
                log::info!(
                    target: "MAIN",
                    "Resume taxonomy profiling with existing mapping output."
                );
                let db = self
                    .config
                    .db
                    .as_ref()
                    .ok_or_else(|| anyhow!("--resume requires the matching --db directory"))?;
                return Ok(PreparedMapping::Metam {
                    output,
                    db: db.clone(),
                });
            }
            log::warn!(
                target: "MAIN",
                "Resume requested, but {}.mapping.tsv.gz was not found; running mapping from scratch.",
                self.config.outprefix
            );
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
        // Invoke the in-process mapper with preset + mode parameters.
        let mapper = Mapper::new(db.clone(), sequences, self.config.outprefix.clone());
        if self.config.verbose || self.config.very_verbose {
            log::info!(
                target: "MAIN",
                "Running mapping with {} threads against {}.",
                self.config.threads,
                db
            );
        } else {
            log::info!(
                target: "MAIN",
                "Running mapping with {} threads.",
                self.config.threads
            );
        }
        let progress_log = MappingProgressLog::streamed_default(!self.config.no_progress);
        if self.config.stream_mappings {
            let profile = new_shared_metam_profile();
            return Ok(PreparedMapping::MetamStream {
                mapper,
                profile,
                db: db.clone(),
            });
        }
        if self.config.keep_mapping_output {
            let output = mapper.run_saved(
                self.config.threads,
                self.config.sequencer,
                self.config.is_paired,
                self.config.mode,
                profiler,
                self.config.very_verbose,
                &self.config.extra_args,
                self.profile_progress_options(),
                progress_log,
            )?;
            return Ok(PreparedMapping::Metam {
                output,
                db: db.clone(),
            });
        }

        let plan = plan_staging_directory(self.config.temp_dir.as_deref(), &self.config.outprefix)?;
        log::info!(
            target: "MAPPER",
            "Staging mappings in '{}'{}; the file will be removed after successful profiling.",
            plan.directory.display(),
            plan.available_bytes
                .map(|bytes| format!(" ({} available)", human_bytes(bytes)))
                .unwrap_or_default()
        );
        let staged = mapper.run_temporary(
            &plan.directory,
            self.config.threads,
            self.config.sequencer,
            self.config.is_paired,
            self.config.mode,
            profiler,
            self.config.very_verbose,
            &self.config.extra_args,
            self.profile_progress_options(),
            progress_log,
        );
        let output = match staged {
            Ok(output) => output,
            Err(error) if plan.fallback_directory.is_some() && is_storage_exhausted(&error) => {
                let fallback = plan.fallback_directory.as_ref().expect("checked above");
                log::warn!(
                    target: "MAPPER",
                    "Temporary mapping storage filled while mapping; retrying in output directory '{}'. The mapping stage will run again.",
                    fallback.display()
                );
                mapper.run_temporary(
                    fallback,
                    self.config.threads,
                    self.config.sequencer,
                    self.config.is_paired,
                    self.config.mode,
                    profiler,
                    self.config.very_verbose,
                    &self.config.extra_args,
                    self.profile_progress_options(),
                    progress_log,
                )?
            }
            Err(error) => return Err(error),
        };
        Ok(PreparedMapping::Metam {
            output,
            db: db.clone(),
        })
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
                flags.push("--low-biomass".to_string());
            }
            if self.config.keep_raw {
                flags.push("--save-raw-profile".to_string());
            }
            if self.config.save_classifications {
                flags.push("--save-classification".to_string());
            }
            if self.config.by_aligned {
                flags.push("--by-mapped".to_string());
            }
            if self.config.strict {
                flags.push("--strict".to_string());
            }
            if self.config.genus_fallback {
                flags.push("--genus-fallback".to_string());
            }
            if self.config.lenient_mapping_input {
                flags.push("--lenient-mapping-input".to_string());
            }
            if self.config.keep_mapping_output {
                flags.push("--save-mapping".to_string());
            }
            if self.config.stream_mappings {
                flags.push("--stream-mappings".to_string());
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
                "    --reuse-mapping {} \\\n",
                "    --temp-dir {} \\\n",
                "    --sequencer {} \\\n",
                "    --mode {} \\\n",
                "    --exact-ties {} \\\n",
                "    --batch-size {} \\\n",
                "    --max-em-iterations {} \\\n",
                "    --min-identity {} \\\n",
                "    -m {} \\\n",
                "    --min-breadth {} \\\n",
                "    --min-chunk-breadth {} \\\n",
                "    --chunk-breadth-scale {} \\\n",
                "    --min-reads {} \\\n",
                "    --min-oebr {} \\\n",
                "    --min-coebr {} \\\n",
                "    -f {} \\\n",
                "    --pathogen-table {} \\\n",
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
            fmt_option(&self.config.reuse_mappings),
            fmt_option(&self.config.temp_dir),
            self.config.sequencer,
            mode,
            self.config.exact_ties,
            fmt_option(&self.config.batch_size),
            self.config.max_em_iterations,
            fmt_option(&self.config.identity),
            fmt_option(&self.config.mapped_len),
            fmt_option(&self.config.breadth),
            fmt_option(&self.config.chunk_breadth),
            self.config.chunk_breadth_scale,
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

    /// Sample mean predicted base accuracy for adaptive long-read community
    /// profiling. Best-effort: FASTA, empty input, or a read error leaves the
    /// fixed preset in effect rather than failing an otherwise valid run.
    fn sample_long_read_base_accuracy(&self) -> Option<f64> {
        let first = self
            .config
            .input_sequences
            .as_deref()?
            .split(',')
            .next()
            .filter(|path| !path.is_empty())?;
        match metam::io::sample_mean_base_accuracy(
            first,
            QUALITY_SAMPLE_READS,
            QUALITY_SAMPLE_MAX_BASES,
        ) {
            Ok(accuracy) => accuracy,
            Err(error) => {
                log::warn!(
                    target: "PROFILER",
                    "Could not sample base quality from {first} for the adaptive long-read identity cutoff ({error:#}); keeping the preset identity cutoff.",
                );
                None
            }
        }
    }

    fn resolve_thresholds(
        &self,
        mapping_metadata: Option<&MappingRunMetadata>,
    ) -> ResolvedProfilerThresholds {
        // Fill in default thresholds based on sequencer and profile mode.
        let sequencer = self.config.sequencer;
        let mut batch_size = self.config.batch_size.unwrap_or(5000);
        let recorded = mapping_metadata
            .filter(|metadata| {
                metadata.mapper_mode == self.config.mode.as_str()
                    && metadata.profiler.strain == self.config.strain
            })
            .map(|metadata| metadata.profiler);
        let mut identity = self
            .config
            .identity
            .or_else(|| recorded.map(|profile| profile.identity));
        let mut mapped_len = self
            .config
            .mapped_len
            .or_else(|| recorded.map(|profile| profile.mapped_len));
        let mut fraction = self
            .config
            .fraction
            .or_else(|| recorded.map(|profile| profile.fraction));
        let mut sampled_mean_base_accuracy =
            mapping_metadata.and_then(|metadata| metadata.profiler.sampled_mean_base_accuracy);

        if recorded.is_some() {
            log::info!(
                target: "PROFILER",
                "Using profiling settings saved in the mapping file; command-line settings take priority.",
            );
        }

        match sequencer {
            Sequencer::Illumina => {
                if batch_size == 0 {
                    batch_size = 5000;
                }
                match self.config.mode {
                    ProfileMode::Recall => {
                        if self.config.strain {
                            identity.get_or_insert(0.96);
                            mapped_len.get_or_insert(50);
                            fraction.get_or_insert(0.65);
                        } else {
                            identity.get_or_insert(0.96);
                            mapped_len.get_or_insert(50);
                            fraction.get_or_insert(0.6);
                        }
                    }
                    ProfileMode::Precision => {
                        if self.config.strain {
                            identity.get_or_insert(0.98);
                            mapped_len.get_or_insert(60);
                            fraction.get_or_insert(0.8);
                        } else {
                            identity.get_or_insert(0.98);
                            mapped_len.get_or_insert(50);
                            fraction.get_or_insert(0.8);
                        }
                    }
                    ProfileMode::Default => {
                        // Strain profiling uses a stricter mapped-fraction
                        // threshold than community profiling.
                        if self.config.strain {
                            identity.get_or_insert(0.97);
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
            Sequencer::Nanopore | Sequencer::PacBio => {
                if batch_size == 0 {
                    batch_size = 5000;
                }
                // Strain profiling keeps its own long-read presets. Length and
                // fraction are combined with OR for long reads.
                if self.config.strain {
                    mapped_len.get_or_insert(500);
                    fraction.get_or_insert(0.30);
                    identity.get_or_insert(match self.config.mode {
                        ProfileMode::Recall => 0.90,
                        ProfileMode::Default => 0.92,
                        ProfileMode::Precision => 0.94,
                    });
                } else {
                    let (default_identity, default_mapped, default_fraction) =
                        match self.config.mode {
                            // Recall relaxes both identity and the
                            // absolute-length gate.
                            ProfileMode::Recall => (0.90, 300, 0.20),
                            ProfileMode::Default => (0.92, 500, 0.20),
                            ProfileMode::Precision => (0.94, 500, 0.30),
                        };
                    // Low-quality long-read samples self-calibrate the identity
                    // cutoff from mean predicted base accuracy. Reused mappings
                    // carry both the original resolved cutoff and sampled
                    // accuracy, so same-mode reuse is exact and a deliberately
                    // changed mode can apply its own margin without the FASTQ.
                    if identity.is_none() {
                        let accuracy_source = if sampled_mean_base_accuracy.is_some() {
                            "saved mapping metadata"
                        } else {
                            sampled_mean_base_accuracy = self.sample_long_read_base_accuracy();
                            "input base qualities"
                        };
                        let resolved_identity = sampled_mean_base_accuracy
                            .and_then(|accuracy| {
                                let cutoff =
                                    adaptive_identity_from_accuracy(accuracy, self.config.mode)?;
                                log::info!(
                                    target: "PROFILER",
                                    "Long-read identity threshold: {cutoff:.4} (auto-calibrated from {accuracy_source}: mean predicted accuracy {accuracy:.4} < {ADAPTIVE_IDENTITY_CEILING:.2}, floor {ADAPTIVE_IDENTITY_FLOOR:.2}).",
                                );
                                Some(cutoff)
                            })
                            .unwrap_or(default_identity);
                        if let Some(accuracy) = sampled_mean_base_accuracy {
                            if accuracy >= ADAPTIVE_IDENTITY_CEILING {
                                log::info!(
                                    target: "PROFILER",
                                    "Long-read identity: mean predicted base accuracy {accuracy:.4} from {accuracy_source} >= {ADAPTIVE_IDENTITY_CEILING:.2}; keeping the preset identity cutoff.",
                                );
                            }
                        }
                        identity = Some(resolved_identity);
                    }
                    mapped_len.get_or_insert(default_mapped);
                    fraction.get_or_insert(default_fraction);
                }
            }
            Sequencer::Assembly => {
                if batch_size == 0 {
                    batch_size = 5000;
                }
                match self.config.mode {
                    ProfileMode::Recall => {
                        identity.get_or_insert(0.92);
                        mapped_len.get_or_insert(400);
                        fraction.get_or_insert(0.55);
                    }
                    ProfileMode::Precision => {
                        identity.get_or_insert(0.95);
                        mapped_len.get_or_insert(600);
                        fraction.get_or_insert(0.70);
                    }
                    ProfileMode::Default => {
                        identity.get_or_insert(0.93);
                        mapped_len.get_or_insert(500);
                        fraction.get_or_insert(0.60);
                    }
                }
            }
        }

        ResolvedProfilerThresholds {
            batch_size,
            identity: identity.unwrap_or(0.0),
            mapped_len: mapped_len.unwrap_or(0),
            fraction: fraction.unwrap_or(0.0),
            sampled_mean_base_accuracy,
        }
    }

    fn build_profiler_config(
        &self,
        mapping: &PreparedMapping,
        thresholds: ResolvedProfilerThresholds,
    ) -> Result<ProfilerConfig> {
        let ResolvedProfilerThresholds {
            batch_size,
            identity,
            mapped_len,
            fraction,
            ..
        } = thresholds;
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
        let suppress_mapping_load_progress = matches!(mapping, PreparedMapping::MetamStream { .. });
        let (sam, mapper_mappings, mapping_profile_memory, db_metadata) = match mapping {
            PreparedMapping::Sam { path, .. } => (path.clone(), None, None, None),
            PreparedMapping::Metam { output, db } => {
                let mapping_path = Path::new(&output.mappings_tsv);
                let database = Arc::new(DatabaseMetadata::load(Path::new(db))?);
                let metadata = load_mapping_run_metadata(mapping_path)?;
                validate_mapping_database_summary(metadata.database, &database)?;
                if metadata.sequencer != self.config.sequencer
                    || metadata.is_paired != self.config.is_paired
                {
                    return Err(anyhow!(
                        "mapping metadata declares sequencer={} paired={}, but profiling resolved sequencer={} paired={}",
                        metadata.sequencer,
                        metadata.is_paired,
                        self.config.sequencer,
                        self.config.is_paired
                    ));
                }
                if metadata.mapper_mode != self.config.mode.as_str() {
                    log::warn!(
                        target: "PROFILER",
                        "Reprofiling in {} mode from mappings generated in {} mode. Profiler thresholds will change, but the saved mapper candidate set cannot be expanded or regenerated.",
                        self.config.mode.as_str(),
                        metadata.mapper_mode
                    );
                }
                (
                    String::new(),
                    Some(output.mappings_tsv.clone()),
                    None,
                    Some(database),
                )
            }
            PreparedMapping::MetamStream { profile, db, .. } => (
                String::new(),
                None,
                Some(profile.clone()),
                Some(Arc::new(DatabaseMetadata::load(Path::new(db))?)),
            ),
        };
        Ok(ProfilerConfig {
            sam,
            mapper_mappings,
            mapping_profile_memory,
            db_metadata,
            sequencer: self.config.sequencer,
            batch_size,
            max_em_iterations: self.config.max_em_iterations.max(1),
            lenient_mapping_input: self.config.lenient_mapping_input,
            is_paired: self.config.is_paired,
            exact_ties: self.config.exact_ties,
            identity,
            mapped_len,
            breadth: self.config.breadth,
            chunk_breadth: self.config.chunk_breadth,
            chunk_breadth_scale: self.config.chunk_breadth_scale,
            min_reads: self.config.min_reads,
            min_oebr: self.config.min_oebr,
            min_coebr: self.config.min_coebr,
            fraction,
            lowbiomass: self.config.lowbiomass,
            keep_raw: self.config.keep_raw,
            save_classifications: self.config.save_classifications,
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
            progress: self.profile_progress_options(),
            suppress_mapping_load_progress,
        })
    }

    fn profile_progress_options(&self) -> ProgressOptions {
        ProgressOptions::new(
            self.config.progress,
            self.config.no_progress,
            self.config.progress_interval,
        )
    }
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
    metam_options: MetamIndexOptions,
) -> Result<()> {
    let params = subsample.map(|fraction| crate::db::SubsampleParams {
        fraction,
        seed,
        segment_length,
        min_length,
    });
    build_database(
        fasta,
        outdir,
        params,
        threads,
        compress_subsampled,
        &metam_options,
    )?;
    Ok(())
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
            reuse_mappings: None,
            sequencer: Sequencer::Illumina,
            is_paired: false,
            exact_ties: ExactTiePolicy::GuardedMap,
            strain: false,
            mode: ProfileMode::Default,
            batch_size: None,
            max_em_iterations: 200,
            lenient_mapping_input: false,
            identity: None,
            mapped_len: None,
            breadth: None,
            chunk_breadth: None,
            chunk_breadth_scale: 1.0,
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            fraction: None,
            lowbiomass: false,
            keep_raw: false,
            save_classifications: false,
            by_aligned: false,
            genus_identity: 0.80,
            low_map_rate_threshold: 0.30,
            genus_fallback: false,
            strict: false,
            keep_mapping_output: false,
            stream_mappings: false,
            temp_dir: None,
            pathogen_host: None,
            host: None,
            verbose: false,
            very_verbose: false,
            progress: false,
            no_progress: false,
            progress_interval: 100_000,
            extra_args: Vec::new(),
        }
    }

    #[test]
    fn library_configuration_validates_chunk_breadth_scaling_invariants() {
        let mut valid = test_app_config();
        valid.chunk_breadth_scale = 0.5;
        MetaxApp { config: valid }
            .validate_chunk_breadth_configuration()
            .unwrap();

        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut cfg = test_app_config();
            cfg.chunk_breadth_scale = invalid;
            assert!(
                MetaxApp { config: cfg }
                    .validate_chunk_breadth_configuration()
                    .is_err(),
                "factor {invalid:?} should be rejected"
            );
        }

        let mut explicit = test_app_config();
        explicit.chunk_breadth_scale = 0.5;
        explicit.chunk_breadth = Some(0.1);
        assert!(MetaxApp { config: explicit }
            .validate_chunk_breadth_configuration()
            .is_err());

        let mut min_reads = test_app_config();
        min_reads.chunk_breadth_scale = 0.5;
        min_reads.min_reads = Some(3);
        assert!(MetaxApp { config: min_reads }
            .validate_chunk_breadth_configuration()
            .is_err());

        let mut low_biomass = test_app_config();
        low_biomass.chunk_breadth_scale = 0.5;
        low_biomass.lowbiomass = true;
        assert!(MetaxApp {
            config: low_biomass
        }
        .validate_chunk_breadth_configuration()
        .is_err());
    }

    fn resolve(cfg: AppConfig) -> (usize, f64, usize, f64) {
        let resolved = MetaxApp { config: cfg }.resolve_thresholds(None);
        (
            resolved.batch_size,
            resolved.identity,
            resolved.mapped_len,
            resolved.fraction,
        )
    }

    fn recorded_long_read_metadata(
        mode: ProfileMode,
        identity: f64,
        mapped_len: usize,
        fraction: f64,
        sampled_mean_base_accuracy: Option<f64>,
    ) -> MappingRunMetadata {
        MappingRunMetadata {
            sequencer: Sequencer::Nanopore,
            is_paired: false,
            mapper_mode: mode.as_str().to_string(),
            profiler: MappingProfilerMetadata {
                strain: false,
                identity,
                mapped_len,
                fraction,
                sampled_mean_base_accuracy,
            },
            database: crate::mapper::tsv::MappingDatabaseSummary {
                is_fractional: false,
                sequences: 1,
                genomes: 1,
                fingerprint: 1,
            },
        }
    }

    #[test]
    fn adaptive_identity_keeps_preset_at_or_above_ceiling() {
        assert!(adaptive_identity_from_accuracy(0.95, ProfileMode::Default).is_none());
        assert!(adaptive_identity_from_accuracy(0.99, ProfileMode::Recall).is_none());
        assert!(adaptive_identity_from_accuracy(f64::NAN, ProfileMode::Default).is_none());
    }

    #[test]
    fn adaptive_identity_applies_mode_margin_below_ceiling() {
        let d = adaptive_identity_from_accuracy(0.90, ProfileMode::Default).unwrap();
        let r = adaptive_identity_from_accuracy(0.90, ProfileMode::Recall).unwrap();
        let p = adaptive_identity_from_accuracy(0.90, ProfileMode::Precision).unwrap();
        assert!((d - 0.87).abs() < 1e-12, "default={d}");
        assert!((r - 0.86).abs() < 1e-12, "recall={r}");
        assert!((p - 0.88).abs() < 1e-12, "precision={p}");
    }

    #[test]
    fn adaptive_identity_never_drops_below_floor() {
        // 0.85 - 0.03 = 0.82 -> floored to 0.84.
        let d = adaptive_identity_from_accuracy(0.85, ProfileMode::Default).unwrap();
        assert!((d - 0.84).abs() < 1e-12, "default={d}");
        // 0.86 - 0.04 = 0.82 -> floored to 0.84.
        let r = adaptive_identity_from_accuracy(0.86, ProfileMode::Recall).unwrap();
        assert!((r - 0.84).abs() < 1e-12, "recall={r}");
    }

    #[test]
    fn long_read_low_quality_sample_lowers_identity_below_preset() {
        // A Nanopore sample whose head reads are all Q10 ('+', accuracy 0.90)
        // self-calibrates identity to 0.90 - 0.03 = 0.87, below the 0.92 preset.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fastq");
        std::fs::write(
            &path,
            "@a\nACGTACGT\n+\n++++++++\n@b\nACGTACGT\n+\n++++++++\n",
        )
        .unwrap();
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Nanopore;
        cfg.mode = ProfileMode::Default;
        cfg.input_sequences = Some(path.to_str().unwrap().to_string());
        let (_, identity, _, _) = resolve(cfg);
        assert!((identity - 0.87).abs() < 1e-9, "identity={identity}");
    }

    #[test]
    fn long_read_high_quality_sample_keeps_preset_identity() {
        // Head reads at Q40 ('I', accuracy 0.9999 >= 0.95) keep the 0.92 preset.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fastq");
        std::fs::write(&path, "@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n").unwrap();
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Nanopore;
        cfg.mode = ProfileMode::Default;
        cfg.input_sequences = Some(path.to_str().unwrap().to_string());
        let (_, identity, _, _) = resolve(cfg);
        assert!((identity - 0.92).abs() < 1e-12, "identity={identity}");
    }

    #[test]
    fn long_read_explicit_identity_overrides_adaptive() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fastq");
        std::fs::write(&path, "@a\nACGT\n+\n++++\n").unwrap();
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Nanopore;
        cfg.mode = ProfileMode::Default;
        cfg.input_sequences = Some(path.to_str().unwrap().to_string());
        cfg.identity = Some(0.91); // explicit --min-identity wins over adaptive
        let (_, identity, _, _) = resolve(cfg);
        assert!((identity - 0.91).abs() < 1e-12, "identity={identity}");
    }

    #[test]
    fn same_mode_reuse_restores_exact_recorded_profiler_gates() {
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Nanopore;
        cfg.mode = ProfileMode::Default;
        let metadata =
            recorded_long_read_metadata(ProfileMode::Default, 0.87321, 417, 0.23, Some(0.90321));

        let resolved = MetaxApp { config: cfg }.resolve_thresholds(Some(&metadata));
        assert_eq!(resolved.identity.to_bits(), 0.87321_f64.to_bits());
        assert_eq!(resolved.mapped_len, 417);
        assert_eq!(resolved.fraction.to_bits(), 0.23_f64.to_bits());
        assert_eq!(
            resolved.sampled_mean_base_accuracy.map(f64::to_bits),
            Some(0.90321_f64.to_bits())
        );
    }

    #[test]
    fn explicit_profiler_gates_override_recorded_values() {
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Nanopore;
        cfg.mode = ProfileMode::Default;
        cfg.identity = Some(0.91);
        cfg.mapped_len = Some(700);
        cfg.fraction = Some(0.40);
        let metadata =
            recorded_long_read_metadata(ProfileMode::Default, 0.87, 400, 0.20, Some(0.90));

        let resolved = MetaxApp { config: cfg }.resolve_thresholds(Some(&metadata));
        assert_eq!(resolved.identity.to_bits(), 0.91_f64.to_bits());
        assert_eq!(resolved.mapped_len, 700);
        assert_eq!(resolved.fraction.to_bits(), 0.40_f64.to_bits());
    }

    #[test]
    fn changed_reuse_mode_recalibrates_from_recorded_mean_accuracy() {
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Nanopore;
        cfg.mode = ProfileMode::Recall;
        let metadata =
            recorded_long_read_metadata(ProfileMode::Default, 0.87, 500, 0.20, Some(0.90));

        let resolved = MetaxApp { config: cfg }.resolve_thresholds(Some(&metadata));
        assert!((resolved.identity - 0.86).abs() < 1e-12);
        assert_eq!(resolved.mapped_len, 300);
        assert!((resolved.fraction - 0.20).abs() < 1e-12);
    }

    #[test]
    fn short_reads_ignore_the_adaptive_long_read_path() {
        // Illumina keeps its fixed preset even with a low-quality FASTQ.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fastq");
        std::fs::write(&path, "@a\nACGT\n+\n++++\n").unwrap();
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Illumina;
        cfg.mode = ProfileMode::Default;
        cfg.input_sequences = Some(path.to_str().unwrap().to_string());
        let (_, identity, _, _) = resolve(cfg);
        assert!((identity - 0.97).abs() < 1e-12, "identity={identity}");
    }

    // ------------------------------------------------------------------
    // Illumina + non-strain + Default mode => identity 0.97 / fraction 0.7
    // ------------------------------------------------------------------

    #[test]
    fn illumina_default_non_strain_uses_new_community_defaults() {
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Illumina;
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
        cfg.sequencer = Sequencer::Illumina;
        cfg.is_paired = true;
        cfg.strain = false;
        cfg.mode = ProfileMode::Default;
        let (_, identity, _, fraction) = resolve(cfg);
        assert!((identity - 0.97).abs() < 1e-12);
        assert!((fraction - 0.70).abs() < 1e-12);
    }

    #[test]
    fn illumina_default_strain_uses_strain_specific_defaults() {
        // Strain mode on the Default preset uses identity 0.97 and
        // fraction 0.8 (mapped_len 50).
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Illumina;
        cfg.strain = true;
        cfg.mode = ProfileMode::Default;
        let (_, identity, mapped_len, fraction) = resolve(cfg);
        assert!((identity - 0.97).abs() < 1e-12);
        assert_eq!(mapped_len, 50);
        assert!((fraction - 0.80).abs() < 1e-12);
    }

    #[test]
    fn illumina_default_strain_paired_also_uses_strain_specific_defaults() {
        // Paired-end short-read strain mode must match single-end.
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Illumina;
        cfg.is_paired = true;
        cfg.strain = true;
        cfg.mode = ProfileMode::Default;
        let (_, identity, _, fraction) = resolve(cfg);
        assert!((identity - 0.97).abs() < 1e-12);
        assert!((fraction - 0.80).abs() < 1e-12);
    }

    #[test]
    fn strain_defaults_apply_only_to_illumina() {
        // The short-read strain thresholds are Illumina-only; a
        // Nanopore strain run uses the long-read length-or-fraction default
        // (mapped_len = 500) and the long-read strain identity (Default = 0.92).
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Nanopore;
        cfg.strain = true;
        cfg.mode = ProfileMode::Default;
        let (_, identity, mapped_len, _) = resolve(cfg);
        assert_eq!(mapped_len, 500);
        assert!(
            (identity - 0.92).abs() < 1e-12,
            "Nanopore strain Default identity should be 0.92, got {identity}"
        );
    }

    #[test]
    fn illumina_community_precision_and_recall_use_community_defaults() {
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Illumina;
        cfg.strain = false;
        cfg.mode = ProfileMode::Precision;
        let (_, identity, mapped_len, fraction) = resolve(cfg.clone());
        assert!((identity - 0.98).abs() < 1e-12);
        assert_eq!(mapped_len, 50);
        assert!((fraction - 0.80).abs() < 1e-12);

        let mut cfg_recall = cfg;
        cfg_recall.mode = ProfileMode::Recall;
        let (_, identity, mapped_len, fraction) = resolve(cfg_recall);
        assert!((identity - 0.96).abs() < 1e-12);
        assert_eq!(mapped_len, 50);
        assert!((fraction - 0.60).abs() < 1e-12);
    }

    #[test]
    fn illumina_strain_precision_and_recall_use_requested_defaults() {
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Illumina;
        cfg.strain = true;
        cfg.mode = ProfileMode::Precision;
        let (_, identity, mapped_len, fraction) = resolve(cfg.clone());
        assert!((identity - 0.98).abs() < 1e-12);
        assert_eq!(mapped_len, 60);
        assert!((fraction - 0.80).abs() < 1e-12);

        let mut cfg_recall = cfg;
        cfg_recall.mode = ProfileMode::Recall;
        let (_, identity, mapped_len, fraction) = resolve(cfg_recall);
        assert!((identity - 0.96).abs() < 1e-12);
        assert_eq!(mapped_len, 50);
        assert!((fraction - 0.65).abs() < 1e-12);
    }

    #[test]
    fn user_overrides_beat_any_default() {
        // Explicit CLI values must pass through untouched regardless of
        // the Illumina-non-strain-Default tightening.
        let mut cfg = test_app_config();
        cfg.sequencer = Sequencer::Illumina;
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
    fn long_read_defaults_use_length_or_fraction_gate() {
        // Community long-read defaults are identity 0.92 with mapped_len >= 500
        // OR mapped fraction >= 0.20.
        for seq in [Sequencer::Nanopore, Sequencer::PacBio] {
            let mut cfg = test_app_config();
            cfg.sequencer = seq;
            cfg.strain = false;
            cfg.mode = ProfileMode::Default;
            let (_, identity, mapped_len, fraction) = resolve(cfg);
            assert_eq!(mapped_len, 500, "{}", seq);
            assert!((fraction - 0.20).abs() < 1e-12, "{}", seq);
            assert!((identity - 0.92).abs() < 1e-12, "{}", seq);
        }
    }

    #[test]
    fn long_read_community_precision_and_recall_use_community_defaults() {
        for seq in [Sequencer::Nanopore, Sequencer::PacBio] {
            let mut cfg = test_app_config();
            cfg.sequencer = seq;
            cfg.strain = false;
            cfg.mode = ProfileMode::Precision;
            let (_, identity, mapped_len, fraction) = resolve(cfg.clone());
            assert!((identity - 0.94).abs() < 1e-12, "{}", seq);
            assert_eq!(mapped_len, 500, "{}", seq);
            assert!((fraction - 0.30).abs() < 1e-12, "{}", seq);

            let mut cfg_recall = cfg;
            cfg_recall.mode = ProfileMode::Recall;
            let (_, identity, mapped_len, fraction) = resolve(cfg_recall);
            assert!((identity - 0.90).abs() < 1e-12, "{}", seq);
            assert_eq!(mapped_len, 300, "{}", seq);
            assert!((fraction - 0.20).abs() < 1e-12, "{}", seq);
        }
    }

    #[test]
    fn long_read_strain_uses_length_or_fraction_defaults() {
        // Strain long-read runs also use the length-or-fraction defaults
        // (mapped_len = 500, fraction = 0.3) across modes; the strain identity
        // is stricter than community by mode (Precision = 0.94, Recall = 0.90).
        for seq in [Sequencer::Nanopore, Sequencer::PacBio] {
            let mut cfg = test_app_config();
            cfg.sequencer = seq;
            cfg.strain = true;
            cfg.mode = ProfileMode::Precision;
            let (_, identity, mapped_len, fraction) = resolve(cfg.clone());
            assert!((identity - 0.94).abs() < 1e-12, "{}", seq);
            assert_eq!(mapped_len, 500, "{}", seq);
            assert!((fraction - 0.30).abs() < 1e-12, "{}", seq);

            let mut cfg_recall = cfg;
            cfg_recall.mode = ProfileMode::Recall;
            let (_, identity, mapped_len, fraction) = resolve(cfg_recall);
            assert!((identity - 0.90).abs() < 1e-12, "{}", seq);
            assert_eq!(mapped_len, 500, "{}", seq);
            assert!((fraction - 0.30).abs() < 1e-12, "{}", seq);
        }
    }
}

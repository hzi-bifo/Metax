use std::cmp::Ordering;
use std::fmt;
use std::fs::{File, Permissions};
use std::io::{BufWriter, Write};
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, SyncSender};
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

use anyhow::{anyhow, bail, Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use tempfile::{Builder, TempPath};

use metam::format::format_mapping;
use metam::index::{read_index_with_cache_progress, IndexCacheBuildProgress, LoadedIndex};
use metam::io::{Record, SequenceReader};
use metam::mapper::{
    apply_pair_awareness, has_proper_pair, limit_mappings_with_workspace, map_read_with_workspace,
    pair_compatibility, rescue_competing_mates_with_workspace, rescue_mate_with_workspace,
    set_mate_status, CandidateBoundary, IdentityMode, LongReadExtension, MapParams, MapProfile,
    Mapping, MappingWorkspace, PairStatus, ReferenceCandidateGroups, RuntimeIndex,
    RuntimeIndexOptions,
};

use crate::db::DatabaseMetadata;
use crate::mapper::tsv::{
    format_sample_metadata_line, store_shared_metam_profile, write_mapping_run_metadata,
    MappingProfilerMetadata, MetamMapping, MetamMappingBatch, MetamMappingReceiver, MetamProfile,
    SharedMetamProfile,
};
use crate::profiler::{parse_reference_taxonomy, ProfileMode};
use crate::progress::{self, ProgressOptions};
use crate::sequencer::Sequencer;

const MAPPING_TSV_GZ_SUFFIX: &str = ".mapping.tsv.gz";
const DEFAULT_ILLUMINA_BATCH_SIZE: usize = 131_072;
const DEFAULT_LONG_READ_BATCH_SIZE: usize = 131_072;
const DEFAULT_LONG_READ_BATCH_BASE_LIMIT: usize = 1 << 30;
// Structured mappings are consumed immediately by the profiler. Keep only a
// few modest chunks in flight so the channel provides pipeline overlap without
// retaining whole mapper work batches (which may contain millions of hits).
const STRUCTURED_MAPPING_CHANNEL_CAPACITY: usize = 4;
const STRUCTURED_MAPPING_CHUNK_RECORDS: usize = 65_536;
const MAPPING_WRITER_BUFFER_BYTES: usize = 1 << 20;
const DEFAULT_STREAMED_MAPPING_PROGRESS_SECS: u64 = 10;
const DEFAULT_CACHE_BUILD_PROGRESS_SECS: u64 = 10;
const DEFAULT_MIN_INSERT: usize = 0;
const DEFAULT_MAX_INSERT: usize = 1000;
/// Taxa placed on one mate only whose missing mate is searched for, per mate,
/// when another taxon already forms a proper pair.
const DEFAULT_MATE_RESCUE_TAXA: usize = 8;
// Rendezvous channels provide double buffering without permitting an
// unbounded queue: while one batch maps, exactly one next input batch may be
// under construction and exactly one previous result batch may be written.
const PIPELINE_CHANNEL_CAPACITY: usize = 0;

pub struct MappingOutput {
    pub mappings_tsv: String,
    temporary_path: Option<TempPath>,
}

impl MappingOutput {
    fn saved(path: PathBuf) -> Self {
        Self {
            mappings_tsv: path.to_string_lossy().into_owned(),
            temporary_path: None,
        }
    }

    fn temporary(path: TempPath) -> Self {
        Self {
            mappings_tsv: path.to_string_lossy().into_owned(),
            temporary_path: Some(path),
        }
    }

    pub fn reused(path: PathBuf) -> Self {
        Self::saved(path)
    }

    pub fn cleanup_after_success(mut self) -> Result<()> {
        if let Some(path) = self.temporary_path.take() {
            let display = path.to_string_lossy().into_owned();
            path.close()
                .with_context(|| format!("failed to remove temporary mapping output {display}"))?;
            log::info!(target: "MAPPER", "Removed temporary mapping output {}.", display);
        }
        Ok(())
    }

    pub fn preserve_after_failure(mut self) -> Result<Option<PathBuf>> {
        let Some(path) = self.temporary_path.take() else {
            return Ok(None);
        };
        path.keep()
            .map(Some)
            .map_err(|error| anyhow!(error).context("failed to preserve temporary mapping output"))
    }
}

pub struct StreamingStructuredMappingOutput {
    pub receiver: MetamMappingReceiver,
    pub profile: SharedMetamProfile,
    pub worker: JoinHandle<Result<()>>,
}

#[derive(Debug)]
struct MappingConsumerClosed;

impl fmt::Display for MappingConsumerClosed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("profiler stopped consuming structured mapping stream")
    }
}

impl std::error::Error for MappingConsumerClosed {}

pub(crate) fn is_mapping_consumer_closed(error: &anyhow::Error) -> bool {
    error.downcast_ref::<MappingConsumerClosed>().is_some()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MappingProgressLog {
    enabled: bool,
    period: Duration,
}

impl MappingProgressLog {
    pub fn streamed_default(enabled: bool) -> Self {
        Self {
            enabled,
            period: Duration::from_secs(DEFAULT_STREAMED_MAPPING_PROGRESS_SECS),
        }
    }
}

struct MappingProgress {
    bar: indicatif::ProgressBar,
    /// Whether the live bar is actually drawn (a TTY, or forced with
    /// `--progress`). When drawn it carries the full running counts, so the
    /// periodic log line is suppressed to avoid garbling it; when hidden
    /// (redirected output) the log line is the only progress signal.
    bar_live: bool,
    log: MappingProgressLog,
    last_log_time: Instant,
    last_logged_count: u64,
    last_label: &'static str,
    start: Instant,
}

struct IndexCachePreparationProgress {
    bar: indicatif::ProgressBar,
    bar_live: bool,
    log_enabled: bool,
    total_bytes: u64,
    start: Instant,
    last_log_time: Instant,
}

struct PerReadMappingOutput {
    read_records: u64,
    mapped_read_records: u64,
    emitted_records: u64,
    profile: MapProfile,
    text: String,
}

struct PerReadStructuredMappingOutput {
    read_records: u64,
    mapped_read_records: u64,
    emitted_records: u64,
    profile: MapProfile,
    mappings: Vec<MetamMapping>,
}

struct ResolvedMappingRun {
    inputs: MetamInputs,
    config: MappingConfig,
    database: DatabaseMetadata,
    runtime_index: Arc<RuntimeIndex>,
    pool: rayon::ThreadPool,
}

impl MappingProgress {
    fn new(progress_options: ProgressOptions, progress_log: MappingProgressLog) -> Self {
        let start = Instant::now();
        Self {
            bar: progress::spinner(progress_options, "loading"),
            bar_live: progress_options.enabled(),
            log: progress_log,
            last_log_time: start,
            last_logged_count: 0,
            last_label: "Mapping reads",
            start,
        }
    }

    fn tick(&mut self, label: &'static str, profile: &MapProfile) {
        self.last_label = label;
        if self.bar_live {
            // The live bar carries the running counts + rate and is the sole
            // progress signal on a TTY, so no periodic log line is emitted.
            let unit = mapping_unit(label);
            progress::tick_message(
                &self.bar,
                format_mapping_progress(unit, profile, self.start),
                profile.total_reads,
            );
            return;
        }
        if !self.log.enabled {
            return;
        }
        let now = Instant::now();
        if now.duration_since(self.last_log_time) < self.log.period {
            return;
        }
        self.last_log_time = now;
        self.last_logged_count = profile.total_reads;
        log_mapping_progress(label, profile, self.start);
    }

    fn finish(&self, profile: &MapProfile) {
        self.bar.finish_and_clear();
        if self.log.enabled
            && self.last_logged_count > 0
            && self.last_logged_count != profile.total_reads
        {
            log_mapping_progress(self.last_label, profile, self.start);
        }
    }
}

impl IndexCachePreparationProgress {
    fn new(options: ProgressOptions, total_bytes: u64) -> Self {
        let start = Instant::now();
        Self {
            bar: progress::bytes_bar(options, total_bytes, "Preparing index cache"),
            bar_live: options.enabled(),
            log_enabled: !options.disable,
            total_bytes,
            start,
            last_log_time: start,
        }
    }

    fn tick(&mut self, update: IndexCacheBuildProgress) {
        let bytes = update.compressed_bytes_read.min(self.total_bytes);
        if self.bar_live {
            progress::refresh_timestamp(&self.bar);
            self.bar.set_position(bytes);
            self.bar.tick();
            return;
        }
        if !self.log_enabled {
            return;
        }
        let now = Instant::now();
        if now.duration_since(self.last_log_time)
            < Duration::from_secs(DEFAULT_CACHE_BUILD_PROGRESS_SECS)
        {
            return;
        }
        self.last_log_time = now;
        self.log_status(bytes, now);
    }

    fn finish(&self) {
        self.bar.finish_and_clear();
    }

    fn log_status(&self, bytes: u64, now: Instant) {
        let elapsed = now.duration_since(self.start).as_secs_f64();
        let bytes_per_second = if elapsed > 0.0 {
            bytes as f64 / elapsed
        } else {
            0.0
        };
        let percent = if self.total_bytes == 0 {
            100.0
        } else {
            bytes as f64 * 100.0 / self.total_bytes as f64
        };
        let eta = if bytes_per_second > 0.0 {
            Duration::from_secs_f64(
                self.total_bytes.saturating_sub(bytes) as f64 / bytes_per_second,
            )
        } else {
            Duration::ZERO
        };
        log::info!(
            target: "MAPPER",
            "Preparing index cache: {}/{} ({:.1}%); {}/s; ETA {}.",
            indicatif::HumanBytes(bytes),
            indicatif::HumanBytes(self.total_bytes),
            percent,
            indicatif::HumanBytes(bytes_per_second.round() as u64),
            indicatif::HumanDuration(eta),
        );
    }
}

fn log_mapping_progress(label: &str, profile: &MapProfile, start: Instant) {
    let elapsed = start.elapsed().as_secs_f64();
    let rate = if elapsed > 0.0 {
        profile.total_reads as f64 / elapsed
    } else {
        0.0
    };
    log::info!(
        target: "MAPPER",
        "{}: {} processed ({:.1}/s); mapped_read_records={}; emitted_records={}",
        label,
        profile.total_reads,
        rate,
        profile.mapped_reads,
        profile.emitted_records
    );
}

/// Short unit noun for the live line, taken from the tick label:
/// "Mapping reads" -> "reads", "Mapping read records" -> "records".
fn mapping_unit(label: &str) -> &str {
    label.rsplit(' ').next().unwrap_or("reads")
}

/// One-line live mapping status, e.g.
/// `1.23M reads · 988k mapped · 3.46M mappings · 14.2k reads/s`. The spinner
/// itself signals the phase, so no label is prefixed. `reads` = pulled from
/// input, `mapped` = reads that produced >= 1 mapping, `mappings` = mapping
/// records emitted to the profiler; the rate is the cumulative throughput
/// (reads / elapsed).
fn format_mapping_progress(unit: &str, profile: &MapProfile, start: Instant) -> String {
    let elapsed = start.elapsed().as_secs_f64();
    let rate = if elapsed > 0.0 {
        profile.total_reads as f64 / elapsed
    } else {
        0.0
    };
    format!(
        "{reads} {unit} · {mapped} mapped · {mappings} mappings · {rate} {unit}/s",
        reads = human_compact(profile.total_reads),
        mapped = human_compact(profile.mapped_reads),
        mappings = human_compact(profile.emitted_records),
        rate = human_compact(rate.round() as u64),
        unit = unit,
    )
}

/// Compact human-readable count with trailing zeros trimmed:
/// `945`, `12.3k`, `988k`, `1.23M`, `2M`, `4.56B`.
fn human_compact(n: u64) -> String {
    if n < 1_000 {
        return n.to_string();
    }
    let (value, suffix, decimals) = if n < 1_000_000 {
        (n as f64 / 1_000.0, "k", 1usize)
    } else if n < 1_000_000_000 {
        (n as f64 / 1_000_000.0, "M", 2usize)
    } else {
        (n as f64 / 1_000_000_000.0, "B", 2usize)
    };
    let mut formatted = format!("{value:.decimals$}");
    if formatted.contains('.') {
        // 988.0 -> 988, 1.20 -> 1.2, 1.23 -> 1.23.
        let trimmed_len = formatted.trim_end_matches('0').trim_end_matches('.').len();
        formatted.truncate(trimmed_len);
    }
    formatted.push_str(suffix);
    formatted
}

fn log_mapping_summary(profile: &MapProfile) {
    log::info!(
        target: "MAPPER",
        "Mapping summary: total_reads={} mapped_read_records={} emitted_records={} elapsed={:.2}s reads_per_second={:.1} query_seeds/read={:.1} seed_hits/read={:.1} candidate_groups/read={:.1} candidates_before_filtering/read={:.1} candidates_after_filtering/read={:.1} identity_comparisons/read={:.1}",
        profile.total_reads,
        profile.mapped_reads,
        profile.emitted_records,
        profile.elapsed_seconds,
        profile.reads_per_second,
        per_read(profile.total_query_seeds, profile.total_reads),
        per_read(profile.total_seed_hits, profile.total_reads),
        per_read(profile.candidate_groups, profile.total_reads),
        per_read(profile.candidates_before_filtering, profile.total_reads),
        per_read(profile.candidates_after_filtering, profile.total_reads),
        per_read(profile.identity_comparisons, profile.total_reads),
    );
}

fn finish_mapping_profile(t0: Instant, profile: &mut MapProfile, progress: &mut MappingProgress) {
    progress.finish(profile);
    let elapsed = t0.elapsed().as_secs_f64();
    profile.elapsed_seconds = elapsed;
    profile.reads_per_second = if elapsed > 0.0 {
        profile.total_reads as f64 / elapsed
    } else {
        0.0
    };
    log_mapping_summary(profile);
}

fn per_read(count: u64, total_reads: u64) -> f64 {
    if total_reads == 0 {
        0.0
    } else {
        count as f64 / total_reads as f64
    }
}

pub struct Mapper {
    db: String,
    sequences: Vec<String>,
    outprefix: String,
}

impl Mapper {
    pub fn new(db: String, sequences: Vec<String>, outprefix: String) -> Self {
        Self {
            db,
            sequences,
            outprefix,
        }
    }

    pub fn run_saved(
        &self,
        threads: usize,
        sequencer: Sequencer,
        is_paired: bool,
        mode: ProfileMode,
        profiler: MappingProfilerMetadata,
        log_command: bool,
        extra_args: &[String],
        progress_options: ProgressOptions,
        progress_log: MappingProgressLog,
    ) -> Result<MappingOutput> {
        let output_path = PathBuf::from(mapping_tsv_gz_path(&self.outprefix));
        let output_directory = output_path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let mut partial = Builder::new()
            .prefix(".metax-mapping-")
            .suffix(".partial.tsv.gz")
            .permissions(Permissions::from_mode(0o666))
            .tempfile_in(output_directory)
            .with_context(|| {
                format!(
                    "failed to create partial mapping output in {}",
                    output_directory.display()
                )
            })?;
        let partial_path = partial.path().to_path_buf();
        self.write_mapping_file(
            partial.as_file_mut(),
            &partial_path,
            threads,
            sequencer,
            is_paired,
            mode,
            profiler,
            log_command,
            extra_args,
            progress_options,
            progress_log,
        )?;
        partial.persist(&output_path).map_err(|error| {
            anyhow!(error.error).context(format!(
                "failed to atomically publish mapping output {}",
                output_path.display()
            ))
        })?;
        log::info!(
            target: "MAPPER",
            "Saved mapping output to {}.",
            output_path.display()
        );
        Ok(MappingOutput::saved(output_path))
    }

    pub fn run_temporary(
        &self,
        directory: &Path,
        threads: usize,
        sequencer: Sequencer,
        is_paired: bool,
        mode: ProfileMode,
        profiler: MappingProfilerMetadata,
        log_command: bool,
        extra_args: &[String],
        progress_options: ProgressOptions,
        progress_log: MappingProgressLog,
    ) -> Result<MappingOutput> {
        let mut temporary = Builder::new()
            .prefix("metax-mapping-")
            .suffix(".tsv.gz")
            .tempfile_in(directory)
            .with_context(|| {
                format!(
                    "failed to create temporary mapping output in {}",
                    directory.display()
                )
            })?;
        let temporary_path = temporary.path().to_path_buf();
        self.write_mapping_file(
            temporary.as_file_mut(),
            &temporary_path,
            threads,
            sequencer,
            is_paired,
            mode,
            profiler,
            log_command,
            extra_args,
            progress_options,
            progress_log,
        )?;
        let path = temporary.into_temp_path();
        log::info!(
            target: "MAPPER",
            "Staged mapping output at {}.",
            path.display()
        );
        Ok(MappingOutput::temporary(path))
    }

    #[allow(clippy::too_many_arguments)]
    fn write_mapping_file(
        &self,
        file: &mut File,
        output_path: &Path,
        threads: usize,
        sequencer: Sequencer,
        is_paired: bool,
        mode: ProfileMode,
        profiler: MappingProfilerMetadata,
        log_command: bool,
        extra_args: &[String],
        progress_options: ProgressOptions,
        progress_log: MappingProgressLog,
    ) -> Result<()> {
        let output_display = output_path.display().to_string();
        let encoder = GzEncoder::new(file, Compression::fast());
        let mut writer = BufWriter::with_capacity(MAPPING_WRITER_BUFFER_BYTES, encoder);
        let profile = self.run_mapping_to_writer(
            &mut writer,
            threads,
            sequencer,
            is_paired,
            mode,
            profiler,
            log_command,
            extra_args,
            progress_options,
            progress_log,
        )?;
        let profile = MetamProfile::from_map_profile(&profile)?;
        writeln!(writer, "{}", format_sample_metadata_line(&profile)).with_context(|| {
            format!("failed to write mapping profile metadata to {output_display}")
        })?;
        writer
            .flush()
            .with_context(|| format!("failed to flush gzipped mapping output {output_display}"))?;
        let encoder = writer.into_inner().map_err(|error| {
            anyhow!(
                "failed to finalize buffered mapping output {}: {}",
                output_display,
                error.error()
            )
        })?;
        encoder
            .finish()
            .with_context(|| format!("failed to finish gzipped mapping output {output_display}"))?;
        log::info!(target: "MAPPER", "In-process mapping finished.");
        Ok(())
    }

    pub fn spawn_structured_stream(
        &self,
        threads: usize,
        sequencer: Sequencer,
        is_paired: bool,
        mode: ProfileMode,
        log_command: bool,
        extra_args: &[String],
        progress_options: ProgressOptions,
        progress_log: MappingProgressLog,
        profile: SharedMetamProfile,
    ) -> Result<StreamingStructuredMappingOutput> {
        let (sender, receiver) =
            mpsc::sync_channel::<MetamMappingBatch>(STRUCTURED_MAPPING_CHANNEL_CAPACITY);
        let mapper = Self {
            db: self.db.clone(),
            sequences: self.sequences.clone(),
            outprefix: self.outprefix.clone(),
        };
        let extra_args = extra_args.to_vec();
        let profile_for_worker = profile.clone();
        let worker = thread::spawn(move || {
            let profile = mapper.run_mapping_to_structured_sender(
                &sender,
                threads,
                sequencer,
                is_paired,
                mode,
                log_command,
                &extra_args,
                progress_options,
                progress_log,
            )?;
            store_shared_metam_profile(
                &profile_for_worker,
                MetamProfile::from_map_profile(&profile)?,
            )?;
            drop(sender);
            Ok(())
        });
        Ok(StreamingStructuredMappingOutput {
            receiver,
            profile,
            worker,
        })
    }

    fn run_mapping_to_writer(
        &self,
        writer: &mut (dyn Write + Send),
        threads: usize,
        sequencer: Sequencer,
        is_paired: bool,
        mode: ProfileMode,
        profiler: MappingProfilerMetadata,
        log_command: bool,
        extra_args: &[String],
        progress_options: ProgressOptions,
        progress_log: MappingProgressLog,
    ) -> Result<MapProfile> {
        let ResolvedMappingRun {
            inputs,
            config,
            database,
            runtime_index,
            pool,
        } = self.resolve_mapping_run(
            sequencer,
            is_paired,
            mode,
            log_command,
            extra_args,
            threads,
            progress_options,
        )?;
        write_mapping_run_metadata(
            writer,
            sequencer,
            is_paired,
            mode.as_str(),
            profiler,
            &database,
        )?;
        drop(database);
        let t0 = Instant::now();
        let mut profile = MapProfile::default();
        let mut progress = MappingProgress::new(progress_options, progress_log);
        if let Some(reads2) = inputs.reads2 {
            map_paired(
                &inputs.reads1,
                &reads2,
                &runtime_index,
                &pool,
                &config,
                writer,
                &mut profile,
                &mut progress,
            )?;
        } else {
            map_single(
                &inputs.reads1,
                &runtime_index,
                &pool,
                &config,
                writer,
                &mut profile,
                &mut progress,
            )?;
        }
        finish_mapping_profile(t0, &mut profile, &mut progress);
        Ok(profile)
    }

    fn run_mapping_to_structured_sender(
        &self,
        sender: &SyncSender<MetamMappingBatch>,
        threads: usize,
        sequencer: Sequencer,
        is_paired: bool,
        mode: ProfileMode,
        log_command: bool,
        extra_args: &[String],
        progress_options: ProgressOptions,
        progress_log: MappingProgressLog,
    ) -> Result<MapProfile> {
        let ResolvedMappingRun {
            inputs,
            config,
            database: _,
            runtime_index,
            pool,
        } = self.resolve_mapping_run(
            sequencer,
            is_paired,
            mode,
            log_command,
            extra_args,
            threads,
            progress_options,
        )?;
        let t0 = Instant::now();
        let mut profile = MapProfile::default();
        let mut progress = MappingProgress::new(progress_options, progress_log);
        if let Some(reads2) = inputs.reads2 {
            map_paired_structured(
                &inputs.reads1,
                &reads2,
                &runtime_index,
                &pool,
                &config,
                sender,
                &mut profile,
                &mut progress,
            )?;
        } else {
            map_single_structured(
                &inputs.reads1,
                &runtime_index,
                &pool,
                &config,
                sender,
                &mut profile,
                &mut progress,
            )?;
        }
        finish_mapping_profile(t0, &mut profile, &mut progress);
        Ok(profile)
    }

    fn resolve_mapping_run(
        &self,
        sequencer: Sequencer,
        is_paired: bool,
        mode: ProfileMode,
        log_command: bool,
        extra_args: &[String],
        threads: usize,
        progress_options: ProgressOptions,
    ) -> Result<ResolvedMappingRun> {
        let db = DatabaseMetadata::load(Path::new(&self.db))?;
        let inputs = self.input_sequences(sequencer, is_paired)?;
        let config = MappingConfig::resolve(sequencer, is_paired, mode, extra_args)?;
        validate_mapping_config(&config, inputs.reads2.as_deref())?;
        let index_path = db.metam_index.to_string_lossy().to_string();
        let cache_path = index_cache_path(&index_path, config.index_cache.as_deref());
        if log_command {
            let reads2 = inputs.reads2.as_deref().unwrap_or("None");
            let index_cache = cache_path.as_deref().unwrap_or("None");
            let extra_args_display = if extra_args.is_empty() {
                "None".to_string()
            } else {
                extra_args.join(" ")
            };
            log::info!(
                target: "MAPPER",
                "Running in-process mapping with resolved parameters:\n\
                 \tindex                     {}\n\
                 \tindex_cache               {}\n\
                 \treads1                    {}\n\
                 \treads2                    {}\n\
                 \tthreads                   {}\n\
                 \tpreset                    {:?}\n\
                 \tbatch_size                {}\n\
                 \tbatch_base_limit_bases    {}\n\
                 \ttop_n                     {}\n\
                 \tscore_ratio               {:.3}\n\
                 \tmin_seeds                 {}\n\
                 \tmin_identity              {:.3}\n\
                 \tmin_mapped_frac           {:.3}\n\
                 \tmin_mapped_len            {}\n\
                 \tno_identity_filter        {}\n\
                 \tignore_n                  {}\n\
                 \tidentity_mode             {:?}\n\
                 \tdiag_bin_size             {}\n\
                 \textension_max_mismatches  {}\n\
                 \tfull_read_interval        {}\n\
                 \tfast_long_read_identity   {}\n\
                 \tlong_read_chaining        {}\n\
                 \tlong_read_extension       {:?}\n\
                 \tquery_seed_stride         {}\n\
                 \tlong_read_rescue          {}\n\
                 \tmax_rescue_candidates     {}\n\
                 \tallow_single_strobe       {}\n\
                 \tpair_aware                {}\n\
                 \tpaired_emission           {:?}\n\
                 \tmin_insert                {}\n\
                 \tmax_insert                {}\n\
                 \tmate_rescue_taxa          {}\n\
                 \textra_mapper_args         {}",
                index_path,
                index_cache,
                inputs.reads1,
                reads2,
                threads.max(1),
                config.preset,
                config.batch_size,
                config
                    .batch_base_limit
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "none".to_string()),
                config.params.top_n,
                config.params.score_ratio,
                config.params.min_seeds,
                config.params.min_identity,
                config.params.min_mapped_frac,
                config
                    .params
                    .min_mapped_len
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "none".to_string()),
                config.params.no_identity_filter,
                config.params.ignore_n,
                config.params.identity_mode,
                config.params.diag_bin_size,
                config.params.extension_max_mismatches,
                config.params.full_read_interval,
                config.params.fast_long_read_identity,
                config.params.long_read_chaining,
                config.params.long_read_extension,
                config.params.query_seed_stride,
                config.params.long_read_rescue,
                config
                    .params
                    .long_read_rescue_limit
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "auto".to_string()),
                config.allow_single_strobe_matches,
                config.pair_aware,
                config.paired_emission,
                config.min_insert,
                config.max_insert,
                config.mate_rescue_taxa,
                extra_args_display,
            );
        }
        let index_load_start = Instant::now();
        let mut cache_progress = None;
        let loaded_index_result = read_index_with_cache_progress(
            &index_path,
            cache_path.as_deref(),
            |update| {
                let display = cache_progress.get_or_insert_with(|| {
                    log::info!(
                        target: "MAPPER",
                        "Preparing the one-time memory-mapped index cache by streaming the compressed index; subsequent runs reuse it."
                    );
                    IndexCachePreparationProgress::new(
                        progress_options,
                        update.compressed_bytes_total,
                    )
                });
                display.tick(update);
            },
        );
        let cache_build_attempted = cache_progress.is_some();
        if let Some(display) = cache_progress.as_ref() {
            display.finish();
        }
        let loaded_index = loaded_index_result?;
        if cache_build_attempted && loaded_index.is_memory_mapped() {
            log::info!(
                target: "MAPPER",
                "Memory-mapped index cache prepared in {:.2}s.",
                index_load_start.elapsed().as_secs_f64()
            );
        } else if cache_path.is_some() && !loaded_index.is_memory_mapped() {
            log::warn!(
                target: "MAPPER",
                "Memory-mapped index cache unavailable after {:.2}s; using the in-memory index. Peak memory can be substantially higher.",
                index_load_start.elapsed().as_secs_f64()
            );
        } else if log_command {
            log::info!(
                target: "MAPPER",
                "Mapping index ready in {:.2}s.",
                index_load_start.elapsed().as_secs_f64()
            );
        }
        let reference_candidate_groups = build_reference_candidate_groups(&loaded_index, &db)?;
        let mut runtime_index = RuntimeIndex::from_loaded_with_options(
            loaded_index,
            RuntimeIndexOptions {
                allow_single_strobe_matches: config.allow_single_strobe_matches,
            },
        );
        runtime_index.set_reference_candidate_groups(reference_candidate_groups)?;
        let runtime_index = Arc::new(runtime_index);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build()
            .context("failed to create mapping thread pool")?;
        Ok(ResolvedMappingRun {
            inputs,
            config,
            database: db,
            runtime_index,
            pool,
        })
    }

    fn input_sequences(&self, sequencer: Sequencer, is_paired: bool) -> Result<MetamInputs> {
        if self.sequences.is_empty() {
            bail!("Please specify input read files.");
        }
        if sequencer.is_long_read() {
            if is_paired || self.sequences.len() != 1 {
                bail!("long-read mapping preset for Nanopore/PacBio is single-end only");
            }
            return Ok(MetamInputs {
                reads1: self.sequences[0].clone(),
                reads2: None,
            });
        }
        if is_paired {
            if self.sequences.len() != 2 {
                bail!("Illumina paired-end mode expects exactly two input read files");
            }
            Ok(MetamInputs {
                reads1: self.sequences[0].clone(),
                reads2: Some(self.sequences[1].clone()),
            })
        } else {
            if self.sequences.len() != 1 {
                bail!(
                    "single-end mapping expects one input file; concatenate inputs explicitly if needed"
                );
            }
            Ok(MetamInputs {
                reads1: self.sequences[0].clone(),
                reads2: None,
            })
        }
    }
}

/// Build the mapper's compact reference-to-taxonomic-endpoint table once.
///
/// Database rows and index references are created in the same order, but an
/// explicit name check prevents a moved or hand-edited metadata file from
/// silently assigning the wrong taxon to a numeric reference ID.
fn build_reference_candidate_groups(
    index: &LoadedIndex,
    database: &DatabaseMetadata,
) -> Result<Vec<ReferenceCandidateGroups>> {
    if index.reference_count() != database.references.len() {
        bail!(
            "mapping index contains {} sequences but database metadata contains {}",
            index.reference_count(),
            database.references.len()
        );
    }

    let mut groups = Vec::with_capacity(database.references.len());
    for (reference_index, reference) in database.references.iter().enumerate() {
        let index_name = index.reference_id(reference_index as u32);
        if index_name != reference.name {
            bail!(
                "mapping index/database reference mismatch at sequence {}: index has '{}', metadata has '{}'",
                reference_index + 1,
                index_name,
                reference.name
            );
        }
        let taxonomy = parse_reference_taxonomy(&reference.taxonomy_field).with_context(|| {
            format!(
                "invalid taxonomy metadata for reference '{}' while preparing mapper candidates",
                reference.name
            )
        })?;
        groups.push(ReferenceCandidateGroups::new(
            taxonomy.profile_taxid,
            taxonomy.primary_taxid.unwrap_or(taxonomy.profile_taxid),
        ));
    }
    Ok(groups)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MappingPreset {
    Illumina,
    Ont,
}

#[derive(Clone, Debug)]
struct MappingConfig {
    preset: MappingPreset,
    params: MapParams,
    batch_size: usize,
    batch_base_limit: Option<usize>,
    allow_single_strobe_matches: bool,
    pair_aware: bool,
    paired_emission: PairedEmissionPolicy,
    min_insert: usize,
    max_insert: usize,
    mate_rescue_taxa: usize,
    index_cache: Option<String>,
}

#[derive(Default)]
struct MapArgOverrides {
    preset: Option<MappingPreset>,
    index_cache: Option<Option<String>>,
    allow_single_strobe_matches: Option<bool>,
    batch_size: Option<usize>,
    top_n: Option<usize>,
    score_ratio: Option<f64>,
    min_seeds: Option<usize>,
    min_identity: Option<f64>,
    min_mapped_frac: Option<f64>,
    no_identity_filter: bool,
    ignore_n: bool,
    identity_mode: Option<IdentityMode>,
    query_seed_stride: Option<usize>,
    no_long_read_rescue: bool,
    long_read_rescue_limit: Option<usize>,
    extend: bool,
    fast_extend: bool,
    pair_aware: bool,
    paired_emission: Option<PairedEmissionPolicy>,
    min_insert: Option<usize>,
    max_insert: Option<usize>,
    mate_rescue_taxa: Option<usize>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PairedEmissionPolicy {
    /// Emit every retained mate mapping for diagnostic output.
    All,
    /// Compatibility mode: prefer concordant mappings, but emit discordant
    /// mappings when no concordant same-reference mapping exists.
    ConcordantFirst,
    /// Normal paired-end profiling mode: emit only concordant mappings.
    ConcordantOnly,
}

impl MappingConfig {
    fn resolve(
        sequencer: Sequencer,
        is_paired: bool,
        mode: ProfileMode,
        extra_args: &[String],
    ) -> Result<Self> {
        let overrides = MapArgOverrides::parse(extra_args)?;
        let inferred_preset = if sequencer.is_long_read() {
            MappingPreset::Ont
        } else {
            MappingPreset::Illumina
        };
        let preset = overrides.preset.unwrap_or(inferred_preset);
        let mut top_n = overrides.top_n;
        let mut score_ratio = overrides.score_ratio;
        let allow_single_strobe_override = overrides.allow_single_strobe_matches;
        let mut allow_single_strobe_matches =
            allow_single_strobe_override.unwrap_or(preset == MappingPreset::Illumina);
        // Paired-end profiling benefits from retaining a broad candidate set:
        // competing species must remain visible so the profiler can preserve
        // conserved-region ambiguity instead of force-assigning a fragment.
        // Keep these independent of profile mode; precision and recall are
        // applied by the profiler's evidence thresholds.
        if preset == MappingPreset::Illumina && is_paired {
            top_n.get_or_insert(50);
            score_ratio.get_or_insert(0.60);
        }
        if mode == ProfileMode::Recall && preset == MappingPreset::Illumina {
            if score_ratio.is_none() {
                score_ratio = Some(0.80);
            }
            if allow_single_strobe_override.is_none() {
                allow_single_strobe_matches = true;
            }
        }
        // Precision is enforced by evidence thresholds, not by hiding close
        // competitors from the profiler. Unless explicitly overridden, it
        // therefore retains the same mapper score ratio as Default (0.90 for
        // unpaired Illumina/assembly input, 0.70 for long reads, and the
        // mode-independent 0.60 paired-end policy above). Illumina recall
        // alone widens candidate acquisition.
        let pair_aware = overrides.pair_aware || (preset == MappingPreset::Illumina && is_paired);
        let mut min_seeds = overrides.min_seeds;
        if min_seeds.is_none() && preset == MappingPreset::Illumina && is_paired {
            min_seeds = Some(1);
        }
        let params = map_params_from_preset(
            preset,
            top_n,
            score_ratio,
            min_seeds,
            overrides.min_identity,
            overrides.min_mapped_frac,
            overrides.no_identity_filter,
            overrides.ignore_n,
            overrides.identity_mode.unwrap_or(IdentityMode::Gapped),
            overrides.extend,
            overrides.fast_extend,
            overrides.query_seed_stride,
            !overrides.no_long_read_rescue,
            overrides.long_read_rescue_limit,
        );
        Ok(Self {
            preset,
            params,
            batch_size: overrides
                .batch_size
                .unwrap_or_else(|| default_batch_size_for_preset(preset)),
            batch_base_limit: match preset {
                MappingPreset::Illumina => None,
                MappingPreset::Ont => Some(DEFAULT_LONG_READ_BATCH_BASE_LIMIT),
            },
            allow_single_strobe_matches,
            pair_aware,
            paired_emission: overrides.paired_emission.unwrap_or(if is_paired {
                PairedEmissionPolicy::ConcordantOnly
            } else {
                PairedEmissionPolicy::All
            }),
            min_insert: overrides.min_insert.unwrap_or(DEFAULT_MIN_INSERT),
            max_insert: overrides.max_insert.unwrap_or(DEFAULT_MAX_INSERT),
            mate_rescue_taxa: overrides
                .mate_rescue_taxa
                .unwrap_or(DEFAULT_MATE_RESCUE_TAXA),
            index_cache: overrides
                .index_cache
                .unwrap_or_else(|| Some("auto".to_string())),
        })
    }
}

/// Render the mapper-only options accepted after `metax profile ... --`.
///
/// The default summary is generated from the same resolver used for real
/// mapping runs, so changing a preset cannot silently leave the help stale.
pub fn mapper_options_help() -> String {
    let defaults = |sequencer, paired, mode| {
        MappingConfig::resolve(sequencer, paired, mode, &[])
            .expect("built-in mapper defaults must always resolve")
    };
    let illumina_se = defaults(Sequencer::Illumina, false, ProfileMode::Default);
    let illumina_pe = defaults(Sequencer::Illumina, true, ProfileMode::Default);
    let long_reads = defaults(Sequencer::Nanopore, false, ProfileMode::Default);
    let illumina_recall = defaults(Sequencer::Illumina, false, ProfileMode::Recall);
    let long_length_gate = long_reads
        .params
        .min_mapped_len
        .map(|value| value.to_string())
        .unwrap_or_else(|| "none".to_string());

    format!(
        r#"Advanced mapper candidate-generation options

Usage:
  metax profile [PROFILE_OPTIONS] -- [MAPPER_OPTIONS]

The separator is significant: options before '--' configure taxonomic
profiling; options after '--' configure which candidate mappings reach the
profiler. Resolved values for a real run are shown by --very-verbose.

Resolved defaults (before mapper overrides):
  Illumina single-end / assembly:
    top-n={se_top_n}, score-ratio={se_score_ratio:.2}, min-seeds={se_min_seeds},
    min-identity={se_identity:.2}, min-mapped-frac={se_fraction:.2},
    single-strobe-matches={se_single_strobe}, batch-size={se_batch_size}
  Illumina paired-end:
    top-n={pe_top_n}, score-ratio={pe_score_ratio:.2}, min-seeds={pe_min_seeds},
    min-identity={pe_identity:.2}, min-mapped-frac={pe_fraction:.2},
    single-strobe-matches={pe_single_strobe}, pair-aware={pe_pair_aware},
    mate-rescue-taxa={mate_rescue_taxa}, paired-emission=concordant-only,
    batch-size={pe_batch_size}
  Nanopore / PacBio:
    top-n={long_top_n}, score-ratio={long_score_ratio:.2}, min-seeds={long_min_seeds},
    min-identity={long_identity:.2}, min-mapped-frac={long_fraction:.2},
    mapped-length OR gate={long_length_gate} bases, single-strobe-matches={long_single_strobe},
    batch-size={long_batch_size} (also capped at 1 GiB of read bases per batch)

Recall mode changes unpaired Illumina score-ratio to {recall_score_ratio:.2};
paired-end and long-read mapper defaults are mode-independent. Precision and
recall otherwise differ through profiler evidence thresholds.

Candidate discovery and retention:
  --top-n <N>
      Maximum retained mappings per read or mate. Larger values preserve more
      ambiguity but increase identity work, output size, and profiler memory.
  --score-ratio <FRACTION>
      Retain candidates scoring at least this fraction of the best seed score.
      Lower values improve recall but can substantially increase runtime.
  --min-seeds <N>
      Minimum supporting query seeds/anchors for a candidate. Lower values are
      more sensitive, especially for short PE reads, but can be much costlier.
  --allow-single-strobe-matches
      Permit lookup by strobe1 when the complete randstrobe pair differs.
      Improves recall but increases seed hits, runtime, and potentially memory.
  --match-both-strobes
      Require the complete randstrobe pair (alias: --require-both-strobes).
      This conflicts with --allow-single-strobe-matches.
  --query-seed-stride <N>
      Use every Nth query seed (alias: --seed-stride; default 1). Values above
      1 reduce work and sensitivity, mainly for very long reads.

Mapper-level acceptance:
  --min-identity <FRACTION>
      Minimum identity for mapper emission. This is distinct from the profiler
      option with the same name before '--'.
  --min-mapped-frac <FRACTION>
      Minimum mapped fraction for mapper emission. Long reads also pass when
      they satisfy the mapped-length OR gate shown above.
  --no-identity-filter
      Disable mapper identity and mapped-coverage rejection. Profiler evidence
      gates still apply, but runtime and mapping output can grow substantially.
  --ignore-n
      Exclude ambiguous N bases from identity comparisons.

Paired-end mapping:
  --pair-aware
      Enable pair-aware ranking and mate rescue (already enabled for Illumina PE).
  --paired-emission <POLICY>
      concordant-only (default), concordant-first, or all. The latter modes are
      mainly useful with --save-mapping for discordance diagnostics.
  --min-insert <BASES>
      Minimum proper-pair template span (default {min_insert}).
  --max-insert <BASES>
      Maximum proper-pair template span (default {max_insert}).
  --mate-rescue-taxa <N>
      When one taxon already forms a proper pair, also search for the missing
      mate of up to N other taxa placed on one mate only, and keep it if that
      pair matches at least as many bases. With 0, mate rescue is limited to
      fragments without any proper pair (default {mate_rescue_taxa}).

Long-read mapping:
  --identity-mode <MODE>
      gapped (default) or ungapped identity.
  --extend
      Extend every long-read candidate before filtering; most accurate and
      potentially expensive.
  --fast-extend
      Extend only retained long-read mappings; cheaper than --extend.
  --no-long-read-rescue
      Disable the additional low-seed long-read candidate rescue pass.
  --max-rescue-candidates <N>
      Limit candidates evaluated by long-read rescue (default: automatic).

Performance and cache:
  --batch-size <RECORDS>
      Mapper input batch size. This is distinct from the profiler --batch-size
      before '--'; reducing it can lower transient read-batch memory.
  --index-cache [PATH]
      Use the automatic memory-mapped cache, or place it at PATH. The default is
      <INDEX>.metamcache beside the mapping index.
  -p, --preset <PRESET>
      illumina or ont. Prefer the profile-level --sequencer option; this is an
      expert compatibility override and can make mapper/profiler presets differ.

Metax manages the index, reads, output, threads, progress, and mapping-file
schema. Set those with profile options before '--', not as mapper overrides.

Example:
  metax profile -d DB -i R1.fq.gz -I R2.fq.gz -o sample -- \
    --top-n 50 --score-ratio 0.60 --min-seeds 1
"#,
        se_top_n = illumina_se.params.top_n,
        se_score_ratio = illumina_se.params.score_ratio,
        se_min_seeds = illumina_se.params.min_seeds,
        se_identity = illumina_se.params.min_identity,
        se_fraction = illumina_se.params.min_mapped_frac,
        se_single_strobe = illumina_se.allow_single_strobe_matches,
        se_batch_size = illumina_se.batch_size,
        pe_top_n = illumina_pe.params.top_n,
        pe_score_ratio = illumina_pe.params.score_ratio,
        pe_min_seeds = illumina_pe.params.min_seeds,
        pe_identity = illumina_pe.params.min_identity,
        pe_fraction = illumina_pe.params.min_mapped_frac,
        pe_single_strobe = illumina_pe.allow_single_strobe_matches,
        pe_pair_aware = illumina_pe.pair_aware,
        pe_batch_size = illumina_pe.batch_size,
        long_top_n = long_reads.params.top_n,
        long_score_ratio = long_reads.params.score_ratio,
        long_min_seeds = long_reads.params.min_seeds,
        long_identity = long_reads.params.min_identity,
        long_fraction = long_reads.params.min_mapped_frac,
        long_single_strobe = long_reads.allow_single_strobe_matches,
        long_batch_size = long_reads.batch_size,
        recall_score_ratio = illumina_recall.params.score_ratio,
        min_insert = illumina_pe.min_insert,
        max_insert = illumina_pe.max_insert,
        mate_rescue_taxa = illumina_pe.mate_rescue_taxa,
    )
}

impl MapArgOverrides {
    fn parse(args: &[String]) -> Result<Self> {
        let mut parsed = Self::default();
        let mut i = 0usize;
        while i < args.len() {
            let arg = &args[i];
            match arg.as_str() {
                "--allow-single-strobe-matches" => {
                    if parsed.allow_single_strobe_matches == Some(false) {
                        bail!("--allow-single-strobe-matches conflicts with --match-both-strobes");
                    }
                    parsed.allow_single_strobe_matches = Some(true);
                }
                "--match-both-strobes" | "--require-both-strobes" => {
                    if parsed.allow_single_strobe_matches == Some(true) {
                        bail!("--match-both-strobes conflicts with --allow-single-strobe-matches");
                    }
                    parsed.allow_single_strobe_matches = Some(false);
                }
                "--no-identity-filter" => parsed.no_identity_filter = true,
                "--ignore-n" => parsed.ignore_n = true,
                "--with-strand" | "--debug-output" => {
                    bail!(
                        "{} is not supported by metax profile because reusable mapping rows have a fixed 11-column schema",
                        arg
                    );
                }
                "--extend" => {
                    if parsed.fast_extend {
                        bail!("--extend conflicts with --fast-extend");
                    }
                    parsed.extend = true;
                }
                "--fast-extend" => {
                    if parsed.extend {
                        bail!("--fast-extend conflicts with --extend");
                    }
                    parsed.fast_extend = true;
                }
                "--pair-aware" => parsed.pair_aware = true,
                "--paired-emission" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.paired_emission = Some(parse_paired_emission_policy(value)?);
                }
                "--quiet" | "--progress" | "--no-progress" => {}
                "--preset" | "-p" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.preset = Some(parse_preset(value)?);
                }
                "--index-cache" => {
                    if i + 1 < args.len() && !args[i + 1].starts_with('-') {
                        i += 1;
                        parsed.index_cache = Some(Some(args[i].clone()));
                    } else {
                        parsed.index_cache = Some(Some("auto".to_string()));
                    }
                }
                "--top-n" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.top_n = Some(parse_positive_usize(value, arg)?);
                }
                "--score-ratio" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.score_ratio = Some(parse_probability(value, arg)?);
                }
                "--min-seeds" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.min_seeds = Some(parse_positive_usize(value, arg)?);
                }
                "--min-identity" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.min_identity = Some(parse_probability(value, arg)?);
                }
                "--min-mapped-frac" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.min_mapped_frac = Some(parse_probability(value, arg)?);
                }
                "--identity-mode" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.identity_mode = Some(parse_identity_mode(value)?);
                }
                "--query-seed-stride" | "--seed-stride" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.query_seed_stride = Some(parse_positive_usize(value, arg)?);
                }
                "--no-long-read-rescue" => parsed.no_long_read_rescue = true,
                "--max-rescue-candidates" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.long_read_rescue_limit = Some(parse_positive_usize(value, arg)?);
                }
                "--batch-size" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.batch_size = Some(parse_positive_usize(value, arg)?);
                }
                "--min-insert" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.min_insert = Some(parse_usize(value, arg)?);
                }
                "--max-insert" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.max_insert = Some(parse_usize(value, arg)?);
                }
                "--mate-rescue-taxa" => {
                    let value = next_arg_value(args, &mut i, arg)?;
                    parsed.mate_rescue_taxa = Some(parse_usize(value, arg)?);
                }
                "--progress-interval" => {
                    let _ = next_arg_value(args, &mut i, arg)?;
                }
                "--index" | "--reads" | "--reads2" | "--out" | "--profile-out" | "--threads"
                | "-i" | "-I" | "-o" | "-t" => {
                    bail!("{arg} is managed by Metax and cannot be passed after --");
                }
                _ => {
                    if let Some((name, value)) = split_long_option(arg) {
                        match name {
                            "--preset" => parsed.preset = Some(parse_preset(value)?),
                            "--index-cache" => parsed.index_cache = Some(Some(value.to_string())),
                            "--top-n" => parsed.top_n = Some(parse_positive_usize(value, name)?),
                            "--score-ratio" => parsed.score_ratio = Some(parse_probability(value, name)?),
                            "--min-seeds" => parsed.min_seeds = Some(parse_positive_usize(value, name)?),
                            "--min-identity" => parsed.min_identity = Some(parse_probability(value, name)?),
                            "--min-mapped-frac" => parsed.min_mapped_frac = Some(parse_probability(value, name)?),
                            "--identity-mode" => parsed.identity_mode = Some(parse_identity_mode(value)?),
                            "--query-seed-stride" | "--seed-stride" => {
                                parsed.query_seed_stride = Some(parse_positive_usize(value, name)?)
                            }
                            "--max-rescue-candidates" => {
                                parsed.long_read_rescue_limit =
                                    Some(parse_positive_usize(value, name)?)
                            }
                            "--paired-emission" => {
                                parsed.paired_emission = Some(parse_paired_emission_policy(value)?)
                            }
                            "--batch-size" => parsed.batch_size = Some(parse_positive_usize(value, name)?),
                            "--min-insert" => parsed.min_insert = Some(parse_usize(value, name)?),
                            "--max-insert" => parsed.max_insert = Some(parse_usize(value, name)?),
                            "--mate-rescue-taxa" => {
                                parsed.mate_rescue_taxa = Some(parse_usize(value, name)?)
                            }
                            "--progress-interval" => {}
                            "--index" | "--reads" | "--reads2" | "--out" | "--profile-out"
                            | "--threads" => bail!(
                                "{name} is managed by Metax and cannot be passed after --"
                            ),
                            _ => bail!(
                                "unsupported mapper option '{}'; run 'metax profile --mapper-help' for the complete list",
                                arg
                            ),
                        }
                    } else {
                        bail!(
                            "unsupported mapper argument '{}'; run 'metax profile --mapper-help' for the complete list",
                            arg
                        );
                    }
                }
            }
            i += 1;
        }
        Ok(parsed)
    }
}

fn read_single_batch(reader: &mut SequenceReader, config: &MappingConfig) -> Result<Vec<Record>> {
    let Some(base_limit) = config.batch_base_limit else {
        return reader.read_batch(config.batch_size);
    };

    let mut batch = Vec::with_capacity(config.batch_size);
    let mut bases = 0usize;
    while batch.len() < config.batch_size {
        let Some(record) = reader.next_record()? else {
            break;
        };
        bases = bases.saturating_add(record.seq.len());
        batch.push(record);
        // The record that crosses the limit stays in this batch. This avoids
        // buffering a look-ahead record and guarantees progress even for one
        // exceptionally long sequence.
        if bases >= base_limit {
            break;
        }
    }
    Ok(batch)
}

fn produce_single_batches(
    reads: &str,
    config: &MappingConfig,
    sender: SyncSender<Vec<Record>>,
) -> Result<()> {
    let mut reader = SequenceReader::from_path(reads)?;
    loop {
        let batch = read_single_batch(&mut reader, config)?;
        if batch.is_empty() {
            return Ok(());
        }
        if sender.send(batch).is_err() {
            // The mapping consumer already failed. Let it report the primary
            // error rather than replacing it with a channel diagnostic.
            return Ok(());
        }
    }
}

fn produce_paired_batches(
    reads1: &str,
    reads2: &str,
    batch_size: usize,
    sender: SyncSender<Vec<(Record, Record)>>,
) -> Result<()> {
    let mut reader1 = SequenceReader::from_path(reads1)?;
    let mut reader2 = SequenceReader::from_path(reads2)?;
    let mut pairs_read = 0u64;
    loop {
        let batch = read_pair_batch(&mut reader1, &mut reader2, batch_size, pairs_read)?;
        if batch.is_empty() {
            return Ok(());
        }
        pairs_read += batch.len() as u64;
        if sender.send(batch).is_err() {
            return Ok(());
        }
    }
}

fn consume_prefetched_batches<B, P, C>(producer: P, mut consumer: C) -> Result<()>
where
    B: Send,
    P: FnOnce(SyncSender<B>) -> Result<()> + Send,
    C: FnMut(B) -> Result<()>,
{
    thread::scope(|scope| {
        let (sender, receiver) = mpsc::sync_channel(PIPELINE_CHANNEL_CAPACITY);
        let producer_handle = scope.spawn(move || producer(sender));

        let consumer_result: Result<()> = (|| {
            while let Ok(batch) = receiver.recv() {
                consumer(batch)?;
            }
            Ok(())
        })();
        // If the consumer failed, disconnect a producer that may currently be
        // blocked handing off its one prefetched batch.
        drop(receiver);
        let producer_result = join_scoped_worker(producer_handle, "mapping input reader")?;
        consumer_result?;
        producer_result
    })
}

fn join_scoped_worker<T>(handle: thread::ScopedJoinHandle<'_, T>, label: &str) -> Result<T> {
    handle
        .join()
        .map_err(|_| anyhow!("{label} thread panicked"))
}

fn add_text_mapping_results(profile: &mut MapProfile, results: &[PerReadMappingOutput]) {
    for result in results {
        profile.total_reads += result.read_records;
        profile.add(&result.profile);
        profile.mapped_reads += result.mapped_read_records;
        profile.emitted_records += result.emitted_records;
    }
}

fn write_text_mapping_results(
    writer: &mut (dyn Write + Send),
    receiver: mpsc::Receiver<Vec<PerReadMappingOutput>>,
) -> Result<()> {
    while let Ok(results) = receiver.recv() {
        for result in results {
            if !result.text.is_empty() {
                writer.write_all(result.text.as_bytes())?;
            }
        }
    }
    Ok(())
}

fn run_text_mapping_pipeline<B, P, M>(
    writer: &mut (dyn Write + Send),
    producer: P,
    mut map_batch: M,
    profile: &mut MapProfile,
    progress: &mut MappingProgress,
    progress_label: &'static str,
) -> Result<()>
where
    B: Send,
    P: FnOnce(SyncSender<B>) -> Result<()> + Send,
    M: FnMut(B) -> Vec<PerReadMappingOutput>,
{
    thread::scope(|scope| {
        let (output_sender, output_receiver) = mpsc::sync_channel(PIPELINE_CHANNEL_CAPACITY);
        let writer_handle =
            scope.spawn(move || write_text_mapping_results(writer, output_receiver));

        let mapping_result = consume_prefetched_batches(producer, |batch| {
            let results = map_batch(batch);
            add_text_mapping_results(profile, &results);
            output_sender
                .send(results)
                .map_err(|_| anyhow!("mapping output writer stopped unexpectedly"))?;
            progress.tick(progress_label, profile);
            Ok(())
        });
        drop(output_sender);
        let writer_result = join_scoped_worker(writer_handle, "mapping output writer")?;

        // A concrete write error is more useful than the channel-closed error
        // observed by the mapper after the writer exits.
        writer_result?;
        mapping_result
    })
}

fn map_single(
    reads: &str,
    index: &Arc<RuntimeIndex>,
    pool: &rayon::ThreadPool,
    config: &MappingConfig,
    writer: &mut (dyn Write + Send),
    profile: &mut MapProfile,
    progress: &mut MappingProgress,
) -> Result<()> {
    run_text_mapping_pipeline(
        writer,
        |sender| produce_single_batches(reads, config, sender),
        |batch: Vec<Record>| {
            let results = pool.install(|| {
                batch
                    .par_iter()
                    .map_init(MappingWorkspace::default, |workspace, record| {
                        let mut local = MapProfile::default();
                        let mappings = map_read_with_workspace(
                            &record.id,
                            &record.seq,
                            index,
                            config.params,
                            &mut local,
                            workspace,
                        );
                        let mapped_read_records = u64::from(!mappings.is_empty());
                        let (text, emitted_records) = format_mapping_lines(mappings);
                        PerReadMappingOutput {
                            read_records: 1,
                            mapped_read_records,
                            emitted_records,
                            profile: local,
                            text,
                        }
                    })
                    .collect::<Vec<PerReadMappingOutput>>()
            });
            drop(batch);
            results
        },
        profile,
        progress,
        "Mapping reads",
    )
}

fn map_single_structured(
    reads: &str,
    index: &Arc<RuntimeIndex>,
    pool: &rayon::ThreadPool,
    config: &MappingConfig,
    sender: &SyncSender<MetamMappingBatch>,
    profile: &mut MapProfile,
    progress: &mut MappingProgress,
) -> Result<()> {
    consume_prefetched_batches(
        |batch_sender| produce_single_batches(reads, config, batch_sender),
        |batch: Vec<Record>| {
            let results = pool.install(|| {
                batch
                    .par_iter()
                    .map_init(MappingWorkspace::default, |workspace, record| {
                        let mut local = MapProfile::default();
                        let mappings = map_read_with_workspace(
                            &record.id,
                            &record.seq,
                            index,
                            config.params,
                            &mut local,
                            workspace,
                        );
                        let mapped_read_records = u64::from(!mappings.is_empty());
                        let emitted_records = mappings.len() as u64;
                        let mappings = mappings.into_iter().map(MetamMapping::from).collect();
                        PerReadStructuredMappingOutput {
                            read_records: 1,
                            mapped_read_records,
                            emitted_records,
                            profile: local,
                            mappings,
                        }
                    })
                    .collect::<Vec<PerReadStructuredMappingOutput>>()
            });
            drop(batch);
            let mut batch_mappings = Vec::with_capacity(STRUCTURED_MAPPING_CHUNK_RECORDS);
            for result in results {
                profile.total_reads += result.read_records;
                profile.add(&result.profile);
                profile.mapped_reads += result.mapped_read_records;
                profile.emitted_records += result.emitted_records;
                append_mapping_chunk(sender, &mut batch_mappings, result.mappings)?;
            }
            send_mapping_batch(sender, batch_mappings)?;
            progress.tick("Mapping reads", profile);
            Ok(())
        },
    )
}

/// Map one read pair once for both the text and structured output paths.
/// Keeping this policy in one function prevents saved mappings and direct
/// profiling from drifting apart.
fn map_record_pair(
    r1: &Record,
    r2: &Record,
    index: &RuntimeIndex,
    config: &MappingConfig,
    workspace: &mut MappingWorkspace,
) -> (Vec<Mapping>, MapProfile) {
    let mut local = MapProfile::default();
    let mut m1 =
        map_read_with_workspace(&r1.id, &r1.seq, index, config.params, &mut local, workspace);
    let mut m2 =
        map_read_with_workspace(&r2.id, &r2.seq, index, config.params, &mut local, workspace);
    // Mates recovered for taxa that compete with a seeded proper pair.
    let mut competing1 = Vec::new();
    let mut competing2 = Vec::new();

    if config.pair_aware {
        let mut rescue = MateRescue {
            index,
            config,
            profile: &mut local,
            workspace: &mut *workspace,
        };
        if m1.is_empty() && !m2.is_empty() {
            m1 = rescue.beside_any(r1, &m2);
        }
        if m2.is_empty() && !m1.is_empty() {
            m2 = rescue.beside_any(r2, &m1);
        }
        if !m1.is_empty() && !m2.is_empty() {
            if has_proper_pair(&m1, &m2, config.min_insert, config.max_insert) {
                // A proper pair on one taxon must not hide a taxon that
                // explains the fragment as well but was seeded on one mate only.
                competing1 = rescue.for_unpaired_taxa(r1, &m2, &m1);
                competing2 = rescue.for_unpaired_taxa(r2, &m1, &m2);
            } else {
                let mut rescued1 = rescue.beside_any(r1, &m2);
                let mut rescued2 = rescue.beside_any(r2, &m1);
                m1.append(&mut rescued1);
                m2.append(&mut rescued2);
            }
        }
    }

    apply_pair_awareness(
        &mut m1,
        &mut m2,
        config.min_insert,
        config.max_insert,
        config.pair_aware,
    );
    limit_mappings_with_workspace(&mut m1, config.params.top_n, workspace);
    limit_mappings_with_workspace(&mut m2, config.params.top_n, workspace);
    // Rescue scores are not on the seed-score scale, so these mates join behind
    // the retained placements instead of competing with them for top-n slots.
    m1.append(&mut competing1);
    m2.append(&mut competing2);
    set_mate_status(&mut m1, &mut m2);

    let mappings = prepare_paired_mappings(
        m1,
        m2,
        config.paired_emission,
        config.min_insert,
        config.max_insert,
    );
    (mappings, local)
}

/// The mate searches of one read pair and the context they share.
struct MateRescue<'a> {
    index: &'a RuntimeIndex,
    config: &'a MappingConfig,
    profile: &'a mut MapProfile,
    workspace: &'a mut MappingWorkspace,
}

impl MateRescue<'_> {
    /// The best placement of `read` beside any placement of its mate.
    fn beside_any(&mut self, read: &Record, mates: &[Mapping]) -> Vec<Mapping> {
        rescue_mate_with_workspace(
            &read.id,
            &read.seq,
            self.index,
            mates,
            self.config.params,
            self.config.min_insert,
            self.config.max_insert,
            self.profile,
            self.workspace,
        )
    }

    /// Placements of `read` that complete a pair, no worse than the best
    /// seeded one, for taxa that `mates` hold without a proper partner in `own`.
    fn for_unpaired_taxa(
        &mut self,
        read: &Record,
        mates: &[Mapping],
        own: &[Mapping],
    ) -> Vec<Mapping> {
        rescue_competing_mates_with_workspace(
            &read.id,
            &read.seq,
            self.index,
            mates,
            own,
            self.config.params,
            self.config.min_insert,
            self.config.max_insert,
            self.config.mate_rescue_taxa,
            self.profile,
            self.workspace,
        )
    }
}

fn map_paired(
    reads1: &str,
    reads2: &str,
    index: &Arc<RuntimeIndex>,
    pool: &rayon::ThreadPool,
    config: &MappingConfig,
    writer: &mut (dyn Write + Send),
    profile: &mut MapProfile,
    progress: &mut MappingProgress,
) -> Result<()> {
    run_text_mapping_pipeline(
        writer,
        |sender| produce_paired_batches(reads1, reads2, config.batch_size, sender),
        |batch: Vec<(Record, Record)>| {
            let results = pool.install(|| {
                batch
                    .par_iter()
                    .map_init(MappingWorkspace::default, |workspace, (r1, r2)| {
                        let (mappings, local) = map_record_pair(r1, r2, index, config, workspace);
                        let mapped_read_records = mapped_mate_count(&mappings);
                        let (text, emitted_records) = format_mapping_lines(mappings);
                        PerReadMappingOutput {
                            read_records: 2,
                            mapped_read_records,
                            emitted_records,
                            profile: local,
                            text,
                        }
                    })
                    .collect::<Vec<_>>()
            });
            drop(batch);
            results
        },
        profile,
        progress,
        "Mapping read records",
    )
}

fn map_paired_structured(
    reads1: &str,
    reads2: &str,
    index: &Arc<RuntimeIndex>,
    pool: &rayon::ThreadPool,
    config: &MappingConfig,
    sender: &SyncSender<MetamMappingBatch>,
    profile: &mut MapProfile,
    progress: &mut MappingProgress,
) -> Result<()> {
    consume_prefetched_batches(
        |batch_sender| produce_paired_batches(reads1, reads2, config.batch_size, batch_sender),
        |batch: Vec<(Record, Record)>| {
            let results = pool.install(|| {
                batch
                    .par_iter()
                    .map_init(MappingWorkspace::default, |workspace, (r1, r2)| {
                        let (mappings, local) = map_record_pair(r1, r2, index, config, workspace);
                        let mapped_read_records = mapped_mate_count(&mappings);
                        let emitted_records = mappings.len() as u64;
                        let mappings = mappings.into_iter().map(MetamMapping::from).collect();
                        PerReadStructuredMappingOutput {
                            read_records: 2,
                            mapped_read_records,
                            emitted_records,
                            profile: local,
                            mappings,
                        }
                    })
                    .collect::<Vec<_>>()
            });
            drop(batch);
            let mut batch_mappings = Vec::with_capacity(STRUCTURED_MAPPING_CHUNK_RECORDS);
            for result in results {
                profile.total_reads += result.read_records;
                profile.add(&result.profile);
                profile.mapped_reads += result.mapped_read_records;
                profile.emitted_records += result.emitted_records;
                append_mapping_chunk(sender, &mut batch_mappings, result.mappings)?;
            }
            send_mapping_batch(sender, batch_mappings)?;
            progress.tick("Mapping read records", profile);
            Ok(())
        },
    )
}

fn format_mapping_lines<I>(mappings: I) -> (String, u64)
where
    I: IntoIterator<Item = Mapping>,
{
    let mappings = mappings.into_iter();
    let (min_records, _) = mappings.size_hint();
    let mut text = String::with_capacity(min_records.saturating_mul(160));
    let mut emitted_records = 0u64;
    for mapping in mappings {
        text.push_str(&format_mapping(&mapping, false, false));
        text.push('\n');
        emitted_records += 1;
    }
    (text, emitted_records)
}

fn send_mapping_batch(
    sender: &SyncSender<MetamMappingBatch>,
    batch: MetamMappingBatch,
) -> Result<()> {
    if batch.is_empty() {
        return Ok(());
    }
    sender
        .send(batch)
        .map_err(|_| anyhow::Error::new(MappingConsumerClosed))
}

fn append_mapping_chunk(
    sender: &SyncSender<MetamMappingBatch>,
    pending: &mut MetamMappingBatch,
    mappings: MetamMappingBatch,
) -> Result<()> {
    for mapping in mappings {
        pending.push(mapping);
        if pending.len() == STRUCTURED_MAPPING_CHUNK_RECORDS {
            send_mapping_batch(sender, std::mem::take(pending))?;
            pending.reserve(STRUCTURED_MAPPING_CHUNK_RECORDS);
        }
    }
    Ok(())
}

/// Select explicit one-to-one R1/R2 placements and emit them in deterministic
/// pair order. A nonzero `pair_id` is the sole pairing contract consumed by
/// the profiler; surplus and discordant mappings retain their mate number but
/// remain diagnostic-only with `pair_id == 0`.
fn prepare_paired_mappings(
    mut mate1: Vec<Mapping>,
    mut mate2: Vec<Mapping>,
    policy: PairedEmissionPolicy,
    min_insert: usize,
    max_insert: usize,
) -> Vec<Mapping> {
    // A per-mate candidate boundary is not evidence that the paired candidate
    // set is incomplete: an omitted placement can compete in PE profiling only
    // if the other mate supports the same reference. Propagating either side's
    // boundary suppressed valid same-reference pairs. Keep the taxon-diverse
    // retained pairs, but do not turn one-sided search truncation into a
    // pair-level coverage veto.
    for mapping in &mut mate1 {
        mapping.mate_number = 1;
        mapping.pair_id = 0;
        mapping.candidate_boundary = CandidateBoundary::COMPLETE;
    }
    for mapping in &mut mate2 {
        mapping.mate_number = 2;
        mapping.pair_id = 0;
        mapping.candidate_boundary = CandidateBoundary::COMPLETE;
    }

    let total = mate1.len().saturating_add(mate2.len());
    let mut mate1 = mate1.into_iter().map(Some).collect::<Vec<_>>();
    let mut mate2 = mate2.into_iter().map(Some).collect::<Vec<_>>();
    let mut paired = Vec::with_capacity(total);

    // R1 mappings are already ordered by pair-adjusted mapper score. For each
    // one, choose the strongest still-unpaired R2 placement on the same
    // reference. This is O(n*m), allocation-free beyond the output vectors,
    // and bounded by the mapper's top_n candidate cap.
    for left_slot in &mut mate1 {
        let best_partner = if let Some(left) = left_slot
            .as_ref()
            .filter(|mapping| mapping.mate_ref_status.is_same())
        {
            let mut best = None;
            for (mate2_index, right) in mate2.iter().enumerate() {
                let Some(right) = right.as_ref() else {
                    continue;
                };
                if !right.mate_ref_status.is_same() || right.ref_id != left.ref_id {
                    continue;
                }
                let compatibility = pair_compatibility(left, right, min_insert, max_insert);
                if partner_is_better(right, mate2_index, compatibility.score, best, &mate2) {
                    best = Some((mate2_index, compatibility));
                }
            }
            best
        } else {
            None
        };

        let Some((mate2_index, compatibility)) = best_partner else {
            continue;
        };
        let pair_id = u32::try_from(paired.len() / 2 + 1)
            .expect("paired mapping count exceeds the u32 output contract");
        let mut left = left_slot.take().expect("selected R1 mapping is present");
        let mut right = mate2[mate2_index]
            .take()
            .expect("selected R2 mapping is present");
        for mapping in [&mut left, &mut right] {
            mapping.pair_id = pair_id;
            mapping.pair_score = compatibility.score.max(0);
            mapping.template_length = compatibility.template_length;
            mapping.pair_status = PairStatus::Paired(compatibility.relation);
        }
        paired.push(left);
        paired.push(right);
    }

    if !paired.is_empty() && policy != PairedEmissionPolicy::All {
        return paired;
    }
    if paired.is_empty() && policy == PairedEmissionPolicy::ConcordantOnly {
        return paired;
    }

    // `all`, or `concordant-first` with no selected pair: retain diagnostics
    // in stable mate-rank order after any explicit pairs.
    let remaining = mate1.len().max(mate2.len());
    for index in 0..remaining {
        if let Some(mapping) = mate1.get_mut(index).and_then(Option::take) {
            paired.push(mapping);
        }
        if let Some(mapping) = mate2.get_mut(index).and_then(Option::take) {
            paired.push(mapping);
        }
    }
    paired
}

fn partner_is_better(
    candidate: &Mapping,
    candidate_index: usize,
    candidate_pair_score: i64,
    current: Option<(usize, metam::mapper::PairCompatibility)>,
    mate2: &[Option<Mapping>],
) -> bool {
    let Some((current_index, current_compatibility)) = current else {
        return true;
    };
    let current = mate2[current_index]
        .as_ref()
        .expect("current R2 partner is present");
    candidate_pair_score
        .cmp(&current_compatibility.score)
        .then_with(|| candidate.adjusted_score.cmp(&current.adjusted_score))
        .then_with(|| candidate.score.cmp(&current.score))
        .then_with(|| candidate.mapped_bases.cmp(&current.mapped_bases))
        .then_with(|| candidate.identity.total_cmp(&current.identity))
        .then_with(|| current_index.cmp(&candidate_index))
        == Ordering::Greater
}

fn mapped_mate_count(mappings: &[Mapping]) -> u64 {
    let mut seen_mate1 = false;
    let mut seen_mate2 = false;
    for mapping in mappings {
        seen_mate1 |= mapping.mate_number == 1;
        seen_mate2 |= mapping.mate_number == 2;
        if seen_mate1 && seen_mate2 {
            return 2;
        }
    }
    u64::from(seen_mate1) + u64::from(seen_mate2)
}

fn read_pair_batch(
    reader1: &mut SequenceReader,
    reader2: &mut SequenceReader,
    batch_size: usize,
    pairs_already_read: u64,
) -> Result<Vec<(Record, Record)>> {
    let mut batch = Vec::with_capacity(batch_size);
    for batch_index in 0..batch_size {
        match (reader1.next_record()?, reader2.next_record()?) {
            (Some(r1), Some(r2)) => {
                validate_pair_ids(&r1, &r2, pairs_already_read + batch_index as u64 + 1)?;
                batch.push((r1, r2));
            }
            (None, None) => break,
            (Some(_), None) => bail!("R1 has more records than R2"),
            (None, Some(_)) => bail!("R2 has more records than R1"),
        }
    }
    Ok(batch)
}

#[inline]
fn validate_pair_ids(r1: &Record, r2: &Record, pair_number: u64) -> Result<()> {
    if r1.id != r2.id {
        bail!(
            "paired FASTQ record {} has mismatched IDs after normalizing mate suffixes: R1='{}', R2='{}'",
            pair_number,
            r1.id,
            r2.id
        );
    }
    Ok(())
}

fn validate_mapping_config(config: &MappingConfig, reads2: Option<&str>) -> Result<()> {
    if !(0.0..=1.0).contains(&config.params.score_ratio) {
        bail!("--score-ratio must be between 0.0 and 1.0");
    }
    if !(0.0..=1.0).contains(&config.params.min_identity) {
        bail!("--min-identity must be between 0.0 and 1.0");
    }
    if !(0.0..=1.0).contains(&config.params.min_mapped_frac) {
        bail!("--min-mapped-frac must be between 0.0 and 1.0");
    }
    if config.params.top_n == 0 {
        bail!("--top-n must be greater than 0");
    }
    if config.params.min_seeds == 0 {
        bail!("--min-seeds must be greater than 0");
    }
    if config.batch_size == 0 {
        bail!("--batch-size must be greater than 0");
    }
    if config.min_insert > config.max_insert {
        bail!("--min-insert must be <= --max-insert");
    }
    if config.preset == MappingPreset::Ont {
        if reads2.is_some() {
            bail!("ONT/PacBio mapping preset is single-end only");
        }
        if config.pair_aware {
            bail!("ONT/PacBio mapping preset does not support --pair-aware");
        }
    } else if matches!(
        config.params.long_read_extension,
        LongReadExtension::Candidate | LongReadExtension::Retained
    ) {
        bail!("--extend and --fast-extend are only supported for ONT/PacBio mapping");
    }
    Ok(())
}

fn map_params_from_preset(
    preset: MappingPreset,
    top_n: Option<usize>,
    score_ratio: Option<f64>,
    min_seeds: Option<usize>,
    min_identity: Option<f64>,
    min_mapped_frac: Option<f64>,
    no_identity_filter: bool,
    ignore_n: bool,
    identity_mode: IdentityMode,
    extend: bool,
    fast_extend: bool,
    query_seed_stride: Option<usize>,
    long_read_rescue: bool,
    long_read_rescue_limit: Option<usize>,
) -> MapParams {
    let mut params = match preset {
        MappingPreset::Illumina => MapParams {
            top_n: 5,
            score_ratio: 0.90,
            min_seeds: 2,
            min_identity: 0.80,
            min_mapped_frac: 0.50,
            min_mapped_len: None,
            no_identity_filter,
            ignore_n,
            identity_mode,
            diag_bin_size: 1,
            extension_max_mismatches: 4,
            full_read_interval: false,
            fast_long_read_identity: false,
            long_read_chaining: false,
            long_read_extension: LongReadExtension::None,
            query_seed_stride: 1,
            long_read_rescue: false,
            long_read_rescue_limit: None,
        },
        MappingPreset::Ont => MapParams {
            top_n: 10,
            // Retaining close long-read candidates consistently removed false
            // species calls in independent sample-0/sample-1 validation
            // without losing true calls. Profile modes still differ through
            // their downstream evidence thresholds.
            score_ratio: 0.70,
            min_seeds: 3,
            min_identity: 0.80,
            min_mapped_frac: 0.30,
            min_mapped_len: Some(300),
            no_identity_filter,
            ignore_n,
            identity_mode,
            diag_bin_size: 128,
            extension_max_mismatches: 40,
            full_read_interval: true,
            fast_long_read_identity: true,
            long_read_chaining: true,
            long_read_extension: if extend {
                LongReadExtension::Candidate
            } else if fast_extend {
                LongReadExtension::Retained
            } else {
                LongReadExtension::None
            },
            query_seed_stride: 1,
            long_read_rescue,
            long_read_rescue_limit,
        },
    };
    if let Some(value) = top_n {
        params.top_n = value;
    }
    if let Some(value) = score_ratio {
        params.score_ratio = value;
    }
    if let Some(value) = min_seeds {
        params.min_seeds = value;
    }
    if let Some(value) = min_identity {
        params.min_identity = value;
    }
    if let Some(value) = min_mapped_frac {
        params.min_mapped_frac = value;
    }
    if let Some(value) = query_seed_stride {
        params.query_seed_stride = value;
    }
    params
}

struct MetamInputs {
    reads1: String,
    reads2: Option<String>,
}

fn default_batch_size_for_preset(preset: MappingPreset) -> usize {
    match preset {
        MappingPreset::Illumina => DEFAULT_ILLUMINA_BATCH_SIZE,
        MappingPreset::Ont => DEFAULT_LONG_READ_BATCH_SIZE,
    }
}

pub fn existing_mapping_output(outprefix: &str) -> Option<MappingOutput> {
    let mappings_tsv = mapping_tsv_gz_path(outprefix);
    if Path::new(&mappings_tsv).exists() {
        return Some(MappingOutput::saved(PathBuf::from(mappings_tsv)));
    }
    None
}

fn mapping_tsv_gz_path(outprefix: &str) -> String {
    format!("{outprefix}{MAPPING_TSV_GZ_SUFFIX}")
}

fn index_cache_path(index_path: &str, cache_arg: Option<&str>) -> Option<String> {
    cache_arg.map(|cache| {
        if cache == "auto" {
            format!("{index_path}.metamcache")
        } else {
            cache.to_string()
        }
    })
}

fn parse_preset(value: &str) -> Result<MappingPreset> {
    match value.to_ascii_lowercase().as_str() {
        "illumina" => Ok(MappingPreset::Illumina),
        "ont" | "nanopore" | "pacbio" => Ok(MappingPreset::Ont),
        _ => bail!("unsupported --preset '{}'; expected illumina or ont", value),
    }
}

fn parse_identity_mode(value: &str) -> Result<IdentityMode> {
    match value.to_ascii_lowercase().as_str() {
        "gapped" => Ok(IdentityMode::Gapped),
        "ungapped" => Ok(IdentityMode::Ungapped),
        _ => bail!(
            "unsupported --identity-mode '{}'; expected gapped or ungapped",
            value
        ),
    }
}

fn parse_paired_emission_policy(value: &str) -> Result<PairedEmissionPolicy> {
    match value.to_ascii_lowercase().as_str() {
        "concordant-first" => Ok(PairedEmissionPolicy::ConcordantFirst),
        "all" => Ok(PairedEmissionPolicy::All),
        "concordant-only" => Ok(PairedEmissionPolicy::ConcordantOnly),
        _ => bail!(
            "unsupported --paired-emission '{}'; expected concordant-first, all, or concordant-only",
            value
        ),
    }
}

fn parse_probability(value: &str, option: &str) -> Result<f64> {
    let parsed: f64 = value
        .parse()
        .with_context(|| format!("{option} value '{value}' is not a number"))?;
    if !parsed.is_finite() || !(0.0..=1.0).contains(&parsed) {
        bail!("{option} must be within [0.0, 1.0], got {parsed}");
    }
    Ok(parsed)
}

fn parse_positive_usize(value: &str, option: &str) -> Result<usize> {
    let parsed = parse_usize(value, option)?;
    if parsed == 0 {
        bail!("{option} must be greater than zero");
    }
    Ok(parsed)
}

fn parse_usize(value: &str, option: &str) -> Result<usize> {
    value
        .parse()
        .with_context(|| format!("{option} value '{value}' is not an unsigned integer"))
}

fn next_arg_value<'a>(args: &'a [String], i: &mut usize, option: &str) -> Result<&'a str> {
    *i += 1;
    args.get(*i)
        .map(String::as_str)
        .ok_or_else(|| anyhow!("{option} requires a value"))
}

fn split_long_option(arg: &str) -> Option<(&str, &str)> {
    let (name, value) = arg.split_once('=')?;
    if name.starts_with("--") {
        Some((name, value))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use metam::mapper::MateRefStatus;

    struct IntentionalWriteFailure;

    impl Write for IntentionalWriteFailure {
        fn write(&mut self, _buffer: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("intentional write failure"))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn prefetched_batches_preserve_order_and_surface_reader_errors() {
        let mut observed = Vec::new();
        consume_prefetched_batches(
            |sender| {
                sender.send(vec![1u8, 2]).unwrap();
                sender.send(vec![3u8, 4]).unwrap();
                Ok(())
            },
            |batch| {
                observed.extend(batch);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(observed, [1, 2, 3, 4]);

        let mut observed_before_error = Vec::new();
        let error = consume_prefetched_batches(
            |sender| {
                sender.send(vec![5u8]).unwrap();
                bail!("intentional input failure")
            },
            |batch| {
                observed_before_error.extend(batch);
                Ok(())
            },
        )
        .unwrap_err();
        assert_eq!(observed_before_error, [5]);
        assert!(error.to_string().contains("intentional input failure"));
    }

    #[test]
    fn text_mapping_pipeline_preserves_order_and_profile_counts() {
        let mut writer = Vec::new();
        let mut profile = MapProfile::default();
        let mut progress = MappingProgress::new(
            ProgressOptions::disabled(),
            MappingProgressLog::streamed_default(false),
        );

        run_text_mapping_pipeline(
            &mut writer,
            |sender| {
                sender.send(vec![1u8, 2]).unwrap();
                sender.send(vec![3u8, 4]).unwrap();
                Ok(())
            },
            |batch| {
                batch
                    .into_iter()
                    .map(|value| PerReadMappingOutput {
                        read_records: 1,
                        mapped_read_records: 1,
                        emitted_records: 1,
                        profile: MapProfile::default(),
                        text: format!("{value}\n"),
                    })
                    .collect()
            },
            &mut profile,
            &mut progress,
            "Mapping reads",
        )
        .unwrap();

        assert_eq!(writer, b"1\n2\n3\n4\n");
        assert_eq!(profile.total_reads, 4);
        assert_eq!(profile.mapped_reads, 4);
        assert_eq!(profile.emitted_records, 4);
    }

    #[test]
    fn text_mapping_pipeline_surfaces_writer_error_without_deadlock() {
        let mut writer = IntentionalWriteFailure;
        let mut profile = MapProfile::default();
        let mut progress = MappingProgress::new(
            ProgressOptions::disabled(),
            MappingProgressLog::streamed_default(false),
        );

        let error = run_text_mapping_pipeline(
            &mut writer,
            |sender| {
                for value in 1u8..=3 {
                    if sender.send(vec![value]).is_err() {
                        return Ok(());
                    }
                }
                Ok(())
            },
            |batch| {
                batch
                    .into_iter()
                    .map(|value| PerReadMappingOutput {
                        read_records: 1,
                        mapped_read_records: 1,
                        emitted_records: 1,
                        profile: MapProfile::default(),
                        text: format!("{value}\n"),
                    })
                    .collect()
            },
            &mut profile,
            &mut progress,
            "Mapping reads",
        )
        .unwrap_err();

        assert!(error.to_string().contains("intentional write failure"));
    }

    #[test]
    fn temporary_mapping_output_is_removed_after_success() {
        let directory = tempfile::tempdir().unwrap();
        let temporary = Builder::new()
            .suffix(".tsv.gz")
            .tempfile_in(directory.path())
            .unwrap();
        let path = temporary.path().to_path_buf();
        let output = MappingOutput::temporary(temporary.into_temp_path());
        assert!(path.exists());
        output.cleanup_after_success().unwrap();
        assert!(!path.exists());
    }

    #[test]
    fn temporary_mapping_output_can_be_preserved_after_failure() {
        let directory = tempfile::tempdir().unwrap();
        let temporary = Builder::new()
            .suffix(".tsv.gz")
            .tempfile_in(directory.path())
            .unwrap();
        let path = temporary.path().to_path_buf();
        let output = MappingOutput::temporary(temporary.into_temp_path());
        let preserved = output.preserve_after_failure().unwrap().unwrap();
        assert_eq!(preserved, path);
        assert!(preserved.exists());
        std::fs::remove_file(preserved).unwrap();
    }

    #[test]
    fn human_compact_scales_by_magnitude() {
        assert_eq!(human_compact(0), "0");
        assert_eq!(human_compact(945), "945");
        assert_eq!(human_compact(1_000), "1k");
        assert_eq!(human_compact(12_345), "12.3k");
        assert_eq!(human_compact(988_000), "988k");
        assert_eq!(human_compact(1_234_567), "1.23M");
        assert_eq!(human_compact(2_000_000), "2M");
        assert_eq!(human_compact(4_560_000_000), "4.56B");
    }

    #[test]
    fn mapping_unit_is_the_label_tail() {
        assert_eq!(mapping_unit("Mapping reads"), "reads");
        assert_eq!(mapping_unit("Mapping read records"), "records");
    }

    #[test]
    fn format_mapping_progress_lists_all_counters() {
        let profile = MapProfile {
            total_reads: 1_234_567,
            mapped_reads: 988_000,
            emitted_records: 3_460_000,
            ..MapProfile::default()
        };
        let line = format_mapping_progress("reads", &profile, Instant::now());
        // Rate is timing-dependent, so pin only the counter prefix and unit suffix.
        assert!(
            line.starts_with("1.23M reads · 988k mapped · 3.46M mappings · "),
            "unexpected line: {line}"
        );
        assert!(line.ends_with(" reads/s"), "unexpected line: {line}");
    }

    #[test]
    fn illumina_uses_large_default_mapping_batch() {
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &[]).unwrap();
        assert_eq!(config.batch_size, DEFAULT_ILLUMINA_BATCH_SIZE);
        assert_eq!(config.batch_base_limit, None);
        assert_eq!(config.params.min_seeds, 1);
        assert_eq!(config.params.top_n, 50);
        assert!((config.params.score_ratio - 0.60).abs() < f64::EPSILON);
        assert!(config.allow_single_strobe_matches);
        assert_eq!(config.paired_emission, PairedEmissionPolicy::ConcordantOnly);
    }

    #[test]
    fn illumina_paired_candidate_defaults_are_mode_independent() {
        for mode in [
            ProfileMode::Recall,
            ProfileMode::Default,
            ProfileMode::Precision,
        ] {
            let config = MappingConfig::resolve(Sequencer::Illumina, true, mode, &[]).unwrap();
            assert_eq!(config.params.top_n, 50);
            assert!((config.params.score_ratio - 0.60).abs() < f64::EPSILON);
        }
    }

    #[test]
    fn illumina_single_end_keeps_two_seed_default() {
        let config =
            MappingConfig::resolve(Sequencer::Illumina, false, ProfileMode::Default, &[]).unwrap();
        assert_eq!(config.params.min_seeds, 2);
        assert_eq!(config.params.top_n, 5);
        assert!((config.params.score_ratio - 0.90).abs() < f64::EPSILON);
        assert_eq!(config.paired_emission, PairedEmissionPolicy::All);
    }

    #[test]
    fn precision_keeps_default_candidate_score_ratios_for_every_input_type() {
        for (sequencer, paired, expected) in [
            (Sequencer::Illumina, false, 0.90),
            (Sequencer::Illumina, true, 0.60),
            (Sequencer::Nanopore, false, 0.70),
            (Sequencer::PacBio, false, 0.70),
            (Sequencer::Assembly, false, 0.90),
        ] {
            let default =
                MappingConfig::resolve(sequencer, paired, ProfileMode::Default, &[]).unwrap();
            let precision =
                MappingConfig::resolve(sequencer, paired, ProfileMode::Precision, &[]).unwrap();
            assert!((default.params.score_ratio - expected).abs() < f64::EPSILON);
            assert_eq!(precision.params.score_ratio, default.params.score_ratio);
        }
    }

    #[test]
    fn recall_widens_unpaired_non_long_read_candidate_ratios() {
        for (sequencer, expected) in [(Sequencer::Illumina, 0.80), (Sequencer::Assembly, 0.80)] {
            let config =
                MappingConfig::resolve(sequencer, false, ProfileMode::Recall, &[]).unwrap();
            assert!((config.params.score_ratio - expected).abs() < f64::EPSILON);
        }
    }

    #[test]
    fn long_read_candidate_ratio_is_mode_independent() {
        for sequencer in [Sequencer::Nanopore, Sequencer::PacBio] {
            for mode in [
                ProfileMode::Recall,
                ProfileMode::Default,
                ProfileMode::Precision,
            ] {
                let config = MappingConfig::resolve(sequencer, false, mode, &[]).unwrap();
                assert!((config.params.score_ratio - 0.70).abs() < f64::EPSILON);
            }
        }
    }

    #[test]
    fn explicit_candidate_options_override_illumina_paired_defaults() {
        let args = vec!["--top-n=7".to_string(), "--score-ratio=0.75".to_string()];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &args).unwrap();
        assert_eq!(config.params.top_n, 7);
        assert!((config.params.score_ratio - 0.75).abs() < f64::EPSILON);
    }

    #[test]
    fn explicit_min_seeds_overrides_illumina_paired_default() {
        let args = vec!["--min-seeds".to_string(), "2".to_string()];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &args).unwrap();
        assert_eq!(config.params.min_seeds, 2);
    }

    #[test]
    fn long_reads_use_large_count_batch_with_byte_safety_limit() {
        let config =
            MappingConfig::resolve(Sequencer::Nanopore, false, ProfileMode::Default, &[]).unwrap();
        assert_eq!(config.batch_size, DEFAULT_LONG_READ_BATCH_SIZE);
        assert_eq!(
            config.batch_base_limit,
            Some(DEFAULT_LONG_READ_BATCH_BASE_LIMIT)
        );
        assert!(!config.allow_single_strobe_matches);
        assert_eq!(config.params.min_seeds, 3);
        // The mapper keeps the conservative 0.30 fraction gate, while its
        // 300-base length rescue exposes substantial partial mappings to the
        // profiler's length-OR-fraction policy.
        assert!((config.params.min_mapped_frac - 0.30).abs() < f64::EPSILON);
        assert_eq!(config.params.min_mapped_len, Some(300));
        assert_eq!(config.params.query_seed_stride, 1);
        assert!(config.params.long_read_rescue);
        assert_eq!(config.params.long_read_rescue_limit, None);
    }

    #[test]
    fn long_read_speed_options_are_parsed() {
        let args = vec![
            "--query-seed-stride".to_string(),
            "2".to_string(),
            "--max-rescue-candidates=64".to_string(),
            "--batch-size".to_string(),
            "32768".to_string(),
        ];
        let config =
            MappingConfig::resolve(Sequencer::Nanopore, false, ProfileMode::Default, &args)
                .unwrap();
        assert_eq!(config.params.query_seed_stride, 2);
        assert!(config.params.long_read_rescue);
        assert_eq!(config.params.long_read_rescue_limit, Some(64));
        assert_eq!(config.batch_size, 32_768);

        let args = vec![
            "--seed-stride=3".to_string(),
            "--no-long-read-rescue".to_string(),
        ];
        let config =
            MappingConfig::resolve(Sequencer::PacBio, false, ProfileMode::Default, &args).unwrap();
        assert_eq!(config.params.query_seed_stride, 3);
        assert!(!config.params.long_read_rescue);
    }

    #[test]
    fn illumina_single_strobe_default_can_be_disabled() {
        let args = vec!["--match-both-strobes".to_string()];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Recall, &args).unwrap();
        assert!(!config.allow_single_strobe_matches);

        let args = vec!["--require-both-strobes".to_string()];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, false, ProfileMode::Default, &args)
                .unwrap();
        assert!(!config.allow_single_strobe_matches);
    }

    #[test]
    fn conflicting_strobe_match_options_are_rejected() {
        let args = vec![
            "--allow-single-strobe-matches".to_string(),
            "--match-both-strobes".to_string(),
        ];
        let err = MappingConfig::resolve(Sequencer::Illumina, false, ProfileMode::Default, &args)
            .expect_err("conflicting strobe options must fail");
        assert!(err.to_string().contains("conflicts"));
    }

    #[test]
    fn extended_mapping_columns_are_rejected_before_mapping() {
        for option in ["--with-strand", "--debug-output"] {
            let error = MappingConfig::resolve(
                Sequencer::Illumina,
                false,
                ProfileMode::Default,
                &[option.to_string()],
            )
            .expect_err("extended output must be rejected by the profiling workflow");
            let message = format!("{error:#}");
            assert!(message.contains(option));
            assert!(message.contains("fixed 11-column schema"));
        }
    }

    #[test]
    fn explicit_mapping_batch_size_overrides_preset_default() {
        let args = vec!["--batch-size".to_string(), "131072".to_string()];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &args).unwrap();
        assert_eq!(config.batch_size, 131_072);
    }

    #[test]
    fn paired_emission_option_accepts_supported_modes() {
        let args = vec![
            "--paired-emission".to_string(),
            "all".to_string(),
            "--top-n=7".to_string(),
        ];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &args).unwrap();
        assert_eq!(config.paired_emission, PairedEmissionPolicy::All);
        assert_eq!(config.params.top_n, 7);

        let args = vec!["--paired-emission=concordant-only".to_string()];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &args).unwrap();
        assert_eq!(config.paired_emission, PairedEmissionPolicy::ConcordantOnly);

        let args = vec!["--paired-emission=concordant-first".to_string()];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &args).unwrap();
        assert_eq!(
            config.paired_emission,
            PairedEmissionPolicy::ConcordantFirst
        );
    }

    #[test]
    fn mate_rescue_taxa_has_a_default_and_can_be_overridden() {
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &[]).unwrap();
        assert_eq!(config.mate_rescue_taxa, DEFAULT_MATE_RESCUE_TAXA);

        let args = vec!["--mate-rescue-taxa".to_string(), "0".to_string()];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &args).unwrap();
        assert_eq!(config.mate_rescue_taxa, 0);

        let args = vec!["--mate-rescue-taxa=3".to_string()];
        let config =
            MappingConfig::resolve(Sequencer::Illumina, true, ProfileMode::Default, &args).unwrap();
        assert_eq!(config.mate_rescue_taxa, 3);
    }

    #[test]
    fn paired_record_ids_must_match_after_reader_normalization() {
        let r1 = Record {
            id: "read-42".to_string(),
            seq: b"ACGT".to_vec(),
        };
        let matching_r2 = Record {
            id: "read-42".to_string(),
            seq: b"TGCA".to_vec(),
        };
        validate_pair_ids(&r1, &matching_r2, 7).unwrap();

        let mismatching_r2 = Record {
            id: "read-43".to_string(),
            seq: b"TGCA".to_vec(),
        };
        let error = validate_pair_ids(&r1, &mismatching_r2, 7).unwrap_err();
        let message = error.to_string();
        assert!(message.contains("record 7"));
        assert!(message.contains("R1='read-42'"));
        assert!(message.contains("R2='read-43'"));
    }

    #[test]
    fn closed_structured_consumer_has_a_typed_error() {
        let (sender, receiver) = mpsc::sync_channel(1);
        drop(receiver);
        let mapping = MetamMapping::from(test_mapping("read", "ref", MateRefStatus::NotApplicable));

        let err = send_mapping_batch(&sender, vec![mapping]).unwrap_err();
        assert!(is_mapping_consumer_closed(&err));
    }

    #[test]
    fn concordant_first_emits_only_concordant_when_present() {
        let mate1 = vec![
            test_mapping("read/1", "ref-a", MateRefStatus::Same),
            test_mapping("read/1", "ref-b", MateRefStatus::Other(Arc::from("ref-c"))),
        ];
        let mate2 = vec![
            test_mapping("read/2", "ref-a", MateRefStatus::Same),
            test_mapping("read/2", "ref-c", MateRefStatus::Other(Arc::from("ref-b"))),
        ];

        let mappings = prepare_paired_mappings(
            mate1,
            mate2,
            PairedEmissionPolicy::ConcordantFirst,
            0,
            1_000,
        );

        assert_eq!(mappings.len(), 2);
        assert_eq!(mappings[0].ref_id, "ref-a");
        assert_eq!(mappings[1].ref_id, "ref-a");
        assert_eq!(mappings[0].mate_number, 1);
        assert_eq!(mappings[1].mate_number, 2);
        assert_eq!(mappings[0].pair_id, 1);
        assert_eq!(mappings[1].pair_id, 1);
    }

    #[test]
    fn concordant_first_keeps_discordant_when_no_concordant_exists() {
        let mate1 = vec![test_mapping("read/1", "ref-a", MateRefStatus::Other(Arc::from("ref-c")))];
        let mate2 = vec![test_mapping("read/2", "ref-c", MateRefStatus::Other(Arc::from("ref-a")))];

        let mappings = prepare_paired_mappings(
            mate1,
            mate2,
            PairedEmissionPolicy::ConcordantFirst,
            0,
            1_000,
        );

        assert_eq!(mappings.len(), 2);
        assert_eq!(mappings[0].ref_id, "ref-a");
        assert_eq!(mappings[1].ref_id, "ref-c");
        assert!(mappings.iter().all(|mapping| mapping.pair_id == 0));
    }

    #[test]
    fn concordant_only_drops_discordant_pairs() {
        let mate1 = vec![test_mapping("read/1", "ref-a", MateRefStatus::Other(Arc::from("ref-c")))];
        let mate2 = vec![test_mapping("read/2", "ref-c", MateRefStatus::Other(Arc::from("ref-a")))];

        let mappings =
            prepare_paired_mappings(mate1, mate2, PairedEmissionPolicy::ConcordantOnly, 0, 1_000);

        assert!(mappings.is_empty());
    }

    #[test]
    fn concordant_only_keeps_concordant_and_drops_other_rows() {
        let mate1 = vec![
            test_mapping("read/1", "ref-a", MateRefStatus::Same),
            test_mapping("read/1", "ref-b", MateRefStatus::Other(Arc::from("ref-c"))),
        ];
        let mate2 = vec![
            test_mapping("read/2", "ref-a", MateRefStatus::Same),
            test_mapping("read/2", "ref-c", MateRefStatus::Other(Arc::from("ref-b"))),
        ];

        let mappings =
            prepare_paired_mappings(mate1, mate2, PairedEmissionPolicy::ConcordantOnly, 0, 1_000);

        assert_eq!(mappings.len(), 2);
        assert_eq!(mappings[0].ref_id, "ref-a");
        assert_eq!(mappings[1].ref_id, "ref-a");
    }

    #[test]
    fn all_emission_keeps_discordant_rows_for_diagnostics() {
        let mate1 = vec![test_mapping("read/1", "ref-a", MateRefStatus::Other(Arc::from("ref-c")))];
        let mate2 = vec![test_mapping("read/2", "ref-c", MateRefStatus::Other(Arc::from("ref-a")))];

        let mappings = prepare_paired_mappings(mate1, mate2, PairedEmissionPolicy::All, 0, 1_000);

        assert_eq!(mappings.len(), 2);
        assert_eq!(mappings[0].ref_id, "ref-a");
        assert_eq!(mappings[1].ref_id, "ref-c");
        assert!(mappings.iter().all(|mapping| mapping.pair_id == 0));
    }

    #[test]
    fn explicit_pairs_are_adjacent_and_identified() {
        let mate1 = vec![
            test_mapping("read/1", "ref-a", MateRefStatus::Same),
            test_mapping("read/1", "ref-b", MateRefStatus::Same),
        ];
        let mate2 = vec![
            test_mapping("read/2", "ref-a", MateRefStatus::Same),
            test_mapping("read/2", "ref-b", MateRefStatus::Same),
        ];

        let ordered = prepare_paired_mappings(mate1, mate2, PairedEmissionPolicy::All, 0, 1_000);

        assert_eq!(
            read_ref_order(&ordered),
            vec![
                ("read/1".to_string(), "ref-a".to_string()),
                ("read/2".to_string(), "ref-a".to_string()),
                ("read/1".to_string(), "ref-b".to_string()),
                ("read/2".to_string(), "ref-b".to_string()),
            ]
        );
        assert_eq!(
            ordered
                .iter()
                .map(|mapping| (mapping.mate_number, mapping.pair_id))
                .collect::<Vec<_>>(),
            vec![(1, 1), (2, 1), (1, 2), (2, 2)]
        );
    }

    #[test]
    fn all_emission_puts_selected_pairs_before_diagnostics() {
        let mate1 = vec![
            test_mapping("read/1", "ref-a", MateRefStatus::Same),
            test_mapping("read/1", "ref-b", MateRefStatus::Other(Arc::from("ref-c"))),
        ];
        let mate2 = vec![
            test_mapping("read/2", "ref-a", MateRefStatus::Same),
            test_mapping("read/2", "ref-c", MateRefStatus::Other(Arc::from("ref-b"))),
        ];

        let ordered = prepare_paired_mappings(mate1, mate2, PairedEmissionPolicy::All, 0, 1_000);

        assert_eq!(
            read_ref_order(&ordered),
            vec![
                ("read/1".to_string(), "ref-a".to_string()),
                ("read/2".to_string(), "ref-a".to_string()),
                ("read/1".to_string(), "ref-b".to_string()),
                ("read/2".to_string(), "ref-c".to_string()),
            ]
        );
        assert_eq!(ordered[0].pair_id, 1);
        assert_eq!(ordered[1].pair_id, 1);
        assert_eq!(ordered[2].pair_id, 0);
        assert_eq!(ordered[3].pair_id, 0);
    }

    #[test]
    fn one_r2_mapping_cannot_be_fabricated_into_multiple_pairs() {
        let mate1 = vec![
            test_mapping("read/1", "ref-a", MateRefStatus::Same),
            test_mapping("read/1", "ref-a", MateRefStatus::Same),
            test_mapping("read/1", "ref-a", MateRefStatus::Same),
        ];
        let mate2 = vec![test_mapping("read/2", "ref-a", MateRefStatus::Same)];

        let concordant = prepare_paired_mappings(
            mate1.clone(),
            mate2.clone(),
            PairedEmissionPolicy::ConcordantOnly,
            0,
            1_000,
        );
        assert_eq!(concordant.len(), 2);
        assert_eq!(concordant[0].mate_number, 1);
        assert_eq!(concordant[1].mate_number, 2);
        assert_eq!(concordant[0].pair_id, concordant[1].pair_id);

        let all = prepare_paired_mappings(mate1, mate2, PairedEmissionPolicy::All, 0, 1_000);
        assert_eq!(all.len(), 4);
        assert_eq!(all.iter().filter(|mapping| mapping.pair_id != 0).count(), 2);
        assert_eq!(all.iter().filter(|mapping| mapping.pair_id == 0).count(), 2);
    }

    #[test]
    fn pair_selection_clears_per_mate_candidate_boundaries() {
        let mut mate1 = test_mapping("read/1", "ref-a", MateRefStatus::Same);
        mate1.candidate_boundary = CandidateBoundary::from_bits(1).unwrap();
        let mut mate2 = test_mapping("read/2", "ref-a", MateRefStatus::Same);
        mate2.candidate_boundary = CandidateBoundary::from_bits(2).unwrap();

        let mappings = prepare_paired_mappings(
            vec![mate1],
            vec![mate2],
            PairedEmissionPolicy::ConcordantOnly,
            0,
            1_000,
        );

        assert_eq!(mappings.len(), 2);
        assert!(mappings
            .iter()
            .all(|mapping| mapping.candidate_boundary.is_complete()));
    }

    fn read_ref_order(mappings: &[Mapping]) -> Vec<(String, String)> {
        mappings
            .iter()
            .map(|mapping| (mapping.read_id.clone(), mapping.ref_id.clone()))
            .collect()
    }

    fn test_mapping(read_id: &str, ref_id: &str, mate_ref_status: MateRefStatus) -> Mapping {
        Mapping {
            read_id: read_id.to_string(),
            ref_id: ref_id.to_string(),
            ref_start: 1,
            ref_end: 101,
            identity: 0.99,
            read_len: 100,
            mapped_bases: 100,
            mate_ref_status,
            strand: 1,
            seed_count: 2,
            score: 100,
            second_score: 90,
            rank: 1,
            num_mappings: 1,
            matched_bases: 99,
            query_start: 0,
            query_end: 100,
            pair_score: 0,
            adjusted_score: 100,
            template_length: 0,
            pair_status: PairStatus::Single,
            mate_number: 0,
            pair_id: 0,
            candidate_groups: ReferenceCandidateGroups::new(1, 1),
            candidate_boundary: CandidateBoundary::COMPLETE,
        }
    }
}

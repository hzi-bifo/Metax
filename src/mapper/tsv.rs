use anyhow::{anyhow, bail, Context, Result};
use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;
use std::sync::mpsc::Receiver;
use std::sync::{Arc, Mutex};

use metam::format::canonical_mapping_identity;
use metam::mapper::{CandidateBoundary, MapProfile, Mapping};

use crate::db::DatabaseMetadata;
use crate::sequencer::Sequencer;

#[derive(Clone, Debug)]
pub struct MetamMapping {
    pub read_id: String,
    pub ref_id: String,
    /// 1-based inclusive start coordinate for profiler coverage accounting.
    pub start_1based: usize,
    /// 1-based exclusive end coordinate for profiler coverage accounting.
    pub end_1based_exclusive: usize,
    pub identity: f64,
    pub read_len: usize,
    pub mapped_bases: usize,
    pub mate_ref_status_same_ref: bool,
    /// Mate number: 0 for an unpaired record, 1 for R1, or 2 for R2.
    pub mate_number: u8,
    /// Read-local identifier shared only by an explicitly selected R1/R2 pair.
    /// Zero marks an unpaired or diagnostic-only mapping.
    pub pair_id: u32,
    /// Taxonomic candidate levels that may be incomplete after configured
    /// mapper boundaries.
    pub candidate_boundary: CandidateBoundary,
}

/// Borrowed view of one coordinate-mapping row.
///
/// Saved mapping files contain the read and reference names inline. Keeping
/// those fields borrowed while the current line is processed avoids allocating
/// two short-lived `String`s for every mapping record. Consumers own a read
/// name only when it first becomes a grouping key; reference names remain
/// lookup-only.
#[derive(Clone, Copy, Debug)]
pub struct MetamMappingRecord<'a> {
    pub read_id: &'a str,
    pub ref_id: &'a str,
    pub start_1based: usize,
    pub end_1based_exclusive: usize,
    pub identity: f64,
    pub read_len: usize,
    pub mapped_bases: usize,
    pub mate_ref_status_same_ref: bool,
    pub mate_number: u8,
    pub pair_id: u32,
    pub candidate_boundary: CandidateBoundary,
}

pub type MetamMappingBatch = Vec<MetamMapping>;
pub type MetamMappingReceiver = Receiver<MetamMappingBatch>;

impl From<Mapping> for MetamMapping {
    fn from(mapping: Mapping) -> Self {
        Self {
            read_id: mapping.read_id,
            ref_id: mapping.ref_id,
            start_1based: (mapping.ref_start as usize).saturating_add(1),
            end_1based_exclusive: (mapping.ref_end as usize).saturating_add(1),
            identity: canonical_mapping_identity(mapping.identity),
            read_len: mapping.read_len,
            mapped_bases: mapping.mapped_bases,
            mate_ref_status_same_ref: mapping.mate_ref_status.is_same(),
            mate_number: mapping.mate_number,
            pair_id: mapping.pair_id,
            candidate_boundary: mapping.candidate_boundary,
        }
    }
}

impl MetamMapping {
    pub fn as_record(&self) -> MetamMappingRecord<'_> {
        MetamMappingRecord {
            read_id: &self.read_id,
            ref_id: &self.ref_id,
            start_1based: self.start_1based,
            end_1based_exclusive: self.end_1based_exclusive,
            identity: self.identity,
            read_len: self.read_len,
            mapped_bases: self.mapped_bases,
            mate_ref_status_same_ref: self.mate_ref_status_same_ref,
            mate_number: self.mate_number,
            pair_id: self.pair_id,
            candidate_boundary: self.candidate_boundary,
        }
    }

    pub fn validate(&self) -> Result<()> {
        self.as_record().validate()
    }
}

impl MetamMappingRecord<'_> {
    pub fn validate(self) -> Result<()> {
        nonempty(self.read_id, "read_id")?;
        nonempty(self.ref_id, "ref_id")?;
        if self.end_1based_exclusive < self.start_1based {
            bail!(
                "mapping end {} is smaller than start {}",
                self.end_1based_exclusive,
                self.start_1based
            );
        }
        if !self.identity.is_finite() || !(0.0..=1.0).contains(&self.identity) {
            bail!(
                "identity must be finite and within [0,1], got {}",
                self.identity
            );
        }
        if self.read_len == 0 {
            bail!("read_len must be greater than zero");
        }
        if self.mapped_bases == 0 {
            bail!("mapped_bases must be greater than zero");
        }
        if self.mapped_bases > self.read_len {
            bail!(
                "mapped_bases {} is greater than read_len {}",
                self.mapped_bases,
                self.read_len
            );
        }
        if self.mate_number > 2 {
            bail!("mate must be 0, 1, or 2, got {}", self.mate_number);
        }
        if self.mate_number == 0 && self.pair_id != 0 {
            bail!("unpaired mapping has nonzero pair_id {}", self.pair_id);
        }
        if self.pair_id != 0 && !self.mate_ref_status_same_ref {
            bail!(
                "mapping in pair {} does not have same-reference mate status",
                self.pair_id
            );
        }
        Ok(())
    }

    #[cfg(test)]
    fn into_owned(self) -> MetamMapping {
        MetamMapping {
            read_id: self.read_id.to_owned(),
            ref_id: self.ref_id.to_owned(),
            start_1based: self.start_1based,
            end_1based_exclusive: self.end_1based_exclusive,
            identity: self.identity,
            read_len: self.read_len,
            mapped_bases: self.mapped_bases,
            mate_ref_status_same_ref: self.mate_ref_status_same_ref,
            mate_number: self.mate_number,
            pair_id: self.pair_id,
            candidate_boundary: self.candidate_boundary,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct MetamProfile {
    pub total_reads: usize,
    pub mapped_reads: usize,
    pub emitted_records: usize,
}

pub type SharedMetamProfile = Arc<Mutex<Option<MetamProfile>>>;

impl MetamProfile {
    pub fn from_map_profile(profile: &MapProfile) -> Result<Self> {
        Ok(Self {
            total_reads: u64_to_usize(profile.total_reads, "total_reads")?,
            mapped_reads: u64_to_usize(profile.mapped_reads, "mapped_reads")?,
            emitted_records: u64_to_usize(profile.emitted_records, "emitted_records")?,
        })
    }
}

pub fn new_shared_metam_profile() -> SharedMetamProfile {
    Arc::new(Mutex::new(None))
}

pub fn store_shared_metam_profile(
    shared: &SharedMetamProfile,
    profile: MetamProfile,
) -> Result<()> {
    let mut guard = shared
        .lock()
        .map_err(|_| anyhow!("in-memory mapping profile lock is poisoned"))?;
    *guard = Some(profile);
    Ok(())
}

pub fn load_shared_metam_profile(shared: &SharedMetamProfile) -> Result<MetamProfile> {
    let guard = shared
        .lock()
        .map_err(|_| anyhow!("in-memory mapping profile lock is poisoned"))?;
    guard
        .as_ref()
        .copied()
        .ok_or_else(|| anyhow!("mapper finished without publishing mapping profile counters"))
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct MetamParseStats {
    pub records: usize,
    pub skipped_malformed: usize,
    pub profile: Option<MetamProfile>,
}

const MAPPING_METADATA_PREFIX: &str = "#metax_mapping";
const DATABASE_SUMMARY_PREFIX: &str = "#metax_mapping_database";
const SAMPLE_METADATA_PREFIX: &str = "#metax_mapping_sample";
const MAPPING_FILE_VERSION: &str = "5";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MappingDatabaseSummary {
    pub is_fractional: bool,
    pub sequences: usize,
    pub genomes: usize,
    /// xxh3-128 content fingerprint of the reference table; identifies the
    /// exact database a reused mapping was produced against.
    pub fingerprint: u128,
}

impl MappingDatabaseSummary {
    pub fn from_database(database: &DatabaseMetadata) -> Self {
        Self {
            is_fractional: database.is_fractional,
            sequences: database.references.len(),
            genomes: database
                .references
                .iter()
                .map(|reference| reference.genome.as_str())
                .collect::<HashSet<_>>()
                .len(),
            fingerprint: database.fingerprint(),
        }
    }
}

/// Profiler settings resolved when a mapping file is created. Persisting these
/// values makes later `--reuse-mapping` runs reproduce the original read-level
/// evidence gates even when a default was calibrated from FASTQ qualities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MappingProfilerMetadata {
    pub strain: bool,
    pub identity: f64,
    pub mapped_len: usize,
    pub fraction: f64,
    /// Mean predicted base accuracy sampled from the input FASTQ, when the
    /// adaptive long-read policy inspected qualities. Keeping the estimate
    /// allows a deliberately changed profile mode to apply its own margin.
    pub sampled_mean_base_accuracy: Option<f64>,
}

impl MappingProfilerMetadata {
    fn validate(self) -> Result<Self> {
        validate_open_unit_interval(self.identity, "profiler identity")?;
        if self.mapped_len == 0 {
            bail!("profiler mapped length must be greater than zero");
        }
        validate_open_unit_interval(self.fraction, "profiler mapped fraction")?;
        if let Some(accuracy) = self.sampled_mean_base_accuracy {
            validate_closed_unit_interval(accuracy, "sampled quality accuracy")?;
        }
        Ok(self)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MappingRunMetadata {
    pub sequencer: Sequencer,
    pub is_paired: bool,
    pub mapper_mode: String,
    pub profiler: MappingProfilerMetadata,
    pub database: MappingDatabaseSummary,
}

/// Write only run-level metadata. Reference metadata remains authoritative in
/// metax_db.tsv.gz and is never copied into a mapping file.
pub fn write_mapping_run_metadata(
    writer: &mut dyn Write,
    sequencer: Sequencer,
    is_paired: bool,
    mapper_mode: &str,
    profiler: MappingProfilerMetadata,
    database: &DatabaseMetadata,
) -> Result<()> {
    validate_mapper_mode(mapper_mode)?;
    let profiler = profiler.validate()?;
    let summary = MappingDatabaseSummary::from_database(database);
    let sampled_mean_base_accuracy = profiler
        .sampled_mean_base_accuracy
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_string());
    writeln!(
        writer,
        "{MAPPING_METADATA_PREFIX}\tversion={MAPPING_FILE_VERSION}\tsequencer={sequencer}\tpaired={is_paired}\tmapper_mode={mapper_mode}\tstrain={}\tprofiler_identity={}\tprofiler_mapped_len={}\tprofiler_fraction={}\tsampled_mean_base_accuracy={sampled_mean_base_accuracy}",
        profiler.strain,
        profiler.identity,
        profiler.mapped_len,
        profiler.fraction,
    )?;
    writeln!(
        writer,
        "{DATABASE_SUMMARY_PREFIX}\tfractional={}\tsequences={}\tgenomes={}\tfingerprint={:032x}",
        summary.is_fractional, summary.sequences, summary.genomes, summary.fingerprint
    )?;
    Ok(())
}

/// Read and validate the two required leading metadata lines.
pub fn load_mapping_run_metadata(path: &Path) -> Result<MappingRunMetadata> {
    if !path.exists() {
        bail!("mapping file {} does not exist", path.display());
    }
    let mut reader = open_metam_mapping_reader(path)?;
    let mut first = String::new();
    loop {
        first.clear();
        if reader.read_line(&mut first)? == 0 {
            bail!("mapping file {} is empty", path.display());
        }
        if !first.trim().is_empty() {
            break;
        }
    }
    let first = first.trim_end_matches(['\r', '\n']);
    if first.split('\t').next() != Some(MAPPING_METADATA_PREFIX) {
        bail!(
            "mapping file {} is missing the required '{}\tversion={}' header; regenerate it with --save-mapping",
            path.display(),
            MAPPING_METADATA_PREFIX,
            MAPPING_FILE_VERSION
        );
    }
    let mapping = parse_key_value_metadata(first, MAPPING_METADATA_PREFIX)?;
    let version = required_metadata(&mapping, "version")?;
    if version != MAPPING_FILE_VERSION {
        bail!(
            "unsupported mapping metadata version '{}' in {}; expected version {}",
            version,
            path.display(),
            MAPPING_FILE_VERSION
        );
    }
    require_exact_metadata_keys(
        &mapping,
        &[
            "version",
            "sequencer",
            "paired",
            "mapper_mode",
            "strain",
            "profiler_identity",
            "profiler_mapped_len",
            "profiler_fraction",
            "sampled_mean_base_accuracy",
        ],
        MAPPING_METADATA_PREFIX,
    )?;
    let sequencer = parse_sequencer(required_metadata(&mapping, "sequencer")?)?;
    let is_paired = parse_bool(required_metadata(&mapping, "paired")?, "paired")?;
    let mapper_mode = required_metadata(&mapping, "mapper_mode")?.to_string();
    validate_mapper_mode(&mapper_mode)?;
    let profiler = MappingProfilerMetadata {
        strain: parse_bool(required_metadata(&mapping, "strain")?, "strain")?,
        identity: parse_open_unit_interval(
            required_metadata(&mapping, "profiler_identity")?,
            "profiler identity",
        )?,
        mapped_len: parse_positive(
            required_metadata(&mapping, "profiler_mapped_len")?,
            "profiler mapped length",
        )?,
        fraction: parse_open_unit_interval(
            required_metadata(&mapping, "profiler_fraction")?,
            "profiler mapped fraction",
        )?,
        sampled_mean_base_accuracy: parse_optional_closed_unit_interval(
            required_metadata(&mapping, "sampled_mean_base_accuracy")?,
            "sampled quality accuracy",
        )?,
    }
    .validate()?;

    let mut second = String::new();
    loop {
        second.clear();
        if reader.read_line(&mut second)? == 0 {
            bail!(
                "mapping file {} is missing {}",
                path.display(),
                DATABASE_SUMMARY_PREFIX
            );
        }
        if !second.trim().is_empty() {
            break;
        }
    }
    let second = second.trim_end_matches(['\r', '\n']);
    let database = parse_key_value_metadata(second, DATABASE_SUMMARY_PREFIX)?;
    require_exact_metadata_keys(
        &database,
        &["fractional", "sequences", "genomes", "fingerprint"],
        DATABASE_SUMMARY_PREFIX,
    )?;
    let database = MappingDatabaseSummary {
        is_fractional: parse_bool(required_metadata(&database, "fractional")?, "fractional")?,
        sequences: parse_positive(required_metadata(&database, "sequences")?, "sequence count")?,
        genomes: parse_positive(required_metadata(&database, "genomes")?, "genome count")?,
        fingerprint: parse_fingerprint(required_metadata(&database, "fingerprint")?)?,
    };
    if database.genomes > database.sequences {
        bail!(
            "mapping metadata declares {} genomes but only {} sequences",
            database.genomes,
            database.sequences
        );
    }
    Ok(MappingRunMetadata {
        sequencer,
        is_paired,
        mapper_mode,
        profiler,
        database,
    })
}

pub fn validate_mapping_database_summary(
    summary: MappingDatabaseSummary,
    database: &DatabaseMetadata,
) -> Result<()> {
    let observed = MappingDatabaseSummary::from_database(database);
    if summary.fingerprint != observed.fingerprint {
        bail!(
            "these mappings were generated against a different database \
             (fingerprint {:032x}: fractional={} sequences={} genomes={}) \
             than the provided --db \
             (fingerprint {:032x}: fractional={} sequences={} genomes={}); \
             provide the matching --db or regenerate the mappings with --save-mapping",
            summary.fingerprint,
            summary.is_fractional,
            summary.sequences,
            summary.genomes,
            observed.fingerprint,
            observed.is_fractional,
            observed.sequences,
            observed.genomes
        );
    }
    Ok(())
}

fn parse_key_value_metadata(line: &str, expected_prefix: &str) -> Result<HashMap<String, String>> {
    let mut tokens = line.split('\t');
    let prefix = tokens.next().unwrap_or_default();
    if prefix != expected_prefix {
        bail!("expected metadata prefix {expected_prefix}, got '{prefix}'");
    }
    let mut fields = HashMap::new();
    for token in tokens {
        let Some((key, value)) = token.split_once('=') else {
            bail!("malformed metadata token '{token}' in {line}");
        };
        if key.is_empty() || value.is_empty() {
            bail!("empty metadata key or value in token '{token}'");
        }
        if fields.insert(key.to_string(), value.to_string()).is_some() {
            bail!("duplicate metadata key '{key}' in {line}");
        }
    }
    Ok(fields)
}

fn parse_fingerprint(value: &str) -> Result<u128> {
    u128::from_str_radix(value, 16)
        .with_context(|| format!("invalid database fingerprint '{value}'"))
}

fn require_exact_metadata_keys(
    fields: &HashMap<String, String>,
    expected: &[&str],
    prefix: &str,
) -> Result<()> {
    for key in expected {
        required_metadata(fields, key)?;
    }
    if fields.len() != expected.len() {
        let mut unexpected: Vec<&str> = fields
            .keys()
            .map(String::as_str)
            .filter(|key| !expected.contains(key))
            .collect();
        unexpected.sort_unstable();
        bail!(
            "{} contains unsupported metadata key(s): {}",
            prefix,
            unexpected.join(", ")
        );
    }
    Ok(())
}

fn required_metadata<'a>(fields: &'a HashMap<String, String>, key: &str) -> Result<&'a str> {
    fields
        .get(key)
        .map(String::as_str)
        .ok_or_else(|| anyhow!("mapping metadata is missing required key '{key}'"))
}

fn parse_sequencer(value: &str) -> Result<Sequencer> {
    match value.to_ascii_lowercase().as_str() {
        "illumina" => Ok(Sequencer::Illumina),
        "nanopore" | "ont" => Ok(Sequencer::Nanopore),
        "pacbio" => Ok(Sequencer::PacBio),
        "assembly" => Ok(Sequencer::Assembly),
        _ => bail!("unsupported sequencer '{value}' in mapping metadata"),
    }
}

fn parse_bool(value: &str, label: &str) -> Result<bool> {
    match value {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => bail!("{label} must be 'true' or 'false', got '{value}'"),
    }
}

fn parse_open_unit_interval(value: &str, label: &str) -> Result<f64> {
    let parsed = value
        .parse::<f64>()
        .with_context(|| format!("invalid {label} '{value}'"))?;
    validate_open_unit_interval(parsed, label)?;
    Ok(parsed)
}

fn parse_optional_closed_unit_interval(value: &str, label: &str) -> Result<Option<f64>> {
    if value == "NA" {
        return Ok(None);
    }
    let parsed = value
        .parse::<f64>()
        .with_context(|| format!("invalid {label} '{value}'"))?;
    validate_closed_unit_interval(parsed, label)?;
    Ok(Some(parsed))
}

fn validate_open_unit_interval(value: f64, label: &str) -> Result<()> {
    if value.is_finite() && value > 0.0 && value <= 1.0 {
        Ok(())
    } else {
        bail!("{label} must be finite and in (0, 1], got {value}")
    }
}

fn validate_closed_unit_interval(value: f64, label: &str) -> Result<()> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(())
    } else {
        bail!("{label} must be finite and in [0, 1], got {value}")
    }
}

fn validate_mapper_mode(value: &str) -> Result<()> {
    if matches!(value, "default" | "recall" | "precision") {
        Ok(())
    } else {
        bail!("unsupported mapper mode '{value}' in mapping metadata")
    }
}

pub fn open_metam_mapping_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open Metam output {}", path.display()))?;
    let reader: Box<dyn Read> = if path.extension().is_some_and(|extension| extension == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(reader)))
}

#[cfg(test)]
pub fn for_each_metam_mapping<R, F>(
    reader: R,
    source: &str,
    lenient: bool,
    mut on_mapping: F,
) -> Result<MetamParseStats>
where
    R: BufRead,
    F: FnMut(MetamMapping) -> Result<()>,
{
    for_each_metam_mapping_record(reader, source, lenient, |mapping| {
        on_mapping(mapping.into_owned())
    })
}

pub fn for_each_metam_mapping_record<R, F>(
    mut reader: R,
    source: &str,
    lenient: bool,
    mut on_mapping: F,
) -> Result<MetamParseStats>
where
    R: BufRead,
    F: for<'a> FnMut(MetamMappingRecord<'a>) -> Result<()>,
{
    let mut stats = MetamParseStats::default();
    let mut line = String::new();
    let mut line_no = 0usize;
    loop {
        line.clear();
        if reader
            .read_line(&mut line)
            .with_context(|| format!("failed reading {source}"))?
            == 0
        {
            break;
        }
        line_no += 1;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with(SAMPLE_METADATA_PREFIX) {
            if stats.profile.is_some() {
                bail!("duplicate mapping sample metadata in {source}");
            }
            stats.profile = Some(parse_sample_metadata_line(line)?);
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let mapping = match parse_mapping_record(line) {
            Ok(mapping) => mapping,
            Err(err) if lenient => {
                stats.skipped_malformed += 1;
                log::warn!(
                    target: "PROFILER",
                    "Skipping malformed mapping line {} in {}: {}",
                    line_no,
                    source,
                    err
                );
                continue;
            }
            Err(err) => {
                return Err(err).with_context(|| {
                    format!(
                        "failed to parse Metam mapping line {} in {}",
                        line_no, source
                    )
                });
            }
        };
        on_mapping(mapping)?;
        stats.records += 1;
    }
    Ok(stats)
}

#[cfg(test)]
pub fn for_each_structured_mapping<F>(
    receiver: MetamMappingReceiver,
    source: &str,
    mut on_mapping: F,
) -> Result<MetamParseStats>
where
    F: FnMut(MetamMapping) -> Result<()>,
{
    let mut stats = MetamParseStats::default();
    for batch in receiver {
        for mapping in batch {
            mapping.validate().with_context(|| {
                format!(
                    "invalid structured mapping record {} from {}",
                    stats.records + 1,
                    source
                )
            })?;
            on_mapping(mapping)?;
            stats.records += 1;
        }
    }
    Ok(stats)
}

pub fn for_each_structured_mapping_record<F>(
    receiver: MetamMappingReceiver,
    source: &str,
    mut on_mapping: F,
) -> Result<MetamParseStats>
where
    F: for<'a> FnMut(MetamMappingRecord<'a>) -> Result<()>,
{
    let mut stats = MetamParseStats::default();
    for batch in receiver {
        for mapping in batch {
            mapping.validate().with_context(|| {
                format!(
                    "invalid structured mapping record {} from {}",
                    stats.records + 1,
                    source
                )
            })?;
            on_mapping(mapping.as_record())?;
            stats.records += 1;
        }
    }
    Ok(stats)
}

pub fn format_sample_metadata_line(profile: &MetamProfile) -> String {
    format!(
        "{SAMPLE_METADATA_PREFIX}\ttotal_reads={}\tmapped_reads={}\temitted_records={}",
        profile.total_reads, profile.mapped_reads, profile.emitted_records
    )
}

fn parse_sample_metadata_line(line: &str) -> Result<MetamProfile> {
    let fields = parse_key_value_metadata(line, SAMPLE_METADATA_PREFIX)?;
    require_exact_metadata_keys(
        &fields,
        &["total_reads", "mapped_reads", "emitted_records"],
        SAMPLE_METADATA_PREFIX,
    )?;
    Ok(MetamProfile {
        total_reads: parse_nonnegative(required_metadata(&fields, "total_reads")?, "total_reads")?,
        mapped_reads: parse_nonnegative(
            required_metadata(&fields, "mapped_reads")?,
            "mapped_reads",
        )?,
        emitted_records: parse_nonnegative(
            required_metadata(&fields, "emitted_records")?,
            "emitted_records",
        )?,
    })
}

#[cfg(test)]
fn parse_mapping_line(line: &str) -> Result<MetamMapping> {
    Ok(parse_mapping_record(line)?.into_owned())
}

fn parse_mapping_record(line: &str) -> Result<MetamMappingRecord<'_>> {
    let mut fields = line.split('\t');
    let read_id = nonempty(next_field(&mut fields, "read_id")?, "read_id")?;
    let ref_id = nonempty(next_field(&mut fields, "ref_id")?, "ref_id")?;
    let ref_start: usize = parse_nonnegative(next_field(&mut fields, "ref_start")?, "ref_start")?;
    let ref_end: usize = parse_nonnegative(next_field(&mut fields, "ref_end")?, "ref_end")?;
    if ref_end < ref_start {
        bail!(
            "ref_end {} is smaller than ref_start {}",
            ref_end,
            ref_start
        );
    }
    let identity_field = next_field(&mut fields, "identity")?;
    let identity: f64 = identity_field
        .parse()
        .with_context(|| format!("invalid identity '{}'", identity_field))?;
    if !identity.is_finite() || !(0.0..=1.0).contains(&identity) {
        bail!("identity must be finite and within [0,1], got {}", identity);
    }
    let identity = canonical_mapping_identity(identity);
    let read_len = parse_positive(next_field(&mut fields, "read_len")?, "read_len")?;
    let mapped_bases = parse_positive(next_field(&mut fields, "mapped_bases")?, "mapped_bases")?;
    if mapped_bases > read_len {
        bail!(
            "mapped_bases {} is greater than read_len {}",
            mapped_bases,
            read_len
        );
    }
    let mate_ref_status = nonempty(
        next_field(&mut fields, "mate_ref_status")?,
        "mate_ref_status",
    )?;
    let mate_ref_status_same_ref = mate_ref_status == "=";
    let mate_number = u8::try_from(parse_nonnegative(next_field(&mut fields, "mate")?, "mate")?)
        .context("mate does not fit in u8")?;
    let pair_id = u32::try_from(parse_nonnegative(
        next_field(&mut fields, "pair_id")?,
        "pair_id",
    )?)
    .context("pair_id does not fit in u32")?;
    let candidate_boundary =
        parse_candidate_boundary(next_field(&mut fields, "candidate_boundary")?)?;
    if fields.next().is_some() {
        bail!("mapping row contains more than the required 11 tab-separated fields");
    }

    let mapping = MetamMappingRecord {
        read_id,
        ref_id,
        start_1based: ref_start.saturating_add(1),
        end_1based_exclusive: ref_end.saturating_add(1),
        identity,
        read_len,
        mapped_bases,
        mate_ref_status_same_ref,
        mate_number,
        pair_id,
        candidate_boundary,
    };
    mapping.validate()?;
    Ok(mapping)
}

fn next_field<'a>(fields: &mut impl Iterator<Item = &'a str>, label: &str) -> Result<&'a str> {
    fields
        .next()
        .ok_or_else(|| anyhow!("missing required field {label}"))
}

fn parse_candidate_boundary(value: &str) -> Result<CandidateBoundary> {
    let bits: u8 = value.parse().with_context(|| {
        format!("candidate_boundary must be an integer from 0 to 3, got '{value}'")
    })?;
    CandidateBoundary::from_bits(bits)
        .ok_or_else(|| anyhow!("candidate_boundary must be an integer from 0 to 3, got '{value}'"))
}

pub fn total_input_units(profile: &MetamProfile, is_paired: bool) -> usize {
    if is_paired {
        profile.total_reads.div_ceil(2)
    } else {
        profile.total_reads
    }
}

fn u64_to_usize(value: u64, label: &str) -> Result<usize> {
    usize::try_from(value).with_context(|| format!("{label} does not fit in usize"))
}

fn nonempty<'a>(value: &'a str, label: &str) -> Result<&'a str> {
    if value.is_empty() {
        Err(anyhow!("{label} is empty"))
    } else {
        Ok(value)
    }
}

fn parse_positive(value: &str, label: &str) -> Result<usize> {
    let parsed = parse_nonnegative(value, label)?;
    if parsed == 0 {
        bail!("{label} must be greater than zero");
    }
    Ok(parsed)
}

fn parse_nonnegative(value: &str, label: &str) -> Result<usize> {
    value
        .parse()
        .with_context(|| format!("invalid {label} '{}'", value))
}

#[cfg(test)]
mod tests {
    use super::*;
    use metam::format::format_mapping;
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

    use crate::db::ReferenceMetadata;

    fn test_database(fractional: bool) -> DatabaseMetadata {
        DatabaseMetadata {
            metam_index: PathBuf::from("metax_db.mxi"),
            references: vec![ReferenceMetadata {
                name: "G|ctg|s:8|100".to_string(),
                length: 100,
                genome: "G".to_string(),
                contig: "ctg".to_string(),
                full_genome_size: 100,
                effective_indexed_size: if fractional { 60 } else { 100 },
                taxonomy_field: "s:8".to_string(),
            }],
            is_fractional: fractional,
        }
    }

    fn test_profiler_metadata() -> MappingProfilerMetadata {
        MappingProfilerMetadata {
            strain: false,
            identity: 0.96,
            mapped_len: 50,
            fraction: 0.60,
            sampled_mean_base_accuracy: Some(0.92),
        }
    }

    fn test_mapping(identity: f64) -> Mapping {
        Mapping {
            read_id: "r1".to_string(),
            ref_id: "ref".to_string(),
            ref_start: 0,
            ref_end: 100,
            identity,
            read_len: 150,
            mapped_bases: 100,
            mate_ref_status: metam::mapper::MateRefStatus::NotApplicable,
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
            pair_status: metam::mapper::PairStatus::Single,
            mate_number: 0,
            pair_id: 0,
            candidate_groups: metam::mapper::ReferenceCandidateGroups::new(1, 1),
            candidate_boundary: CandidateBoundary::COMPLETE,
        }
    }

    #[test]
    fn parses_metam_mapping_line() {
        let m = parse_mapping_line("r1\tref\t0\t100\t0.9700\t150\t100\tNA\t0\t0\t0").unwrap();
        assert_eq!(m.read_id, "r1");
        assert_eq!(m.ref_id, "ref");
        assert_eq!(m.start_1based, 1);
        assert_eq!(m.end_1based_exclusive, 101);
        assert_eq!(m.mapped_bases, 100);
        assert!(!m.mate_ref_status_same_ref);
        assert_eq!(m.mate_number, 0);
        assert_eq!(m.pair_id, 0);
        assert!(m.candidate_boundary.is_complete());
        let paired = parse_mapping_line("r1/1\tref\t0\t100\t0.9700\t150\t100\t=\t1\t7\t3").unwrap();
        assert!(paired.mate_ref_status_same_ref);
        assert_eq!(paired.mate_number, 1);
        assert_eq!(paired.pair_id, 7);
        assert!(paired.candidate_boundary.profile_incomplete());
        assert!(paired.candidate_boundary.terminal_incomplete());
    }

    #[test]
    fn mapping_record_borrows_identifiers_from_the_input_line() {
        let line = "sample-read\tGCF_1|contig|s:7\t0\t100\t0.9700\t150\t100\tNA\t0\t0\t0";
        let record = parse_mapping_record(line).unwrap();
        assert_eq!(record.read_id, "sample-read");
        assert_eq!(record.ref_id, "GCF_1|contig|s:7");

        let input_start = line.as_ptr() as usize;
        let input_end = input_start + line.len();
        assert!((input_start..input_end).contains(&(record.read_id.as_ptr() as usize)));
        assert!((input_start..input_end).contains(&(record.ref_id.as_ptr() as usize)));

        let owned = record.into_owned();
        assert_eq!(owned.read_id, record.read_id);
        assert_eq!(owned.ref_id, record.ref_id);
        assert_eq!(owned.identity.to_bits(), record.identity.to_bits());
    }

    #[test]
    fn saved_and_structured_identity_are_bit_identical() {
        for identity in [0.8, 0.949_994, 0.949_996, 0.950_004, 0.993_333_333, 1.0] {
            let mapping = test_mapping(identity);
            let structured = MetamMapping::from(mapping.clone());
            let serialized = format_mapping(&mapping, false, false);
            let parsed = parse_mapping_line(&serialized).unwrap();
            assert_eq!(parsed.identity.to_bits(), structured.identity.to_bits());
        }
    }

    #[test]
    fn canonical_identity_has_stable_threshold_behavior() {
        let below = MetamMapping::from(test_mapping(0.949_994));
        let above = MetamMapping::from(test_mapping(0.949_996));
        assert!(below.identity < 0.95);
        assert_eq!(above.identity.to_bits(), 0.95_f64.to_bits());
    }

    #[test]
    fn near_tied_identities_canonicalize_identically() {
        let first = MetamMapping::from(test_mapping(0.950_001));
        let second = MetamMapping::from(test_mapping(0.950_004));
        assert_eq!(first.identity.to_bits(), second.identity.to_bits());
        assert_eq!(
            format_mapping(&test_mapping(0.950_001), false, false),
            format_mapping(&test_mapping(0.950_004), false, false)
        );
    }

    #[test]
    fn roundtrips_sample_metadata() {
        let profile = MetamProfile {
            total_reads: 123,
            mapped_reads: 45,
            emitted_records: 67,
        };
        let line = format_sample_metadata_line(&profile);
        assert!(line.starts_with("#metax_mapping_sample\t"));
        assert_eq!(parse_sample_metadata_line(&line).unwrap(), profile);
    }

    #[test]
    fn roundtrips_clean_run_metadata_without_reference_rows() {
        let database = test_database(false);
        let mut file = NamedTempFile::new().unwrap();
        write_mapping_run_metadata(
            file.as_file_mut(),
            Sequencer::Illumina,
            true,
            "recall",
            test_profiler_metadata(),
            &database,
        )
        .unwrap();
        writeln!(file, "r1\tG|ctg|s:8|100\t0\t100\t1\t100\t100\t=\t1\t1\t0").unwrap();
        file.flush().unwrap();

        let text = std::fs::read_to_string(file.path()).unwrap();
        assert_eq!(text.lines().count(), 3);
        let fingerprint = test_database(false).fingerprint();
        assert!(text.starts_with(&format!(
            "#metax_mapping\tversion=5\tsequencer=Illumina\tpaired=true\tmapper_mode=recall\tstrain=false\tprofiler_identity=0.96\tprofiler_mapped_len=50\tprofiler_fraction=0.6\tsampled_mean_base_accuracy=0.92\n#metax_mapping_database\tfractional=false\tsequences=1\tgenomes=1\tfingerprint={fingerprint:032x}\n"
        )));
        assert!(!text.contains("#metax_mapping_reference"));

        let metadata = load_mapping_run_metadata(file.path()).unwrap();
        assert_eq!(metadata.sequencer, Sequencer::Illumina);
        assert!(metadata.is_paired);
        assert_eq!(metadata.mapper_mode, "recall");
        assert_eq!(metadata.profiler, test_profiler_metadata());
        assert_eq!(
            metadata.database,
            MappingDatabaseSummary {
                is_fractional: false,
                sequences: 1,
                genomes: 1,
                fingerprint,
            }
        );
    }

    #[test]
    fn database_summary_rejects_wrong_database() {
        let summary = MappingDatabaseSummary::from_database(&test_database(false));
        let error = validate_mapping_database_summary(summary, &test_database(true)).unwrap_err();
        assert!(error.to_string().contains("different database"));
    }

    #[test]
    fn mapping_without_required_run_metadata_is_rejected() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "r1\tref\t0\t100\t1\t100\t100\t=\t1\t1\t0").unwrap();
        file.flush().unwrap();
        let error = load_mapping_run_metadata(file.path()).unwrap_err();
        assert!(error.to_string().contains("missing the required"));
        assert!(error.to_string().contains("--save-mapping"));
    }

    #[test]
    fn obsolete_mapping_version_is_reported_before_newer_required_keys() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "#metax_mapping\tversion=4\tsequencer=Illumina\tpaired=true\tmapper_mode=default"
        )
        .unwrap();
        file.flush().unwrap();

        let error = load_mapping_run_metadata(file.path()).unwrap_err();
        let message = error.to_string();
        assert!(message.contains("unsupported mapping metadata version '4'"));
        assert!(message.contains("expected version 5"));
    }

    #[test]
    fn obsolete_embedded_reference_metadata_is_rejected() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "#metax_mapping\tversion=5\tsequencer=Illumina\tpaired=false\tmapper_mode=default\tstrain=false\tprofiler_identity=0.97\tprofiler_mapped_len=50\tprofiler_fraction=0.7\tsampled_mean_base_accuracy=NA"
        )
        .unwrap();
        writeln!(
            file,
            "#metax_mapping_database\tversion=3\tfractional=false\tsequences=1\tgenomes=1\tfingerprint=00000000000000000000000000000000"
        )
        .unwrap();
        file.flush().unwrap();

        let error = load_mapping_run_metadata(file.path()).unwrap_err();
        assert!(error.to_string().contains("unsupported metadata key"));
        assert!(error.to_string().contains("version"));
    }

    #[test]
    fn invalid_profiler_metadata_is_rejected_before_writing() {
        let mut output = Vec::new();
        let error = write_mapping_run_metadata(
            &mut output,
            Sequencer::Nanopore,
            false,
            "default",
            MappingProfilerMetadata {
                fraction: 0.0,
                ..test_profiler_metadata()
            },
            &test_database(false),
        )
        .unwrap_err();
        assert!(error.to_string().contains("profiler mapped fraction"));
        assert!(output.is_empty());
    }

    #[test]
    fn malformed_sampled_accuracy_is_rejected_on_load() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "#metax_mapping\tversion=5\tsequencer=Nanopore\tpaired=false\tmapper_mode=default\tstrain=false\tprofiler_identity=0.9\tprofiler_mapped_len=500\tprofiler_fraction=0.2\tsampled_mean_base_accuracy=1.1"
        )
        .unwrap();
        file.flush().unwrap();

        let error = load_mapping_run_metadata(file.path()).unwrap_err();
        assert!(error.to_string().contains("sampled quality accuracy"));
    }

    #[test]
    fn parser_skips_comments_and_returns_embedded_profile() {
        let input = concat!(
            "#arbitrary comment\n",
            "r1\tref\t0\t100\t0.9700\t150\t100\tNA\t0\t0\t0\n",
            "#metax_mapping_sample\ttotal_reads=10\tmapped_reads=1\temitted_records=1\n",
        );
        let mut seen = 0usize;
        let stats = for_each_metam_mapping(input.as_bytes(), "test", false, |_mapping| {
            seen += 1;
            Ok(())
        })
        .unwrap();
        assert_eq!(seen, 1);
        assert_eq!(stats.records, 1);
        assert_eq!(
            stats.profile,
            Some(MetamProfile {
                total_reads: 10,
                mapped_reads: 1,
                emitted_records: 1,
            })
        );
    }

    #[test]
    fn rejects_obsolete_optional_counters_and_extra_mapping_fields() {
        let error = parse_sample_metadata_line(
            "#metax_mapping_sample\ttotal_reads=NA\tmapped_reads=1\temitted_records=1",
        )
        .unwrap_err();
        assert!(error.to_string().contains("total_reads"));

        let error =
            parse_mapping_line("r1\tref\t0\t100\t1\t100\t100\tNA\t0\t0\t0\textra").unwrap_err();
        assert!(error.to_string().contains("required 11"));
    }

    #[test]
    fn rejects_mapped_bases_larger_than_read_len() {
        let err = parse_mapping_line("r1\tref\t0\t100\t0.9700\t90\t100\tNA\t0\t0\t0").unwrap_err();
        assert!(err.to_string().contains("mapped_bases"));
    }

    #[test]
    fn rejects_inconsistent_explicit_pair_metadata() {
        let unpaired_with_pair =
            parse_mapping_line("r1\tref\t0\t100\t1\t100\t100\t=\t0\t1\t0").unwrap_err();
        assert!(unpaired_with_pair.to_string().contains("nonzero pair_id"));

        let discordant_pair =
            parse_mapping_line("r1\tref\t0\t100\t1\t100\t100\tother\t1\t1\t0").unwrap_err();
        assert!(discordant_pair
            .to_string()
            .contains("does not have same-reference mate status"));

        let invalid_mate =
            parse_mapping_line("r1\tref\t0\t100\t1\t100\t100\t=\t3\t1\t0").unwrap_err();
        assert!(invalid_mate.to_string().contains("mate must be 0, 1, or 2"));

        let invalid_candidate_boundary =
            parse_mapping_line("r1\tref\t0\t100\t1\t100\t100\tNA\t0\t0\t4").unwrap_err();
        assert!(invalid_candidate_boundary
            .to_string()
            .contains("candidate_boundary must be an integer from 0 to 3"));
    }

    #[test]
    fn lenient_parser_skips_bad_mapping_lines() {
        let input = b"good\tref\t0\t10\t1.0\t10\t10\tNA\t0\t0\t0\n\
                      bad\tref\t0\t10\t1.0\t5\t10\tNA\t0\t0\t0\n\
                      good2\tref\t0\t10\t1.0\t10\t10\t=\t1\t1\t0\n";
        let mut seen = Vec::new();
        let stats = for_each_metam_mapping(&input[..], "test", true, |mapping| {
            seen.push(mapping.read_id);
            Ok(())
        })
        .unwrap();
        assert_eq!(stats.records, 2);
        assert_eq!(stats.skipped_malformed, 1);
        assert_eq!(seen, vec!["good", "good2"]);
    }

    #[test]
    fn structured_stream_rejects_invalid_mapper_records() {
        let (sender, receiver) = std::sync::mpsc::channel();
        sender
            .send(vec![MetamMapping {
                read_id: "r1".to_string(),
                ref_id: "ref".to_string(),
                start_1based: 1,
                end_1based_exclusive: 101,
                identity: 0.99,
                read_len: 90,
                mapped_bases: 100,
                mate_ref_status_same_ref: false,
                mate_number: 0,
                pair_id: 0,
                candidate_boundary: CandidateBoundary::COMPLETE,
            }])
            .unwrap();
        drop(sender);

        let err = for_each_structured_mapping(receiver, "test stream", |_| Ok(())).unwrap_err();
        let message = format!("{err:#}");
        assert!(message.contains("invalid structured mapping record 1"));
        assert!(message.contains("mapped_bases"));
    }
}

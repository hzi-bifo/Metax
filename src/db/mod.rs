use anyhow::{anyhow, bail, Context, Result};
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use indexmap::IndexMap;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::progress::{self, ProgressOptions};
use metam::{
    index::{
        build_index_file_with_progress, BuildProgress, ReferenceOrderDigest,
        DEFAULT_INDEX_SEED_CHUNK_MIB,
    },
    progress::log_message,
    seeds::{SeedMode, StrobeStartDistance},
};

pub const DB_TSV_FILE: &str = "metax_db.tsv.gz";
const METAM_INDEX_FILE: &str = "metax_db.mxi";
const MAPPING_REFERENCE_FILE: &str = "metax_reference.fasta";
/// Metadata layout for databases whose references are all genomes. It is
/// still written for such databases so that earlier Metax 2 releases can read
/// them.
const DB_VERSION_GENOMES: &str = "2";
/// Metadata layout with a per-reference `reference_type` column, written when
/// at least one genome is not of type `genome`.
const DB_VERSION_TYPED: &str = "3";
/// Typed layout plus a `declared_genome_size` column, written when at least
/// one marker genome declares the size of its whole genome.
const DB_VERSION_SIZED: &str = "4";
const COLUMNS_GENOMES: &str =
    "reference_name\treference_length\tfull_genome_size\teffective_indexed_size\tgenome\tcontig\ttaxonomy";
pub const DEFAULT_SEED_CHUNK_MIB: usize = DEFAULT_INDEX_SEED_CHUNK_MIB;

/// What the reference sequences of one genome represent.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub enum ReferenceType {
    /// Assembled genome sequences: chromosomes, plasmids, scaffolds or contigs.
    #[default]
    Genome,
    /// A set of marker genes, one reference sequence per gene. The genome size
    /// is the summed marker length, so read counts are reported as mapped.
    Markers,
}

impl ReferenceType {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Genome => "genome",
            Self::Markers => "markers",
        }
    }
}

impl std::str::FromStr for ReferenceType {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value.to_ascii_lowercase().as_str() {
            "genome" => Ok(Self::Genome),
            "markers" => Ok(Self::Markers),
            _ => bail!("unsupported reference type '{value}'; expected genome or markers"),
        }
    }
}

impl std::fmt::Display for ReferenceType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Reference type of every genome in an index build: a default and the
/// genomes listed in a `genome_id<TAB>type` table.
#[derive(Clone, Debug, Default)]
pub struct ReferenceTypes {
    default: ReferenceType,
    by_genome: HashMap<String, ReferenceType>,
}

impl ReferenceTypes {
    pub fn new(default: ReferenceType, table: Option<&Path>) -> Result<Self> {
        let by_genome = match table {
            Some(path) => read_reference_type_table(path)?,
            None => HashMap::new(),
        };
        Ok(Self { default, by_genome })
    }

    /// Whether a table sets the type of single genomes.
    pub fn has_table(&self) -> bool {
        !self.by_genome.is_empty()
    }

    fn of(&self, genome: &str) -> ReferenceType {
        self.by_genome.get(genome).copied().unwrap_or(self.default)
    }
}

fn read_reference_type_table(path: &Path) -> Result<HashMap<String, ReferenceType>> {
    let reader = BufReader::new(open_text_reader(path)?);
    let mut types = HashMap::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("failed reading {}", path.display()))?;
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (genome, kind) = line
            .split_once('\t')
            .filter(|(_, kind)| !kind.contains('\t'))
            .ok_or_else(|| {
                anyhow!(
                    "line {} of {} must be genome_id<TAB>reference_type",
                    line_no + 1,
                    path.display()
                )
            })?;
        let kind: ReferenceType = kind
            .trim()
            .parse()
            .with_context(|| format!("line {} of {}", line_no + 1, path.display()))?;
        if types
            .insert(genome.trim().to_string(), kind)
            .is_some_and(|old| old != kind)
        {
            bail!(
                "genome {} has two reference types in {}",
                genome.trim(),
                path.display()
            );
        }
    }
    Ok(types)
}

#[derive(Clone, Debug)]
pub struct DatabaseMetadata {
    pub metam_index: PathBuf,
    pub references: Vec<ReferenceMetadata>,
    pub is_fractional: bool,
}

impl DatabaseMetadata {
    /// `Some(type)` when every reference has that type, `None` when mixed.
    pub fn uniform_reference_type(&self) -> Option<ReferenceType> {
        let first = self.references.first()?.reference_type;
        self.references
            .iter()
            .all(|reference| reference.reference_type == first)
            .then_some(first)
    }

    /// `genome`, `markers`, or `mixed`.
    pub fn reference_type_label(&self) -> &'static str {
        self.uniform_reference_type()
            .map_or("mixed", ReferenceType::as_str)
    }

    /// Stable content fingerprint of the reference table. Reused mappings
    /// (`--reuse-mapping` / `--resume`) record this so profiling can verify
    /// they were produced against this exact database. It hashes the parsed
    /// reference identities in load order — not the compressed file bytes — so
    /// recompressing `metax_db.tsv.gz` never changes it, while any real
    /// reference, length, taxonomy, or reference-type difference does. Genome
    /// rows hash as before reference types existed, so databases of genomes
    /// keep their fingerprint.
    pub fn fingerprint(&self) -> u128 {
        use xxhash_rust::xxh3::Xxh3;

        let mut hasher = Xxh3::new();
        hasher.update(if self.is_fractional {
            b"fractional\n"
        } else {
            b"absolute\n"
        });
        let mut row = Vec::new();
        for reference in &self.references {
            row.clear();
            // Tab-delimited canonical identity, mirroring metax_db.tsv.gz columns.
            row.extend_from_slice(reference.name.as_bytes());
            for value in [
                reference.length,
                reference.full_genome_size,
                reference.effective_indexed_size,
            ] {
                row.push(b'\t');
                push_decimal(&mut row, value);
            }
            for text in [
                &reference.genome,
                &reference.contig,
                &reference.taxonomy_field,
            ] {
                row.push(b'\t');
                row.extend_from_slice(text.as_bytes());
            }
            if reference.reference_type != ReferenceType::Genome {
                row.push(b'\t');
                row.extend_from_slice(reference.reference_type.as_str().as_bytes());
            }
            if let Some(size) = reference.declared_genome_size {
                row.push(b'\t');
                push_decimal(&mut row, size);
            }
            row.push(b'\n');
            hasher.update(&row);
        }
        hasher.digest128()
    }
}

/// Append `value` in decimal, as `Display` writes it.
fn push_decimal(out: &mut Vec<u8>, mut value: usize) {
    let mut digits = [0u8; 20];
    let mut start = digits.len();
    loop {
        start -= 1;
        digits[start] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    out.extend_from_slice(&digits[start..]);
}

#[derive(Clone, Debug)]
pub struct ReferenceMetadata {
    pub name: String,
    pub length: usize,
    pub genome: String,
    pub contig: String,
    pub full_genome_size: usize,
    pub effective_indexed_size: usize,
    pub taxonomy_field: String,
    pub reference_type: ReferenceType,
    /// Size of the whole genome as declared in the reference FASTA, kept for
    /// marker genomes only, whose `full_genome_size` is their summed marker
    /// length. Used for the unclassified estimate, never for scaling.
    pub declared_genome_size: Option<usize>,
}

impl ReferenceMetadata {
    /// Taxid of the profile endpoint the reference supports.
    pub fn profile_taxid(&self) -> Result<u32> {
        Ok(parse_taxonomy_field(&self.taxonomy_field)?.profile_taxid)
    }
}

/// Genome-type species of every genus that also holds marker-type species,
/// keyed by genus. Reads from the part of a marker-only genome that is not
/// indexed can land on these genomes of its relatives.
pub fn genome_species_next_to_markers(
    species: impl IntoIterator<Item = (ReferenceType, u32)>,
    genus_of: impl Fn(u32) -> Option<u32>,
) -> HashMap<u32, HashSet<u32>> {
    let mut genera: HashMap<u32, (bool, HashSet<u32>)> = HashMap::new();
    for (kind, species) in species {
        let Some(genus) = genus_of(species) else {
            continue;
        };
        let (has_markers, genomes) = genera.entry(genus).or_default();
        match kind {
            ReferenceType::Markers => *has_markers = true,
            ReferenceType::Genome => {
                genomes.insert(species);
            }
        }
    }
    genera
        .into_iter()
        .filter(|(_, (has_markers, genomes))| *has_markers && !genomes.is_empty())
        .map(|(genus, (_, genomes))| (genus, genomes))
        .collect()
}

#[derive(Clone, Debug)]
pub struct MetamIndexOptions {
    pub k: usize,
    pub w: usize,
    pub max_occ: usize,
    pub hash_bits: u8,
    pub seed_mode: String,
    pub match_both_strobes: bool,
    pub strobe_start_dist_min: u16,
    pub strobe_start_dist_max: u16,
    pub seed_chunk_mib: usize,
    pub temp_dir: Option<PathBuf>,
    pub quiet: bool,
    pub progress: bool,
    pub no_progress: bool,
    pub progress_interval: u64,
}

#[derive(Clone, Copy)]
pub struct SubsampleParams {
    pub fraction: f64,
    pub seed: u64,
    pub segment_length: usize,
    pub min_length: usize,
}

impl DatabaseMetadata {
    pub fn load(path: &Path) -> Result<Self> {
        let resolved = resolve_db_path(path)?;
        load_db_tsv(&resolved)
    }

    pub fn load_resolved(path: &Path) -> Result<(PathBuf, Self)> {
        let resolved = resolve_db_path(path)?;
        let metadata = load_db_tsv(&resolved)?;
        Ok((resolved, metadata))
    }
}

pub fn build_database(
    fasta: &str,
    outdir: &str,
    subsample: Option<SubsampleParams>,
    threads: Option<usize>,
    _compress_reference: bool,
    types: &ReferenceTypes,
    options: &MetamIndexOptions,
) -> Result<PathBuf> {
    let outdir_path = Path::new(outdir);
    match fs::symlink_metadata(outdir_path) {
        Ok(_) => bail!(
            "output database directory {} already exists; choose a new path or remove it before indexing",
            outdir_path.display()
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(error).with_context(|| {
                format!(
                    "failed to inspect output database path {}",
                    outdir_path.display()
                )
            });
        }
    }

    let fasta_path = Path::new(fasta);
    if !fasta_path.exists() {
        bail!("reference FASTA {} does not exist", fasta);
    }
    if let Some(params) = subsample {
        validate_subsample_params(params, threads)?;
    }

    let input = read_reference_fasta(fasta_path, types)?;
    let selected = select_index_records(&input, subsample, threads)?;
    drop(input);
    if selected.is_empty() {
        bail!("reference FASTA produced no indexable records");
    }

    let staging_path = create_database_staging_directory(outdir_path)?;
    let mut cleanup_staging = RemoveDirOnDrop::new(staging_path.clone());
    let force_sampled_field =
        requires_effective_size_reference_names(&selected, subsample.is_some());
    let index_path = staging_path.join(METAM_INDEX_FILE);

    // Whenever effective indexed size differs from full genome size, rewrite
    // reference names so every mapping row carries both values. This includes
    // explicit fractional subsampling and references whose declared full size
    // is larger than the supplied indexed sequence.
    let (references, indexed) = if force_sampled_field {
        let reference_path = staging_path.join(MAPPING_REFERENCE_FILE);
        let _cleanup_reference = RemoveFileOnDrop::new(reference_path.clone());
        let references =
            write_metam_reference(&reference_path, fasta_path, selected, force_sampled_field)?;
        let indexed = run_metam_index(
            &reference_path,
            &index_path,
            threads.unwrap_or(1).max(1),
            options,
        )?;
        gzip_file(&reference_path)
            .with_context(|| format!("failed to compress {}", reference_path.display()))?;
        (references, indexed)
    } else {
        let references = build_direct_reference_metadata(selected)?;
        let indexed = run_metam_index(
            fasta_path,
            &index_path,
            threads.unwrap_or(1).max(1),
            options,
        )?;
        (references, indexed)
    };
    verify_index_reference_order(&references, indexed)?;

    let metadata = DatabaseMetadata {
        metam_index: PathBuf::from(METAM_INDEX_FILE),
        is_fractional: force_sampled_field,
        references,
    };
    let tsv_path = staging_path.join(DB_TSV_FILE);
    write_db_tsv(&tsv_path, &metadata)?;
    validate_staged_database(&tsv_path, &metadata)?;
    fs::rename(&staging_path, outdir_path).with_context(|| {
        format!(
            "failed to publish completed database at {}",
            outdir_path.display()
        )
    })?;
    cleanup_staging.disarm();
    Ok(outdir_path.join(DB_TSV_FILE))
}

fn create_database_staging_directory(target: &Path) -> Result<PathBuf> {
    let parent = target
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = target.file_name().ok_or_else(|| {
        anyhow!(
            "output database path {} must name a directory",
            target.display()
        )
    })?;
    fs::create_dir_all(parent)
        .with_context(|| format!("failed to create output parent {}", parent.display()))?;
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or(0);
    for attempt in 0..100u32 {
        let mut staging_name = std::ffi::OsString::from(".");
        staging_name.push(name);
        staging_name.push(format!(
            ".metax-staging-{}-{nanos}-{attempt}",
            process::id()
        ));
        let path = parent.join(staging_name);
        match fs::create_dir(&path) {
            Ok(()) => return Ok(path),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(error).with_context(|| {
                    format!("failed to create database staging area {}", path.display())
                });
            }
        }
    }
    bail!(
        "failed to create a unique database staging area beside {}",
        target.display()
    )
}

fn validate_staged_database(path: &Path, expected: &DatabaseMetadata) -> Result<()> {
    let loaded = load_db_tsv(path).context("failed to validate staged database metadata")?;
    if loaded.is_fractional != expected.is_fractional
        || loaded.references.len() != expected.references.len()
        || loaded.fingerprint() != expected.fingerprint()
    {
        bail!("staged database metadata did not round-trip exactly");
    }
    let index = fs::metadata(&loaded.metam_index).with_context(|| {
        format!(
            "staged database mapping index {} is missing",
            loaded.metam_index.display()
        )
    })?;
    if !index.is_file() || index.len() == 0 {
        bail!(
            "staged database mapping index {} is not a non-empty file",
            loaded.metam_index.display()
        );
    }
    Ok(())
}

fn requires_effective_size_reference_names(
    records: &[SelectedRecord],
    subsampling_requested: bool,
) -> bool {
    subsampling_requested
        || records
            .iter()
            .any(|record| record.effective_indexed_size < record.full_genome_size)
}

struct RemoveFileOnDrop {
    path: PathBuf,
    armed: bool,
}

struct RemoveDirOnDrop {
    path: PathBuf,
    armed: bool,
}

impl RemoveDirOnDrop {
    fn new(path: PathBuf) -> Self {
        Self { path, armed: true }
    }

    fn disarm(&mut self) {
        self.armed = false;
    }
}

impl Drop for RemoveDirOnDrop {
    fn drop(&mut self) {
        if self.armed {
            let _ = fs::remove_dir_all(&self.path);
        }
    }
}

impl RemoveFileOnDrop {
    fn new(path: PathBuf) -> Self {
        Self { path, armed: true }
    }

    fn disarm(&mut self) {
        self.armed = false;
    }
}

impl Drop for RemoveFileOnDrop {
    fn drop(&mut self) {
        if self.armed {
            let _ = fs::remove_file(&self.path);
        }
    }
}

fn resolve_db_path(path: &Path) -> Result<PathBuf> {
    if path.is_dir() {
        let gz = path.join(DB_TSV_FILE);
        if gz.exists() {
            return Ok(gz);
        }
        bail!("{} does not contain {}", path.display(), DB_TSV_FILE);
    }
    if path.file_name().and_then(|name| name.to_str()) != Some(DB_TSV_FILE) {
        bail!(
            "database metadata must be a {} file or a directory containing it: {}",
            DB_TSV_FILE,
            path.display()
        );
    }
    Ok(path.to_path_buf())
}

fn load_db_tsv(path: &Path) -> Result<DatabaseMetadata> {
    let reader = open_metadata_reader(path)?;
    let base = path.parent().unwrap_or_else(|| Path::new("."));
    // Columns after the genome columns: none (version 2), reference_type
    // (version 3), reference_type and declared_genome_size (version 4).
    let mut extra_columns: Option<usize> = None;
    let mut declared_type: Option<String> = None;
    let mut metam_index: Option<PathBuf> = None;
    let mut is_fractional = false;
    let mut references = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("failed reading {}", path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.as_slice() {
            ["#metax_db_version", version] => {
                extra_columns = Some(match *version {
                    DB_VERSION_GENOMES => 0,
                    DB_VERSION_TYPED => 1,
                    DB_VERSION_SIZED => 2,
                    _ => bail!(
                        "unsupported Metax DB version {} in {}",
                        version,
                        path.display()
                    ),
                });
            }
            ["#mapping_index", value] => {
                metam_index = Some(resolve_relative_path(base, value));
            }
            ["#fractional", value] => {
                is_fractional = parse_bool(value)
                    .with_context(|| format!("invalid fractional flag on line {}", line_no + 1))?;
            }
            ["#reference_type", value] if extra_columns.is_some_and(|extra| extra > 0) => {
                declared_type = Some((*value).to_string());
            }
            ["#columns", columns @ ..] => {
                if columns.join("\t") != column_names(extra_columns.unwrap_or(0)) {
                    bail!(
                        "unsupported database metadata columns on line {} in {}: {}",
                        line_no + 1,
                        path.display(),
                        line
                    );
                }
            }
            [name, length, full, effective, genome, contig, taxonomy, rest @ ..]
                if Some(rest.len()) == extra_columns =>
            {
                let reference_type = match rest.first() {
                    Some(value) => value.parse().with_context(|| {
                        format!("invalid reference type on line {}", line_no + 1)
                    })?,
                    None => ReferenceType::Genome,
                };
                let declared_genome_size = match rest.get(1) {
                    None | Some(&"NA") => None,
                    Some(_) if reference_type != ReferenceType::Markers => bail!(
                        "a declared genome size is kept only for marker genomes (line {})",
                        line_no + 1
                    ),
                    Some(value) => Some(
                        parse_positive_usize(value, "declared genome size").with_context(|| {
                            format!("invalid declared genome size on line {}", line_no + 1)
                        })?,
                    ),
                };
                references.push(ReferenceMetadata {
                    name: (*name).to_string(),
                    length: parse_positive_usize(length, "reference length").with_context(
                        || format!("invalid reference length on line {}", line_no + 1),
                    )?,
                    genome: (*genome).to_string(),
                    contig: (*contig).to_string(),
                    full_genome_size: parse_positive_usize(full, "full genome size").with_context(
                        || format!("invalid full genome size on line {}", line_no + 1),
                    )?,
                    effective_indexed_size: parse_positive_usize(
                        effective,
                        "effective indexed size",
                    )
                    .with_context(|| {
                        format!("invalid effective indexed size on line {}", line_no + 1)
                    })?,
                    taxonomy_field: (*taxonomy).to_string(),
                    reference_type,
                    declared_genome_size,
                });
            }
            _ => bail!(
                "malformed database metadata line {} in {}: {}",
                line_no + 1,
                path.display(),
                line
            ),
        }
    }

    if extra_columns.is_none() {
        bail!("{} is not a Metax database metadata file", path.display());
    }
    let metam_index =
        metam_index.ok_or_else(|| anyhow!("{} is missing mapping_index", path.display()))?;
    if references.is_empty() {
        bail!("{} contains no reference metadata", path.display());
    }
    let metadata = DatabaseMetadata {
        metam_index,
        references,
        is_fractional,
    };
    if let Some(declared) = declared_type {
        if declared != metadata.reference_type_label() {
            bail!(
                "{} declares reference_type {} but its rows are {}",
                path.display(),
                declared,
                metadata.reference_type_label()
            );
        }
    }
    Ok(metadata)
}

fn open_metadata_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file =
        File::open(path).with_context(|| format!("failed to open database {}", path.display()))?;
    if path.extension().is_some_and(|extension| extension == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn resolve_relative_path(base: &Path, value: &str) -> PathBuf {
    let path = Path::new(value);
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base.join(path)
    }
}

fn write_db_tsv(path: &Path, metadata: &DatabaseMetadata) -> Result<()> {
    if path.extension().is_some_and(|extension| extension == "gz") {
        let file =
            File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut encoder = GzEncoder::new(file, Compression::fast());
        write_db_tsv_body(&mut encoder, metadata)?;
        encoder
            .finish()
            .with_context(|| format!("failed to finalize {}", path.display()))?;
        return Ok(());
    }
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    write_db_tsv_body(&mut writer, metadata)?;
    writer.flush()?;
    Ok(())
}

/// The `#columns` line of a metadata layout with `extra` columns after the
/// genome columns.
fn column_names(extra: usize) -> String {
    let mut columns = COLUMNS_GENOMES.to_string();
    for name in ["reference_type", "declared_genome_size"]
        .iter()
        .take(extra)
    {
        columns.push('\t');
        columns.push_str(name);
    }
    columns
}

fn write_db_tsv_body(writer: &mut dyn Write, metadata: &DatabaseMetadata) -> Result<()> {
    let sized = metadata
        .references
        .iter()
        .any(|reference| reference.declared_genome_size.is_some());
    let typed = sized || metadata.uniform_reference_type() != Some(ReferenceType::Genome);
    let (version, extra_columns) = if sized {
        (DB_VERSION_SIZED, 2)
    } else if typed {
        (DB_VERSION_TYPED, 1)
    } else {
        (DB_VERSION_GENOMES, 0)
    };
    writeln!(writer, "#metax_db_version\t{version}")?;
    writeln!(
        writer,
        "#mapping_index\t{}",
        metadata.metam_index.to_string_lossy()
    )?;
    writeln!(writer, "#fractional\t{}", metadata.is_fractional)?;
    if typed {
        writeln!(
            writer,
            "#reference_type\t{}",
            metadata.reference_type_label()
        )?;
    }
    writeln!(writer, "#columns\t{}", column_names(extra_columns))?;
    for r in &metadata.references {
        ensure_tsv_safe(&r.name, "reference name")?;
        ensure_tsv_safe(&r.genome, "genome accession")?;
        ensure_tsv_safe(&r.contig, "contig id")?;
        ensure_tsv_safe(&r.taxonomy_field, "taxonomy field")?;
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.name,
            r.length,
            r.full_genome_size,
            r.effective_indexed_size,
            r.genome,
            r.contig,
            r.taxonomy_field
        )?;
        if typed {
            write!(writer, "\t{}", r.reference_type)?;
        }
        if sized {
            match r.declared_genome_size {
                Some(size) => write!(writer, "\t{size}")?,
                None => write!(writer, "\tNA")?,
            }
        }
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct TaxonomySpec {
    raw: String,
    primary_taxid: Option<u32>,
    profile_taxid: u32,
}

#[derive(Clone, Debug)]
struct InputContig {
    contig: String,
    reference_name: String,
    length: usize,
}

#[derive(Clone, Debug)]
struct GenomeInput {
    genome: String,
    taxonomy: TaxonomySpec,
    provided_genome_size: Option<usize>,
    contigs: Vec<InputContig>,
    full_genome_size: usize,
    reference_type: ReferenceType,
    /// The provided size of a marker genome, kept beside its marker length.
    declared_genome_size: Option<usize>,
}

/// Position of one reference FASTA record: its genome in
/// `ReferenceInput::genomes` and its contig within that genome.
#[derive(Clone, Copy, Debug)]
struct RecordPosition {
    genome: u32,
    contig: u32,
}

/// Reference FASTA contents: genomes in order of first appearance, and the
/// position of every record in input order. The records of a genome need not
/// be adjacent. The index stores references in input order, so the database
/// metadata follows `record_order` rather than the genome order.
struct ReferenceInput {
    genomes: Vec<GenomeInput>,
    record_order: Vec<RecordPosition>,
}

/// Reference count and `ReferenceOrderDigest` reported by the index build.
#[derive(Clone, Copy, Debug)]
struct IndexedReferences {
    count: usize,
    order_digest: u64,
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

#[derive(Clone)]
struct SelectedRecord {
    genome: String,
    contig: String,
    source_contig: String,
    source_name: String,
    taxonomy: TaxonomySpec,
    reference_type: ReferenceType,
    full_genome_size: usize,
    effective_indexed_size: usize,
    declared_genome_size: Option<usize>,
    source_start: usize,
    source_end: usize,
}

impl SelectedRecord {
    fn len(&self) -> usize {
        self.source_end
            .saturating_sub(self.source_start)
            .saturating_add(1)
    }
}

fn read_reference_fasta(path: &Path, types: &ReferenceTypes) -> Result<ReferenceInput> {
    let reader = open_text_reader(path)?;
    let reader = BufReader::new(reader);
    let mut genomes: IndexMap<String, GenomeInput> = IndexMap::new();
    let mut record_order: Vec<RecordPosition> = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_len = 0usize;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("failed reading FASTA {}", path.display()))?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(header) = current_header.take() {
                record_order.push(store_input_record(&mut genomes, &header, current_len)?);
                current_len = 0;
            }
            current_header = Some(rest.trim().to_string());
        } else {
            let sequence_line = line.trim();
            if sequence_line.is_empty() {
                continue;
            }
            let header = current_header.as_deref().ok_or_else(|| {
                anyhow!(
                    "FASTA sequence encountered before first header at line {} in {}",
                    line_no + 1,
                    path.display()
                )
            })?;
            validate_sequence_line(sequence_line, header, current_len)?;
            current_len += sequence_line.len();
        }
    }
    if let Some(header) = current_header.take() {
        record_order.push(store_input_record(&mut genomes, &header, current_len)?);
    }
    if genomes.is_empty() {
        bail!("reference FASTA contains no records");
    }

    let mut declared_marker_sizes = 0usize;
    for genome in genomes.values_mut() {
        genome.reference_type = types.of(&genome.genome);
        let sequence_sum: usize = genome.contigs.iter().map(|c| c.length).sum();
        if sequence_sum == 0 {
            bail!("genome {} has no sequence bases", genome.genome);
        }
        genome.full_genome_size = match genome.provided_genome_size {
            Some(size) if size < sequence_sum => bail!(
                "provided genome size {} for {} is smaller than summed contig length {}",
                size,
                genome.genome,
                sequence_sum
            ),
            // Marker read counts are reported as mapped, so a marker genome's
            // size is its summed marker length; the declared size of the whole
            // genome is kept beside it.
            Some(size)
                if size > sequence_sum && genome.reference_type == ReferenceType::Markers =>
            {
                declared_marker_sizes += 1;
                genome.declared_genome_size = Some(size);
                sequence_sum
            }
            Some(size) => {
                if size > sequence_sum {
                    eprintln!(
                        "[index] warning: provided genome size {} for {} exceeds summed contig length {}; profiling will treat this as a partially indexed genome",
                        size,
                        genome.genome,
                        sequence_sum
                    );
                }
                size
            }
            None => sequence_sum,
        };
    }
    if declared_marker_sizes > 0 {
        eprintln!(
            "[index] note: {declared_marker_sizes} marker genomes declare a genome size larger than their marker genes; it is kept for the unclassified estimate, and their read counts and depth use the marker genes"
        );
    }
    let unlisted = types
        .by_genome
        .keys()
        .filter(|genome| !genomes.contains_key(genome.as_str()))
        .count();
    if unlisted > 0 {
        eprintln!(
            "[index] note: {unlisted} genomes of the reference type table are not in the reference FASTA"
        );
    }

    Ok(ReferenceInput {
        genomes: genomes.into_values().collect(),
        record_order,
    })
}

fn open_text_reader(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    if path.extension().map(|ext| ext == "gz").unwrap_or(false) {
        Ok(Box::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

fn store_input_record(
    genomes: &mut IndexMap<String, GenomeInput>,
    header: &str,
    length: usize,
) -> Result<RecordPosition> {
    if header.is_empty() {
        bail!("empty FASTA header");
    }
    let parsed = parse_input_header(header)?;
    if length == 0 {
        bail!(
            "reference {}|{} has empty sequence",
            parsed.genome,
            parsed.contig
        );
    }

    let genome_index = match genomes.entry(parsed.genome.clone()) {
        indexmap::map::Entry::Occupied(entry) => entry.index(),
        indexmap::map::Entry::Vacant(entry) => {
            let index = entry.index();
            entry.insert(GenomeInput {
                genome: parsed.genome.clone(),
                taxonomy: parsed.taxonomy.clone(),
                provided_genome_size: parsed.genome_size,
                contigs: Vec::new(),
                full_genome_size: 0,
                reference_type: ReferenceType::Genome,
                declared_genome_size: None,
            });
            index
        }
    };
    let entry = &mut genomes[genome_index];
    if entry.taxonomy != parsed.taxonomy {
        bail!(
            "inconsistent taxonomy for genome {}: '{}' vs '{}'",
            parsed.genome,
            entry.taxonomy.raw,
            parsed.taxonomy.raw
        );
    }
    match (entry.provided_genome_size, parsed.genome_size) {
        (Some(a), Some(b)) if a != b => bail!(
            "inconsistent genome size for {}: {} vs {}",
            parsed.genome,
            a,
            b
        ),
        (None, Some(size)) => entry.provided_genome_size = Some(size),
        _ => {}
    }
    if entry.contigs.iter().any(|c| c.contig == parsed.contig) {
        bail!(
            "duplicate contig id '{}' within genome '{}'",
            parsed.contig,
            parsed.genome
        );
    }
    entry.contigs.push(InputContig {
        contig: parsed.contig,
        reference_name: header.to_string(),
        length,
    });
    Ok(RecordPosition {
        genome: u32::try_from(genome_index)
            .context("reference FASTA contains more than 2^32 genomes")?,
        contig: u32::try_from(entry.contigs.len() - 1)
            .with_context(|| format!("genome {} has more than 2^32 records", entry.genome))?,
    })
}

struct ParsedHeader {
    genome: String,
    contig: String,
    taxonomy: TaxonomySpec,
    genome_size: Option<usize>,
}

fn parse_input_header(header: &str) -> Result<ParsedHeader> {
    let fields: Vec<&str> = header.split('|').collect();
    if fields.len() != 3 && fields.len() != 4 {
        bail!(
            "malformed header '{}'; expected genome|contig|taxonomy or genome|contig|taxonomy|genome_size",
            header
        );
    }
    let genome = required_id_field(fields[0], "genome accession", header)?;
    let contig = required_id_field(fields[1], "contig id", header)?;
    let taxonomy = parse_taxonomy_field(fields[2])
        .with_context(|| format!("malformed taxonomy field in header '{}'", header))?;
    let genome_size = if fields.len() == 4 && !fields[3].trim().is_empty() {
        Some(parse_positive_usize(fields[3].trim(), "genome size")?)
    } else {
        None
    };
    Ok(ParsedHeader {
        genome,
        contig,
        taxonomy,
        genome_size,
    })
}

fn required_id_field(value: &str, label: &str, header: &str) -> Result<String> {
    let value = value.trim();
    if value.is_empty() {
        bail!("empty {label} in header '{}'", header);
    }
    if value.contains('|') || value.contains('\t') || value.chars().any(char::is_whitespace) {
        bail!(
            "{label} '{}' in header '{}' must not contain whitespace, tab, or pipe",
            value,
            header
        );
    }
    Ok(value.to_string())
}

fn parse_taxonomy_field(value: &str) -> Result<TaxonomySpec> {
    let raw = value.trim();
    if raw.is_empty() {
        bail!("taxonomy field is empty");
    }
    if raw.contains('\t') || raw.chars().any(char::is_whitespace) {
        bail!("taxonomy field must not contain whitespace or tabs");
    }
    let mut seen: HashMap<char, u32> = HashMap::new();
    for token in raw.split(';') {
        if token.is_empty() {
            bail!("empty taxonomy token in '{}'", raw);
        }
        let (prefix, taxid) = token
            .split_once(':')
            .ok_or_else(|| anyhow!("taxonomy token '{}' is missing ':'", token))?;
        let mut chars = prefix.chars();
        let rank = chars
            .next()
            .ok_or_else(|| anyhow!("empty taxonomy prefix in '{}'", token))?;
        if chars.next().is_some() || !matches!(rank, 't' | 's' | 'g' | 'f' | 'o' | 'c' | 'p' | 'k')
        {
            bail!("unsupported taxonomy prefix '{}' in '{}'", prefix, raw);
        }
        if seen.contains_key(&rank) {
            bail!("duplicate taxonomy prefix '{}' in '{}'", rank, raw);
        }
        seen.insert(rank, parse_positive_u32(taxid, "taxid")?);
    }
    if seen.contains_key(&'t') && !seen.contains_key(&'s') {
        bail!(
            "taxonomy field '{}' has t:<strain_taxid> but no s:<species_taxid>",
            raw
        );
    }
    let profile_taxid = if let Some(species) = seen.get(&'s') {
        *species
    } else {
        ['g', 'f', 'o', 'c', 'p', 'k']
            .iter()
            .find_map(|rank| seen.get(rank).copied())
            .ok_or_else(|| anyhow!("taxonomy field '{}' contains no supported taxid", raw))?
    };
    Ok(TaxonomySpec {
        raw: raw.to_string(),
        primary_taxid: seen.get(&'t').copied(),
        profile_taxid,
    })
}

fn parse_positive_u32(value: &str, label: &str) -> Result<u32> {
    let parsed: u32 = value
        .parse()
        .with_context(|| format!("{label} '{}' is not an unsigned integer", value))?;
    if parsed == 0 {
        bail!("{label} must be greater than zero");
    }
    Ok(parsed)
}

fn parse_positive_usize(value: &str, label: &str) -> Result<usize> {
    let parsed: usize = value
        .parse()
        .with_context(|| format!("{label} '{}' is not an unsigned integer", value))?;
    if parsed == 0 {
        bail!("{label} must be greater than zero");
    }
    Ok(parsed)
}

fn validate_sequence_line(sequence: &str, header: &str, offset: usize) -> Result<()> {
    for (idx, base) in sequence.bytes().enumerate() {
        normalize_base(base).with_context(|| {
            format!(
                "invalid base '{}' at offset {} in reference '{}'",
                base as char,
                offset + idx,
                header
            )
        })?;
    }
    Ok(())
}

fn normalize_base(base: u8) -> Result<u8> {
    match base.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' | b'N' => Ok(base.to_ascii_uppercase()),
        b'U' => Ok(b'T'),
        b'R' | b'Y' | b'K' | b'M' | b'S' | b'W' | b'B' | b'D' | b'H' | b'V' => Ok(b'N'),
        _ => bail!("not a DNA/IUPAC base"),
    }
}

fn validate_subsample_params(params: SubsampleParams, threads: Option<usize>) -> Result<()> {
    if !(0.0 < params.fraction && params.fraction < 1.0) {
        bail!("subsample fraction must be between 0 and 1 (exclusive)");
    }
    if params.segment_length == 0 {
        bail!("segment length must be greater than zero");
    }
    if params.min_length == 0 {
        bail!("minimum length must be greater than zero");
    }
    if matches!(threads, Some(0)) {
        bail!("thread count must be at least 1");
    }
    Ok(())
}

fn select_index_records(
    input: &ReferenceInput,
    subsample: Option<SubsampleParams>,
    threads: Option<usize>,
) -> Result<Vec<SelectedRecord>> {
    let genomes = input.genomes.as_slice();
    if genomes.is_empty() {
        return Ok(Vec::new());
    }
    let params = match subsample {
        Some(params) => params,
        None => {
            // Whole records in input order: the index is built from the input
            // FASTA as it is, so its references come in this order.
            let effective: Vec<usize> = genomes
                .iter()
                .map(|genome| genome.contigs.iter().map(|c| c.length).sum())
                .collect();
            return Ok(input
                .record_order
                .iter()
                .map(|position| {
                    let genome = &genomes[position.genome as usize];
                    let contig = &genome.contigs[position.contig as usize];
                    SelectedRecord {
                        genome: genome.genome.clone(),
                        contig: contig.contig.clone(),
                        source_contig: contig.contig.clone(),
                        source_name: contig.reference_name.clone(),
                        taxonomy: genome.taxonomy.clone(),
                        reference_type: genome.reference_type,
                        full_genome_size: genome.full_genome_size,
                        effective_indexed_size: effective[position.genome as usize],
                        declared_genome_size: genome.declared_genome_size,
                        source_start: 1,
                        source_end: contig.length,
                    }
                })
                .collect());
        }
    };

    let mut seed_rng = StdRng::seed_from_u64(params.seed);
    let genome_seeds: Vec<u64> = genomes.iter().map(|_| seed_rng.gen()).collect();
    let sequential = || -> Result<Vec<Vec<SelectedRecord>>> {
        genomes
            .iter()
            .zip(genome_seeds.iter())
            .map(|(genome, &seed)| {
                let mut rng = StdRng::seed_from_u64(seed);
                subsample_genome(genome, params, &mut rng)
            })
            .collect()
    };
    let per_genome = if let Some(thread_count) = threads {
        if thread_count > 1 {
            let pool = ThreadPoolBuilder::new()
                .num_threads(thread_count)
                .build()
                .with_context(|| {
                    format!("failed to build thread pool with {thread_count} threads")
                })?;
            let mut results = pool.install(|| -> Result<Vec<(usize, Vec<SelectedRecord>)>> {
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
            results.into_iter().map(|(_, records)| records).collect()
        } else {
            sequential()?
        }
    } else {
        sequential()?
    };
    Ok(per_genome.into_iter().flatten().collect())
}

fn subsample_genome<R: rand::Rng + ?Sized>(
    genome: &GenomeInput,
    params: SubsampleParams,
    rng: &mut R,
) -> Result<Vec<SelectedRecord>> {
    const SHORT_GENOME_THRESHOLD: usize = 500_000;

    let mut segments = Vec::new();
    // Short genomes, and marker genes of any size, are indexed whole.
    if genome.reference_type == ReferenceType::Markers
        || genome.full_genome_size < SHORT_GENOME_THRESHOLD
    {
        for (idx, contig) in genome.contigs.iter().enumerate() {
            let len = contig.length;
            if len > 0 {
                segments.push(Segment {
                    contig_index: idx,
                    start: 1,
                    end: len,
                });
            }
        }
    } else {
        let target_len =
            ((params.fraction * genome.full_genome_size as f64).ceil() as usize).max(1);
        let mut contig_candidates: Vec<Vec<Segment>> = Vec::new();
        let mut total_slots = 0usize;
        for (idx, contig) in genome.contigs.iter().enumerate() {
            let len = contig.length;
            if len < params.segment_length {
                continue;
            }
            let max_segments = len / params.segment_length;
            if max_segments == 0 {
                continue;
            }
            let partition_size = len / max_segments;
            let mut contig_segments = Vec::new();
            let mut partition_start = 0usize;
            for partition in 0..max_segments {
                let partition_end = if partition == max_segments - 1 {
                    len
                } else {
                    (partition + 1) * partition_size
                };
                if partition_end <= partition_start {
                    continue;
                }
                let available = partition_end.saturating_sub(partition_start);
                if available < params.segment_length {
                    partition_start = partition_end;
                    continue;
                }
                let max_start = partition_end - params.segment_length;
                let start = rng.gen_range(partition_start..=max_start) + 1;
                let end = start + params.segment_length - 1;
                contig_segments.push(Segment {
                    contig_index: idx,
                    start,
                    end,
                });
                partition_start = partition_end;
            }
            if !contig_segments.is_empty() {
                contig_segments.sort_by_key(|seg| seg.start);
                total_slots += contig_segments.len();
                contig_candidates.push(contig_segments);
            }
        }
        if total_slots == 0 || total_slots * params.segment_length < target_len {
            segments = select_longest_contigs(genome, target_len, params);
        } else {
            let target_segments = target_len.div_ceil(params.segment_length);
            let mut allocations = vec![0usize; contig_candidates.len()];
            if target_segments >= total_slots {
                for (i, segs) in contig_candidates.iter().enumerate() {
                    allocations[i] = segs.len();
                }
            } else {
                let total_slots_f = total_slots as f64;
                for (i, segs) in contig_candidates.iter().enumerate() {
                    let share = segs.len() as f64 / total_slots_f;
                    allocations[i] =
                        ((share * target_segments as f64).floor() as usize).min(segs.len());
                }
                let mut remainder = target_segments.saturating_sub(allocations.iter().sum());
                let mut fractional: Vec<(usize, f64)> = contig_candidates
                    .iter()
                    .enumerate()
                    .map(|(i, segs)| {
                        let share = segs.len() as f64 / total_slots_f;
                        let exact = share * target_segments as f64;
                        (i, exact - allocations[i] as f64)
                    })
                    .collect();
                fractional
                    .sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
                for (idx, _) in fractional {
                    if remainder == 0 {
                        break;
                    }
                    if allocations[idx] < contig_candidates[idx].len() {
                        allocations[idx] += 1;
                        remainder -= 1;
                    }
                }
                if remainder > 0 {
                    for (i, segs) in contig_candidates.iter().enumerate() {
                        while allocations[i] < segs.len() && remainder > 0 {
                            allocations[i] += 1;
                            remainder -= 1;
                        }
                    }
                }
            }
            for (allocation, segs) in allocations.into_iter().zip(contig_candidates.iter()) {
                if allocation == 0 {
                    continue;
                }
                if allocation >= segs.len() {
                    segments.extend(segs.iter().cloned());
                } else {
                    for j in 0..allocation {
                        let idx = ((j + 1) * segs.len()) / allocation - 1;
                        segments.push(segs[idx].clone());
                    }
                }
            }
        }
    }

    segments.sort_by_key(|seg| (seg.contig_index, seg.start));
    let segments = enforce_min_segment_length(segments, genome, params.min_length);
    let effective_indexed_size: usize = segments.iter().map(Segment::len).sum();
    if effective_indexed_size == 0 {
        bail!("subsampling selected no bases for genome {}", genome.genome);
    }
    let mut records = Vec::new();
    for segment in segments {
        let contig = &genome.contigs[segment.contig_index];
        let start = segment.start;
        let end = segment.end.min(contig.length);
        if start == 0 || start > end {
            continue;
        }
        records.push(SelectedRecord {
            genome: genome.genome.clone(),
            contig: format!("{}:{}-{}", contig.contig, segment.start, segment.end),
            source_contig: contig.contig.clone(),
            source_name: contig.reference_name.clone(),
            taxonomy: genome.taxonomy.clone(),
            reference_type: genome.reference_type,
            full_genome_size: genome.full_genome_size,
            effective_indexed_size,
            declared_genome_size: genome.declared_genome_size,
            source_start: start,
            source_end: end,
        });
    }
    Ok(records)
}

fn select_longest_contigs(
    genome: &GenomeInput,
    target_len: usize,
    params: SubsampleParams,
) -> Vec<Segment> {
    let mut indices: Vec<usize> = (0..genome.contigs.len()).collect();
    indices.sort_by_key(|&idx| Reverse(genome.contigs[idx].length));
    let apply_filter = genome.full_genome_size > params.min_length.saturating_mul(10);
    let mut accumulated = 0usize;
    let mut segments = Vec::new();
    let mut fallback = Vec::new();
    for idx in indices {
        if accumulated >= target_len {
            break;
        }
        let len = genome.contigs[idx].length;
        if len == 0 {
            continue;
        }
        let remaining = target_len.saturating_sub(accumulated);
        let end = len.min(remaining.max(1));
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
        for segment in fallback {
            if accumulated >= target_len {
                break;
            }
            accumulated += segment.len();
            segments.push(segment);
        }
    }
    segments
}

fn enforce_min_segment_length(
    segments: Vec<Segment>,
    genome: &GenomeInput,
    min_length: usize,
) -> Vec<Segment> {
    if segments.is_empty() || genome.full_genome_size <= min_length.saturating_mul(10) {
        return segments;
    }
    let filtered: Vec<Segment> = segments
        .iter()
        .filter(|segment| segment.len() >= min_length)
        .cloned()
        .collect();
    if filtered.is_empty() {
        segments
    } else {
        filtered
    }
}

fn build_direct_reference_metadata(records: Vec<SelectedRecord>) -> Result<Vec<ReferenceMetadata>> {
    let mut seen_names = HashSet::new();
    let mut references = Vec::with_capacity(records.len());
    for record in records {
        let length = record.len();
        if !seen_names.insert(record.source_name.clone()) {
            bail!(
                "reference FASTA contains duplicate sequence ID {}",
                record.source_name
            );
        }
        references.push(ReferenceMetadata {
            name: record.source_name,
            length,
            genome: record.genome,
            contig: record.contig,
            full_genome_size: record.full_genome_size,
            effective_indexed_size: record.effective_indexed_size,
            taxonomy_field: record.taxonomy.raw,
            reference_type: record.reference_type,
            declared_genome_size: record.declared_genome_size,
        });
    }
    Ok(references)
}

fn write_metam_reference(
    path: &Path,
    source_fasta: &Path,
    records: Vec<SelectedRecord>,
    force_sampled_field: bool,
) -> Result<Vec<ReferenceMetadata>> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    let mut seen_names = HashSet::new();
    let record_count = records.len();
    let mut selected_by_source: HashMap<(String, String), Vec<PreparedRecord>> = HashMap::new();
    for record in records {
        let name = metax_reference_name(&record, force_sampled_field);
        if !seen_names.insert(name.clone()) {
            bail!("generated duplicate reference name '{}'", name);
        }
        selected_by_source
            .entry((record.genome.clone(), record.source_contig.clone()))
            .or_default()
            .push(PreparedRecord { name, record });
    }
    drop(seen_names);

    // Metadata rows are added as the records are written, so they follow the
    // written FASTA (source order), which is the order the index stores.
    let mut references = Vec::with_capacity(record_count);
    stream_selected_reference_records(
        source_fasta,
        &mut selected_by_source,
        &mut writer,
        &mut references,
    )?;
    writer.flush()?;
    if !selected_by_source.is_empty() {
        let missing = selected_by_source
            .keys()
            .take(5)
            .map(|(genome, contig)| format!("{genome}|{contig}"))
            .collect::<Vec<_>>()
            .join(", ");
        bail!(
            "selected reference records were not found while streaming source FASTA: {}",
            missing
        );
    }
    Ok(references)
}

#[derive(Clone)]
struct PreparedRecord {
    name: String,
    record: SelectedRecord,
}

fn stream_selected_reference_records(
    source_fasta: &Path,
    selected_by_source: &mut HashMap<(String, String), Vec<PreparedRecord>>,
    writer: &mut dyn Write,
    references: &mut Vec<ReferenceMetadata>,
) -> Result<()> {
    let reader = open_text_reader(source_fasta)?;
    let reader = BufReader::new(reader);
    let mut current_key: Option<(String, String)> = None;
    let mut current_header: Option<String> = None;
    let mut capture = false;
    let mut sequence = String::new();

    for line in reader.lines() {
        let line =
            line.with_context(|| format!("failed reading FASTA {}", source_fasta.display()))?;
        if let Some(rest) = line.strip_prefix('>') {
            flush_selected_source_record(
                selected_by_source,
                writer,
                references,
                current_key.take(),
                current_header.take(),
                &sequence,
            )?;
            sequence.clear();
            let header = rest.trim().to_string();
            let parsed = parse_input_header(&header)?;
            let key = (parsed.genome, parsed.contig);
            capture = selected_by_source.contains_key(&key);
            current_key = Some(key);
            current_header = Some(header);
        } else if capture {
            sequence.push_str(line.trim());
        }
    }
    flush_selected_source_record(
        selected_by_source,
        writer,
        references,
        current_key.take(),
        current_header.take(),
        &sequence,
    )?;
    Ok(())
}

fn flush_selected_source_record(
    selected_by_source: &mut HashMap<(String, String), Vec<PreparedRecord>>,
    writer: &mut dyn Write,
    references: &mut Vec<ReferenceMetadata>,
    key: Option<(String, String)>,
    header: Option<String>,
    sequence: &str,
) -> Result<()> {
    let Some(key) = key else {
        return Ok(());
    };
    let Some(records) = selected_by_source.remove(&key) else {
        return Ok(());
    };
    let header = header.unwrap_or_else(|| format!("{}|{}", key.0, key.1));
    for prepared in records {
        if prepared.record.source_start == 0
            || prepared.record.source_end < prepared.record.source_start
        {
            bail!(
                "invalid selected interval {}:{}-{} for reference '{}'",
                prepared.record.source_contig,
                prepared.record.source_start,
                prepared.record.source_end,
                prepared.name
            );
        }
        if prepared.record.source_end > sequence.len() {
            bail!(
                "selected interval {}:{}-{} exceeds source length {} for FASTA record '{}'",
                prepared.record.source_contig,
                prepared.record.source_start,
                prepared.record.source_end,
                sequence.len(),
                header
            );
        }
        let start = prepared.record.source_start - 1;
        let end = prepared.record.source_end;
        writeln!(writer, ">{}", prepared.name)?;
        write_normalized_sequence(writer, &sequence.as_bytes()[start..end])?;
        let length = prepared.record.len();
        let record = prepared.record;
        references.push(ReferenceMetadata {
            name: prepared.name,
            length,
            genome: record.genome,
            contig: record.contig,
            full_genome_size: record.full_genome_size,
            effective_indexed_size: record.effective_indexed_size,
            taxonomy_field: record.taxonomy.raw,
            reference_type: record.reference_type,
            declared_genome_size: record.declared_genome_size,
        });
    }
    Ok(())
}

fn write_normalized_sequence(writer: &mut dyn Write, sequence: &[u8]) -> Result<()> {
    let mut line_len = 0usize;
    for &base in sequence {
        let normalized = normalize_base(base)?;
        writer.write_all(&[normalized])?;
        line_len += 1;
        if line_len == 60 {
            writer.write_all(b"\n")?;
            line_len = 0;
        }
    }
    if line_len > 0 {
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn metax_reference_name(record: &SelectedRecord, force_sampled_field: bool) -> String {
    let primary = record
        .taxonomy
        .primary_taxid
        .map(|v| v.to_string())
        .unwrap_or_else(|| "0".to_string());
    let mut name = format!(
        "{}|{}|{}|{}|{}",
        record.genome,
        primary,
        record.taxonomy.profile_taxid,
        record.contig,
        record.full_genome_size
    );
    if force_sampled_field {
        name.push('|');
        name.push_str(&record.effective_indexed_size.to_string());
    }
    name
}

fn run_metam_index(
    reference_path: &Path,
    index_path: &Path,
    threads: usize,
    options: &MetamIndexOptions,
) -> Result<IndexedReferences> {
    let seed_mode = parse_seed_mode(&options.seed_mode)?;
    let strobe_start_dist =
        StrobeStartDistance::new(options.strobe_start_dist_min, options.strobe_start_dist_max)?;
    let seed_chunk_bases = seed_chunk_mib_to_bases(options.seed_chunk_mib)?;
    let progress_options = ProgressOptions::new(
        options.progress,
        options.no_progress || options.quiet,
        options.progress_interval,
    );
    let mut progress = IndexProgress::new(progress_options);
    let reference = reference_path.to_string_lossy().to_string();
    let index = index_path.to_string_lossy().to_string();
    let default_temp_parent = index_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let temp_parent = options.temp_dir.as_deref().unwrap_or(default_temp_parent);
    if !options.quiet {
        log_message(
            "index",
            format_args!("temporary seed workspace parent={}", temp_parent.display()),
        );
    }
    let (header, stats) = build_index_file_with_progress(
        &reference,
        &index,
        options.k,
        options.w,
        options.max_occ,
        options.hash_bits,
        seed_mode,
        options.match_both_strobes,
        Some(strobe_start_dist),
        threads,
        seed_chunk_bases,
        Some(temp_parent),
        |event| match event {
            BuildProgress::Phase(phase) => progress.phase(phase),
            BuildProgress::ReferencesRead(count) => progress.references_read(count as u64),
            BuildProgress::ReferencesSeeded(done, total) => {
                progress.references_seeded(done as u64, total as u64)
            }
            BuildProgress::Writing => progress.writing(),
        },
    )?;
    progress.finish(header.reference_count as u64);
    if !options.quiet {
        log_message(
            "index",
            format_args!(
                "raw_seed_occurrences={} indexed_seed_occurrences={} filtered_hashes={} max_partition_seed_occurrences={} temp_partitions={} seed_chunk_mib={}",
                stats.raw_seed_occurrences,
                header.indexed_seed_occurrences,
                stats.filtered_hashes,
                stats.max_partition_seed_occurrences,
                stats.temp_partitions,
                stats.seed_chunk_bases / (1024 * 1024)
            ),
        );
    }
    Ok(IndexedReferences {
        count: header.reference_count,
        order_digest: stats.reference_order_digest,
    })
}

/// `metax profile` addresses references by their position in the index and
/// requires the database metadata to list them in the same order. Checking it
/// here makes a mismatch fail the build, before a database is published.
fn verify_index_reference_order(
    references: &[ReferenceMetadata],
    indexed: IndexedReferences,
) -> Result<()> {
    let mut digest = ReferenceOrderDigest::new();
    for reference in references {
        digest.push(&reference.name);
    }
    if references.len() != indexed.count || digest.finish() != indexed.order_digest {
        bail!(
            "internal error: the mapping index ({} references) and the database metadata ({} rows) do not list the same references in the same order; the database was not written",
            indexed.count,
            references.len()
        );
    }
    Ok(())
}

fn seed_chunk_mib_to_bases(seed_chunk_mib: usize) -> Result<usize> {
    if seed_chunk_mib == 0 {
        bail!("--seed-chunk-mib must be greater than 0");
    }
    seed_chunk_mib
        .checked_mul(1024 * 1024)
        .context("--seed-chunk-mib is too large for this platform")
}

fn parse_seed_mode(value: &str) -> Result<SeedMode> {
    match value.to_ascii_lowercase().as_str() {
        "minimizer" => Ok(SeedMode::Minimizer),
        "randstrobe" => Ok(SeedMode::Randstrobe),
        _ => bail!(
            "unsupported seed mode '{}'; expected minimizer or randstrobe",
            value
        ),
    }
}

struct IndexProgress {
    options: ProgressOptions,
    spinner: indicatif::ProgressBar,
    counting: bool,
}

impl IndexProgress {
    fn new(options: ProgressOptions) -> Self {
        Self {
            options,
            spinner: progress::spinner(options, "Starting index build"),
            counting: false,
        }
    }

    fn phase(&mut self, phase: &'static str) {
        self.set_spinner(phase);
    }

    fn references_read(&mut self, count: u64) {
        self.count("Reading references", count);
    }

    fn references_seeded(&mut self, done: u64, _total: u64) {
        self.count("Seeding references", done);
    }

    fn writing(&mut self) {
        self.set_spinner("Writing index");
    }

    fn finish(&mut self, count: u64) {
        self.spinner
            .set_message(format!("Index build complete: {count} references"));
        self.spinner.finish_and_clear();
    }

    fn set_spinner(&mut self, message: &'static str) {
        if self.counting {
            self.spinner.finish_and_clear();
            self.spinner = progress::spinner(self.options, message);
            self.counting = false;
        } else {
            self.spinner.set_message(message);
            self.spinner.tick();
        }
    }

    fn count(&mut self, label: &'static str, count: u64) {
        if !self.counting {
            self.spinner.finish_and_clear();
            self.spinner = progress::counter(self.options, label);
            self.counting = true;
        }
        progress::tick_count(&self.spinner, label, count, self.options.interval());
    }
}

fn gzip_file(path: &Path) -> Result<PathBuf> {
    let gz_path = path.with_extension(format!(
        "{}.gz",
        path.extension()
            .map(|ext| ext.to_string_lossy().to_string())
            .unwrap_or_else(|| "gz".to_string())
    ));
    let mut cleanup_output = RemoveFileOnDrop::new(gz_path.clone());
    let input =
        File::open(path).with_context(|| format!("failed to open {} for gzip", path.display()))?;
    let mut reader = BufReader::new(input);
    let output = File::create(&gz_path)
        .with_context(|| format!("failed to create {}", gz_path.display()))?;
    let mut encoder = GzEncoder::new(output, Compression::fast());
    std::io::copy(&mut reader, &mut encoder)?;
    encoder.finish()?;
    fs::remove_file(path)
        .with_context(|| format!("failed to remove {} after gzip", path.display()))?;
    cleanup_output.disarm();
    Ok(gz_path)
}

fn parse_bool(value: &str) -> Result<bool> {
    match value {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => bail!("expected true or false"),
    }
}

fn ensure_tsv_safe(value: &str, label: &str) -> Result<()> {
    if value.contains('\t') || value.contains('\n') || value.contains('\r') {
        bail!("{label} contains a tab or newline: '{}'", value);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(header: &str) -> ParsedHeader {
        parse_input_header(header).unwrap()
    }

    fn test_index_options(temp_dir: Option<PathBuf>) -> MetamIndexOptions {
        MetamIndexOptions {
            k: 11,
            w: 11,
            max_occ: 100,
            hash_bits: 32,
            seed_mode: "randstrobe".into(),
            match_both_strobes: true,
            strobe_start_dist_min: 20,
            strobe_start_dist_max: 60,
            seed_chunk_mib: 1,
            temp_dir,
            quiet: true,
            progress: false,
            no_progress: true,
            progress_interval: 1_000,
        }
    }

    #[test]
    fn parses_supported_taxonomy_shapes() {
        let h = parse("GCF_1|contig1|t:10;s:20|1000");
        assert_eq!(h.genome, "GCF_1");
        assert_eq!(h.contig, "contig1");
        assert_eq!(h.taxonomy.primary_taxid, Some(10));
        assert_eq!(h.taxonomy.profile_taxid, 20);
        assert_eq!(h.genome_size, Some(1000));

        let h = parse("GCF_2|c|s:562");
        assert_eq!(h.taxonomy.primary_taxid, None);
        assert_eq!(h.taxonomy.profile_taxid, 562);

        let h = parse("GCF_3|c|g:561|");
        assert_eq!(h.taxonomy.primary_taxid, None);
        assert_eq!(h.taxonomy.profile_taxid, 561);
    }

    #[test]
    fn rejects_malformed_headers_and_taxonomy() {
        assert!(parse_input_header("G|c").is_err());
        assert!(parse_input_header("G|c|x:1").is_err());
        assert!(parse_input_header("G|c|t:1").is_err());
        assert!(parse_input_header("G|c|s:0").is_err());
        assert!(parse_input_header("G|c|s:1|abc").is_err());
    }

    #[test]
    fn generated_reference_names_match_profiler_layout() {
        let record = SelectedRecord {
            genome: "G".into(),
            contig: "ctg:1-10".into(),
            source_contig: "ctg".into(),
            source_name: "G|ctg|t:7;s:8|100".into(),
            taxonomy: TaxonomySpec {
                raw: "t:7;s:8".into(),
                primary_taxid: Some(7),
                profile_taxid: 8,
            },
            reference_type: ReferenceType::Genome,
            full_genome_size: 100,
            effective_indexed_size: 10,
            declared_genome_size: None,
            source_start: 1,
            source_end: 10,
        };
        assert_eq!(metax_reference_name(&record, true), "G|7|8|ctg:1-10|100|10");
    }

    #[test]
    fn partial_input_without_explicit_subsampling_still_records_effective_size() {
        let mut record = SelectedRecord {
            genome: "G".into(),
            contig: "ctg".into(),
            source_contig: "ctg".into(),
            source_name: "G|ctg|s:8|100".into(),
            taxonomy: TaxonomySpec {
                raw: "s:8".into(),
                primary_taxid: None,
                profile_taxid: 8,
            },
            reference_type: ReferenceType::Genome,
            full_genome_size: 100,
            effective_indexed_size: 60,
            declared_genome_size: None,
            source_start: 1,
            source_end: 60,
        };
        assert!(requires_effective_size_reference_names(
            std::slice::from_ref(&record),
            false
        ));
        assert_eq!(metax_reference_name(&record, true), "G|0|8|ctg|100|60");

        record.effective_indexed_size = 100;
        assert!(!requires_effective_size_reference_names(
            std::slice::from_ref(&record),
            false
        ));
        assert!(requires_effective_size_reference_names(&[record], true));
    }

    #[test]
    fn db_tsv_v2_has_comment_metadata_and_plain_reference_rows() {
        let metadata = DatabaseMetadata {
            metam_index: PathBuf::from("metax_db.mxi"),
            is_fractional: false,
            references: vec![ReferenceMetadata {
                name: "GCF_1|ctg1|s:20|1000".into(),
                length: 500,
                genome: "GCF_1".into(),
                contig: "ctg1".into(),
                full_genome_size: 1000,
                effective_indexed_size: 1000,
                taxonomy_field: "s:20".into(),
                reference_type: ReferenceType::Genome,
                declared_genome_size: None,
            }],
        };
        let mut out = Vec::new();
        write_db_tsv_body(&mut out, &metadata).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines = text.lines().collect::<Vec<_>>();

        assert_eq!(lines[0], "#metax_db_version\t2");
        assert_eq!(lines[1], "#mapping_index\tmetax_db.mxi");
        assert_eq!(lines[2], "#fractional\tfalse");
        assert_eq!(
            lines[3],
            "#columns\treference_name\treference_length\tfull_genome_size\teffective_indexed_size\tgenome\tcontig\ttaxonomy"
        );
        assert_eq!(
            lines[4],
            "GCF_1|ctg1|s:20|1000\t500\t1000\t1000\tGCF_1\tctg1\ts:20"
        );
        assert!(!lines[4].starts_with("reference\t"));
    }

    #[test]
    fn rejects_legacy_database_version() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("legacy.tsv");
        std::fs::write(
            &path,
            "#metax_db_version\t1\nmapping_index\tmetax_db.mxi\nfractional\tfalse\nreference\tG|c|s:20\t100\tG\tc\t100\t100\ts:20\n",
        )
        .unwrap();

        let error = load_db_tsv(&path).unwrap_err();
        assert!(error.to_string().contains("unsupported Metax DB version 1"));
    }

    #[test]
    fn rejects_unprefixed_legacy_metadata_keys() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("legacy-keys.tsv");
        std::fs::write(
            &path,
            "#metax_db_version\t2\nmapping_index\tmetax_db.mxi\n#fractional\tfalse\n",
        )
        .unwrap();

        let error = load_db_tsv(&path).unwrap_err();
        assert!(error
            .to_string()
            .contains("malformed database metadata line 2"));
    }

    #[test]
    fn rejects_prefixed_legacy_reference_rows() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("legacy-reference.tsv");
        std::fs::write(
            &path,
            "#metax_db_version\t2\n#mapping_index\tmetax_db.mxi\n#fractional\tfalse\nreference\tG|c|s:20\t100\tG\tc\t100\t100\ts:20\n",
        )
        .unwrap();

        let error = load_db_tsv(&path).unwrap_err();
        assert!(error
            .to_string()
            .contains("malformed database metadata line 4"));
    }

    #[test]
    fn failed_index_build_removes_database_staging_area() {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("reference.fa");
        let output = directory.path().join("database");
        std::fs::write(
            &fasta,
            format!(
                ">G1|contig1|s:20\n{}\n",
                "ACGTGCAATGCCGTAAGTCCGATCGTACGATGCTAGCTAGGCTA".repeat(10)
            ),
        )
        .unwrap();
        let mut options = test_index_options(Some(directory.path().to_path_buf()));
        options.max_occ = 0;

        let error = build_database(
            fasta.to_str().unwrap(),
            output.to_str().unwrap(),
            None,
            Some(2),
            false,
            &ReferenceTypes::default(),
            &options,
        )
        .unwrap_err();

        assert!(error.to_string().contains("--max-occ"));
        assert!(!output.exists());
        let names = std::fs::read_dir(directory.path())
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .collect::<Vec<_>>();
        assert!(names
            .iter()
            .all(|name| !name.to_string_lossy().contains("metax-staging")));
    }

    #[test]
    fn existing_output_directory_is_never_overwritten() {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("reference.fa");
        let output = directory.path().join("database");
        std::fs::write(&fasta, ">G1|contig1|s:20\nACGTACGTACGTACGT\n").unwrap();
        std::fs::create_dir(&output).unwrap();
        let sentinel = output.join("user-file");
        std::fs::write(&sentinel, b"keep me").unwrap();
        let options = test_index_options(Some(directory.path().to_path_buf()));

        let error = build_database(
            fasta.to_str().unwrap(),
            output.to_str().unwrap(),
            None,
            Some(1),
            false,
            &ReferenceTypes::default(),
            &options,
        )
        .unwrap_err();

        assert!(error.to_string().contains("already exists"));
        assert_eq!(std::fs::read(sentinel).unwrap(), b"keep me");
    }

    #[test]
    fn successful_database_publish_contains_only_final_artifacts() {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("reference.fa");
        let output = directory.path().join("database");
        std::fs::write(
            &fasta,
            format!(
                ">G1|contig1|s:20\n{}\n",
                "ACGTGCAATGCCGTAAGTCCGATCGTACGATGCTAGCTAGGCTA".repeat(20)
            ),
        )
        .unwrap();
        let options = test_index_options(None);

        let metadata_path = build_database(
            fasta.to_str().unwrap(),
            output.to_str().unwrap(),
            None,
            Some(2),
            false,
            &ReferenceTypes::default(),
            &options,
        )
        .unwrap();

        assert_eq!(metadata_path, output.join(DB_TSV_FILE));
        assert!(output.join(DB_TSV_FILE).is_file());
        assert!(output.join(METAM_INDEX_FILE).is_file());
        let mut names = std::fs::read_dir(&output)
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .collect::<Vec<_>>();
        names.sort();
        assert_eq!(
            names,
            vec![
                std::ffi::OsString::from(METAM_INDEX_FILE),
                std::ffi::OsString::from(DB_TSV_FILE)
            ]
        );
        let loaded = DatabaseMetadata::load(&output).unwrap();
        assert_eq!(loaded.references.len(), 1);
    }

    /// Deterministic pseudo-random DNA, so every record yields seeds.
    fn test_sequence(seed: u64, len: usize) -> String {
        let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state % 4) as usize] as char
            })
            .collect()
    }

    /// Genome G1 has two records with a G2 record between them.
    fn write_interleaved_fasta(path: &Path, g1_size: &str) {
        let records = [
            format!("G1|c1|t:11;s:20{g1_size}"),
            "G2|c1|s:30".to_string(),
            format!("G1|c2|t:11;s:20{g1_size}"),
        ];
        let fasta = records
            .iter()
            .enumerate()
            .map(|(i, header)| format!(">{header}\n{}\n", test_sequence(i as u64 + 1, 600)))
            .collect::<String>();
        std::fs::write(path, fasta).unwrap();
    }

    fn index_reference_ids(database: &Path) -> Vec<String> {
        let index =
            metam::index::read_index(database.join(METAM_INDEX_FILE).to_str().unwrap()).unwrap();
        index.references.into_iter().map(|r| r.id).collect()
    }

    fn build_interleaved(subsample: Option<SubsampleParams>, g1_size: &str) -> DatabaseMetadata {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("reference.fa");
        let output = directory.path().join("database");
        write_interleaved_fasta(&fasta, g1_size);
        build_database(
            fasta.to_str().unwrap(),
            output.to_str().unwrap(),
            subsample,
            Some(2),
            false,
            &ReferenceTypes::default(),
            &test_index_options(None),
        )
        .unwrap();
        let loaded = DatabaseMetadata::load(&output).unwrap();
        let names: Vec<String> = loaded.references.iter().map(|r| r.name.clone()).collect();
        assert_eq!(index_reference_ids(&output), names);
        let genomes: Vec<&str> = loaded
            .references
            .iter()
            .map(|r| r.genome.as_str())
            .collect();
        assert_eq!(genomes, ["G1", "G2", "G1"]);
        loaded
    }

    #[test]
    fn interleaved_genome_records_keep_input_order() {
        let loaded = build_interleaved(None, "");
        let names: Vec<&str> = loaded.references.iter().map(|r| r.name.as_str()).collect();
        assert_eq!(names, ["G1|c1|t:11;s:20", "G2|c1|s:30", "G1|c2|t:11;s:20"]);
        let g1 = &loaded.references[0];
        assert_eq!(
            (g1.full_genome_size, g1.effective_indexed_size),
            (1200, 1200)
        );
        assert_eq!(loaded.references[2].full_genome_size, 1200);
    }

    #[test]
    fn interleaved_genome_records_keep_input_order_when_fractional() {
        let params = SubsampleParams {
            fraction: 0.5,
            seed: 42,
            segment_length: 50_000,
            min_length: 3_000,
        };
        let loaded = build_interleaved(Some(params), "");
        assert!(loaded.is_fractional);
        let contigs: Vec<&str> = loaded
            .references
            .iter()
            .map(|r| r.contig.as_str())
            .collect();
        assert_eq!(contigs, ["c1:1-600", "c1:1-600", "c2:1-600"]);
    }

    #[test]
    fn interleaved_genome_records_keep_input_order_with_declared_size() {
        // A declared size above the summed length takes the rewritten-reference path.
        let loaded = build_interleaved(None, "|5000");
        assert!(loaded.is_fractional);
        let g1 = &loaded.references[0];
        assert_eq!(
            (g1.full_genome_size, g1.effective_indexed_size),
            (5000, 1200)
        );
    }

    #[test]
    fn index_order_check_rejects_a_different_order() {
        let reference = |name: &str| ReferenceMetadata {
            name: name.to_string(),
            length: 1,
            genome: "G".into(),
            contig: name.to_string(),
            full_genome_size: 1,
            effective_indexed_size: 1,
            taxonomy_field: "s:1".into(),
            reference_type: ReferenceType::Genome,
            declared_genome_size: None,
        };
        let references = vec![reference("a"), reference("b")];
        let digest = |names: &[&str]| {
            let mut digest = ReferenceOrderDigest::new();
            names.iter().for_each(|name| digest.push(name));
            digest.finish()
        };
        let indexed = |names: &[&str]| IndexedReferences {
            count: names.len(),
            order_digest: digest(names),
        };
        assert!(verify_index_reference_order(&references, indexed(&["a", "b"])).is_ok());
        assert!(verify_index_reference_order(&references, indexed(&["b", "a"])).is_err());
        assert!(verify_index_reference_order(&references, indexed(&["a"])).is_err());
    }

    fn build_typed(
        types: &str,
        subsample: Option<SubsampleParams>,
        g1_size: &str,
    ) -> (DatabaseMetadata, String) {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("reference.fa");
        let table = directory.path().join("types.tsv");
        let output = directory.path().join("database");
        write_interleaved_fasta(&fasta, g1_size);
        std::fs::write(&table, types).unwrap();
        let built = build_database(
            fasta.to_str().unwrap(),
            output.to_str().unwrap(),
            subsample,
            Some(2),
            false,
            &ReferenceTypes::new(ReferenceType::Genome, Some(&table)).unwrap(),
            &test_index_options(None),
        )
        .unwrap();
        let mut text = String::new();
        MultiGzDecoder::new(File::open(built).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        (DatabaseMetadata::load(&output).unwrap(), text)
    }

    #[test]
    fn genome_databases_keep_the_version_2_layout_and_fingerprint() {
        let (loaded, text) = build_typed("# no overrides\n", None, "");
        assert!(text.starts_with("#metax_db_version\t2\n"));
        assert!(!text.contains("reference_type"));
        assert_eq!(loaded.reference_type_label(), "genome");

        use std::fmt::Write as _;
        let mut hasher = xxhash_rust::xxh3::Xxh3::new();
        hasher.update(b"absolute\n");
        for r in &loaded.references {
            let mut row = String::new();
            writeln!(
                row,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.name,
                r.length,
                r.full_genome_size,
                r.effective_indexed_size,
                r.genome,
                r.contig,
                r.taxonomy_field
            )
            .unwrap();
            hasher.update(row.as_bytes());
        }
        assert_eq!(loaded.fingerprint(), hasher.digest128());
    }

    #[test]
    fn fingerprint_hashes_the_formatted_rows() {
        use std::fmt::Write as _;
        let reference =
            |name: &str, length, reference_type, declared_genome_size| ReferenceMetadata {
                name: name.to_string(),
                length,
                genome: "G1".to_string(),
                contig: format!("{name}_contig"),
                full_genome_size: 1_000_000_007,
                effective_indexed_size: 10,
                taxonomy_field: "t:436646856;s:2056026598".to_string(),
                reference_type,
                declared_genome_size,
            };
        for is_fractional in [false, true] {
            let database = DatabaseMetadata {
                metam_index: PathBuf::from("index.mxi"),
                references: vec![
                    reference("a", 1, ReferenceType::Genome, None),
                    reference("b", 987_654_321, ReferenceType::Markers, Some(2_540_361)),
                    reference("c", 9, ReferenceType::Markers, None),
                    reference("d", usize::MAX, ReferenceType::Genome, None),
                ],
                is_fractional,
            };
            let mut hasher = xxhash_rust::xxh3::Xxh3::new();
            hasher.update(if is_fractional {
                b"fractional\n"
            } else {
                b"absolute\n"
            });
            for r in &database.references {
                let mut row = String::new();
                write!(
                    row,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    r.name,
                    r.length,
                    r.full_genome_size,
                    r.effective_indexed_size,
                    r.genome,
                    r.contig,
                    r.taxonomy_field
                )
                .unwrap();
                if r.reference_type != ReferenceType::Genome {
                    write!(row, "\t{}", r.reference_type).unwrap();
                }
                if let Some(size) = r.declared_genome_size {
                    write!(row, "\t{size}").unwrap();
                }
                row.push('\n');
                hasher.update(row.as_bytes());
            }
            assert_eq!(database.fingerprint(), hasher.digest128());
        }
    }

    #[test]
    fn mixed_databases_record_reference_types() {
        // G1 declares a genome size: a marker genome keeps it beside its marker
        // length, which its read counts and depth use.
        let params = SubsampleParams {
            fraction: 0.5,
            seed: 7,
            segment_length: 100,
            min_length: 50,
        };
        let (loaded, text) = build_typed("G1\tMarkers\n", Some(params), "|5000");
        assert!(text.starts_with("#metax_db_version\t4\n"));
        assert!(text.contains("\treference_type\tdeclared_genome_size\n"));
        assert!(text.contains("\n#reference_type\tmixed\n"));
        assert_eq!(loaded.reference_type_label(), "mixed");
        for r in &loaded.references {
            let expected = if r.genome == "G1" {
                ReferenceType::Markers
            } else {
                ReferenceType::Genome
            };
            assert_eq!(r.reference_type, expected);
            if r.genome == "G1" {
                assert_eq!((r.full_genome_size, r.effective_indexed_size), (1200, 1200));
                assert_eq!(r.declared_genome_size, Some(5000));
            } else {
                assert_eq!(r.declared_genome_size, None);
            }
        }
        let mut without_sizes = loaded.clone();
        for r in &mut without_sizes.references {
            r.declared_genome_size = None;
        }
        assert_ne!(loaded.fingerprint(), without_sizes.fingerprint());
        let mut genomes_only = loaded.clone();
        for r in &mut genomes_only.references {
            r.reference_type = ReferenceType::Genome;
        }
        assert_ne!(loaded.fingerprint(), genomes_only.fingerprint());
    }

    #[test]
    fn a_marker_size_equal_to_the_marker_length_is_not_kept() {
        // Marker FASTAs that write the summed marker length as genome size
        // keep the version 3 layout.
        let (loaded, text) = build_typed("G1\tMarkers\n", None, "|1200");
        assert!(text.starts_with("#metax_db_version\t3\n"));
        assert!(loaded
            .references
            .iter()
            .all(|r| r.declared_genome_size.is_none()));
    }

    #[test]
    fn marker_genomes_are_never_subsampled() {
        let genome = GenomeInput {
            genome: "M".into(),
            taxonomy: parse_taxonomy_field("s:5").unwrap(),
            provided_genome_size: None,
            contigs: (0..600)
                .map(|i| InputContig {
                    contig: format!("gene{i}"),
                    reference_name: format!("M|gene{i}|s:5"),
                    length: 1_000,
                })
                .collect(),
            full_genome_size: 600_000,
            reference_type: ReferenceType::Markers,
            declared_genome_size: None,
        };
        let params = SubsampleParams {
            fraction: 0.1,
            seed: 1,
            segment_length: 50_000,
            min_length: 3_000,
        };
        let records = subsample_genome(&genome, params, &mut StdRng::seed_from_u64(1)).unwrap();
        assert_eq!(records.len(), 600);
        assert!(records
            .iter()
            .all(|r| r.len() == 1_000 && r.effective_indexed_size == 600_000));
    }

    #[test]
    fn genome_species_next_to_markers_are_found_per_genus() {
        let genus_of = |species: u32| (species < 100).then_some(species / 10);
        let species = [
            (ReferenceType::Genome, 11),
            (ReferenceType::Markers, 12),
            (ReferenceType::Genome, 21),
            (ReferenceType::Genome, 22),
            (ReferenceType::Markers, 31),
            (ReferenceType::Markers, 500),
            (ReferenceType::Genome, 600),
        ];
        let genera = genome_species_next_to_markers(species, genus_of);
        assert_eq!(genera.len(), 1);
        assert_eq!(genera[&1], HashSet::from([11]));
    }

    #[test]
    fn reference_type_tables_are_validated() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("types.tsv");
        std::fs::write(
            &path,
            "# genome\ttype\nA\tmarkers\nB\tGENOME\n\nA\tmarkers\n",
        )
        .unwrap();
        let types = ReferenceTypes::new(ReferenceType::Genome, Some(&path)).unwrap();
        assert_eq!(types.of("A"), ReferenceType::Markers);
        assert_eq!(types.of("B"), ReferenceType::Genome);
        assert_eq!(types.of("C"), ReferenceType::Genome);

        for bad in [
            "A\tmarkers\nA\tgenome\n",
            "A\tplasmid\n",
            "A markers\n",
            "A\tmarkers\textra\n",
        ] {
            std::fs::write(&path, bad).unwrap();
            assert!(
                ReferenceTypes::new(ReferenceType::Genome, Some(&path)).is_err(),
                "{bad:?}"
            );
        }
    }

    #[test]
    fn metadata_layouts_are_checked() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("metax_db.tsv");
        let header = |version: &str, extra: &str, columns: &str| {
            format!("#metax_db_version\t{version}\n#mapping_index\tmetax_db.mxi\n#fractional\tfalse\n{extra}#columns\t{columns}\n")
        };
        let row = "G|c|s:20\t100\t100\t100\tG\tc\ts:20";
        let typed_columns = format!("{COLUMNS_GENOMES}\treference_type");
        let cases = [
            // version 2 with a type column
            (header("2", "", &typed_columns) + row + "\tmarkers\n", false),
            // version 3 without one
            (
                header("3", "#reference_type\tmarkers\n", COLUMNS_GENOMES) + row + "\n",
                false,
            ),
            // declared type differs from the rows
            (
                header("3", "#reference_type\tmixed\n", &typed_columns) + row + "\tmarkers\n",
                false,
            ),
            (
                header("3", "#reference_type\tmarkers\n", &typed_columns) + row + "\tmarkers\n",
                true,
            ),
            (header("2", "", COLUMNS_GENOMES) + row + "\n", true),
        ];
        for (text, ok) in cases {
            std::fs::write(&path, &text).unwrap();
            assert_eq!(load_db_tsv(&path).is_ok(), ok, "{text}");
        }
    }
}

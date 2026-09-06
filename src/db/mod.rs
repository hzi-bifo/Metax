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
    index::{build_index_file_with_progress, BuildProgress, DEFAULT_INDEX_SEED_CHUNK_MIB},
    progress::log_message,
    seeds::{SeedMode, StrobeStartDistance},
};

pub const DB_TSV_FILE: &str = "metax_db.tsv.gz";
const METAM_INDEX_FILE: &str = "metax_db.mxi";
const MAPPING_REFERENCE_FILE: &str = "metax_reference.fasta";
const DB_VERSION: &str = "2";
pub const DEFAULT_SEED_CHUNK_MIB: usize = DEFAULT_INDEX_SEED_CHUNK_MIB;

#[derive(Clone, Debug)]
pub struct DatabaseMetadata {
    pub metam_index: PathBuf,
    pub references: Vec<ReferenceMetadata>,
    pub is_fractional: bool,
}

impl DatabaseMetadata {
    /// Stable content fingerprint of the reference table. Reused mappings
    /// (`--reuse-mapping` / `--resume`) record this so profiling can verify
    /// they were produced against this exact database. It hashes the parsed
    /// reference identities in load order — not the compressed file bytes — so
    /// recompressing `metax_db.tsv.gz` never changes it, while any real
    /// reference, length, or taxonomy difference does.
    pub fn fingerprint(&self) -> u128 {
        use std::fmt::Write as _;
        use xxhash_rust::xxh3::Xxh3;

        let mut hasher = Xxh3::new();
        hasher.update(if self.is_fractional {
            b"fractional\n"
        } else {
            b"absolute\n"
        });
        let mut row = String::new();
        for reference in &self.references {
            row.clear();
            // Tab-delimited canonical identity, mirroring metax_db.tsv.gz columns.
            let _ = writeln!(
                row,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                reference.name,
                reference.length,
                reference.full_genome_size,
                reference.effective_indexed_size,
                reference.genome,
                reference.contig,
                reference.taxonomy_field,
            );
            hasher.update(row.as_bytes());
        }
        hasher.digest128()
    }
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

    let genomes = read_reference_fasta(fasta_path)?;
    let selected = select_index_records(&genomes, subsample, threads)?;
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
    let references = if force_sampled_field {
        let reference_path = staging_path.join(MAPPING_REFERENCE_FILE);
        let _cleanup_reference = RemoveFileOnDrop::new(reference_path.clone());
        let references =
            write_metam_reference(&reference_path, fasta_path, selected, force_sampled_field)?;
        run_metam_index(
            &reference_path,
            &index_path,
            threads.unwrap_or(1).max(1),
            options,
        )?;
        gzip_file(&reference_path)
            .with_context(|| format!("failed to compress {}", reference_path.display()))?;
        references
    } else {
        let references = build_direct_reference_metadata(selected)?;
        run_metam_index(
            fasta_path,
            &index_path,
            threads.unwrap_or(1).max(1),
            options,
        )?;
        references
    };

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
    let mut version_ok = false;
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
                if *version != DB_VERSION {
                    bail!(
                        "unsupported Metax DB version {} in {}",
                        version,
                        path.display()
                    );
                }
                version_ok = true;
            }
            ["#mapping_index", value] => {
                metam_index = Some(resolve_relative_path(base, value));
            }
            ["#fractional", value] => {
                is_fractional = parse_bool(value)
                    .with_context(|| format!("invalid fractional flag on line {}", line_no + 1))?;
            }
            ["#columns", "reference_name", "reference_length", "full_genome_size", "effective_indexed_size", "genome", "contig", "taxonomy"] =>
                {}
            ["#columns", ..] => bail!(
                "unsupported database metadata columns on line {} in {}: {}",
                line_no + 1,
                path.display(),
                line
            ),
            [name, length, full, effective, genome, contig, taxonomy] => {
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

    if !version_ok {
        bail!("{} is not a Metax database metadata file", path.display());
    }
    let metam_index =
        metam_index.ok_or_else(|| anyhow!("{} is missing mapping_index", path.display()))?;
    if references.is_empty() {
        bail!("{} contains no reference metadata", path.display());
    }

    Ok(DatabaseMetadata {
        metam_index,
        references,
        is_fractional,
    })
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

fn write_db_tsv_body(writer: &mut dyn Write, metadata: &DatabaseMetadata) -> Result<()> {
    writeln!(writer, "#metax_db_version\t{DB_VERSION}")?;
    writeln!(
        writer,
        "#mapping_index\t{}",
        metadata.metam_index.to_string_lossy()
    )?;
    writeln!(writer, "#fractional\t{}", metadata.is_fractional)?;
    writeln!(
        writer,
        "#columns\treference_name\treference_length\tfull_genome_size\teffective_indexed_size\tgenome\tcontig\ttaxonomy"
    )?;
    for r in &metadata.references {
        ensure_tsv_safe(&r.name, "reference name")?;
        ensure_tsv_safe(&r.genome, "genome accession")?;
        ensure_tsv_safe(&r.contig, "contig id")?;
        ensure_tsv_safe(&r.taxonomy_field, "taxonomy field")?;
        writeln!(
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
    full_genome_size: usize,
    effective_indexed_size: usize,
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

fn read_reference_fasta(path: &Path) -> Result<Vec<GenomeInput>> {
    let reader = open_text_reader(path)?;
    let reader = BufReader::new(reader);
    let mut genomes: IndexMap<String, GenomeInput> = IndexMap::new();
    let mut current_header: Option<String> = None;
    let mut current_len = 0usize;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("failed reading FASTA {}", path.display()))?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(header) = current_header.take() {
                store_input_record(&mut genomes, &header, current_len)?;
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
        store_input_record(&mut genomes, &header, current_len)?;
    }
    if genomes.is_empty() {
        bail!("reference FASTA contains no records");
    }

    for genome in genomes.values_mut() {
        let sequence_sum: usize = genome.contigs.iter().map(|c| c.length).sum();
        if sequence_sum == 0 {
            bail!("genome {} has no sequence bases", genome.genome);
        }
        genome.full_genome_size = match genome.provided_genome_size {
            Some(size) => {
                if size < sequence_sum {
                    bail!(
                        "provided genome size {} for {} is smaller than summed contig length {}",
                        size,
                        genome.genome,
                        sequence_sum
                    );
                }
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

    Ok(genomes.into_values().collect())
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
) -> Result<()> {
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

    let entry = genomes
        .entry(parsed.genome.clone())
        .or_insert_with(|| GenomeInput {
            genome: parsed.genome.clone(),
            taxonomy: parsed.taxonomy.clone(),
            provided_genome_size: parsed.genome_size,
            contigs: Vec::new(),
            full_genome_size: 0,
        });
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
    Ok(())
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
    genomes: &[GenomeInput],
    subsample: Option<SubsampleParams>,
    threads: Option<usize>,
) -> Result<Vec<SelectedRecord>> {
    if genomes.is_empty() {
        return Ok(Vec::new());
    }
    let params = match subsample {
        Some(params) => params,
        None => {
            return Ok(genomes
                .iter()
                .flat_map(|genome| {
                    let effective: usize = genome.contigs.iter().map(|c| c.length).sum();
                    genome.contigs.iter().map(move |contig| SelectedRecord {
                        genome: genome.genome.clone(),
                        contig: contig.contig.clone(),
                        source_contig: contig.contig.clone(),
                        source_name: contig.reference_name.clone(),
                        taxonomy: genome.taxonomy.clone(),
                        full_genome_size: genome.full_genome_size,
                        effective_indexed_size: effective,
                        source_start: 1,
                        source_end: contig.length,
                    })
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
    if genome.full_genome_size < SHORT_GENOME_THRESHOLD {
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
            full_genome_size: genome.full_genome_size,
            effective_indexed_size,
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
    let mut references = Vec::with_capacity(records.len());
    let mut selected_by_source: HashMap<(String, String), Vec<PreparedRecord>> = HashMap::new();
    for record in records {
        let name = metax_reference_name(&record, force_sampled_field);
        if !seen_names.insert(name.clone()) {
            bail!("generated duplicate reference name '{}'", name);
        }
        references.push(ReferenceMetadata {
            name: name.clone(),
            length: record.len(),
            genome: record.genome.clone(),
            contig: record.contig.clone(),
            full_genome_size: record.full_genome_size,
            effective_indexed_size: record.effective_indexed_size,
            taxonomy_field: record.taxonomy.raw.clone(),
        });
        selected_by_source
            .entry((record.genome.clone(), record.source_contig.clone()))
            .or_default()
            .push(PreparedRecord { name, record });
    }

    stream_selected_reference_records(source_fasta, &mut selected_by_source, &mut writer)?;
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
        current_key.take(),
        current_header.take(),
        &sequence,
    )?;
    Ok(())
}

fn flush_selected_source_record(
    selected_by_source: &mut HashMap<(String, String), Vec<PreparedRecord>>,
    writer: &mut dyn Write,
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
) -> Result<()> {
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
            full_genome_size: 100,
            effective_indexed_size: 10,
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
            full_genome_size: 100,
            effective_indexed_size: 60,
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
}

use std::collections::HashSet;
use std::fs::File;
use std::io::{copy, Read};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::db::{database_directory, DatabaseMetadata, ReferenceType};
use crate::taxonomy::{TaxonomyRecord, DATABASE_TAXDMP_DIR, TAXONOMY_RECORD_FILE};

#[derive(Debug)]
struct IndexInfo {
    version: u32,
    metam_version: String,
    k: usize,
    w: usize,
    seed_mode: &'static str,
    max_occ: usize,
    hash_bits: u8,
    match_both_strobes: bool,
    strobe_start_dist_min: String,
    strobe_start_dist_max: String,
    created_unix_seconds: u64,
    reference_count: usize,
    total_reference_bases: u64,
    indexed_seed_hashes: u64,
    indexed_seed_occurrences: u64,
    filtered_seed_hashes: u64,
    reference_storage_bytes: u64,
    file_size_bytes: u64,
}

pub fn run_info(db: Option<&str>, index: Option<&str>) -> Result<()> {
    match (db, index) {
        (Some(db), None) => run_db_info(Path::new(db)),
        (None, Some(index)) => run_index_info(Path::new(index)),
        (Some(_), Some(_)) => bail!("use either --db or --index, not both"),
        (None, None) => bail!("please provide --db or --index"),
    }
}

fn run_db_info(path: &Path) -> Result<()> {
    let (metadata_path, metadata) = DatabaseMetadata::load_resolved(path)?;
    print_kv("info_type", "metax_database");
    print_kv("metadata_path", metadata_path.display());
    print_kv("mapping_index", metadata.metam_index.display());
    print_kv("fractional", metadata.is_fractional);
    print_metadata_summary(&metadata);
    print_taxonomy_summary(&database_directory(path)?)?;
    if metadata.metam_index.exists() {
        print_index_summary(&metadata.metam_index)?;
    } else {
        print_kv("index_present", false);
        eprintln!(
            "warning: mapping index {} does not exist; showing metadata only",
            metadata.metam_index.display()
        );
    }
    Ok(())
}

fn run_index_info(path: &Path) -> Result<()> {
    print_kv("info_type", "mapping_index");
    print_kv("index_path", path.display());
    print_index_summary(path)
}

fn print_metadata_summary(metadata: &DatabaseMetadata) {
    let mut genomes = HashSet::new();
    let mut marker_genomes = 0usize;
    let mut sized_marker_genomes = 0usize;
    let mut reference_bases = 0u64;
    let mut full_genome_bases = 0u64;
    let mut effective_indexed_bases = 0u64;
    for reference in &metadata.references {
        reference_bases += reference.length as u64;
        // Genome sizes repeat on every row of a genome; count each genome once.
        if genomes.insert(&*reference.genome) {
            full_genome_bases += reference.full_genome_size as u64;
            effective_indexed_bases += reference.effective_indexed_size as u64;
            marker_genomes += usize::from(reference.reference_type == ReferenceType::Markers);
            sized_marker_genomes += usize::from(reference.declared_genome_size.is_some());
        }
    }
    print_kv("reference_type", metadata.reference_type_label());
    print_kv("metadata_reference_count", metadata.references.len());
    print_kv("metadata_genome_count", genomes.len());
    print_kv("metadata_marker_genome_count", marker_genomes);
    if marker_genomes > 0 {
        print_kv("metadata_marker_genomes_with_size", sized_marker_genomes);
    }
    print_kv("metadata_reference_bases", reference_bases);
    print_kv("metadata_full_genome_bases", full_genome_bases);
    print_kv("metadata_effective_indexed_bases", effective_indexed_bases);
}

/// Whether the database stores its taxonomy dump, and how many taxa its
/// taxonomy record holds.
fn print_taxonomy_summary(database: &Path) -> Result<()> {
    print_kv(
        "taxonomy_dump_stored",
        database.join(DATABASE_TAXDMP_DIR).is_dir(),
    );
    let record = database.join(TAXONOMY_RECORD_FILE);
    if record.exists() {
        print_kv("taxonomy_record_taxa", TaxonomyRecord::load(&record)?.len());
    } else {
        print_kv("taxonomy_record_taxa", "NA");
    }
    Ok(())
}

fn print_index_summary(path: &Path) -> Result<()> {
    let info = read_index_info(path)?;
    print_kv("index_present", true);
    print_kv("index_path", path.display());
    print_kv("index_version", info.version);
    print_kv("metam_version", info.metam_version);
    print_kv("k", info.k);
    print_kv("w", info.w);
    print_kv("seed_mode", info.seed_mode);
    print_kv("max_occ", info.max_occ);
    print_kv("hash_bits", info.hash_bits);
    print_kv("match_both_strobes", info.match_both_strobes);
    print_kv("strobe_start_dist_min", info.strobe_start_dist_min);
    print_kv("strobe_start_dist_max", info.strobe_start_dist_max);
    print_kv("created_unix_seconds", info.created_unix_seconds);
    print_kv("reference_count", info.reference_count);
    print_kv("total_reference_bases", info.total_reference_bases);
    print_kv("indexed_seed_hashes", info.indexed_seed_hashes);
    print_kv("indexed_seed_occurrences", info.indexed_seed_occurrences);
    print_kv("filtered_seed_hashes", info.filtered_seed_hashes);
    print_kv("reference_storage_bytes", info.reference_storage_bytes);
    print_kv("index_file_size_bytes", info.file_size_bytes);
    print_kv("index_file_size_human", human_bytes(info.file_size_bytes));
    if info.indexed_seed_hashes > 0 {
        let avg = info.indexed_seed_occurrences as f64 / info.indexed_seed_hashes as f64;
        print_kv("avg_occurrences_per_key", format!("{avg:.3}"));
    }
    Ok(())
}

fn read_index_info(path: &Path) -> Result<IndexInfo> {
    let file =
        File::open(path).with_context(|| format!("failed to open index {}", path.display()))?;
    let file_size_bytes = file.metadata()?.len();
    let mut reader = zstd::stream::Decoder::new(file)
        .with_context(|| format!("failed to initialize zstd decoder for {}", path.display()))?;

    let mut magic = [0u8; 8];
    reader.read_exact(&mut magic)?;
    let version = read_varint(&mut reader)? as u32;
    metam::index::validate_index_identity(magic, version)
        .with_context(|| format!("cannot read index {}", path.display()))?;
    let metam_version = read_string(&mut reader)?;
    let k = read_varint(&mut reader)? as usize;
    let w = read_varint(&mut reader)? as usize;
    let seed_mode_code = read_varint(&mut reader)?;
    let seed_mode = match seed_mode_code {
        0 => "minimizer",
        1 => "randstrobe",
        _ => bail!(
            "invalid seed mode code {seed_mode_code} in {}",
            path.display()
        ),
    };
    let max_occ = read_varint(&mut reader)? as usize;
    let mut hash_bits = [0u8; 1];
    reader.read_exact(&mut hash_bits)?;
    let hash_bits = hash_bits[0];
    let match_both_strobes = if version >= 5 {
        let mut value = [0u8; 1];
        reader.read_exact(&mut value)?;
        match value[0] {
            0 => false,
            1 => true,
            _ => bail!("invalid match-both-strobes flag in {}", path.display()),
        }
    } else {
        false
    };
    let (strobe_start_dist_min, strobe_start_dist_max) = if version >= 6 {
        let min = read_varint(&mut reader)?;
        let max = read_varint(&mut reader)?;
        if min == 0 && max == 0 {
            ("NA".to_string(), "NA".to_string())
        } else {
            (min.to_string(), max.to_string())
        }
    } else if seed_mode == "randstrobe" {
        ("legacy".to_string(), "legacy".to_string())
    } else {
        ("NA".to_string(), "NA".to_string())
    };
    let created_unix_seconds = read_varint(&mut reader)?;
    let filtered_seed_hashes = read_varint(&mut reader)?;
    let reference_count = read_varint(&mut reader)? as usize;

    let mut total_reference_bases = 0u64;
    let mut reference_storage_bytes = 0u64;
    for _ in 0..reference_count {
        skip_len_prefixed_bytes(&mut reader)?;
        let len = read_varint(&mut reader)?;
        let packed_len = skip_len_prefixed_bytes(&mut reader)?;
        let ambiguous_len = skip_len_prefixed_bytes(&mut reader)?;
        total_reference_bases = total_reference_bases.saturating_add(len);
        reference_storage_bytes =
            reference_storage_bytes.saturating_add(packed_len.saturating_add(ambiguous_len));
    }
    let indexed_seed_occurrences = read_varint(&mut reader)?;
    let indexed_seed_hashes = read_varint(&mut reader)?;

    Ok(IndexInfo {
        version,
        metam_version,
        k,
        w,
        seed_mode,
        max_occ,
        hash_bits,
        match_both_strobes,
        strobe_start_dist_min,
        strobe_start_dist_max,
        created_unix_seconds,
        reference_count,
        total_reference_bases,
        indexed_seed_hashes,
        indexed_seed_occurrences,
        filtered_seed_hashes,
        reference_storage_bytes,
        file_size_bytes,
    })
}

fn read_string<R: Read>(reader: &mut R) -> Result<String> {
    let len = read_varint(reader)?;
    let len = usize::try_from(len).context("string length does not fit in usize")?;
    let mut bytes = vec![0u8; len];
    reader.read_exact(&mut bytes)?;
    String::from_utf8(bytes).context("index string is not valid UTF-8")
}

fn skip_len_prefixed_bytes<R: Read>(reader: &mut R) -> Result<u64> {
    let len = read_varint(reader)?;
    let mut limited = reader.take(len);
    let skipped = copy(&mut limited, &mut std::io::sink())?;
    if skipped != len {
        bail!("truncated length-prefixed byte field in index");
    }
    Ok(len)
}

fn read_varint<R: Read>(reader: &mut R) -> Result<u64> {
    let mut value = 0u64;
    let mut shift = 0u32;
    loop {
        let mut byte = [0u8; 1];
        reader.read_exact(&mut byte)?;
        value |= u64::from(byte[0] & 0x7f) << shift;
        if byte[0] & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
        if shift >= 64 {
            bail!("invalid varint in index");
        }
    }
}

fn human_bytes(bytes: u64) -> String {
    let units = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0usize;
    while value >= 1024.0 && unit + 1 < units.len() {
        value /= 1024.0;
        unit += 1;
    }
    format!("{value:.2} {}", units[unit])
}

fn print_kv(key: &str, value: impl std::fmt::Display) {
    println!("{key}\t{value}");
}

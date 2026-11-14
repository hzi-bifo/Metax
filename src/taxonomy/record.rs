//! Record of a database's taxa in the taxonomy dump it was built with.

use super::Taxonomy;
use anyhow::{anyhow, bail, Context, Result};
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use xxhash_rust::xxh3::xxh3_64;

/// File of a database directory that holds its taxonomy record.
pub const TAXONOMY_RECORD_FILE: &str = "metax_taxonomy.tsv.gz";
const RECORD_VERSION: &str = "1";
const COLUMNS: &str = "taxid\tlineage_hash";
/// Taxids named in a message before it is cut short.
const LISTED_TAXIDS: usize = 10;

/// A hash of the lineage of every taxid the references of a database name,
/// taken in the dump the database was built with. Another dump gives the
/// same output for these taxa exactly when it gives the same hashes.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TaxonomyRecord {
    /// (taxid, lineage hash), sorted by taxid.
    lineages: Vec<(u32, u64)>,
    provenance: Option<DumpProvenance>,
}

/// Where a record's dump came from: the dump files stored with the database
/// and the Metax version and time that recorded them.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DumpProvenance {
    pub metax_version: String,
    /// UTC time, as `YYYY-MM-DDThh:mm:ssZ`.
    pub created: String,
    /// (file name, bytes, SHA-256) of each stored dump file.
    pub files: Vec<(String, u64, String)>,
}

impl DumpProvenance {
    /// Provenance of `files` recorded now by this Metax version.
    pub fn now(files: Vec<(String, u64, String)>) -> Self {
        Self {
            metax_version: env!("CARGO_PKG_VERSION").to_string(),
            created: chrono::Utc::now().format("%Y-%m-%dT%H:%M:%SZ").to_string(),
            files,
        }
    }
}

impl TaxonomyRecord {
    /// Record `taxids` in `taxonomy`. Fails when the dump lacks any of them.
    pub fn build(taxonomy: &Taxonomy, taxids: &BTreeSet<u32>) -> Result<Self> {
        let mut lineages = Vec::with_capacity(taxids.len());
        let mut missing = Vec::new();
        for &taxid in taxids {
            match lineage_hash(taxonomy, taxid) {
                Some(hash) => lineages.push((taxid, hash)),
                None => missing.push(taxid),
            }
        }
        if !missing.is_empty() {
            bail!(
                "{} of {} reference taxids are not in the taxonomy dump ({})",
                missing.len(),
                taxids.len(),
                list_taxids(&missing)
            );
        }
        Ok(Self {
            lineages,
            provenance: None,
        })
    }

    pub fn with_provenance(self, provenance: Option<DumpProvenance>) -> Self {
        Self { provenance, ..self }
    }

    pub fn provenance(&self) -> Option<&DumpProvenance> {
        self.provenance.as_ref()
    }

    /// Number of recorded taxa.
    pub fn len(&self) -> usize {
        self.lineages.len()
    }

    pub fn is_empty(&self) -> bool {
        self.lineages.is_empty()
    }

    pub fn write(&self, path: &Path) -> Result<()> {
        let file =
            File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut writer = GzEncoder::new(BufWriter::new(file), Compression::default());
        writeln!(writer, "#metax_taxonomy_version\t{RECORD_VERSION}")?;
        if let Some(provenance) = &self.provenance {
            writeln!(writer, "#metax_version\t{}", provenance.metax_version)?;
            writeln!(writer, "#created\t{}", provenance.created)?;
            for (name, bytes, sha256) in &provenance.files {
                writeln!(writer, "#dump_file\t{name}\t{bytes}\t{sha256}")?;
            }
        }
        writeln!(writer, "#columns\t{COLUMNS}")?;
        for (taxid, hash) in &self.lineages {
            writeln!(writer, "{taxid}\t{hash:016x}")?;
        }
        writer
            .finish()
            .and_then(|mut inner| inner.flush())
            .with_context(|| format!("failed to finish {}", path.display()))
    }

    /// The record of the database in `database`, or `None` for a database
    /// built without one.
    pub fn load_from_database(database: &Path) -> Result<Option<Self>> {
        let path = database.join(TAXONOMY_RECORD_FILE);
        if !path.exists() {
            return Ok(None);
        }
        Self::load(&path).map(Some)
    }

    /// The recorded taxa whose lineage differs in `taxonomy`, split into
    /// those it merges into another taxid and the others, and the recorded
    /// taxa it lacks.
    fn compare(&self, taxonomy: &Taxonomy) -> Comparison {
        let mut comparison = Comparison::default();
        for &(taxid, hash) in &self.lineages {
            match lineage_hash(taxonomy, taxid) {
                Some(current) if current == hash => {}
                Some(_) if taxonomy.resolve_current_taxid(taxid) != Some(taxid) => {
                    comparison.merged.push(taxid)
                }
                Some(_) => comparison.changed.push(taxid),
                None => comparison.missing.push(taxid),
            }
        }
        comparison
    }

    /// Compare `taxonomy`, loaded from `dump`, with the record. Merged or
    /// otherwise changed taxa are reported as a warning; taxa the dump lacks
    /// stop the run unless `allow_missing`.
    pub fn check(&self, taxonomy: &Taxonomy, dump: &Path, allow_missing: bool) -> Result<()> {
        let Comparison {
            merged,
            changed,
            missing,
        } = self.compare(taxonomy);
        if merged.is_empty() && changed.is_empty() && missing.is_empty() {
            log::info!(
                target: "TAXONOMY",
                "The taxonomy dump in {} matches the database ({} taxa{}).",
                dump.display(),
                self.len(),
                self.provenance.as_ref().map_or(String::new(), |provenance| format!(
                    "; recorded by Metax {} on {}",
                    provenance.metax_version, provenance.created
                ))
            );
            return Ok(());
        }
        let mut differences = Vec::new();
        for (taxids, what) in [
            (&merged, "are merged into other taxids in it"),
            (&changed, "have another name, rank or lineage"),
            (&missing, "are not in it"),
        ] {
            if !taxids.is_empty() {
                differences.push(format!("{} {what} ({})", taxids.len(), list_taxids(taxids)));
            }
        }
        let detail = format!(
            "taxonomy dump in {} differs from the one the database was built with: of the {} database taxa, {}",
            dump.display(),
            self.len(),
            differences.join("; ")
        );
        if !missing.is_empty() && !allow_missing {
            bail!(
                "the {detail}. Use the dump stored with the database (used when --dmp-dir is omitted), or pass --skip-missing-taxa to continue without the references of these taxa"
            );
        }
        log::warn!(
            target: "TAXONOMY",
            "The {detail}. The profile uses the names and lineages of this dump."
        );
        Ok(())
    }

    pub fn load(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        let reader = BufReader::new(MultiGzDecoder::new(file));
        let mut version = None;
        let mut lineages: Vec<(u32, u64)> = Vec::new();
        let mut metax_version = None;
        let mut created = None;
        let mut files = Vec::new();
        for (index, line) in reader.lines().enumerate() {
            let line = line.with_context(|| format!("failed to read {}", path.display()))?;
            let parsed = match line.split_once('\t') {
                Some(("#metax_taxonomy_version", value)) => {
                    if value != RECORD_VERSION {
                        bail!("unsupported taxonomy record version {value}");
                    }
                    version = Some(value.to_string());
                    Ok(())
                }
                Some(("#metax_version", value)) if version.is_some() => {
                    metax_version = Some(value.to_string());
                    Ok(())
                }
                Some(("#created", value)) if version.is_some() => {
                    created = Some(value.to_string());
                    Ok(())
                }
                Some(("#dump_file", value)) if version.is_some() => {
                    match value.split('\t').collect::<Vec<_>>()[..] {
                        [name, bytes, sha256] if bytes.parse::<u64>().is_ok() => {
                            files.push((
                                name.to_string(),
                                bytes.parse().expect("checked"),
                                sha256.to_string(),
                            ));
                            Ok(())
                        }
                        _ => Err(anyhow!("malformed dump file line")),
                    }
                }
                Some(("#columns", value)) if value == COLUMNS => Ok(()),
                Some((taxid, hash)) if version.is_some() && !taxid.starts_with('#') => {
                    let taxid = taxid.parse::<u32>().ok();
                    let hash = u64::from_str_radix(hash, 16).ok();
                    match (taxid, hash) {
                        (Some(taxid), Some(hash))
                            if lineages.last().is_none_or(|(last, _)| *last < taxid) =>
                        {
                            lineages.push((taxid, hash));
                            Ok(())
                        }
                        _ => Err(anyhow!("malformed or unsorted row")),
                    }
                }
                _ => Err(anyhow!("unexpected line")),
            };
            parsed.with_context(|| format!("{} line {}: {}", path.display(), index + 1, line))?;
        }
        if version.is_none() {
            bail!("{} is not a Metax taxonomy record", path.display());
        }
        let provenance = match (metax_version, created) {
            (Some(metax_version), Some(created)) => Some(DumpProvenance {
                metax_version,
                created,
                files,
            }),
            (None, None) if files.is_empty() => None,
            _ => bail!("{} has incomplete provenance lines", path.display()),
        };
        Ok(Self {
            lineages,
            provenance,
        })
    }
}

/// Recorded taxa that another dump treats differently.
#[derive(Debug, Default)]
struct Comparison {
    merged: Vec<u32>,
    changed: Vec<u32>,
    missing: Vec<u32>,
}

/// Hash of the scientific name and rank of `taxid` and of each of its
/// ancestors, or `None` when the dump lacks `taxid`. A merged taxid hashes
/// the lineage of the taxid it was merged into, which is what the profile
/// reports for it. Changing this layout invalidates stored records.
fn lineage_hash(taxonomy: &Taxonomy, taxid: u32) -> Option<u64> {
    let lineage = taxonomy.get_parents(taxid);
    if lineage.is_empty() {
        return None;
    }
    let mut text = Vec::new();
    for (node, name, rank) in &lineage {
        writeln!(text, "{node}\t{rank}\t{name}").expect("writing to a vector cannot fail");
    }
    Some(xxh3_64(&text))
}

fn list_taxids(taxids: &[u32]) -> String {
    let mut listed = taxids
        .iter()
        .take(LISTED_TAXIDS)
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    if taxids.len() > LISTED_TAXIDS {
        listed.push_str(", ...");
    }
    format!("taxid {listed}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashMap, HashSet};

    fn taxonomy(species_name: &str) -> Taxonomy {
        Taxonomy::from_parts(
            HashMap::from([
                (1, (1, "no rank".to_string())),
                (10, (1, "genus".to_string())),
                (20, (10, "species".to_string())),
                (30, (10, "species".to_string())),
            ]),
            HashMap::from([
                (1, "root".to_string()),
                (10, "Genus ten".to_string()),
                (20, species_name.to_string()),
                (30, "Species thirty".to_string()),
            ]),
            HashMap::from([(25, 20)]),
            HashSet::new(),
        )
    }

    #[test]
    fn record_round_trips_and_follows_lineage_content() {
        let taxids = BTreeSet::from([20, 25, 30]);
        let record = TaxonomyRecord::build(&taxonomy("Species twenty"), &taxids).unwrap();
        assert_eq!(record.len(), 3);
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join(TAXONOMY_RECORD_FILE);
        record.write(&path).unwrap();
        assert_eq!(TaxonomyRecord::load(&path).unwrap(), record);
        let with_provenance = record
            .clone()
            .with_provenance(Some(DumpProvenance::now(vec![(
                "nodes.dmp".to_string(),
                12,
                "ab".repeat(32),
            )])));
        with_provenance.write(&path).unwrap();
        assert_eq!(TaxonomyRecord::load(&path).unwrap(), with_provenance);

        let renamed = TaxonomyRecord::build(&taxonomy("Species renamed"), &taxids).unwrap();
        assert_ne!(renamed, record);
        assert_eq!(renamed.lineages[2], record.lineages[2]);
        // The merged taxid hashes the lineage of the taxid it was merged into.
        assert_eq!(record.lineages[0].1, record.lineages[1].1);
    }

    #[test]
    fn changed_taxa_warn_and_missing_taxa_stop_unless_allowed() {
        let taxids = BTreeSet::from([20, 25, 30]);
        let record = TaxonomyRecord::build(&taxonomy("Species twenty"), &taxids).unwrap();
        let dump = Path::new("dump");
        record
            .check(&taxonomy("Species twenty"), dump, false)
            .unwrap();

        // Renaming species 20 changes taxid 20 and taxid 25, which the dump
        // merges into 20.
        let renamed = taxonomy("Species renamed");
        let comparison = record.compare(&renamed);
        assert_eq!(comparison.merged, [25]);
        assert_eq!(comparison.changed, [20]);
        assert!(comparison.missing.is_empty());
        record.check(&renamed, dump, false).unwrap();

        let mut without_thirty = taxonomy("Species twenty");
        without_thirty.nodes.remove(&30);
        without_thirty.names.remove(&30);
        assert_eq!(record.compare(&without_thirty).missing, [30]);
        let message = record
            .check(&without_thirty, dump, false)
            .unwrap_err()
            .to_string();
        assert!(message.contains("of the 3 database taxa, 1 are not in it (taxid 30)"));
        assert!(message.contains("--skip-missing-taxa"));
        record.check(&without_thirty, dump, true).unwrap();
    }

    #[test]
    fn taxids_missing_from_the_dump_are_an_error() {
        let error =
            TaxonomyRecord::build(&taxonomy("Species twenty"), &BTreeSet::from([20, 40, 50]))
                .unwrap_err();
        let message = error.to_string();
        assert!(message.contains("2 of 3 reference taxids"));
        assert!(message.contains("taxid 40, 50"));
    }

    #[test]
    fn malformed_records_are_rejected() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join(TAXONOMY_RECORD_FILE);
        for body in [
            "20\t00000000000000ff\n",
            "#metax_taxonomy_version\t2\n",
            "#metax_taxonomy_version\t1\n#columns\ttaxid\tlineage_hash\n30\t1\n20\t2\n",
            "#metax_taxonomy_version\t1\n20\tnot-hex\n",
        ] {
            let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::fast());
            encoder.write_all(body.as_bytes()).unwrap();
            encoder.finish().unwrap();
            assert!(TaxonomyRecord::load(&path).is_err(), "{body:?}");
        }
    }
}

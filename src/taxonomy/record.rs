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
        Ok(Self { lineages })
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
        writeln!(writer, "#columns\t{COLUMNS}")?;
        for (taxid, hash) in &self.lineages {
            writeln!(writer, "{taxid}\t{hash:016x}")?;
        }
        writer
            .finish()
            .and_then(|mut inner| inner.flush())
            .with_context(|| format!("failed to finish {}", path.display()))
    }

    pub fn load(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        let reader = BufReader::new(MultiGzDecoder::new(file));
        let mut version = None;
        let mut lineages: Vec<(u32, u64)> = Vec::new();
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
        Ok(Self { lineages })
    }
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

        let renamed = TaxonomyRecord::build(&taxonomy("Species renamed"), &taxids).unwrap();
        assert_ne!(renamed, record);
        assert_eq!(renamed.lineages[2], record.lineages[2]);
        // The merged taxid hashes the lineage of the taxid it was merged into.
        assert_eq!(record.lineages[0].1, record.lineages[1].1);
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

mod dump;
mod record;

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

pub use dump::TaxonomySource;
pub(crate) use dump::{store_dump, DATABASE_TAXDMP_ARCHIVE, DATABASE_TAXDMP_DIR};
pub use record::{TaxonomyRecord, TAXONOMY_RECORD_FILE};

#[derive(Clone, Debug)]
pub struct Taxonomy {
    // tax_id -> (parent_id, rank)
    nodes: HashMap<u32, (u32, String)>,
    // tax_id -> scientific name
    names: HashMap<u32, String>,
    // deprecated tax_id -> current tax_id
    merged: HashMap<u32, u32>,
    // deleted tax_ids to suppress
    del_nodes: HashSet<u32>,
}

impl Taxonomy {
    /// Locate the taxonomy dump: `dmp_dir` when given, else the dump stored
    /// in the `database` directory, else the managed NCBI dump.
    pub(crate) fn prepare_dmp(
        dmp_dir: Option<&str>,
        database: Option<&Path>,
    ) -> Result<TaxonomySource> {
        TaxonomySource::prepare(dmp_dir, database)
    }

    pub(super) fn from_dir(dir: &Path) -> Result<Self> {
        let nodes = dump_file_path(dir, "nodes.dmp");
        let names = dump_file_path(dir, "names.dmp");
        let merged = dump_file_path(dir, "merged.dmp");
        let delnodes = dump_file_path(dir, "delnodes.dmp");
        Self::from_files(&nodes, &names, &merged, Some(&delnodes))
    }

    /// Load the dump files from a `.tar.gz` archive without unpacking it.
    /// Entries are matched by file name in any folder of the archive, and
    /// may themselves be gzipped; other entries are skipped.
    pub(crate) fn from_archive(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        let mut archive = tar::Archive::new(MultiGzDecoder::new(BufReader::new(file)));
        let mut nodes = None;
        let mut names = None;
        let mut merged = None;
        let mut del_nodes = None;
        for entry in archive
            .entries()
            .with_context(|| format!("failed to read {}", path.display()))?
        {
            let entry = entry.with_context(|| format!("failed to read {}", path.display()))?;
            if !entry.header().entry_type().is_file() {
                continue;
            }
            let entry_path = entry
                .path()
                .with_context(|| format!("invalid entry path in {}", path.display()))?
                .into_owned();
            let Some(file_name) = entry_path.file_name().and_then(|name| name.to_str()) else {
                continue;
            };
            let (name, gzipped) = match file_name.strip_suffix(".gz") {
                Some(name) => (name, true),
                None => (file_name, false),
            };
            if !dump::REQUIRED_DUMP_FILES.contains(&name) {
                continue;
            }
            let source = format!("{}:{}", path.display(), entry_path.display());
            let reader: Box<dyn BufRead> = if gzipped {
                Box::new(BufReader::new(MultiGzDecoder::new(entry)))
            } else {
                Box::new(BufReader::new(entry))
            };
            let duplicate = match name {
                "nodes.dmp" => nodes.replace(load_nodes(reader, &source)?).is_some(),
                "names.dmp" => names.replace(load_names(reader, &source)?).is_some(),
                "merged.dmp" => merged.replace(load_merged(reader, &source)?).is_some(),
                _ => del_nodes.replace(load_delnodes(reader, &source)?).is_some(),
            };
            if duplicate {
                bail!("{} holds {name} more than once", path.display());
            }
        }
        let missing = [
            ("nodes.dmp", nodes.is_none()),
            ("names.dmp", names.is_none()),
            ("merged.dmp", merged.is_none()),
            ("delnodes.dmp", del_nodes.is_none()),
        ]
        .into_iter()
        .filter_map(|(name, absent)| absent.then_some(name))
        .collect::<Vec<_>>();
        if !missing.is_empty() {
            bail!("{} lacks {}", path.display(), missing.join(", "));
        }
        let taxonomy = Self {
            nodes: nodes.expect("checked above"),
            names: names.expect("checked above"),
            merged: merged.expect("checked above"),
            del_nodes: del_nodes.expect("checked above"),
        };
        taxonomy.validate()?;
        Ok(taxonomy)
    }

    pub fn from_files(
        nodes: &Path,
        names: &Path,
        merged: &Path,
        delnodes: Option<&Path>,
    ) -> Result<Self> {
        // Load each taxonomy table into memory for fast lookups.
        let label = |path: &Path| path.display().to_string();
        let nodes_iter = load_nodes(open_dump_file(nodes)?, &label(nodes))?;
        let names_iter = load_names(open_dump_file(names)?, &label(names))?;
        let merged_iter = load_merged(open_dump_file(merged)?, &label(merged))?;
        let del_nodes = if let Some(path) = delnodes {
            load_delnodes(open_dump_file(path)?, &label(path))?
        } else {
            HashSet::new()
        };
        let taxonomy = Self {
            nodes: nodes_iter,
            names: names_iter,
            merged: merged_iter,
            del_nodes,
        };
        taxonomy.validate()?;
        Ok(taxonomy)
    }

    fn validate(&self) -> Result<()> {
        if self.nodes.is_empty() {
            bail!("nodes.dmp contains no taxonomy nodes");
        }
        if self.names.is_empty() {
            bail!("names.dmp contains no scientific names");
        }
        match self.nodes.get(&1) {
            Some((1, _)) => {}
            Some((parent, _)) => bail!("invalid taxonomy root: taxid 1 has parent {parent}"),
            None => bail!("nodes.dmp does not contain the NCBI taxonomy root (taxid 1)"),
        }
        if !self.names.contains_key(&1) {
            bail!("names.dmp does not contain a scientific name for taxid 1");
        }
        if self.names.len() != self.nodes.len() {
            bail!(
                "taxonomy node/name count mismatch (nodes={}, scientific_names={}); the dump may be incomplete",
                self.nodes.len(),
                self.names.len()
            );
        }
        Ok(())
    }

    /// Resolve a possibly deprecated taxid to the current taxonomy node.
    ///
    /// NCBI normally emits direct old-to-current entries in `merged.dmp`, but
    /// following the complete chain also makes mixed-age dump snapshots safe.
    /// Cyclic merge data, deleted terminal nodes, and unknown nodes are
    /// rejected instead of being allowed to masquerade as current taxids.
    pub fn resolve_current_taxid(&self, tax_id: u32) -> Option<u32> {
        // Almost every database taxid is already current. Keep that path to
        // one merged-table probe and two validation probes.
        if !self.merged.contains_key(&tax_id) {
            return (!self.del_nodes.contains(&tax_id) && self.nodes.contains_key(&tax_id))
                .then_some(tax_id);
        }

        // Detect a malformed merge cycle without allocating a visited set on
        // every taxonomy lookup.
        let mut slow = tax_id;
        let mut fast = tax_id;
        while let Some(next_slow) = self.merged.get(&slow).copied() {
            slow = next_slow;

            let Some(next_fast) = self.merged.get(&fast).copied() else {
                break;
            };
            let Some(next_fast) = self.merged.get(&next_fast).copied() else {
                break;
            };
            fast = next_fast;
            if slow == fast {
                return None;
            }
        }

        let mut current = tax_id;
        while let Some(next) = self.merged.get(&current).copied() {
            current = next;
        }
        if self.del_nodes.contains(&current) || !self.nodes.contains_key(&current) {
            None
        } else {
            Some(current)
        }
    }

    pub fn get_name(&self, tax_id: u32) -> Option<&str> {
        let tax_id = self.resolve_current_taxid(tax_id)?;
        self.names.get(&tax_id).map(|s| s.as_str())
    }

    pub fn get_rank(&self, tax_id: u32) -> Option<&str> {
        let tax_id = self.resolve_current_taxid(tax_id)?;
        self.nodes.get(&tax_id).map(|(_, rank)| rank.as_str())
    }

    pub fn get_parents(&self, tax_id: u32) -> Vec<(u32, String, String)> {
        // Walk up the taxonomy to root, returning the lineage.
        let mut parents = Vec::new();
        let Some(mut current) = self.resolve_current_taxid(tax_id) else {
            return parents;
        };
        let mut seen = HashSet::new();
        loop {
            if !seen.insert(current) {
                break;
            }
            let Some((parent_id, rank)) = self.nodes.get(&current) else {
                break;
            };
            let name = self.names.get(&current).cloned().unwrap_or_default();
            parents.push((current, name, rank.clone()));
            if *parent_id == current {
                break;
            }
            let Some(parent) = self.resolve_current_taxid(*parent_id) else {
                break;
            };
            current = parent;
        }
        parents
    }

    #[cfg(test)]
    pub(crate) fn from_parts(
        nodes: HashMap<u32, (u32, String)>,
        names: HashMap<u32, String>,
        merged: HashMap<u32, u32>,
        del_nodes: HashSet<u32>,
    ) -> Self {
        Self {
            nodes,
            names,
            merged,
            del_nodes,
        }
    }

    pub fn get_majority_lca(
        &self,
        tax_ids: &HashSet<u32>,
        cutoff_fraction: f64,
    ) -> Option<(u32, String, String)> {
        // Find the deepest node shared by at least cutoff_fraction of paths.
        if tax_ids.is_empty() {
            return None;
        }
        // Build root-to-leaf paths for each taxid, then scan for consensus.
        let paths: Vec<Vec<(u32, String, String)>> = tax_ids
            .iter()
            .map(|tax_id| {
                let mut parents = self.get_parents(*tax_id);
                parents.reverse();
                parents
            })
            .collect();
        let min_count = (cutoff_fraction * paths.len() as f64).ceil() as usize;
        let mut lca: Option<(u32, String, String)> = None;
        let mut level = 0;
        loop {
            let mut counter: IndexMap<u32, usize> = IndexMap::new();
            let mut exhausted = true;
            for path in &paths {
                if level < path.len() {
                    exhausted = false;
                    let (taxid, name, rank) = &path[level];
                    *counter.entry(*taxid).or_insert(0) += 1;
                    if counter[taxid] >= min_count {
                        lca = Some((*taxid, name.clone(), rank.clone()));
                    }
                }
            }
            if exhausted {
                break;
            }
            level += 1;
        }
        lca
    }
}

/// Path of the dump file `name` in `dir`: the plain file, or its gzipped
/// form `name.gz` when only that exists.
pub(crate) fn dump_file_path(dir: &Path, name: &str) -> PathBuf {
    let plain = dir.join(name);
    if !plain.exists() {
        let gzipped = dir.join(format!("{name}.gz"));
        if gzipped.exists() {
            return gzipped;
        }
    }
    plain
}

fn open_dump_file(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    if path.extension().is_some_and(|extension| extension == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn load_nodes(reader: impl BufRead, source: &str) -> Result<HashMap<u32, (u32, String)>> {
    // nodes.dmp: tax_id | parent_tax_id | rank
    let mut nodes = HashMap::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line = line.with_context(|| format!("failed to read {} line {line_number}", source))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('|').map(str::trim);
        let tax_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing taxid in {} line {line_number}", source))?
            .parse::<u32>()
            .with_context(|| format!("invalid taxid in {} line {line_number}", source))?;
        let parent_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing parent taxid in {} line {line_number}", source))?
            .parse::<u32>()
            .with_context(|| format!("invalid parent taxid in {} line {line_number}", source))?;
        let rank = fields
            .next()
            .with_context(|| format!("missing rank in {} line {line_number}", source))?
            .to_string();
        if nodes.insert(tax_id, (parent_id, rank)).is_some() {
            bail!("duplicate taxid {tax_id} in {} line {line_number}", source);
        }
    }
    Ok(nodes)
}

fn load_names(reader: impl BufRead, source: &str) -> Result<HashMap<u32, String>> {
    // names.dmp: tax_id | name | ... | name_class (keep scientific name only)
    let mut names = HashMap::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line = line.with_context(|| format!("failed to read {} line {line_number}", source))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('|').map(str::trim);
        let tax_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing taxid in {} line {line_number}", source))?
            .parse::<u32>()
            .with_context(|| format!("invalid taxid in {} line {line_number}", source))?;
        let name = fields
            .next()
            .with_context(|| format!("missing name in {} line {line_number}", source))?;
        let _unique_name = fields.next().with_context(|| {
            format!("missing unique-name field in {} line {line_number}", source)
        })?;
        let name_class = fields
            .next()
            .with_context(|| format!("missing name class in {} line {line_number}", source))?;
        if name_class == "scientific name" && names.insert(tax_id, name.to_string()).is_some() {
            bail!(
                "duplicate scientific name for taxid {tax_id} in {} line {line_number}",
                source
            );
        }
    }
    Ok(names)
}

fn load_merged(reader: impl BufRead, source: &str) -> Result<HashMap<u32, u32>> {
    // merged.dmp: old_tax_id | new_tax_id
    let mut merged = HashMap::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line = line.with_context(|| format!("failed to read {} line {line_number}", source))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('|').map(str::trim);
        let old_tax_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing old taxid in {} line {line_number}", source))?
            .parse::<u32>()
            .with_context(|| format!("invalid old taxid in {} line {line_number}", source))?;
        let new_tax_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing new taxid in {} line {line_number}", source))?
            .parse::<u32>()
            .with_context(|| format!("invalid new taxid in {} line {line_number}", source))?;
        if merged.insert(old_tax_id, new_tax_id).is_some() {
            bail!(
                "duplicate merged taxid {old_tax_id} in {} line {line_number}",
                source
            );
        }
    }
    Ok(merged)
}

fn load_delnodes(reader: impl BufRead, source: &str) -> Result<HashSet<u32>> {
    // delnodes.dmp: list of deleted tax_ids
    let mut delnodes = HashSet::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line = line.with_context(|| format!("failed to read {} line {line_number}", source))?;
        if line.trim().is_empty() {
            continue;
        }
        let tax_id = line
            .split('|')
            .next()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing taxid in {} line {line_number}", source))?
            .parse::<u32>()
            .with_context(|| format!("invalid taxid in {} line {line_number}", source))?;
        if !delnodes.insert(tax_id) {
            bail!(
                "duplicate deleted taxid {tax_id} in {} line {line_number}",
                source
            );
        }
    }
    Ok(delnodes)
}

#[cfg(test)]
mod tests {
    use super::Taxonomy;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::collections::{HashMap, HashSet};
    use std::io::Write;

    fn taxonomy(merged: HashMap<u32, u32>, deleted: HashSet<u32>) -> Taxonomy {
        Taxonomy::from_parts(
            HashMap::from([
                (1, (1, "no rank".to_string())),
                (10, (1, "genus".to_string())),
                (20, (10, "species".to_string())),
                (30, (20, "subspecies".to_string())),
            ]),
            HashMap::from([
                (1, "root".to_string()),
                (10, "Genus ten".to_string()),
                (20, "Species twenty".to_string()),
                (30, "Subspecies thirty".to_string()),
            ]),
            merged,
            deleted,
        )
    }

    #[test]
    fn resolves_complete_merged_taxid_chains() {
        let taxonomy = taxonomy(HashMap::from([(28, 29), (29, 30)]), HashSet::new());
        assert_eq!(taxonomy.resolve_current_taxid(28), Some(30));
        assert_eq!(taxonomy.get_name(28), Some("Subspecies thirty"));
        assert_eq!(taxonomy.get_rank(28), Some("subspecies"));
        assert_eq!(
            taxonomy
                .get_parents(28)
                .into_iter()
                .map(|(taxid, _, _)| taxid)
                .collect::<Vec<_>>(),
            vec![30, 20, 10, 1]
        );
    }

    #[test]
    fn rejects_deleted_missing_and_cyclic_taxids() {
        let taxonomy = taxonomy(
            HashMap::from([(40, 41), (41, 40), (50, 20)]),
            HashSet::from([20]),
        );
        assert_eq!(taxonomy.resolve_current_taxid(40), None);
        assert_eq!(taxonomy.resolve_current_taxid(50), None);
        assert_eq!(taxonomy.resolve_current_taxid(999), None);
        assert!(taxonomy.get_parents(40).is_empty());
    }

    #[test]
    fn loads_gzipped_and_plain_dump_files() {
        let temp = tempfile::tempdir().unwrap();
        let gzip = |name: &str, contents: &str| {
            let file = std::fs::File::create(temp.path().join(name)).unwrap();
            let mut encoder = GzEncoder::new(file, Compression::fast());
            encoder.write_all(contents.as_bytes()).unwrap();
            encoder.finish().unwrap();
        };
        gzip(
            "nodes.dmp.gz",
            "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tspecies\t|\n",
        );
        gzip(
            "names.dmp.gz",
            "1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tSpecies two\t|\t\t|\tscientific name\t|\n",
        );
        gzip("merged.dmp.gz", "3\t|\t2\t|\n");
        std::fs::write(temp.path().join("delnodes.dmp"), "4\t|\n").unwrap();

        let taxonomy = Taxonomy::from_dir(temp.path()).unwrap();
        assert_eq!(taxonomy.get_name(3), Some("Species two"));
        assert_eq!(taxonomy.resolve_current_taxid(4), None);
    }
}

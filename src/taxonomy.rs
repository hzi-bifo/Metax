mod dump;

use anyhow::{bail, Context, Result};
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub(crate) use dump::TaxonomySource;

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
    pub(crate) fn prepare_dmp(dmp_dir: Option<&str>) -> Result<TaxonomySource> {
        TaxonomySource::prepare(dmp_dir)
    }

    pub(super) fn from_dir(dir: &Path) -> Result<Self> {
        let nodes = dir.join("nodes.dmp");
        let names = dir.join("names.dmp");
        let merged = dir.join("merged.dmp");
        let delnodes = dir.join("delnodes.dmp");
        Self::from_files(&nodes, &names, &merged, Some(&delnodes))
    }

    pub fn from_files(
        nodes: &Path,
        names: &Path,
        merged: &Path,
        delnodes: Option<&Path>,
    ) -> Result<Self> {
        // Load each taxonomy table into memory for fast lookups.
        let nodes_iter = load_nodes(nodes)?;
        let names_iter = load_names(names)?;
        let merged_iter = load_merged(merged)?;
        let del_nodes = if let Some(path) = delnodes {
            load_delnodes(path)?
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

fn load_nodes(path: &Path) -> Result<HashMap<u32, (u32, String)>> {
    // nodes.dmp: tax_id | parent_tax_id | rank
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut nodes = HashMap::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line =
            line.with_context(|| format!("failed to read {} line {line_number}", path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('|').map(str::trim);
        let tax_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing taxid in {} line {line_number}", path.display()))?
            .parse::<u32>()
            .with_context(|| format!("invalid taxid in {} line {line_number}", path.display()))?;
        let parent_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| {
                format!(
                    "missing parent taxid in {} line {line_number}",
                    path.display()
                )
            })?
            .parse::<u32>()
            .with_context(|| {
                format!(
                    "invalid parent taxid in {} line {line_number}",
                    path.display()
                )
            })?;
        let rank = fields
            .next()
            .with_context(|| format!("missing rank in {} line {line_number}", path.display()))?
            .to_string();
        if nodes.insert(tax_id, (parent_id, rank)).is_some() {
            bail!(
                "duplicate taxid {tax_id} in {} line {line_number}",
                path.display()
            );
        }
    }
    Ok(nodes)
}

fn load_names(path: &Path) -> Result<HashMap<u32, String>> {
    // names.dmp: tax_id | name | ... | name_class (keep scientific name only)
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut names = HashMap::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line =
            line.with_context(|| format!("failed to read {} line {line_number}", path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('|').map(str::trim);
        let tax_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing taxid in {} line {line_number}", path.display()))?
            .parse::<u32>()
            .with_context(|| format!("invalid taxid in {} line {line_number}", path.display()))?;
        let name = fields
            .next()
            .with_context(|| format!("missing name in {} line {line_number}", path.display()))?;
        let _unique_name = fields.next().with_context(|| {
            format!(
                "missing unique-name field in {} line {line_number}",
                path.display()
            )
        })?;
        let name_class = fields.next().with_context(|| {
            format!(
                "missing name class in {} line {line_number}",
                path.display()
            )
        })?;
        if name_class == "scientific name" && names.insert(tax_id, name.to_string()).is_some() {
            bail!(
                "duplicate scientific name for taxid {tax_id} in {} line {line_number}",
                path.display()
            );
        }
    }
    Ok(names)
}

fn load_merged(path: &Path) -> Result<HashMap<u32, u32>> {
    // merged.dmp: old_tax_id | new_tax_id
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut merged = HashMap::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line =
            line.with_context(|| format!("failed to read {} line {line_number}", path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('|').map(str::trim);
        let old_tax_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing old taxid in {} line {line_number}", path.display()))?
            .parse::<u32>()
            .with_context(|| {
                format!("invalid old taxid in {} line {line_number}", path.display())
            })?;
        let new_tax_id = fields
            .next()
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing new taxid in {} line {line_number}", path.display()))?
            .parse::<u32>()
            .with_context(|| {
                format!("invalid new taxid in {} line {line_number}", path.display())
            })?;
        if merged.insert(old_tax_id, new_tax_id).is_some() {
            bail!(
                "duplicate merged taxid {old_tax_id} in {} line {line_number}",
                path.display()
            );
        }
    }
    Ok(merged)
}

fn load_delnodes(path: &Path) -> Result<HashSet<u32>> {
    // delnodes.dmp: list of deleted tax_ids
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut delnodes = HashSet::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line =
            line.with_context(|| format!("failed to read {} line {line_number}", path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let tax_id = line
            .split('|')
            .next()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .with_context(|| format!("missing taxid in {} line {line_number}", path.display()))?
            .parse::<u32>()
            .with_context(|| format!("invalid taxid in {} line {line_number}", path.display()))?;
        if !delnodes.insert(tax_id) {
            bail!(
                "duplicate deleted taxid {tax_id} in {} line {line_number}",
                path.display()
            );
        }
    }
    Ok(delnodes)
}

#[cfg(test)]
mod tests {
    use super::Taxonomy;
    use std::collections::{HashMap, HashSet};

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
}

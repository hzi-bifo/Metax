use anyhow::{Context, Result};
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

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
    pub fn load_from_dmp(dmp_dir: Option<&str>) -> Result<Self> {
        // Default to a bundled "data" directory next to the executable.
        let dir = if let Some(dir) = dmp_dir {
            PathBuf::from(dir)
        } else {
            let exe = std::env::current_exe().context("failed to locate executable")?;
            exe.parent()
                .map(|p| p.join("data"))
                .context("failed to derive default taxonomy path")?
        };
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
        Ok(Self {
            nodes: nodes_iter,
            names: names_iter,
            merged: merged_iter,
            del_nodes,
        })
    }

    pub fn get_name(&self, tax_id: u32) -> Option<&str> {
        // Apply merged/deleted mappings before returning a display name.
        let tax_id = self.merged.get(&tax_id).copied().unwrap_or(tax_id);
        if self.del_nodes.contains(&tax_id) {
            return None;
        }
        self.names.get(&tax_id).map(|s| s.as_str())
    }

    pub fn get_rank(&self, tax_id: u32) -> Option<&str> {
        // Apply merged/deleted mappings before returning a rank.
        let tax_id = self.merged.get(&tax_id).copied().unwrap_or(tax_id);
        if self.del_nodes.contains(&tax_id) {
            return None;
        }
        self.nodes.get(&tax_id).map(|(_, rank)| rank.as_str())
    }

    pub fn get_parents(&self, mut tax_id: u32) -> Vec<(u32, String, String)> {
        // Walk up the taxonomy to root, returning the lineage.
        let mut parents = Vec::new();
        loop {
            let current = self.merged.get(&tax_id).copied().unwrap_or(tax_id);
            if let Some((parent_id, rank)) = self.nodes.get(&current) {
                let name = self.names.get(&current).cloned().unwrap_or_default();
                parents.push((current, name, rank.clone()));
                if *parent_id == current {
                    break;
                }
                tax_id = *parent_id;
                if tax_id == 1 {
                    parents.push((1, "root".to_string(), "no rank".to_string()));
                    break;
                }
            } else {
                break;
            }
        }
        parents
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
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('|').collect();
        if parts.len() < 3 {
            continue;
        }
        let tax_id = parts[0].trim().parse::<u32>()?;
        let parent_id = parts[1].trim().parse::<u32>()?;
        let rank = parts[2].trim().to_string();
        nodes.insert(tax_id, (parent_id, rank));
    }
    Ok(nodes)
}

fn load_names(path: &Path) -> Result<HashMap<u32, String>> {
    // names.dmp: tax_id | name | ... | name_class (keep scientific name only)
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut names = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('|').collect();
        if parts.len() < 4 {
            continue;
        }
        if parts[3].contains("scientific name") {
            let tax_id = parts[0].trim().parse::<u32>()?;
            let name = parts[1].trim().to_string();
            names.insert(tax_id, name);
        }
    }
    Ok(names)
}

fn load_merged(path: &Path) -> Result<HashMap<u32, u32>> {
    // merged.dmp: old_tax_id | new_tax_id
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut merged = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('|').collect();
        if parts.len() < 2 {
            continue;
        }
        let old_tax_id = parts[0].trim().parse::<u32>()?;
        let new_tax_id = parts[1].trim().parse::<u32>()?;
        merged.insert(old_tax_id, new_tax_id);
    }
    Ok(merged)
}

fn load_delnodes(path: &Path) -> Result<HashSet<u32>> {
    // delnodes.dmp: list of deleted tax_ids
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut delnodes = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('|').collect();
        if let Some(id) = parts.first() {
            if let Ok(tax_id) = id.trim().parse::<u32>() {
                delnodes.insert(tax_id);
            }
        }
    }
    Ok(delnodes)
}

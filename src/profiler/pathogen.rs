//! Pathogen-host annotation table shared by both profilers.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{Context, Result};

/// Pathogen metadata from the pathogen-host mapping file.
///
/// Annotates results with host and disease information when pathogen
/// filtering is enabled. Each field holds a semicolon-separated list.
#[derive(Clone, Debug)]
pub(crate) struct PathogenEntry {
    pub(crate) host_taxids: String,
    pub(crate) host_names: String,
    pub(crate) diseases: String,
}

/// Load the pathogen-host mapping from a TSV file keyed by pathogen taxid.
///
/// The first line is a header. Rows with fewer than four columns are skipped,
/// so a table carrying extra trailing columns still loads.
pub(crate) fn load_pathogen_table(path: &str) -> Result<HashMap<u32, PathogenEntry>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open pathogen host table: {}", path))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for (idx, line) in reader.lines().enumerate() {
        let line = line?;
        if idx == 0 {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }
        let pathogen_taxid: u32 = parts[0].parse()?;
        map.insert(
            pathogen_taxid,
            PathogenEntry {
                host_taxids: parts[1].to_string(),
                host_names: parts[2].to_string(),
                diseases: parts[3].to_string(),
            },
        );
    }
    Ok(map)
}

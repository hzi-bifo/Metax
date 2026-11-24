//! `--host`: a profile row is judged by its own taxid in the pathogen table
//! (`--pathogen-table`). A taxon the table lists is kept when one of its hosts
//! is the given host or lies below it in the taxonomy, and is annotated with
//! the entry's host names, host taxids and diseases; a taxon the table lists
//! for other hosts only is dropped. A taxon the table does not list, or lists
//! without hosts (an empty or `unknown` host field), is kept and annotated
//! `unknown`; a strain row with such a taxid follows its species. Rows are
//! dropped before abundance is normalised; their reads go to no other taxon.
//! The profile gains the three annotation columns, and its files (profile,
//! decisions) take the prefix `<prefix>.pathogen`.

use std::collections::HashMap;

use super::higher::Lineages;
use super::{Dropped, Removal, TaxonBlock};
use crate::profiler::output::{ProfileAnnotations, ProfileRow};
use crate::profiler::pathogen::PathogenEntry;

/// What the pathogen table says about a taxon for the host.
#[derive(Clone, Debug, PartialEq)]
enum Verdict {
    /// Not in the table: kept, annotated `unknown`.
    Unlisted,
    /// A pathogen of the host: kept with its entry's annotations.
    Pathogen(ProfileAnnotations),
    /// Listed for other hosts only: dropped.
    OtherHosts,
}

/// The host, the pathogen table, and the taxa's lineages.
pub(super) struct HostFilter<'a> {
    pub(super) host: u32,
    pub(super) table: &'a HashMap<u32, PathogenEntry>,
    pub(super) lineages: &'a Lineages<'a>,
}

/// The distinct items of a `;`-separated field, in their order; `None` when
/// the field is empty or `unknown`.
fn items(field: &str) -> Option<Vec<&str>> {
    if field.is_empty() || field == "unknown" {
        return None;
    }
    let mut items: Vec<&str> = Vec::new();
    for item in field.split(';') {
        if !items.contains(&item) {
            items.push(item);
        }
    }
    Some(items)
}

impl HostFilter<'_> {
    fn verdict(&self, taxid: u32) -> Verdict {
        let Some(entry) = self.table.get(&taxid) else {
            return Verdict::Unlisted;
        };
        // A taxon listed without hosts is judged as an unlisted one.
        let Some(listed) = items(&entry.host_taxids) else {
            return Verdict::Unlisted;
        };
        let hosts = listed
            .into_iter()
            .filter_map(|item| item.parse::<u32>().ok())
            .collect::<Vec<_>>();
        let of_host = hosts.iter().any(|&listed| {
            listed == self.host
                || self
                    .lineages
                    .of(listed)
                    .iter()
                    .any(|node| node.0 == self.host)
        });
        if !of_host {
            return Verdict::OtherHosts;
        }
        let joined =
            |items: Option<Vec<&str>>| items.map_or_else(|| "NA".to_string(), |i| i.join(";"));
        let host_taxids = hosts.iter().map(u32::to_string).collect::<Vec<_>>();
        Verdict::Pathogen(ProfileAnnotations {
            host_names: joined(items(&entry.host_names)),
            host_taxids: host_taxids.join(";"),
            diseases: joined(items(&entry.diseases)),
        })
    }

    /// Annotate `row` by `verdict`; `false` when the row is dropped.
    fn keep(row: &mut ProfileRow, verdict: &Verdict) -> bool {
        row.annotations = Some(match verdict {
            Verdict::Unlisted => ProfileAnnotations {
                host_names: "unknown".to_string(),
                host_taxids: "unknown".to_string(),
                diseases: "unknown".to_string(),
            },
            Verdict::Pathogen(annotations) => annotations.clone(),
            Verdict::OtherHosts => return false,
        });
        true
    }

    /// Annotate the rows of `blocks` and drop those that are not pathogens of
    /// the host, and the blocks left without a row. Records the taxa dropped
    /// whole, and every strain row dropped by its genome.
    pub(super) fn apply(&self, blocks: &mut Vec<TaxonBlock>, dropped: &mut Dropped) {
        let removal = Removal::NotPathogen(self.host);
        blocks.retain_mut(|block| {
            let own = self.verdict(block.taxid);
            if let Some(row) = &mut block.taxon {
                if !Self::keep(row, &own) {
                    block.taxon = None;
                }
            }
            block.genomes.retain_mut(|row| {
                let verdict = match self.verdict(row.taxid) {
                    Verdict::Unlisted => own.clone(),
                    listed => listed,
                };
                let kept = Self::keep(row, &verdict);
                if !kept {
                    dropped.genomes.insert(row.best_fit.assembly_id, removal);
                }
                kept
            });
            let kept = block.taxon.is_some() || !block.genomes.is_empty();
            if !kept {
                dropped.taxa.insert(block.taxid, removal);
            }
            kept
        });
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::super::profile_row;
    use std::sync::Arc;

    use super::*;
    use crate::profiler::evidence::EvidenceCoverage;
    use crate::taxonomy::Taxonomy;

    /// Primates 9443 > Homo 9605 > human 9606; mouse 10090 elsewhere.
    fn taxonomy() -> Taxonomy {
        let node = |parent: u32, rank: &str| (parent, rank.to_string());
        Taxonomy::from_parts(
            HashMap::from([
                (1, node(1, "no rank")),
                (9443, node(1, "order")),
                (9605, node(9443, "genus")),
                (9606, node(9605, "species")),
                (10090, node(1, "species")),
            ]),
            HashMap::new(),
            HashMap::new(),
            HashSet::new(),
        )
    }

    fn entry(hosts: &str, names: &str, diseases: &str) -> PathogenEntry {
        PathogenEntry {
            host_taxids: hosts.to_string(),
            host_names: names.to_string(),
            diseases: diseases.to_string(),
        }
    }

    /// Species 1 (a human pathogen, with strain rows: one under its own taxid, one under strain 11, one under
    /// strain 12 listed for mice only), species 2 (mice only), species 3 (not listed), genus 100 (not listed),
    /// family 1000 (human), family 2000 (mice only).
    fn table() -> HashMap<u32, PathogenEntry> {
        HashMap::from([
            (
                1,
                entry(
                    "9606;9606",
                    "Homo sapiens;Homo sapiens",
                    "meningitis;sepsis",
                ),
            ),
            (12, entry("10090", "Mus musculus", "unknown")),
            (2, entry("10090", "Mus musculus", "unknown")),
            (1000, entry("9606", "Homo sapiens", "unknown")),
            (2000, entry("10090", "Mus musculus", "unknown")),
        ])
    }

    fn row(taxid: u32, rank: &str, genome: u32) -> ProfileRow {
        profile_row(
            format!("Taxon {taxid}"),
            taxid,
            rank,
            (10.0, 1.0),
            &EvidenceCoverage::default(),
            (genome, Arc::from(format!("GCF_{genome}"))),
            (f64::NAN, f64::NAN),
        )
    }

    fn block(taxid: u32, taxon: Option<ProfileRow>, genomes: Vec<ProfileRow>) -> TaxonBlock {
        TaxonBlock {
            taxid,
            unit: taxon.as_ref().is_none_or(|row| row.rank == "species"),
            taxon,
            genomes,
            breadth: 1.0,
        }
    }

    #[test]
    fn a_listed_taxon_is_kept_for_its_hosts_and_their_ancestors_only() {
        let table = table();
        let taxonomy = taxonomy();
        let lineages = Lineages::new(&taxonomy);
        let filter = |host| HostFilter {
            host,
            table: &table,
            lineages: &lineages,
        };
        let pathogen = filter(9606).verdict(1);
        assert_eq!(
            pathogen,
            Verdict::Pathogen(ProfileAnnotations {
                host_names: "Homo sapiens".to_string(),
                host_taxids: "9606".to_string(),
                diseases: "meningitis;sepsis".to_string(),
            })
        );
        // The genus of the listed host counts as the host.
        assert!(matches!(filter(9605).verdict(1), Verdict::Pathogen(_)));
        assert_eq!(filter(9606).verdict(2), Verdict::OtherHosts);
        assert_eq!(filter(10090).verdict(1), Verdict::OtherHosts);
        assert_eq!(filter(9606).verdict(3), Verdict::Unlisted);
    }

    #[test]
    fn rows_are_annotated_or_dropped_and_strains_follow_their_species() {
        let table = table();
        let taxonomy = taxonomy();
        let lineages = Lineages::new(&taxonomy);
        let mut blocks = vec![
            block(
                1,
                Some(row(1, "species", 10)),
                vec![
                    row(1, "strain", 10),
                    row(11, "strain", 11),
                    row(12, "strain", 12),
                ],
            ),
            block(2, Some(row(2, "species", 20)), vec![row(21, "strain", 21)]),
            block(3, Some(row(3, "species", 30)), Vec::new()),
            block(100, Some(row(100, "genus", 40)), Vec::new()),
            block(1000, Some(row(1000, "family", 50)), Vec::new()),
            block(2000, Some(row(2000, "family", 60)), Vec::new()),
        ];
        let mut dropped = Dropped::default();
        HostFilter {
            host: 9606,
            table: &table,
            lineages: &lineages,
        }
        .apply(&mut blocks, &mut dropped);
        let kept = blocks
            .iter()
            .flat_map(|block| block.taxon.iter().chain(&block.genomes))
            .map(|row| {
                let annotations = row.annotations.as_ref().unwrap();
                (row.taxid, row.rank.as_str(), annotations.diseases.as_str())
            })
            .collect::<Vec<_>>();
        assert_eq!(
            kept,
            [
                (1, "species", "meningitis;sepsis"),
                (1, "strain", "meningitis;sepsis"),
                // Strain 11 is not listed: it follows its species.
                (11, "strain", "meningitis;sepsis"),
                (3, "species", "unknown"),
                (100, "genus", "unknown"),
                (1000, "family", "NA"),
            ]
        );
        // Strain 12 is listed for mice only; species 2 goes with its strain row; family 2000.
        let mut taxa = dropped.taxa.keys().copied().collect::<Vec<_>>();
        taxa.sort_unstable();
        assert_eq!(taxa, [2, 2000]);
        let mut genomes = dropped.genomes.keys().copied().collect::<Vec<_>>();
        genomes.sort_unstable();
        assert_eq!(genomes, [12, 21]);
        assert_eq!(
            dropped.taxa[&2].reason(),
            "removed by --host: not a pathogen of host 9606"
        );
    }

    #[test]
    fn a_taxon_listed_without_hosts_is_judged_as_unlisted() {
        let table = HashMap::from([
            (4, entry("", "Homo sapiens", "sepsis")),
            (5, entry("unknown", "unknown", "unknown")),
            (41, entry("", "", "")),
        ]);
        let taxonomy = taxonomy();
        let lineages = Lineages::new(&taxonomy);
        let filter = HostFilter {
            host: 9606,
            table: &table,
            lineages: &lineages,
        };
        assert_eq!(filter.verdict(4), Verdict::Unlisted);
        assert_eq!(filter.verdict(5), Verdict::Unlisted);
        let mut dropped = Dropped::default();
        let mut blocks = vec![block(
            4,
            Some(row(4, "species", 40)),
            vec![row(41, "strain", 41)],
        )];
        filter.apply(&mut blocks, &mut dropped);
        let host_taxids = |row: &ProfileRow| row.annotations.as_ref().unwrap().host_taxids.clone();
        assert_eq!(host_taxids(blocks[0].taxon.as_ref().unwrap()), "unknown");
        assert_eq!(host_taxids(&blocks[0].genomes[0]), "unknown");
        // A strain row listed without hosts follows its species, here one
        // listed for mice only.
        let mice = HashMap::from([
            (2, entry("10090", "Mus musculus", "unknown")),
            (41, entry("", "", "")),
        ]);
        let mut blocks = vec![block(
            2,
            Some(row(2, "species", 20)),
            vec![row(41, "strain", 42)],
        )];
        HostFilter {
            host: 9606,
            table: &mice,
            lineages: &lineages,
        }
        .apply(&mut blocks, &mut dropped);
        assert!(blocks.is_empty());
        assert_eq!(dropped.genomes.keys().copied().collect::<Vec<_>>(), [42]);
    }

    #[test]
    fn abundance_is_normalised_over_the_rows_kept() {
        let table = table();
        let taxonomy = taxonomy();
        let lineages = Lineages::new(&taxonomy);
        let mut blocks = vec![
            block(1, Some(row(1, "species", 10)), vec![row(11, "strain", 11)]),
            block(2, Some(row(2, "species", 20)), Vec::new()),
            block(3, Some(row(3, "species", 30)), Vec::new()),
        ];
        HostFilter {
            host: 9606,
            table: &table,
            lineages: &lineages,
        }
        .apply(&mut blocks, &mut Dropped::default());
        super::super::set_abundance(&mut blocks);
        // Three rows of depth 1 are left: a third each.
        let abundances = blocks
            .iter()
            .flat_map(|block| block.taxon.iter().chain(&block.genomes))
            .map(|row| row.abundance)
            .collect::<Vec<_>>();
        assert_eq!(abundances.len(), 3);
        assert!(abundances
            .iter()
            .all(|&abundance| (abundance - 100.0 / 3.0).abs() < 1e-9));
    }
}

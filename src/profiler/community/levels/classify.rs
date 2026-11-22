//! Per-read classifications (`--save-classification`,
//! `<prefix>.classify.txt.gz`): one row per fragment (a read pair counts once)
//! with its read id, the taxid and rank of its largest share, the taxids of its
//! shares and the shares, largest first (equal shares by taxid); a fragment
//! with a single share has share `1`.
//!
//! A fragment's shares are the ones the profile's rows count it with:
//!
//! - a fragment with candidate species at the species identity is shared
//!   among the kept species by the second EM (its read group's posteriors). A
//!   species reported at species counts its share in its species row, and a
//!   unit whose references are labelled above species in its row at that
//!   rank; with `--strain`, a fragment that the genome EM shares among the
//!   genomes reported on their own splits its species share among their strain
//!   rows by that EM's posteriors. A species reported at its genus counts its share
//!   in the genus row; a species left out, or dropped by `--host`,
//!   `--min-reads` or `--min-breadth`, counts it nowhere;
//! - a fragment none of whose candidate species is kept goes whole to the
//!   roll-up's genus row of its best candidate, when that species' reads were
//!   rolled up;
//! - a fragment without a species-level candidate is, on noisy long reads,
//!   credited to the kept species tied with its best match, shared by their
//!   weights and counted as above; otherwise the family level's EM shares it
//!   among the families (its posteriors), counted in the family rows.
//!
//! Shares on taxa without a row in the profile are dropped and the others
//! renormalised; a fragment without any is not written, as are those of a
//! database of sampled genomes (no file is written). Rows follow the mapping
//! file's read order.
//!
//! Summed over the file, a taxon's shares equal its reads in the profile
//! (column 4), as the profile counts each fragment the same way, except where
//! shares were renormalised and where a species row is left out for holding
//! less than one read.

use std::collections::HashMap;
use std::io::Write;
use std::sync::Arc;

use anyhow::Result;

use super::genomes::SpeciesGenomes;
use super::{Fragment, TaxonBlock, UnitShares};
use crate::profiler::community::MappingData;
use crate::profiler::output::{create_classification_writer, finish_classification_writer};
use crate::profiler::rollup::GenusRow;

/// A fragment's shares: taxid and rank of a profile row, and the share.
type Shares<'r> = Vec<(u32, &'r str, f64)>;

/// A read group's posteriors over units.
type Posteriors = [(u32, f64)];

/// The profile's rows a share can go to.
struct ProfileTaxa<'r> {
    /// The rank of each taxon's row: a unit's row (species, or the taxon its
    /// references are labelled with), a genus or a family row.
    rows: HashMap<u32, &'r str>,
    /// Taxid of each genome's strain row, by assembly id.
    strain_rows: HashMap<u32, u32>,
    /// Genus row of each species reported at its genus.
    at_genus: HashMap<u32, u32>,
}

impl<'r> ProfileTaxa<'r> {
    fn of(
        blocks: &'r [TaxonBlock],
        genome_level: &HashMap<u32, SpeciesGenomes>,
        genus_members: &[(u32, u32, f64)],
        assembly_names: &[Arc<str>],
    ) -> Self {
        let mut taxa = Self {
            rows: HashMap::new(),
            strain_rows: HashMap::new(),
            at_genus: HashMap::new(),
        };
        for block in blocks {
            if let Some(row) = &block.taxon {
                taxa.rows.insert(block.taxid, row.rank.as_str());
            }
            let strains = genome_level
                .get(&block.taxid)
                .map_or(&[][..], |g| &g.strains);
            for row in &block.genomes {
                let genome = row.best_fit.as_ref().map(|b| &*b.genome);
                if let Some(strain) = strains
                    .iter()
                    .find(|s| Some(&*assembly_names[s.assembly_id as usize]) == genome)
                {
                    taxa.strain_rows.insert(strain.assembly_id, row.taxid);
                }
            }
        }
        for &(genus, species, _) in genus_members {
            if taxa.rows.contains_key(&genus) {
                taxa.at_genus.insert(species, genus);
            }
        }
        taxa
    }
}

/// How the engine's levels share the fragments.
pub(super) struct Classifications<'a> {
    /// The fragments with candidate species, each with its read group in the
    /// second EM over the kept species (`None` without a kept candidate), and
    /// the groups' posteriors.
    pub(super) fragments: &'a [Fragment],
    pub(super) group_of: &'a [Option<usize>],
    pub(super) posteriors: &'a [Vec<(u32, f64)>],
    pub(super) genome_level: &'a HashMap<u32, SpeciesGenomes>,
    /// Species reported at their genus: (genus, species, reads).
    pub(super) genus_members: &'a [(u32, u32, f64)],
    pub(super) rollup_rows: &'a [GenusRow],
    /// The fragments offered to the roll-up: read and best candidate.
    pub(super) orphans: &'a [(usize, u32)],
    /// Credited fragments: read and shares of the kept species.
    pub(super) credited: &'a [(usize, Vec<(u32, f64)>)],
    pub(super) families: Option<&'a UnitShares>,
}

impl Classifications<'_> {
    /// Write `path` for the profile's `blocks`; returns the fragments written.
    pub(super) fn write(
        &self,
        path: &str,
        blocks: &[TaxonBlock],
        data: &MappingData,
    ) -> Result<usize> {
        let rows = self.shares(blocks, &data.assembly_names);
        let mut writer = create_classification_writer(path)?;
        for (read, shares) in &rows {
            writeln!(
                writer,
                "{}",
                classification_row(&data.read_alignments[*read].read_name, shares)
            )?;
        }
        finish_classification_writer(writer)?;
        Ok(rows.len())
    }

    /// Each written fragment's read and its shares on the profile's rows,
    /// renormalised, in read order.
    fn shares<'r>(
        &self,
        blocks: &'r [TaxonBlock],
        assembly_names: &[Arc<str>],
    ) -> Vec<(usize, Shares<'r>)> {
        let taxa = ProfileTaxa::of(
            blocks,
            self.genome_level,
            self.genus_members,
            assembly_names,
        );
        // The genome EM's posteriors of each fragment it shares, by read.
        let mut genome_shares: HashMap<usize, (u32, &Posteriors)> = HashMap::new();
        for (&species, genomes) in self.genome_level {
            for &(read, group) in &genomes.shares.members {
                genome_shares.insert(read, (species, &genomes.shares.posteriors[group]));
            }
        }
        let species_share = |shares: &mut Shares<'r>, read: usize, species: u32, share: f64| {
            if let Some(&genus) = taxa.at_genus.get(&species) {
                shares.push((genus, "genus", share));
            } else if let Some(&(_, posteriors)) =
                genome_shares.get(&read).filter(|(of, _)| *of == species)
            {
                for &(genome, posterior) in posteriors {
                    if let Some(&taxid) = taxa.strain_rows.get(&genome) {
                        shares.push((taxid, "strain", share * posterior));
                    }
                }
            } else if let Some(&rank) = taxa.rows.get(&species) {
                shares.push((species, rank, share));
            }
        };
        let mut rows = Vec::new();
        for (fragment, group) in self.fragments.iter().zip(self.group_of) {
            let Some(group) = group else {
                continue;
            };
            let mut shares = Vec::new();
            for &(species, posterior) in &self.posteriors[*group] {
                species_share(&mut shares, fragment.read, species, posterior);
            }
            rows.extend(renormalised(shares).map(|shares| (fragment.read, shares)));
        }
        for (read, credited) in self.credited {
            let mut shares = Vec::new();
            for &(species, share) in credited {
                species_share(&mut shares, *read, species, share);
            }
            rows.extend(renormalised(shares).map(|shares| (*read, shares)));
        }
        let genus_of_donor: HashMap<u32, u32> = self
            .rollup_rows
            .iter()
            .flat_map(|row| row.donors.iter().map(|&donor| (donor, row.genus)))
            .collect();
        for &(read, donor) in self.orphans {
            if let Some(&genus) = genus_of_donor
                .get(&donor)
                .filter(|genus| taxa.rows.contains_key(genus))
            {
                rows.push((read, vec![(genus, "genus", 1.0)]));
            }
        }
        if let Some(families) = self.families {
            for (&read, group) in families.reads.iter().zip(&families.group_of) {
                let Some(group) = group else {
                    continue;
                };
                let shares = families.posteriors[*group]
                    .iter()
                    .filter(|(family, _)| taxa.rows.contains_key(family))
                    .map(|&(family, posterior)| (family, "family", posterior))
                    .collect();
                rows.extend(renormalised(shares).map(|shares| (read, shares)));
            }
        }
        rows.sort_unstable_by_key(|row| row.0);
        rows
    }
}

/// `shares` with those on one row summed, renormalised to 1, largest first
/// (equal shares by taxid, then rank); `None` without a positive share.
fn renormalised(mut shares: Shares) -> Option<Shares> {
    shares.retain(|share| share.2.is_finite() && share.2 > 0.0);
    shares.sort_by_key(|share| (share.0, share.1));
    shares.dedup_by(|later, kept| {
        let same = (later.0, later.1) == (kept.0, kept.1);
        if same {
            kept.2 += later.2;
        }
        same
    });
    let total: f64 = shares.iter().map(|share| share.2).sum();
    if total <= 0.0 {
        return None;
    }
    for share in &mut shares {
        share.2 /= total;
    }
    shares.sort_by(|a, b| b.2.total_cmp(&a.2).then(a.0.cmp(&b.0)).then(a.1.cmp(b.1)));
    Some(shares)
}

/// One row of the file: read id, taxid and rank of the largest share, the
/// taxids, and the shares (`1` for a single one, else five decimals).
fn classification_row(read: &str, shares: &Shares) -> String {
    let (taxid, rank, _) = shares[0];
    let taxids = shares
        .iter()
        .map(|share| share.0.to_string())
        .collect::<Vec<_>>()
        .join(";");
    let fractions = if shares.len() == 1 {
        "1".to_string()
    } else {
        shares
            .iter()
            .map(|share| format!("{:.5}", share.2))
            .collect::<Vec<_>>()
            .join(";")
    };
    format!("{read}\t{taxid}\t{rank}\t{taxids}\t{fractions}")
}

#[cfg(test)]
mod tests {
    use super::super::genomes::{StrainGenome, StrainShares};
    use super::super::profile_row;
    use super::*;
    use crate::profiler::evidence::EvidenceCoverage;

    /// A profile row of `rank` for `taxid` on genome `GCF_<genome>`.
    fn row(taxid: u32, rank: &str, genome: u32) -> crate::profiler::output::ProfileRow {
        profile_row(
            format!("Taxon {taxid}"),
            taxid,
            rank,
            (1.0, 1.0),
            &EvidenceCoverage::default(),
            Arc::from(format!("GCF_{genome}")),
            Some((f64::NAN, f64::NAN)),
        )
    }

    fn block(
        taxid: u32,
        taxon: Option<crate::profiler::output::ProfileRow>,
        genomes: Vec<crate::profiler::output::ProfileRow>,
    ) -> TaxonBlock {
        TaxonBlock {
            taxid,
            unit: taxon.as_ref().is_none_or(|row| row.rank == "species"),
            taxon,
            genomes,
            breadth: 1.0,
        }
    }

    fn fragment(read: usize) -> Fragment {
        Fragment {
            read,
            candidates: Vec::new(),
            bases: 100.0,
        }
    }

    fn strain(assembly_id: u32) -> StrainGenome {
        StrainGenome {
            assembly_id,
            reads: 1.0,
            bases: 1.0,
            mean_identity: 1.0,
            coverage: EvidenceCoverage::default(),
            asm_len: 1.0,
        }
    }

    /// Accessions `GCF_<id>` of assembly ids 0-99.
    fn accessions() -> Vec<Arc<str>> {
        (0..100).map(|id| Arc::from(format!("GCF_{id}"))).collect()
    }

    /// The profile: species 1 (species row), species 4 with a species row
    /// and strain rows on genomes 40 (strain 400) and 41 (under the species'
    /// taxid), species 2 at genus 100, genus 100 (with the roll-up of species
    /// 5), and family 1000. Species 3 is kept but not reported (left out).
    fn blocks() -> Vec<TaxonBlock> {
        vec![
            block(1, Some(row(1, "species", 10)), Vec::new()),
            block(
                4,
                Some(row(4, "species", 40)),
                vec![row(400, "strain", 40), row(4, "strain", 41)],
            ),
            block(100, Some(row(100, "genus", 20)), Vec::new()),
            block(1000, Some(row(1000, "family", 80)), Vec::new()),
        ]
    }

    #[test]
    fn fragments_are_labelled_with_the_profile_rows_that_count_them() {
        let fragments = (0..6).map(fragment).collect::<Vec<_>>();
        // Read 0: species 1 only. Read 1: species 1 and species 3 (left out). Read 2: species 2 (at genus 100).
        // Read 3: species 4, a fragment the genome EM shares among its genomes. Read 4: species 4, not shared by the
        // genome EM. Read 5: no kept candidate (roll-up).
        let group_of = [Some(0), Some(1), Some(2), Some(3), Some(3), None];
        let posteriors = vec![
            vec![(1, 1.0)],
            vec![(1, 0.6), (3, 0.4)],
            vec![(2, 1.0)],
            vec![(4, 0.9), (1, 0.1)],
        ];
        let genome_level = HashMap::from([(
            4,
            SpeciesGenomes {
                best_fit: 40,
                mean_identity: 1.0,
                strains: vec![strain(40), strain(41)],
                shares: StrainShares {
                    members: vec![(3, 0)],
                    posteriors: vec![vec![(40, 0.75), (41, 0.25)]],
                },
            },
        )]);
        let rollup_rows = vec![GenusRow {
            genus: 100,
            fragments: 1.0,
            equivalent: 1.0,
            depth: 0.0,
            mapped_len: 0.0,
            witness: EvidenceCoverage::default(),
            donors: vec![5],
            reads: None,
        }];
        let families = UnitShares {
            reads: vec![6, 7],
            group_of: vec![Some(0), Some(1)],
            // Read 6: family 1000 and 2000 (no row). Read 7: family 2000 only.
            posteriors: vec![vec![(1000, 0.5), (2000, 0.5)], vec![(2000, 1.0)]],
        };
        let classifications = Classifications {
            fragments: &fragments,
            group_of: &group_of,
            posteriors: &posteriors,
            genome_level: &genome_level,
            genus_members: &[(100, 2, 1.0)],
            rollup_rows: &rollup_rows,
            orphans: &[(5, 5)],
            // Read 8: credited to species 1 and 4.
            credited: &[(8, vec![(1, 0.25), (4, 0.75)])],
            families: Some(&families),
        };
        let rows = classifications
            .shares(&blocks(), &accessions())
            .into_iter()
            .map(|(read, shares)| classification_row(&format!("read{read}"), &shares))
            .collect::<Vec<_>>();
        assert_eq!(
            rows,
            [
                "read0\t1\tspecies\t1\t1",
                // The left-out species' share is dropped and the rest renormalised.
                "read1\t1\tspecies\t1\t1",
                "read2\t100\tgenus\t100\t1",
                // The species share split among the strain rows.
                "read3\t400\tstrain\t400;4;1\t0.67500;0.22500;0.10000",
                "read4\t4\tspecies\t4;1\t0.90000;0.10000",
                "read5\t100\tgenus\t100\t1",
                "read6\t1000\tfamily\t1000\t1",
                "read8\t4\tspecies\t4;1\t0.75000;0.25000",
            ]
        );
        // The genome row on genome 41, under the species' taxid, keeps its rank.
        let blocks = blocks();
        let shares = classifications.shares(&blocks, &accessions());
        let read3 = &shares.iter().find(|row| row.0 == 3).unwrap().1;
        assert_eq!((read3[1].0, read3[1].1), (4, "strain"));
        assert!((read3[1].2 - 0.225).abs() < 1e-12);
    }

    #[test]
    fn shares_on_one_row_are_summed_and_renormalised() {
        let shares = vec![
            (100, "genus", 0.25),
            (7, "species", 0.5),
            (100, "genus", 0.25),
            (8, "species", 0.0),
        ];
        let shares = renormalised(shares).unwrap();
        assert_eq!(shares.len(), 2);
        assert_eq!((shares[0].0, shares[0].1), (7, "species"));
        assert!((shares[0].2 - 0.5).abs() < 1e-12 && (shares[1].2 - 0.5).abs() < 1e-12);
        // Equal shares come by taxid.
        assert_eq!(shares[1].0, 100);
        assert!(renormalised(vec![(1, "species", 0.0)]).is_none());
        assert_eq!(
            classification_row("r", &shares),
            "r\t7\tspecies\t7;100\t0.50000;0.50000"
        );
    }

    #[test]
    fn a_fragment_counted_in_a_labelled_unit_s_row_carries_its_rank() {
        // Unit 100 is labelled at genus; read 0 is shared between it and species 1.
        let fragments = vec![fragment(0)];
        let posteriors = vec![vec![(100, 0.7), (1, 0.3)]];
        let blocks = vec![
            block(1, Some(row(1, "species", 10)), Vec::new()),
            TaxonBlock {
                unit: true,
                ..block(100, Some(row(100, "genus", 20)), Vec::new())
            },
        ];
        let classifications = Classifications {
            fragments: &fragments,
            group_of: &[Some(0)],
            posteriors: &posteriors,
            genome_level: &HashMap::new(),
            genus_members: &[],
            rollup_rows: &[],
            orphans: &[],
            credited: &[],
            families: None,
        };
        let rows = classifications
            .shares(&blocks, &accessions())
            .into_iter()
            .map(|(read, shares)| classification_row(&format!("read{read}"), &shares))
            .collect::<Vec<_>>();
        assert_eq!(rows, ["read0\t100\tgenus\t100;1\t0.70000;0.30000"]);
    }
}

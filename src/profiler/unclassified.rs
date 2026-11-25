//! Share of the input that the final profile does not explain.
//!
//! Every profile row explains the fragments in its read-count column, which
//! for a fractional genome is already scaled to the whole genome. A species
//! represented by marker genes counts only the fragments on its genes there;
//! its whole genome offers more placements, by the ratio of the placements of
//! a sequence of its declared genome size to those its genes offer. Those
//! placements count a fragment on every gene it lies on, so the ratio is
//! multiplied by the gene hits per fragment, which exceed 1 where the two
//! mates of a pair lie on two genes. Long reads differ in length and can
//! cover several genes, so for them the ratio is measured on the reads, from
//! their hits on the genes that need all the bases the gates ask of a read.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::Write;

use csv::Writer;

use super::markers::MarkerModel;
use super::output::{ProfileRow, PROFILE_COLUMNS};
use crate::db::ReferenceType;

/// Whole-genome fragments per fragment on marker genes, by species; `None` for
/// a species with a marker genome that declares no genome size or whose genes
/// cannot hold a read of the sample's length.
pub(crate) type MarkerScales = HashMap<u32, Option<f64>>;

/// The marker genomes of the species represented only by marker genomes, by
/// species: their assembly ids.
pub(crate) type MarkerSpecies = HashMap<u32, BTreeSet<u32>>;

/// The fragments resolved to a marker genome: for short reads their hits on
/// its genes, for long reads the whole-genome reads they stand for, `per_base
/// G + fixed` for a genome of `G` bases ([`MarkerModel::whole_genome_reads`]
/// summed over the reads).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct GenomeReads {
    reads: u64,
    gene_hits: u64,
    per_base: f64,
    fixed: f64,
}

impl GenomeReads {
    /// One short-read fragment with `gene_hits` hits on the genome's genes.
    pub(crate) fn add_fragment(&mut self, gene_hits: u32) {
        self.reads += 1;
        self.gene_hits += u64::from(gene_hits);
    }

    /// One long read with `gene_hits` hits of `(per_base, fixed)` each.
    pub(crate) fn add(&mut self, gene_hits: u32, (per_base, fixed): (f64, f64)) {
        let hits = f64::from(gene_hits);
        self.reads += 1;
        self.per_base += hits * per_base;
        self.fixed += hits * fixed;
    }

    #[cfg(test)]
    pub(crate) fn reads(&self) -> u64 {
        self.reads
    }

    /// The whole-genome reads of long reads on a genome of `size` bases.
    fn whole_genome(&self, size: usize) -> f64 {
        self.per_base * size as f64 + self.fixed
    }

    /// The whole-genome reads per read of long reads on a genome of `size`
    /// bases; `None` without a read.
    pub(crate) fn whole_genome_per_read(&self, size: usize) -> Option<f64> {
        (self.reads > 0).then(|| self.whole_genome(size) / self.reads as f64)
    }

    /// The whole-genome fragments per short-read fragment on a genome whose
    /// fragment `scale` is [`MarkerModel::genome_scale`]: the scale times the
    /// gene hits, over the fragments.
    pub(crate) fn short_read_scale(&self, scale: f64) -> f64 {
        scale * self.gene_hits as f64 / self.reads as f64
    }

    /// Gene hits per short-read fragment, 1 without fragments.
    fn hits_per_fragment(&self) -> f64 {
        if self.reads == 0 {
            1.0
        } else {
            self.gene_hits as f64 / self.reads as f64
        }
    }
}

/// One reference as the estimate needs it.
pub(crate) struct ReferenceKind {
    pub(crate) taxid: u32,
    pub(crate) assembly_id: u32,
    pub(crate) reference_type: ReferenceType,
}

/// The species represented only by marker genomes, with their marker
/// genomes. A species that also has genome references counts its fragments as
/// they are.
pub(crate) fn marker_species(references: impl IntoIterator<Item = ReferenceKind>) -> MarkerSpecies {
    let mut with_genomes = HashSet::new();
    let mut marker_genomes: MarkerSpecies = HashMap::new();
    for reference in references {
        match reference.reference_type {
            ReferenceType::Genome => {
                with_genomes.insert(reference.taxid);
            }
            ReferenceType::Markers => {
                marker_genomes
                    .entry(reference.taxid)
                    .or_default()
                    .insert(reference.assembly_id);
            }
        }
    }
    marker_genomes.retain(|taxid, _| !with_genomes.contains(taxid));
    marker_genomes
}

/// Scale of every marker species, at the genome sizes its marker genomes
/// declare to the marker `model`: averaged over its marker genomes, the
/// model's times the gene hits per fragment of the fragments `resolved` to the
/// genome, or for long reads the whole-genome reads of its resolved reads over
/// their number.
pub(crate) fn marker_scales(
    species: &MarkerSpecies,
    model: Option<&MarkerModel>,
    resolved: &HashMap<u32, GenomeReads>,
    long_reads: bool,
) -> MarkerScales {
    species
        .iter()
        .map(|(&taxid, genomes)| {
            let scale = || {
                let model = model?;
                if long_reads {
                    let (mut reads, mut whole_genome) = (0, 0.0);
                    for &assembly_id in genomes {
                        if let Some(genome) = resolved.get(&assembly_id) {
                            reads += genome.reads;
                            whole_genome += genome.whole_genome(model.declared_size(assembly_id)?);
                        }
                    }
                    (reads > 0).then(|| whole_genome / reads as f64)
                } else {
                    let mut sum = 0.0;
                    for &assembly_id in genomes {
                        let hits = resolved
                            .get(&assembly_id)
                            .map_or(1.0, GenomeReads::hits_per_fragment);
                        sum += model.declared_genome_scale(assembly_id)? * hits;
                    }
                    Some(sum / genomes.len() as f64)
                }
            };
            (taxid, scale())
        })
        .collect()
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct UnclassifiedEstimate {
    /// Input fragments that no row explains.
    pub(crate) fragments: f64,
    /// Their share of the input, between 0 and 1.
    pub(crate) share: f64,
}

/// The estimate for `rows`, or the first species (by taxid) without a marker
/// scale.
pub(crate) fn estimate(
    rows: &[ProfileRow],
    scales: &MarkerScales,
    input_fragments: usize,
) -> Result<UnclassifiedEstimate, u32> {
    let mut explained = 0.0;
    let mut missing: Option<u32> = None;
    for row in rows {
        match scales.get(&row.taxid) {
            None => explained += row.equivalent_mapcount,
            Some(Some(scale)) => explained += row.equivalent_mapcount * scale,
            Some(None) => missing = Some(missing.map_or(row.taxid, |taxid| taxid.min(row.taxid))),
        }
    }
    if let Some(taxid) = missing {
        return Err(taxid);
    }
    let input = input_fragments as f64;
    let explained_share = if input > 0.0 {
        (explained / input).min(1.0)
    } else {
        1.0
    };
    Ok(UnclassifiedEstimate {
        fragments: (input - explained).max(0.0),
        share: 1.0 - explained_share,
    })
}

/// The `unclassified` profile row; NA throughout when the estimate is
/// unavailable, and NA in the columns of a row's coverage, best-fit genome and
/// identity.
pub(crate) fn write_row<W: Write>(
    writer: &mut Writer<W>,
    estimate: Option<UnclassifiedEstimate>,
) -> csv::Result<()> {
    let mut fields = vec![
        "unclassified".to_string(),
        "0".to_string(),
        "unclassified".to_string(),
    ];
    let (fragments, share) = match estimate {
        Some(estimate) => (
            format!("{:.5}", estimate.fragments),
            format!("{:.5}", estimate.share * 100.0),
        ),
        None => ("NA".to_string(), "NA".to_string()),
    };
    fields.push(fragments);
    fields.push("NA".to_string());
    fields.push(share);
    fields.resize(PROFILE_COLUMNS.len(), "NA".to_string());
    writer.write_record(fields)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profiler::markers::{FragmentGeometry, FragmentLayout};
    use crate::profiler::output::BestFitReport;
    use std::sync::Arc;

    fn row(taxid: u32, fragments: f64) -> ProfileRow {
        ProfileRow {
            taxname: taxid.to_string(),
            taxid,
            rank: "species".to_string(),
            equivalent_mapcount: fragments,
            depth: 0.0,
            abundance: 0.0,
            breadth: 0.0,
            expected_breadth: 0.0,
            breadth_pvalue: None,
            fixed_chunk: 0.0,
            flex_chunk: 0.0,
            expected_flex_chunk: 0.0,
            chunk_pvalue: None,
            best_fit: BestFitReport {
                assembly_id: 1,
                genome: Arc::from("GCF_1"),
                evidence_depth: 0.0,
                identity: (f64::NAN, f64::NAN),
            },
            annotations: None,
        }
    }

    #[test]
    fn marker_species_count_as_their_whole_genome() {
        // Species 1 has a genome reference, species 2 marker genes only.
        let scales = MarkerScales::from([(2, Some(50.0))]);
        let rows = [row(1, 600.0), row(2, 4.0)];
        let estimate = estimate(&rows, &scales, 1_000).unwrap();
        assert!((estimate.fragments - 200.0).abs() < 1e-9);
        assert!((estimate.share - 0.2).abs() < 1e-12);
        // More explained than the input: nothing is unclassified.
        let full = super::estimate(&[row(1, 1_200.0)], &scales, 1_000).unwrap();
        assert_eq!((full.fragments, full.share), (0.0, 0.0));
    }

    #[test]
    fn a_marker_species_without_a_genome_size_makes_the_estimate_unavailable() {
        let scales = MarkerScales::from([(2, None), (3, None)]);
        assert_eq!(
            estimate(&[row(3, 1.0), row(2, 1.0), row(1, 1.0)], &scales, 10),
            Err(2)
        );
    }

    #[test]
    fn genome_references_decide_and_marker_genomes_are_averaged() {
        let kind = |taxid, assembly_id, reference_type| ReferenceKind {
            taxid,
            assembly_id,
            reference_type,
        };
        let references = [
            kind(1, 0, ReferenceType::Genome),
            kind(1, 1, ReferenceType::Markers),
            kind(2, 2, ReferenceType::Markers),
            kind(3, 3, ReferenceType::Markers),
        ];
        let species = marker_species(references);
        assert_eq!(
            species,
            MarkerSpecies::from([(2, BTreeSet::from([2])), (3, BTreeSet::from([3]))])
        );
        // Without a marker model no scale can be computed.
        let scales = marker_scales(&species, None, &HashMap::new(), false);
        assert_eq!(scales.get(&1), None);
        assert_eq!(scales.get(&2), Some(&None));
        assert_eq!(scales.get(&3), Some(&None));
    }

    #[test]
    fn pairs_on_two_genes_raise_the_scale_by_their_gene_hits() {
        let geometry = FragmentGeometry::new(FragmentLayout::PairsAndLoneMates, 150, 400, 50, 0.7);
        let model =
            MarkerModel::from_gene_lengths(&[(1, &[1_000, 1_000]), (2, &[1_000])], geometry)
                .with_declared(HashMap::from([(1, 1_000_000), (2, 500_000)]));
        let species = MarkerSpecies::from([(7, BTreeSet::from([1, 2]))]);
        // Genome 1: four fragments, one with its two mates on the two genes;
        // none resolved to genome 2.
        let mut resolved = HashMap::new();
        for hits in [1, 2, 1, 1] {
            resolved
                .entry(1)
                .or_insert_with(GenomeReads::default)
                .add_fragment(hits);
        }
        let scales = marker_scales(&species, Some(&model), &resolved, false);
        let expected = (model.genome_scale(1, 1_000_000).unwrap() * 1.25
            + model.genome_scale(2, 500_000).unwrap())
            / 2.0;
        assert_eq!(scales[&7], Some(expected));
    }

    #[test]
    fn long_reads_measure_the_scale_of_a_marker_species() {
        // Genome 3 declares no size.
        let geometry = FragmentGeometry::new(FragmentLayout::LongReads, 3_000, 3_000, 500, 0.2);
        let model = MarkerModel::from_gene_lengths(&[], geometry).with_declared(HashMap::from([
            (2, 1_000),
            (4, 3_000),
            (5, 1_000),
        ]));
        let species = MarkerSpecies::from([
            (2, BTreeSet::from([2, 4])),
            (3, BTreeSet::from([3])),
            (5, BTreeSet::from([5])),
        ]);
        let mut reads = HashMap::new();
        // Genome 2: one read hitting two genes, each standing for 0.001 G + 1.
        reads
            .entry(2)
            .or_insert_with(GenomeReads::default)
            .add(2, (0.001, 1.0));
        // Genome 4: one read with one hit of 0.002 G.
        reads
            .entry(4)
            .or_insert_with(GenomeReads::default)
            .add(1, (0.002, 0.0));
        reads
            .entry(3)
            .or_insert_with(GenomeReads::default)
            .add(1, (0.001, 0.0));
        let scales = marker_scales(&species, Some(&model), &reads, true);
        // (2 x (1 + 1) + 6) whole-genome reads over 2 reads.
        assert_eq!(scales[&2], Some(5.0));
        // No declared size, or no read to measure on.
        assert_eq!(scales[&3], None);
        assert_eq!(scales[&5], None);
    }

    #[test]
    fn the_row_shows_fragments_and_share_or_na() {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());
        write_row(
            &mut writer,
            Some(UnclassifiedEstimate {
                fragments: 200.0,
                share: 0.2,
            }),
        )
        .unwrap();
        write_row(&mut writer, None).unwrap();
        let text = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        let na = ["NA"; 11].join("\t");
        assert_eq!(
            text,
            format!(
                "unclassified\t0\tunclassified\t200.00000\tNA\t20.00000\t{na}\n\
                 unclassified\t0\tunclassified\tNA\tNA\tNA\t{na}\n"
            )
        );
        assert!(text.lines().all(|line| line.split('\t').count() == 17));
    }
}

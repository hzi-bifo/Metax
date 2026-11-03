//! Share of the input that the final profile does not explain.
//!
//! Every profile row explains the fragments in its read-count column, which
//! for a fractional genome is already scaled to the whole genome. A species
//! represented by marker genes counts only the fragments on its genes there;
//! its whole genome offers more placements, by the ratio of the placements of
//! a sequence of its declared genome size to those its genes offer.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Write;

use csv::Writer;

use super::markers::MarkerModel;
use super::output::ProfileRow;
use crate::db::ReferenceType;

/// Whole-genome fragments per fragment on marker genes, by species; `None` for
/// a species with a marker genome that declares no genome size or whose genes
/// cannot hold a read of the sample's length.
pub(crate) type MarkerScales = HashMap<u32, Option<f64>>;

/// One reference as the estimate needs it.
pub(crate) struct ReferenceKind {
    pub(crate) taxid: u32,
    pub(crate) assembly_id: u32,
    pub(crate) reference_type: ReferenceType,
    pub(crate) declared_genome_size: Option<usize>,
}

/// Scales of the species represented only by marker genomes, averaged over
/// their marker genomes. A species that also has genome references counts its
/// fragments as they are.
pub(crate) fn marker_scales(
    references: impl IntoIterator<Item = ReferenceKind>,
    model: Option<&MarkerModel>,
) -> MarkerScales {
    let mut with_genomes = HashSet::new();
    let mut marker_genomes: HashMap<u32, BTreeMap<u32, Option<usize>>> = HashMap::new();
    for reference in references {
        match reference.reference_type {
            ReferenceType::Genome => {
                with_genomes.insert(reference.taxid);
            }
            ReferenceType::Markers => {
                marker_genomes
                    .entry(reference.taxid)
                    .or_default()
                    .insert(reference.assembly_id, reference.declared_genome_size);
            }
        }
    }
    marker_genomes
        .into_iter()
        .filter(|(taxid, _)| !with_genomes.contains(taxid))
        .map(|(taxid, genomes)| {
            let mean = || {
                let model = model?;
                let mut sum = 0.0;
                for (&assembly_id, &size) in &genomes {
                    sum += model.genome_scale(assembly_id, size?)?;
                }
                Some(sum / genomes.len() as f64)
            };
            (taxid, mean())
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
/// unavailable.
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
    fields.extend(std::iter::repeat_n("NA".to_string(), 7));
    writer.write_record(fields)
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let kind = |taxid, assembly_id, reference_type, size| ReferenceKind {
            taxid,
            assembly_id,
            reference_type,
            declared_genome_size: size,
        };
        let references = [
            kind(1, 0, ReferenceType::Genome, None),
            kind(1, 1, ReferenceType::Markers, Some(1_000)),
            kind(2, 2, ReferenceType::Markers, Some(1_000)),
            kind(3, 3, ReferenceType::Markers, None),
        ];
        // Without a marker model no scale can be computed.
        let scales = marker_scales(references, None);
        assert_eq!(scales.get(&1), None);
        assert_eq!(scales.get(&2), Some(&None));
        assert_eq!(scales.get(&3), Some(&None));
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
        let na = ["NA"; 7].join("\t");
        assert_eq!(
            text,
            format!(
                "unclassified\t0\tunclassified\t200.00000\tNA\t20.00000\t{na}\n\
                 unclassified\t0\tunclassified\tNA\tNA\tNA\t{na}\n"
            )
        );
    }
}

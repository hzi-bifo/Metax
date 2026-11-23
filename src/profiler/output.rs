use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;

use anyhow::{Context, Result};
use csv::Writer;
use flate2::write::GzEncoder;
use flate2::Compression;

/// Optional host/pathogen annotation columns appended to a profile row.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ProfileAnnotations {
    pub(crate) host_names: String,
    pub(crate) host_taxids: String,
    pub(crate) diseases: String,
}

/// Columns of the profile. Every row has all of them, the unclassified row
/// included.
pub(crate) const PROFILE_COLUMNS: [&str; 17] = [
    "name",
    "taxid",
    "rank",
    "reads",
    "depth",
    "abundance",
    "breadth",
    "expected_breadth",
    "breadth_pvalue",
    "fixed_chunk_breadth",
    "flex_chunk_breadth",
    "expected_flex_chunk_breadth",
    "chunk_pvalue",
    "evidence_depth",
    "best_fit_genome",
    "mean_identity",
    "estimated_ani",
];

/// Columns that `--host` appends to every row of the profile.
pub(crate) const HOST_COLUMNS: [&str; 3] = ["host_names", "host_taxids", "diseases"];

/// The genome a taxon kept by the coverage evidence is reported on, and the
/// depth of that evidence there before abundance estimation.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BestFitReport {
    pub(crate) genome: Arc<str>,
    pub(crate) evidence_depth: f64,
    /// Mean identity of the row's fragments to the genome and the estimated
    /// ANI of the sample's strain to it; NaN, written `NA`, for a value that
    /// is not known.
    pub(crate) identity: (f64, f64),
}

/// Typed profile row retained at full precision until the final write.
///
/// Filtering and abundance normalization must use these numeric fields, not
/// their five-decimal display representation. This prevents values near a
/// cutoff from crossing it merely because they were formatted and reparsed.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ProfileRow {
    pub(crate) taxname: String,
    pub(crate) taxid: u32,
    pub(crate) rank: String,
    pub(crate) equivalent_mapcount: f64,
    pub(crate) depth: f64,
    pub(crate) abundance: f64,
    pub(crate) breadth: f64,
    pub(crate) expected_breadth: f64,
    pub(crate) breadth_pvalue: Option<f64>,
    pub(crate) fixed_chunk: f64,
    pub(crate) flex_chunk: f64,
    pub(crate) expected_flex_chunk: f64,
    pub(crate) chunk_pvalue: Option<f64>,
    pub(crate) best_fit: BestFitReport,
    pub(crate) annotations: Option<ProfileAnnotations>,
}

/// Gzip writer of the per-read classifications, `<prefix>.classify.txt.gz`.
pub(crate) type ClassificationWriter = BufWriter<GzEncoder<File>>;

/// Create the per-read classification file at `path`. Gzip level 2: the
/// backend's level 1 is a quick mode that compresses much less for little
/// gain in speed.
pub(crate) fn create_classification_writer(path: &str) -> Result<ClassificationWriter> {
    let file = File::create(path).with_context(|| format!("failed to create {path}"))?;
    Ok(BufWriter::new(GzEncoder::new(file, Compression::new(2))))
}

/// Write the gzip trailer. Dropping the writer would also close the stream,
/// but would hide a failed write.
pub(crate) fn finish_classification_writer(writer: ClassificationWriter) -> Result<()> {
    writer
        .into_inner()
        .map_err(|error| error.into_error())?
        .finish()?;
    Ok(())
}

impl ProfileRow {
    pub(crate) fn write<W: Write>(&self, writer: &mut Writer<W>) -> csv::Result<()> {
        let mut fields = Vec::with_capacity(
            PROFILE_COLUMNS.len() + self.annotations.as_ref().map_or(0, |_| HOST_COLUMNS.len()),
        );
        fields.push(self.taxname.clone());
        fields.push(self.taxid.to_string());
        fields.push(self.rank.clone());
        fields.push(format_five(self.equivalent_mapcount));
        fields.push(format_five(self.depth));
        fields.push(format_five(self.abundance));
        fields.push(format_five(self.breadth));
        fields.push(format_five(self.expected_breadth));
        fields.push(format_probability(self.breadth_pvalue));
        fields.push(format_five(self.fixed_chunk));
        fields.push(format_five(self.flex_chunk));
        fields.push(format_five(self.expected_flex_chunk));
        fields.push(format_probability(self.chunk_pvalue));
        fields.push(format_five(self.best_fit.evidence_depth));
        fields.push(self.best_fit.genome.to_string());
        fields.push(format_known(self.best_fit.identity.0));
        fields.push(format_known(self.best_fit.identity.1));
        if let Some(annotations) = &self.annotations {
            fields.push(annotations.host_names.clone());
            fields.push(annotations.host_taxids.clone());
            fields.push(annotations.diseases.clone());
        }
        writer.write_record(fields)
    }
}

#[inline]
pub(crate) fn format_probability(value: Option<f64>) -> String {
    match value {
        Some(value) if value.is_finite() => {
            let adjusted = if value.abs() < 1e-12 { 0.0 } else { value };
            if adjusted.abs() > 0.0 && adjusted.abs() < 1e-4 {
                format!("{adjusted:.6e}")
            } else {
                format_five(adjusted)
            }
        }
        _ => "NA".to_string(),
    }
}

#[inline]
pub(crate) fn format_five(value: f64) -> String {
    format!("{value:.5}")
}

/// A value that may be unknown (not finite): `NA` then.
fn format_known(value: f64) -> String {
    if value.is_finite() {
        format_five(value)
    } else {
        "NA".to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::{BestFitReport, ProfileAnnotations, ProfileRow, HOST_COLUMNS, PROFILE_COLUMNS};
    use std::sync::Arc;

    fn row() -> ProfileRow {
        ProfileRow {
            taxname: "Species A".to_string(),
            taxid: 101,
            rank: "species".to_string(),
            equivalent_mapcount: 3.0,
            depth: 0.2,
            abundance: 0.25,
            breadth: 0.1,
            expected_breadth: 0.11,
            breadth_pvalue: Some(1.0),
            fixed_chunk: 0.2,
            flex_chunk: 0.12,
            expected_flex_chunk: 0.13,
            chunk_pvalue: None,
            best_fit: BestFitReport {
                genome: Arc::from("GCF_000000001.1"),
                evidence_depth: 1.5,
                identity: (f64::NAN, 0.99),
            },
            annotations: None,
        }
    }

    #[test]
    fn writer_formats_only_at_the_output_boundary() {
        let mut row = row();
        row.annotations = Some(ProfileAnnotations {
            host_names: "host".to_string(),
            host_taxids: "7".to_string(),
            diseases: "condition".to_string(),
        });
        let mut writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(Vec::new());
        row.write(&mut writer).unwrap();
        let bytes = writer.into_inner().unwrap();
        let text = String::from_utf8(bytes).unwrap();
        let fields = text.trim_end().split('\t').collect::<Vec<_>>();
        assert_eq!(fields.len(), PROFILE_COLUMNS.len() + HOST_COLUMNS.len());
        assert_eq!(fields[5], "0.25000");
        assert_eq!(fields[8], "1.00000");
        assert_eq!(fields[12], "NA");
    }

    #[test]
    fn writer_reports_the_best_fit_genome_before_host_columns() {
        let write = |row: &ProfileRow| {
            let mut writer = csv::WriterBuilder::new()
                .has_headers(false)
                .delimiter(b'\t')
                .from_writer(Vec::new());
            row.write(&mut writer).unwrap();
            let text = String::from_utf8(writer.into_inner().unwrap()).unwrap();
            text.trim_end()
                .split('\t')
                .map(str::to_string)
                .collect::<Vec<_>>()
        };
        let mut row = row();
        assert_eq!(
            write(&row)[13..],
            ["1.50000", "GCF_000000001.1", "NA", "0.99000"]
        );
        row.annotations = Some(ProfileAnnotations {
            host_names: "host".to_string(),
            host_taxids: "7".to_string(),
            diseases: "condition".to_string(),
        });
        let fields = write(&row);
        assert_eq!(
            fields[13..],
            [
                "1.50000",
                "GCF_000000001.1",
                "NA",
                "0.99000",
                "host",
                "7",
                "condition"
            ]
        );
    }

    #[test]
    fn writer_preserves_scientific_notation_for_small_probabilities() {
        let mut row = row();
        row.breadth_pvalue = Some(1e-5);
        let mut writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(Vec::new());
        row.write(&mut writer).unwrap();
        let text = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        assert_eq!(text.trim_end().split('\t').nth(8), Some("1.000000e-5"));
    }
}

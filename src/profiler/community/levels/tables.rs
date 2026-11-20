//! Tables of the level-wise engine: the lines that open its profile, two lines
//! starting with `#`: the run's metadata, in the style of the mapping file's
//! header, then the column names.

use crate::db::DatabaseMetadata;
use crate::profiler::ProfilerConfig;
use crate::sequencer::Sequencer;

/// Columns of the level-wise engine's profile. Every row has all of them, the
/// unclassified row included.
pub(super) const PROFILE_COLUMNS: [&str; 17] = [
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

/// The run named on the first line of each table.
pub(super) struct TableHeader {
    /// Fingerprint of the database, as in the mapping file's header; `None`
    /// without database metadata.
    pub(super) fingerprint: Option<u128>,
    pub(super) sequencer: Sequencer,
    pub(super) paired: bool,
}

impl TableHeader {
    pub(super) fn of(cfg: &ProfilerConfig) -> Self {
        Self {
            fingerprint: cfg
                .db_metadata
                .as_deref()
                .map(DatabaseMetadata::fingerprint),
            sequencer: cfg.sequencer,
            paired: cfg.is_paired,
        }
    }

    /// The two header lines of table `table`: `#metax_<table>` with the run's
    /// metadata, then `#` and the `columns`.
    pub(super) fn lines(&self, table: &str, columns: &[&str]) -> String {
        let fingerprint = self
            .fingerprint
            .map_or_else(|| "NA".to_string(), |f| format!("{f:032x}"));
        format!(
            "#metax_{table}\tversion={}\tengine=level\tdatabase_fingerprint={fingerprint}\tsequencer={}\tpaired={}\n#{}\n",
            env!("CARGO_PKG_VERSION"),
            self.sequencer,
            self.paired,
            columns.join("\t")
        )
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::super::profile_row;
    use super::*;
    use crate::profiler::evidence::EvidenceCoverage;
    use crate::profiler::unclassified;

    fn header() -> TableHeader {
        TableHeader {
            fingerprint: None,
            sequencer: Sequencer::Illumina,
            paired: true,
        }
    }

    #[test]
    fn the_profile_opens_with_its_run_and_the_name_of_every_column() {
        let text = header().lines("profile", &PROFILE_COLUMNS);
        let lines = text.lines().collect::<Vec<_>>();
        assert_eq!(lines.len(), 2);
        assert_eq!(
            lines[0],
            format!(
                "#metax_profile\tversion={}\tengine=level\tdatabase_fingerprint=NA\tsequencer=Illumina\tpaired=true",
                env!("CARGO_PKG_VERSION")
            )
        );
        assert!(lines[1].starts_with("#name\ttaxid\trank\treads\tdepth\tabundance\t"));
        assert!(
            lines[1].ends_with("\tevidence_depth\tbest_fit_genome\tmean_identity\testimated_ani")
        );
        let with_database = TableHeader {
            fingerprint: Some(0xabc),
            ..header()
        };
        assert!(with_database
            .lines("profile", &PROFILE_COLUMNS)
            .contains("\tdatabase_fingerprint=00000000000000000000000000000abc\t"));
        // Every row has a column per name: a species row and the unclassified row, also when no row precedes it.
        let mut writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(Vec::new());
        profile_row(
            "Species".to_string(),
            1,
            "species",
            (5.0, 1.0),
            &EvidenceCoverage::default(),
            Arc::from("GCF_1"),
            Some((f64::NAN, f64::NAN)),
        )
        .write(&mut writer)
        .unwrap();
        unclassified::write_row(&mut writer, None, true).unwrap();
        let text = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        let rows = text.lines().collect::<Vec<_>>();
        assert_eq!(rows.len(), 2);
        assert!(rows
            .iter()
            .all(|row| row.split('\t').count() == PROFILE_COLUMNS.len()));
    }
}

use std::collections::HashMap;
use std::io::Write;

use csv::Writer;

/// Optional host/pathogen annotation columns appended to a profile row.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ProfileAnnotations {
    pub(crate) host_names: String,
    pub(crate) host_taxids: String,
    pub(crate) diseases: String,
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
    pub(crate) annotations: Option<ProfileAnnotations>,
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum ClassificationSource {
    Direct,
    Candidates(Vec<u32>),
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct ResolvedClassificationKey {
    pub(crate) taxid: u32,
    pub(crate) source: ClassificationSource,
}

#[derive(Debug)]
pub(crate) struct ResolvedClassificationGroup {
    pub(crate) key: ResolvedClassificationKey,
    pub(crate) count: usize,
    pub(crate) reads: Vec<String>,
}

impl ResolvedClassificationGroup {
    /// Borrow read names in deterministic lexical order without cloning their
    /// strings. Called only when optional classification output is requested.
    pub(crate) fn sorted_reads(&self) -> Vec<&str> {
        let mut reads = self.reads.iter().map(String::as_str).collect::<Vec<_>>();
        reads.sort_unstable();
        reads
    }
}

/// Grouped classifications retained until final taxon filtering.
///
/// Read names are held only when classification output was explicitly
/// requested. Counts are always kept so progress summaries can describe the
/// same retained assignments as `profile.txt` without imposing per-read
/// memory on normal runs.
#[derive(Debug, Default)]
pub(crate) struct ResolvedClassifications {
    indices: HashMap<ResolvedClassificationKey, usize>,
    pub(crate) groups: Vec<ResolvedClassificationGroup>,
    retain_read_names: bool,
}

impl ResolvedClassifications {
    pub(crate) fn new(retain_read_names: bool) -> Self {
        Self {
            retain_read_names,
            ..Self::default()
        }
    }

    /// Borrow classification groups in a stable output order. The temporary
    /// vector contains pointers only and is allocated solely for `classify.txt`.
    pub(crate) fn sorted_groups(&self) -> Vec<&ResolvedClassificationGroup> {
        let mut groups = self.groups.iter().collect::<Vec<_>>();
        groups.sort_unstable_by(|left, right| left.key.cmp(&right.key));
        groups
    }

    pub(crate) fn record_read(&mut self, taxid: u32, read: &str) {
        let reads = self.retain_read_names.then(|| vec![read.to_string()]);
        self.record(taxid, ClassificationSource::Direct, 1, reads);
    }

    pub(crate) fn record_direct_group(
        &mut self,
        taxid: u32,
        count: usize,
        reads: Option<Vec<String>>,
    ) {
        self.record(taxid, ClassificationSource::Direct, count, reads);
    }

    pub(crate) fn record_group(
        &mut self,
        taxid: u32,
        source_taxids: Vec<u32>,
        count: usize,
        reads: Option<Vec<String>>,
    ) {
        self.record(
            taxid,
            ClassificationSource::Candidates(source_taxids),
            count,
            reads,
        );
    }

    fn record(
        &mut self,
        taxid: u32,
        source: ClassificationSource,
        count: usize,
        reads: Option<Vec<String>>,
    ) {
        if count == 0 {
            return;
        }
        let key = ResolvedClassificationKey { taxid, source };
        if let Some(&index) = self.indices.get(&key) {
            let group = &mut self.groups[index];
            group.count = group.count.saturating_add(count);
            if self.retain_read_names {
                if let Some(mut reads) = reads {
                    group.reads.append(&mut reads);
                }
            }
            return;
        }
        let index = self.groups.len();
        self.indices.insert(key.clone(), index);
        self.groups.push(ResolvedClassificationGroup {
            key,
            count,
            reads: if self.retain_read_names {
                reads.unwrap_or_default()
            } else {
                Vec::new()
            },
        });
    }
}

impl ProfileRow {
    pub(crate) fn normalize_abundance(&mut self, total_abundance: f64) {
        self.abundance = if total_abundance > 0.0 && total_abundance.is_finite() {
            self.abundance / total_abundance * 100.0
        } else {
            0.0
        };
    }

    pub(crate) fn write<W: Write>(&self, writer: &mut Writer<W>) -> csv::Result<()> {
        let mut fields = Vec::with_capacity(if self.annotations.is_some() { 16 } else { 13 });
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
fn format_five(value: f64) -> String {
    format!("{value:.5}")
}

#[cfg(test)]
mod tests {
    use super::{
        ClassificationSource, ProfileAnnotations, ProfileRow, ResolvedClassifications,
    };

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
            annotations: None,
        }
    }

    #[test]
    fn normalization_uses_unrounded_abundance() {
        let mut row = row();
        row.normalize_abundance(0.5);
        assert_eq!(row.abundance, 50.0);
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
        assert_eq!(fields.len(), 16);
        assert_eq!(fields[5], "0.25000");
        assert_eq!(fields[8], "1.00000");
        assert_eq!(fields[12], "NA");
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

    #[test]
    fn count_only_classifications_do_not_retain_read_names() {
        let mut classifications = ResolvedClassifications::new(false);
        classifications.record_read(101, "read-one");
        classifications.record_read(101, "read-two");
        assert_eq!(classifications.groups.len(), 1);
        assert_eq!(classifications.groups[0].count, 2);
        assert!(classifications.groups[0].reads.is_empty());
    }

    #[test]
    fn classification_views_are_sorted_without_mutating_storage() {
        let mut classifications = ResolvedClassifications::new(true);
        classifications.record(
            202,
            ClassificationSource::Candidates(vec![202, 303]),
            1,
            Some(vec!["read-z".to_string(), "read-a".to_string()]),
        );
        classifications.record_read(101, "read-c");

        let groups = classifications.sorted_groups();
        assert_eq!(groups[0].key.taxid, 101);
        assert_eq!(groups[1].key.taxid, 202);
        assert_eq!(groups[1].sorted_reads(), vec!["read-a", "read-z"]);
        assert_eq!(classifications.groups[0].reads[0], "read-z");
    }
}

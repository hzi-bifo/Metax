//! Coverage evidence for genomes represented by marker genes.
//!
//! Every reference sequence of a marker genome is one gene. A fragment counts
//! on a gene only when enough of a read lies inside it, so a short gene offers
//! fewer placements than its length suggests and bases near gene ends are
//! covered less often than bases in the middle. The genome model of
//! [`calculate_coverage`] treats a genome as one long sequence cut into equal
//! chunks; for a marker genome the genes are the chunks, each weighted by the
//! placements it offers, so that the expected values match how reads can
//! actually land:
//!
//! - amount: the fraction of marker genes hit (reported as both chunk
//!   breadths), and the number of genes hit for `--min-marker-genes`;
//! - evenness: genes hit against the genes expected to be hit by the same
//!   fragments, with an exact lower-tail p-value (reported as the expected
//!   chunk breadth, chunk ratio and chunk p-value);
//! - base breadth: covered marker bases against the expected covered bases,
//!   which accounts for the lower coverage of bases near gene ends.
//!
//! [`calculate_coverage`]: super::evidence::calculate_coverage

use std::collections::HashMap;

use super::evidence::{
    adjusted_deviation_pvalue, merged_interval_bases, EvidenceCoverage, PresenceAssembly,
};
use crate::db::{DatabaseMetadata, ReferenceType};

/// Marker genes of the marker genomes and how this sample's fragments land on them.
#[derive(Debug)]
pub(crate) struct MarkerModel {
    /// Genes of every marker genome by assembly id, in reference order.
    genes: HashMap<u32, Vec<Gene>>,
    geometry: FragmentGeometry,
}

#[derive(Clone, Copy, Debug)]
struct Gene {
    reference_id: u32,
    len: u32,
}

impl MarkerModel {
    /// The marker genes of `db`, or `None` when it has no marker genomes.
    /// `assembly_id` gives the profiler's id of a genome.
    pub(crate) fn from_database(
        db: &DatabaseMetadata,
        assembly_id: impl Fn(&str) -> u32,
        geometry: FragmentGeometry,
    ) -> Option<Self> {
        let mut genes: HashMap<u32, Vec<Gene>> = HashMap::new();
        for (reference_id, reference) in db.references.iter().enumerate() {
            if reference.reference_type == ReferenceType::Markers {
                genes
                    .entry(assembly_id(&reference.genome))
                    .or_default()
                    .push(Gene {
                        reference_id: u32::try_from(reference_id).ok()?,
                        len: u32::try_from(reference.length).ok()?,
                    });
            }
        }
        (!genes.is_empty()).then_some(Self { genes, geometry })
    }

    pub(crate) fn is_marker_genome(&self, assembly_id: u32) -> bool {
        self.genes.contains_key(&assembly_id)
    }

    /// Coverage of a marker genome, or `None` for a genome of another type or
    /// when no gene can hold a read of this length.
    pub(crate) fn coverage(
        &self,
        assembly_id: u32,
        assembly: &mut PresenceAssembly,
        tests: usize,
    ) -> Option<EvidenceCoverage> {
        let genes = self.genes.get(&assembly_id)?;
        marker_coverage(genes, &self.geometry, assembly, tests)
    }
}

/// Which fragments the profiler counts on a gene.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum FragmentLayout {
    /// Single-end or long reads.
    Reads,
    /// Read pairs with both mates on the gene.
    Pairs,
    /// Read pairs, and lone mates whose mate lies beyond the gene end.
    PairsAndLoneMates,
}

/// Read and fragment lengths of the sample, and the bases a read needs inside
/// a gene to count.
#[derive(Clone, Copy, Debug)]
pub(crate) struct FragmentGeometry {
    layout: FragmentLayout,
    read_len: i64,
    template_len: i64,
    min_inside: i64,
}

impl FragmentGeometry {
    /// `min_mapped_len` and `min_fraction` are the profiler's read-level gates.
    pub(crate) fn new(
        layout: FragmentLayout,
        read_len: usize,
        template_len: usize,
        min_mapped_len: usize,
        min_fraction: f64,
    ) -> Self {
        let read_len = i64::try_from(read_len.max(1)).unwrap_or(i64::MAX);
        let by_fraction = (min_fraction.clamp(0.0, 1.0) * read_len as f64).ceil() as i64;
        let min_inside = i64::try_from(min_mapped_len)
            .unwrap_or(i64::MAX)
            .max(by_fraction)
            .clamp(1, read_len);
        let template_len = i64::try_from(template_len)
            .unwrap_or(i64::MAX)
            .max(read_len);
        Self {
            layout,
            read_len,
            template_len,
            min_inside,
        }
    }

    /// Start positions of a read with at least `min_inside` bases on a gene.
    fn read_starts(&self, len: i64) -> i64 {
        if len < self.min_inside {
            return 0;
        }
        len - 2 * self.min_inside + self.read_len + 1
    }

    /// Fragment placements that count on a gene of `len` bases.
    fn placements(&self, len: i64) -> f64 {
        let reads = self.read_starts(len);
        let placements = match self.layout {
            FragmentLayout::Reads => reads,
            FragmentLayout::Pairs => (len - self.template_len + 1).max(0),
            // Either mate places the fragment: the start ranges of the two
            // mates, offset by template minus read length, joined.
            FragmentLayout::PairsAndLoneMates => {
                (reads + self.template_len - self.read_len).min(2 * reads)
            }
        };
        placements as f64
    }

    /// Gene coordinates `[first, last]` of the start of the reported read.
    fn reported_read_starts(&self, len: i64) -> (i64, i64) {
        match self.layout {
            // Both mates lie on the gene, so the reported mate does too.
            FragmentLayout::Pairs => (0, len - self.read_len),
            _ => (self.min_inside - self.read_len, len - self.min_inside),
        }
    }

    /// Expected covered bases of a gene of `len` bases that receives each of
    /// `fragments` fragments with probability `share`, each reporting one read
    /// placed uniformly over its possible starts.
    fn expected_covered_bases(&self, len: i64, share: f64, fragments: f64) -> f64 {
        let (first, last) = self.reported_read_starts(len);
        let starts = last - first + 1;
        if starts <= 0 || share <= 0.0 {
            return 0.0;
        }
        // A base covered by `c` of the possible starts is covered by one
        // fragment with probability share * c / starts.
        let covered: Vec<f64> = (0..=self.read_len)
            .map(|c| 1.0 - (1.0 - share * c as f64 / starts as f64).powf(fragments))
            .collect();
        (0..len)
            .map(|base| {
                let c = (base.min(last) - (base - self.read_len + 1).max(first) + 1)
                    .clamp(0, self.read_len);
                covered[c as usize]
            })
            .sum()
    }
}

impl std::fmt::Display for FragmentGeometry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let layout = match self.layout {
            FragmentLayout::Reads => "reads",
            FragmentLayout::Pairs => "read pairs",
            FragmentLayout::PairsAndLoneMates => "read pairs and lone mates",
        };
        write!(
            f,
            "{layout}, read length {}, template length {}, at least {} bases of a read inside a gene",
            self.read_len, self.template_len, self.min_inside
        )
    }
}

fn marker_coverage(
    genes: &[Gene],
    geometry: &FragmentGeometry,
    assembly: &mut PresenceAssembly,
    tests: usize,
) -> Option<EvidenceCoverage> {
    let placements: Vec<f64> = genes
        .iter()
        .map(|gene| geometry.placements(i64::from(gene.len)))
        .collect();
    let total_placements: f64 = placements.iter().sum();
    if total_placements <= 0.0 {
        return None;
    }
    let fragments = assembly.count as f64;
    let marker_bases: f64 = genes.iter().map(|gene| f64::from(gene.len)).sum();

    let mut hit = vec![false; genes.len()];
    for interval in &assembly.intervals {
        if let Ok(index) = genes.binary_search_by_key(&interval.reference_id, |g| g.reference_id) {
            hit[index] = true;
        }
    }
    let genes_hit = hit.iter().filter(|&&hit| hit).count();
    let hit_probabilities: Vec<f64> = placements
        .iter()
        .map(|placements| 1.0 - (1.0 - placements / total_placements).powf(fragments))
        .collect();
    let expected_genes: f64 = hit_probabilities.iter().sum();
    let expected_covered: f64 = genes
        .iter()
        .zip(&placements)
        .map(|(gene, placements)| {
            geometry.expected_covered_bases(
                i64::from(gene.len),
                placements / total_placements,
                fragments,
            )
        })
        .sum();

    let count = usize::try_from(assembly.count).unwrap_or(usize::MAX);
    let covered = merged_interval_bases(&mut assembly.intervals, |record| record.span) as f64;
    let depth = assembly.reference_span_bases as f64 / marker_bases;
    let breadth = (covered / marker_bases).min(1.0);
    let expected_breadth = expected_covered / marker_bases;
    let gene_fraction = genes_hit as f64 / genes.len() as f64;
    let oebr = ratio(breadth, expected_breadth);
    let coebr = ratio(genes_hit as f64, expected_genes);
    let tests = tests.max(1) as f64;
    Some(EvidenceCoverage {
        count,
        depth,
        breadth,
        expected_breadth,
        fixed_chunk: gene_fraction,
        flex_chunk: gene_fraction,
        expected_flex_chunk: expected_genes / genes.len() as f64,
        oebr,
        coebr,
        breadth_pvalue: adjusted_deviation_pvalue(
            count,
            Some(depth),
            breadth,
            expected_breadth,
            oebr,
            tests,
        ),
        chunk_pvalue: Some((at_most(&hit_probabilities, genes_hit) * tests).min(1.0)),
        marker_genes_hit: Some(genes_hit),
    })
}

fn ratio(observed: f64, expected: f64) -> f64 {
    if expected > 0.0 {
        observed / expected
    } else {
        f64::NAN
    }
}

/// Probability that at most `k` of independent events with the given
/// probabilities occur (Poisson-binomial lower tail).
fn at_most(probabilities: &[f64], k: usize) -> f64 {
    // `dist[j]`: probability of exactly `j` events so far; mass above `k` is
    // never needed and is dropped.
    let mut dist = vec![0.0; k + 1];
    dist[0] = 1.0;
    for &p in probabilities {
        for j in (0..=k).rev() {
            let from_below = if j > 0 { dist[j - 1] * p } else { 0.0 };
            dist[j] = dist[j] * (1.0 - p) + from_below;
        }
    }
    dist.iter().sum::<f64>().min(1.0)
}

/// Lengths observed in the mappings, for a median.
#[derive(Debug, Default)]
pub(crate) struct LengthHistogram {
    counts: Vec<u64>,
}

impl LengthHistogram {
    const MAX_LEN: usize = 100_000;

    pub(crate) fn add(&mut self, len: usize) {
        let len = len.min(Self::MAX_LEN);
        if self.counts.len() <= len {
            self.counts.resize(len + 1, 0);
        }
        self.counts[len] += 1;
    }

    pub(crate) fn median(&self) -> Option<usize> {
        let total: u64 = self.counts.iter().sum();
        let mut seen = 0u64;
        self.counts.iter().position(|&count| {
            seen += count;
            total > 0 && 2 * seen >= total
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profiler::evidence::PresenceInterval;

    fn geometry(layout: FragmentLayout) -> FragmentGeometry {
        // 150 bp reads, 400 bp fragments, 105 bases inside a gene.
        FragmentGeometry::new(layout, 150, 400, 50, 0.7)
    }

    #[test]
    fn placements_follow_the_fragment_layout() {
        let reads = geometry(FragmentLayout::Reads);
        assert_eq!(reads.min_inside, 105);
        assert_eq!(reads.placements(1_000), 941.0);
        assert_eq!(reads.placements(100), 0.0);
        assert_eq!(geometry(FragmentLayout::Pairs).placements(1_000), 601.0);
        assert_eq!(geometry(FragmentLayout::Pairs).placements(300), 0.0);
        let lone = geometry(FragmentLayout::PairsAndLoneMates);
        assert_eq!(lone.placements(1_000), 941.0 + 250.0);
        // On a short gene the two mates' start ranges do not overlap.
        assert_eq!(lone.placements(250), 2.0 * 191.0);
    }

    #[test]
    fn expected_coverage_is_lower_near_gene_ends() {
        let reads = geometry(FragmentLayout::Reads);
        let many = reads.expected_covered_bases(1_000, 1.0, 10_000.0);
        assert!((many - 1_000.0).abs() < 1e-6);
        // One read covers 150 bases, fewer when it overhangs a gene end.
        let one = reads.expected_covered_bases(1_000, 1.0, 1.0);
        assert!(one < 150.0 && one > 130.0, "{one}");
    }

    #[test]
    fn poisson_binomial_lower_tail_is_exact() {
        assert!((at_most(&[0.5, 0.5], 0) - 0.25).abs() < 1e-12);
        assert!((at_most(&[0.5, 0.5], 1) - 0.75).abs() < 1e-12);
        assert!((at_most(&[0.2, 0.9, 0.5], 3) - 1.0).abs() < 1e-12);
        assert!((at_most(&[0.2, 0.9, 0.5], 0) - 0.8 * 0.1 * 0.5).abs() < 1e-12);
    }

    fn assembly(hits: &[(u32, u32)]) -> PresenceAssembly {
        PresenceAssembly {
            asm_len: 0.0,
            count: hits.len() as u64,
            reference_span_bases: 150 * hits.len() as u64,
            intervals: hits
                .iter()
                .map(|&(reference_id, start)| PresenceInterval {
                    reference_id,
                    span: (start, start + 150),
                    fixed_chunk: (start, start + 150),
                    flex_chunk: (start, start + 150),
                })
                .collect(),
        }
    }

    #[test]
    fn reads_spread_over_genes_pass_and_clumped_reads_do_not() {
        let genes: Vec<Gene> = (0..20)
            .map(|reference_id| Gene {
                reference_id,
                len: 1_000,
            })
            .collect();
        let geometry = geometry(FragmentLayout::PairsAndLoneMates);
        let spread: Vec<(u32, u32)> = (0..20).map(|gene| (gene, 400)).collect();
        let even = marker_coverage(&genes, &geometry, &mut assembly(&spread), 1).unwrap();
        assert_eq!(even.marker_genes_hit, Some(20));
        assert!(
            even.coebr > 1.0 && even.chunk_pvalue == Some(1.0),
            "{even:?}"
        );

        let clumped: Vec<(u32, u32)> = (0..20).map(|read| (read % 2, 400)).collect();
        let clumped = marker_coverage(&genes, &geometry, &mut assembly(&clumped), 1).unwrap();
        assert_eq!(clumped.marker_genes_hit, Some(2));
        assert!(clumped.coebr < 0.2, "{clumped:?}");
        assert!(clumped.chunk_pvalue.unwrap() < 1e-5, "{clumped:?}");
        assert_eq!(clumped.fixed_chunk, 0.1);
    }

    #[test]
    fn length_histogram_median() {
        let mut lengths = LengthHistogram::default();
        assert_eq!(lengths.median(), None);
        for len in [150, 150, 151, 400, 1_000_000] {
            lengths.add(len);
        }
        assert_eq!(lengths.median(), Some(151));
    }
}

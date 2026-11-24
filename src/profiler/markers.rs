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
use super::statistics::NEUTRAL_RATIO;
use crate::db::{DatabaseMetadata, ReferenceType};

/// Marker genes of the marker genomes and how this sample's fragments land on them.
#[derive(Debug)]
pub(crate) struct MarkerModel {
    /// Genes of every marker genome by assembly id, in reference order.
    genes: HashMap<u32, Vec<Gene>>,
    /// Fragment placements summed over the genes of every marker genome
    /// that can hold a read, by assembly id.
    placements: HashMap<u32, f64>,
    /// Gene lengths of every marker genome by assembly id, for the
    /// placements of each long read's own length; empty for other reads.
    long_read_genes: HashMap<u32, GeneLengths>,
    /// The whole-genome size that each marker genome declares, by assembly
    /// id.
    declared: HashMap<u32, usize>,
    geometry: FragmentGeometry,
}

/// The gene lengths of one genome for long reads, sorted, and their sum.
#[derive(Debug)]
struct GeneLengths {
    lens: Box<[u32]>,
    total: i64,
}

impl GeneLengths {
    /// Index of the first gene that needs `inside` bases of a read; the
    /// shorter genes before it need more than half of themselves,
    /// `len / 2 + 1` bases.
    fn first_full_gate(&self, inside: i64) -> usize {
        self.lens
            .partition_point(|&len| i64::from(len / 2 + 1) < inside)
    }

    /// Read starts that put `min(inside, len / 2 + 1)` bases of a read of
    /// `read_len` bases inside each gene of `len` bases, summed over the
    /// genes: `sum(len + read_len + 1 - 2 min(inside, len / 2 + 1))`.
    fn placements(&self, inside: i64, read_len: i64) -> f64 {
        let first = self.first_full_gate(inside);
        let below_bases: i64 = self.lens[..first]
            .iter()
            .map(|&len| i64::from(len / 2 + 1))
            .sum();
        let genes = self.lens.len() as i64;
        let inside_bases = below_bases + inside * (genes - first as i64);
        (self.total + genes * (read_len + 1) - 2 * inside_bases) as f64
    }

    /// [`Self::placements`] over the genes that need all `inside` bases.
    fn full_gate_placements(&self, inside: i64, read_len: i64) -> f64 {
        let genes = &self.lens[self.first_full_gate(inside)..];
        let bases: i64 = genes.iter().map(|&len| i64::from(len)).sum();
        (bases + genes.len() as i64 * (read_len + 1 - 2 * inside)) as f64
    }
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
        assembly_of: &[u32],
        geometry: FragmentGeometry,
    ) -> Option<Self> {
        let mut genes: HashMap<u32, Vec<Gene>> = HashMap::new();
        let mut declared = HashMap::new();
        for (reference_id, reference) in db.references.iter().enumerate() {
            if reference.reference_type == ReferenceType::Markers {
                let assembly_id = assembly_of[reference_id];
                genes.entry(assembly_id).or_default().push(Gene {
                    reference_id: u32::try_from(reference_id).ok()?,
                    len: u32::try_from(reference.length).ok()?,
                });
                if let Some(size) = reference.declared_genome_size {
                    declared.insert(assembly_id, size);
                }
            }
        }
        (!genes.is_empty()).then(|| Self::new(genes, declared, geometry))
    }

    fn new(
        genes: HashMap<u32, Vec<Gene>>,
        declared: HashMap<u32, usize>,
        geometry: FragmentGeometry,
    ) -> Self {
        let placements = genes
            .iter()
            .filter_map(|(&assembly_id, genes)| {
                let placements: f64 = genes
                    .iter()
                    .map(|gene| geometry.placements(i64::from(gene.len)))
                    .sum();
                (placements > 0.0).then_some((assembly_id, placements))
            })
            .collect();
        let long_read_genes = if geometry.layout == FragmentLayout::LongReads {
            genes
                .iter()
                .map(|(&assembly_id, genes)| {
                    let mut lens: Box<[u32]> = genes.iter().map(|gene| gene.len).collect();
                    lens.sort_unstable();
                    let total = lens.iter().map(|&len| i64::from(len)).sum();
                    (assembly_id, GeneLengths { lens, total })
                })
                .collect()
        } else {
            HashMap::new()
        };
        Self {
            genes,
            placements,
            long_read_genes,
            declared,
            geometry,
        }
    }

    /// A model of marker genomes given as `(assembly id, gene lengths)`,
    /// without declared genome sizes.
    #[cfg(test)]
    pub(crate) fn from_gene_lengths(genomes: &[(u32, &[u32])], geometry: FragmentGeometry) -> Self {
        let mut reference_id = 0;
        let genes = genomes
            .iter()
            .map(|&(assembly_id, lens)| {
                let genes = lens
                    .iter()
                    .map(|&len| {
                        reference_id += 1;
                        Gene { reference_id, len }
                    })
                    .collect();
                (assembly_id, genes)
            })
            .collect();
        Self::new(genes, HashMap::new(), geometry)
    }

    /// The model with the whole-genome sizes `declared` by assembly id.
    #[cfg(test)]
    pub(crate) fn with_declared(self, declared: HashMap<u32, usize>) -> Self {
        Self { declared, ..self }
    }

    pub(crate) fn is_marker_genome(&self, assembly_id: u32) -> bool {
        self.genes.contains_key(&assembly_id)
    }

    /// Fragment placements that the marker genes of a genome offer, summed
    /// over its genes; `None` for a genome without marker genes or when no
    /// gene can hold a read of this length. A genome of `G` bases offers
    /// about `G` placements, so the fragments on its marker genes give the
    /// depth of the whole genome over this length.
    pub(crate) fn fragment_placements(&self, assembly_id: u32) -> Option<f64> {
        self.placements.get(&assembly_id).copied()
    }

    /// Placements that the marker genes of a genome offer the fragment of a
    /// read of `read_len` bases: for a long read, those of its own length,
    /// since long reads differ widely in length; otherwise
    /// [`Self::fragment_placements`].
    pub(crate) fn placements_for_read(&self, assembly_id: u32, read_len: u32) -> Option<f64> {
        if self.geometry.layout != FragmentLayout::LongReads {
            return self.fragment_placements(assembly_id);
        }
        let read_len = i64::from(read_len.max(1));
        self.long_read_genes
            .get(&assembly_id)
            .map(|genes| genes.placements(self.geometry.gate_bases(read_len), read_len))
    }

    /// Whole-genome reads that a long read of `read_len` bases with hits on
    /// the genes `reference_ids` of a genome stands for: the hits that count
    /// and, per hit, `(a, b)` with the reads `a G + b` for a genome of `G`
    /// bases. Only genes that need all `m` bases the gates ask of the read
    /// count, as a read's hits on them are found reliably, while hits on
    /// shorter genes, which need more than half of the gene, are often
    /// missed. A hit stands for the placements of the whole genome,
    /// `G + read_len + 1 - 2 m`, over those of these genes. `None` for other
    /// reads and for a genome without such genes.
    pub(crate) fn whole_genome_reads(
        &self,
        assembly_id: u32,
        read_len: u32,
        reference_ids: impl IntoIterator<Item = u32>,
    ) -> Option<(u32, (f64, f64))> {
        if self.geometry.layout != FragmentLayout::LongReads {
            return None;
        }
        let read_len = i64::from(read_len.max(1));
        let inside = self.geometry.gate_bases(read_len);
        let on_genes = self
            .long_read_genes
            .get(&assembly_id)?
            .full_gate_placements(inside, read_len);
        if on_genes <= 0.0 {
            return None;
        }
        let genes = self.genes.get(&assembly_id)?;
        let hits = reference_ids
            .into_iter()
            .filter(|&reference_id| {
                genes
                    .binary_search_by_key(&reference_id, |gene| gene.reference_id)
                    .is_ok_and(|index| i64::from(genes[index].len / 2 + 1) >= inside)
            })
            .count();
        let outside = read_len + 1 - 2 * inside;
        Some((
            u32::try_from(hits).ok()?,
            (1.0 / on_genes, outside as f64 / on_genes),
        ))
    }

    /// Fragments that a whole genome of `genome_size` bases offers per
    /// fragment that its marker genes offer; `None` for a genome without
    /// marker genes or when no gene can hold a read of this length.
    pub(crate) fn genome_scale(&self, assembly_id: u32, genome_size: usize) -> Option<f64> {
        let on_genes = self.fragment_placements(assembly_id)?;
        let genome = self.geometry.placements(i64::try_from(genome_size).ok()?);
        Some(genome / on_genes)
    }

    /// [`Self::genome_scale`] at the genome size that genome `assembly_id`
    /// declares; `None` also without a declared size.
    pub(crate) fn declared_genome_scale(&self, assembly_id: u32) -> Option<f64> {
        self.genome_scale(assembly_id, self.declared_size(assembly_id)?)
    }

    /// The whole-genome size that genome `assembly_id` declares.
    pub(crate) fn declared_size(&self, assembly_id: u32) -> Option<usize> {
        self.declared.get(&assembly_id).copied()
    }

    /// Whether the sample's reads are long reads, whose placements on the
    /// genes depend on each read's own length.
    pub(crate) fn long_reads(&self) -> bool {
        self.geometry.layout == FragmentLayout::LongReads
    }

    /// The depth that a read of `read_len` bases on `genes` genes of genome
    /// `assembly_id` adds to the whole genome: its length on each gene over
    /// the placements that the genome's genes offer such a read
    /// ([`Self::placements_for_read`]); `None` when no gene can hold it.
    pub(crate) fn read_depth(&self, assembly_id: u32, read_len: u32, genes: u32) -> Option<f64> {
        self.placements_for_read(assembly_id, read_len)
            .map(|placements| f64::from(genes) * f64::from(read_len) / placements)
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
    /// Single-end short reads.
    Reads,
    /// Long reads, which count on a gene with the length or the read fraction
    /// of the length-or-fraction gate inside it, or more than half of the
    /// gene when that is shorter.
    LongReads,
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
    min_mapped_len: i64,
    min_fraction: f64,
    /// [`Self::gate_bases`] of the sample's read length.
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
        let template_len = i64::try_from(template_len)
            .unwrap_or(i64::MAX)
            .max(read_len);
        let mut geometry = Self {
            layout,
            read_len,
            template_len,
            min_mapped_len: i64::try_from(min_mapped_len).unwrap_or(i64::MAX),
            min_fraction: min_fraction.clamp(0.0, 1.0),
            min_inside: 0,
        };
        geometry.min_inside = geometry.gate_bases(read_len);
        geometry
    }

    /// Bases of a read of `read_len` bases that the read-level gates need on
    /// a reference: both the mapped length and the read fraction for short
    /// reads, either of them for long reads.
    fn gate_bases(&self, read_len: i64) -> i64 {
        let by_fraction = (self.min_fraction * read_len as f64).ceil() as i64;
        let bases = if self.layout == FragmentLayout::LongReads {
            self.min_mapped_len.min(by_fraction)
        } else {
            self.min_mapped_len.max(by_fraction)
        };
        bases.clamp(1, read_len)
    }

    /// Bases of a read needed inside a gene of `len` bases: for long reads at
    /// most more than half of the gene.
    fn min_inside(&self, len: i64) -> i64 {
        if self.layout == FragmentLayout::LongReads {
            self.min_inside.min(len / 2 + 1)
        } else {
            self.min_inside
        }
    }

    /// Start positions of a read with enough bases on a gene.
    fn read_starts(&self, len: i64) -> i64 {
        let inside = self.min_inside(len);
        if len < inside {
            return 0;
        }
        len - 2 * inside + self.read_len + 1
    }

    /// Fragment placements that count on a gene of `len` bases.
    fn placements(&self, len: i64) -> f64 {
        let reads = self.read_starts(len);
        let placements = match self.layout {
            FragmentLayout::Reads | FragmentLayout::LongReads => reads,
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
            _ => (
                self.min_inside(len) - self.read_len,
                len - self.min_inside(len),
            ),
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
        let covered = |c: i64| 1.0 - (1.0 - share * c as f64 / starts as f64).powf(fragments);
        let starts_over = |base: i64| {
            (base.min(last) - (base - self.read_len + 1).max(first) + 1).clamp(0, self.read_len)
        };
        if self.read_len < len {
            // Fewer start counts than bases: work each out once.
            let by_starts: Vec<f64> = (0..=self.read_len).map(covered).collect();
            (0..len)
                .map(|base| by_starts[starts_over(base) as usize])
                .sum()
        } else {
            (0..len).map(|base| covered(starts_over(base))).sum()
        }
    }
}

impl std::fmt::Display for FragmentGeometry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let layout = match self.layout {
            FragmentLayout::Reads => "reads",
            FragmentLayout::LongReads => {
                return write!(
                    f,
                    "long reads, median read length {}, at least {} bases of a read or more than half of a gene inside it",
                    self.read_len, self.min_inside
                );
            }
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
    // Gene hits: a long read, or the two mates of a pair, can cover several
    // genes.
    let fragments = (assembly.count + assembly.other_gene_hits) as f64;
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
    // Too few genes hit, tested only outside the window the genome p-values use.
    let chunk_pvalue = if NEUTRAL_RATIO.contains(&coebr) {
        1.0
    } else {
        (at_most(&hit_probabilities, genes_hit) * tests).min(1.0)
    };
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
        chunk_pvalue: Some(chunk_pvalue),
        marker_genes_hit: Some(genes_hit),
        window_spread: None,
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
            other_gene_hits: 0,
            reference_span_bases: 150 * hits.len() as u64,
            intervals: hits
                .iter()
                .map(|&(reference_id, start)| PresenceInterval {
                    reference_id,
                    span: (start, start + 150),
                    mate_span: (0, 0),
                    fixed_chunk: (start, start + 150),
                    flex_chunk: (start, start + 150),
                })
                .collect(),
            fits: 0,
            sampling_scale: 1.0,
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
    fn a_small_gene_shortfall_with_many_reads_is_not_tested() {
        let genes: Vec<Gene> = (0..200)
            .map(|reference_id| Gene {
                reference_id,
                len: 1_000,
            })
            .collect();
        let geometry = geometry(FragmentLayout::PairsAndLoneMates);
        let coverage = |genes_hit: u32| {
            let hits: Vec<(u32, u32)> = (0..4_000).map(|read| (read % genes_hit, 400)).collect();
            marker_coverage(&genes, &geometry, &mut assembly(&hits), 1).unwrap()
        };
        // 4,000 fragments reach nearly every gene, so the exact test alone
        // would reject 190 of 200 genes.
        let p = 1.0 - (1.0 - 1.0 / 200.0f64).powf(4_000.0);
        assert!(at_most(&[p; 200], 190) < 1e-5);
        let near = coverage(190);
        assert!(NEUTRAL_RATIO.contains(&near.coebr), "{near:?}");
        assert_eq!(near.chunk_pvalue, Some(1.0));
        let short = coverage(170);
        assert!(short.coebr < 0.9, "{short:?}");
        assert!(short.chunk_pvalue.unwrap() < 1e-5, "{short:?}");
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

    #[test]
    fn long_reads_count_on_genes_shorter_than_the_gates() {
        // Either gate for long reads: 500 bases, as 20 % of 10,000 is more.
        let long = FragmentGeometry::new(FragmentLayout::LongReads, 10_000, 10_000, 500, 0.2);
        assert_eq!(long.min_inside, 500);
        assert_eq!(long.placements(1_000), 10_001.0);
        // A 300-base gene counts a read with more than half of it, 151 bases.
        assert_eq!(long.placements(300), (300 - 2 * 151 + 10_001) as f64);
        // Every base of a short gene can be covered.
        assert!((long.expected_covered_bases(300, 1.0, 1e6) - 300.0).abs() < 1e-9);
        // Short single-end reads need both gates.
        let short = FragmentGeometry::new(FragmentLayout::Reads, 150, 150, 50, 0.7);
        assert_eq!(short.min_inside, 105);
        assert_eq!(short.placements(100), 0.0);
    }

    #[test]
    fn a_long_read_is_placed_by_its_own_length() {
        let lens = [300u32, 800, 2_000];
        let at = |read_len: usize| {
            FragmentGeometry::new(FragmentLayout::LongReads, read_len, read_len, 500, 0.2)
        };
        let model = MarkerModel::from_gene_lengths(&[(0, &lens)], at(3_000));
        for read_len in [200u32, 1_000, 3_000, 10_000] {
            let geometry = at(read_len as usize);
            let expected: f64 = lens
                .iter()
                .map(|&len| geometry.placements(i64::from(len)))
                .sum();
            assert_eq!(model.placements_for_read(0, read_len), Some(expected));
        }
        assert_eq!(model.placements_for_read(1, 3_000), None);
        // Other layouts keep the placements of the sample's lengths.
        let pairs = MarkerModel::from_gene_lengths(
            &[(0, &lens)],
            FragmentGeometry::new(FragmentLayout::Pairs, 150, 400, 50, 0.7),
        );
        assert_eq!(
            pairs.placements_for_read(0, 10_000),
            pairs.fragment_placements(0)
        );
    }

    #[test]
    fn the_estimate_counts_hits_on_genes_that_need_the_whole_gate() {
        // Genes of 300, 800 and 2,000 bases: reference ids 1, 2 and 3.
        let lens = [300u32, 800, 2_000];
        let geometry = FragmentGeometry::new(FragmentLayout::LongReads, 3_000, 3_000, 500, 0.2);
        let model = MarkerModel::from_gene_lengths(&[(0, &lens)], geometry);
        // A 3,000-base read needs 500 bases: only the 2,000-base gene needs all.
        let (hits, per_hit) = model.whole_genome_reads(0, 3_000, [1, 2, 3]).unwrap();
        let on_genes = f64::from(2_000 + 3_001 - 2 * 500);
        assert_eq!(hits, 1);
        assert_eq!(
            per_hit,
            (1.0 / on_genes, f64::from(3_001 - 2 * 500) / on_genes)
        );
        // A 1,000-base read needs 200 bases: the 800- and 2,000-base genes count.
        let (hits, (per_base, _)) = model.whole_genome_reads(0, 1_000, [1, 2]).unwrap();
        assert_eq!(hits, 1);
        assert_eq!(
            per_base,
            1.0 / f64::from(800 + 2_000 + 2 * (1_001 - 2 * 200))
        );
        // A read with hits on short genes only still counts, with no hit.
        assert_eq!(model.whole_genome_reads(0, 3_000, [1]).unwrap().0, 0);
        assert!(model.whole_genome_reads(1, 3_000, [1]).is_none());
    }

    #[test]
    fn a_whole_genome_offers_its_placements_over_those_of_its_genes() {
        // Pairs of 400 bases: a gene of L bases offers L - 400 + 1 placements.
        let geometry = FragmentGeometry::new(FragmentLayout::Pairs, 150, 400, 50, 0.7);
        let model = MarkerModel::from_gene_lengths(&[(0, &[1_000, 500])], geometry);
        let scale = model.genome_scale(0, 10_000).unwrap();
        assert!((scale - 9_601.0 / 702.0).abs() < 1e-9);
        assert_eq!(model.genome_scale(1, 10_000), None);
        // The same at the size the genome declares; none without one.
        assert_eq!(model.declared_genome_scale(0), None);
        let model = model.with_declared(HashMap::from([(0, 10_000)]));
        assert_eq!(model.declared_genome_scale(0), Some(scale));
    }

    #[test]
    fn marker_genes_offer_their_summed_placements() {
        // 150 bp pairs of 400 bases with lone mates: a gene of 1,000 bases
        // offers 941 + 250 placements, one of 250 bases 2 x 191, and one of
        // 100 bases none.
        let genes = |lens: &[u32]| -> Vec<Gene> {
            lens.iter()
                .enumerate()
                .map(|(reference_id, &len)| Gene {
                    reference_id: reference_id as u32,
                    len,
                })
                .collect()
        };
        let model = MarkerModel::new(
            HashMap::from([(0, genes(&[1_000, 250, 100])), (1, genes(&[100]))]),
            HashMap::new(),
            geometry(FragmentLayout::PairsAndLoneMates),
        );
        assert_eq!(model.fragment_placements(0), Some(1_191.0 + 382.0));
        assert_eq!(model.fragment_placements(1), None);
        assert_eq!(model.fragment_placements(2), None);
    }
}

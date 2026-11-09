//! Community Profiler - terminal-rank taxonomic abundance estimation
//!
//! This module implements the default profiling mode. It reports each
//! reference at the lowest taxonomy rank that reference's metadata
//! confidently supports, rather than rolling every assignment up to species.
//!
//! # Algorithm Overview
//!
//! 1. **Load mappings**: Read Metam coordinate mappings, either streamed from
//!    the mapper or replayed from a saved mapping file.
//!
//! 2. **Reference Metadata**: Take genome, taxonomy, size, and reference type
//!    of every reference from the database metadata. MANUAL.md describes the
//!    reference format.
//!
//! 3. **Coverage Calculation**: Compute coverage metrics:
//!    - Breadth: Fraction of genome covered by at least one read
//!    - Fixed and flex chunk: Coverage of genome chunks by at least one read;
//!      for marker genomes the genes are the chunks
//!
//! 4. **Taxonomy Assignment**: Assign reads to taxa based on alignment
//!    quality, filtering by identity, mapped length, and fraction thresholds.
//!
//! 5. **Ambiguity Resolution**: Multi-mapping reads are resolved using an
//!    Expectation-Maximization (EM) algorithm that iteratively redistributes
//!    reads based on relative taxon abundances.
//!
//! 6. **Statistical Filtering**: Final taxa are filtered based on:
//!    - Breadth/expected breadth ratio
//!    - Chunk breadth/expected ratio
//!    - P-values from coverage distribution tests
//!
//! # Key Data Structures
//!
//! - `ReadBundle`: Groups all alignments for a single read
//! - `AlignmentEntry`: Information about a single read-to-reference alignment
//! - `TaxonCoverage`: Coverage metrics aggregated at the terminal taxon
//! - `ReadAssignment`: Final taxonomic assignment for a read

use std::cmp::{Ordering, Reverse};
use std::collections::{hash_map::Entry, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;
use std::str;
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::taxonomy::Taxonomy;
use csv::WriterBuilder;

use crate::db::{DatabaseMetadata, ReferenceType};
use crate::mapper::tsv::{
    for_each_metam_mapping_record, for_each_structured_mapping_record, load_shared_metam_profile,
    open_metam_mapping_reader, total_input_units, MetamMappingReceiver, MetamMappingRecord,
    MetamParseStats, MetamProfile,
};
use crate::progress;

use super::evidence::{
    calculate_coverage, evaluate_ranked_witness, guarded_map_winner, resolved_candidate_index,
    standard_ratio_lower_bounds, witness_is_better, CandidatePlacement, CandidateScratch,
    EvidenceAlignment, EvidenceCoverage, EvidenceMode, EvidencePolicy, EvidenceThresholds,
    ExactTiePlacement, ExactTies, MappingRole, PresenceAssembly, PresenceInterval, PresenceKey,
    TaxonEvidence, WitnessDecision, STANDARD_RATIO_UPPER,
};
#[cfg(test)]
use super::evidence::{PremapChunkGuard, RescueKinds};
use super::identity::{
    identity_baseline, move_rows_to_genus, BestPlacementScratch, BestPlacements,
};
use super::markers::{FragmentGeometry, FragmentLayout, LengthHistogram, MarkerModel};
use super::output::{
    BestFitReport, ClassificationSource, ProfileAnnotations, ProfileRow, ResolvedClassifications,
};
use super::pathogen::{load_pathogen_table, PathogenEntry};
use super::rollup::{GenomeLayout, GenusRollup, GenusRow, OrphanedFragments};
use super::statistics::{coverage_deviation_pvalue, round5};
use super::unclassified::{self, GenomeReads, MarkerSpecies, ReferenceKind};
use super::{
    coverage_passes, equivalent_read_count_passes_minimum, log_reference_taxonomy_adjustments,
    log_resolved_profiler_parameters, parse_reference_taxonomy, push_string_keyed_value,
    record_genome_sampling_scale, record_reference_taxonomy_adjustment, sampling_scale_from_sizes,
    scale_automatic_chunk_breadth, strip_pair_suffix, ProfileMode, ProfileRank, ProfilerConfig,
    ReadSupport, ReferenceTaxonomyAdjustments, ReferenceTaxonomyResolver,
    ResolvedProfilerParameters, StableReadMap,
};

/// Fixed coordinate radius around each accepted read span for fixed-chunk
/// coverage.
const NATIVE_FIXED_FLANK_LEN: usize = 1_000;

fn record_terminal_rank(
    ranks: &mut HashMap<u32, ProfileRank>,
    taxid: u32,
    rank: ProfileRank,
    reference: &str,
) -> Result<()> {
    match ranks.entry(taxid) {
        Entry::Vacant(entry) => {
            entry.insert(rank);
            Ok(())
        }
        Entry::Occupied(entry) if *entry.get() == rank => Ok(()),
        Entry::Occupied(entry) => bail!(
            "conflicting terminal ranks for taxid {}: '{}' and '{}' (reference '{}')",
            taxid,
            entry.get().as_str(),
            rank.as_str(),
            reference,
        ),
    }
}

/// Information about a single alignment of a read to a reference.
///
/// Derived fields are computed once where the alignment enters the profiler --
/// in [`CommunityProfiler::entry_from_metam_mapping`] -- and stored here ready
/// to use. They
/// are invariant for a given alignment, so the taxonomy-assignment loop is
/// purely arithmetic and touches no `HashMap` or `String::split` per
/// alignment.
type AlignmentEntry = EvidenceAlignment;

/// Narrow a genomic coordinate/length to `u32` for the packed
/// [`AlignmentEntry`] fields. Coordinates are bounded by `ref_len + 1`, and the
/// index rejects references `>= 2^32` bases, so the `.min` only guards the
/// unreachable 4 Gbp-contig edge — every real value is preserved exactly.
#[inline]
fn clamp_u32(value: usize) -> u32 {
    value.min(u32::MAX as usize) as u32
}

/// Widen a packed `(u32, u32)` coordinate pair back to `(usize, usize)` when
/// feeding the (much smaller) per-taxon interval maps and coverage merge.
#[inline]
fn widen_pair(pair: (u32, u32)) -> (usize, usize) {
    (pair.0 as usize, pair.1 as usize)
}

/// Bundle of all alignments for a single read.
///
/// A read may align to multiple reference sequences (multi-mapping).
/// This structure groups all alignments for downstream processing.
#[derive(Clone, Debug)]
struct ReadBundle {
    /// Read identifier (pair suffix stripped for paired-end)
    read_name: Arc<str>,

    /// All alignments for this read, ordered by quality
    alignments: Vec<AlignmentEntry>,
}

/// Coverage metrics aggregated at a terminal taxonomy endpoint.
#[derive(Clone, Debug, Default)]
struct TaxonCoverage {
    /// Breadth of coverage: fraction of genome with ≥1 read
    breadth: f64,

    /// Fixed chunk coverage: uses fixed-size flanking regions
    fixed_chunk: f64,

    /// Flexible chunk coverage: uses √(genome_size)-sized regions
    flex_chunk: f64,
}

/// Unique key identifying a specific accession within an assembly.
///
/// Used as a HashMap key for aggregating coverage intervals.
#[derive(Clone, Debug, Hash, Ord, PartialEq, Eq, PartialOrd)]
struct AccKey {
    /// Species taxonomy ID
    taxid: u32,

    /// Assembly identifier (first field of reference name)
    asm: Arc<str>,

    /// Accession identifier (fourth field of reference name)
    acc: Arc<str>,
}

/// Reference metadata pre-parsed for the Metam coordinate path.
///
/// Metam emits the reference name on every mapping row. Splitting that
/// reference name and re-looking-up assembly lengths in the hot loop creates
/// avoidable allocator and hash-table work. This cache moves that work to DB
/// load time, where each reference is parsed once.
#[derive(Clone, Debug)]
struct MetamReferenceEntry {
    ref_len: usize,
    asm: Arc<str>,
    acc: Arc<str>,
    profile_taxid: u32,
    profile_rank: ProfileRank,
    reference_id: u32,
    assembly_id: u32,
    sampling_scale: f32,
    asm_len: f64,
    flex_len: usize,
    is_marker: bool,
}

/// Statistics about alignment filtering.
///
/// Tracks how many alignments were evaluated, passed filters,
/// and why alignments were rejected.
#[derive(Clone, Copy, Debug, Default)]
struct FilterStats {
    evaluated: usize,

    passed: usize,

    /// Filtered due to insufficient mapped length
    filtered_mapped_len: usize,

    /// Filtered due to low sequence identity
    filtered_identity: usize,

    /// Filtered due to low aligned fraction
    filtered_fraction: usize,

    /// Skipped due to missing reference metadata
    skipped_reference: usize,
}

impl FilterStats {
    /// Merge statistics from a parallel worker.
    fn merge(&mut self, other: FilterStats) {
        self.evaluated += other.evaluated;
        self.passed += other.passed;
        self.filtered_mapped_len += other.filtered_mapped_len;
        self.filtered_identity += other.filtered_identity;
        self.filtered_fraction += other.filtered_fraction;
        self.skipped_reference += other.skipped_reference;
    }
}

/// Results from the taxonomy assignment phase.
///
/// Contains all data needed for abundance estimation and output generation.
#[derive(Default)]
struct AssignmentResult {
    /// Per-read taxonomic assignments
    read_assignments: Vec<ReadAssignment>,

    /// Alignment intervals per accession (for breadth calculation)
    acc_intervals: HashMap<AccKey, Vec<(usize, usize)>>,

    /// Fixed chunk intervals per accession
    acc_fixed_chunk_intervals: HashMap<AccKey, Vec<(usize, usize)>>,

    /// Flexible chunk intervals per accession
    acc_flex_chunk_intervals: HashMap<AccKey, Vec<(usize, usize)>>,

    /// Terminal-rank, assembly-coherent coordinate evidence for native short-
    /// and long-read profiling. Different assemblies are always independent
    /// witnesses and are never unioned.
    presence: HashMap<PresenceKey, PresenceAssembly>,

    /// Candidate, resolved-best, and direct-clear fragment counts used by the
    /// evidence-policy rescue tiers.
    taxon_evidence: HashMap<u32, TaxonEvidence>,

    /// Fragments whose best placement lies on each taxon, for
    /// `--identity-margin`.
    best_placements: HashMap<u32, BestPlacements>,

    /// Compact coordinate payload for exact top-score short-read ties.
    /// Populated only when guarded MAP is active.
    exact_ties: ExactTies,

    /// Terminal taxa observed as the only valid candidate for at least one
    /// fragment. This is the multiple-testing basis for the evidence policy.
    uniquely_supported_taxa: HashSet<u32>,

    /// Fragments for which a mapper boundary may hide an unrepresented
    /// profile endpoint. Emitted candidates remain in abundance EM, but
    /// cannot establish direct coordinate evidence from an incomplete set.
    boundary_ambiguous_fragments: usize,

    filter_stats: FilterStats,

    /// The fragments resolved to each marker genome, their gene hits and, for
    /// long reads, the whole-genome reads they stand for, by assembly id;
    /// filled only for `--estimate-unclassified`.
    marker_genome_reads: HashMap<u32, GenomeReads>,
}

/// Rank at which a read is taxonomically assigned.
///
/// Reads are normally assigned at the `Species` rank. When the genus
/// fallback mode activates (low-mapping Illumina, pathogen mode off),
/// reads that fail species thresholds but meet the fallback criteria
/// are assigned at the `Genus` rank instead.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum AssignRank {
    Species,
    Genus,
    Family,
    Order,
    Class,
    Phylum,
    Kingdom,
    /// Candidate endpoints span more than one terminal rank.
    Mixed,
}

impl AssignRank {
    fn from_profile_rank(rank: ProfileRank) -> Self {
        match rank {
            ProfileRank::Species => Self::Species,
            ProfileRank::Genus => Self::Genus,
            ProfileRank::Family => Self::Family,
            ProfileRank::Order => Self::Order,
            ProfileRank::Class => Self::Class,
            ProfileRank::Phylum => Self::Phylum,
            ProfileRank::Kingdom => Self::Kingdom,
        }
    }
}

/// Taxonomic assignment for a single read.
///
/// A read may be assigned to one taxon (unambiguous) or multiple taxa
/// (ambiguous multi-mapping). Ambiguous assignments are resolved in
/// the EM phase.
///
/// Native candidates can span terminal ranks; those assignments use
/// [`AssignRank::Mixed`]. Legacy and fallback assignments remain homogeneous.
#[derive(Clone, Debug)]
struct ReadAssignment {
    read: Arc<str>,

    /// Homogeneous terminal rank, genus fallback, or `Mixed`.
    rank: AssignRank,

    /// `true` when this assignment came from the native terminal-evidence
    /// resolver and its candidates must pass the assembly-coherent presence
    /// policy. Genus-fallback assignments use the pooled coverage path even
    /// when their taxid is also a direct database endpoint.
    requires_terminal_evidence: bool,

    /// Candidate terminal taxids (may be >1 for multi-mappers).
    taxids: Vec<u32>,

    /// Depth contribution for each candidate (mapped_len / genome_size)
    depths: Vec<f64>,

    /// Mapped lengths for each candidate
    mapped_lens: Vec<f64>,

    /// Genome-specific multiplier for each representative candidate mapping.
    sampling_scales: Vec<f64>,

    /// Index in `taxids` of the candidate with the best placement. Native
    /// assignments only; the genus roll-up credits an orphaned fragment to
    /// that species.
    best: u32,
}

// `best` fills alignment padding: the per-fragment record must not grow.
#[cfg(target_pointer_width = "64")]
const _: () = assert!(std::mem::size_of::<ReadAssignment>() == 120);

#[derive(Clone, Copy, Debug, Default)]
struct TaxonTotals {
    support: ReadSupport,
    depth: f64,
    mapped_len: f64,
}

#[derive(Clone, Debug, Default)]
struct AmbiguousTotals {
    observed_count: f64,
    depths: Vec<f64>,
    mapped_lens: Vec<f64>,
    /// Sum of the candidate genome scales over all reads in this group.
    equivalent_support_sums: Vec<f64>,
}

/// Optional, diagnostic abundance estimate built from every quality-qualified
/// read assignment before taxon-level coverage/evidence selection.
///
/// This is intentionally separate from the production EM state: enabling
/// `--save-raw-profile` must not add candidates to `profile.txt` or
/// `classify.txt`, or otherwise alter the established profiling policy.
struct RawProfileEstimate {
    totals: HashMap<u32, TaxonTotals>,
    coverage: HashMap<u32, TaxonCoverage>,
}

/// Community-level taxonomy profiler.
///
/// Implements terminal-rank profiling by grouping native mappings at each
/// reference-supported endpoint and resolving ambiguous assignments with EM.
///
/// # Pipeline
///
/// 1. Load coordinate mappings and extract reference metadata
/// 2. Compute coverage intervals for each reference sequence
/// 3. Assign reads to taxa based on alignment quality
/// 4. Resolve multi-mapping reads with EM algorithm
/// 5. Filter taxa based on coverage statistics
/// 6. Output abundance profile and classification files
pub struct CommunityProfiler {
    cfg: ProfilerConfig,

    /// Taxonomy database (shared across threads)
    taxonomy: Arc<Taxonomy>,
}

impl CommunityProfiler {
    /// Create a new community profiler with the given configuration.
    ///
    /// # Arguments
    /// * `cfg` - Profiler configuration
    /// * `taxonomy` - NCBI taxonomy database (wrapped in Arc)
    pub fn new(cfg: ProfilerConfig, taxonomy: Arc<Taxonomy>) -> Result<Self> {
        Ok(Self { cfg, taxonomy })
    }

    /// Run the community profiling pipeline.
    ///
    /// This is the main entry point that orchestrates the full pipeline:
    /// 1. Load and parse alignments
    /// 2. Compute coverage metrics
    /// 3. Assign taxonomy
    /// 4. Write output files
    pub fn run(mut self) -> Result<()> {
        log::info!(target: "PROFILER", "Initializing the profiler!");

        let mappings = self
            .cfg
            .mapper_mappings
            .clone()
            .ok_or_else(|| anyhow!("missing Metax mapping file"))?;
        let db = Arc::clone(
            self.cfg
                .db_metadata
                .as_ref()
                .ok_or_else(|| anyhow!("missing reference metadata for native mappings"))?,
        );
        let mapping_data = self.load_metam(Path::new(&mappings), &db)?;
        self.profile(mapping_data)
    }

    pub fn run_with_metam_mappings(
        mut self,
        receiver: MetamMappingReceiver,
        source: &str,
    ) -> Result<()> {
        log::info!(target: "PROFILER", "Initializing the profiler!");
        let db = Arc::clone(
            self.cfg
                .db_metadata
                .as_ref()
                .ok_or_else(|| anyhow!("missing Metax database metadata for streamed mappings"))?,
        );
        let mapping_data = self.load_metam_mappings(receiver, source, &db)?;
        self.profile(mapping_data)
    }

    /// Lower an automatic long-read identity cutoff when the reads align worse
    /// than their base qualities claim, as when a basecaller over-states its
    /// accuracy. The upper quartile of the reads' best alignment identities
    /// then takes the claimed accuracy's place: reads from well-represented
    /// genomes align at their sequencing accuracy, while divergent or poorly
    /// represented organisms lower only the lower quantiles.
    fn lower_identity_to_alignments(&mut self, reads: &[ReadBundle]) {
        let cfg = &self.cfg;
        let Some(claimed) = cfg.identity_check_accuracy else {
            return;
        };
        let score = |aln: &&AlignmentEntry| f64::from(aln.identity) * f64::from(aln.mapped_len);
        let mut best: Vec<f32> = reads
            .iter()
            .filter_map(|read| {
                read.alignments
                    .iter()
                    .filter(|aln| {
                        coverage_passes(
                            aln.mapped_len as usize >= cfg.mapped_len,
                            f64::from(aln.fraction) >= cfg.fraction,
                            cfg.coverage_gate_is_or(),
                        )
                    })
                    .max_by(|a, b| score(a).total_cmp(&score(b)))
                    .map(|aln| aln.identity)
            })
            .collect();
        if best.is_empty() {
            return;
        }
        let upper = best.len() * 3 / 4;
        let aligned = f64::from(*best.select_nth_unstable_by(upper, f32::total_cmp).1);
        let lowered = (aligned < claimed)
            .then(|| crate::adaptive_identity_from_accuracy(aligned, cfg.mode))
            .flatten()
            .filter(|&cutoff| cutoff < cfg.identity);
        if let Some(cutoff) = lowered {
            log::info!(
                target: "PROFILER",
                "Long-read identity threshold lowered to {cutoff:.4}: reads align worse than their base qualities claim (upper quartile of read identities {aligned:.4} < claimed accuracy {claimed:.4}).",
            );
            self.cfg.identity = cutoff;
        }
    }

    fn load_metam_profile(&self, embedded: Option<MetamProfile>) -> Result<MetamProfile> {
        if let Some(profile) = embedded {
            return Ok(profile);
        }
        if let Some(shared) = &self.cfg.mapping_profile_memory {
            return load_shared_metam_profile(shared);
        }
        bail!("mapping file is missing embedded profile counters; rerun mapping with this Metax version")
    }

    fn load_metam(&mut self, mappings_path: &Path, db: &DatabaseMetadata) -> Result<MappingData> {
        if !mappings_path.exists() {
            bail!(
                "Metam mapping file {} does not exist",
                mappings_path.display()
            );
        }
        let reader = open_metam_mapping_reader(mappings_path)?;
        self.load_metam_reader(reader, &mappings_path.display().to_string(), db)
    }

    fn load_metam_reader<R: BufRead>(
        &mut self,
        reader: R,
        source: &str,
        db: &DatabaseMetadata,
    ) -> Result<MappingData> {
        let lenient = self.cfg.lenient_mapping_input;
        self.load_metam_records(source, db, |on_mapping| {
            for_each_metam_mapping_record(reader, source, lenient, on_mapping)
        })
    }

    fn load_metam_mappings(
        &mut self,
        receiver: MetamMappingReceiver,
        source: &str,
        db: &DatabaseMetadata,
    ) -> Result<MappingData> {
        self.load_metam_records(source, db, |on_mapping| {
            for_each_structured_mapping_record(receiver, source, on_mapping)
        })
    }

    fn load_metam_records<F>(
        &mut self,
        source: &str,
        db: &DatabaseMetadata,
        feed: F,
    ) -> Result<MappingData>
    where
        F: FnOnce(
            &mut dyn for<'a> FnMut(MetamMappingRecord<'a>) -> Result<()>,
        ) -> Result<MetamParseStats>,
    {
        log::info!(target: "PROFILER", "Loading Metam coordinate mappings ...");
        let fixed_flank_len = NATIVE_FIXED_FLANK_LEN;
        let mut asm_lengths_raw: HashMap<String, f64> = HashMap::new();
        let mut genome_sampling_scales: HashMap<String, f64> = HashMap::new();
        let mut profile_taxids: HashSet<u32> = HashSet::new();
        let mut terminal_ranks: HashMap<u32, ProfileRank> = HashMap::new();
        let mut reference_taxonomies = Vec::with_capacity(db.references.len());
        let mut taxonomy_adjustments = ReferenceTaxonomyAdjustments::new();
        let mut taxonomy_resolver = ReferenceTaxonomyResolver::new(&self.taxonomy);
        let has_marker_genomes = db.uniform_reference_type() != Some(ReferenceType::Genome);
        let mut taxids_with_genomes: HashSet<u32> = HashSet::new();
        for reference in &db.references {
            let declared_taxonomy = parse_reference_taxonomy(&reference.taxonomy_field)
                .with_context(|| {
                    format!(
                        "invalid taxonomy metadata for reference '{}'",
                        reference.name
                    )
                })?;
            let taxonomy = taxonomy_resolver.resolve(declared_taxonomy, &reference.name)?;
            record_reference_taxonomy_adjustment(
                &mut taxonomy_adjustments,
                declared_taxonomy,
                taxonomy,
            );
            let asm = reference.genome.clone();
            *asm_lengths_raw.entry(asm).or_insert(0.0) += reference.length as f64;
            profile_taxids.insert(taxonomy.profile_taxid);
            if reference.reference_type == ReferenceType::Genome {
                taxids_with_genomes.insert(taxonomy.profile_taxid);
            }
            record_terminal_rank(
                &mut terminal_ranks,
                taxonomy.profile_taxid,
                taxonomy.profile_rank,
                &reference.name,
            )?;
            let scale = sampling_scale_from_sizes(
                reference.full_genome_size,
                reference.effective_indexed_size,
            );
            record_genome_sampling_scale(&mut genome_sampling_scales, &reference.genome, scale)?;
            reference_taxonomies.push(taxonomy);
        }
        log_reference_taxonomy_adjustments(
            &taxonomy_adjustments,
            self.cfg.verbose || self.cfg.very_verbose,
        );
        let asm_flex_lengths: HashMap<String, usize> = asm_lengths_raw
            .iter()
            .map(|(asm, len)| {
                let flex = len.sqrt().round() as usize;
                (asm.clone(), flex.max(1))
            })
            .collect();
        let (species_normalize_raw, mut species_to_genus_raw) =
            self.build_taxonomy_lookups(&profile_taxids);
        let normalize = |taxid: u32| species_normalize_raw.get(&taxid).copied().unwrap_or(taxid);
        let mut marker_only_species: HashMap<u32, Option<u32>> = HashMap::new();
        if has_marker_genomes {
            // Genus fallback places reads that no species explains at the genus
            // of their best hit. A marker gene holds too little of a genome for
            // that: reads failing every species there come mostly from the
            // homologous genes of relatives. Only species with genome-type
            // references take part.
            let species_with_genomes: HashSet<u32> = taxids_with_genomes
                .iter()
                .map(|&taxid| normalize(taxid))
                .collect();
            if self.cfg.removed_to_genus {
                marker_only_species = profile_taxids
                    .iter()
                    .map(|&taxid| normalize(taxid))
                    .filter(|species| !species_with_genomes.contains(species))
                    .map(|species| (species, species_to_genus_raw.get(&species).copied()))
                    .collect();
            }
            species_to_genus_raw.retain(|species, _| species_with_genomes.contains(species));
        }
        // The roll-up needs every genus, also on a fractional database, where
        // genus fallback is off.
        let rollup_genus_of = self
            .cfg
            .removed_to_genus
            .then(|| species_to_genus_raw.clone());
        let db_is_subsampled = db.is_fractional;
        let species_normalize = Arc::new(species_normalize_raw);
        let species_to_genus = Arc::new(if db_is_subsampled && !self.cfg.genus_fallback {
            HashMap::new()
        } else {
            species_to_genus_raw
        });
        let mut assembly_ids: HashMap<&str, u32> = HashMap::new();
        let mut assembly_names: Vec<Arc<str>> = Vec::new();
        for reference in &db.references {
            if !assembly_ids.contains_key(reference.genome.as_str()) {
                let next = u32::try_from(assembly_ids.len())
                    .context("database contains more than 2^32 assemblies")?;
                assembly_ids.insert(reference.genome.as_str(), next);
                assembly_names.push(Arc::from(reference.genome.as_str()));
            }
        }
        let genus_rollup = rollup_genus_of.map(|genus_of| {
            let layout = GenomeLayout::new(db.references.iter().map(|reference| {
                (
                    assembly_ids[reference.genome.as_str()],
                    reference.length as u64,
                )
            }));
            Arc::new(GenusRollup {
                genus_of,
                marker_only: marker_only_species,
                layout,
                fractional: db_is_subsampled,
            })
        });
        let mut metam_refs: HashMap<String, MetamReferenceEntry> =
            HashMap::with_capacity(db.references.len());
        for (reference_index, (reference, taxonomy)) in
            db.references.iter().zip(&reference_taxonomies).enumerate()
        {
            let asm = reference.genome.as_str();
            let raw_profile_taxid = taxonomy.profile_taxid;
            let acc = reference.contig.as_str();
            let asm_len = match asm_lengths_raw.get(asm).copied() {
                Some(value) if value > 0.0 => value,
                _ => continue,
            };
            let profile_taxid = species_normalize
                .get(&raw_profile_taxid)
                .copied()
                .unwrap_or(raw_profile_taxid);
            let flex_len = asm_flex_lengths
                .get(asm)
                .copied()
                .unwrap_or(fixed_flank_len);
            let sampling_scale = genome_sampling_scales.get(asm).copied().unwrap_or(1.0) as f32;
            let reference_id = u32::try_from(reference_index)
                .context("database contains more than 2^32 references")?;
            let assembly_id = assembly_ids[asm];
            if metam_refs.contains_key(&reference.name) {
                bail!(
                    "database metadata contains duplicate reference name '{}'",
                    reference.name
                );
            }
            metam_refs.insert(
                reference.name.clone(),
                MetamReferenceEntry {
                    ref_len: reference.length,
                    asm: Arc::<str>::from(asm),
                    acc: Arc::<str>::from(acc),
                    profile_taxid,
                    profile_rank: taxonomy.profile_rank,
                    reference_id,
                    assembly_id,
                    sampling_scale,
                    asm_len,
                    flex_len,
                    is_marker: reference.reference_type == ReferenceType::Markers,
                },
            );
        }

        let mut read_map: StableReadMap<AlignmentEntry> = StableReadMap::default();
        let mut num_skipped_reference = 0usize;
        let mut num_alignments = 0usize;
        let mut read_lengths = LengthHistogram::default();
        let mut lone_mates = false;
        let load_progress_options = self.cfg.mapping_load_progress_options();
        let load_progress = progress::counter(load_progress_options, "Loading mapping records");
        let parse_result = feed(&mut |mapping| {
            num_alignments += 1;
            progress::tick_count(
                &load_progress,
                "Loading mapping records",
                num_alignments as u64,
                load_progress_options.interval(),
            );
            if self.cfg.is_paired && !mapping.is_paired_evidence() {
                return Ok(());
            }
            if has_marker_genomes {
                read_lengths.add(mapping.read_len);
                lone_mates |= mapping.is_mate_beyond_reference_end();
            }
            match self.entry_from_metam_mapping(mapping, &metam_refs, fixed_flank_len) {
                Some((read_id, entry)) => {
                    push_string_keyed_value(&mut read_map, read_id, entry);
                }
                None => {
                    num_skipped_reference += 1;
                }
            }
            Ok(())
        });
        load_progress.finish_and_clear();
        let parse_stats = parse_result?;
        if parse_stats.skipped_malformed > 0 {
            log::warn!(
                target: "PROFILER",
                "Skipped {} malformed mapping records from {}.",
                parse_stats.skipped_malformed,
                source
            );
        }

        let aligned_reads = read_map.len();
        let mapped_unit = if self.cfg.is_paired {
            "fragments"
        } else {
            "reads"
        };
        let profile = self.load_metam_profile(parse_stats.profile)?;
        let total_reads = total_input_units(&profile, self.cfg.is_paired).max(aligned_reads);
        let read_alignments = read_map
            .into_iter()
            .map(|(read_name, alignments)| ReadBundle {
                read_name: Arc::<str>::from(read_name),
                alignments,
            })
            .collect::<Vec<_>>();
        self.lower_identity_to_alignments(&read_alignments);
        let marker_model = if has_marker_genomes {
            let geometry =
                sample_fragment_geometry(&self.cfg, &read_alignments, &read_lengths, lone_mates);
            log::info!(target: "PROFILER", "Marker genes: {geometry}.");
            MarkerModel::from_database(db, |genome| assembly_ids[genome], geometry).map(Arc::new)
        } else {
            None
        };
        let unclassified_species = self.cfg.estimate_unclassified.then(|| {
            unclassified::marker_species(db.references.iter().zip(&reference_taxonomies).map(
                |(reference, taxonomy)| {
                    let taxid = taxonomy.profile_taxid;
                    ReferenceKind {
                        taxid: species_normalize.get(&taxid).copied().unwrap_or(taxid),
                        assembly_id: assembly_ids[reference.genome.as_str()],
                        reference_type: reference.reference_type,
                        declared_genome_size: reference.declared_genome_size,
                    }
                },
            ))
        });

        let mapping_rate = if total_reads == 0 {
            0.0
        } else {
            aligned_reads as f64 / total_reads as f64
        };
        // The low-map gate asks how much of the sample the reference explains, so
        // fragments count as the genome share they represent.
        let gate_rate =
            genome_equivalent_mapping_rate(&read_alignments, total_reads, marker_model.as_deref());
        if gate_rate != mapping_rate {
            log::info!(
                target: "PROFILER",
                "Genome-equivalent mapping rate for the low-map gate: {:.5}.",
                gate_rate
            );
        }
        let gate_active = fallback_gate_active(&self.cfg, gate_rate);
        let fallback_mode = genus_fallback_mode_active(&self.cfg, gate_active, db_is_subsampled);
        let chunk_plan = plan_chunk_breadth(
            &self.cfg,
            gate_rate,
            total_reads,
            aligned_reads,
            db_is_subsampled,
            gate_active,
        );
        let min_breadth = self.cfg.breadth.unwrap_or(0.0);
        let chunk_thresholds =
            resolve_community_chunk_breadth(&self.cfg, chunk_plan.estimate_reads);
        let min_chunk_breadth_value = chunk_thresholds.standard;
        let chunk_breadth_source =
            community_chunk_breadth_source(&self.cfg, chunk_plan, chunk_thresholds);

        log::info!(target: "PROFILER", "Total number of reads: {}", total_reads);
        log::info!(target: "PROFILER", "Number of mapping records: {}", num_alignments);
        log::info!(
            target: "PROFILER",
            "Number of mapped {}: {}",
            mapped_unit,
            aligned_reads
        );
        log::info!(target: "PROFILER", "Metam reported mapped read records: {}", profile.mapped_reads);
        log::info!(target: "PROFILER", "Metam emitted mapping records: {}", profile.emitted_records);
        if self.cfg.lowbiomass && self.cfg.chunk_breadth.is_none() {
            log::info!(
                target: "PROFILER",
                "Automatic minimum chunk breadth filter disabled by --low-biomass."
            );
        } else if self.cfg.chunk_breadth.is_none() {
            if chunk_plan.use_aligned_basis {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from mapped {} ({}).",
                    mapped_unit,
                    aligned_reads
                );
            } else {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from total reads ({}).",
                    total_reads
                );
            }
        }
        log_chunk_breadth_scaling(chunk_thresholds);
        log_resolved_profiler_parameters(
            &self.cfg,
            ResolvedProfilerParameters {
                min_breadth,
                min_chunk_breadth: min_chunk_breadth_value,
                chunk_breadth_source,
                chunk_breadth_basis_reads: chunk_plan.estimate_reads,
                total_reads,
                mapped_reads: aligned_reads,
                mapping_rate,
            },
        );

        Ok(MappingData {
            read_alignments,
            asm_lengths: Arc::new(asm_lengths_raw),
            min_breadth,
            min_chunk_breadth: min_chunk_breadth_value,
            safeguard_min_chunk_breadth: chunk_thresholds.safeguard,
            has_scaled_genomes: genome_sampling_scales.values().any(|scale| *scale > 1.0),
            total_reads,
            aligned_reads,
            mapping_rate,
            fallback_mode,
            species_to_genus,
            terminal_ranks: Arc::new(terminal_ranks),
            skipped_reference_at_load: num_skipped_reference,
            marker_model,
            genus_rollup,
            unclassified_species,
            assembly_names,
        })
    }

    fn entry_from_metam_mapping<'a>(
        &self,
        mapping: MetamMappingRecord<'a>,
        metam_refs: &HashMap<String, MetamReferenceEntry>,
        fixed_flank_len: usize,
    ) -> Option<(&'a str, AlignmentEntry)> {
        if self.cfg.is_paired && !mapping.is_paired_evidence() {
            return None;
        }
        let reference = metam_refs.get(mapping.ref_id)?;
        if reference.ref_len == 0 || mapping.read_len == 0 || mapping.mapped_bases == 0 {
            return None;
        }
        let start = mapping.start_1based.clamp(1, reference.ref_len.max(1));
        let mut end = mapping.end_1based_exclusive.max(start.saturating_add(1));
        end = end.min(reference.ref_len.saturating_add(1));
        let midpoint = start + ((end.saturating_sub(start)) / 2);
        let fixed_chunk_start = start.saturating_sub(fixed_flank_len).max(1);
        let fixed_chunk_end = end
            .saturating_add(fixed_flank_len)
            .min(reference.ref_len.saturating_add(1));
        let flex_chunk_start = midpoint.saturating_sub(reference.flex_len).max(1);
        let flex_chunk_end = midpoint
            .saturating_add(reference.flex_len)
            .min(reference.ref_len.saturating_add(1));
        let fraction = (mapping.mapped_bases as f64 / mapping.read_len as f64).min(1.0);
        let read_name: &'a str = if self.cfg.is_paired {
            strip_pair_suffix(mapping.read_id)
        } else {
            mapping.read_id
        };
        Some((
            read_name,
            AlignmentEntry {
                asm: Arc::clone(&reference.asm),
                acc: Arc::clone(&reference.acc),
                profile_taxid: reference.profile_taxid,
                profile_rank: reference.profile_rank,
                reference_id: reference.reference_id,
                assembly_id: reference.assembly_id,
                role: if self.cfg.is_paired {
                    match mapping.mate_number {
                        1 => MappingRole::Mate1,
                        2 => MappingRole::Mate2,
                        _ => return None,
                    }
                } else {
                    MappingRole::Unpaired
                },
                pair_id: mapping.pair_id,
                candidate_boundary: mapping.candidate_boundary,
                sampling_scale: reference.sampling_scale,
                asm_len: reference.asm_len,
                identity: mapping.identity as f32,
                fraction: fraction as f32,
                mapped_len: clamp_u32(mapping.mapped_bases),
                read_len: clamp_u32(mapping.read_len),
                span: (clamp_u32(start), clamp_u32(end)),
                fixed_chunk: (clamp_u32(fixed_chunk_start), clamp_u32(fixed_chunk_end)),
                flex_chunk: (clamp_u32(flex_chunk_start), clamp_u32(flex_chunk_end)),
                marker_gene_len: if reference.is_marker {
                    clamp_u32(reference.ref_len)
                } else {
                    0
                },
            },
        ))
    }

    /// Precompute taxonomy lookups used during per-read assignment.
    ///
    /// Given the set of raw taxids of the database references, this produces:
    ///
    /// - `species_normalize`: sparse `raw -> species` map. Only contains
    ///   entries where the raw taxid is **below** species rank and a
    ///   species ancestor exists. Entries where raw is already species
    ///   (or has no species ancestor) are omitted; callers treat a
    ///   missing entry as "use the raw taxid unchanged". This is the
    ///   mechanism that merges multiple sub-species / strains belonging
    ///   to the same species into a single species-level assignment.
    ///
    /// - `species_to_genus`: `species -> genus` map, keyed by the
    ///   **post-normalization** species taxid so that the genus
    ///   fallback path looks up the genus once per species regardless
    ///   of how many sub-taxon flavours share that species.
    ///
    /// Taxids without a resolvable genus are simply absent from
    /// `species_to_genus`; the genus fallback silently skips them.
    fn build_taxonomy_lookups(
        &self,
        raw_taxids: &HashSet<u32>,
    ) -> (HashMap<u32, u32>, HashMap<u32, u32>) {
        let tax = self.taxonomy.as_ref();
        let mut species_normalize: HashMap<u32, u32> = HashMap::new();
        let mut species_to_genus: HashMap<u32, u32> = HashMap::with_capacity(raw_taxids.len());
        for &raw in raw_taxids {
            // Determine the species anchor (or fall back to the raw id
            // when no species ancestor exists, to preserve behaviour for
            // references that already point above species rank).
            let species = match resolve_species_taxid(tax, raw) {
                Some(s) => {
                    if s != raw {
                        species_normalize.insert(raw, s);
                    }
                    s
                }
                None => raw,
            };
            if let Entry::Vacant(slot) = species_to_genus.entry(species) {
                if let Some(g) = resolve_genus_taxid(tax, species) {
                    slot.insert(g);
                }
            }
        }
        (species_normalize, species_to_genus)
    }

    /// Run the main profiling algorithm on loaded alignment data.
    ///
    /// Steps:
    /// 1. Assign taxonomy to reads based on alignment quality
    /// 2. Compute coverage metrics for each taxon
    /// 3. Resolve ambiguous multi-mapping reads with EM
    /// 4. Filter taxa based on coverage statistics
    /// 5. Write output files
    fn profile(&self, mut mapping_data: MappingData) -> Result<()> {
        let MappingData {
            read_alignments,
            fallback_mode,
            species_to_genus,
            mapping_rate,
            total_reads,
            aligned_reads,
            skipped_reference_at_load,
            marker_model,
            ..
        } = &mapping_data;
        let verbose_logging = self.cfg.verbose || self.cfg.very_verbose;
        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "Assigning taxonomy for reads with {} threads ...",
                self.cfg.threads
            );
            log::info!(
                target: "PROFILER",
                "Run mapping rate carried into assignment: {:.5} ({}/{}). Genus fallback = {}.",
                mapping_rate,
                aligned_reads,
                total_reads,
                if *fallback_mode { "enabled" } else { "disabled" },
            );
        }
        // Derived from the final configuration, after any identity check.
        let evidence_policy = EvidencePolicy::from_config(&self.cfg).ok_or_else(|| {
            anyhow!("community profiling takes Illumina reads or single-end long reads")
        })?;
        let mut assignment = self.assign_taxonomy(
            read_alignments,
            *fallback_mode,
            Arc::clone(species_to_genus),
            evidence_policy,
            marker_model.as_deref(),
        )?;
        // The `skipped_reference` diagnostic counts the mapping records
        // rejected at load (unparsable species taxid or missing assembly
        // length).
        assignment.filter_stats.skipped_reference += *skipped_reference_at_load;
        if verbose_logging {
            log::info!(target: "PROFILER", "Assigning taxonomy finished!");
        }
        let candidate_reads = assignment.read_assignments.len();
        if assignment.boundary_ambiguous_fragments > 0 {
            log::info!(
                target: "PROFILER",
                "Candidate-boundary safeguard kept {} {} in abundance ambiguity but excluded them from direct coordinate evidence.",
                assignment.boundary_ambiguous_fragments,
                if self.cfg.is_paired { "fragments" } else { "reads" },
            );
        }
        let assignment_unit = if self.cfg.is_paired {
            "fragments"
        } else {
            "reads"
        };
        // Assignment owns every value required by EM and presence filtering.
        // Release the dominant all-mappings buffer before coverage reduction
        // and output construction so mapper/profiler peak RSS does not carry
        // both representations at once.
        drop(std::mem::take(&mut mapping_data.read_alignments));
        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "{} with candidate taxon assignments before coverage filtering: {}",
                if self.cfg.is_paired { "Fragments" } else { "Reads" },
                candidate_reads
            );
            log::info!(
                target: "PROFILER",
                "Profiling taxonomy abundance ..."
            );
        }
        let final_classified = self.write_outputs(mapping_data, assignment, evidence_policy)?;
        log::info!(
            target: "PROFILER",
            "Number of {} assigned to retained taxa: {}",
            assignment_unit,
            final_classified
        );
        log::info!(target: "PROFILER", "Taxonomy profiling finished.");
        Ok(())
    }

    /// Assign taxonomy to reads based on alignment quality.
    ///
    /// Uses rayon for parallel processing of read batches. Each read's
    /// alignments are filtered by identity, mapped length, and fraction
    /// thresholds, then the best alignment per species is selected.
    ///
    /// # Algorithm
    ///
    /// For each read:
    /// 1. Filter alignments by quality thresholds
    /// 2. Parse reference names to extract taxonomy info
    /// 3. For each species, keep only the highest-identity alignment
    /// 4. Record coverage intervals for breadth calculation
    /// 5. Store taxonomic assignment (may be multi-species)
    ///
    /// # Returns
    ///
    /// AssignmentResult containing per-read assignments and coverage intervals
    fn assign_taxonomy(
        &self,
        read_alignments: &[ReadBundle],
        fallback_mode: bool,
        species_to_genus: Arc<HashMap<u32, u32>>,
        evidence_policy: EvidencePolicy,
        markers: Option<&MarkerModel>,
    ) -> Result<AssignmentResult> {
        let batch_size = self.cfg.batch_size.max(1);
        let cfg = Arc::new(self.cfg.clone());

        // Create rayon thread pool for parallel batch processing
        let pool = ThreadPoolBuilder::new()
            .num_threads(self.cfg.threads.max(1))
            .build()?;

        let pb = progress::bar(
            self.cfg.progress,
            read_alignments.len() as u64,
            "Assigning taxonomy",
        );
        let progress = pb.clone();

        // Process read batches in parallel
        let mut results: Vec<AssignmentResult> = pool.install(|| {
            read_alignments
                .par_chunks(batch_size)
                .map_with(progress, |progress_bar, chunk| {
                    let result = process_evidence_batch(
                        chunk,
                        &cfg,
                        evidence_policy,
                        fallback_mode,
                        species_to_genus.as_ref(),
                        markers,
                    );
                    progress_bar.inc(chunk.len() as u64);
                    result
                })
                .collect()
        });
        pb.finish_and_clear();

        // Aggregate results from all batches
        let mut final_result = AssignmentResult::default();
        for result in results.drain(..) {
            // Merge read assignments
            final_result
                .read_assignments
                .extend(result.read_assignments);

            // Merge coverage intervals
            for (key, intervals) in result.acc_intervals {
                final_result
                    .acc_intervals
                    .entry(key)
                    .or_default()
                    .extend(intervals);
            }
            for (key, intervals) in result.acc_fixed_chunk_intervals {
                final_result
                    .acc_fixed_chunk_intervals
                    .entry(key)
                    .or_default()
                    .extend(intervals);
            }
            for (key, intervals) in result.acc_flex_chunk_intervals {
                final_result
                    .acc_flex_chunk_intervals
                    .entry(key)
                    .or_default()
                    .extend(intervals);
            }
            for (key, presence) in result.presence {
                final_result
                    .presence
                    .entry(key)
                    .or_default()
                    .merge(presence);
            }
            for (assembly_id, reads) in result.marker_genome_reads {
                final_result
                    .marker_genome_reads
                    .entry(assembly_id)
                    .or_default()
                    .merge(reads);
            }
            for (taxid, evidence) in result.taxon_evidence {
                let target = final_result.taxon_evidence.entry(taxid).or_default();
                target.candidate_fragments = target
                    .candidate_fragments
                    .saturating_add(evidence.candidate_fragments);
                target.resolved_candidate_fragments = target
                    .resolved_candidate_fragments
                    .saturating_add(evidence.resolved_candidate_fragments);
                target.best_fragments = target
                    .best_fragments
                    .saturating_add(evidence.best_fragments);
                target.clear_best_fragments = target
                    .clear_best_fragments
                    .saturating_add(evidence.clear_best_fragments);
                target.mapped_fragments = target
                    .mapped_fragments
                    .saturating_add(evidence.mapped_fragments);
                target.identity_sum += evidence.identity_sum;
            }
            for (taxid, placements) in result.best_placements {
                final_result
                    .best_placements
                    .entry(taxid)
                    .or_default()
                    .add(placements);
            }
            final_result.exact_ties.append(result.exact_ties);
            final_result
                .uniquely_supported_taxa
                .extend(result.uniquely_supported_taxa);
            final_result.boundary_ambiguous_fragments = final_result
                .boundary_ambiguous_fragments
                .saturating_add(result.boundary_ambiguous_fragments);

            // Merge filter statistics
            final_result.filter_stats.merge(result.filter_stats);
        }
        Ok(final_result)
    }

    /// Write profiling results to output files.
    ///
    /// This is a large function that:
    /// 1. Computes final coverage metrics per terminal taxon
    /// 2. Performs EM algorithm for ambiguous read redistribution
    /// 3. Calculates abundance estimates and statistical p-values
    /// 4. Filters taxa based on coverage criteria
    /// 5. Writes profile.txt and, when requested, classify.txt output files
    ///
    /// # Output Files
    ///
    /// - `<prefix>.profile.txt`: Filtered abundance profile with statistics
    /// - `<prefix>.rprofile.txt`: Post-EM profile before taxon-level evidence
    ///   and reporting filters (if keep_raw=true)
    /// - `<prefix>.classify.txt`: Optional per-read assignments
    ///   (`--save-classification`)
    ///
    /// # EM Algorithm
    ///
    /// For multi-mapping reads, the EM algorithm iteratively:
    /// 1. E-step: Redistribute reads based on current abundance estimates
    /// 2. M-step: Recalculate abundances from redistributed reads
    /// 3. Converge when relative change < 1e-5
    ///
    /// # Returns
    ///
    /// Number of taxonomically classified reads
    fn write_outputs(
        &self,
        mapping_data: MappingData,
        mut assignment: AssignmentResult,
        evidence_policy: EvidencePolicy,
    ) -> Result<usize> {
        let MappingData {
            asm_lengths,
            min_breadth,
            min_chunk_breadth,
            safeguard_min_chunk_breadth,
            has_scaled_genomes,
            terminal_ranks,
            marker_model,
            genus_rollup,
            total_reads,
            unclassified_species,
            assembly_names,
            ..
        } = mapping_data;
        let verbose_logging = self.cfg.verbose || self.cfg.very_verbose;
        // Majority-LCA rows are evaluated later than native coordinate
        // witnesses. Resolve the same platform mode here so their default
        // ratio bounds come from the standard evidence policy rather than the
        // older output-only genus/species defaults.
        let lca_ratio_mode = evidence_policy.mode;

        let mut taxon_coverage: HashMap<u32, TaxonCoverage> = HashMap::new();
        let mut coverage_keys = assignment.acc_intervals.keys().collect::<Vec<_>>();
        coverage_keys.sort_unstable();
        for key in coverage_keys {
            let intervals = &assignment.acc_intervals[key];
            let merged = merge_intervals(intervals);
            let merged_fixed = merge_intervals(
                assignment
                    .acc_fixed_chunk_intervals
                    .get(key)
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]),
            );
            let merged_flex = merge_intervals(
                assignment
                    .acc_flex_chunk_intervals
                    .get(key)
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]),
            );
            let asm_len = asm_lengths
                .get(key.asm.as_ref())
                .copied()
                .unwrap_or(1.0)
                .max(1.0);
            let total_span: usize = merged.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let total_fixed: usize = merged_fixed.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let total_flex: usize = merged_flex.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
            let coverage = taxon_coverage.entry(key.taxid).or_default();
            coverage.breadth += total_span as f64 / asm_len;
            coverage.fixed_chunk += total_fixed as f64 / asm_len;
            coverage.flex_chunk += total_flex as f64 / asm_len;
        }

        let resolved = resolve_coordinate_evidence(
            &self.cfg,
            evidence_policy,
            &mut assignment,
            min_breadth,
            min_chunk_breadth,
            safeguard_min_chunk_breadth,
            CoverageModel {
                markers: marker_model.as_deref(),
                layout: genus_rollup.as_deref().map(|rollup| &rollup.layout),
            },
        )?;
        rescale_to_best_fit(
            &mut assignment.read_assignments,
            &resolved.best_fit,
            marker_model.as_deref(),
        );
        let best_fit = resolved.best_fit;
        let mut evidence_decisions = resolved.decisions;
        let mut evidence_selected_taxa = resolved.selected_taxa;
        if let Some(rollup) = genus_rollup.as_deref() {
            let even = rollup.even_coverage_species(
                &evidence_decisions,
                |taxid| terminal_ranks.get(&taxid).copied(),
                min_chunk_breadth,
            );
            log::info!(
                target: "PROFILER",
                "Species reported for even coverage of their genome: {}.",
                even.len()
            );
            for taxid in even {
                if let Some(decision) = evidence_decisions.get_mut(&taxid) {
                    decision.rescues.even_coverage = true;
                }
                evidence_selected_taxa.insert(taxid);
                if self.cfg.verbose || self.cfg.very_verbose {
                    log::info!(
                        target: "PROFILER",
                        "  {} ({}) has even coverage.",
                        self.taxonomy.get_name(taxid).unwrap_or("unknown"),
                        taxid,
                    );
                }
            }
        }

        // Preserve only compact coordinate measurements for the optional raw
        // profile. Its separate EM runs later, after every production value is
        // finalized, so diagnostic work cannot influence production
        // aggregation or filtering. The large coordinate stores are never
        // cloned.
        let raw_profile_coverage = self.cfg.keep_raw.then(|| taxon_coverage.clone());

        let classification_path = format!("{}.classify.txt", self.cfg.outprefix);
        let database_is_fractional = self
            .cfg
            .db_metadata
            .as_ref()
            .map(|metadata| metadata.is_fractional)
            .unwrap_or(has_scaled_genomes);
        let write_classifications = self.cfg.save_classifications && !database_is_fractional;
        let mut resolved_classifications = ResolvedClassifications::new(write_classifications);
        if self.cfg.save_classifications {
            if database_is_fractional {
                log::warn!(
                    target: "PROFILER",
                    "Cannot write classify.txt for a fractional database because per-read assignments cannot be expanded safely to full-index equivalents."
                );
            } else {
                log::info!(
                    target: "PROFILER",
                    "Writing per-read classifications to {}",
                    classification_path
                );
            }
        }

        let mut unambiguous: HashMap<u32, TaxonTotals> = HashMap::new();
        let mut ambiguous_counts: HashMap<Vec<u32>, AmbiguousTotals> = HashMap::new();
        let mut ambiguous_reads: HashMap<Vec<u32>, Vec<String>> = HashMap::new();
        let mut orphans = genus_rollup
            .as_deref()
            .map(|rollup| OrphanedFragments::new(rollup, write_classifications));

        for assignment in &assignment.read_assignments {
            let mut final_taxids = Vec::new();
            let mut final_depths = Vec::new();
            let mut final_maplens = Vec::new();
            let mut final_scales = Vec::new();
            for (((taxid, depth), maplen), sampling_scale) in assignment
                .taxids
                .iter()
                .zip(&assignment.depths)
                .zip(&assignment.mapped_lens)
                .zip(&assignment.sampling_scales)
            {
                let coverage = taxon_coverage.get(taxid).cloned().unwrap_or_default();
                let evidence_selected =
                    assignment.requires_terminal_evidence && evidence_selected_taxa.contains(taxid);
                let legacy_selected = !assignment.requires_terminal_evidence;
                if evidence_selected || (legacy_selected && coverage.breadth >= min_breadth) {
                    final_taxids.push(*taxid);
                    final_depths.push(*depth);
                    final_maplens.push(*maplen);
                    final_scales.push(*sampling_scale);
                }
            }
            match final_taxids.len() {
                0 => {
                    // The coverage evidence rejected every candidate.
                    if let Some(orphans) = orphans
                        .as_mut()
                        .filter(|_| assignment.requires_terminal_evidence)
                    {
                        let best = assignment.best as usize;
                        orphans.offer(
                            &assignment.taxids,
                            best,
                            assignment.depths[best],
                            assignment.mapped_lens[best],
                            assignment.sampling_scales[best],
                            &assignment.read,
                        );
                    }
                    continue;
                }
                1 => {
                    let taxid = final_taxids[0];
                    let entry = unambiguous.entry(taxid).or_default();
                    entry.support.add_weighted(1.0, final_scales[0], 1.0);
                    entry.depth += final_depths[0];
                    entry.mapped_len += final_maplens[0];
                    resolved_classifications.record_read(taxid, &assignment.read);
                }
                _ => {
                    if write_classifications {
                        ambiguous_reads
                            .entry(final_taxids.clone())
                            .or_default()
                            .push(assignment.read.to_string());
                    }
                    ambiguous_counts
                        .entry(final_taxids.clone())
                        .and_modify(|entry| {
                            entry.observed_count += 1.0;
                            accumulate(&mut entry.depths, &final_depths);
                            accumulate(&mut entry.mapped_lens, &final_maplens);
                            accumulate(&mut entry.equivalent_support_sums, &final_scales);
                        })
                        .or_insert_with(|| AmbiguousTotals {
                            observed_count: 1.0,
                            depths: final_depths.clone(),
                            mapped_lens: final_maplens.clone(),
                            equivalent_support_sums: final_scales.clone(),
                        });
                }
            }
        }

        let direct_taxa_before_lca: HashSet<u32> = unambiguous.keys().copied().collect();
        let mut taxids_with_lca: Vec<Vec<u32>> = Vec::new();
        let mut lca_ratio_lower_bounds: HashMap<u32, (f64, f64)> = HashMap::new();
        // Multiple ambiguous candidate groups may fold into the same LCA. Keep
        // one globally strongest coherent descendant witness for that LCA
        // instead of letting HashMap iteration order choose the last group.
        let mut coherent_lca_witnesses: HashMap<u32, (u32, WitnessDecision)> = HashMap::new();
        // Groups of genus-fallback reads have no evidence decisions. Keep their
        // per-field maxima in a separate map so an LCA group cannot overwrite
        // coverage already accumulated for the taxon itself.
        let mut legacy_lca_coverage: HashMap<u32, TaxonCoverage> = HashMap::new();
        let mut taxids_with_any_ambiguous: Vec<Vec<u32>> = Vec::new();
        let mut ambiguous_removed: HashMap<Vec<u32>, AmbiguousTotals> = HashMap::new();

        // Iterate compact candidate groups in lexical taxid order. The values
        // stay borrowed, so reproducibility costs one pointer vector rather
        // than cloning every candidate/value vector.
        let mut initial_group_keys = ambiguous_counts.keys().collect::<Vec<_>>();
        initial_group_keys.sort_unstable();
        for taxids in initial_group_keys {
            let totals = &ambiguous_counts[taxids];
            let presence: Vec<bool> = taxids
                .iter()
                .map(|taxid| unambiguous.contains_key(taxid))
                .collect();
            if presence.iter().all(|&p| p) {
                continue;
            } else if presence.iter().any(|&p| p) {
                taxids_with_any_ambiguous.push(taxids.clone());
                let mut filtered_taxids = Vec::new();
                let mut filtered_depths = Vec::new();
                let mut filtered_maplens = Vec::new();
                let mut filtered_scale_sums = Vec::new();
                for (idx, taxid) in taxids.iter().enumerate() {
                    if presence[idx] {
                        filtered_taxids.push(*taxid);
                        filtered_depths.push(totals.depths[idx]);
                        filtered_maplens.push(totals.mapped_lens[idx]);
                        filtered_scale_sums.push(totals.equivalent_support_sums[idx]);
                    }
                }
                if filtered_taxids.len() == 1 {
                    let taxid = filtered_taxids[0];
                    let entry = unambiguous.entry(taxid).or_default();
                    entry
                        .support
                        .add_weighted(totals.observed_count, filtered_scale_sums[0], 1.0);
                    entry.depth += filtered_depths[0];
                    entry.mapped_len += filtered_maplens[0];
                    resolved_classifications.record_direct_group(
                        taxid,
                        totals.observed_count.round() as usize,
                        ambiguous_reads.remove(taxids),
                    );
                } else if !filtered_taxids.is_empty() {
                    ambiguous_removed
                        .entry(filtered_taxids.clone())
                        .and_modify(|entry| {
                            entry.observed_count += totals.observed_count;
                            accumulate(&mut entry.depths, &filtered_depths);
                            accumulate(&mut entry.mapped_lens, &filtered_maplens);
                            accumulate(&mut entry.equivalent_support_sums, &filtered_scale_sums);
                        })
                        .or_insert_with(|| AmbiguousTotals {
                            observed_count: totals.observed_count,
                            depths: filtered_depths.clone(),
                            mapped_lens: filtered_maplens.clone(),
                            equivalent_support_sums: filtered_scale_sums.clone(),
                        });
                    if let Some(mut reads) = ambiguous_reads.remove(taxids) {
                        ambiguous_reads
                            .entry(filtered_taxids)
                            .or_default()
                            .append(&mut reads);
                    }
                }
            } else {
                taxids_with_lca.push(taxids.clone());
                let taxid_set: HashSet<u32> = taxids.iter().copied().collect();
                if let Some((lca_taxid, _lca_name, lca_rank)) =
                    self.taxonomy.get_majority_lca(&taxid_set, 0.7)
                {
                    if let Some(bounds) = standard_lca_ratio_lower_bounds(&lca_rank, lca_ratio_mode)
                    {
                        lca_ratio_lower_bounds.insert(lca_taxid, bounds);
                    }
                    // Carry the coverage of a single coherent descendant
                    // witness. Maximizing breadth, fixed chunk and flex chunk
                    // independently could take each from a different genome,
                    // and the reporting path then divided those by an
                    // expectation rebuilt from the group's pooled depth and
                    // count -- so OEBR/COEBR compared one genome's observation
                    // against another population's expectation. Choosing one
                    // witness keeps observed, expected and both p-values on the
                    // same coordinate system.
                    let coherent_witness = coherent_lca_witness(taxids, &evidence_decisions);
                    if let Some(witness) = coherent_witness {
                        retain_coherent_lca_witness(
                            &mut coherent_lca_witnesses,
                            lca_taxid,
                            witness,
                        );
                    } else {
                        // No member has coordinate evidence (a group of
                        // genus-fallback reads). Retain the maxima over their
                        // pooled coverage.
                        let pick = |select: fn(&TaxonCoverage) -> f64| {
                            taxids
                                .iter()
                                .filter_map(|taxid| taxon_coverage.get(taxid).map(select))
                                .fold(0.0, f64::max)
                        };
                        let legacy = TaxonCoverage {
                            breadth: pick(|cov| cov.breadth),
                            fixed_chunk: pick(|cov| cov.fixed_chunk),
                            flex_chunk: pick(|cov| cov.flex_chunk),
                        };
                        merge_legacy_lca_coverage(&mut legacy_lca_coverage, lca_taxid, legacy);
                    }
                    let depth_mean = totals.depths.iter().sum::<f64>() / totals.depths.len() as f64;
                    let maplen_mean =
                        totals.mapped_lens.iter().sum::<f64>() / totals.mapped_lens.len() as f64;
                    let equivalent_mean = totals.equivalent_support_sums.iter().sum::<f64>()
                        / totals.equivalent_support_sums.len() as f64;
                    let entry = unambiguous.entry(lca_taxid).or_default();
                    entry
                        .support
                        .add_weighted(totals.observed_count, equivalent_mean, 1.0);
                    entry.depth += depth_mean;
                    entry.mapped_len += maplen_mean;
                    resolved_classifications.record_group(
                        lca_taxid,
                        taxids.clone(),
                        totals.observed_count.round() as usize,
                        ambiguous_reads.remove(taxids),
                    );
                }
            }
        }

        let lca_metrics: HashMap<u32, EvidenceCoverage> = coherent_lca_witnesses
            .into_iter()
            .map(|(lca_taxid, (_, decision))| (lca_taxid, decision.metrics))
            .collect();

        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "Majority-LCA folding: {} ambiguous candidate group(s) folded; {} distinct LCA taxa carried a coherent descendant witness.",
                taxids_with_lca.len(),
                lca_metrics.len(),
            );
        }
        for taxids in taxids_with_lca {
            ambiguous_counts.remove(&taxids);
            ambiguous_reads.remove(&taxids);
        }
        for taxids in taxids_with_any_ambiguous {
            ambiguous_counts.remove(&taxids);
            ambiguous_reads.remove(&taxids);
        }
        let mut ambiguous_removed = ambiguous_removed.into_iter().collect::<Vec<_>>();
        ambiguous_removed.sort_unstable_by(|left, right| left.0.cmp(&right.0));
        for (taxids, values) in ambiguous_removed {
            ambiguous_counts
                .entry(taxids)
                .and_modify(|entry| {
                    entry.observed_count += values.observed_count;
                    accumulate(&mut entry.depths, &values.depths);
                    accumulate(&mut entry.mapped_lens, &values.mapped_lens);
                    accumulate(
                        &mut entry.equivalent_support_sums,
                        &values.equivalent_support_sums,
                    );
                })
                .or_insert(values);
        }

        let rollup_metrics = match orphans {
            Some(orphans) => self.add_genus_rows(
                orphans.into_genus_rows(&evidence_decisions, &evidence_selected_taxa),
                &mut unambiguous,
                &mut resolved_classifications,
            ),
            None => HashMap::new(),
        };
        let rolled_up_genera: HashSet<u32> = rollup_metrics.keys().copied().collect();

        let mut em_groups = ambiguous_counts.iter().collect::<Vec<_>>();
        em_groups.sort_unstable_by(|left, right| left.0.cmp(right.0));
        let mut em_taxids = unambiguous.keys().copied().collect::<Vec<_>>();
        em_taxids.sort_unstable();
        let em = shared_read_em(
            &unambiguous,
            &em_groups,
            self.cfg.max_em_iterations,
            "Community",
            verbose_logging,
        );
        let mut ambiguous_fraction = group_shares(&em_groups, &em.previous);
        let species_count_depth = em.totals;

        let norm_factor: f64 = 0.01
            * em_taxids
                .iter()
                .map(|taxid| species_count_depth.get(taxid).map_or(0.0, |totals| totals.depth))
                .sum::<f64>();
        let num_taxa = species_count_depth.len();
        if verbose_logging {
            log::info!(
                target: "PROFILER",
                "Legacy/LCA p-value hypothesis family: {} post-EM taxa.",
                num_taxa,
            );
        }

        let mut out_taxa_list = Vec::new();
        let host_taxid = self.cfg.host.as_deref().and_then(|h| h.parse::<u32>().ok());
        let pathogen_map = if self.cfg.host.is_some() {
            self.cfg
                .pathogen_host
                .as_deref()
                .map(load_pathogen_table)
                .transpose()?
        } else {
            None
        };
        let mut entries: Vec<(u32, TaxonTotals)> = species_count_depth
            .iter()
            .map(|(taxid, totals)| (*taxid, *totals))
            .collect();
        entries.sort_by(|left, right| {
            right
                .1
                .depth
                .total_cmp(&left.1.depth)
                .then_with(|| left.0.cmp(&right.0))
        });

        for (taxid, totals) in entries {
            let mapcount = totals.support.observed;
            if mapcount < 0.5 {
                continue;
            }
            let equivalent_mapcount = totals.support.equivalent;
            let depth = totals.depth;
            let maplen = totals.mapped_len;
            let abundance = if norm_factor <= f64::EPSILON {
                0.0
            } else {
                depth / norm_factor
            };
            let rank = terminal_ranks
                .get(&taxid)
                .copied()
                .map(ProfileRank::as_str)
                .or_else(|| self.taxonomy.get_rank(taxid))
                .unwrap_or("no rank")
                .to_string();
            let taxname = self
                .taxonomy
                .get_name(taxid)
                .unwrap_or("Taxonomy deprecated");
            if taxname == "Taxonomy deprecated" {
                continue;
            }
            let evidence_metrics = if evidence_selected_taxa.contains(&taxid) {
                evidence_decisions
                    .get(&taxid)
                    .map(|decision| decision.metrics)
            } else {
                // A majority-LCA row reports the coherent coverage of its chosen
                // descendant witness instead of dividing pooled observations by
                // an expectation rebuilt from pooled depth and count. Report
                // filtering still applies to LCA rows; only the numbers change.
                // A roll-up genus row reports its strongest donor likewise.
                lca_metrics
                    .get(&taxid)
                    .or_else(|| rollup_metrics.get(&taxid))
                    .copied()
            };
            let (
                breadth,
                fixed_chunk,
                flex_chunk,
                expected_breadth,
                cov_prob,
                expected_flex_chunk,
                chunkcov_prob,
            ) = if let Some(metrics) = evidence_metrics {
                // Presence statistics are all derived from the same resolved,
                // assembly-coherent coordinate set. EM depth remains solely
                // the abundance estimate and must not redefine its expected
                // breadth or p-values.
                (
                    metrics.breadth,
                    metrics.fixed_chunk,
                    metrics.flex_chunk,
                    metrics.expected_breadth,
                    metrics.breadth_pvalue,
                    metrics.expected_flex_chunk,
                    metrics.chunk_pvalue,
                )
            } else {
                let coverage = legacy_output_coverage(
                    taxid,
                    &direct_taxa_before_lca,
                    &taxon_coverage,
                    &legacy_lca_coverage,
                );
                let breadth = coverage.breadth.min(1.0);
                let fixed_chunk = coverage.fixed_chunk.min(1.0);
                let flex_chunk = coverage.flex_chunk.min(1.0);
                let expected_breadth = crate::profiler::evidence::expected_breadth(depth, mapcount);
                let cov_prob = if rank == "species" {
                    coverage_deviation_pvalue(breadth, expected_breadth, Some(depth), mapcount, num_taxa)
                } else {
                    None
                };
                let genome_len = if depth <= f64::EPSILON {
                    0.0
                } else {
                    maplen / depth
                };
                let num_flex_chunks = if genome_len <= 0.0 {
                    1
                } else {
                    (genome_len.sqrt() / 2.0).floor().max(1.0) as usize
                };
                let expected_flex_chunk = 1.0 - (1.0 - 1.0 / num_flex_chunks as f64).powf(mapcount);
                let chunkcov_prob = if rank == "species" {
                    coverage_deviation_pvalue(flex_chunk, expected_flex_chunk, None, mapcount, num_taxa)
                } else {
                    None
                };
                (
                    breadth,
                    fixed_chunk,
                    flex_chunk,
                    expected_breadth,
                    cov_prob,
                    expected_flex_chunk,
                    chunkcov_prob,
                )
            };
            let Some(annotations) = self.row_annotations(taxid, host_taxid, pathogen_map.as_ref())
            else {
                continue;
            };
            // A taxon kept by the coverage evidence also reports its best-fit
            // genome and the depth of that evidence there before abundance.
            let best_fit_report = best_fit
                .get(&taxid)
                .filter(|_| evidence_selected_taxa.contains(&taxid))
                .zip(evidence_metrics)
                .and_then(|(genome, metrics)| {
                    Some(BestFitReport {
                        genome: Arc::clone(assembly_names.get(genome.assembly_id as usize)?),
                        evidence_depth: metrics.depth,
                    })
                });
            out_taxa_list.push(ProfileRow {
                taxname: taxname.to_string(),
                taxid,
                rank,
                equivalent_mapcount,
                depth,
                abundance,
                breadth,
                expected_breadth,
                breadth_pvalue: cov_prob,
                fixed_chunk,
                flex_chunk,
                expected_flex_chunk,
                chunk_pvalue: chunkcov_prob,
                best_fit: best_fit_report,
                annotations,
            });
        }

        let mut normalized: Vec<ProfileRow> = out_taxa_list
            .into_iter()
            .filter(|row| {
                community_profile_row_passes(
                    row,
                    &self.cfg,
                    min_chunk_breadth,
                    &evidence_selected_taxa,
                    &rolled_up_genera,
                    &lca_ratio_lower_bounds,
                )
            })
            .collect();
        if let Some(gap) = self.cfg.relative_identity_gap {
            let removed = remove_species_explained_by_relatives(
                &mut normalized,
                &self.taxonomy,
                &assignment.taxon_evidence,
                gap,
            );
            log::info!(
                target: "PROFILER",
                "Relative identity guard removed {} species.",
                removed.len()
            );
            if verbose_logging {
                for (taxid, relative) in removed {
                    log::info!(
                        target: "PROFILER",
                        "  {taxid} is explained by {relative}, which has more depth and more identical reads."
                    );
                }
            }
        }
        // Before species move to their genus, while marker species still
        // carry their own scale.
        let unclassified = unclassified_species.as_ref().map(|species| {
            // Long reads differ in length and in the genes each covers, so
            // their scales come from the reads themselves.
            let scales = unclassified::marker_scales(
                species,
                marker_model.as_deref(),
                &assignment.marker_genome_reads,
                self.cfg.coverage_gate_is_or(),
            );
            let estimate = unclassified::estimate(&normalized, &scales, total_reads);
            match estimate {
                Ok(estimate) => log::info!(
                    target: "PROFILER",
                    "Unclassified estimate: {:.2}% of {} input {} are not explained by the profile.",
                    estimate.share * 100.0,
                    total_reads,
                    if self.cfg.is_paired { "fragments" } else { "reads" }
                ),
                Err(taxid) => log::warn!(
                    target: "PROFILER",
                    "Unclassified estimate unavailable: species {taxid} is represented by marker genes that declare no genome size, cannot hold a read of this length, or have no resolved long read."
                ),
            }
            estimate.ok()
        });
        if let Some(margin) = self.cfg.identity_margin {
            let moved = self.report_low_identity_species_at_genus(
                &mut normalized,
                &assignment.best_placements,
                margin,
                host_taxid,
                pathogen_map.as_ref(),
            );
            resolved_classifications.relabel(&moved);
            relabel_ambiguous_fractions(&mut ambiguous_fraction, &moved);
        }
        let total_abundance = normalized.iter().map(|row| row.abundance).sum::<f64>();
        for row in &mut normalized {
            row.normalize_abundance(total_abundance);
        }
        let retained_taxa = normalized
            .iter()
            .map(|row| row.taxid)
            .collect::<HashSet<_>>();
        let final_classified_reads = write_community_classifications(
            write_classifications.then_some(classification_path.as_str()),
            &resolved_classifications,
            &ambiguous_counts,
            &ambiguous_reads,
            &ambiguous_fraction,
            &retained_taxa,
            &terminal_ranks,
            &self.taxonomy,
        )?;

        // Production assignment, EM, report filtering, and abundance
        // normalization are complete. Release their grouped working sets
        // before the opt-in shadow EM to avoid overlapping the two memory
        // peaks. Keep only per-read candidate vectors needed by raw EM, then
        // release the remaining coordinate-heavy AssignmentResult fields.
        drop(ambiguous_counts);
        drop(ambiguous_reads);
        drop(ambiguous_fraction);
        drop(unambiguous);
        drop(species_count_depth);
        let raw_assignments = raw_profile_coverage
            .as_ref()
            .map(|_| std::mem::take(&mut assignment.read_assignments));
        drop(assignment);

        let raw_taxa_list = if let (Some(mut coverage), Some(assignments)) =
            (raw_profile_coverage, raw_assignments)
        {
            let totals = estimate_unfiltered_taxon_totals(
                &assignments,
                &self.taxonomy,
                &mut coverage,
                self.cfg.max_em_iterations,
                verbose_logging,
            )?;
            self.build_raw_profile_rows(
                &RawProfileEstimate { totals, coverage },
                &evidence_decisions,
                &terminal_ranks,
                host_taxid,
                pathogen_map.as_ref(),
            )
        } else {
            Vec::new()
        };

        let out_prefix = if self.cfg.host.is_some() {
            format!("{}.pathogen", self.cfg.outprefix)
        } else {
            self.cfg.outprefix.clone()
        };
        let profile_path = format!("{}.profile.txt", out_prefix);
        let mut profile_writer = WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(File::create(&profile_path)?);
        for row in &normalized {
            row.write(&mut profile_writer)?;
        }
        if let Some(estimate) = unclassified {
            unclassified::write_row(&mut profile_writer, estimate)?;
        }
        profile_writer.flush()?;

        if self.cfg.keep_raw {
            let raw_path = format!("{}.rprofile.txt", out_prefix);
            let mut raw_writer = WriterBuilder::new()
                .has_headers(false)
                .delimiter(b'\t')
                .from_writer(File::create(&raw_path)?);
            for row in raw_taxa_list {
                row.write(&mut raw_writer)?;
            }
            raw_writer.flush()?;
        }

        Ok(final_classified_reads)
    }

    /// Add the genus rows of the roll-up to the unambiguous totals
    /// that enter the EM, and their fragments to the classifications. Returns
    /// the coverage that each row reports.
    fn add_genus_rows(
        &self,
        rows: Vec<GenusRow>,
        unambiguous: &mut HashMap<u32, TaxonTotals>,
        classifications: &mut ResolvedClassifications,
    ) -> HashMap<u32, EvidenceCoverage> {
        let unit = if self.cfg.is_paired {
            "fragments"
        } else {
            "reads"
        };
        let name = |taxid: u32| self.taxonomy.get_name(taxid).unwrap_or("unknown");
        log::info!(
            target: "PROFILER",
            "Reads of removed species reported at genus: {} genus rows, {:.0} {} from {} species.",
            rows.len(),
            // A float sum of no rows is -0.0, printed "-0".
            rows.iter().fold(0.0, |total, row| total + row.fragments),
            unit,
            rows.iter().map(|row| row.donors.len()).sum::<usize>(),
        );
        let mut metrics = HashMap::with_capacity(rows.len());
        for row in rows {
            if self.cfg.verbose || self.cfg.very_verbose {
                let donors = row
                    .donors
                    .iter()
                    .map(|&taxid| name(taxid))
                    .collect::<Vec<_>>();
                log::info!(
                    target: "PROFILER",
                    "  {} ({}): {:.0} {} of {}",
                    name(row.genus),
                    row.genus,
                    row.fragments,
                    unit,
                    donors.join(", "),
                );
            }
            let totals = unambiguous.entry(row.genus).or_default();
            totals
                .support
                .add_weighted(row.fragments, row.equivalent, 1.0);
            totals.depth += row.depth;
            totals.mapped_len += row.mapped_len;
            classifications.record_group(row.genus, row.donors, row.fragments as usize, row.reads);
            metrics.insert(row.genus, row.witness);
        }
        metrics
    }

    /// Format the optional raw community profile from its independent EM
    /// estimate.
    ///
    /// Unlike the production table, these rows are not filtered by minimum
    /// read support, breadth/chunk amount, OEBR/COEBR, p-values, or the
    /// coordinate-evidence accept/reject decision. Mapping-quality gates,
    /// taxonomy normalization, ambiguity/LCA handling, and EM refinement have
    /// already been applied. This makes the table useful for reproducible
    /// downstream threshold exploration without changing production calls.
    fn build_raw_profile_rows(
        &self,
        estimate: &RawProfileEstimate,
        evidence_decisions: &HashMap<u32, WitnessDecision>,
        terminal_ranks: &HashMap<u32, ProfileRank>,
        host_taxid: Option<u32>,
        pathogen_map: Option<&HashMap<u32, PathogenEntry>>,
    ) -> Vec<ProfileRow> {
        let norm_factor = 0.01
            * estimate
                .totals
                .values()
                .map(|totals| totals.depth)
                .sum::<f64>();
        let num_taxa = estimate.totals.len().max(1);
        let mut entries = estimate
            .totals
            .iter()
            .map(|(&taxid, &totals)| (taxid, totals))
            .collect::<Vec<_>>();
        entries.sort_by(|left, right| {
            right
                .1
                .depth
                .partial_cmp(&left.1.depth)
                .unwrap_or(Ordering::Equal)
                .then_with(|| left.0.cmp(&right.0))
        });

        let mut rows = Vec::with_capacity(entries.len());
        for (taxid, totals) in entries {
            let mapcount = totals.support.observed;
            let equivalent_mapcount = totals.support.equivalent;
            let depth = totals.depth;
            if mapcount <= f64::EPSILON
                && equivalent_mapcount <= f64::EPSILON
                && depth <= f64::EPSILON
            {
                continue;
            }
            let maplen = totals.mapped_len;
            let abundance = if norm_factor <= f64::EPSILON {
                0.0
            } else {
                depth / norm_factor
            };
            let rank = terminal_ranks
                .get(&taxid)
                .copied()
                .map(ProfileRank::as_str)
                .or_else(|| self.taxonomy.get_rank(taxid))
                .unwrap_or("no rank")
                .to_string();
            let Some(taxname) = self.taxonomy.get_name(taxid) else {
                continue;
            };

            let evidence_metrics = evidence_decisions
                .get(&taxid)
                .map(|decision| decision.metrics);
            let (
                breadth,
                fixed_chunk,
                flex_chunk,
                expected_breadth,
                cov_prob,
                expected_flex_chunk,
                chunkcov_prob,
            ) = if let Some(metrics) = evidence_metrics {
                (
                    metrics.breadth,
                    metrics.fixed_chunk,
                    metrics.flex_chunk,
                    metrics.expected_breadth,
                    metrics.breadth_pvalue,
                    metrics.expected_flex_chunk,
                    metrics.chunk_pvalue,
                )
            } else {
                let coverage = estimate.coverage.get(&taxid).cloned().unwrap_or_default();
                let breadth = coverage.breadth.min(1.0);
                let fixed_chunk = coverage.fixed_chunk.min(1.0);
                let flex_chunk = coverage.flex_chunk.min(1.0);
                let expected_breadth = crate::profiler::evidence::expected_breadth(depth, mapcount);
                let cov_prob = if rank == "species" {
                    coverage_deviation_pvalue(breadth, expected_breadth, Some(depth), mapcount, num_taxa)
                } else {
                    None
                };
                let genome_len = if depth <= f64::EPSILON {
                    0.0
                } else {
                    maplen / depth
                };
                let num_flex_chunks = if genome_len <= 0.0 {
                    1
                } else {
                    (genome_len.sqrt() / 2.0).floor().max(1.0) as usize
                };
                let expected_flex_chunk = 1.0 - (1.0 - 1.0 / num_flex_chunks as f64).powf(mapcount);
                let chunkcov_prob = if rank == "species" {
                    coverage_deviation_pvalue(flex_chunk, expected_flex_chunk, None, mapcount, num_taxa)
                } else {
                    None
                };
                (
                    breadth,
                    fixed_chunk,
                    flex_chunk,
                    expected_breadth,
                    cov_prob,
                    expected_flex_chunk,
                    chunkcov_prob,
                )
            };

            let Some(annotations) = self.row_annotations(taxid, host_taxid, pathogen_map) else {
                continue;
            };
            rows.push(ProfileRow {
                taxname: taxname.to_string(),
                taxid,
                rank,
                equivalent_mapcount,
                depth,
                abundance,
                breadth,
                expected_breadth,
                breadth_pvalue: cov_prob,
                fixed_chunk,
                flex_chunk,
                expected_flex_chunk,
                chunk_pvalue: chunkcov_prob,
                best_fit: None,
                annotations,
            });
        }
        rows
    }

    /// `--identity-margin`: report at their genus the reported species whose
    /// best placements are clearly less identical than those of the sample's
    /// typical species. Returns the species moved, each with its genus.
    fn report_low_identity_species_at_genus(
        &self,
        rows: &mut Vec<ProfileRow>,
        placements: &HashMap<u32, BestPlacements>,
        margin: f64,
        host_taxid: Option<u32>,
        pathogen_map: Option<&HashMap<u32, PathogenEntry>>,
    ) -> HashMap<u32, u32> {
        let species = rows
            .iter()
            .filter(|row| row.rank == "species")
            .map(|row| row.taxid);
        let Some(judged) = identity_baseline(species, placements) else {
            log::info!(
                target: "PROFILER",
                "Identity margin not applied: fewer than 3 reported species have best placements."
            );
            return HashMap::new();
        };
        let moves: HashMap<u32, u32> = judged
            .below(margin)
            .filter_map(|species| {
                let genus = resolve_genus_taxid(&self.taxonomy, species.taxid)?;
                Some((species.taxid, genus))
            })
            .collect();
        let moved = move_rows_to_genus(rows, &moves, |genus| {
            let name = self.taxonomy.get_name(genus)?.to_string();
            Some((name, self.row_annotations(genus, host_taxid, pathogen_map)?))
        });
        log::info!(
            target: "PROFILER",
            "Identity baseline {:.5} over {} species; {} species more than {} below it reported at genus.",
            judged.baseline,
            judged.species.len(),
            moved.len(),
            margin,
        );
        if self.cfg.verbose || self.cfg.very_verbose {
            // --verbose lists the species below the baseline, --very-verbose
            // every judged species.
            let name = |taxid: u32| self.taxonomy.get_name(taxid).unwrap_or("unknown");
            for species in &judged.species {
                let below = judged.is_below(species, margin);
                if !below && !self.cfg.very_verbose {
                    continue;
                }
                let outcome = match moved.get(&species.taxid) {
                    Some(&genus) => format!("reported at {} ({genus})", name(genus)),
                    None if below => "kept: no genus to report".to_string(),
                    None => "kept".to_string(),
                };
                log::info!(
                    target: "PROFILER",
                    "  {} ({}): identity {:.5} over {} fragments, {}",
                    name(species.taxid),
                    species.taxid,
                    species.identity,
                    species.fragments,
                    outcome,
                );
            }
        }
        moved
    }

    /// Host and disease columns of a profile row; rows have none without
    /// `--host`. `None` when the host filter excludes the taxon.
    fn row_annotations(
        &self,
        taxid: u32,
        host_taxid: Option<u32>,
        pathogen_map: Option<&HashMap<u32, PathogenEntry>>,
    ) -> Option<Option<ProfileAnnotations>> {
        let Some(host_taxid) = host_taxid else {
            return Some(None);
        };
        let annotations = match pathogen_map.and_then(|map| map.get(&taxid)) {
            Some(entry) => {
                let (report_it, host_names, host_taxids, diseases) =
                    self.annotate_host_disease(entry, host_taxid);
                if !report_it {
                    return None;
                }
                let or_na = |value: Option<String>| value.unwrap_or_else(|| "NA".to_string());
                ProfileAnnotations {
                    host_names: or_na(host_names),
                    host_taxids: or_na(host_taxids),
                    diseases: or_na(diseases),
                }
            }
            None => ProfileAnnotations {
                host_names: "unknown".to_string(),
                host_taxids: "unknown".to_string(),
                diseases: "unknown".to_string(),
            },
        };
        Some(Some(annotations))
    }

    /// Annotate a pathogen with host and disease information.
    ///
    /// Checks if the pathogen can infect the specified host by examining
    /// the taxonomy lineage of the pathogen's known hosts.
    ///
    /// # Arguments
    /// * `entry` - Pathogen metadata from the pathogen-host mapping file
    /// * `host_taxid` - NCBI taxonomy ID of the host organism
    ///
    /// # Returns
    ///
    /// Tuple of (should_report, host_names, host_taxids, diseases)
    /// - should_report: true if this pathogen can infect the specified host
    /// - host_names: Semicolon-separated host names (or None)
    /// - host_taxids: Semicolon-separated host taxids (or None)
    /// - diseases: Semicolon-separated disease names (or None)
    fn annotate_host_disease(
        &self,
        entry: &PathogenEntry,
        host_taxid: u32,
    ) -> (bool, Option<String>, Option<String>, Option<String>) {
        let mut host_taxids_set: HashSet<u32> = HashSet::new();
        let mut host_names_set: HashSet<String> = HashSet::new();
        let mut diseases_set: HashSet<String> = HashSet::new();
        if !entry.host_taxids.is_empty() && entry.host_taxids != "unknown" {
            for tax in entry.host_taxids.split(';') {
                if let Ok(val) = tax.parse::<u32>() {
                    host_taxids_set.insert(val);
                }
            }
        }
        if !entry.host_names.is_empty() && entry.host_names != "unknown" {
            for name in entry.host_names.split(';') {
                host_names_set.insert(name.to_string());
            }
        }
        if !entry.diseases.is_empty() && entry.diseases != "unknown" {
            for disease in entry.diseases.split(';') {
                diseases_set.insert(disease.to_string());
            }
        }
        let host_taxids = if host_taxids_set.is_empty() {
            None
        } else {
            Some(
                host_taxids_set
                    .iter()
                    .map(|t| t.to_string())
                    .collect::<Vec<_>>()
                    .join(";"),
            )
        };
        let host_names = if host_names_set.is_empty() {
            None
        } else {
            Some(host_names_set.iter().cloned().collect::<Vec<_>>().join(";"))
        };
        let diseases = if diseases_set.is_empty() {
            None
        } else {
            Some(diseases_set.iter().cloned().collect::<Vec<_>>().join(";"))
        };
        let mut report_it = false;
        if host_taxids_set.contains(&host_taxid) {
            report_it = true;
        } else if !host_taxids_set.is_empty() {
            let mut lineage: HashSet<u32> = HashSet::new();
            for taxid_val in &host_taxids_set {
                lineage.extend(
                    self.taxonomy
                        .get_parents(*taxid_val)
                        .into_iter()
                        .map(|(t, _, _)| t),
                );
            }
            if lineage.contains(&host_taxid) {
                report_it = true;
            }
        }
        (report_it, host_names, host_taxids, diseases)
    }
}

/// Process a batch of reads for taxonomy assignment.
///
/// For each read in the batch:
/// 1. Filter alignments by mapped_len, identity, and fraction thresholds
/// 2. Parse reference names to extract taxonomy information
/// 3. Keep the best (highest identity) alignment per species
/// 4. Record coverage intervals for breadth calculation
///
/// # Arguments
/// * `chunk` - Slice of ReadBundle to process
/// * `cfg` - Profiler configuration with filtering thresholds
/// * `asm_lengths` - Assembly lengths for depth normalization
///
/// # Returns
///
/// AssignmentResult with read assignments and coverage intervals for this batch
/// Per-taxid best-alignment bookkeeping used while scanning a read's
/// alignments. `best` holds `(identity, depth, mapped_len, sampling_scale)`;
/// `order` keeps first-seen order for deterministic output; `intervals` captures the
/// coverage interval triple from the currently-best alignment.
type BestMap = HashMap<u32, (f64, f64, f64, f64)>;
type IntervalMap = HashMap<u32, (AccKey, (usize, usize), (usize, usize), (usize, usize))>;

struct CoordinateEvidenceResult {
    decisions: HashMap<u32, WitnessDecision>,
    selected_taxa: HashSet<u32>,
    /// The genome each present taxon is reported on.
    best_fit: HashMap<u32, BestFitGenome>,
}

/// A present taxon's best-fit genome, which also sets its abundance.
#[derive(Clone, Copy, Debug)]
struct BestFitGenome {
    assembly_id: u32,
    asm_len: f64,
    sampling_scale: f64,
}

/// A taxon's abundance comes from the genome it is reported on: each
/// fragment's depth for a present taxon is its mapped bases over that
/// genome's length, with that genome's sampling scale, whichever of the
/// taxon's genomes the fragment placed best on. Depths on marker genomes keep
/// their gene model.
fn rescale_to_best_fit(
    assignments: &mut [ReadAssignment],
    best_fit: &HashMap<u32, BestFitGenome>,
    markers: Option<&MarkerModel>,
) {
    for assignment in assignments
        .iter_mut()
        .filter(|assignment| assignment.requires_terminal_evidence)
    {
        for (index, taxid) in assignment.taxids.iter().enumerate() {
            let Some(genome) = best_fit.get(taxid) else {
                continue;
            };
            if markers.is_some_and(|model| model.is_marker_genome(genome.assembly_id)) {
                continue;
            }
            assignment.depths[index] = assignment.mapped_lens[index] / genome.asm_len.max(1.0);
            assignment.sampling_scales[index] = genome.sampling_scale;
        }
    }
}

/// Choose the one descendant witness whose coverage a majority-LCA row should
/// report.
///
/// A majority-LCA row has no coordinates of its own. Maximizing breadth, fixed
/// chunk and flex chunk independently across descendants could take each from a
/// different genome, and the report path then divides those by an expectation
/// rebuilt from the group's pooled depth and count -- comparing one genome's
/// observation against another population's expectation. Reporting a single
/// witness keeps observed, expected, OEBR, COEBR and both p-values on one
/// coordinate system.
///
/// Selection is the codebase's existing witness preference, made total with a
/// lowest-taxid tiebreak so the choice cannot depend on candidate order.
/// Returns `None` when no member has coordinate evidence at all, as for a
/// group of genus-fallback reads.
fn coherent_lca_witness(
    taxids: &[u32],
    evidence_decisions: &HashMap<u32, WitnessDecision>,
) -> Option<(u32, WitnessDecision)> {
    taxids
        .iter()
        .filter_map(|taxid| {
            evidence_decisions
                .get(taxid)
                .map(|decision| (*taxid, *decision))
        })
        .reduce(coherent_witness_max)
}

/// Return the preferred coherent witness using a total, order-independent
/// comparison. The taxid tiebreak is required because `witness_is_better`
/// intentionally compares evidence strength only.
#[inline]
fn coherent_witness_max(
    left: (u32, WitnessDecision),
    right: (u32, WitnessDecision),
) -> (u32, WitnessDecision) {
    let right_is_better = witness_is_better(right.1, left.1);
    let equivalent = !right_is_better && !witness_is_better(left.1, right.1);
    if right_is_better || (equivalent && right.0 < left.0) {
        right
    } else {
        left
    }
}

/// Retain one coherent coordinate witness across every ambiguous group that
/// resolves to the same LCA. This is O(1) per folded group and cannot depend on
/// the randomized iteration order of `ambiguous_counts`.
fn retain_coherent_lca_witness(
    witnesses: &mut HashMap<u32, (u32, WitnessDecision)>,
    lca_taxid: u32,
    candidate: (u32, WitnessDecision),
) {
    match witnesses.entry(lca_taxid) {
        Entry::Vacant(entry) => {
            entry.insert(candidate);
        }
        Entry::Occupied(mut entry) => {
            let retained = coherent_witness_max(*entry.get(), candidate);
            entry.insert(retained);
        }
    }
}

/// Merge the pooled LCA coverage of genus-fallback groups without mutating
/// direct taxon coverage. This preserves the established component-wise maxima while making
/// the result independent of ambiguous-group iteration order.
fn merge_legacy_lca_coverage(
    coverage: &mut HashMap<u32, TaxonCoverage>,
    lca_taxid: u32,
    candidate: TaxonCoverage,
) {
    let retained = coverage.entry(lca_taxid).or_default();
    retained.breadth = retained.breadth.max(candidate.breadth);
    retained.fixed_chunk = retained.fixed_chunk.max(candidate.fixed_chunk);
    retained.flex_chunk = retained.flex_chunk.max(candidate.flex_chunk);
}

/// Select output coverage for the legacy path. Evidence directly assigned to
/// a taxon takes precedence; otherwise an LCA row uses the deterministic
/// descendant fallback, then ordinary taxon coverage.
fn legacy_output_coverage(
    taxid: u32,
    direct_taxa_before_lca: &HashSet<u32>,
    direct: &HashMap<u32, TaxonCoverage>,
    lca: &HashMap<u32, TaxonCoverage>,
) -> TaxonCoverage {
    if direct_taxa_before_lca.contains(&taxid) {
        return direct.get(&taxid).cloned().unwrap_or_default();
    }
    lca.get(&taxid)
        .or_else(|| direct.get(&taxid))
        .cloned()
        .unwrap_or_default()
}

fn strongest_witnesses(
    evaluations: &[(PresenceKey, WitnessDecision)],
) -> HashMap<u32, WitnessDecision> {
    let mut strongest: HashMap<u32, (PresenceKey, WitnessDecision)> = HashMap::new();
    for &(key, decision) in evaluations {
        match strongest.entry(key.taxid) {
            Entry::Vacant(entry) => {
                entry.insert((key, decision));
            }
            Entry::Occupied(mut entry) => {
                if keyed_witness_is_better((key, decision), *entry.get()) {
                    entry.insert((key, decision));
                }
            }
        }
    }
    strongest
        .into_iter()
        .map(|(taxid, (_, decision))| (taxid, decision))
        .collect()
}

/// Compare coherent assembly witnesses with a stable final key. The policy
/// comparison remains authoritative; assembly identity is consulted only when
/// neither witness is stronger. This is constant-time and avoids sorting all
/// witnesses merely to make equal-evidence output reproducible.
#[inline]
fn keyed_witness_is_better(
    left: (PresenceKey, WitnessDecision),
    right: (PresenceKey, WitnessDecision),
) -> bool {
    if witness_is_better(left.1, right.1) {
        return true;
    }
    if witness_is_better(right.1, left.1) {
        return false;
    }
    left.0 < right.0
}

/// What choosing and using a best-fit genome needs of each witness genome:
/// fits and fragments from its resolved fragments, length and sampling scale.
#[derive(Clone, Copy, Debug, Default)]
struct WitnessGenome {
    fits: u64,
    count: u64,
    asm_len: f64,
    sampling_scale: f64,
}

type GenomeSupport = HashMap<PresenceKey, WitnessGenome>;

/// Report each present taxon on its best-fit genome: of its genomes among
/// `witnesses` that pass, the one on which most of its resolved fragments
/// have their best placement; equal: more fragments placed on it, the larger
/// fixed chunk breadth, then the lower assembly id. Which taxa are present is
/// unchanged.
fn report_best_fit(
    decisions: &mut HashMap<u32, WitnessDecision>,
    witnesses: impl IntoIterator<Item = (PresenceKey, WitnessDecision)>,
    support: &GenomeSupport,
) -> HashMap<u32, BestFitGenome> {
    let rank = |(key, decision): &(PresenceKey, WitnessDecision)| {
        let genome = support.get(key).copied().unwrap_or_default();
        (
            genome.fits,
            genome.count,
            decision.metrics.fixed_chunk,
            Reverse(key.assembly_id),
        )
    };
    let mut best: HashMap<u32, (PresenceKey, WitnessDecision)> = HashMap::new();
    for witness in witnesses
        .into_iter()
        .filter(|(_, decision)| decision.accepted())
    {
        match best.entry(witness.0.taxid) {
            Entry::Vacant(entry) => {
                entry.insert(witness);
            }
            Entry::Occupied(mut entry) => {
                let (left, right) = (rank(&witness), rank(entry.get()));
                let better = (left.0, left.1)
                    .cmp(&(right.0, right.1))
                    .then_with(|| left.2.total_cmp(&right.2))
                    .then_with(|| left.3.cmp(&right.3))
                    .is_gt();
                if better {
                    entry.insert(witness);
                }
            }
        }
    }
    best.into_iter()
        .map(|(taxid, (key, decision))| {
            decisions.insert(taxid, decision);
            let genome = support.get(&key).copied().unwrap_or_default();
            let best_fit = BestFitGenome {
                assembly_id: key.assembly_id,
                asm_len: genome.asm_len,
                sampling_scale: genome.sampling_scale,
            };
            (taxid, best_fit)
        })
        .collect()
}

/// A species absent before MAP may be added on a witness that passes after MAP
/// and already met the pre-MAP fixed-chunk guard on the same genome.
fn map_addition_qualifies(
    key: PresenceKey,
    post: WitnessDecision,
    baseline_selected: &HashSet<u32>,
    baseline_by_key: &HashMap<PresenceKey, WitnessDecision>,
    min_chunk_breadth: f64,
) -> bool {
    if !key.profile_rank.is_species() || baseline_selected.contains(&key.taxid) {
        return false;
    }
    let Some(premap_fraction) = post.required_premap_guard.fraction() else {
        return false;
    };
    let Some(pre) = baseline_by_key.get(&key) else {
        return false;
    };
    let premap_threshold = min_chunk_breadth * premap_fraction;
    pre.metrics.fixed_chunk.partial_cmp(&premap_threshold) == Some(Ordering::Greater)
}

/// Select MAP-only species additions without allowing MAP evidence to satisfy
/// its own independent-coverage guard.
///
/// The lookup is keyed by both species and assembly. This is intentional: a
/// broad pre-MAP witness on one assembly cannot authorize MAP-created
/// coverage on another assembly of the same species.
fn guarded_map_additions(
    baseline_selected: &HashSet<u32>,
    baseline_by_key: &HashMap<PresenceKey, WitnessDecision>,
    post_evaluations: &[(PresenceKey, WitnessDecision)],
    min_chunk_breadth: f64,
) -> HashMap<u32, WitnessDecision> {
    let mut additions: HashMap<u32, (PresenceKey, WitnessDecision)> = HashMap::new();
    for &(key, post) in post_evaluations {
        if !map_addition_qualifies(
            key,
            post,
            baseline_selected,
            baseline_by_key,
            min_chunk_breadth,
        ) {
            continue;
        }
        match additions.entry(key.taxid) {
            Entry::Vacant(entry) => {
                entry.insert((key, post));
            }
            Entry::Occupied(mut entry) => {
                if keyed_witness_is_better((key, post), *entry.get()) {
                    entry.insert((key, post));
                }
            }
        }
    }
    additions
        .into_iter()
        .map(|(taxid, (_, decision))| (taxid, decision))
        .collect()
}

/// Build the seeded prior used by guarded short-read exact-tie MAP.
///
/// Only species with at least one unique-candidate fragment are eligible and
/// only their unique-candidate read/mate-1 depth seeds the prior. Every ordinary
/// ambiguous candidate group then updates that prior by fractional EM. This
/// is deliberately separate from the final abundance EM: MAP can contribute
/// coordinate evidence, but it never hard-assigns abundance.
fn seeded_short_read_map_priors(
    assignments: &[ReadAssignment],
    eligible_species: &HashSet<u32>,
    max_iterations: usize,
) -> (HashMap<u32, f64>, usize) {
    if eligible_species.is_empty() {
        return (HashMap::new(), 0);
    }
    let mut species = eligible_species.iter().copied().collect::<Vec<_>>();
    species.sort_unstable();
    let indices = species
        .iter()
        .enumerate()
        .map(|(index, &taxid)| (taxid, index))
        .collect::<HashMap<_, _>>();
    let mut base = vec![0.0f64; species.len()];
    let mut grouped: HashMap<Vec<usize>, Vec<f64>> = HashMap::new();

    for assignment in assignments {
        if assignment.rank != AssignRank::Species || assignment.taxids.is_empty() {
            continue;
        }
        if assignment.taxids.len() == 1 {
            if let Some(&index) = indices.get(&assignment.taxids[0]) {
                let depth = assignment.depths.first().copied().unwrap_or(0.0);
                if depth.is_finite() && depth > 0.0 {
                    base[index] += depth;
                }
            }
            continue;
        }

        let mut key = Vec::with_capacity(assignment.taxids.len());
        let mut depths = Vec::with_capacity(assignment.taxids.len());
        for (&taxid, &depth) in assignment.taxids.iter().zip(&assignment.depths) {
            if let Some(&index) = indices.get(&taxid) {
                key.push(index);
                depths.push(if depth.is_finite() && depth > 0.0 {
                    depth
                } else {
                    0.0
                });
            }
        }
        if key.is_empty() {
            continue;
        }
        if let Some(total) = grouped.get_mut(key.as_slice()) {
            for (target, value) in total.iter_mut().zip(depths) {
                *target += value;
            }
        } else {
            grouped.insert(key, depths);
        }
    }

    let mut groups = grouped.into_iter().collect::<Vec<_>>();
    groups.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    let mut current = base.clone();
    let mut iterations = 0usize;
    for iteration in 0..max_iterations.max(1) {
        let mut next = base.clone();
        for (candidate_indices, depths) in &groups {
            if candidate_indices.len() == 1 {
                next[candidate_indices[0]] += depths[0];
                continue;
            }
            let denominator = candidate_indices
                .iter()
                .map(|&index| current[index])
                .sum::<f64>();
            if denominator <= f64::EPSILON || !denominator.is_finite() {
                continue;
            }
            for (&index, &depth) in candidate_indices.iter().zip(depths) {
                next[index] += depth * (current[index] / denominator);
            }
        }
        let numerator = next
            .iter()
            .zip(&current)
            .map(|(next, current)| (next - current).abs())
            .sum::<f64>();
        let denominator = current.iter().sum::<f64>().max(1e-12);
        current = next;
        iterations = iteration + 1;
        if numerator / denominator <= 1e-5 {
            break;
        }
    }

    let priors = species
        .into_iter()
        .zip(current)
        .filter(|(_, depth)| depth.is_finite() && *depth > 0.0)
        .collect();
    (priors, iterations)
}

/// Add one guarded-MAP placement and report the `presence` index it landed on,
/// so the caller can mark that witness's coverage as needing recomputation.
fn add_exact_tie_placement(
    presence: &mut Vec<(PresenceKey, PresenceAssembly)>,
    presence_indices: &mut HashMap<PresenceKey, usize>,
    taxid: u32,
    placement: ExactTiePlacement,
) -> Result<usize> {
    let key = PresenceKey {
        taxid,
        profile_rank: ProfileRank::Species,
        assembly_id: placement.assembly_id,
    };
    let index = if let Some(&index) = presence_indices.get(&key) {
        index
    } else {
        let index = presence.len();
        presence.push((key, PresenceAssembly::default()));
        presence_indices.insert(key, index);
        index
    };
    let assembly = &mut presence[index].1;
    if assembly.asm_len <= 0.0 {
        assembly.asm_len = placement.asm_len;
    } else {
        let tolerance = 1e-12 * assembly.asm_len.abs().max(placement.asm_len.abs()).max(1.0);
        if (assembly.asm_len - placement.asm_len).abs() > tolerance {
            bail!(
                "inconsistent effective assembly size for species {} assembly {}: {} and {}",
                taxid,
                placement.assembly_id,
                assembly.asm_len,
                placement.asm_len
            );
        }
    }
    assembly.count = assembly.count.saturating_add(1);
    assembly.reference_span_bases = assembly
        .reference_span_bases
        .saturating_add(u64::from(placement.span.1.saturating_sub(placement.span.0)));
    assembly.intervals.push(PresenceInterval {
        reference_id: placement.reference_id,
        span: placement.span,
        mate_span: placement.mate_span,
        fixed_chunk: placement.fixed_chunk,
        flex_chunk: placement.flex_chunk,
    });
    Ok(index)
}

/// How the coverage of an assembly-coherent witness is measured: a marker
/// genome by the genes its reads hit, any other genome by its intervals, with
/// the window spread when the genus roll-up needs it.
#[derive(Clone, Copy)]
struct CoverageModel<'a> {
    markers: Option<&'a MarkerModel>,
    layout: Option<&'a GenomeLayout>,
}

impl CoverageModel<'_> {
    fn measure(
        self,
        key: &PresenceKey,
        assembly: &mut PresenceAssembly,
        multiple_testing_taxa: usize,
    ) -> EvidenceCoverage {
        if let Some(markers) = self
            .markers
            .and_then(|model| model.coverage(key.assembly_id, assembly, multiple_testing_taxa))
        {
            return markers;
        }
        let mut metrics = calculate_coverage(assembly, multiple_testing_taxa);
        metrics.window_spread = self
            .layout
            .and_then(|layout| layout.window_spread(key.assembly_id, &assembly.intervals));
        metrics
    }
}

fn resolve_coordinate_evidence(
    cfg: &ProfilerConfig,
    policy: EvidencePolicy,
    assignment: &mut AssignmentResult,
    min_breadth: f64,
    standard_min_chunk_breadth: f64,
    safeguard_min_chunk_breadth: f64,
    coverage_model: CoverageModel<'_>,
) -> Result<CoordinateEvidenceResult> {
    let thresholds = EvidenceThresholds {
        mode: policy.mode,
        min_breadth,
        standard_min_chunk_breadth,
        safeguard_min_chunk_breadth,
        min_oebr: cfg.min_oebr,
        min_coebr: cfg.min_coebr,
        strict: cfg.strict,
        min_marker_genes: cfg.min_marker_genes,
    };
    let mut presence = std::mem::take(&mut assignment.presence)
        .into_iter()
        .collect::<Vec<_>>();
    let baseline_witness_count = presence.len();
    // Correct for every hypothesis actually tested: one breadth test and one
    // chunk test run per assembly-coherent witness. `uniquely_supported_taxa`
    // was a proxy that matched neither the tests performed nor the taxa
    // reported, and it under-corrected -- leaving adjusted p-values smaller
    // than the testing burden justifies, so the gate rejected coherent taxa.
    // Both the baseline and post-MAP passes use this same family so the two
    // remain directly comparable.
    let multiple_testing_taxa = baseline_witness_count.max(1);
    let coverage = |key: &PresenceKey, assembly: &mut PresenceAssembly| {
        coverage_model.measure(key, assembly, multiple_testing_taxa)
    };
    if cfg.verbose || cfg.very_verbose {
        log::info!(
            target: "PROFILER",
            "Coverage p-value hypothesis family: {} assembly-coherent witnesses tested (uniquely supported taxa: {}).",
            multiple_testing_taxa,
            assignment.uniquely_supported_taxa.len(),
        );
    }
    let pool = ThreadPoolBuilder::new()
        .num_threads(cfg.threads.max(1))
        .build()?;
    let baseline_evaluations = pool.install(|| {
        presence
            .par_iter_mut()
            .map(|(key, assembly)| {
                let metrics = coverage(key, assembly);
                let evidence = assignment
                    .taxon_evidence
                    .get(&key.taxid)
                    .copied()
                    .unwrap_or_default();
                (
                    *key,
                    evaluate_ranked_witness(metrics, evidence, thresholds, key.profile_rank),
                )
            })
            .collect::<Vec<_>>()
    });
    let mut baseline_decisions = strongest_witnesses(&baseline_evaluations);
    let baseline_selected = baseline_decisions
        .iter()
        .filter_map(|(&taxid, decision)| decision.accepted().then_some(taxid))
        .collect::<HashSet<_>>();
    let support: GenomeSupport = presence
        .iter()
        .map(|(key, assembly)| {
            let genome = WitnessGenome {
                fits: assembly.fits,
                count: assembly.count,
                asm_len: assembly.asm_len,
                sampling_scale: assembly.sampling_scale,
            };
            (*key, genome)
        })
        .collect();
    let baseline_best_fit = report_best_fit(
        &mut baseline_decisions,
        baseline_evaluations.iter().copied(),
        &support,
    );

    if !policy.collect_exact_ties {
        if cfg.verbose || cfg.very_verbose {
            let standard = baseline_decisions
                .values()
                .filter(|decision| decision.standard)
                .count();
            let rescued = baseline_decisions
                .values()
                .filter(|decision| decision.rescued())
                .count();
            log::info!(
                target: "PROFILER",
                "Resolved coordinate evidence ({}) selected {} terminal taxa ({} standard, {} rescue) from {} assembly-coherent witnesses.",
                policy.mode.label(),
                baseline_selected.len(),
                standard,
                rescued,
                baseline_witness_count,
            );
        }
        return Ok(CoordinateEvidenceResult {
            decisions: baseline_decisions,
            selected_taxa: baseline_selected,
            best_fit: baseline_best_fit,
        });
    }

    // This is also checked by the CLI, but library callers receive a safe
    // baseline rather than an unguarded MAP pass if the amount guard is absent.
    if !policy.mode.supports_guarded_map() || safeguard_min_chunk_breadth <= 0.0 {
        log::warn!(
            target: "PROFILER",
            "Guarded short-read MAP requires native Illumina evidence and a positive active chunk-breadth cutoff; exact ties remain ambiguous."
        );
        return Ok(CoordinateEvidenceResult {
            decisions: baseline_decisions,
            selected_taxa: baseline_selected,
            best_fit: baseline_best_fit,
        });
    }

    let exact_ties = std::mem::take(&mut assignment.exact_ties);
    let exact_tie_count = exact_ties.fragment_count();
    if exact_tie_count == 0 {
        return Ok(CoordinateEvidenceResult {
            decisions: baseline_decisions,
            selected_taxa: baseline_selected,
            best_fit: baseline_best_fit,
        });
    }
    // Guarded MAP is a calibrated species-only policy. Higher-rank terminal
    // endpoints remain in abundance EM but must never seed or receive MAP
    // coordinate evidence.
    let map_eligible_species = presence
        .iter()
        .filter_map(|(key, _)| {
            (key.profile_rank.is_species()
                && assignment.uniquely_supported_taxa.contains(&key.taxid))
            .then_some(key.taxid)
        })
        .collect::<HashSet<_>>();
    let (priors, prior_iterations) = seeded_short_read_map_priors(
        &assignment.read_assignments,
        &map_eligible_species,
        cfg.max_em_iterations,
    );
    if priors.is_empty() {
        return Ok(CoordinateEvidenceResult {
            decisions: baseline_decisions,
            selected_taxa: baseline_selected,
            best_fit: baseline_best_fit,
        });
    }
    let baseline_by_key = baseline_evaluations
        .iter()
        .copied()
        .collect::<HashMap<_, _>>();
    let mut presence_indices = presence
        .iter()
        .enumerate()
        .map(|(index, (key, _))| (*key, index))
        .collect::<HashMap<_, _>>();
    // Update the cloned counters directly. This avoids two temporary per-taxon
    // hash maps while preserving the immutable baseline evidence.
    let mut mapped_taxon_evidence = assignment.taxon_evidence.clone();
    let mut resolved_ties = 0usize;
    let mut posterior_sum = 0.0f64;
    // Coverage is a pure function of an assembly's intervals, count and length,
    // so only witnesses that actually receive a MAP placement need their
    // intervals merged and sorted again. Track those; witnesses appended past
    // the baseline are new and always need a first computation.
    let mut coverage_dirty = vec![false; baseline_witness_count];
    for fragment in exact_ties.fragments() {
        let Some((winner, posterior)) = guarded_map_winner(fragment.tied(), &priors) else {
            continue;
        };
        resolved_ties += 1;
        posterior_sum += posterior;
        // Every candidate of this fragment was compared against the winner, so
        // the whole fragment now belongs to a resolved comparison — not only
        // the tied prefix. Crediting just the prefix left a taxon that lost the
        // tie outright with a too-small denominator, inflating the
        // resolved-best fraction that platform and low-abundance rescue read.
        for candidate in fragment.all() {
            let evidence = mapped_taxon_evidence.entry(candidate.taxid).or_default();
            evidence.resolved_candidate_fragments =
                evidence.resolved_candidate_fragments.saturating_add(1);
        }
        let winner_evidence = mapped_taxon_evidence.entry(winner.taxid).or_default();
        winner_evidence.best_fragments = winner_evidence.best_fragments.saturating_add(1);
        for &placement in exact_ties.placements(winner) {
            let index = add_exact_tie_placement(
                &mut presence,
                &mut presence_indices,
                winner.taxid,
                placement,
            )?;
            if let Some(dirty) = coverage_dirty.get_mut(index) {
                *dirty = true;
            }
        }
    }
    let post_evaluations = pool.install(|| {
        presence
            .par_iter_mut()
            .enumerate()
            .map(|(index, (key, assembly))| {
                // Reuse the baseline coverage unless this witness received new
                // MAP coordinates. `calculate_coverage` depends only on the
                // assembly contents and the (unchanged) hypothesis count, so
                // the reused value is exactly what recomputing would produce.
                let metrics = match coverage_dirty.get(index) {
                    Some(false) => baseline_evaluations[index].1.metrics,
                    _ => coverage(key, assembly),
                };
                // The decision itself must always be re-derived: guarded MAP
                // updates evidence denominators for every candidate of a
                // resolved fragment, so a witness whose coordinates were never
                // touched can still change tier.
                let evidence = mapped_taxon_evidence
                    .get(&key.taxid)
                    .copied()
                    .unwrap_or_default();
                (
                    *key,
                    evaluate_ranked_witness(metrics, evidence, thresholds, key.profile_rank),
                )
            })
            .collect::<Vec<_>>()
    });
    let mut post_decisions = strongest_witnesses(&post_evaluations);
    let mut selected_taxa = baseline_selected.clone();
    let additions = guarded_map_additions(
        &baseline_selected,
        &baseline_by_key,
        &post_evaluations,
        safeguard_min_chunk_breadth,
    );
    let guarded_additions = additions.len();
    let low_abundance_additions = additions
        .values()
        .filter(|decision| decision.rescues.low_abundance_pe)
        .count();
    for (taxid, decision) in additions {
        selected_taxa.insert(taxid);
        // The globally strongest post-MAP witness may belong to an assembly
        // that did not satisfy the pre-MAP guard. Report the strongest
        // qualifying same-assembly witness instead.
        post_decisions.insert(taxid, decision);
    }
    // MAP is monotonic. If added conserved intervals make an upper-ratio or
    // p-value statistic fail, retain the baseline witness and its metrics.
    for &taxid in &baseline_selected {
        let post_still_passes = post_decisions
            .get(&taxid)
            .is_some_and(|decision| decision.accepted());
        if !post_still_passes {
            if let Some(&baseline) = baseline_decisions.get(&taxid) {
                post_decisions.insert(taxid, baseline);
            }
        }
    }
    // Each present species is reported on its best-fit genome among the
    // witnesses that decided it: those passing after MAP (for an addition,
    // also meeting the pre-MAP guard), or the baseline ones it was retained on.
    let passes_after_map = post_evaluations
        .iter()
        .filter(|(_, decision)| decision.accepted())
        .map(|(key, _)| key.taxid)
        .collect::<HashSet<_>>();
    let deciding = post_evaluations
        .iter()
        .copied()
        .filter(|&(key, post)| {
            baseline_selected.contains(&key.taxid)
                || map_addition_qualifies(
                    key,
                    post,
                    &baseline_selected,
                    &baseline_by_key,
                    safeguard_min_chunk_breadth,
                )
        })
        .chain(baseline_evaluations.iter().copied().filter(|(key, _)| {
            baseline_selected.contains(&key.taxid) && !passes_after_map.contains(&key.taxid)
        }));
    let best_fit = report_best_fit(&mut post_decisions, deciding, &support);

    if cfg.verbose || cfg.very_verbose {
        let posterior_mean = if resolved_ties == 0 {
            0.0
        } else {
            posterior_sum / resolved_ties as f64
        };
        let premap_guard = match policy.mode {
            EvidenceMode::Paired => 0.5 * safeguard_min_chunk_breadth,
            EvidenceMode::Single => safeguard_min_chunk_breadth,
            EvidenceMode::Long => unreachable!("long-read MAP is disabled"),
        };
        let post_requirement = match policy.mode {
            EvidenceMode::Paired => "standard-or-rescue",
            EvidenceMode::Single => "standard-only",
            EvidenceMode::Long => unreachable!("long-read MAP is disabled"),
        };
        log::info!(
            target: "PROFILER",
            "Guarded {} exact-tie MAP: eligible_species={} prior_iterations={} exact_ties={} resolved={} mean_posterior={:.5}; baseline_species={} guarded_additions={} (low_abundance={}) final_species={} (ordinary same-assembly pre-MAP fixed_chunk > {:.6}; PE low-abundance guard > {:.6}; post-MAP {}).",
            policy.mode.label(),
            priors.len(),
            prior_iterations,
            exact_tie_count,
            resolved_ties,
            posterior_mean,
            baseline_selected.len(),
            guarded_additions,
            low_abundance_additions,
            selected_taxa.len(),
            premap_guard,
            0.05 * safeguard_min_chunk_breadth,
            post_requirement,
        );
    }

    Ok(CoordinateEvidenceResult {
        decisions: post_decisions,
        selected_taxa,
        best_fit,
    })
}

/// Mapped fragments as a share of the input fragments, each counting as the
/// share of a genome it represents: the sampling scale of the genome of its
/// best placement, so that an index of sampled genomes is not mistaken for a
/// sample of unknown or host DNA. Fragments placed only on marker genes, whose
/// genome size is unknown, are left out. Equals the plain mapping rate for a
/// database of whole genomes.
fn genome_equivalent_mapping_rate(
    bundles: &[ReadBundle],
    total_reads: usize,
    markers: Option<&MarkerModel>,
) -> f64 {
    if total_reads == 0 {
        return 0.0;
    }
    let mapped: f64 = bundles
        .iter()
        .filter_map(|bundle| {
            bundle
                .alignments
                .iter()
                .find(|row| markers.is_none_or(|model| !model.is_marker_genome(row.assembly_id)))
                .map(|row| f64::from(row.sampling_scale))
        })
        .fold(0.0, |total, scale| total + scale);
    (mapped / total_reads as f64).min(1.0)
}

/// Read length, template length and read layout of this sample, measured
/// on the mappings, for the coverage expected on marker genes.
fn sample_fragment_geometry(
    cfg: &ProfilerConfig,
    bundles: &[ReadBundle],
    read_lengths: &LengthHistogram,
    lone_mates: bool,
) -> FragmentGeometry {
    let read_len = read_lengths.median().unwrap_or(1);
    let mut templates = LengthHistogram::default();
    if cfg.is_paired {
        for rows in bundles.iter().map(|bundle| &bundle.alignments) {
            for mate1 in rows
                .iter()
                .filter(|row| row.role == MappingRole::Mate1 && row.pair_id != 0)
            {
                if let Some(mate2) = rows
                    .iter()
                    .find(|row| row.role == MappingRole::Mate2 && row.pair_id == mate1.pair_id)
                {
                    let start = mate1.span.0.min(mate2.span.0);
                    let end = mate1.span.1.max(mate2.span.1);
                    templates.add(end.saturating_sub(start) as usize);
                }
            }
        }
    }
    let layout = match (cfg.is_paired, lone_mates) {
        (false, _) if cfg.coverage_gate_is_or() => FragmentLayout::LongReads,
        (false, _) => FragmentLayout::Reads,
        (true, false) => FragmentLayout::Pairs,
        (true, true) => FragmentLayout::PairsAndLoneMates,
    };
    let template_len = templates.median().unwrap_or(read_len);
    FragmentGeometry::new(layout, read_len, template_len, cfg.mapped_len, cfg.fraction)
}

/// Depth that one fragment adds to its candidate: its mapped bases over the
/// genome length. A fragment on a marker gene stands for one fragment of the
/// whole genome: it adds its whole read length, also when the read is clipped
/// at a gene end, over the placements that the genome's genes offer a read of
/// its length, so that marker and genome taxa have comparable depths; a long
/// read adds this for every gene of the genome it covers, as the placements
/// count every gene. Marker genes that cannot hold a read keep mapped bases
/// over their summed length.
fn fragment_depth(
    placement: CandidatePlacement,
    mapped: u32,
    markers: Option<&MarkerModel>,
) -> f64 {
    match markers
        .and_then(|model| model.placements_for_read(placement.assembly_id, placement.read_len))
    {
        Some(placements) => {
            f64::from(1 + placement.other_genes) * f64::from(placement.read_len) / placements
        }
        None => f64::from(mapped) / placement.asm_len.max(1.0),
    }
}

/// Native short-/long-read assignment policy.
///
/// All valid terminal-taxon candidates are retained through evidence
/// resolution; endpoints with independent presence support remain available
/// to abundance EM. Baseline coordinate evidence is deliberately stricter:
/// only a platform-specific, margin-resolved best endpoint contributes
/// intervals. Guarded MAP can later add an exact short-read tie only when its
/// posterior guard passes. Each represented assembly remains an independent
/// witness; the one of a fragment's best placement is credited with a fit.
fn process_evidence_batch(
    chunk: &[ReadBundle],
    cfg: &ProfilerConfig,
    policy: EvidencePolicy,
    fallback_mode: bool,
    species_to_genus: &HashMap<u32, u32>,
    markers: Option<&MarkerModel>,
) -> AssignmentResult {
    let mut result = AssignmentResult::default();
    let mut scratch = CandidateScratch::default();
    let mut indexed_candidates: Vec<usize> = Vec::new();
    let mut best_identities: Vec<(u32, f32)> = Vec::new();
    let mut best_placements = BestPlacementScratch::default();

    // Genus fallback remains a separate, lower-rank path and runs only when
    // no terminal candidate survives the native policy.
    let genus_fb_mapped_min = genus_fallback_min_mapped_len(cfg);
    let mut genus_best: BestMap = HashMap::new();
    let mut genus_order: Vec<u32> = Vec::new();
    let mut genus_intervals: IntervalMap = HashMap::new();

    for read in chunk {
        if policy.collect_identities {
            // Best mapping identity of the fragment to each taxon, before any gate.
            best_identities.clear();
            for row in &read.alignments {
                match best_identities
                    .iter_mut()
                    .find(|(taxid, _)| *taxid == row.profile_taxid)
                {
                    Some((_, best)) => *best = best.max(row.identity),
                    None => best_identities.push((row.profile_taxid, row.identity)),
                }
            }
            for &(taxid, identity) in &best_identities {
                let evidence = result.taxon_evidence.entry(taxid).or_default();
                evidence.mapped_fragments = evidence.mapped_fragments.saturating_add(1);
                evidence.identity_sum += f64::from(identity);
            }
        }
        if policy.collect_best_placements {
            best_placements.credit(&read.alignments, &mut result.best_placements);
        }

        let candidates = scratch.build(&read.alignments, policy);
        result.filter_stats.evaluated = result
            .filter_stats
            .evaluated
            .saturating_add(read.alignments.len());

        if !candidates.is_empty() {
            let resolved_index = resolved_candidate_index(candidates, policy.mode);
            let profile_boundary_incomplete = candidates
                .iter()
                .any(|candidate| candidate.candidate_boundary.profile_incomplete());
            if profile_boundary_incomplete {
                result.boundary_ambiguous_fragments =
                    result.boundary_ambiguous_fragments.saturating_add(1);
            }
            result.filter_stats.passed =
                result.filter_stats.passed.saturating_add(candidates.len());
            for candidate in candidates {
                let evidence = result.taxon_evidence.entry(candidate.taxid).or_default();
                evidence.candidate_fragments = evidence.candidate_fragments.saturating_add(1);
                if resolved_index.is_some() {
                    evidence.resolved_candidate_fragments =
                        evidence.resolved_candidate_fragments.saturating_add(1);
                }
            }
            if candidates.len() == 1 && !profile_boundary_incomplete {
                result.uniquely_supported_taxa.insert(candidates[0].taxid);
            }

            // Stable taxid order canonicalizes ambiguous EM keys independent
            // of mapping row order and candidate score ties.
            indexed_candidates.clear();
            indexed_candidates.extend(0..candidates.len());
            indexed_candidates.sort_unstable_by_key(|&index| {
                (candidates[index].taxid, candidates[index].profile_rank)
            });
            let mut taxids = Vec::with_capacity(candidates.len());
            let mut depths = Vec::with_capacity(candidates.len());
            let mut mapped_lens = Vec::with_capacity(candidates.len());
            let mut sampling_scales = Vec::with_capacity(candidates.len());
            // Candidates arrive best placement first.
            let best = indexed_candidates
                .iter()
                .position(|&index| index == 0)
                .unwrap_or(0) as u32;
            for &index in &indexed_candidates {
                let candidate = &candidates[index];
                let placement = candidate.best;
                let mapped = match policy.mode {
                    EvidenceMode::Paired | EvidenceMode::Single => placement.r1_mapped,
                    EvidenceMode::Long => placement.both_mapped,
                };
                taxids.push(candidate.taxid);
                depths.push(fragment_depth(placement, mapped, markers));
                mapped_lens.push(f64::from(mapped));
                sampling_scales.push(placement.sampling_scale);
            }
            result.read_assignments.push(ReadAssignment {
                read: Arc::clone(&read.read_name),
                rank: homogeneous_candidate_rank(candidates),
                requires_terminal_evidence: true,
                taxids,
                depths,
                mapped_lens,
                sampling_scales,
                best,
            });

            if policy.collect_exact_ties
                && candidates
                    .iter()
                    .all(|candidate| candidate.profile_rank.is_species())
            {
                result.exact_ties.record(candidates);
            }

            if let Some(index) = resolved_index {
                // Direct unique or margin-resolved winners count as both
                // resolved-best and clear-best support. An unresolved tie
                // counts for neither taxon; guarded MAP may later add only its
                // successfully resolved winner to `best_fragments`.
                let evidence = result
                    .taxon_evidence
                    .entry(candidates[index].taxid)
                    .or_default();
                evidence.best_fragments = evidence.best_fragments.saturating_add(1);
                evidence.clear_best_fragments = evidence.clear_best_fragments.saturating_add(1);
                let best_assembly = candidates[index].best.assembly_id;
                for placement in &candidates[index].placements {
                    let key = PresenceKey {
                        taxid: placement.taxid,
                        profile_rank: placement.profile_rank,
                        assembly_id: placement.assembly_id,
                    };
                    let presence = result.presence.entry(key).or_default();
                    if presence.asm_len <= 0.0 {
                        presence.asm_len = placement.asm_len;
                        presence.sampling_scale = placement.sampling_scale;
                    }
                    presence.count = presence.count.saturating_add(1);
                    presence.fits += u64::from(placement.assembly_id == best_assembly);
                    presence.other_gene_hits = presence
                        .other_gene_hits
                        .saturating_add(u64::from(placement.other_genes));
                    presence.reference_span_bases = presence.reference_span_bases.saturating_add(
                        u64::from(placement.span.1.saturating_sub(placement.span.0)),
                    );
                    presence.intervals.push(PresenceInterval {
                        reference_id: placement.reference_id,
                        span: placement.span,
                        mate_span: placement.mate_span,
                        fixed_chunk: placement.fixed_chunk,
                        flex_chunk: placement.flex_chunk,
                    });
                    if cfg.estimate_unclassified {
                        if cfg.coverage_gate_is_or() {
                            let genes = std::iter::once(placement.reference_id).chain(
                                candidates[index]
                                    .other_gene_hits
                                    .iter()
                                    .filter(|hit| hit.assembly_id == placement.assembly_id)
                                    .map(|hit| hit.interval.reference_id),
                            );
                            if let Some((hits, per_hit)) = markers.and_then(|model| {
                                model.whole_genome_reads(
                                    placement.assembly_id,
                                    placement.read_len,
                                    genes,
                                )
                            }) {
                                result
                                    .marker_genome_reads
                                    .entry(placement.assembly_id)
                                    .or_default()
                                    .add(hits, per_hit);
                            }
                        } else if markers
                            .is_some_and(|model| model.is_marker_genome(placement.assembly_id))
                        {
                            result
                                .marker_genome_reads
                                .entry(placement.assembly_id)
                                .or_default()
                                .add_fragment(1 + placement.other_genes);
                        }
                    }
                }
                let resolved = &candidates[index];
                for hit in &resolved.other_gene_hits {
                    let key = PresenceKey {
                        taxid: resolved.taxid,
                        profile_rank: resolved.profile_rank,
                        assembly_id: hit.assembly_id,
                    };
                    let presence = result.presence.entry(key).or_default();
                    presence.reference_span_bases = presence.reference_span_bases.saturating_add(
                        u64::from(hit.interval.span.1.saturating_sub(hit.interval.span.0)),
                    );
                    presence.intervals.push(hit.interval);
                }
            }
            continue;
        }

        if !fallback_mode {
            continue;
        }
        genus_best.clear();
        genus_order.clear();
        genus_intervals.clear();
        for aln in &read.alignments {
            // One row per fragment: the first mate of a pair, a lone mate whose
            // mate lies beyond its reference end, or a single-end/long read.
            if aln.role == MappingRole::Mate2 && aln.pair_id != 0 {
                continue;
            }
            if !aln.profile_rank.is_species() {
                continue;
            }
            let identity = f64::from(aln.identity);
            let fraction = f64::from(aln.fraction);
            let length_ok = aln.mapped_len as usize >= genus_fb_mapped_min;
            let fraction_ok = fraction >= GENUS_FALLBACK_MIN_FRACTION;
            if identity < cfg.genus_identity
                || !coverage_passes(length_ok, fraction_ok, cfg.coverage_gate_is_or())
            {
                continue;
            }
            let Some(genus_taxid) = species_to_genus.get(&aln.profile_taxid).copied() else {
                continue;
            };
            update_best(
                &mut genus_best,
                &mut genus_order,
                &mut genus_intervals,
                genus_taxid,
                identity,
                aln.mapped_len as f64 / aln.asm_len,
                aln,
            );
        }
        if !genus_order.is_empty() {
            commit_assignment(
                &mut result,
                &read.read_name,
                AssignRank::Genus,
                &genus_order,
                &genus_best,
                &mut genus_intervals,
            );
        }
    }
    result
}

fn homogeneous_candidate_rank(candidates: &[super::evidence::TerminalCandidate]) -> AssignRank {
    let Some(first) = candidates.first() else {
        return AssignRank::Mixed;
    };
    if candidates
        .iter()
        .all(|candidate| candidate.profile_rank == first.profile_rank)
    {
        AssignRank::from_profile_rank(first.profile_rank)
    } else {
        AssignRank::Mixed
    }
}

/// Update per-taxid best-alignment bookkeeping, keeping the alignment with
/// the highest identity as representative for coverage intervals.
///
/// `AccKey` construction is deferred into the branches that actually
/// update the representative — the two `String` clones (`asm`, `acc`)
/// are skipped entirely when the alignment doesn't beat the current
/// best for its taxid.
#[inline]
fn update_best(
    best: &mut BestMap,
    order: &mut Vec<u32>,
    intervals: &mut IntervalMap,
    taxid: u32,
    identity: f64,
    depth: f64,
    aln: &AlignmentEntry,
) {
    match best.entry(taxid) {
        Entry::Vacant(entry) => {
            order.push(taxid);
            entry.insert((
                identity,
                depth,
                aln.mapped_len as f64,
                f64::from(aln.sampling_scale),
            ));
            let acc_key = AccKey {
                taxid,
                asm: aln.asm.clone(),
                acc: aln.acc.clone(),
            };
            intervals.insert(
                taxid,
                (
                    acc_key,
                    widen_pair(aln.span),
                    widen_pair(aln.fixed_chunk),
                    widen_pair(aln.flex_chunk),
                ),
            );
        }
        Entry::Occupied(mut entry) => {
            if identity > entry.get().0 {
                *entry.get_mut() = (
                    identity,
                    depth,
                    aln.mapped_len as f64,
                    f64::from(aln.sampling_scale),
                );
                let acc_key = AccKey {
                    taxid,
                    asm: aln.asm.clone(),
                    acc: aln.acc.clone(),
                };
                intervals.insert(
                    taxid,
                    (
                        acc_key,
                        widen_pair(aln.span),
                        widen_pair(aln.fixed_chunk),
                        widen_pair(aln.flex_chunk),
                    ),
                );
            }
        }
    }
}

/// Commit the best-per-taxid results for one read into the shared
/// `AssignmentResult`: appends a `ReadAssignment` and flushes coverage
/// intervals for the chosen rank.
///
/// The caller passes borrowed views of the per-read bookkeeping
/// (`order`, `best`, `intervals`) because those are hoisted outside
/// the read loop in [`process_evidence_batch`] and reused across reads. The
/// intervals map is drained (`remove(&taxid)`) so the next read can
/// `.clear()` it cheaply. To keep the three interval maps keyed by
/// the same `AccKey`, the key is cloned into the first two maps and moved
/// into the third.
fn commit_assignment(
    result: &mut AssignmentResult,
    read_name: &str,
    rank: AssignRank,
    order: &[u32],
    best: &BestMap,
    intervals: &mut IntervalMap,
) {
    let mut taxids = Vec::with_capacity(order.len());
    let mut depths = Vec::with_capacity(order.len());
    let mut maplens = Vec::with_capacity(order.len());
    let mut sampling_scales = Vec::with_capacity(order.len());
    for &taxid in order {
        if let Some((_, depth, maplen, sampling_scale)) = best.get(&taxid) {
            taxids.push(taxid);
            depths.push(*depth);
            maplens.push(*maplen);
            sampling_scales.push(*sampling_scale);
        }
        if let Some((key, span, fixed, flex)) = intervals.remove(&taxid) {
            result
                .acc_intervals
                .entry(key.clone())
                .or_default()
                .push(span);
            result
                .acc_fixed_chunk_intervals
                .entry(key.clone())
                .or_default()
                .push(fixed);
            result
                .acc_flex_chunk_intervals
                .entry(key)
                .or_default()
                .push(flex);
        }
    }
    if taxids.is_empty() {
        return;
    }
    result.read_assignments.push(ReadAssignment {
        read: Arc::<str>::from(read_name),
        rank,
        requires_terminal_evidence: false,
        taxids,
        depths,
        mapped_lens: maplens,
        sampling_scales,
        best: 0,
    });
}

/// Merge overlapping genomic intervals.
///
/// Takes a list of (start, end) coordinate pairs and returns a list
/// where overlapping intervals have been merged into single intervals.
///
/// Used to calculate breadth of coverage - the total span of unique
/// positions covered by alignments.
///
/// # Algorithm
///
/// 1. Sort intervals by start position
/// 2. Iterate through, extending current interval or starting new one
/// 3. Intervals overlap if interval.start < current.end
///
/// # Example
///
/// Input: [(1,5), (3,7), (10,15)]
/// Output: [(1,7), (10,15)]
fn merge_intervals(intervals: &[(usize, usize)]) -> Vec<(usize, usize)> {
    if intervals.len() <= 1 {
        return intervals.to_vec();
    }
    let mut sorted = intervals.to_vec();
    sorted.sort_by_key(|i| i.0);
    let mut merged = Vec::with_capacity(sorted.len());
    merged.push(sorted[0]);
    for interval in sorted.into_iter().skip(1) {
        if let Some(last) = merged.last_mut() {
            if interval.0 < last.1 {
                last.1 = last.1.max(interval.1);
            } else {
                merged.push(interval);
            }
        }
    }
    merged
}

// ============================================================================
// Genus-level fallback helpers
// ============================================================================
//
// Genus fallback can assign a read at genus rank when no species hit qualifies.
// Automatic activation is limited to low-map, non-pathogen Illumina runs on
// full databases; `--genus-fallback` enables it explicitly.

/// Minimum mapped length (inclusive) for a genus-fallback candidate on
/// **Illumina** (short-read) data. Long-read platforms use `cfg.mapped_len`
/// instead (see [`genus_fallback_min_mapped_len`]).
const GENUS_FALLBACK_MIN_MAPPED_LEN: usize = 50;
/// Minimum aligned fraction (inclusive) for a genus-fallback candidate.
const GENUS_FALLBACK_MIN_FRACTION: f64 = 0.60;

/// Mapped-length floor for genus-fallback candidates: Illumina uses a fixed
/// 50 bp minimum; all other sequencers use the same cutoff as species
/// (`cfg.mapped_len`, from defaults or `-m`).
#[inline]
fn genus_fallback_min_mapped_len(cfg: &ProfilerConfig) -> usize {
    if cfg.sequencer.is_illumina() {
        GENUS_FALLBACK_MIN_MAPPED_LEN
    } else {
        cfg.mapped_len
    }
}

/// Return `true` when the automatic genus-fallback gate is active.
fn fallback_gate_active(cfg: &ProfilerConfig, mapping_rate: f64) -> bool {
    cfg.is_low_map_run(mapping_rate)
}

/// Return whether genus-level assignment is active for this run.
#[inline]
fn genus_fallback_mode_active(
    cfg: &ProfilerConfig,
    low_map_gate_active: bool,
    db_is_subsampled: bool,
) -> bool {
    cfg.genus_fallback || (low_map_gate_active && !db_is_subsampled)
}

/// Run-level decision for the read count fed into the automatic chunk-breadth
/// curve.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ChunkBreadthPlan {
    /// `true` when the unscaled basis is mapped reads; `false` for total reads.
    use_aligned_basis: bool,
    /// `true` when aligned basis was selected by the automatic low-map rule.
    by_aligned_auto: bool,
    /// `true` when subsampled fallback uses total reads for threshold estimation.
    force_total_for_subsampled_fallback: bool,
    /// Unscaled read count selected as the chunk-breadth basis.
    basis_reads: usize,
    /// Read count actually passed to the automatic chunk-breadth curve.
    estimate_reads: usize,
    /// Whether `basis_reads` was multiplied by `FALLBACK_CHUNK_BASIS_SCALE`.
    scaled: bool,
}

/// Choose the read-count basis used to estimate chunk breadth.
#[inline]
fn plan_chunk_breadth(
    cfg: &ProfilerConfig,
    mapping_rate: f64,
    total_reads: usize,
    aligned_reads: usize,
    db_is_subsampled: bool,
    low_map_gate_active: bool,
) -> ChunkBreadthPlan {
    let chunk_breadth_pinned = cfg.chunk_breadth.is_some();
    let force_total_for_subsampled_fallback =
        db_is_subsampled && cfg.genus_fallback && !chunk_breadth_pinned;
    let by_aligned_effective =
        cfg.resolve_by_aligned(mapping_rate) && !force_total_for_subsampled_fallback;
    let basis_reads = if by_aligned_effective {
        aligned_reads
    } else {
        total_reads
    };

    let low_map_fallback_active = low_map_gate_active && !db_is_subsampled;
    let scaled = !chunk_breadth_pinned && by_aligned_effective && low_map_fallback_active;
    let estimate_reads = fallback_scale_chunk_basis(basis_reads, scaled);

    ChunkBreadthPlan {
        use_aligned_basis: by_aligned_effective,
        by_aligned_auto: by_aligned_effective && !cfg.by_aligned && low_map_gate_active,
        force_total_for_subsampled_fallback,
        basis_reads,
        estimate_reads,
        scaled,
    }
}

/// Multiplier used for automatic low-map genus fallback.
const FALLBACK_CHUNK_BASIS_SCALE: f64 = 1.5;

/// Recall-mode scaling starts only once the automatic cutoff is based on more
/// than one million effective reads/fragments. This leaves shallow samples on
/// the default curve, where lowering an already-small amount gate
/// mostly adds noise.
const RECALL_CHUNK_SCALE_MIN_BASIS: usize = 1_000_000;
const SHORT_RECALL_CHUNK_SCALE_FLOOR: f64 = 0.70;
const LONG_RECALL_CHUNK_SCALE: f64 = 0.80;

/// Standard and safeguard cutoffs are deliberately distinct. A user scale
/// first adjusts the shared automatic safeguard; recall mode may then lower
/// only the ordinary standard-evidence gate.
#[derive(Clone, Copy, Debug, PartialEq)]
struct CommunityChunkBreadthThresholds {
    /// Calibrated read-count curve before any user or mode adjustment.
    automatic: f64,
    /// Automatic cutoff after the explicit user scale. Rescue and guarded-MAP
    /// safeguards use this value.
    standard: f64,
    safeguard: f64,
    user_scale: f64,
    recall_scale: f64,
}

impl CommunityChunkBreadthThresholds {
    #[inline]
    fn is_user_scaled(self) -> bool {
        self.user_scale != 1.0 && self.automatic > 0.0
    }

    #[inline]
    fn is_recall_scaled(self) -> bool {
        self.recall_scale < 1.0 && self.safeguard > 0.0
    }

    #[inline]
    fn is_scaled(self) -> bool {
        self.is_user_scaled() || self.is_recall_scaled()
    }
}

/// Resolve community chunk-breadth thresholds once.
fn resolve_community_chunk_breadth(
    cfg: &ProfilerConfig,
    basis_reads: usize,
) -> CommunityChunkBreadthThresholds {
    // Explicit user policy is authoritative in every mode, including
    // low-biomass mode and recall mode.
    if let Some(cutoff) = cfg.chunk_breadth {
        return CommunityChunkBreadthThresholds {
            automatic: cutoff,
            standard: cutoff,
            safeguard: cutoff,
            user_scale: 1.0,
            recall_scale: 1.0,
        };
    }

    // Low-biomass mode has no automatic chunk-breadth filtering, and
    // therefore nothing to scale.
    if cfg.lowbiomass {
        return CommunityChunkBreadthThresholds {
            automatic: 0.0,
            standard: 0.0,
            safeguard: 0.0,
            user_scale: 1.0,
            recall_scale: 1.0,
        };
    }

    let automatic = community_chunk_breadth_cutoff(basis_reads, cfg.sequencer.is_illumina());
    let safeguard = scale_automatic_chunk_breadth(automatic, cfg.chunk_breadth_scale);
    let recall_scale = automatic_recall_chunk_breadth_scale(cfg, basis_reads);
    let standard = if recall_scale == 1.0 {
        safeguard
    } else {
        safeguard * recall_scale
    };
    CommunityChunkBreadthThresholds {
        automatic,
        standard,
        safeguard,
        user_scale: cfg.chunk_breadth_scale,
        recall_scale,
    }
}

/// Scale only an automatic community cutoff in recall mode. The short-read
/// curve is continuous at one million reads and approaches 0.70 smoothly:
/// `scale = 0.70 + 0.30 / m`, where `m = basis_reads / 1e6`. Long reads use
/// a fixed 0.80 multiplier above the same depth gate.
#[inline]
fn automatic_recall_chunk_breadth_scale(cfg: &ProfilerConfig, basis_reads: usize) -> f64 {
    if cfg.mode != ProfileMode::Recall
        || cfg.strain
        || cfg.lowbiomass
        || cfg.chunk_breadth.is_some()
        || basis_reads <= RECALL_CHUNK_SCALE_MIN_BASIS
    {
        return 1.0;
    }

    if cfg.sequencer.is_illumina() {
        let millions = basis_reads as f64 / 1_000_000.0;
        (SHORT_RECALL_CHUNK_SCALE_FLOOR + (1.0 - SHORT_RECALL_CHUNK_SCALE_FLOOR) / millions)
            .clamp(SHORT_RECALL_CHUNK_SCALE_FLOOR, 1.0)
    } else {
        LONG_RECALL_CHUNK_SCALE
    }
}

fn community_chunk_breadth_source(
    cfg: &ProfilerConfig,
    plan: ChunkBreadthPlan,
    thresholds: CommunityChunkBreadthThresholds,
) -> &'static str {
    if cfg.chunk_breadth.is_some() {
        "user-supplied --min-chunk-breadth"
    } else if cfg.lowbiomass {
        "disabled by --low-biomass"
    } else {
        match (thresholds.is_scaled(), plan.scaled, plan.use_aligned_basis) {
            (true, true, _) => "recall-scaled auto-estimate from scaled mapped reads",
            (true, false, true) => "recall-scaled auto-estimate from mapped reads",
            (true, false, false) => "recall-scaled auto-estimate from total reads",
            (false, true, _) => "auto-estimated from scaled mapped reads",
            (false, false, true) => "auto-estimated from mapped reads",
            (false, false, false) => "auto-estimated from total reads",
        }
    }
}

fn log_chunk_breadth_scaling(thresholds: CommunityChunkBreadthThresholds) {
    if thresholds.is_user_scaled() {
        let uncapped = thresholds.automatic * thresholds.user_scale;
        if uncapped > 1.0 {
            log::info!(
                target: "PROFILER",
                "User chunk-breadth scale adjusted the automatic safeguard: {:.6} x {} = {:.6}, capped at 1.000000. Rescue and guarded-MAP safeguards use the adjusted value.",
                thresholds.automatic,
                thresholds.user_scale,
                thresholds.safeguard,
            );
        } else {
            log::info!(
                target: "PROFILER",
                "User chunk-breadth scale adjusted the automatic safeguard: {:.6} x {} = {:.6}. Rescue and guarded-MAP safeguards use the adjusted value.",
                thresholds.automatic,
                thresholds.user_scale,
                thresholds.safeguard,
            );
        }
    }
    if thresholds.is_recall_scaled() {
        log::info!(
            target: "PROFILER",
            "Recall mode scaled only the automatic standard chunk-breadth threshold: {:.6} x {:.4} = {:.6}; rescue and guarded-MAP safeguards retain {:.6}, and OEBR/COEBR thresholds are unchanged.",
            thresholds.safeguard,
            thresholds.recall_scale,
            thresholds.standard,
            thresholds.safeguard,
        );
    }
}

/// Scale the read count used as input to the automatic chunk-breadth curve.
#[inline]
fn fallback_scale_chunk_basis(chunk_basis_reads: usize, scale: bool) -> usize {
    if scale {
        ((chunk_basis_reads as f64) * FALLBACK_CHUNK_BASIS_SCALE).round() as usize
    } else {
        chunk_basis_reads
    }
}

/// Derive the genus taxid for a species taxid using the NCBI taxonomy
/// lineage (no header parsing).
///
/// Returns `None` if the species has no ancestor at rank "genus".
fn resolve_genus_taxid(taxonomy: &Taxonomy, species_taxid: u32) -> Option<u32> {
    taxonomy
        .get_parents(species_taxid)
        .into_iter()
        .find(|(_, _, rank)| rank == "genus")
        .map(|(taxid, _, _)| taxid)
}

/// Resolve a reference taxid to its species ancestor when one is available.
fn resolve_species_taxid(taxonomy: &Taxonomy, taxid: u32) -> Option<u32> {
    if matches!(taxonomy.get_rank(taxid), Some("species")) {
        return Some(taxid);
    }
    taxonomy
        .get_parents(taxid)
        .into_iter()
        .find(|(_, _, rank)| rank == "species")
        .map(|(id, _, _)| id)
}

/// Depth-scaled minimum chunk-breadth cutoff calibrated for community
/// profiling.
///
/// Both platforms use a bounded, saturating Michaelis-Menten form in
/// `m = num_reads / 1e6`:
///
/// - Illumina: `0.66*m/(m+1) + 0.24*m/(m+100)`.
/// - Nanopore/PacBio: `0.85*m/(m+0.78)`.
///
/// The cutoff rises smoothly with sequencing depth and then plateaus. The
/// 0.95 clamp is defensive; both fitted plateaus are lower.
fn community_chunk_breadth_cutoff(num_reads: usize, illumina: bool) -> f64 {
    let m = num_reads as f64 / 1_000_000.0;
    let value = if illumina {
        0.66 * m / (m + 1.0) + 0.24 * m / (m + 100.0)
    } else {
        0.85 * m / (m + 0.78)
    };
    value.min(0.95)
}

/// Historical default lower ratio bound for non-LCA output rows.
fn rank_filter_lower_bound(rank: &str) -> f64 {
    if rank == "genus" {
        0.65
    } else {
        0.75
    }
}

/// Resolve the standard evidence-tier ratio defaults for a majority-LCA rank.
/// Unsupported taxonomy ranks retain the established filter rather than
/// guessing a new policy.
#[inline]
fn standard_lca_ratio_lower_bounds(rank: &str, mode: EvidenceMode) -> Option<(f64, f64)> {
    Some(standard_ratio_lower_bounds(
        ProfileRank::from_taxonomy_rank(rank)?,
        mode,
    ))
}

/// Established ratio filter for non-LCA output rows.
fn ratio_filter_passes(
    rank: &str,
    breadth: f64,
    expected_breadth: f64,
    flex_chunk: f64,
    expected_flex: f64,
    min_oebr: Option<f64>,
    min_coebr: Option<f64>,
) -> bool {
    let default_lower = rank_filter_lower_bound(rank);
    ratio_filter_passes_with_defaults(
        breadth,
        expected_breadth,
        flex_chunk,
        expected_flex,
        (default_lower, default_lower),
        min_oebr,
        min_coebr,
    )
}

/// Apply OEBR/COEBR filtering with independently selected default lower
/// bounds. Explicit CLI overrides remain authoritative. The inclusive 1.50
/// upper bound is shared by the established output path and standard evidence.
fn ratio_filter_passes_with_defaults(
    breadth: f64,
    expected_breadth: f64,
    flex_chunk: f64,
    expected_flex: f64,
    default_lowers: (f64, f64),
    min_oebr: Option<f64>,
    min_coebr: Option<f64>,
) -> bool {
    if expected_breadth <= 0.0 || expected_flex <= 0.0 {
        return false;
    }
    let min_oebr = min_oebr.unwrap_or(default_lowers.0);
    let min_coebr = min_coebr.unwrap_or(default_lowers.1);
    let br_ratio = breadth / expected_breadth;
    let fl_ratio = flex_chunk / expected_flex;
    br_ratio.is_finite()
        && fl_ratio.is_finite()
        && min_oebr.is_finite()
        && min_coebr.is_finite()
        && br_ratio >= min_oebr
        && br_ratio <= STANDARD_RATIO_UPPER
        && fl_ratio >= min_coebr
        && fl_ratio <= STANDARD_RATIO_UPPER
}

/// Estimate abundance from every mapping-quality-qualified assignment without
/// applying taxon-level evidence/report filters.
///
/// This mirrors the established majority-LCA and EM mechanics used by the
/// production branch, but intentionally has no classification side effects and
/// never changes the production assignment state. It is invoked only for
/// `--save-raw-profile`.
fn estimate_unfiltered_taxon_totals(
    assignments: &[ReadAssignment],
    taxonomy: &Taxonomy,
    coverage: &mut HashMap<u32, TaxonCoverage>,
    max_iterations: usize,
    verbose_logging: bool,
) -> Result<HashMap<u32, TaxonTotals>> {
    let mut unambiguous: HashMap<u32, TaxonTotals> = HashMap::new();
    let mut ambiguous: HashMap<Vec<u32>, AmbiguousTotals> = HashMap::new();

    for assignment in assignments {
        let candidate_count = assignment.taxids.len();
        if assignment.depths.len() != candidate_count
            || assignment.mapped_lens.len() != candidate_count
            || assignment.sampling_scales.len() != candidate_count
        {
            bail!(
                "inconsistent raw-profile candidate vectors for read '{}': taxids={}, depths={}, mapped_lens={}, sampling_scales={}",
                assignment.read,
                candidate_count,
                assignment.depths.len(),
                assignment.mapped_lens.len(),
                assignment.sampling_scales.len(),
            );
        }
        if candidate_count == 0 {
            continue;
        }
        for index in 0..candidate_count {
            let depth = assignment.depths[index];
            let mapped_len = assignment.mapped_lens[index];
            let sampling_scale = assignment.sampling_scales[index];
            if !depth.is_finite()
                || depth < 0.0
                || !mapped_len.is_finite()
                || mapped_len < 0.0
                || !sampling_scale.is_finite()
                || sampling_scale <= 0.0
            {
                bail!(
                    "invalid raw-profile candidate values for read '{}' taxid {}: depth={}, mapped_len={}, sampling_scale={}",
                    assignment.read,
                    assignment.taxids[index],
                    depth,
                    mapped_len,
                    sampling_scale,
                );
            }
        }

        if candidate_count == 1 {
            let entry = unambiguous.entry(assignment.taxids[0]).or_default();
            entry
                .support
                .add_weighted(1.0, assignment.sampling_scales[0], 1.0);
            entry.depth += assignment.depths[0];
            entry.mapped_len += assignment.mapped_lens[0];
            continue;
        }

        if let Some(entry) = ambiguous.get_mut(assignment.taxids.as_slice()) {
            entry.observed_count += 1.0;
            accumulate(&mut entry.depths, &assignment.depths);
            accumulate(&mut entry.mapped_lens, &assignment.mapped_lens);
            accumulate(
                &mut entry.equivalent_support_sums,
                &assignment.sampling_scales,
            );
        } else {
            ambiguous.insert(
                assignment.taxids.clone(),
                AmbiguousTotals {
                    observed_count: 1.0,
                    depths: assignment.depths.clone(),
                    mapped_lens: assignment.mapped_lens.clone(),
                    equivalent_support_sums: assignment.sampling_scales.clone(),
                },
            );
        }
    }

    // Remove candidates with no independent support exactly as the production
    // EM does. If an entire candidate set lacks such support, fold it to the
    // established 70%-majority LCA instead of inventing a species prior.
    let mut em_groups: HashMap<Vec<u32>, AmbiguousTotals> = HashMap::new();
    let mut ambiguous = ambiguous.into_iter().collect::<Vec<_>>();
    ambiguous.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    for (taxids, totals) in ambiguous {
        let independently_supported = taxids
            .iter()
            .map(|taxid| unambiguous.contains_key(taxid))
            .collect::<Vec<_>>();
        if independently_supported.iter().all(|&present| present) {
            merge_ambiguous_totals(&mut em_groups, taxids, totals);
            continue;
        }
        if independently_supported.iter().any(|&present| present) {
            let mut filtered_taxids = Vec::with_capacity(taxids.len());
            let mut filtered_depths = Vec::with_capacity(taxids.len());
            let mut filtered_maplens = Vec::with_capacity(taxids.len());
            let mut filtered_scale_sums = Vec::with_capacity(taxids.len());
            for (index, taxid) in taxids.iter().enumerate() {
                if independently_supported[index] {
                    filtered_taxids.push(*taxid);
                    filtered_depths.push(totals.depths[index]);
                    filtered_maplens.push(totals.mapped_lens[index]);
                    filtered_scale_sums.push(totals.equivalent_support_sums[index]);
                }
            }
            if filtered_taxids.len() == 1 {
                let entry = unambiguous.entry(filtered_taxids[0]).or_default();
                entry
                    .support
                    .add_weighted(totals.observed_count, filtered_scale_sums[0], 1.0);
                entry.depth += filtered_depths[0];
                entry.mapped_len += filtered_maplens[0];
            } else if !filtered_taxids.is_empty() {
                merge_ambiguous_totals(
                    &mut em_groups,
                    filtered_taxids,
                    AmbiguousTotals {
                        observed_count: totals.observed_count,
                        depths: filtered_depths,
                        mapped_lens: filtered_maplens,
                        equivalent_support_sums: filtered_scale_sums,
                    },
                );
            }
            continue;
        }

        let taxid_set = taxids.iter().copied().collect::<HashSet<_>>();
        let Some((lca_taxid, _lca_name, _lca_rank)) = taxonomy.get_majority_lca(&taxid_set, 0.7)
        else {
            continue;
        };
        let lca_breadth = taxids
            .iter()
            .filter_map(|taxid| coverage.get(taxid))
            .map(|value| value.breadth)
            .fold(0.0, f64::max);
        let lca_fixed = taxids
            .iter()
            .filter_map(|taxid| coverage.get(taxid))
            .map(|value| value.fixed_chunk)
            .fold(0.0, f64::max);
        let lca_flex = taxids
            .iter()
            .filter_map(|taxid| coverage.get(taxid))
            .map(|value| value.flex_chunk)
            .fold(0.0, f64::max);
        coverage.insert(
            lca_taxid,
            TaxonCoverage {
                breadth: lca_breadth,
                fixed_chunk: lca_fixed,
                flex_chunk: lca_flex,
            },
        );
        let candidate_count = totals.depths.len() as f64;
        let depth_mean = totals.depths.iter().sum::<f64>() / candidate_count;
        let maplen_mean = totals.mapped_lens.iter().sum::<f64>() / candidate_count;
        let equivalent_mean = totals.equivalent_support_sums.iter().sum::<f64>() / candidate_count;
        let entry = unambiguous.entry(lca_taxid).or_default();
        entry
            .support
            .add_weighted(totals.observed_count, equivalent_mean, 1.0);
        entry.depth += depth_mean;
        entry.mapped_len += maplen_mean;
    }

    let mut em_groups = em_groups.iter().collect::<Vec<_>>();
    em_groups.sort_unstable_by(|left, right| left.0.cmp(right.0));
    Ok(shared_read_em(
        &unambiguous,
        &em_groups,
        max_iterations,
        "Raw-profile",
        verbose_logging,
    )
    .totals)
}

/// Result of [`shared_read_em`]: the final totals and the totals its last
/// round started from, whose shares classify each group's reads.
struct SharedReadEm {
    totals: HashMap<u32, TaxonTotals>,
    previous: HashMap<u32, TaxonTotals>,
}

/// How a group of shared reads is divided among its candidates.
enum Shares {
    /// One candidate holds at least 0.99 of the group's depth and takes it all.
    Whole(usize),
    /// Each candidate takes the share of its depth.
    Split(Vec<f64>),
}

/// Shares of a group's candidates in proportion to their depth in `totals`,
/// or `None` when they have no depth.
fn depth_shares(candidates: &[u32], totals: &HashMap<u32, TaxonTotals>) -> Option<Shares> {
    let depth = |taxid: &u32| totals.get(taxid).map_or(0.0, |totals| totals.depth);
    let sum_depth = candidates.iter().map(depth).sum::<f64>();
    if sum_depth <= f64::EPSILON || !sum_depth.is_finite() {
        return None;
    }
    let shares = candidates
        .iter()
        .map(|taxid| depth(taxid) / sum_depth)
        .collect::<Vec<_>>();
    let (top, &largest) = shares
        .iter()
        .enumerate()
        .max_by(|left, right| left.1.total_cmp(right.1))?;
    Some(if largest >= 0.99 {
        Shares::Whole(top)
    } else {
        Shares::Split(shares)
    })
}

/// Add candidate `index` of a group, weighted by `share`, to its totals.
fn add_share(
    totals: &mut HashMap<u32, TaxonTotals>,
    taxid: u32,
    group: &AmbiguousTotals,
    index: usize,
    share: f64,
) {
    let entry = totals.entry(taxid).or_default();
    entry.support.add_weighted(
        group.observed_count,
        group.equivalent_support_sums[index],
        share,
    );
    entry.depth += group.depths[index] * share;
    entry.mapped_len += group.mapped_lens[index] * share;
}

/// Divide the reads of each group shared by several taxa among its candidates
/// in proportion to their depth, starting from the reads each taxon has alone
/// (`base`), until the summed depth changes by at most 1e-5 (relative) or
/// `max_iterations` rounds have run.
fn shared_read_em(
    base: &HashMap<u32, TaxonTotals>,
    groups: &[(&Vec<u32>, &AmbiguousTotals)],
    max_iterations: usize,
    label: &str,
    verbose: bool,
) -> SharedReadEm {
    let mut taxids = base.keys().copied().collect::<Vec<_>>();
    for (candidates, _) in groups {
        taxids.extend(candidates.iter().copied());
    }
    taxids.sort_unstable();
    taxids.dedup();
    let depth = |totals: &HashMap<u32, TaxonTotals>, taxid: &u32| {
        totals.get(taxid).map_or(0.0, |totals| totals.depth)
    };
    let max_iterations = max_iterations.max(1);
    let mut current = base.clone();
    let mut previous = HashMap::new();
    for iteration in 0..max_iterations {
        let mut next = base.clone();
        for &(candidates, group) in groups {
            match depth_shares(candidates, &current) {
                Some(Shares::Whole(top)) => add_share(&mut next, candidates[top], group, top, 1.0),
                Some(Shares::Split(shares)) => {
                    for (index, share) in shares.into_iter().enumerate() {
                        add_share(&mut next, candidates[index], group, index, share);
                    }
                }
                None => {}
            }
        }
        if verbose && iteration % 5 == 0 {
            log::info!(target: "PROFILER", "{label} EM iteration {iteration} ...");
        }
        let change = taxids
            .iter()
            .map(|taxid| (depth(&next, taxid) - depth(&current, taxid)).abs())
            .sum::<f64>();
        let scale = taxids
            .iter()
            .map(|taxid| depth(&current, taxid))
            .sum::<f64>()
            .max(1e-12);
        previous = std::mem::replace(&mut current, next);
        if change / scale <= 1e-5 {
            if verbose {
                log::info!(target: "PROFILER", "{label} EM finished after {iteration} iterations.");
            }
            return SharedReadEm {
                totals: current,
                previous,
            };
        }
    }
    log::warn!(
        target: "PROFILER",
        "{label} EM reached the maximum of {max_iterations} iterations before convergence."
    );
    SharedReadEm {
        totals: current,
        previous,
    }
}

/// Shares of each group's reads among its candidates, from the totals the
/// last EM round started from, for the read classifications.
fn group_shares(
    groups: &[(&Vec<u32>, &AmbiguousTotals)],
    totals: &HashMap<u32, TaxonTotals>,
) -> HashMap<Vec<u32>, HashMap<u32, f64>> {
    let mut fractions = HashMap::with_capacity(groups.len());
    for &(candidates, _) in groups {
        let shares = match depth_shares(candidates, totals) {
            Some(Shares::Whole(top)) => HashMap::from([(candidates[top], 1.0)]),
            Some(Shares::Split(shares)) => {
                let mut map = candidates
                    .iter()
                    .zip(shares)
                    .map(|(&taxid, share)| (taxid, round5(share)))
                    .collect();
                normalize_map(&mut map);
                map
            }
            None => continue,
        };
        fractions.insert(candidates.clone(), shares);
    }
    fractions
}

fn merge_ambiguous_totals(
    groups: &mut HashMap<Vec<u32>, AmbiguousTotals>,
    taxids: Vec<u32>,
    totals: AmbiguousTotals,
) {
    match groups.entry(taxids) {
        Entry::Vacant(entry) => {
            entry.insert(totals);
        }
        Entry::Occupied(mut entry) => {
            let current = entry.get_mut();
            current.observed_count += totals.observed_count;
            accumulate(&mut current.depths, &totals.depths);
            accumulate(&mut current.mapped_lens, &totals.mapped_lens);
            accumulate(
                &mut current.equivalent_support_sums,
                &totals.equivalent_support_sums,
            );
        }
    }
}

const TAXON_PROB_CUTOFF: f64 = 1e-5;

/// Probability filter for final taxa rows.
fn pvalue_filter_passes(
    rank: &str,
    cov_prob: Option<f64>,
    chunk_prob: Option<f64>,
    strict: bool,
) -> bool {
    match (cov_prob, chunk_prob) {
        (Some(cov), Some(chunk)) => {
            if strict {
                cov >= TAXON_PROB_CUTOFF && chunk >= TAXON_PROB_CUTOFF
            } else {
                !(cov < TAXON_PROB_CUTOFF && chunk < TAXON_PROB_CUTOFF)
            }
        }
        // Legacy coverage probabilities are only computed for species rows, so
        // for every higher rank their absence is expected rather than evidence
        // of anything, and strict mode must not reject a test that was never
        // run. A species row missing them is a different matter and still
        // fails.
        (None, None) if strict && rank != "species" => true,
        _ => !strict,
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Element-wise accumulation: target[i] += values[i].
/// Extends target vector if values is longer.
fn accumulate(target: &mut Vec<f64>, values: &[f64]) {
    if target.len() < values.len() {
        target.resize(values.len(), 0.0);
    }
    for (t, v) in target.iter_mut().zip(values.iter()) {
        *t += *v;
    }
}

/// Normalize map values to sum to 1.0.
/// Used for converting counts to fractions.
fn normalize_map(map: &mut HashMap<u32, f64>) {
    let mut taxids = map.keys().copied().collect::<Vec<_>>();
    taxids.sort_unstable();
    let sum: f64 = taxids.iter().map(|taxid| map[taxid]).sum::<f64>();
    if sum <= f64::EPSILON {
        return;
    }
    for value in map.values_mut() {
        *value = round5(*value / sum);
    }
}

fn community_profile_row_passes(
    row: &ProfileRow,
    cfg: &ProfilerConfig,
    min_chunk_breadth: f64,
    evidence_selected_taxa: &HashSet<u32>,
    rolled_up_genera: &HashSet<u32>,
    lca_ratio_lower_bounds: &HashMap<u32, (f64, f64)>,
) -> bool {
    let min_reads_ok = equivalent_read_count_passes_minimum(row.equivalent_mapcount, cfg.min_reads);
    if evidence_selected_taxa.contains(&row.taxid) || rolled_up_genera.contains(&row.taxid) {
        // The evidence tiers, or the genus roll-up, chose this taxon; cutoffs
        // the user set explicitly still apply to it, including to rescued taxa.
        return min_reads_ok && explicit_cutoffs_pass(row, cfg);
    }
    let fixed_chunk_ok = row.fixed_chunk > min_chunk_breadth;
    let ratio_ok = if let Some(default_lowers) = lca_ratio_lower_bounds.get(&row.taxid).copied() {
        ratio_filter_passes_with_defaults(
            row.breadth,
            row.expected_breadth,
            row.flex_chunk,
            row.expected_flex_chunk,
            default_lowers,
            cfg.min_oebr,
            cfg.min_coebr,
        )
    } else {
        ratio_filter_passes(
            &row.rank,
            row.breadth,
            row.expected_breadth,
            row.flex_chunk,
            row.expected_flex_chunk,
            cfg.min_oebr,
            cfg.min_coebr,
        )
    };
    min_reads_ok
        && fixed_chunk_ok
        && ratio_ok
        && pvalue_filter_passes(&row.rank, row.breadth_pvalue, row.chunk_pvalue, cfg.strict)
}

/// Remove reported species whose reads a stronger relative explains better:
/// a reported species of the same genus with more depth whose candidate reads
/// are at least `gap` more identical. Reads from genome regions missing from
/// the index land on relatives at lower identity, so such a species mostly
/// consists of the relative's reads. A species truly present next to a much
/// more abundant relative can be removed as well, because the relative's reads
/// also land on it and lower its mean identity. Returns the removed species,
/// each with the relative that explains it.
fn remove_species_explained_by_relatives(
    rows: &mut Vec<ProfileRow>,
    taxonomy: &Taxonomy,
    evidence: &HashMap<u32, TaxonEvidence>,
    gap: f64,
) -> Vec<(u32, u32)> {
    // (taxid, depth, mean identity) of the reported species of each genus.
    // Species without a genus in the taxonomy are not checked.
    let mut genera: HashMap<u32, Vec<(u32, f64, f64)>> = HashMap::new();
    for row in rows.iter().filter(|row| row.rank == "species") {
        let genus = taxonomy
            .get_parents(row.taxid)
            .into_iter()
            .find(|(_, _, rank)| rank == "genus")
            .map(|(genus, _, _)| genus);
        let identity = evidence
            .get(&row.taxid)
            .and_then(TaxonEvidence::mean_identity);
        if let (Some(genus), Some(identity)) = (genus, identity) {
            genera
                .entry(genus)
                .or_default()
                .push((row.taxid, row.depth, identity));
        }
    }
    let mut removed = Vec::new();
    for species in genera.values() {
        for &(taxid, depth, identity) in species {
            let explaining = species.iter().find(|&&(_, other_depth, other_identity)| {
                other_depth > depth && other_identity - identity >= gap
            });
            if let Some(&(relative, ..)) = explaining {
                removed.push((taxid, relative));
            }
        }
    }
    let explained: HashSet<u32> = removed.iter().map(|&(taxid, _)| taxid).collect();
    rows.retain(|row| !explained.contains(&row.taxid));
    removed.sort_unstable();
    removed
}

/// Breadth, chunk breadth and ratio cutoffs given on the command line.
fn explicit_cutoffs_pass(row: &ProfileRow, cfg: &ProfilerConfig) -> bool {
    let at_least = |cutoff: Option<f64>, value: f64| cutoff.is_none_or(|cutoff| value >= cutoff);
    at_least(cfg.breadth, row.breadth)
        && cfg
            .chunk_breadth
            .is_none_or(|cutoff| row.fixed_chunk > cutoff)
        && at_least(cfg.min_oebr, row.breadth / row.expected_breadth)
        && at_least(cfg.min_coebr, row.flex_chunk / row.expected_flex_chunk)
}

/// Move the shares of relabelled taxa in ambiguous classifications to their
/// new taxid.
fn relabel_ambiguous_fractions(
    fractions: &mut HashMap<Vec<u32>, HashMap<u32, f64>>,
    relabel: &HashMap<u32, u32>,
) {
    if relabel.is_empty() {
        return;
    }
    for shares in fractions.values_mut() {
        for (&from, &to) in relabel {
            if let Some(share) = shares.remove(&from) {
                *shares.entry(to).or_default() += share;
            }
        }
    }
}

fn write_community_classifications(
    path: Option<&str>,
    resolved: &ResolvedClassifications,
    ambiguous_counts: &HashMap<Vec<u32>, AmbiguousTotals>,
    ambiguous_reads: &HashMap<Vec<u32>, Vec<String>>,
    ambiguous_fraction: &HashMap<Vec<u32>, HashMap<u32, f64>>,
    retained_taxa: &HashSet<u32>,
    terminal_ranks: &HashMap<u32, ProfileRank>,
    taxonomy: &Taxonomy,
) -> Result<usize> {
    let mut writer = path
        .map(|path| {
            File::create(path)
                .with_context(|| format!("failed to create {path}"))
                .map(BufWriter::new)
        })
        .transpose()?;
    let mut classified = 0usize;

    for group in &resolved.groups {
        if retained_taxa.contains(&group.key.taxid) {
            classified = classified.saturating_add(group.count);
        }
    }

    // Classification rows are optional and potentially numerous. Sort only
    // when they are requested, and sort borrowed pointers rather than cloning
    // read names. The normal count-only path above remains allocation-free.
    if let Some(writer) = writer.as_mut() {
        for group in resolved.sorted_groups() {
            if !retained_taxa.contains(&group.key.taxid) {
                continue;
            }
            let rank = terminal_ranks
                .get(&group.key.taxid)
                .copied()
                .map(ProfileRank::as_str)
                .or_else(|| taxonomy.get_rank(group.key.taxid))
                .unwrap_or("no rank");
            let sources = match &group.key.source {
                ClassificationSource::Direct => group.key.taxid.to_string(),
                ClassificationSource::Candidates(taxids) => taxids
                    .iter()
                    .map(u32::to_string)
                    .collect::<Vec<_>>()
                    .join(";"),
            };
            for read in group.sorted_reads() {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t1",
                    read, group.key.taxid, rank, sources
                )?;
            }
        }
    }

    let mut keys = ambiguous_counts.keys().collect::<Vec<_>>();
    keys.sort_unstable();
    for taxids in keys {
        let Some(fractions) = ambiguous_fraction.get(taxids) else {
            continue;
        };
        let mut retained = fractions
            .iter()
            .filter_map(|(&taxid, &fraction)| {
                (retained_taxa.contains(&taxid) && fraction.is_finite() && fraction > 0.0)
                    .then_some((taxid, fraction))
            })
            .collect::<Vec<_>>();
        retained.sort_unstable_by_key(|(taxid, _)| *taxid);
        let total = retained.iter().map(|(_, fraction)| *fraction).sum::<f64>();
        if total <= f64::EPSILON || !total.is_finite() {
            continue;
        }
        for (_, fraction) in &mut retained {
            *fraction /= total;
        }
        retained.sort_by(|left, right| {
            right
                .1
                .total_cmp(&left.1)
                .then_with(|| left.0.cmp(&right.0))
        });
        let count = ambiguous_counts[taxids].observed_count.round() as usize;
        classified = classified.saturating_add(count);

        let Some(writer) = writer.as_mut() else {
            continue;
        };
        let final_taxid = retained[0].0;
        let rank = terminal_ranks
            .get(&final_taxid)
            .copied()
            .map(ProfileRank::as_str)
            .or_else(|| taxonomy.get_rank(final_taxid))
            .unwrap_or("no rank");
        let candidate_taxids = retained
            .iter()
            .map(|(taxid, _)| taxid.to_string())
            .collect::<Vec<_>>()
            .join(";");
        let candidate_fractions = retained
            .iter()
            .map(|(_, fraction)| format!("{fraction:.5}"))
            .collect::<Vec<_>>()
            .join(";");
        if let Some(reads) = ambiguous_reads.get(taxids) {
            let mut reads = reads.iter().map(String::as_str).collect::<Vec<_>>();
            reads.sort_unstable();
            for read in reads {
                writeln!(
                    writer,
                    "{read}\t{final_taxid}\t{rank}\t{candidate_taxids}\t{candidate_fractions}"
                )?;
            }
        }
    }

    if let Some(writer) = writer.as_mut() {
        writer.flush()?;
    }
    Ok(classified)
}

/// Aggregated data from mapping loading.
///
/// Contains all information needed for the profiling phase.
struct MappingData {
    /// Per-read alignment bundles
    read_alignments: Vec<ReadBundle>,

    /// Total assembly lengths: asm_id -> total_length. Shared with the
    /// worker context via `Arc` so no clone is needed when threading
    /// it through the profile / assign_taxonomy / write_outputs chain.
    asm_lengths: Arc<HashMap<String, f64>>,

    /// Minimum breadth threshold for reporting taxa
    min_breadth: f64,

    /// Effective standard chunk-breadth threshold. In deep recall-mode
    /// community runs this may be lower than the original automatic cutoff.
    min_chunk_breadth: f64,

    /// Original unscaled cutoff retained by rescue and guarded-MAP safeguards.
    safeguard_min_chunk_breadth: f64,

    /// Whether at least one genome has a full/effective size multiplier > 1.
    has_scaled_genomes: bool,

    /// Total reads in the input.
    total_reads: usize,

    /// Distinct reads that produced at least one retained alignment.
    aligned_reads: usize,

    /// Run-level mapping rate = aligned_reads / total_reads.
    /// Used to gate the genus-level fallback mode.
    mapping_rate: f64,

    /// Whether genus-level fallback assignment is active for this run.
    /// Automatic activation requires low-map Illumina non-pathogen data on a
    /// non-subsampled DB; explicit `--genus-fallback` activates it directly.
    fallback_mode: bool,

    /// Mapping from species taxid to its enclosing genus taxid, derived from
    /// the taxonomy lineage (not parsed from reference headers). Keyed by
    /// the **post-normalization** species taxid so that sub-species
    /// aliases all map through the same entry.
    /// Species without a resolvable genus are simply absent from this map.
    species_to_genus: Arc<HashMap<u32, u32>>,

    /// Most-specific reference-supported rank for every terminal taxid. Other
    /// rows (genus fallback, LCA) take their rank from the taxonomy.
    terminal_ranks: Arc<HashMap<u32, ProfileRank>>,

    /// Mapping records discarded at load because their reference had no
    /// parsable species taxid or assembly length. Seeds
    /// [`FilterStats::skipped_reference`].
    skipped_reference_at_load: usize,

    /// Marker genes and fragment geometry when the database has marker genomes.
    marker_model: Option<Arc<MarkerModel>>,

    /// Inputs of the genus roll-up; `None` when it is off (`--strict-coverage`).
    genus_rollup: Option<Arc<GenusRollup>>,

    /// Inputs of `--estimate-unclassified`; `None` when it is off.
    unclassified_species: Option<MarkerSpecies>,

    /// Genome accession of each assembly id.
    assembly_names: Vec<Arc<str>>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profiler::evidence::EvidenceCoverage;
    use crate::profiler::{ExactTiePolicy, ProfileMode};
    use crate::sequencer::Sequencer;

    fn approx_equal(lhs: f64, rhs: f64) {
        assert!(
            (lhs - rhs).abs() < 1e-9,
            "left {} not approximately equal to right {}",
            lhs,
            rhs
        );
    }

    fn placement(assembly_id: u32, mapped: u32, read_len: u32) -> CandidatePlacement {
        CandidatePlacement {
            taxid: 1,
            profile_rank: ProfileRank::Species,
            reference_id: 0,
            assembly_id,
            sampling_scale: 1.0,
            asm_len: 2_000.0,
            r1_mapped: mapped,
            both_mapped: mapped,
            read_len,
            other_genes: 0,
            pair_errors: 0.0,
            long_score: 0.0,
            span: (1, 1 + mapped),
            mate_span: (0, 0),
            fixed_chunk: (1, 1 + mapped),
            flex_chunk: (1, 1 + mapped),
            input_order: 0,
        }
    }

    #[test]
    fn a_long_read_adds_depth_for_every_marker_gene_it_covers() {
        // Long reads of 3,000 bases on genes of 1,000 bases: 500 bases inside
        // a gene, 1,000 - 1,000 + 3,001 = 3,001 placements per gene.
        let geometry = FragmentGeometry::new(FragmentLayout::LongReads, 3_000, 3_000, 500, 0.2);
        let markers = MarkerModel::from_gene_lengths(&[(7, &[1_000, 1_000])], geometry);
        let mut read = placement(7, 900, 3_000);
        approx_equal(fragment_depth(read, 900, Some(&markers)), 3_000.0 / 6_002.0);
        read.other_genes = 1;
        approx_equal(fragment_depth(read, 900, Some(&markers)), 6_000.0 / 6_002.0);
        // A read of another length is placed by its own length.
        let short = placement(7, 900, 1_000);
        approx_equal(
            fragment_depth(short, 900, Some(&markers)),
            1_000.0 / (2.0 * (1_000.0 + 1_001.0 - 400.0)),
        );
    }

    #[test]
    fn a_fragment_on_a_marker_gene_adds_its_read_over_the_genes_placements() {
        // 150 bp pairs of 400 bases with lone mates: genes of 1,000 bases
        // offer 1,191 placements each; a gene of 100 bases holds no read.
        let geometry = FragmentGeometry::new(FragmentLayout::PairsAndLoneMates, 150, 400, 50, 0.7);
        let markers =
            MarkerModel::from_gene_lengths(&[(7, &[1_000, 1_000]), (8, &[100])], geometry);
        let depth = |assembly, mapped, markers| {
            fragment_depth(placement(assembly, mapped, 150), mapped, markers)
        };
        approx_equal(depth(7, 150, Some(&markers)), 150.0 / 2_382.0);
        // A read clipped at a gene end stands for its whole read.
        approx_equal(depth(7, 110, Some(&markers)), 150.0 / 2_382.0);
        // Genome references and marker genes without placements keep mapped
        // bases over the summed length.
        approx_equal(depth(9, 140, Some(&markers)), 140.0 / 2_000.0);
        approx_equal(depth(8, 90, Some(&markers)), 90.0 / 2_000.0);
        approx_equal(depth(7, 140, None), 140.0 / 2_000.0);
    }

    #[test]
    fn expected_breadth_is_clamped_for_high_depths() {
        let expected = crate::profiler::evidence::expected_breadth(10.0, 5.0);
        approx_equal(expected, 1.0);
        assert!(expected.is_finite());

        let zero = crate::profiler::evidence::expected_breadth(10.0, 0.0);
        approx_equal(zero, 0.0);
    }

    // ------------------------------------------------------------------
    // Helpers for fallback / rank-aware tests.
    // ------------------------------------------------------------------

    /// Minimal `ProfilerConfig` for unit tests. Callers mutate only the
    /// fields relevant to the scenario under test.
    fn make_cfg(sequencer: Sequencer, host: Option<&str>) -> ProfilerConfig {
        ProfilerConfig {
            mapper_mappings: None,
            mapping_profile_memory: None,
            db_metadata: None,
            sequencer,
            batch_size: 1,
            max_em_iterations: 200,
            lenient_mapping_input: false,
            is_paired: false,
            exact_ties: ExactTiePolicy::Ambiguous,
            identity: 0.95,
            identity_check_accuracy: None,
            mapped_len: 50,
            breadth: Some(0.0),
            chunk_breadth: Some(0.0),
            chunk_breadth_scale: 1.0,
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            min_marker_genes: 3,
            relative_identity_gap: None,
            removed_to_genus: false,
            estimate_unclassified: false,
            identity_margin: None,
            fraction: 0.6,
            lowbiomass: false,
            keep_raw: false,
            save_classifications: false,
            by_aligned: false,
            genus_identity: 0.80,
            low_map_rate_threshold: 0.30,
            genus_fallback: false,
            strict: false,
            pathogen_host: None,
            host: host.map(|s| s.to_string()),
            threads: 1,
            verbose: false,
            very_verbose: false,
            progress: crate::progress::ProgressOptions::default(),
            suppress_mapping_load_progress: false,
            outprefix: "test".to_string(),
            mode: ProfileMode::Default,
            dmp_dir: None,
            strain: false,
        }
    }

    fn test_taxonomy(directory: &Path) -> Arc<Taxonomy> {
        let nodes = directory.join("nodes.dmp");
        let names = directory.join("names.dmp");
        let merged = directory.join("merged.dmp");
        std::fs::write(
            &nodes,
            "1\t|\t1\t|\tno rank\t|\n101\t|\t1\t|\tspecies\t|\n102\t|\t1\t|\tspecies\t|\n",
        )
        .unwrap();
        std::fs::write(
            &names,
            "1\t|\troot\t|\t\t|\tscientific name\t|\n101\t|\tSpecies one\t|\t\t|\tscientific name\t|\n102\t|\tSpecies two\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();
        std::fs::write(&merged, "").unwrap();
        Arc::new(Taxonomy::from_files(&nodes, &names, &merged, None).unwrap())
    }

    fn raw_output_fixture() -> (MappingData, AssignmentResult) {
        let good_key = AccKey {
            taxid: 101,
            asm: Arc::<str>::from("asm101"),
            acc: Arc::<str>::from("ref101"),
        };
        let rejected_key = AccKey {
            taxid: 102,
            asm: Arc::<str>::from("asm102"),
            acc: Arc::<str>::from("ref102"),
        };
        let assignment = AssignmentResult {
            read_assignments: vec![
                ReadAssignment {
                    read: "accepted-read".into(),
                    rank: AssignRank::Species,
                    requires_terminal_evidence: false,
                    taxids: vec![101],
                    depths: vec![0.1],
                    mapped_lens: vec![100.0],
                    sampling_scales: vec![1.0],
                    best: 0,
                },
                ReadAssignment {
                    read: "rejected-read".into(),
                    rank: AssignRank::Species,
                    requires_terminal_evidence: false,
                    taxids: vec![102],
                    depths: vec![0.1],
                    mapped_lens: vec![100.0],
                    sampling_scales: vec![1.0],
                    best: 0,
                },
            ],
            acc_intervals: HashMap::from([
                (good_key.clone(), vec![(0, 100)]),
                (rejected_key.clone(), vec![(0, 10)]),
            ]),
            acc_fixed_chunk_intervals: HashMap::from([
                (good_key.clone(), vec![(0, 100)]),
                (rejected_key.clone(), vec![(0, 10)]),
            ]),
            acc_flex_chunk_intervals: HashMap::from([
                (good_key, vec![(0, 67)]),
                (rejected_key, vec![(0, 10)]),
            ]),
            ..AssignmentResult::default()
        };
        let mapping_data = MappingData {
            read_alignments: Vec::new(),
            asm_lengths: Arc::new(HashMap::from([
                ("asm101".to_string(), 1_000.0),
                ("asm102".to_string(), 1_000.0),
            ])),
            min_breadth: 0.05,
            min_chunk_breadth: 0.05,
            safeguard_min_chunk_breadth: 0.05,
            has_scaled_genomes: false,
            total_reads: 2,
            aligned_reads: 2,
            mapping_rate: 1.0,
            fallback_mode: false,
            species_to_genus: Arc::new(HashMap::new()),
            terminal_ranks: Arc::new(HashMap::new()),
            skipped_reference_at_load: 0,
            marker_model: None,
            genus_rollup: None,
            unclassified_species: None,
            assembly_names: Vec::new(),
        };
        (mapping_data, assignment)
    }

    #[test]
    fn raw_em_refines_all_quality_qualified_candidates() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let assignments = vec![
            ReadAssignment {
                read: "unique-one".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: false,
                taxids: vec![101],
                depths: vec![1.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
                best: 0,
            },
            ReadAssignment {
                read: "unique-two".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: false,
                taxids: vec![102],
                depths: vec![3.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
                best: 0,
            },
            ReadAssignment {
                read: "ambiguous".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: false,
                taxids: vec![101, 102],
                depths: vec![2.0, 2.0],
                mapped_lens: vec![100.0, 100.0],
                sampling_scales: vec![1.0, 1.0],
                best: 0,
            },
        ];
        let totals = estimate_unfiltered_taxon_totals(
            &assignments,
            &taxonomy,
            &mut HashMap::new(),
            200,
            false,
        )
        .unwrap();

        approx_equal(totals[&101].support.observed, 1.25);
        approx_equal(totals[&102].support.observed, 1.75);
        approx_equal(totals[&101].depth, 1.5);
        approx_equal(totals[&102].depth, 4.5);
    }

    #[test]
    fn raw_em_distributes_across_mixed_terminal_ranks_without_parent_rollup() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let mut assignments = vec![ReadAssignment {
            read: "unique-genus".into(),
            rank: AssignRank::Genus,
            requires_terminal_evidence: true,
            taxids: vec![111],
            depths: vec![2.0],
            mapped_lens: vec![100.0],
            sampling_scales: vec![1.0],
            best: 0,
        }];
        for taxid in [1111, 1112, 1113] {
            assignments.push(ReadAssignment {
                read: format!("unique-species-{taxid}").into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![taxid],
                depths: vec![1.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
                best: 0,
            });
        }
        assignments.push(ReadAssignment {
            read: "mixed-rank-ambiguous".into(),
            rank: AssignRank::Mixed,
            requires_terminal_evidence: true,
            taxids: vec![111, 1111, 1112, 1113],
            depths: vec![1.0; 4],
            mapped_lens: vec![100.0; 4],
            sampling_scales: vec![1.0; 4],
            best: 0,
        });

        let totals = estimate_unfiltered_taxon_totals(
            &assignments,
            &taxonomy,
            &mut HashMap::new(),
            200,
            false,
        )
        .unwrap();

        approx_equal(totals[&111].depth, 2.4);
        for taxid in [1111, 1112, 1113] {
            approx_equal(totals[&taxid].depth, 1.2);
        }
    }

    #[test]
    fn raw_em_rejects_inconsistent_candidate_vectors() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let assignments = vec![ReadAssignment {
            read: "broken".into(),
            rank: AssignRank::Species,
            requires_terminal_evidence: false,
            taxids: vec![101, 102],
            depths: vec![1.0],
            mapped_lens: vec![100.0, 100.0],
            sampling_scales: vec![1.0, 1.0],
            best: 0,
        }];
        let error = estimate_unfiltered_taxon_totals(
            &assignments,
            &taxonomy,
            &mut HashMap::new(),
            200,
            false,
        )
        .unwrap_err();
        assert!(error
            .to_string()
            .contains("inconsistent raw-profile candidate vectors"));
    }

    #[test]
    fn raw_rows_retain_rejected_native_witness_measurements() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let profiler =
            CommunityProfiler::new(make_cfg(Sequencer::Illumina, None), taxonomy).unwrap();
        let estimate = RawProfileEstimate {
            totals: HashMap::from([(
                102,
                TaxonTotals {
                    support: ReadSupport {
                        observed: 0.25,
                        equivalent: 0.25,
                    },
                    depth: 0.02,
                    mapped_len: 20.0,
                },
            )]),
            coverage: HashMap::new(),
        };
        let decisions = HashMap::from([(
            102,
            WitnessDecision {
                metrics: EvidenceCoverage {
                    breadth: 0.01234,
                    expected_breadth: 0.05678,
                    fixed_chunk: 0.02345,
                    flex_chunk: 0.03456,
                    expected_flex_chunk: 0.06789,
                    breadth_pvalue: Some(0.0002),
                    chunk_pvalue: Some(0.25),
                    ..EvidenceCoverage::default()
                },
                standard: false,
                rescues: RescueKinds::default(),
                required_premap_guard: PremapChunkGuard::None,
            },
        )]);

        let rows =
            profiler.build_raw_profile_rows(&estimate, &decisions, &HashMap::new(), None, None);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].taxid, 102);
        assert_eq!(rows[0].equivalent_mapcount, 0.25);
        assert_eq!(rows[0].abundance, 100.0);
        assert_eq!(rows[0].breadth, 0.01234);
        assert_eq!(rows[0].expected_breadth, 0.05678);
        assert_eq!(rows[0].fixed_chunk, 0.02345);
        assert_eq!(rows[0].flex_chunk, 0.03456);
        assert_eq!(rows[0].expected_flex_chunk, 0.06789);
    }

    #[test]
    fn enabling_raw_profile_does_not_change_production_outputs() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let baseline_prefix = directory.path().join("baseline");
        let raw_prefix = directory.path().join("with-raw");

        let mut baseline_cfg = make_cfg(Sequencer::Illumina, None);
        baseline_cfg.outprefix = baseline_prefix.to_string_lossy().into_owned();
        baseline_cfg.save_classifications = true;
        baseline_cfg.breadth = Some(0.05);
        baseline_cfg.chunk_breadth = Some(0.05);
        let baseline = CommunityProfiler::new(baseline_cfg, Arc::clone(&taxonomy)).unwrap();
        let (mapping_data, assignment) = raw_output_fixture();
        assert_eq!(
            baseline
                .write_outputs(
                    mapping_data,
                    assignment,
                    EvidencePolicy::from_config(&baseline.cfg).unwrap()
                )
                .unwrap(),
            1
        );

        let mut raw_cfg = make_cfg(Sequencer::Illumina, None);
        raw_cfg.outprefix = raw_prefix.to_string_lossy().into_owned();
        raw_cfg.save_classifications = true;
        raw_cfg.keep_raw = true;
        raw_cfg.breadth = Some(0.05);
        raw_cfg.chunk_breadth = Some(0.05);
        let with_raw = CommunityProfiler::new(raw_cfg, taxonomy).unwrap();
        let (mapping_data, assignment) = raw_output_fixture();
        assert_eq!(
            with_raw
                .write_outputs(
                    mapping_data,
                    assignment,
                    EvidencePolicy::from_config(&with_raw.cfg).unwrap()
                )
                .unwrap(),
            1
        );

        let baseline_profile =
            std::fs::read(format!("{}.profile.txt", baseline_prefix.to_string_lossy())).unwrap();
        let raw_profile =
            std::fs::read(format!("{}.profile.txt", raw_prefix.to_string_lossy())).unwrap();
        assert_eq!(baseline_profile, raw_profile);

        let baseline_classify = std::fs::read(format!(
            "{}.classify.txt",
            baseline_prefix.to_string_lossy()
        ))
        .unwrap();
        let raw_classify =
            std::fs::read(format!("{}.classify.txt", raw_prefix.to_string_lossy())).unwrap();
        assert_eq!(baseline_classify, raw_classify);
        let classify_text = String::from_utf8(raw_classify).unwrap();
        assert!(classify_text.contains("accepted-read\t101"));
        assert!(!classify_text.contains("rejected-read"));

        let raw_table =
            std::fs::read_to_string(format!("{}.rprofile.txt", raw_prefix.to_string_lossy()))
                .unwrap();
        assert!(raw_table.contains("Species one\t101\t"));
        assert!(raw_table.contains("Species two\t102\t"));
    }

    #[test]
    fn final_minimum_support_controls_profile_classification_and_count_together() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let output = directory.path().join("minimum-support");
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.outprefix = output.to_string_lossy().into_owned();
        cfg.save_classifications = true;
        cfg.min_reads = Some(2);
        let profiler = CommunityProfiler::new(cfg, taxonomy).unwrap();
        let (mapping_data, assignment) = raw_output_fixture();

        assert_eq!(
            profiler
                .write_outputs(
                    mapping_data,
                    assignment,
                    EvidencePolicy::from_config(&profiler.cfg).unwrap()
                )
                .unwrap(),
            0
        );
        assert!(
            std::fs::read(format!("{}.profile.txt", output.to_string_lossy()))
                .unwrap()
                .is_empty()
        );
        assert!(
            std::fs::read(format!("{}.classify.txt", output.to_string_lossy()))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn species_explained_by_a_more_identical_relative_are_removed() {
        // Family 100 holds genera 10 (species 1, 2, 3) and 11 (species 4);
        // species 5 is in family 101.
        let nodes = [
            (1u32, (1u32, "no rank")),
            (100, (1, "family")),
            (101, (1, "family")),
            (10, (100, "genus")),
            (11, (100, "genus")),
            (12, (101, "genus")),
            (1001, (10, "species")),
            (1002, (10, "species")),
            (1003, (10, "species")),
            (1004, (11, "species")),
            (1005, (12, "species")),
        ]
        .into_iter()
        .map(|(id, (parent, rank))| (id, (parent, rank.to_string())))
        .collect();
        let taxonomy = Taxonomy::from_parts(nodes, HashMap::new(), HashMap::new(), HashSet::new());
        let row = |taxid: u32, depth: f64| ProfileRow {
            taxname: taxid.to_string(),
            taxid,
            rank: "species".to_string(),
            equivalent_mapcount: 10.0,
            depth,
            abundance: depth,
            breadth: 0.1,
            expected_breadth: 0.1,
            breadth_pvalue: Some(1.0),
            fixed_chunk: 0.1,
            flex_chunk: 0.1,
            expected_flex_chunk: 0.1,
            chunk_pvalue: Some(1.0),
            best_fit: None,
            annotations: None,
        };
        let evidence: HashMap<u32, TaxonEvidence> = [
            (1001, 0.998),
            (1002, 0.95),
            (1003, 0.993),
            (1004, 0.96),
            (1005, 0.90),
        ]
        .into_iter()
        .map(|(taxid, identity)| {
            let evidence = TaxonEvidence {
                mapped_fragments: 10,
                identity_sum: 10.0 * identity,
                ..TaxonEvidence::default()
            };
            (taxid, evidence)
        })
        .collect();
        let rows = || {
            vec![
                row(1001, 1.0),
                row(1002, 0.3),
                row(1003, 0.2),
                row(1004, 0.1),
                row(1005, 0.05),
            ]
        };

        // Species 1002 is explained by 1001 of the same genus; 1003 is within
        // the gap, 1004 is in another genus of the family and 1005 is in
        // another family.
        let mut all = rows();
        let removed = remove_species_explained_by_relatives(&mut all, &taxonomy, &evidence, 0.01);
        assert_eq!(removed, [(1002, 1001)]);
        let kept: Vec<u32> = all.iter().map(|row| row.taxid).collect();
        assert_eq!(kept, [1001, 1003, 1004, 1005]);

        // With a gap of 0.05, 1002 (0.048 below 1001) stays as well.
        let mut wide = rows();
        let removed = remove_species_explained_by_relatives(&mut wide, &taxonomy, &evidence, 0.05);
        assert!(removed.is_empty());
    }

    #[test]
    fn community_filter_uses_unrounded_chunk_breadth() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        let row = |fixed_chunk| ProfileRow {
            taxname: "Species one".to_string(),
            taxid: 101,
            rank: "species".to_string(),
            equivalent_mapcount: 10.0,
            depth: 0.1,
            abundance: 1.0,
            breadth: 0.1,
            expected_breadth: 0.1,
            breadth_pvalue: Some(1.0),
            fixed_chunk,
            flex_chunk: 0.1,
            expected_flex_chunk: 0.1,
            chunk_pvalue: Some(1.0),
            best_fit: None,
            annotations: None,
        };
        let selected = HashSet::new();
        let lca = HashMap::new();
        assert!(!community_profile_row_passes(
            &row(0.499_996),
            &cfg,
            0.5,
            &selected,
            &HashSet::new(),
            &lca,
        ));
        assert!(community_profile_row_passes(
            &row(0.500_004),
            &cfg,
            0.5,
            &selected,
            &HashSet::new(),
            &lca,
        ));
    }

    #[test]
    fn unresolved_native_terminal_does_not_enter_legacy_zero_breadth_path() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let output = directory.path().join("unresolved-native");
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.outprefix = output.to_string_lossy().into_owned();
        cfg.breadth = Some(0.0);
        cfg.chunk_breadth = Some(0.0);
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let profiler = CommunityProfiler::new(cfg, taxonomy).unwrap();
        let mapping_data = MappingData {
            read_alignments: Vec::new(),
            asm_lengths: Arc::new(HashMap::from([("asm101".to_string(), 1_000.0)])),
            min_breadth: 0.0,
            min_chunk_breadth: 0.0,
            safeguard_min_chunk_breadth: 0.0,
            has_scaled_genomes: false,
            total_reads: 1,
            aligned_reads: 1,
            mapping_rate: 1.0,
            fallback_mode: false,
            species_to_genus: Arc::new(HashMap::new()),
            terminal_ranks: Arc::new(HashMap::from([(101, ProfileRank::Species)])),
            skipped_reference_at_load: 0,
            marker_model: None,
            genus_rollup: None,
            unclassified_species: None,
            assembly_names: Vec::new(),
        };
        let assignment = AssignmentResult {
            read_assignments: vec![ReadAssignment {
                read: "unresolved".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![101],
                depths: vec![0.1],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
                best: 0,
            }],
            ..AssignmentResult::default()
        };

        assert_eq!(
            profiler
                .write_outputs(mapping_data, assignment, policy)
                .unwrap(),
            0
        );
        assert!(
            std::fs::read(format!("{}.profile.txt", output.to_string_lossy()))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn matching_terminal_taxid_does_not_reclassify_legacy_evidence() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        let output = directory.path().join("legacy-taxid-collision");
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.outprefix = output.to_string_lossy().into_owned();
        cfg.breadth = Some(0.0);
        cfg.chunk_breadth = Some(0.0);
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let profiler = CommunityProfiler::new(cfg, taxonomy).unwrap();
        let key = AccKey {
            taxid: 101,
            asm: Arc::<str>::from("asm101"),
            acc: Arc::<str>::from("ref101"),
        };
        let mapping_data = MappingData {
            read_alignments: Vec::new(),
            asm_lengths: Arc::new(HashMap::from([("asm101".to_string(), 1_000.0)])),
            min_breadth: 0.0,
            min_chunk_breadth: 0.0,
            safeguard_min_chunk_breadth: 0.0,
            has_scaled_genomes: false,
            total_reads: 1,
            aligned_reads: 1,
            mapping_rate: 1.0,
            fallback_mode: true,
            species_to_genus: Arc::new(HashMap::new()),
            terminal_ranks: Arc::new(HashMap::from([(101, ProfileRank::Species)])),
            skipped_reference_at_load: 0,
            marker_model: None,
            genus_rollup: None,
            unclassified_species: None,
            assembly_names: Vec::new(),
        };
        let assignment = AssignmentResult {
            read_assignments: vec![ReadAssignment {
                read: "legacy-fallback".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: false,
                taxids: vec![101],
                depths: vec![0.1],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
                best: 0,
            }],
            acc_intervals: HashMap::from([(key.clone(), vec![(0, 100)])]),
            acc_fixed_chunk_intervals: HashMap::from([(key.clone(), vec![(0, 100)])]),
            acc_flex_chunk_intervals: HashMap::from([(key, vec![(0, 67)])]),
            ..AssignmentResult::default()
        };

        assert_eq!(
            profiler
                .write_outputs(mapping_data, assignment, policy)
                .unwrap(),
            1
        );
        let profile =
            std::fs::read_to_string(format!("{}.profile.txt", output.to_string_lossy())).unwrap();
        assert_eq!(profile.lines().count(), 1);
        assert!(profile.contains("Species one\t101\tspecies\t"));
    }

    #[test]
    fn conflicting_terminal_rank_for_one_taxid_is_rejected() {
        let mut ranks = HashMap::new();
        record_terminal_rank(&mut ranks, 111, ProfileRank::Genus, "genus-ref").unwrap();
        record_terminal_rank(&mut ranks, 111, ProfileRank::Genus, "second-genus-ref").unwrap();
        let error =
            record_terminal_rank(&mut ranks, 111, ProfileRank::Species, "species-ref").unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("conflicting terminal ranks for taxid 111"));
        assert!(message.contains("genus"));
        assert!(message.contains("species"));
    }

    #[test]
    fn equivalent_mapping_rate_counts_the_genome_share_of_fragments() {
        let bundle = |scale: f32| {
            let mut row = native_mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Unpaired);
            row.sampling_scale = scale;
            ReadBundle {
                read_name: Arc::from("read"),
                alignments: vec![row],
            }
        };
        assert_eq!(
            genome_equivalent_mapping_rate(&[bundle(1.0), bundle(1.0)], 4, None),
            0.5
        );
        assert_eq!(
            genome_equivalent_mapping_rate(&[bundle(10.0)], 20, None),
            0.5
        );
        assert_eq!(
            genome_equivalent_mapping_rate(&[bundle(10.0)], 4, None),
            1.0
        );
        assert_eq!(genome_equivalent_mapping_rate(&[], 0, None), 0.0);
    }

    fn native_mapping(
        species: u32,
        reference_id: u32,
        assembly_id: u32,
        identity: f32,
        mapped_len: u32,
        read_len: u32,
        role: MappingRole,
    ) -> AlignmentEntry {
        native_mapping_at_rank(
            species,
            ProfileRank::Species,
            reference_id,
            assembly_id,
            identity,
            mapped_len,
            read_len,
            role,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn native_mapping_at_rank(
        species: u32,
        profile_rank: ProfileRank,
        reference_id: u32,
        assembly_id: u32,
        identity: f32,
        mapped_len: u32,
        read_len: u32,
        role: MappingRole,
    ) -> AlignmentEntry {
        AlignmentEntry {
            asm: Arc::<str>::from(format!("asm{assembly_id}")),
            acc: Arc::<str>::from(format!("ref{reference_id}")),
            profile_taxid: species,
            profile_rank,
            reference_id,
            assembly_id,
            role,
            pair_id: if matches!(role, MappingRole::Mate1 | MappingRole::Mate2) {
                reference_id.max(1)
            } else {
                0
            },
            candidate_boundary: metam::mapper::CandidateBoundary::COMPLETE,
            sampling_scale: 1.0,
            asm_len: 1_000.0,
            identity,
            fraction: mapped_len as f32 / read_len as f32,
            mapped_len,
            read_len,
            span: (100, 100 + mapped_len),
            fixed_chunk: (1, 300 + mapped_len),
            flex_chunk: (1, 500 + mapped_len),
            marker_gene_len: 0,
        }
    }

    fn witness(
        fixed_chunk: f64,
        standard: bool,
        rescue: bool,
        required_premap_guard: PremapChunkGuard,
    ) -> WitnessDecision {
        WitnessDecision {
            metrics: EvidenceCoverage {
                fixed_chunk,
                ..EvidenceCoverage::default()
            },
            standard,
            rescues: RescueKinds {
                platform: rescue,
                ..RescueKinds::default()
            },
            required_premap_guard,
        }
    }

    fn low_abundance_witness(fixed_chunk: f64) -> WitnessDecision {
        WitnessDecision {
            metrics: EvidenceCoverage {
                fixed_chunk,
                ..EvidenceCoverage::default()
            },
            standard: false,
            rescues: RescueKinds {
                low_abundance_pe: true,
                ..RescueKinds::default()
            },
            required_premap_guard: PremapChunkGuard::FivePercent,
        }
    }

    #[test]
    fn guarded_map_requires_independent_coverage_on_the_same_assembly() {
        let species = 101;
        let assembly_one = PresenceKey {
            taxid: species,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        };
        let assembly_two = PresenceKey {
            taxid: species,
            profile_rank: ProfileRank::Species,
            assembly_id: 2,
        };
        let baseline_selected = HashSet::new();
        let mut baseline = HashMap::from([
            (
                assembly_one,
                witness(0.70, false, false, PremapChunkGuard::None),
            ),
            (
                assembly_two,
                witness(0.20, false, false, PremapChunkGuard::None),
            ),
        ]);
        let post = vec![(
            assembly_two,
            witness(0.70, true, false, PremapChunkGuard::Full),
        )];

        let additions = guarded_map_additions(&baseline_selected, &baseline, &post, 0.60);
        assert!(
            additions.is_empty(),
            "coverage on assembly one must not authorize MAP evidence on assembly two"
        );

        baseline.insert(
            assembly_two,
            witness(0.61, false, false, PremapChunkGuard::None),
        );
        let additions = guarded_map_additions(&baseline_selected, &baseline, &post, 0.60);
        assert!(additions.contains_key(&species));
    }

    #[test]
    fn guarded_map_uses_mode_specific_post_map_and_amount_guards() {
        let key = PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        };
        let baseline_selected = HashSet::new();
        let se_rescue_only = vec![(key, witness(0.70, false, true, PremapChunkGuard::None))];

        let se_baseline =
            HashMap::from([(key, witness(0.61, false, false, PremapChunkGuard::None))]);
        assert!(
            guarded_map_additions(&baseline_selected, &se_baseline, &se_rescue_only, 0.60,)
                .is_empty(),
            "SE MAP-only calls must pass the full post-MAP standard tier"
        );

        let pe_rescue_only = vec![(key, witness(0.70, false, true, PremapChunkGuard::Half))];
        let pe_baseline =
            HashMap::from([(key, witness(0.31, false, false, PremapChunkGuard::None))]);
        assert!(
            guarded_map_additions(&baseline_selected, &pe_baseline, &pe_rescue_only, 0.60,)
                .contains_key(&101),
            "PE retains its validated half-cutoff plus rescue-tier guard"
        );
    }

    #[test]
    fn equivalent_assembly_witnesses_use_stable_key_tiebreak() {
        let low_key = PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 3,
        };
        let high_key = PresenceKey {
            assembly_id: 9,
            ..low_key
        };
        let mut low = witness(0.50, true, false, PremapChunkGuard::Half);
        low.metrics.flex_chunk = 0.30;
        let mut high = low;
        high.metrics.flex_chunk = 0.90;

        let forward = strongest_witnesses(&[(high_key, high), (low_key, low)]);
        let reverse = strongest_witnesses(&[(low_key, low), (high_key, high)]);
        assert_eq!(forward[&101].metrics.flex_chunk, 0.30);
        assert_eq!(reverse[&101].metrics.flex_chunk, 0.30);
    }

    #[test]
    fn paired_low_abundance_map_uses_relaxed_same_assembly_guard_only_for_its_tier() {
        let species = 101;
        let assembly_one = PresenceKey {
            taxid: species,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        };
        let assembly_two = PresenceKey {
            assembly_id: 2,
            ..assembly_one
        };
        let selected = HashSet::new();
        let low_post = vec![(assembly_one, low_abundance_witness(0.29))];

        let above_guard = HashMap::from([(
            assembly_one,
            witness(0.031, false, false, PremapChunkGuard::None),
        )]);
        assert!(
            guarded_map_additions(&selected, &above_guard, &low_post, 0.60)
                .contains_key(&species),
            "the calibrated PE low-abundance tier accepts independent pre-MAP coverage above 0.05*C"
        );

        let at_guard = HashMap::from([(
            assembly_one,
            witness(0.030, false, false, PremapChunkGuard::None),
        )]);
        assert!(
            guarded_map_additions(&selected, &at_guard, &low_post, 0.60).is_empty(),
            "the independent-coverage comparison is strict"
        );

        let wrong_assembly = HashMap::from([(
            assembly_two,
            witness(0.50, false, false, PremapChunkGuard::None),
        )]);
        assert!(
            guarded_map_additions(&selected, &wrong_assembly, &low_post, 0.60).is_empty(),
            "coverage from another assembly must not authorize the rescue"
        );

        let ordinary_post = vec![(
            assembly_one,
            witness(0.70, false, true, PremapChunkGuard::Half),
        )];
        assert!(
            guarded_map_additions(&selected, &above_guard, &ordinary_post, 0.60).is_empty(),
            "ordinary PE standard/platform evidence retains its 0.50*C pre-MAP guard"
        );
    }

    #[test]
    fn evidence_policy_is_limited_to_validated_community_inputs() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        assert_eq!(
            EvidencePolicy::from_config(&cfg).map(|policy| policy.mode),
            Some(EvidenceMode::Single)
        );

        cfg.is_paired = true;
        assert_eq!(
            EvidencePolicy::from_config(&cfg).map(|policy| policy.mode),
            Some(EvidenceMode::Paired)
        );

        cfg.strain = true;
        assert!(EvidencePolicy::from_config(&cfg).is_none());

        cfg.strain = false;
        cfg.is_paired = false;
        cfg.sequencer = Sequencer::Nanopore;
        assert_eq!(
            EvidencePolicy::from_config(&cfg).map(|policy| policy.mode),
            Some(EvidenceMode::Long)
        );

        cfg.is_paired = true;
        assert!(EvidencePolicy::from_config(&cfg).is_none());
    }

    #[test]
    fn a_candidate_with_nearly_all_depth_takes_the_whole_group() {
        let totals = |depths: &[(u32, f64)]| {
            depths
                .iter()
                .map(|&(taxid, depth)| {
                    (
                        taxid,
                        TaxonTotals {
                            depth,
                            ..TaxonTotals::default()
                        },
                    )
                })
                .collect::<HashMap<_, _>>()
        };
        let whole = depth_shares(&[1, 2], &totals(&[(1, 1.0), (2, 99.0)]));
        assert!(matches!(whole, Some(Shares::Whole(1))));
        let Some(Shares::Split(shares)) = depth_shares(&[1, 2], &totals(&[(1, 1.0), (2, 3.0)]))
        else {
            panic!("expected a split");
        };
        assert_eq!(shares, vec![0.25, 0.75]);
        assert!(depth_shares(&[1, 2], &totals(&[])).is_none());
    }

    #[test]
    fn fraction_normalization_is_independent_of_hash_insertion_order() {
        let mut forward = HashMap::new();
        forward.insert(1, 1.0e16);
        forward.insert(2, 1.0);
        forward.insert(3, 1.0);

        let mut reverse = HashMap::new();
        reverse.insert(3, 1.0);
        reverse.insert(2, 1.0);
        reverse.insert(1, 1.0e16);

        normalize_map(&mut forward);
        normalize_map(&mut reverse);
        for taxid in [1, 2, 3] {
            assert_eq!(forward[&taxid].to_bits(), reverse[&taxid].to_bits());
        }
    }

    #[test]
    fn exact_tie_storage_follows_the_resolved_short_read_policy() {
        let rows = vec![
            native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let bundle = || ReadBundle {
            read_name: "tie".into(),
            alignments: rows.clone(),
        };
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = true;

        let default_policy = EvidencePolicy::from_config(&cfg).unwrap();
        let default_result = process_evidence_batch(
            &[bundle()],
            &cfg,
            default_policy,
            false,
            &HashMap::new(),
            None,
        );
        assert_eq!(default_result.exact_ties.fragment_count(), 0);

        cfg.exact_ties = ExactTiePolicy::GuardedMap;
        let guarded_policy = EvidencePolicy::from_config(&cfg).unwrap();
        assert!(guarded_policy.collect_exact_ties);
        let guarded_result = process_evidence_batch(
            &[bundle()],
            &cfg,
            guarded_policy,
            false,
            &HashMap::new(),
            None,
        );
        assert_eq!(guarded_result.exact_ties.fragment_count(), 1);
        assert!(guarded_result.presence.is_empty());
        for taxid in [101, 102] {
            assert_eq!(
                guarded_result.taxon_evidence[&taxid].resolved_candidate_fragments,
                0
            );
            assert_eq!(guarded_result.taxon_evidence[&taxid].best_fragments, 0);
            assert_eq!(
                guarded_result.taxon_evidence[&taxid].clear_best_fragments,
                0
            );
        }
    }

    #[test]
    fn single_end_policy_keeps_em_candidates_and_uses_half_error_for_coverage() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let bundle = ReadBundle {
            read_name: "single".into(),
            alignments: vec![
                native_mapping(101, 10, 1, 1.000, 100, 100, MappingRole::Unpaired),
                native_mapping(102, 20, 2, 0.995, 100, 100, MappingRole::Unpaired),
            ],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        assert_eq!(policy.mode, EvidenceMode::Single);
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new(), None);

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].taxids, vec![101, 102]);
        assert!(result.presence.contains_key(&PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        }));
        assert!(!result.presence.contains_key(&PresenceKey {
            taxid: 102,
            profile_rank: ProfileRank::Species,
            assembly_id: 2,
        }));
        assert_eq!(result.taxon_evidence[&101].best_fragments, 1);
        assert_eq!(result.taxon_evidence[&101].clear_best_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].clear_best_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].resolved_candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].resolved_candidate_fragments, 1);
    }

    #[test]
    fn candidate_boundary_singleton_stays_in_em_without_unique_breadth() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        let mut mapping = native_mapping(101, 10, 1, 1.0, 100, 100, MappingRole::Unpaired);
        mapping.candidate_boundary = metam::mapper::CandidateBoundary::from_bits(1).unwrap();
        let bundle = ReadBundle {
            read_name: "candidate-boundary".into(),
            alignments: vec![mapping],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new(), None);

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].taxids, vec![101]);
        assert_eq!(result.boundary_ambiguous_fragments, 1);
        assert!(result.uniquely_supported_taxa.is_empty());
        assert!(result.presence.is_empty());
        assert_eq!(result.taxon_evidence[&101].candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&101].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].clear_best_fragments, 0);
    }

    #[test]
    fn mixed_rank_exact_tie_stays_in_em_without_manufacturing_breadth() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.exact_ties = ExactTiePolicy::GuardedMap;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let bundle = ReadBundle {
            read_name: "mixed-rank".into(),
            alignments: vec![
                native_mapping_at_rank(
                    111,
                    ProfileRank::Genus,
                    10,
                    1,
                    0.99,
                    100,
                    100,
                    MappingRole::Unpaired,
                ),
                native_mapping_at_rank(
                    1111,
                    ProfileRank::Species,
                    20,
                    2,
                    0.99,
                    100,
                    100,
                    MappingRole::Unpaired,
                ),
            ],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new(), None);

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].rank, AssignRank::Mixed);
        assert_eq!(result.read_assignments[0].taxids, vec![111, 1111]);
        assert!(result.presence.is_empty());
        assert_eq!(result.exact_ties.fragment_count(), 0);
    }

    #[test]
    fn long_read_mixed_rank_tie_stays_in_em_without_manufacturing_breadth() {
        let mut cfg = make_cfg(Sequencer::Nanopore, None);
        cfg.identity = 0.92;
        cfg.mapped_len = 500;
        cfg.fraction = 0.20;
        let bundle = ReadBundle {
            read_name: "mixed-rank-long".into(),
            alignments: vec![
                native_mapping_at_rank(
                    111,
                    ProfileRank::Genus,
                    10,
                    1,
                    0.98,
                    1_000,
                    2_000,
                    MappingRole::Unpaired,
                ),
                native_mapping_at_rank(
                    1111,
                    ProfileRank::Species,
                    20,
                    2,
                    0.98,
                    1_000,
                    2_000,
                    MappingRole::Unpaired,
                ),
            ],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        assert_eq!(policy.mode, EvidenceMode::Long);
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new(), None);

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].rank, AssignRank::Mixed);
        assert_eq!(result.read_assignments[0].taxids, vec![111, 1111]);
        assert!(result.presence.is_empty());
        assert_eq!(result.exact_ties.fragment_count(), 0);
    }

    #[test]
    fn direct_genus_reference_keeps_genus_rank_and_coordinates() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let bundle = ReadBundle {
            read_name: "genus-only".into(),
            alignments: vec![native_mapping_at_rank(
                111,
                ProfileRank::Genus,
                10,
                1,
                0.99,
                100,
                100,
                MappingRole::Unpaired,
            )],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new(), None);

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].rank, AssignRank::Genus);
        assert!(result.presence.contains_key(&PresenceKey {
            taxid: 111,
            profile_rank: ProfileRank::Genus,
            assembly_id: 1,
        }));
    }

    #[test]
    fn guarded_map_collects_single_end_exact_ties_without_premature_coverage() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.exact_ties = ExactTiePolicy::GuardedMap;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let bundle = ReadBundle {
            read_name: "single-tie".into(),
            alignments: vec![
                native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Unpaired),
                native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Unpaired),
            ],
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new(), None);

        assert_eq!(result.exact_ties.fragment_count(), 1);
        assert!(result.presence.is_empty());
        assert_eq!(result.read_assignments[0].taxids, vec![101, 102]);
    }

    #[test]
    fn seeded_map_prior_uses_unique_depth_and_fractional_ambiguous_em() {
        let assignments = vec![
            ReadAssignment {
                read: "unique-1".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![101],
                depths: vec![95.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
                best: 0,
            },
            ReadAssignment {
                read: "unique-2".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![102],
                depths: vec![5.0],
                mapped_lens: vec![100.0],
                sampling_scales: vec![1.0],
                best: 0,
            },
            ReadAssignment {
                read: "ambiguous".into(),
                rank: AssignRank::Species,
                requires_terminal_evidence: true,
                taxids: vec![101, 102, 103],
                depths: vec![100.0, 100.0, 100.0],
                mapped_lens: vec![100.0; 3],
                sampling_scales: vec![1.0; 3],
                best: 0,
            },
        ];
        let eligible = HashSet::from([101, 102]);
        let (priors, iterations) = seeded_short_read_map_priors(&assignments, &eligible, 200);
        assert!(iterations > 0);
        assert!(!priors.contains_key(&103));
        let total = priors[&101] + priors[&102];
        assert!((priors[&101] / total - 0.95).abs() < 1e-12);
        assert!((priors[&102] / total - 0.05).abs() < 1e-12);
    }

    #[test]
    fn paired_exact_tie_stays_in_em_but_creates_no_coordinate_evidence() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = true;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let rows = vec![
            native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(101, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(102, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let bundle = ReadBundle {
            read_name: "tie".into(),
            alignments: rows,
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new(), None);

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].taxids, vec![101, 102]);
        assert_eq!(result.read_assignments[0].mapped_lens, vec![100.0, 100.0]);
        assert!(result.presence.is_empty());
        assert_eq!(result.taxon_evidence[&101].candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&101].resolved_candidate_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].resolved_candidate_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].clear_best_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].clear_best_fragments, 0);
        assert_eq!(result.exact_ties.fragment_count(), 0);
    }

    #[test]
    fn automatic_long_read_identity_follows_reads_that_align_worse() {
        let unpaired = MappingRole::Unpaired;
        let reads: Vec<ReadBundle> = [0.80_f32, 0.84, 0.88, 0.89]
            .iter()
            .enumerate()
            .map(|(i, &identity)| ReadBundle {
                read_name: format!("r{i}").into(),
                alignments: vec![
                    native_mapping(101, 10, 1, identity, 1_000, 2_000, unpaired),
                    // Higher identity over fewer bases is not the read's best record.
                    native_mapping(102, 20, 2, identity + 0.05, 600, 2_000, unpaired),
                ],
            })
            .collect();
        let resolved = |claimed: Option<f64>, identity: f64| {
            let mut cfg = make_cfg(Sequencer::Nanopore, None);
            cfg.mapped_len = 500;
            cfg.fraction = 0.20;
            cfg.identity = identity;
            cfg.identity_check_accuracy = claimed;
            let taxonomy = Taxonomy::from_parts(
                HashMap::new(),
                HashMap::new(),
                HashMap::new(),
                HashSet::new(),
            );
            let mut profiler = CommunityProfiler::new(cfg, Arc::new(taxonomy)).unwrap();
            profiler.lower_identity_to_alignments(&reads);
            profiler.cfg.identity
        };
        // Upper quartile 0.89 below the claimed 0.94: cutoff 0.89 - 0.03.
        assert!((resolved(Some(0.94), 0.91) - 0.86).abs() < 1e-6);
        // Reads align as claimed, the cutoff would rise, or it is not automatic.
        assert_eq!(resolved(Some(0.88), 0.85), 0.85);
        assert_eq!(resolved(Some(0.99), 0.80), 0.80);
        assert_eq!(resolved(None, 0.91), 0.91);
    }

    #[test]
    fn paired_margin_resolution_uses_mate1_for_abundance_and_coverage() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = true;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let rows = vec![
            native_mapping(101, 10, 1, 1.00, 90, 100, MappingRole::Mate1),
            native_mapping(101, 10, 1, 1.00, 100, 100, MappingRole::Mate2),
            native_mapping(102, 20, 2, 0.98, 100, 100, MappingRole::Mate1),
            native_mapping(102, 20, 2, 0.98, 100, 100, MappingRole::Mate2),
        ];
        let bundle = ReadBundle {
            read_name: "resolved".into(),
            alignments: rows,
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new(), None);

        assert_eq!(result.read_assignments.len(), 1);
        assert_eq!(result.read_assignments[0].taxids, vec![101, 102]);
        assert_eq!(result.read_assignments[0].mapped_lens, vec![90.0, 100.0]);
        let presence = &result.presence[&PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        }];
        assert_eq!(presence.count, 1);
        assert_eq!(presence.reference_span_bases, 90);
        assert_eq!(presence.intervals.len(), 1);
        assert_eq!(result.taxon_evidence[&101].best_fragments, 1);
        assert_eq!(result.taxon_evidence[&101].clear_best_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].best_fragments, 0);
        assert_eq!(result.taxon_evidence[&102].clear_best_fragments, 0);
        assert_eq!(result.taxon_evidence[&101].resolved_candidate_fragments, 1);
        assert_eq!(result.taxon_evidence[&102].resolved_candidate_fragments, 1);
        assert!(!result.presence.contains_key(&PresenceKey {
            taxid: 102,
            profile_rank: ProfileRank::Species,
            assembly_id: 2,
        }));
    }

    #[test]
    fn resolved_species_keeps_assemblies_as_independent_witnesses() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = true;
        cfg.identity = 0.97;
        cfg.mapped_len = 50;
        cfg.fraction = 0.70;
        let rows = vec![
            native_mapping(101, 10, 1, 1.00, 100, 100, MappingRole::Mate1),
            native_mapping(101, 10, 1, 1.00, 100, 100, MappingRole::Mate2),
            native_mapping(101, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            native_mapping(101, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let bundle = ReadBundle {
            read_name: "two-assemblies".into(),
            alignments: rows,
        };
        let policy = EvidencePolicy::from_config(&cfg).unwrap();
        let result = process_evidence_batch(&[bundle], &cfg, policy, false, &HashMap::new(), None);

        assert_eq!(result.presence.len(), 2);
        assert!(result.presence.contains_key(&PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 1,
        }));
        assert!(result.presence.contains_key(&PresenceKey {
            taxid: 101,
            profile_rank: ProfileRank::Species,
            assembly_id: 2,
        }));
        // The fragment fits genome 1 best.
        let fits = |assembly_id| {
            result.presence[&PresenceKey {
                taxid: 101,
                profile_rank: ProfileRank::Species,
                assembly_id,
            }]
                .fits
        };
        assert_eq!((fits(1), fits(2)), (1, 0));
    }

    #[test]
    fn abundance_uses_the_best_fit_genome_of_each_present_taxon() {
        let mut assignments = vec![ReadAssignment {
            read: Arc::from("r"),
            rank: AssignRank::Mixed,
            requires_terminal_evidence: true,
            taxids: vec![1, 2],
            depths: vec![150.0 / 1_000.0, 150.0 / 1_000.0],
            mapped_lens: vec![150.0, 150.0],
            sampling_scales: vec![1.0, 1.0],
            best: 0,
        }];
        let best_fit = HashMap::from([(
            1,
            BestFitGenome {
                assembly_id: 7,
                asm_len: 1_500.0,
                sampling_scale: 2.0,
            },
        )]);
        rescale_to_best_fit(&mut assignments, &best_fit, None);
        // Taxon 1 is measured on its 1,500-base genome; taxon 2, not present,
        // keeps the genome of its own placement.
        assert_eq!(assignments[0].depths, vec![0.1, 0.15]);
        assert_eq!(assignments[0].sampling_scales, vec![2.0, 1.0]);
    }

    #[test]
    fn a_present_taxon_is_reported_on_its_best_fit_passing_genome() {
        let key = |taxid, assembly_id| PresenceKey {
            taxid,
            profile_rank: ProfileRank::Species,
            assembly_id,
        };
        // Depth tags the genome a decision was measured on.
        let witness = |assembly_id: u32, passes: bool, fixed_chunk: f64| WitnessDecision {
            metrics: EvidenceCoverage {
                depth: f64::from(assembly_id),
                fixed_chunk,
                ..EvidenceCoverage::default()
            },
            standard: passes,
            rescues: RescueKinds::default(),
            required_premap_guard: PremapChunkGuard::None,
        };
        let fit = |fits, count| WitnessGenome {
            fits,
            count,
            ..WitnessGenome::default()
        };
        let support = GenomeSupport::from([
            (key(1, 5), fit(3, 9)),
            (key(1, 6), fit(4, 4)),
            (key(1, 7), fit(9, 9)),
            (key(2, 8), fit(2, 5)),
            (key(2, 9), fit(2, 6)),
            (key(3, 1), fit(1, 1)),
            (key(3, 2), fit(1, 1)),
            (key(4, 4), fit(1, 1)),
            (key(4, 3), fit(1, 1)),
        ]);
        let mut decisions = HashMap::new();
        report_best_fit(
            &mut decisions,
            [
                // Most fits among the passing genomes; genome 7 fits more
                // fragments best but fails.
                (key(1, 5), witness(5, true, 0.9)),
                (key(1, 6), witness(6, true, 0.1)),
                (key(1, 7), witness(7, false, 0.9)),
                // Equal fits: more fragments placed.
                (key(2, 8), witness(8, true, 0.9)),
                (key(2, 9), witness(9, true, 0.1)),
                // Equal fits and fragments: the larger fixed chunk breadth.
                (key(3, 1), witness(1, true, 0.2)),
                (key(3, 2), witness(2, true, 0.3)),
                // All equal: the lower assembly id.
                (key(4, 4), witness(4, true, 0.5)),
                (key(4, 3), witness(3, true, 0.5)),
                // No passing genome: nothing to report.
                (key(5, 1), witness(1, false, 0.9)),
            ],
            &support,
        );
        let genome = |taxid: u32| decisions[&taxid].metrics.depth;
        assert_eq!(decisions.len(), 4);
        assert_eq!(
            (genome(1), genome(2), genome(3), genome(4)),
            (6.0, 9.0, 2.0, 3.0)
        );
    }

    /// Intermediate synthetic alignment used inside tests.
    ///
    /// `AlignmentEntry` stores the pre-parsed
    /// `asm`/`acc`/`species_taxid` rather than the raw `ref_name`. Tests
    /// describe alignments by reference name and (identity, mapped fraction),
    /// and let [`run_single_with_normalize`] derive the rest.
    #[derive(Clone, Debug)]
    struct RawAln {
        ref_name: String,
        quality: (f64, f64),
        mapped_len: usize,
        span: (usize, usize),
    }

    /// Construct a synthetic alignment descriptor suitable for
    /// [`run_single`] / [`run_single_with_normalize`].
    ///
    /// `ref_name` must match the expected `asm|taxid|species_taxid|acc|...`
    /// layout so the profiler can extract the species taxid.
    fn aln_entry(
        ref_name: &str,
        quality: (f64, f64),
        mapped_len: usize,
        span: (usize, usize),
    ) -> RawAln {
        RawAln {
            ref_name: ref_name.to_string(),
            quality,
            mapped_len,
            span,
        }
    }

    /// Convert a [`RawAln`] (reference name `assembly|…|taxid|accession`)
    /// into an [`AlignmentEntry`] with the derived fields the profiler uses.
    /// Returns `None` when the synthetic alignment has no parsable species
    /// taxid or no known assembly length, so tests can cover rejected
    /// records.
    fn build_aln(
        raw: &RawAln,
        asm_lengths: &HashMap<String, f64>,
        species_normalize: &HashMap<u32, u32>,
    ) -> Option<AlignmentEntry> {
        let mut parts = raw.ref_name.split('|');
        let asm = parts.next()?;
        if asm.is_empty() {
            return None;
        }
        let _ = parts.next();
        let raw_taxid = parts.next()?.parse::<u32>().ok()?;
        let species_taxid = species_normalize
            .get(&raw_taxid)
            .copied()
            .unwrap_or(raw_taxid);
        let acc = parts.next().unwrap_or("");
        let asm_len = asm_lengths.get(asm).copied().unwrap_or(0.0);
        if asm_len <= 0.0 {
            return None;
        }
        let (identity, fraction) = raw.quality;
        Some(AlignmentEntry {
            asm: Arc::<str>::from(asm),
            acc: Arc::<str>::from(acc),
            profile_taxid: species_taxid,
            profile_rank: ProfileRank::Species,
            reference_id: u32::MAX,
            assembly_id: u32::MAX,
            role: MappingRole::Unpaired,
            pair_id: 0,
            candidate_boundary: metam::mapper::CandidateBoundary::COMPLETE,
            sampling_scale: 1.0,
            asm_len,
            identity: identity as f32,
            fraction: fraction as f32,
            mapped_len: clamp_u32(raw.mapped_len),
            read_len: clamp_u32(raw.mapped_len),
            span: (clamp_u32(raw.span.0), clamp_u32(raw.span.1)),
            fixed_chunk: (clamp_u32(raw.span.0), clamp_u32(raw.span.1)),
            flex_chunk: (clamp_u32(raw.span.0), clamp_u32(raw.span.1)),
            marker_gene_len: 0,
        })
    }

    // ------------------------------------------------------------------
    // Fallback gate tests
    // ------------------------------------------------------------------

    #[test]
    fn fallback_gate_activates_only_below_mapping_rate_threshold() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        assert!(fallback_gate_active(&cfg, 0.0));
        assert!(fallback_gate_active(&cfg, 0.29));
        assert!(fallback_gate_active(&cfg, 0.2999999));
        // Exactly at the 0.30 threshold is NOT active (strict `<`).
        assert!(!fallback_gate_active(&cfg, 0.30));
        assert!(!fallback_gate_active(&cfg, 0.50));
        assert!(!fallback_gate_active(&cfg, 1.0));
    }

    #[test]
    fn fallback_gate_disabled_when_host_is_set() {
        // Any host string (e.g. human 9606) disables fallback regardless
        // of sequencer or mapping rate.
        let cfg_host = make_cfg(Sequencer::Illumina, Some("9606"));
        assert!(!fallback_gate_active(&cfg_host, 0.0));
        assert!(!fallback_gate_active(&cfg_host, 0.1));
    }

    #[test]
    fn fallback_gate_covers_long_reads_but_not_assembly() {
        // Nanopore and PacBio honour the low-map auto-trigger, so a low
        // mapping rate enables genus fallback for them too.
        for seq in [Sequencer::Nanopore, Sequencer::PacBio] {
            let cfg = make_cfg(seq, None);
            assert!(
                fallback_gate_active(&cfg, 0.0),
                "expected fallback enabled for low-map {seq}"
            );
            assert!(
                !fallback_gate_active(&cfg, 0.99),
                "expected fallback disabled for high-map {seq}"
            );
        }
    }

    #[test]
    fn fallback_gate_respects_custom_low_map_rate_threshold() {
        // With the default 0.30 threshold, 0.40 is NOT low-map.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        assert!(!fallback_gate_active(&cfg, 0.40));
        // Raise the threshold: now 0.40 IS low-map (strict `<`).
        cfg.low_map_rate_threshold = 0.50;
        assert!(fallback_gate_active(&cfg, 0.40));
        assert!(fallback_gate_active(&cfg, 0.4999));
        assert!(!fallback_gate_active(&cfg, 0.50));
        // Lower the threshold: even 0.05 no longer triggers when = threshold.
        cfg.low_map_rate_threshold = 0.05;
        assert!(!fallback_gate_active(&cfg, 0.05));
        assert!(fallback_gate_active(&cfg, 0.0499));
    }

    #[test]
    fn resolve_by_aligned_honors_explicit_flag() {
        // --by-mapped flag present: mapped basis is forced regardless of
        // sequencer, host, or mapping rate.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.by_aligned = true;
        assert!(cfg.resolve_by_aligned(0.99));
        cfg.sequencer = Sequencer::Nanopore;
        assert!(cfg.resolve_by_aligned(0.99));
        cfg.host = Some("9606".to_string());
        assert!(cfg.resolve_by_aligned(0.99));
    }

    fn pvalue_filter(rank: &str, cov: Option<f64>, chunk: Option<f64>, strict: bool) -> bool {
        pvalue_filter_passes(rank, cov, chunk, strict)
    }

    #[test]
    fn strict_filter_requires_both_probs_above_cutoff() {
        // Both probabilities at or above the 1e-5 cutoff => keep in both modes.
        assert!(pvalue_filter("species", Some(1e-5), Some(1e-3), false));
        assert!(pvalue_filter("species", Some(1e-5), Some(1e-3), true));

        // Exactly one below cutoff:
        //   non-strict (OR): kept
        //   strict    (AND): dropped
        assert!(pvalue_filter("species", Some(1e-6), Some(1e-3), false));
        assert!(!pvalue_filter("species", Some(1e-6), Some(1e-3), true));
        assert!(pvalue_filter("species", Some(1e-3), Some(1e-6), false));
        assert!(!pvalue_filter("species", Some(1e-3), Some(1e-6), true));

        // Both below cutoff => dropped in both modes.
        assert!(!pvalue_filter("species", Some(1e-6), Some(1e-7), false));
        assert!(!pvalue_filter("species", Some(1e-6), Some(1e-7), true));
    }

    #[test]
    fn strict_filter_handles_missing_probabilities() {
        // Non-strict mode keeps rows with missing probabilities.
        assert!(pvalue_filter("species", None, Some(1e-3), false));
        assert!(pvalue_filter("species", Some(1e-3), None, false));
        assert!(pvalue_filter("species", None, None, false));

        // Species rows have p-values, so strict mode drops rows with missing
        // probability evidence.
        assert!(!pvalue_filter("species", None, Some(1e-3), true));
        assert!(!pvalue_filter("species", Some(1e-3), None, true));
        assert!(!pvalue_filter("species", None, None, true));

        // Legacy probabilities are computed for species only, so every higher
        // rank is kept under strict when both are absent -- not just genus.
        for rank in ["genus", "family", "order", "class", "phylum"] {
            assert!(
                pvalue_filter(rank, None, None, true),
                "{rank} rows must not be rejected for a test that was never run"
            );
            assert!(!pvalue_filter(rank, None, Some(1e-3), true));
            assert!(!pvalue_filter(rank, Some(1e-3), None, true));
        }
    }

    #[test]
    fn resolve_by_aligned_auto_forces_for_low_map_read_platforms_non_pathogen() {
        // Default cfg has --by-mapped unset (false); auto-rule applies.
        let cfg = make_cfg(Sequencer::Illumina, None);
        assert!(!cfg.by_aligned);
        // Below default 0.30 threshold => aligned basis is auto-forced.
        assert!(cfg.resolve_by_aligned(0.10));
        // At or above threshold => total-read basis.
        assert!(!cfg.resolve_by_aligned(0.30));
        assert!(!cfg.resolve_by_aligned(0.50));
        // Pathogen mode on disables auto-force.
        let cfg_host = make_cfg(Sequencer::Illumina, Some("9606"));
        assert!(!cfg_host.resolve_by_aligned(0.10));
        // Long-read platforms also auto-force on low map.
        let cfg_np = make_cfg(Sequencer::Nanopore, None);
        assert!(cfg_np.resolve_by_aligned(0.10));
    }

    #[test]
    fn fallback_scale_chunk_basis_applies_only_when_scale_true() {
        // Pass-through when scaling is off.
        assert_eq!(fallback_scale_chunk_basis(0, false), 0);
        assert_eq!(fallback_scale_chunk_basis(1, false), 1);
        assert_eq!(fallback_scale_chunk_basis(100_000, false), 100_000);

        // When scaling is on, count is multiplied by FALLBACK_CHUNK_BASIS_SCALE
        // and rounded to the nearest integer.
        assert_eq!(fallback_scale_chunk_basis(0, true), 0);
        assert_eq!(fallback_scale_chunk_basis(1, true), 2); // 1.5 -> 2
        assert_eq!(fallback_scale_chunk_basis(2, true), 3);
        assert_eq!(fallback_scale_chunk_basis(10, true), 15);
        assert_eq!(fallback_scale_chunk_basis(100, true), 150);
        assert_eq!(fallback_scale_chunk_basis(100_000, true), 150_000);
        // Round-half-to-even-ish: 3 * 1.5 = 4.5 => rounds to 5.
        assert_eq!(fallback_scale_chunk_basis(3, true), 5);
    }

    #[test]
    fn fallback_chunk_basis_scale_constant_is_one_and_a_half() {
        assert!((FALLBACK_CHUNK_BASIS_SCALE - 1.5).abs() < 1e-12);
    }

    #[test]
    fn chunk_breadth_plan_scales_only_low_map_aligned_basis() {
        let total_reads = 100_000usize;
        let aligned_reads = 10_000usize;

        // --by-mapped alone: mapped basis, same equation, no 1.5x.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.by_aligned = true;
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // Automatic low-map fallback on a full DB: aligned basis with 1.5x.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        let gate = fallback_gate_active(&cfg, 0.10);
        assert!(genus_fallback_mode_active(&cfg, gate, false));
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, false, gate);
        assert!(plan.use_aligned_basis);
        assert!(plan.by_aligned_auto);
        assert!(plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, 15_000);

        // Explicit --genus-fallback on non-low-map data without --by-mapped:
        // total basis, no scale.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.genus_fallback = true;
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(!plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, total_reads);
        assert_eq!(plan.estimate_reads, total_reads);

        // Explicit --genus-fallback plus --by-mapped on non-low-map data:
        // mapped basis, same equation, no scale.
        cfg.by_aligned = true;
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // Explicit --genus-fallback on a subsampled DB still enables fallback
        // assignment, but chunk-breadth estimation uses total reads.
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, true, false);
        assert!(genus_fallback_mode_active(&cfg, false, true));
        assert!(!plan.use_aligned_basis);
        assert!(plan.force_total_for_subsampled_fallback);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, total_reads);
        assert_eq!(plan.estimate_reads, total_reads);

        // A low-map subsampled run without explicit --genus-fallback keeps the
        // automatic fallback assignment disabled, so it must not scale.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        let gate = fallback_gate_active(&cfg, 0.10);
        assert!(!genus_fallback_mode_active(&cfg, gate, true));
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, true, gate);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // User-pinned --min-chunk-breadth is never rescaled.
        cfg.chunk_breadth = Some(0.25);
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, false, gate);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.estimate_reads, aligned_reads);
    }

    #[test]
    fn fallback_gate_active_for_illumina_single_and_paired() {
        // Single-end
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = false;
        assert!(fallback_gate_active(&cfg, 0.1));
        // Paired-end
        cfg.is_paired = true;
        assert!(fallback_gate_active(&cfg, 0.1));
    }

    // ------------------------------------------------------------------
    // Legacy non-LCA rank filter bounds
    // ------------------------------------------------------------------

    #[test]
    fn legacy_rank_filter_lower_bounds_remain_unchanged() {
        assert_eq!(rank_filter_lower_bound("species"), 0.75);
        assert_eq!(rank_filter_lower_bound("genus"), 0.65);
        // Unknown ranks default to species-level stringency.
        assert_eq!(rank_filter_lower_bound("family"), 0.75);
    }

    #[test]
    fn community_chunk_breadth_curve_matches_spec() {
        // Short (Illumina): 0.66*m/(m+1) + 0.24*m/(m+100), plateau 0.90.
        let short = |n| community_chunk_breadth_cutoff(n, true);
        assert!((short(0) - 0.0).abs() < 1e-12);
        assert!((short(10_000) - 0.00656).abs() < 1e-4);
        assert!((short(100_000) - 0.06024).abs() < 1e-4);
        assert!((short(1_000_000) - 0.33238).abs() < 1e-4); // C = 0.33 at 1M reads
        assert!((short(3_300_000) - 0.51418).abs() < 1e-4);
        assert!((short(16_600_000) - 0.65667).abs() < 1e-4);
        // approaches the 0.90 plateau at extreme depth, always below the clamp
        let short_deep = short(100_000_000_000);
        assert!(short_deep > 0.899 && short_deep < 0.95);

        // Long (Nanopore/PacBio): 0.85*m/(m+0.78), plateau 0.85.
        let long = |n| community_chunk_breadth_cutoff(n, false);
        assert!((long(0) - 0.0).abs() < 1e-12);
        assert!((long(10_000) - 0.01076).abs() < 1e-4);
        assert!((long(100_000) - 0.09659).abs() < 1e-4);
        assert!((long(1_000_000) - 0.47753).abs() < 1e-4);
        assert!((long(3_300_000) - 0.68750).abs() < 1e-4);
        assert!((long(16_600_000) - 0.81185).abs() < 1e-4);
        let long_deep = long(100_000_000_000);
        assert!(long_deep > 0.849 && long_deep < 0.85);

        // Monotonic non-decreasing in read count on both platforms.
        for &(a, b) in &[
            (1_000usize, 10_000usize),
            (100_000, 1_000_000),
            (1_000_000, 16_600_000),
        ] {
            assert!(
                community_chunk_breadth_cutoff(b, true) >= community_chunk_breadth_cutoff(a, true)
            );
            assert!(
                community_chunk_breadth_cutoff(b, false)
                    >= community_chunk_breadth_cutoff(a, false)
            );
        }
    }

    #[test]
    fn recall_scales_only_deep_automatic_community_chunk_cutoffs() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.mode = ProfileMode::Recall;

        let at_gate = resolve_community_chunk_breadth(&cfg, 1_000_000);
        approx_equal(at_gate.recall_scale, 1.0);
        approx_equal(at_gate.standard, at_gate.safeguard);
        assert_eq!(at_gate.standard.to_bits(), at_gate.safeguard.to_bits());

        let two_million = resolve_community_chunk_breadth(&cfg, 2_000_000);
        approx_equal(two_million.recall_scale, 0.85);
        approx_equal(
            two_million.standard,
            two_million.safeguard * two_million.recall_scale,
        );

        let five_million = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(five_million.recall_scale, 0.76);
        let deep = resolve_community_chunk_breadth(&cfg, 16_600_000);
        approx_equal(deep.recall_scale, 0.70 + 0.30 / 16.6);
        assert!(deep.recall_scale > SHORT_RECALL_CHUNK_SCALE_FLOOR);

        cfg.sequencer = Sequencer::Nanopore;
        let long_at_gate = resolve_community_chunk_breadth(&cfg, 1_000_000);
        approx_equal(long_at_gate.recall_scale, 1.0);
        let long_deep = resolve_community_chunk_breadth(&cfg, 1_000_001);
        approx_equal(long_deep.recall_scale, LONG_RECALL_CHUNK_SCALE);
        approx_equal(
            long_deep.standard,
            long_deep.safeguard * LONG_RECALL_CHUNK_SCALE,
        );

        cfg.sequencer = Sequencer::PacBio;
        approx_equal(
            resolve_community_chunk_breadth(&cfg, 5_000_000).recall_scale,
            LONG_RECALL_CHUNK_SCALE,
        );
    }

    #[test]
    fn recall_chunk_scaling_respects_user_and_nonstandard_policy_boundaries() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.mode = ProfileMode::Recall;
        cfg.chunk_breadth = Some(0.37);
        let explicit = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(explicit.standard, 0.37);
        approx_equal(explicit.safeguard, 0.37);
        approx_equal(explicit.recall_scale, 1.0);
        assert_eq!(explicit.standard.to_bits(), explicit.safeguard.to_bits());

        cfg.chunk_breadth = None;
        cfg.lowbiomass = true;
        let low_biomass = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(low_biomass.standard, 0.0);
        approx_equal(low_biomass.safeguard, 0.0);
        approx_equal(low_biomass.recall_scale, 1.0);

        // A read minimum adds to the automatic cutoff and leaves it unchanged.
        cfg.lowbiomass = false;
        let without_min_reads = resolve_community_chunk_breadth(&cfg, 20_000_000);
        cfg.min_reads = Some(3);
        assert_eq!(
            resolve_community_chunk_breadth(&cfg, 20_000_000),
            without_min_reads
        );

        cfg.min_reads = None;
        cfg.mode = ProfileMode::Default;
        let default = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(default.recall_scale, 1.0);
        assert_eq!(default.standard.to_bits(), default.safeguard.to_bits());
        cfg.mode = ProfileMode::Precision;
        let precision = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(precision.recall_scale, 1.0);
        assert_eq!(precision.standard.to_bits(), precision.safeguard.to_bits());

        cfg.mode = ProfileMode::Recall;
        cfg.strain = true;
        approx_equal(
            resolve_community_chunk_breadth(&cfg, 20_000_000).recall_scale,
            1.0,
        );
    }

    #[test]
    fn user_chunk_scale_adjusts_automatic_standard_and_safeguard_cutoffs() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;

        let unchanged = resolve_community_chunk_breadth(&cfg, 5_000_000);
        assert_eq!(unchanged.automatic.to_bits(), unchanged.safeguard.to_bits());
        assert_eq!(unchanged.automatic.to_bits(), unchanged.standard.to_bits());

        cfg.chunk_breadth_scale = 0.5;
        let relaxed = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(relaxed.safeguard, relaxed.automatic * 0.5);
        assert_eq!(relaxed.standard.to_bits(), relaxed.safeguard.to_bits());

        cfg.chunk_breadth_scale = 10.0;
        let capped = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(capped.safeguard, 1.0);
        approx_equal(capped.standard, 1.0);
    }

    #[test]
    fn recall_scaling_composes_after_user_chunk_scaling() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.chunk_breadth_scale = 0.5;
        cfg.mode = ProfileMode::Recall;

        let thresholds = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(thresholds.user_scale, 0.5);
        approx_equal(thresholds.recall_scale, 0.76);
        approx_equal(thresholds.safeguard, thresholds.automatic * 0.5);
        approx_equal(
            thresholds.standard,
            thresholds.safeguard * thresholds.recall_scale,
        );
    }

    #[test]
    fn ratio_filter_genus_passes_at_066_species_fails_at_066() {
        // Both ratios = 0.66 => genus passes (>= 0.65), species fails (< 0.75).
        let br_obs = 0.66;
        let br_exp = 1.0;
        let fx_obs = 0.66;
        let fx_exp = 1.0;
        assert!(ratio_filter_passes(
            "genus", br_obs, br_exp, fx_obs, fx_exp, None, None
        ));
        assert!(!ratio_filter_passes(
            "species", br_obs, br_exp, fx_obs, fx_exp, None, None
        ));
    }

    #[test]
    fn coherent_lca_witness_is_single_strongest_and_order_independent() {
        // Weak-but-standard beats strong-but-rescued, matching witness_is_better.
        let decisions = HashMap::from([
            (10u32, witness(0.90, false, true, PremapChunkGuard::Half)),
            (20u32, witness(0.10, true, false, PremapChunkGuard::Half)),
        ]);
        let (taxid, chosen) = coherent_lca_witness(&[10, 20], &decisions).unwrap();
        assert_eq!(taxid, 20, "a standard witness outranks a rescued one");
        assert!(chosen.standard);
        // Candidate order must not matter.
        assert_eq!(
            coherent_lca_witness(&[20, 10], &decisions).unwrap().0,
            taxid
        );

        // Fully equivalent witnesses resolve to the lowest taxid either way.
        let tied = HashMap::from([
            (77u32, witness(0.50, true, false, PremapChunkGuard::Half)),
            (33u32, witness(0.50, true, false, PremapChunkGuard::Half)),
        ]);
        assert_eq!(coherent_lca_witness(&[77, 33], &tied).unwrap().0, 33);
        assert_eq!(coherent_lca_witness(&[33, 77], &tied).unwrap().0, 33);

        // The whole coverage set comes from that one witness, never a mixture.
        let mixed = HashMap::from([
            (
                1u32,
                WitnessDecision {
                    metrics: EvidenceCoverage {
                        fixed_chunk: 0.80,
                        breadth: 0.10,
                        flex_chunk: 0.10,
                        ..EvidenceCoverage::default()
                    },
                    standard: true,
                    rescues: RescueKinds::default(),
                    required_premap_guard: PremapChunkGuard::Half,
                },
            ),
            (
                2u32,
                WitnessDecision {
                    metrics: EvidenceCoverage {
                        fixed_chunk: 0.20,
                        breadth: 0.99,
                        flex_chunk: 0.99,
                        ..EvidenceCoverage::default()
                    },
                    standard: true,
                    rescues: RescueKinds::default(),
                    required_premap_guard: PremapChunkGuard::Half,
                },
            ),
        ]);
        let (_, chosen) = coherent_lca_witness(&[1, 2], &mixed).unwrap();
        assert_eq!(chosen.metrics.fixed_chunk, 0.80);
        assert_eq!(
            chosen.metrics.breadth, 0.10,
            "breadth must come from the chosen witness, not the per-field maximum"
        );
        assert_eq!(chosen.metrics.flex_chunk, 0.10);

        // No member with coordinate evidence -> no coherent witness.
        assert!(coherent_lca_witness(&[404, 405], &decisions).is_none());
        assert!(coherent_lca_witness(&[], &decisions).is_none());
    }

    #[test]
    fn repeated_lca_groups_retain_one_global_coherent_witness() {
        let weak = (20u32, witness(0.20, true, false, PremapChunkGuard::Half));
        let strong = (10u32, witness(0.80, true, false, PremapChunkGuard::Half));

        let mut forward = HashMap::new();
        retain_coherent_lca_witness(&mut forward, 100, weak);
        retain_coherent_lca_witness(&mut forward, 100, strong);

        let mut reverse = HashMap::new();
        retain_coherent_lca_witness(&mut reverse, 100, strong);
        retain_coherent_lca_witness(&mut reverse, 100, weak);

        assert_eq!(forward[&100].0, 10);
        assert_eq!(reverse[&100].0, 10);
        assert_eq!(forward[&100].1.metrics.fixed_chunk, 0.80);
        assert_eq!(reverse[&100].1.metrics.fixed_chunk, 0.80);

        // Equivalent witnesses use the source taxid as a deterministic final
        // tiebreak, including when they arrive from different LCA groups.
        let tied_high = (77u32, witness(0.50, true, false, PremapChunkGuard::Half));
        let tied_low = (33u32, witness(0.50, true, false, PremapChunkGuard::Half));
        retain_coherent_lca_witness(&mut forward, 200, tied_high);
        retain_coherent_lca_witness(&mut forward, 200, tied_low);
        assert_eq!(forward[&200].0, 33);
    }

    #[test]
    fn legacy_lca_coverage_is_merged_without_replacing_direct_coverage() {
        let direct = HashMap::from([(
            100u32,
            TaxonCoverage {
                breadth: 0.40,
                fixed_chunk: 0.50,
                flex_chunk: 0.60,
            },
        )]);
        let direct_taxa = HashSet::from([100u32]);
        let mut lca = HashMap::new();
        merge_legacy_lca_coverage(
            &mut lca,
            100,
            TaxonCoverage {
                breadth: 0.90,
                fixed_chunk: 0.10,
                flex_chunk: 0.20,
            },
        );
        merge_legacy_lca_coverage(
            &mut lca,
            100,
            TaxonCoverage {
                breadth: 0.30,
                fixed_chunk: 0.80,
                flex_chunk: 0.70,
            },
        );

        assert_eq!(lca[&100].breadth, 0.90);
        assert_eq!(lca[&100].fixed_chunk, 0.80);
        assert_eq!(lca[&100].flex_chunk, 0.70);

        let retained_direct = legacy_output_coverage(100, &direct_taxa, &direct, &lca);
        assert_eq!(retained_direct.breadth, 0.40);
        assert_eq!(retained_direct.fixed_chunk, 0.50);
        assert_eq!(retained_direct.flex_chunk, 0.60);

        let lca_only = legacy_output_coverage(100, &HashSet::new(), &direct, &lca);
        assert_eq!(lca_only.breadth, 0.90);
        assert_eq!(lca_only.fixed_chunk, 0.80);
        assert_eq!(lca_only.flex_chunk, 0.70);
    }

    #[test]
    fn majority_lca_uses_standard_rank_and_platform_ratio_defaults() {
        let short_species =
            standard_lca_ratio_lower_bounds("species", EvidenceMode::Paired).unwrap();
        assert_eq!(short_species, (0.70, 0.75));
        assert!(ratio_filter_passes_with_defaults(
            0.70,
            1.0,
            0.75,
            1.0,
            short_species,
            None,
            None,
        ));
        assert!(!ratio_filter_passes_with_defaults(
            0.70,
            1.0,
            0.74,
            1.0,
            short_species,
            None,
            None,
        ));

        let long_species = standard_lca_ratio_lower_bounds("species", EvidenceMode::Long).unwrap();
        assert_eq!(long_species, (0.70, 0.70));
        assert!(ratio_filter_passes_with_defaults(
            0.70,
            1.0,
            0.70,
            1.0,
            long_species,
            None,
            None,
        ));

        for rank in [
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "superkingdom",
        ] {
            let defaults = standard_lca_ratio_lower_bounds(rank, EvidenceMode::Single).unwrap();
            assert_eq!(defaults, (0.65, 0.65));
            assert!(ratio_filter_passes_with_defaults(
                0.65, 1.0, 0.65, 1.0, defaults, None, None,
            ));
        }

        assert_eq!(
            standard_lca_ratio_lower_bounds("no rank", EvidenceMode::Long),
            None
        );
    }

    #[test]
    fn majority_lca_ratio_defaults_still_honor_cli_overrides_and_upper_bound() {
        let defaults = standard_lca_ratio_lower_bounds("family", EvidenceMode::Paired).unwrap();
        assert!(!ratio_filter_passes_with_defaults(
            0.79,
            1.0,
            0.65,
            1.0,
            defaults,
            Some(0.80),
            None,
        ));
        assert!(!ratio_filter_passes_with_defaults(
            1.51, 1.0, 0.65, 1.0, defaults, None, None,
        ));
    }

    #[test]
    fn ratio_filter_both_fail_above_15() {
        // Ratio slightly above 1.5 => both ranks fail.
        let br_obs = 1.51;
        let br_exp = 1.0;
        let fx_obs = 1.51;
        let fx_exp = 1.0;
        assert!(!ratio_filter_passes(
            "species", br_obs, br_exp, fx_obs, fx_exp, None, None
        ));
        assert!(!ratio_filter_passes(
            "genus", br_obs, br_exp, fx_obs, fx_exp, None, None
        ));
        // At exactly 1.5 both pass (inclusive upper bound per spec).
        assert!(ratio_filter_passes(
            "species", 1.5, 1.0, 1.5, 1.0, None, None
        ));
        assert!(ratio_filter_passes("genus", 1.5, 1.0, 1.5, 1.0, None, None));
    }

    #[test]
    fn ratio_filter_requires_positive_expectations() {
        // Non-positive expected values must fail regardless of observed.
        assert!(!ratio_filter_passes(
            "species", 1.0, 0.0, 1.0, 1.0, None, None
        ));
        assert!(!ratio_filter_passes(
            "genus", 1.0, 1.0, 1.0, 0.0, None, None
        ));
    }

    #[test]
    fn ratio_filter_honours_independent_lower_bound_overrides() {
        assert!(ratio_filter_passes(
            "species",
            0.80,
            1.0,
            0.78,
            1.0,
            Some(0.80),
            Some(0.78),
        ));
        assert!(!ratio_filter_passes(
            "species",
            0.79,
            1.0,
            0.78,
            1.0,
            Some(0.80),
            Some(0.78),
        ));
        assert!(!ratio_filter_passes(
            "species",
            0.80,
            1.0,
            0.77,
            1.0,
            Some(0.80),
            Some(0.78),
        ));
    }

    // ------------------------------------------------------------------
    // Batch assignment: species precedence & genus fallback thresholds
    // ------------------------------------------------------------------

    /// Helper: run the batch assignment on a single read and return the
    /// resulting `ReadAssignment` (if any). Uses an empty
    /// species-normalization map (no sub-species remapping).
    fn run_single(
        cfg: &ProfilerConfig,
        asm_lengths: &HashMap<String, f64>,
        s2g: &HashMap<u32, u32>,
        fallback_mode: bool,
        read_name: &str,
        alns: Vec<RawAln>,
    ) -> Option<ReadAssignment> {
        run_single_with_normalize(
            cfg,
            asm_lengths,
            &HashMap::new(),
            s2g,
            fallback_mode,
            read_name,
            alns,
        )
    }

    /// Helper that additionally accepts a `species_normalize` map so
    /// tests can exercise the sub-species → species remapping path.
    fn run_single_with_normalize(
        cfg: &ProfilerConfig,
        asm_lengths: &HashMap<String, f64>,
        species_normalize: &HashMap<u32, u32>,
        s2g: &HashMap<u32, u32>,
        fallback_mode: bool,
        read_name: &str,
        alns: Vec<RawAln>,
    ) -> Option<ReadAssignment> {
        // Pre-parse each synthetic alignment as the mapping loader would —
        // rejects are filtered out, mirroring the real pipeline where
        // unparsable records never reach the batch assignment.
        let alignments: Vec<AlignmentEntry> = alns
            .iter()
            .filter_map(|raw| build_aln(raw, asm_lengths, species_normalize))
            .collect();
        let bundle = ReadBundle {
            read_name: Arc::<str>::from(read_name),
            alignments,
        };
        let policy = EvidencePolicy::from_config(cfg).expect("read input");
        let mut result = process_evidence_batch(&[bundle], cfg, policy, fallback_mode, s2g, None);
        result.read_assignments.pop()
    }

    fn mk_asm_map() -> HashMap<String, f64> {
        let mut m = HashMap::new();
        m.insert("asmA".to_string(), 1000.0);
        m.insert("asmB".to_string(), 1000.0);
        m
    }

    fn mk_genus_map() -> HashMap<u32, u32> {
        let mut m = HashMap::new();
        m.insert(101, 5001); // species 101 -> genus 5001
        m.insert(102, 5001); // species 102 -> genus 5001
        m.insert(103, 5002); // species 103 -> genus 5002
        m
    }

    /// (identity, mapped fraction) of a species-qualified hit.
    const PERFECT: (f64, f64) = (1.0, 1.0);
    /// (identity, mapped fraction) of a hit for genus fallback.
    const NEAR_MISS: (f64, f64) = (0.92, 1.0);
    /// Mapped fraction 0.5, below GENUS_FALLBACK_MIN_FRACTION.
    const LOW_FRACTION: (f64, f64) = (1.0, 0.5);

    // ------------------------------------------------------------------
    // Sub-species -> species normalization
    // ------------------------------------------------------------------

    #[test]
    fn subspecies_taxid_is_normalized_to_species() {
        // Raw parts[2] = 9001 is a sub-species of species 101.
        // With `species_normalize = { 9001 -> 101 }` the read must be
        // assigned to 101, not 9001. This is independent of fallback.
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        let mut norm = HashMap::new();
        norm.insert(9001u32, 101u32);
        let alns = vec![aln_entry(
            "asmA|9001|9001|accA|1000",
            PERFECT,
            100,
            (1, 100),
        )];
        // With fallback OFF (previous-logic path): still normalizes.
        let no_fb = run_single_with_normalize(
            &cfg,
            &asm,
            &norm,
            &s2g,
            /*fallback*/ false,
            "r1",
            alns.clone(),
        )
        .expect("species assignment");
        assert_eq!(no_fb.rank, AssignRank::Species);
        assert_eq!(no_fb.taxids, vec![101]);

        // With fallback ON: species precedence still picks species 101
        // (not the raw sub-species id and not genus).
        let with_fb = run_single_with_normalize(&cfg, &asm, &norm, &s2g, true, "r2", alns)
            .expect("species assignment");
        assert_eq!(with_fb.rank, AssignRank::Species);
        assert_eq!(with_fb.taxids, vec![101]);
    }

    #[test]
    fn multiple_subspecies_of_same_species_merge_into_one() {
        // Two alignments hit different sub-species (9001 and 9002) of
        // the same species 101. They must collapse into a single
        // species-level assignment with taxid 101.
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        let mut norm = HashMap::new();
        norm.insert(9001u32, 101u32);
        norm.insert(9002u32, 101u32);
        let alns = vec![
            aln_entry("asmA|9001|9001|accA|1000", PERFECT, 100, (1, 100)),
            aln_entry("asmB|9002|9002|accB|1000", PERFECT, 100, (1, 100)),
        ];
        let assignment = run_single_with_normalize(&cfg, &asm, &norm, &s2g, false, "r1", alns)
            .expect("species assignment");
        assert_eq!(assignment.rank, AssignRank::Species);
        // One unique taxid (101), not [9001, 9002] and not [101, 101].
        assert_eq!(assignment.taxids, vec![101]);
    }

    #[test]
    fn missing_normalize_entry_keeps_raw_taxid() {
        // When the raw taxid has no entry in `species_normalize` the
        // raw id is used as-is. This preserves the original behaviour
        // for references whose third field is already a species taxid
        // (or is above species with no species ancestor).
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        let norm: HashMap<u32, u32> = HashMap::new();
        let alns = vec![aln_entry("asmA|1|101|accA|1000", PERFECT, 100, (1, 100))];
        let assignment = run_single_with_normalize(&cfg, &asm, &norm, &s2g, false, "r1", alns)
            .expect("species assignment");
        assert_eq!(assignment.rank, AssignRank::Species);
        assert_eq!(assignment.taxids, vec![101]);
    }

    #[test]
    fn subspecies_normalization_applies_to_genus_fallback_key() {
        // Only a near-miss alignment (genus-fallback candidate) against
        // sub-species 9001, which normalizes to species 101 whose genus
        // is 5001. The genus fallback path must therefore emit genus
        // 5001 (derived via the normalized species), not whatever
        // s2g[9001] would have been (s2g is keyed by species here).
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map(); // only 101, 102, 103 are species keys
        let mut norm = HashMap::new();
        norm.insert(9001u32, 101u32);
        let alns = vec![aln_entry(
            "asmA|9001|9001|accA|1000",
            NEAR_MISS,
            100,
            (1, 100),
        )];
        let assignment = run_single_with_normalize(&cfg, &asm, &norm, &s2g, true, "r1", alns)
            .expect("genus assignment");
        assert_eq!(assignment.rank, AssignRank::Genus);
        assert_eq!(assignment.taxids, vec![5001]);
    }

    #[test]
    fn species_first_precedence_overrides_genus_fallback() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // Two alignments: first qualifies as species (perfect), second only
        // as genus fallback. Species must win even with fallback_mode=true.
        let alns = vec![
            aln_entry("asmA|1|101|accA|1000", PERFECT, 100, (1, 100)),
            aln_entry("asmB|1|103|accB|1000", NEAR_MISS, 100, (1, 100)),
        ];
        let assignment =
            run_single(&cfg, &asm, &s2g, /*fallback*/ true, "r1", alns).expect("assignment");
        assert_eq!(assignment.rank, AssignRank::Species);
        assert_eq!(assignment.taxids, vec![101]);
    }

    #[test]
    fn genus_fallback_only_when_no_species_qualifies() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // Only a near-miss alignment: passes genus thresholds, not species.
        let alns = vec![aln_entry("asmA|1|101|accA|1000", NEAR_MISS, 100, (1, 100))];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns).expect("assignment");
        assert_eq!(assignment.rank, AssignRank::Genus);
        assert_eq!(assignment.taxids, vec![5001]);
    }

    #[test]
    fn genus_fallback_skipped_when_gate_inactive() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // Same near-miss alignment; fallback_mode = false => no assignment.
        let alns = vec![aln_entry("asmA|1|101|accA|1000", NEAR_MISS, 100, (1, 100))];
        let assignment = run_single(&cfg, &asm, &s2g, /*fallback*/ false, "r1", alns);
        assert!(
            assignment.is_none(),
            "expected no assignment when fallback is disabled"
        );
    }

    #[test]
    fn genus_fallback_uses_fixed_mapped_len_and_fraction_constants() {
        // Even when the user sets cfg.mapped_len and cfg.fraction higher
        // than the fallback constants, the genus path must still accept
        // mapped_len >= 50 and fraction >= 0.60 alignments.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.mapped_len = 200; // aggressive user threshold
        cfg.fraction = 0.95; //  aggressive user threshold
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // Identity 0.92, fraction 1.0, mapped_len 60 => species disqualified
        // (mapped_len < cfg.mapped_len), but genus-qualified (mapped_len >= 50,
        // fraction >= 0.6, identity >= cfg.genus_identity).
        let alns = vec![aln_entry("asmA|1|101|accA|1000", NEAR_MISS, 60, (1, 60))];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns).expect("assignment");
        assert_eq!(assignment.rank, AssignRank::Genus);
        assert_eq!(assignment.taxids, vec![5001]);
    }

    #[test]
    fn genus_fallback_accepts_high_identity_alignment_that_fails_species_fraction() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.identity = 0.95;
        cfg.fraction = 0.80;
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // identity = 1.0 passes species identity, but fraction = 0.70 fails
        // the species fraction cutoff while still passing the genus fallback
        // fraction cutoff. High identity must not be rejected by an upper-bound
        // species-identity check in the fallback path.
        let alns = vec![aln_entry("asmA|1|101|accA|1000", (1.0, 0.7), 70, (1, 70))];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns)
            .expect("high-identity partial alignment should fall back to genus");
        assert_eq!(assignment.rank, AssignRank::Genus);
        assert_eq!(assignment.taxids, vec![5001]);
    }

    #[test]
    fn genus_fallback_rejected_when_mapped_len_below_50() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // mapped_len = 49 violates the fixed 50 floor for genus fallback.
        let alns = vec![aln_entry("asmA|1|101|accA|1000", NEAR_MISS, 49, (1, 49))];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns);
        assert!(
            assignment.is_none(),
            "expected no assignment when mapped_len < 50"
        );
    }

    #[test]
    fn genus_fallback_long_reads_pass_on_length_or_fraction() {
        // For long reads the coverage gate is length OR fraction. A Nanopore
        // read whose identity misses the species cutoff falls to the genus
        // fallback, which now accepts it on EITHER the mapped-length floor
        // (cfg.mapped_len) or the fraction floor (GENUS_FALLBACK_MIN_FRACTION),
        // and rejects it only when it meets neither.
        let mut cfg = make_cfg(Sequencer::Nanopore, None);
        cfg.mapped_len = 100;
        cfg.fraction = 0.6;
        cfg.identity = 0.95; // NEAR_MISS identity 0.92 < 0.95 -> falls to genus
        let asm = mk_asm_map();
        let s2g = mk_genus_map();

        // mapped_len 80 < 100 (length fails) but NEAR_MISS fraction = 1.0:
        // accepted via fraction.
        let via_fraction = vec![aln_entry("asmA|1|101|accA|1000", NEAR_MISS, 80, (1, 80))];
        assert_eq!(
            run_single(&cfg, &asm, &s2g, true, "r1", via_fraction)
                .expect("accepted on fraction")
                .rank,
            AssignRank::Genus
        );

        // mapped_len 100 >= 100: accepted via length.
        let via_length = vec![aln_entry("asmA|1|101|accA|1000", NEAR_MISS, 100, (1, 100))];
        assert_eq!(
            run_single(&cfg, &asm, &s2g, true, "r1", via_length)
                .expect("accepted on length")
                .rank,
            AssignRank::Genus
        );

        // mapped_len 80 < 100 AND LOW_FRACTION fraction 0.5 < 0.60: both fail ->
        // rejected.
        let neither = vec![aln_entry("asmA|1|101|accA|1000", LOW_FRACTION, 80, (1, 80))];
        assert!(
            run_single(&cfg, &asm, &s2g, true, "r1", neither).is_none(),
            "long reads: genus fallback rejects an alignment failing both length and fraction"
        );
    }

    #[test]
    fn genus_fallback_rejected_when_fraction_below_060() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        let s2g = mk_genus_map();
        // Fraction = 50/100 = 0.5 (<0.60) but identity 1.0 and
        // mapped_len above species threshold. Note this actually qualifies
        // for species (identity >= cfg.identity=0.95, mapped_len >= 50) *except*
        // fraction < cfg.fraction => species disqualified. Then fallback:
        // fraction < 0.60 => also disqualified.
        let alns = vec![aln_entry("asmA|1|101|accA|1000", (1.0, 0.5), 50, (1, 50))];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns);
        assert!(
            assignment.is_none(),
            "expected no assignment when fraction < 0.60"
        );
    }

    #[test]
    fn genus_fallback_identity_cutoff_is_cfg_genus_identity() {
        // identity of a NEAR_MISS = 92/100 = 0.92.
        let asm = mk_asm_map();
        let s2g = mk_genus_map();

        // With the default genus_identity=0.80, a 0.92 near-miss qualifies.
        let cfg_default = make_cfg(Sequencer::Illumina, None);
        assert!((cfg_default.genus_identity - 0.80).abs() < 1e-12);
        let alns = vec![aln_entry("asmA|1|101|accA|1000", NEAR_MISS, 100, (1, 100))];
        let assignment = run_single(&cfg_default, &asm, &s2g, true, "r1", alns.clone())
            .expect("default genus_identity should accept 0.92 near-miss");
        assert_eq!(assignment.rank, AssignRank::Genus);

        // Raising genus_identity above the alignment identity rejects it.
        let mut cfg_strict = make_cfg(Sequencer::Illumina, None);
        cfg_strict.genus_identity = 0.93;
        let assignment = run_single(&cfg_strict, &asm, &s2g, true, "r1", alns.clone());
        assert!(
            assignment.is_none(),
            "identity 0.92 must be rejected when genus_identity=0.93"
        );

        // A 0.82 alignment: ACCEPTED at the new default 0.80 cutoff,
        // REJECTED if the cutoff is raised to 0.85. This exercises the
        // boundary around the old (0.85) and new (0.80) defaults.
        let boundary = vec![aln_entry(
            "asmA|1|101|accA|1000",
            (0.82, 1.0),
            100,
            (1, 100),
        )];
        let assignment = run_single(&cfg_default, &asm, &s2g, true, "r1", boundary.clone())
            .expect("default genus_identity=0.80 should accept 0.82 alignment");
        assert_eq!(assignment.rank, AssignRank::Genus);

        let mut cfg_old = make_cfg(Sequencer::Illumina, None);
        cfg_old.genus_identity = 0.85;
        let assignment = run_single(&cfg_old, &asm, &s2g, true, "r1", boundary);
        assert!(
            assignment.is_none(),
            "0.82 identity must be rejected when genus_identity is raised to 0.85"
        );

        // Truly low-identity alignment (0.78) must still be rejected at
        // the new default 0.80 cutoff.
        let cfg_default = make_cfg(Sequencer::Illumina, None);
        let below_default = vec![aln_entry(
            "asmA|1|101|accA|1000",
            (0.78, 1.0),
            100,
            (1, 100),
        )];
        let assignment = run_single(&cfg_default, &asm, &s2g, true, "r1", below_default);
        assert!(
            assignment.is_none(),
            "0.78 identity must be rejected at default genus_identity=0.80"
        );
    }

    #[test]
    fn genus_fallback_skips_alignment_without_resolvable_genus() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        let asm = mk_asm_map();
        // species 999 is missing from the s2g map (no genus resolvable).
        let s2g: HashMap<u32, u32> = HashMap::new();
        let alns = vec![aln_entry("asmA|1|999|accA|1000", NEAR_MISS, 100, (1, 100))];
        let assignment = run_single(&cfg, &asm, &s2g, true, "r1", alns);
        assert!(
            assignment.is_none(),
            "expected no assignment when genus cannot be resolved"
        );
    }

    #[test]
    fn representative_mapping_carries_its_genome_scale() {
        let lower_identity = AlignmentEntry {
            asm: Arc::<str>::from("full-genome"),
            acc: Arc::<str>::from("contig-a"),
            profile_taxid: 100,
            profile_rank: ProfileRank::Species,
            reference_id: u32::MAX,
            assembly_id: u32::MAX,
            role: MappingRole::Unpaired,
            pair_id: 0,
            candidate_boundary: metam::mapper::CandidateBoundary::COMPLETE,
            sampling_scale: 1.0,
            asm_len: 1_000.0,
            identity: 0.98,
            fraction: 1.0,
            mapped_len: 100,
            read_len: 100,
            span: (1, 101),
            fixed_chunk: (1, 101),
            flex_chunk: (1, 101),
            marker_gene_len: 0,
        };
        let mut higher_identity = lower_identity.clone();
        higher_identity.asm = Arc::<str>::from("fractional-genome");
        higher_identity.sampling_scale = 10.0;
        higher_identity.identity = 0.99;

        let mut best = BestMap::new();
        let mut order = Vec::new();
        let mut intervals = IntervalMap::new();
        update_best(
            &mut best,
            &mut order,
            &mut intervals,
            100,
            f64::from(lower_identity.identity),
            0.1,
            &lower_identity,
        );
        update_best(
            &mut best,
            &mut order,
            &mut intervals,
            100,
            f64::from(higher_identity.identity),
            0.1,
            &higher_identity,
        );

        assert_eq!(best.get(&100).map(|entry| entry.3), Some(10.0));
        assert_eq!(
            intervals.get(&100).map(|entry| entry.0.asm.as_ref()),
            Some("fractional-genome")
        );
    }
}

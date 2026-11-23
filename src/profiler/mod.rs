//! Profiler Module - Taxonomic abundance estimation from mapping data
//!
//! The mappings of a sample are loaded with the reference metadata
//! (`community.rs`) and profiled by the level-wise engine
//! (`community/levels.rs`): species first, then the genomes within each kept
//! species (`--strain`), and the fragments that no kept species takes at the
//! lowest rank above species that their matches and evidence support. Each
//! reference is reported at the most specific rank its taxonomy label
//! supports.
//!
//! # Coverage Metrics
//!
//! The profiler computes several coverage metrics:
//!
//! - **Breadth**: Fraction of genome covered by at least one read
//! - **Depth**: Average number of reads covering each position
//! - **Chunk breadth**: Coverage of fixed-size genome chunks (for QC)
//! - **Flex chunk breadth**: Coverage of variable-size chunks (based on √genome_size)

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::hash::BuildHasher;
use std::sync::Arc;

use anyhow::{Context, Result};
use clap::ValueEnum;
use xxhash_rust::xxh3::Xxh3DefaultBuilder;

use crate::db::{DatabaseMetadata, ReferenceType};
use crate::mapper::tsv::{MetamMappingReceiver, SharedMetamProfile};
use crate::progress::ProgressOptions;
use crate::sequencer::Sequencer;
use crate::taxonomy::Taxonomy;

pub(crate) mod community;
pub(crate) mod evidence;
mod likelihood_em;
mod markers;
mod output;
mod pathogen;
mod rollup;
mod statistics;
mod unclassified;

/// Read-name aggregation map with a stable, high-throughput hash builder.
///
/// Its hasher has a fixed seed, so that the reads of a mapping stream come out
/// in the same order in every run, and the sums taken in that order are the
/// same. XXH3 is also faster than the default keyed hasher on short read
/// names.
pub(crate) type StableReadMap<T> = HashMap<String, Vec<T>, Xxh3DefaultBuilder>;

#[repr(u8)]
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum ProfileRank {
    Species,
    Genus,
    Family,
    Order,
    Class,
    Phylum,
    Kingdom,
}

impl ProfileRank {
    pub(crate) fn from_prefix(prefix: char) -> Option<Self> {
        match prefix {
            's' => Some(Self::Species),
            'g' => Some(Self::Genus),
            'f' => Some(Self::Family),
            'o' => Some(Self::Order),
            'c' => Some(Self::Class),
            'p' => Some(Self::Phylum),
            'k' => Some(Self::Kingdom),
            _ => None,
        }
    }

    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Species => "species",
            Self::Genus => "genus",
            Self::Family => "family",
            Self::Order => "order",
            Self::Class => "class",
            Self::Phylum => "phylum",
            Self::Kingdom => "kingdom",
        }
    }

    pub(crate) fn from_taxonomy_rank(rank: &str) -> Option<Self> {
        match rank {
            "species" => Some(Self::Species),
            "genus" => Some(Self::Genus),
            "family" => Some(Self::Family),
            "order" => Some(Self::Order),
            "class" => Some(Self::Class),
            "phylum" => Some(Self::Phylum),
            "kingdom" | "superkingdom" => Some(Self::Kingdom),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct ReferenceTaxonomy {
    pub primary_taxid: Option<u32>,
    pub profile_taxid: u32,
    pub profile_rank: ProfileRank,
}

pub(crate) fn parse_reference_taxonomy(value: &str) -> Result<ReferenceTaxonomy> {
    let mut primary_taxid = None;
    let mut by_rank: HashMap<char, u32> = HashMap::new();
    for token in value.split(';') {
        let Some((prefix, taxid)) = token.split_once(':') else {
            anyhow::bail!("taxonomy token '{}' is missing ':'", token);
        };
        let mut chars = prefix.chars();
        let Some(rank) = chars.next() else {
            anyhow::bail!("empty taxonomy prefix in '{}'", token);
        };
        if chars.next().is_some() {
            anyhow::bail!("unsupported taxonomy prefix '{}' in '{}'", prefix, value);
        }
        if !matches!(rank, 't' | 's' | 'g' | 'f' | 'o' | 'c' | 'p' | 'k') {
            anyhow::bail!("unsupported taxonomy prefix '{}' in '{}'", prefix, value);
        }
        let taxid = taxid
            .parse::<u32>()
            .with_context(|| format!("invalid taxid '{}' in taxonomy '{}'", taxid, value))?;
        if taxid == 0 {
            anyhow::bail!("taxid must be greater than zero in taxonomy '{}'", value);
        }
        if by_rank.insert(rank, taxid).is_some() {
            anyhow::bail!("duplicate taxonomy prefix '{}' in '{}'", rank, value);
        }
        if rank == 't' {
            primary_taxid = Some(taxid);
        }
    }
    if primary_taxid.is_some() && !by_rank.contains_key(&'s') {
        anyhow::bail!(
            "taxonomy '{}' has t:<strain_taxid> but no s:<species_taxid>",
            value
        );
    }
    let (profile_rank, profile_taxid) = if let Some(species) = by_rank.get(&'s') {
        (ProfileRank::Species, *species)
    } else {
        ['g', 'f', 'o', 'c', 'p', 'k']
            .iter()
            .find_map(|rank| {
                by_rank.get(rank).copied().map(|taxid| {
                    (
                        ProfileRank::from_prefix(*rank)
                            .expect("supported profile-rank prefix must have an enum value"),
                        taxid,
                    )
                })
            })
            .ok_or_else(|| anyhow::anyhow!("taxonomy '{}' contains no supported taxid", value))?
    };
    Ok(ReferenceTaxonomy {
        primary_taxid,
        profile_taxid,
        profile_rank,
    })
}

#[derive(Clone, Copy, Debug)]
struct ResolvedProfileEndpoint {
    current_taxid: u32,
    profile_taxid: u32,
    profile_rank: ProfileRank,
    strain_like: bool,
}

/// Cached reconciliation of reference-header taxonomy against the dump used
/// for this profiling run.
///
/// Header prefixes describe taxonomy at index-build time. Taxids can be
/// merged or change rank in later NCBI dumps, so the current dump is the
/// authority. Each distinct taxid and lineage is resolved once; repeated
/// contigs use only in-memory hash lookups.
pub(crate) struct ReferenceTaxonomyResolver<'a> {
    taxonomy: &'a Taxonomy,
    current_taxids: HashMap<u32, u32>,
    profile_endpoints: HashMap<u32, ResolvedProfileEndpoint>,
    resolved_headers: HashMap<ReferenceTaxonomy, ReferenceTaxonomy>,
}

impl<'a> ReferenceTaxonomyResolver<'a> {
    pub(crate) fn new(taxonomy: &'a Taxonomy) -> Self {
        Self {
            taxonomy,
            current_taxids: HashMap::new(),
            profile_endpoints: HashMap::new(),
            resolved_headers: HashMap::new(),
        }
    }

    #[inline]
    fn current_taxid(&mut self, taxid: u32) -> Option<u32> {
        if let Some(current) = self.current_taxids.get(&taxid).copied() {
            return Some(current);
        }
        let current = self.taxonomy.resolve_current_taxid(taxid)?;
        self.current_taxids.insert(taxid, current);
        Some(current)
    }

    fn profile_endpoint(
        &mut self,
        current_taxid: u32,
        declared_taxid: u32,
        reference: &str,
    ) -> Result<ResolvedProfileEndpoint> {
        if let Some(endpoint) = self.profile_endpoints.get(&current_taxid).copied() {
            return Ok(endpoint);
        }
        let lineage = self.taxonomy.get_parents(current_taxid);
        let current_rank = lineage
            .first()
            .map(|(_, _, rank)| rank.as_str())
            .expect("a resolved taxonomy node must start a lineage");
        let (profile_taxid, profile_rank) = lineage
            .iter()
            .find_map(|(taxid, _, rank)| {
                ProfileRank::from_taxonomy_rank(rank)
                    .map(|profile_rank| (*taxid, profile_rank))
            })
            .with_context(|| {
                format!(
                    "reference '{}' taxid {} resolves to current taxid {} at rank '{}', but neither it nor its ancestors has a supported profile rank",
                    reference, declared_taxid, current_taxid, current_rank,
                )
            })?;
        let endpoint = ResolvedProfileEndpoint {
            current_taxid,
            profile_taxid,
            profile_rank,
            strain_like: matches!(current_rank, "subspecies" | "strain"),
        };
        self.profile_endpoints.insert(current_taxid, endpoint);
        Ok(endpoint)
    }

    pub(crate) fn resolve(
        &mut self,
        declared: ReferenceTaxonomy,
        reference: &str,
    ) -> Result<ReferenceTaxonomy> {
        if let Some(resolved) = self.resolved_headers.get(&declared).copied() {
            return Ok(resolved);
        }

        // A stale strain/subspecies taxid does not invalidate an otherwise
        // current species endpoint. Conversely, when a stale species taxid
        // was deleted without a merge, a current `t:` lineage can recover its
        // current species ancestor.
        let current_primary_taxid = declared
            .primary_taxid
            .and_then(|taxid| self.current_taxid(taxid));
        let current_profile_taxid = self.current_taxid(declared.profile_taxid);
        let endpoint_source = current_profile_taxid
            .or(current_primary_taxid)
            .with_context(|| {
                format!(
                    "reference '{}' taxonomy cannot be resolved in the current dump: terminal {}:{} and primary {:?} are absent, deleted, or have invalid merged-taxid chains",
                    reference,
                    declared.profile_rank.as_str(),
                    declared.profile_taxid,
                    declared.primary_taxid,
                )
            })?;
        let endpoint =
            match self.profile_endpoint(endpoint_source, declared.profile_taxid, reference) {
                Ok(endpoint) => endpoint,
                Err(profile_error)
                    if current_profile_taxid.is_some()
                        && current_primary_taxid.is_some()
                        && current_primary_taxid != current_profile_taxid =>
                {
                    self.profile_endpoint(
                        current_primary_taxid.expect("guarded as present"),
                        declared
                            .primary_taxid
                            .expect("current primary requires declared primary"),
                        reference,
                    )
                    .with_context(|| {
                        format!("terminal-taxid lineage also failed: {profile_error:#}")
                    })?
                }
                Err(error) => return Err(error),
            };

        let primary_taxid = current_primary_taxid.or_else(|| {
            (declared.primary_taxid.is_none() && endpoint.strain_like)
                .then_some(endpoint.current_taxid)
        });
        let resolved = ReferenceTaxonomy {
            primary_taxid,
            profile_taxid: endpoint.profile_taxid,
            profile_rank: endpoint.profile_rank,
        };
        self.resolved_headers.insert(declared, resolved);
        Ok(resolved)
    }
}

pub(crate) type ReferenceTaxonomyAdjustments =
    BTreeMap<(ReferenceTaxonomy, ReferenceTaxonomy), usize>;

pub(crate) fn record_reference_taxonomy_adjustment(
    adjustments: &mut ReferenceTaxonomyAdjustments,
    declared: ReferenceTaxonomy,
    resolved: ReferenceTaxonomy,
) {
    if declared != resolved {
        *adjustments.entry((declared, resolved)).or_default() += 1;
    }
}

fn taxonomy_adjustment_description(
    declared: ReferenceTaxonomy,
    resolved: ReferenceTaxonomy,
) -> String {
    let primary = if declared.primary_taxid == resolved.primary_taxid {
        format!(
            "primary remains {}",
            declared
                .primary_taxid
                .map_or_else(|| "none".to_string(), |taxid| taxid.to_string())
        )
    } else {
        format!(
            "primary {} → {}",
            declared
                .primary_taxid
                .map_or_else(|| "none".to_string(), |taxid| taxid.to_string()),
            resolved
                .primary_taxid
                .map_or_else(|| "none".to_string(), |taxid| taxid.to_string())
        )
    };
    let terminal = if declared.profile_taxid == resolved.profile_taxid
        && declared.profile_rank == resolved.profile_rank
    {
        format!(
            "terminal remains {}:{}",
            declared.profile_rank.as_str(),
            declared.profile_taxid
        )
    } else {
        format!(
            "terminal {}:{} → {}:{}",
            declared.profile_rank.as_str(),
            declared.profile_taxid,
            resolved.profile_rank.as_str(),
            resolved.profile_taxid
        )
    };
    format!("{primary}; {terminal}")
}

/// References skipped because the taxonomy dump cannot resolve their
/// taxids, which is allowed only with --skip-missing-taxa.
#[derive(Default)]
pub(crate) struct UnresolvedReferences {
    references: usize,
    genomes: BTreeSet<String>,
}

impl UnresolvedReferences {
    /// Skip a reference of `genome` whose taxonomy failed to resolve with
    /// `error`, or return the error when skipping is not allowed.
    pub(crate) fn skip(&mut self, genome: &str, error: anyhow::Error, allowed: bool) -> Result<()> {
        if !allowed {
            return Err(error);
        }
        self.references += 1;
        if !self.genomes.contains(genome) {
            self.genomes.insert(genome.to_string());
        }
        Ok(())
    }

    pub(crate) fn log(&self) {
        if self.references == 0 {
            return;
        }
        const LISTED: usize = 5;
        let mut listed = self
            .genomes
            .iter()
            .take(LISTED)
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join(", ");
        if self.genomes.len() > LISTED {
            listed.push_str(", ...");
        }
        log::warn!(
            target: "PROFILER",
            "Ignored {} reference sequence(s) of {} genome(s) whose taxids are not in the taxonomy dump ({}); reads that map only to them stay unassigned, and reads of these organisms may be assigned to relatives.",
            self.references,
            self.genomes.len(),
            listed
        );
    }
}

pub(crate) fn log_reference_taxonomy_adjustments(
    adjustments: &ReferenceTaxonomyAdjustments,
    verbose: bool,
) {
    if adjustments.is_empty() {
        return;
    }
    let sequence_count: usize = adjustments.values().sum();
    log::warn!(
        target: "PROFILER",
        "The current taxonomy dump reconciled taxonomy metadata for {} database sequence(s) across {} transition(s); dump taxonomy is authoritative.",
        sequence_count,
        adjustments.len(),
    );

    if !verbose {
        return;
    }

    let mut transitions: Vec<_> = adjustments.iter().collect();
    transitions
        .sort_unstable_by(|left, right| right.1.cmp(left.1).then_with(|| left.0.cmp(right.0)));
    const MAX_LOGGED_TRANSITIONS: usize = 8;
    for ((declared, resolved), count) in transitions.iter().take(MAX_LOGGED_TRANSITIONS) {
        log::info!(
            target: "PROFILER",
            "Taxonomy updated for {} database sequence(s): {}.",
            count,
            taxonomy_adjustment_description(*declared, *resolved),
        );
    }
    if transitions.len() > MAX_LOGGED_TRANSITIONS {
        log::info!(
            target: "PROFILER",
            "{} additional taxonomy transition(s) were summarized but not logged individually.",
            transitions.len() - MAX_LOGGED_TRANSITIONS,
        );
    }
}

/// Return the full-index-equivalent support multiplier for one genome.
/// Fully indexed genomes deliberately receive exactly `1.0`.
pub(crate) fn sampling_scale_from_sizes(
    full_genome_size: usize,
    effective_indexed_size: usize,
) -> f64 {
    if full_genome_size == 0
        || effective_indexed_size == 0
        || effective_indexed_size >= full_genome_size
    {
        1.0
    } else {
        full_genome_size as f64 / effective_indexed_size as f64
    }
}

/// Record one scale per genome and reject internally inconsistent metadata.
pub(crate) fn record_genome_sampling_scale(
    scales: &mut HashMap<String, f64>,
    genome: &str,
    scale: f64,
) -> Result<()> {
    if genome.is_empty() || !scale.is_finite() || scale < 1.0 {
        anyhow::bail!("invalid sampling scale {scale} for genome '{genome}'");
    }
    if let Some(previous) = scales.get(genome).copied() {
        let tolerance = 1e-12 * previous.abs().max(scale.abs()).max(1.0);
        if (previous - scale).abs() > tolerance {
            anyhow::bail!(
                "inconsistent sampling scales for genome '{}': {} and {}",
                genome,
                previous,
                scale
            );
        }
    } else {
        scales.insert(genome.to_owned(), scale);
    }
    Ok(())
}

pub(crate) fn equivalent_read_count_passes_minimum(
    equivalent_read_count: f64,
    minimum: Option<usize>,
) -> bool {
    match minimum {
        Some(minimum) => {
            equivalent_read_count.is_finite() && equivalent_read_count >= minimum as f64
        }
        None => true,
    }
}

pub(crate) fn strip_pair_suffix(read_id: &str) -> &str {
    match read_id.rsplit_once('/') {
        Some((base, "1" | "2")) => base,
        _ => read_id,
    }
}

pub(crate) fn push_string_keyed_value<T, S>(
    map: &mut HashMap<String, Vec<T>, S>,
    key: &str,
    value: T,
) where
    S: BuildHasher,
{
    if let Some(values) = map.get_mut(key) {
        values.push(value);
    } else {
        map.insert(key.to_owned(), vec![value]);
    }
}

/// Profiling mode that controls the sensitivity/specificity trade-off.
///
/// Each mode sets different default thresholds for mapping filtering:
///
/// - **Recall**: Lower thresholds for higher sensitivity
///   - Accepts more mappings
///
/// - **Precision**: Higher thresholds for higher specificity
///   - Stricter filtering
///
/// - **Default**: Balanced settings for general use
///   - Good trade-off between sensitivity and specificity
///   - Recommended for most applications
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ProfileMode {
    /// High sensitivity mode - more permissive filtering
    Recall,
    /// High specificity mode - stricter filtering
    Precision,
    /// Balanced mode - recommended for general use
    Default,
}

impl ProfileMode {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Recall => "recall",
            Self::Precision => "precision",
            Self::Default => "default",
        }
    }
}

/// Configuration for the taxonomy profiler.
///
/// This structure holds all parameters that control profiling behavior,
/// from input/output paths to filtering thresholds to algorithm options.
/// Values here are the resolved settings of the run.
///
/// # Filtering Thresholds
///
/// The key thresholds that determine which mappings are retained:
///
/// - `identity`: Minimum mapper-reported per-base identity for native mappings
/// - `mapped_len`: Minimum number of mapped read bases
/// - `fraction`: Minimum fraction of read that must be mapped
/// - `breadth`: Minimum genome breadth of coverage to report a taxon
/// - `chunk_breadth`: Minimum chunked coverage
/// - `min_reads`: Minimum profile read count to report a taxon
/// - `min_oebr`: Minimum observed/expected breadth ratio
/// - `min_coebr`: Minimum observed/expected chunk breadth ratio
#[derive(Debug, Clone)]
pub struct ProfilerConfig {
    /// Also profile the genomes within each kept species (`--strain`).
    pub genome_level: bool,

    /// Path to a Metax mapping TSV. Unset when mappings are streamed from the
    /// mapper.
    pub mapper_mappings: Option<String>,

    /// In-memory mapping counters for the default streamed mapper path.
    pub mapping_profile_memory: Option<SharedMetamProfile>,

    /// Metax database metadata used to interpret native mapping reference IDs.
    /// File-backed mapping output can provide this metadata directly.
    pub db_metadata: Option<Arc<DatabaseMetadata>>,

    /// Sequencing platform (Illumina, Nanopore, PacBio)
    /// Affects default threshold selection
    pub sequencer: Sequencer,

    /// Maximum EM iterations for ambiguous-read redistribution.
    pub max_em_iterations: usize,

    /// Skip malformed Metam mapping rows instead of aborting the run.
    pub lenient_mapping_input: bool,

    /// Whether input is paired-end Illumina data
    /// Affects read name parsing and coverage calculation
    pub is_paired: bool,

    /// Minimum mapping identity threshold (0.0-1.0), against the mapper's
    /// canonical per-base identity.
    pub identity: f64,

    /// Claimed base accuracy behind an automatic long-read identity cutoff.
    /// When set, the cutoff is lowered if the reads align worse than this.
    pub identity_check_accuracy: Option<f64>,

    /// Long reads whose sampled mean predicted base accuracy is below the
    /// noisy-read limit (as the mapper decides): many reads of a present
    /// species then align below the species identity cutoff.
    pub noisy_long_reads: bool,

    /// Minimum number of mapped bases to keep a mapping
    pub mapped_len: usize,

    /// Minimum genome breadth of coverage to report a taxon
    /// If None, defaults to 0.0 (report all)
    pub breadth: Option<f64>,

    /// Minimum chunk breadth threshold
    /// If None, estimated automatically from read count
    pub chunk_breadth: Option<f64>,

    /// Multiplicative adjustment for an automatically estimated chunk-breadth
    /// cutoff. The CLI guarantees a finite value greater than zero.
    pub chunk_breadth_scale: f64,

    /// Minimum profile read count to report a taxon.
    pub min_reads: Option<usize>,

    /// Minimum observed/expected breadth ratio.
    pub min_oebr: Option<f64>,

    /// Minimum observed/expected chunk breadth ratio.
    pub min_coebr: Option<f64>,

    /// Minimum number of distinct marker genes with assigned reads for a
    /// species represented by marker genes.
    pub min_marker_genes: usize,

    /// Drop a species when a same-genus species with more depth has a mean
    /// read identity at least this much higher. Off when `None`.
    pub relative_identity_gap: Option<f64>,

    /// Report species removed by the species tests whose reads cover their
    /// genome evenly, and the fragments that no kept species takes at the
    /// lowest rank above species that their matches and their evidence
    /// support. On unless `--strict-coverage`.
    pub above_species: bool,

    /// Estimate the share of the input that the profile does not explain and
    /// append it as an `unclassified` row.
    pub estimate_unclassified: bool,

    /// On a fractional database, leave out a species, or report it at its
    /// genus, when the mean identity of its best placements lies more than
    /// this below the sample's baseline. Off when `None`, and the loader sets
    /// it off unless some genome is indexed in part.
    pub identity_margin: Option<f64>,

    /// Minimum fraction of read that must map (0.0-1.0)
    pub fraction: f64,

    /// Enable low-biomass sample heuristics
    /// Disables chunk breadth threshold
    pub lowbiomass: bool,

    /// Write how every taxon judged was decided (`<prefix>.decisions.tsv`).
    pub save_decisions: bool,

    /// Write optional per-read taxonomic assignments.
    pub save_classifications: bool,

    /// Force mapped-read basis for the chunk-breadth estimate.
    pub by_aligned: bool,

    /// Minimum identity of a read without a species-level candidate for the
    /// levels above species. Comparison is inclusive (`identity >= family_identity`).
    /// Default 0.80.
    pub family_identity: f64,

    /// Mapping-rate cut-off for the low-map path. A run is considered
    /// low-map (auto `--by-mapped` forced on) when
    /// `mapped_reads / total_reads < low_map_rate_threshold`. The comparison
    /// is strict. Default 0.30.
    pub low_map_rate_threshold: f64,

    /// `--host` with its pathogen table.
    pub host: Option<Host>,

    /// Number of threads for parallel processing
    pub threads: usize,

    /// Enable verbose logging
    pub verbose: bool,

    /// Enable very verbose logging with full commands
    pub very_verbose: bool,

    /// Live progress-bar behavior for profiling phases.
    pub progress: ProgressOptions,

    /// Suppress the live mapping-record loading counter. This is used for
    /// streamed in-process mapping, where the mapper already emits periodic
    /// progress logs for the same pipeline stage.
    pub suppress_mapping_load_progress: bool,

    /// Output file prefix
    pub outprefix: String,

    /// Profiling mode preset (recall, precision, default)
    pub mode: ProfileMode,

    /// Directory containing taxonomy dump files (optional)
    pub dmp_dir: Option<String>,

    /// Skip references whose taxids the taxonomy dump cannot resolve instead
    /// of failing (--skip-missing-taxa).
    pub ignore_unresolved_references: bool,
}

/// `--host`: the host's NCBI taxid, and the pathogen table
/// (`--pathogen-table`) that says which taxa are its pathogens.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Host {
    pub taxid: u32,
    pub pathogen_table: String,
}

/// Combine the mapped-length and mapped-fraction coverage checks per the
/// platform gate: an alignment passes on EITHER cutoff for long reads
/// (`or_gate` true), and only on BOTH otherwise. See
/// [`ProfilerConfig::coverage_gate_is_or`].
pub(crate) fn coverage_passes(length_ok: bool, fraction_ok: bool, or_gate: bool) -> bool {
    if or_gate {
        length_ok || fraction_ok
    } else {
        length_ok && fraction_ok
    }
}

impl ProfilerConfig {
    /// Whether the coverage gate accepts an alignment on EITHER the
    /// mapped-length or the mapped-fraction cutoff (long-read platforms:
    /// Nanopore/PacBio) instead of requiring BOTH (Illumina). Long-read lengths
    /// vary a lot, so a large absolute alignment on a very long read is
    /// confident even at low fraction, and a mostly-aligned shorter read is
    /// confident even below the length cutoff.
    pub fn coverage_gate_is_or(&self) -> bool {
        self.sequencer.is_long_read()
    }

    /// Return whether the current run should be treated as "low-map" for
    /// the chunk-breadth basis.
    ///
    /// Requires all of:
    /// 1. pathogen-detection mode is off (`host` is `None`),
    /// 2. `mapping_rate < low_map_rate_threshold` (strict `<`).
    ///
    /// Non-finite `mapping_rate` values always return `false`.
    pub fn is_low_map_run(&self, mapping_rate: f64) -> bool {
        mapping_rate.is_finite()
            && self.host.is_none()
            && mapping_rate < self.low_map_rate_threshold
    }

    pub fn mapping_load_progress_options(&self) -> ProgressOptions {
        if self.suppress_mapping_load_progress {
            ProgressOptions::disabled()
        } else {
            self.progress
        }
    }

    /// Whether `classify.txt.gz` is written: when requested, for any database
    /// except a fractional one, whose scaled counts do not match single reads.
    /// Warns when the database has marker-gene references, because the file
    /// then lists only the reads placed on those genes for their species.
    /// `has_scaled_genomes` stands in for the fractional flag when no database
    /// metadata is loaded.
    pub(crate) fn writes_classifications(&self, has_scaled_genomes: bool) -> bool {
        if !self.save_classifications {
            return false;
        }
        let database = self.db_metadata.as_deref();
        if database.map_or(has_scaled_genomes, |database| database.is_fractional) {
            log::warn!(
                target: "PROFILER",
                "Cannot write classify.txt.gz for a fractional database because per-read assignments cannot be expanded safely to full-index equivalents."
            );
            return false;
        }
        if database.is_some_and(|database| {
            database
                .references
                .iter()
                .any(|reference| reference.reference_type == ReferenceType::Markers)
        }) {
            log::warn!(
                target: "PROFILER",
                "The database has marker-gene references: for species represented by marker genes, classify.txt.gz lists only the reads placed on those genes, a small part of each genome, and not the reads from the rest of their genomes. Use profile.txt for abundances."
            );
        }
        true
    }
}

/// A configuration for unit tests: Illumina single-end reads in the default
/// mode, nothing optional; tests set the fields they need.
#[cfg(test)]
impl ProfilerConfig {
    pub(crate) fn for_tests() -> Self {
        Self {
            genome_level: false,
            mapper_mappings: None,
            mapping_profile_memory: None,
            db_metadata: None,
            sequencer: Sequencer::Illumina,
            max_em_iterations: 200,
            lenient_mapping_input: false,
            is_paired: false,
            identity: 0.95,
            identity_check_accuracy: None,
            noisy_long_reads: false,
            mapped_len: 50,
            breadth: None,
            chunk_breadth: None,
            chunk_breadth_scale: 1.0,
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            min_marker_genes: 3,
            relative_identity_gap: None,
            above_species: false,
            estimate_unclassified: false,
            identity_margin: None,
            fraction: 0.6,
            lowbiomass: false,
            save_decisions: false,
            save_classifications: false,
            by_aligned: false,
            family_identity: 0.80,
            low_map_rate_threshold: 0.30,
            host: None,
            threads: 1,
            verbose: false,
            very_verbose: false,
            progress: ProgressOptions::default(),
            suppress_mapping_load_progress: false,
            outprefix: "test".to_string(),
            mode: ProfileMode::Default,
            dmp_dir: None,
            ignore_unresolved_references: false,
        }
    }
}

/// Apply a user-selected multiplier to an automatically estimated
/// chunk-breadth cutoff. Preserve the original value bit-for-bit at the
/// default factor, and cap stricter adjustments at a valid fraction of 1.0.
#[inline]
pub(crate) fn scale_automatic_chunk_breadth(cutoff: f64, factor: f64) -> f64 {
    debug_assert!(cutoff.is_finite() && cutoff >= 0.0);
    debug_assert!(factor.is_finite() && factor > 0.0);
    if factor == 1.0 {
        cutoff
    } else {
        (cutoff * factor).min(1.0)
    }
}

pub(crate) struct ResolvedProfilerParameters<'a> {
    pub min_breadth: f64,
    pub min_chunk_breadth: f64,
    pub chunk_breadth_source: &'a str,
    pub chunk_breadth_basis_reads: usize,
    pub total_reads: usize,
    pub mapped_reads: usize,
    pub mapping_rate: f64,
}

pub(crate) fn log_resolved_profiler_parameters(
    cfg: &ProfilerConfig,
    resolved: ResolvedProfilerParameters<'_>,
) {
    log::info!(
        target: "PROFILER",
        "Minimum chunk breadth threshold: {:.6} ({}, basis_reads={}).",
        resolved.min_chunk_breadth,
        resolved.chunk_breadth_source,
        resolved.chunk_breadth_basis_reads
    );
    if !cfg.very_verbose {
        return;
    }
    log::info!(
        target: "PROFILER",
        "Resolved profiler parameters:\n\
         \tsequencer                    {}\n\
         \tmode                         {:?}\n\
         \tis_paired                    {}\n\
         \tthreads                      {}\n\
         \tmax_em_iterations            {}\n\
         \tidentity                     {:.3}\n\
         \tmapped_len                   {}\n\
         \tfraction                     {:.3}\n\
         \tmin_breadth                  {:.6}\n\
         \tmin_chunk_breadth            {:.6}\n\
         \tchunk_breadth_scale          {}\n\
         \tchunk_breadth_source         {}\n\
         \tchunk_breadth_basis_reads    {}\n\
         \tmin_reads                    {}\n\
         \tmin_oebr                     {}\n\
         \tmin_coebr                    {}\n\
         \tmin_marker_genes             {}\n\
         \trelative_identity_gap        {}\n\
         \tabove_species                {}\n\
         \testimate_unclassified        {}\n\
         \tidentity_margin              {}\n\
         \tby_mapped                    {}\n\
         \tfamily_identity              {:.3}\n\
         \tlow_map_rate_threshold       {:.3}\n\
         \tlowbiomass                   {}\n\
         \tsave_decisions               {}\n\
         \tsave_classifications         {}\n\
         \tlenient_mapping_input        {}\n\
         \tsuppress_mapping_load_progress {}\n\
         \ttotal_reads                  {}\n\
         \tmapped_reads                 {}\n\
         \tmapping_rate                 {:.5}\n\
         \tpathogen_host                {}\n\
         \thost                         {}\n\
         \tdmp_dir                      {}\n\
         \toutprefix                    {}",
        cfg.sequencer,
        cfg.mode,
        cfg.is_paired,
        cfg.threads,
        cfg.max_em_iterations,
        cfg.identity,
        cfg.mapped_len,
        cfg.fraction,
        resolved.min_breadth,
        resolved.min_chunk_breadth,
        cfg.chunk_breadth_scale,
        resolved.chunk_breadth_source,
        resolved.chunk_breadth_basis_reads,
        display_option(cfg.min_reads),
        display_option(cfg.min_oebr),
        display_option(cfg.min_coebr),
        cfg.min_marker_genes,
        display_guard(cfg.relative_identity_gap),
        cfg.above_species,
        cfg.estimate_unclassified,
        display_guard(cfg.identity_margin),
        cfg.by_aligned,
        cfg.family_identity,
        cfg.low_map_rate_threshold,
        cfg.lowbiomass,
        cfg.save_decisions,
        cfg.save_classifications,
        cfg.lenient_mapping_input,
        cfg.suppress_mapping_load_progress,
        resolved.total_reads,
        resolved.mapped_reads,
        resolved.mapping_rate,
        display_option(cfg.host.as_ref().map(|host| &host.pathogen_table)),
        display_option(cfg.host.as_ref().map(|host| host.taxid)),
        display_option(cfg.dmp_dir.as_deref()),
        cfg.outprefix,
    );
}

/// An identity guard's value, `off` when it is not set.
fn display_guard(value: Option<f64>) -> String {
    value.map_or_else(|| "off".to_string(), |value| value.to_string())
}

fn display_option<T: ToString>(value: Option<T>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "None".to_string())
}

/// The profiler of a sample: loads its mappings ([`community`]) and profiles
/// them with the level-wise engine.
pub struct Profiler {
    cfg: ProfilerConfig,
    taxonomy: Taxonomy,
}

impl Profiler {
    pub fn new(cfg: ProfilerConfig, taxonomy: Taxonomy) -> Self {
        Self { cfg, taxonomy }
    }

    /// Profile the saved or staged mapping file.
    pub fn run(mut self) -> Result<()> {
        log::info!(target: "PROFILER", "Initializing the profiler!");
        let data = self.load_mapping_file()?;
        self.profile(data)
    }

    /// Profile structured in-process Metam mappings, without writing and
    /// parsing a mapping file.
    pub fn run_with_metam_mappings(
        mut self,
        receiver: MetamMappingReceiver,
        source: &str,
    ) -> Result<()> {
        log::info!(target: "PROFILER", "Initializing the profiler!");
        let data = self.load_mapping_stream(receiver, source)?;
        self.profile(data)
    }
}

#[cfg(test)]
mod reference_taxonomy_tests {
    use super::{
        parse_reference_taxonomy, taxonomy_adjustment_description, ProfileRank, ReferenceTaxonomy,
        ReferenceTaxonomyResolver, UnresolvedReferences,
    };
    use crate::taxonomy::Taxonomy;
    use std::collections::{HashMap, HashSet};

    fn current_taxonomy() -> Taxonomy {
        Taxonomy::from_parts(
            HashMap::from([
                (1, (1, "no rank".to_string())),
                (10, (1, "genus".to_string())),
                (20, (10, "species".to_string())),
                (30, (20, "subspecies".to_string())),
                (40, (1, "genus".to_string())),
                (50, (20, "no rank".to_string())),
            ]),
            HashMap::new(),
            HashMap::from([(29, 30), (39, 40), (49, 50)]),
            HashSet::from([48, 999]),
        )
    }

    #[test]
    fn taxonomy_adjustment_is_human_readable() {
        let declared = ReferenceTaxonomy {
            primary_taxid: Some(572_307),
            profile_taxid: 29_176,
            profile_rank: ProfileRank::Species,
        };
        let resolved = ReferenceTaxonomy {
            primary_taxid: Some(29_176),
            ..declared
        };
        assert_eq!(
            taxonomy_adjustment_description(declared, resolved),
            "primary 572307 → 29176; terminal remains species:29176"
        );
    }

    #[test]
    fn taxonomy_adjustment_formats_missing_and_changed_endpoints() {
        let declared = ReferenceTaxonomy {
            primary_taxid: None,
            profile_taxid: 100,
            profile_rank: ProfileRank::Species,
        };
        let resolved = ReferenceTaxonomy {
            primary_taxid: Some(200),
            profile_taxid: 20,
            profile_rank: ProfileRank::Genus,
        };
        assert_eq!(
            taxonomy_adjustment_description(declared, resolved),
            "primary none → 200; terminal species:100 → genus:20"
        );
    }

    #[test]
    fn preserves_the_most_specific_supported_terminal_rank() {
        let strain = parse_reference_taxonomy("t:101;s:100").unwrap();
        assert_eq!(strain.primary_taxid, Some(101));
        assert_eq!(strain.profile_taxid, 100);
        assert_eq!(strain.profile_rank, ProfileRank::Species);

        for (field, taxid, rank) in [
            ("g:300", 300, ProfileRank::Genus),
            ("f:400", 400, ProfileRank::Family),
            ("o:500", 500, ProfileRank::Order),
            ("c:600", 600, ProfileRank::Class),
            ("p:700", 700, ProfileRank::Phylum),
            ("k:2", 2, ProfileRank::Kingdom),
        ] {
            let taxonomy = parse_reference_taxonomy(field).unwrap();
            assert_eq!(taxonomy.primary_taxid, None);
            assert_eq!(taxonomy.profile_taxid, taxid);
            assert_eq!(taxonomy.profile_rank, rank);
        }
    }

    #[test]
    fn species_remains_the_terminal_rank_when_higher_ancestors_are_present() {
        let taxonomy = parse_reference_taxonomy("s:100;g:10;f:5").unwrap();
        assert_eq!(taxonomy.profile_taxid, 100);
        assert_eq!(taxonomy.profile_rank, ProfileRank::Species);
    }

    #[test]
    fn community_endpoint_skips_strain_and_subspecies_taxids() {
        let declared = parse_reference_taxonomy("t:29;s:20").unwrap();
        let taxonomy = current_taxonomy();
        let resolved = ReferenceTaxonomyResolver::new(&taxonomy)
            .resolve(declared, "strain-reference")
            .unwrap();
        assert_eq!(resolved.primary_taxid, Some(30));
        assert_eq!(resolved.profile_taxid, 20);
        assert_eq!(resolved.profile_rank, ProfileRank::Species);
    }

    #[test]
    fn stale_species_prefix_at_a_subspecies_node_uses_species_fallback() {
        let declared = parse_reference_taxonomy("s:30").unwrap();
        let taxonomy = current_taxonomy();
        let resolved = ReferenceTaxonomyResolver::new(&taxonomy)
            .resolve(declared, "reclassified-reference")
            .unwrap();
        assert_eq!(resolved.primary_taxid, Some(30));
        assert_eq!(resolved.profile_taxid, 20);
        assert_eq!(resolved.profile_rank, ProfileRank::Species);
    }

    #[test]
    fn current_dump_rank_replaces_a_stale_header_rank() {
        let declared = parse_reference_taxonomy("s:39").unwrap();
        let taxonomy = current_taxonomy();
        let resolved = ReferenceTaxonomyResolver::new(&taxonomy)
            .resolve(declared, "merged-reference")
            .unwrap();
        assert_eq!(resolved.primary_taxid, None);
        assert_eq!(resolved.profile_taxid, 40);
        assert_eq!(resolved.profile_rank, ProfileRank::Genus);
    }

    #[test]
    fn unresolved_references_are_skipped_only_when_allowed() {
        let mut unresolved = UnresolvedReferences::default();
        assert!(unresolved
            .skip("G1", anyhow::anyhow!("cannot be resolved"), false)
            .is_err());
        assert_eq!(unresolved.references, 0);
        for genome in ["G1", "G1", "G2"] {
            unresolved
                .skip(genome, anyhow::anyhow!("cannot be resolved"), true)
                .unwrap();
        }
        assert_eq!(unresolved.references, 3);
        assert_eq!(unresolved.genomes.len(), 2);
    }

    #[test]
    fn unresolved_header_taxid_is_a_clear_error() {
        let declared = parse_reference_taxonomy("s:999").unwrap();
        let taxonomy = current_taxonomy();
        let error = ReferenceTaxonomyResolver::new(&taxonomy)
            .resolve(declared, "missing-reference")
            .unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("missing-reference"));
        assert!(message.contains("cannot be resolved"));
    }

    #[test]
    fn current_primary_lineage_recovers_a_deleted_species_taxid() {
        let declared = parse_reference_taxonomy("t:49;s:999").unwrap();
        let taxonomy = current_taxonomy();
        let resolved = ReferenceTaxonomyResolver::new(&taxonomy)
            .resolve(declared, "recoverable-reference")
            .unwrap();
        assert_eq!(resolved.primary_taxid, Some(50));
        assert_eq!(resolved.profile_taxid, 20);
        assert_eq!(resolved.profile_rank, ProfileRank::Species);
    }

    #[test]
    fn deleted_primary_falls_back_to_a_current_species() {
        let declared = parse_reference_taxonomy("t:48;s:20").unwrap();
        let taxonomy = current_taxonomy();
        let resolved = ReferenceTaxonomyResolver::new(&taxonomy)
            .resolve(declared, "species-only-reference")
            .unwrap();
        assert_eq!(resolved.primary_taxid, None);
        assert_eq!(resolved.profile_taxid, 20);
        assert_eq!(resolved.profile_rank, ProfileRank::Species);
    }
}

#[cfg(test)]
mod subsample_scaling_tests {
    use super::{
        push_string_keyed_value, record_genome_sampling_scale, sampling_scale_from_sizes,
        strip_pair_suffix,
    };
    use std::collections::HashMap;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-12,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn scaling_from_sizes_uses_typed_full_and_effective_lengths() {
        assert_close(sampling_scale_from_sizes(1_000, 250), 4.0);
    }

    #[test]
    fn scaling_from_sizes_ignores_full_indexes() {
        assert_close(sampling_scale_from_sizes(1_000, 1_000), 1.0);
        assert_close(sampling_scale_from_sizes(250, 1_000), 1.0);
    }

    #[test]
    fn scaling_is_independent_per_genome_and_rejects_inconsistency() {
        let mut scaling = HashMap::new();
        record_genome_sampling_scale(&mut scaling, "full", 1.0).unwrap();
        record_genome_sampling_scale(&mut scaling, "sampled", 4.0).unwrap();
        record_genome_sampling_scale(&mut scaling, "sampled", 4.0).unwrap();

        assert_close(*scaling.get("full").unwrap(), 1.0);
        assert_close(*scaling.get("sampled").unwrap(), 4.0);
        assert!(record_genome_sampling_scale(&mut scaling, "sampled", 2.0).is_err());
    }

    #[test]
    fn borrowed_mapping_keys_are_owned_only_once() {
        let mut values = HashMap::new();
        push_string_keyed_value(&mut values, "read", 1usize);
        push_string_keyed_value(&mut values, "read", 2usize);

        assert_eq!(values.len(), 1);
        assert_eq!(values.get("read"), Some(&vec![1, 2]));
        assert_eq!(strip_pair_suffix("read/1"), "read");
        assert_eq!(strip_pair_suffix("read/2"), "read");
        assert_eq!(strip_pair_suffix("read/x"), "read/x");
    }
}

#[cfg(test)]
mod chunk_basis_tests {
    use std::path::PathBuf;
    use std::sync::Arc;

    use crate::db::{DatabaseMetadata, ReferenceMetadata, ReferenceType};
    use crate::sequencer::Sequencer;

    use super::{scale_automatic_chunk_breadth, Host, ProfilerConfig};

    #[test]
    fn automatic_chunk_scale_preserves_default_and_caps_at_one() {
        let cutoff = 0.123456789;
        assert_eq!(
            scale_automatic_chunk_breadth(cutoff, 1.0).to_bits(),
            cutoff.to_bits()
        );
        assert_eq!(scale_automatic_chunk_breadth(cutoff, 0.5), cutoff * 0.5);
        assert_eq!(scale_automatic_chunk_breadth(0.75, 2.0), 1.0);
    }

    #[test]
    fn a_low_map_run_maps_below_the_threshold_without_a_host() {
        let mut cfg = ProfilerConfig::for_tests();
        // Strictly below the threshold, on every read platform, single-end or
        // paired.
        for sequencer in [Sequencer::Illumina, Sequencer::Nanopore, Sequencer::PacBio] {
            for paired in [false, true] {
                cfg.sequencer = sequencer;
                cfg.is_paired = paired;
                assert!(cfg.is_low_map_run(0.0), "{sequencer}");
                assert!(cfg.is_low_map_run(0.2999999), "{sequencer}");
                assert!(!cfg.is_low_map_run(0.30), "{sequencer}");
                assert!(!cfg.is_low_map_run(0.99), "{sequencer}");
            }
        }
        assert!(!cfg.is_low_map_run(f64::NAN));
        // Another threshold.
        cfg.low_map_rate_threshold = 0.50;
        assert!(cfg.is_low_map_run(0.4999));
        assert!(!cfg.is_low_map_run(0.50));
        // A host switches it off.
        cfg.host = Some(Host {
            taxid: 9606,
            pathogen_table: "pathogens.tsv".to_string(),
        });
        assert!(!cfg.is_low_map_run(0.0));
    }

    #[test]
    fn classifications_are_written_for_every_database_except_fractional_ones() {
        let reference = |reference_type| ReferenceMetadata {
            name: "ref".to_string(),
            length: 1_000,
            genome: "genome".into(),
            contig: "ref".to_string(),
            full_genome_size: 1_000,
            effective_indexed_size: 1_000,
            taxonomy_field: "s:8".into(),
            reference_type,
            declared_genome_size: None,
        };
        let database = |types: &[ReferenceType], is_fractional| {
            Some(Arc::new(DatabaseMetadata {
                metam_index: PathBuf::from("metax_db.mxi"),
                references: types.iter().copied().map(reference).collect(),
                is_fractional,
            }))
        };
        let mut cfg = ProfilerConfig::for_tests();
        cfg.db_metadata = database(&[ReferenceType::Genome], false);
        assert!(!cfg.writes_classifications(false));

        cfg.save_classifications = true;
        assert!(cfg.writes_classifications(false));
        cfg.db_metadata = database(&[ReferenceType::Markers], false);
        assert!(cfg.writes_classifications(false));
        cfg.db_metadata = database(&[ReferenceType::Genome, ReferenceType::Markers], false);
        assert!(cfg.writes_classifications(false));
        cfg.db_metadata = database(&[ReferenceType::Genome], true);
        assert!(!cfg.writes_classifications(false));

        cfg.db_metadata = None;
        assert!(cfg.writes_classifications(false));
        assert!(!cfg.writes_classifications(true));
    }
}

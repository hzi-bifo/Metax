//! Profiler Module - Taxonomic abundance estimation from mapping data
//!
//! This module contains the core profiling logic for Metax. It provides two
//! profiling strategies:
//!
//! - **Community Profiler** (`community.rs`): Terminal-rank profiling that
//!   reports each reference at its lowest confidently known taxonomy rank.
//!   Best for typical metagenomic samples.
//!
//! - **Strain Profiler** (`strain.rs`): Strain-level profiling that can
//!   distinguish between closely related strains within a species. More
//!   computationally intensive but provides finer-grained resolution.
//!   Experimental: it runs its own filter chain with its own presets rather
//!   than the community evidence policy.
//!
//! # Profiling Pipeline
//!
//! Both profilers follow a similar pipeline:
//!
//! 1. **Load mappings**: Parse Metam mappings or SAM/BAM records and extract reference metadata
//! 2. **Compute coverage**: Calculate breadth and depth for each reference
//! 3. **Assign taxonomy**: Map read mappings to taxonomic units
//! 4. **Resolve ambiguity**: Use EM algorithm for multi-mapping reads
//! 5. **Filter and output**: Apply coverage filters and write results
//!
//! # Coverage Metrics
//!
//! The profiler computes several coverage metrics:
//!
//! - **Breadth**: Fraction of genome covered by at least one read
//! - **Depth**: Average number of reads covering each position
//! - **Chunk breadth**: Coverage of fixed-size genome chunks (for QC)
//! - **Flex chunk breadth**: Coverage of variable-size chunks (based on √genome_size)

use std::collections::{BTreeMap, HashMap};
use std::hash::BuildHasher;
use std::path::Path;
use std::sync::Arc;
use std::thread::JoinHandle;

use anyhow::{bail, Context, Result};
use clap::ValueEnum;
use xxhash_rust::xxh3::Xxh3DefaultBuilder;

use crate::db::DatabaseMetadata;
use crate::mapper::tsv::{MetamMappingReceiver, SharedMetamProfile};
use crate::progress::ProgressOptions;
use crate::sequencer::Sequencer;
use crate::taxonomy::Taxonomy;

// Sub-modules for different profiling strategies
mod alignment;
pub(crate) mod community;
pub(crate) mod evidence;
mod output;
mod pathogen;
mod statistics;
pub(crate) mod strain;

/// Read-name aggregation map with a stable, high-throughput hash builder.
///
/// Native mapping records arrive in a repeatable stream, but `HashMap`'s
/// per-process random seed made the later `into_iter()` order vary and changed
/// floating-point reduction order. XXH3 is already a project dependency and
/// is faster than the standard keyed hasher for short sequencing identifiers.
pub(crate) type StableReadMap<T> = HashMap<String, Vec<T>, Xxh3DefaultBuilder>;

/// Parse one SAM `@SQ` dictionary record and require the fields needed for
/// coordinate normalization. The compatibility profiler cannot safely infer a
/// missing reference name or length from alignment rows.
fn parse_sam_sequence_record(line: &str) -> Result<(String, usize)> {
    let mut name = None;
    let mut length = None;
    for field in line.split('\t').skip(1) {
        if let Some(value) = field.strip_prefix("SN:") {
            name = Some(value);
        } else if let Some(value) = field.strip_prefix("LN:") {
            length = Some(value);
        }
    }
    let name = name
        .filter(|value| !value.is_empty())
        .ok_or_else(|| anyhow::anyhow!("SAM @SQ record is missing a non-empty SN field: {line}"))?;
    let raw_length = length
        .filter(|value| !value.is_empty())
        .ok_or_else(|| anyhow::anyhow!("SAM @SQ record is missing an LN field: {line}"))?;
    let length = raw_length
        .parse::<usize>()
        .with_context(|| format!("invalid SAM @SQ length '{raw_length}' for reference '{name}'"))?;
    if length == 0 {
        bail!("SAM @SQ reference '{}' has zero length", name);
    }
    Ok((name.to_string(), length))
}

/// Require a usable sequence dictionary before parsing SAM/BAM mappings.
/// Metax needs reference lengths for every coverage statistic, so accepting a
/// headerless alignment file would silently discard all mapped records.
fn require_alignment_reference_dictionary(path: &Path, references: usize) -> Result<()> {
    if references == 0 {
        bail!(
            "alignment file '{}' has no usable reference sequence dictionary; include @SQ records with SN and positive LN fields",
            path.display()
        );
    }
    Ok(())
}

/// Join a SAM parsing worker without discarding its panic. Worker panics are
/// exceptional, so extracting the payload has no cost on successful runs.
fn join_alignment_worker(handle: JoinHandle<()>, label: &str) -> Result<()> {
    handle.join().map_err(|payload| {
        let message = if let Some(message) = payload.downcast_ref::<&str>() {
            (*message).to_string()
        } else if let Some(message) = payload.downcast_ref::<String>() {
            message.clone()
        } else {
            "unknown panic payload".to_string()
        };
        anyhow::anyhow!("{label} panicked: {message}")
    })
}

#[cfg(test)]
mod sam_header_tests {
    use std::path::Path;

    use super::{
        join_alignment_worker, parse_sam_sequence_record, require_alignment_reference_dictionary,
    };

    #[test]
    fn sam_sequence_record_requires_name_and_positive_length() {
        assert_eq!(
            parse_sam_sequence_record("@SQ\tSN:ref-a\tLN:123").unwrap(),
            ("ref-a".to_string(), 123)
        );
        for line in [
            "@SQ\tLN:123",
            "@SQ\tSN:ref-a",
            "@SQ\tSN:ref-a\tLN:0",
            "@SQ\tSN:ref-a\tLN:not-a-number",
        ] {
            assert!(parse_sam_sequence_record(line).is_err(), "line={line}");
        }
    }

    #[test]
    fn alignment_input_requires_a_reference_dictionary() {
        let error = require_alignment_reference_dictionary(Path::new("headerless.sam"), 0)
            .expect_err("a missing sequence dictionary must fail");
        let message = format!("{error:#}");
        assert!(message.contains("headerless.sam"));
        assert!(message.contains("no usable reference sequence dictionary"));
        require_alignment_reference_dictionary(Path::new("valid.sam"), 1).unwrap();
    }

    #[test]
    fn alignment_worker_panic_is_reported() {
        let handle = std::thread::spawn(|| panic!("synthetic worker failure"));
        let error = join_alignment_worker(handle, "SAM parsing worker")
            .expect_err("worker panic must be propagated");
        let message = format!("{error:#}");
        assert!(message.contains("SAM parsing worker panicked"));
        assert!(message.contains("synthetic worker failure"));
    }
}

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

    #[inline]
    pub(crate) fn is_species(self) -> bool {
        self == Self::Species
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

/// Read a genome-specific sampling scale from the internal SAM reference-name
/// layout `...|full_genome_size|effective_indexed_size`. Missing or invalid
/// size fields are treated as an unscaled, fully indexed reference.
pub(crate) fn sampling_scale_from_name_parts(parts: &[&str]) -> f64 {
    if parts.len() < 6 {
        return 1.0;
    }
    let full = match parts[4].parse::<usize>() {
        Ok(value) if value > 0 => value,
        _ => return 1.0,
    };
    let effective = match parts[5].parse::<usize>() {
        Ok(value) if value > 0 => value,
        _ => return 1.0,
    };
    sampling_scale_from_sizes(full, effective)
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

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct ReadSupport {
    pub observed: f64,
    pub equivalent: f64,
}

impl ReadSupport {
    pub fn add_weighted(&mut self, observed: f64, equivalent_sum: f64, weight: f64) {
        self.observed += observed * weight;
        self.equivalent += equivalent_sum * weight;
    }
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

/// Policy for exact, equal-score cross-species ties in short-read evidence.
///
/// `GuardedMap` uses a seeded EM prior only for exact ties, requires a
/// uniquely dominant posterior, and cannot create the independent coverage
/// needed for a new species call. It is the default for native single-end and
/// paired-end Illumina community profiling.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, ValueEnum)]
pub enum ExactTiePolicy {
    /// Keep exact cross-species ties ambiguous for coordinate evidence.
    Ambiguous,
    /// Resolve only strongly supported exact ties with the guarded MAP rule.
    #[default]
    GuardedMap,
}

impl std::fmt::Display for ExactTiePolicy {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Ambiguous => formatter.write_str("ambiguous"),
            Self::GuardedMap => formatter.write_str("guarded-map"),
        }
    }
}

/// Configuration for the taxonomy profiler.
///
/// This structure holds all parameters that control profiling behavior,
/// from input/output paths to filtering thresholds to algorithm options.
/// Values here are the resolved defaults used by the profilers.
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
    /// Path to the input SAM/BAM mapping file
    pub sam: String,

    /// Path to Metam coordinate TSV output for native Metam-backed profiling.
    /// When set, the profiler reads coordinates directly instead of parsing SAM.
    pub mapper_mappings: Option<String>,

    /// In-memory mapping counters for the default streamed mapper path.
    pub mapping_profile_memory: Option<SharedMetamProfile>,

    /// Metax database metadata used to interpret native mapping reference IDs.
    /// File-backed mapping output can provide this metadata directly.
    pub db_metadata: Option<Arc<DatabaseMetadata>>,

    /// Sequencing platform (Illumina, Nanopore, PacBio)
    /// Affects default threshold selection
    pub sequencer: Sequencer,

    /// Number of reads to process per batch
    /// Larger batches improve throughput but use more memory
    pub batch_size: usize,

    /// Maximum EM iterations for ambiguous-read redistribution.
    pub max_em_iterations: usize,

    /// Skip malformed Metam mapping rows instead of aborting the run.
    pub lenient_mapping_input: bool,

    /// Whether input is paired-end Illumina data
    /// Affects read name parsing and coverage calculation
    pub is_paired: bool,

    /// How exact cross-species short-read score ties contribute coordinate
    /// evidence.
    pub exact_ties: ExactTiePolicy,

    /// Minimum mapping identity threshold (0.0-1.0). Native mappings use the
    /// mapper's canonical per-base identity; legacy SAM/BAM retains its
    /// gap-compressed CIGAR identity.
    pub identity: f64,

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

    /// Minimum fraction of read that must map (0.0-1.0)
    pub fraction: f64,

    /// Enable low-biomass sample heuristics
    /// Disables chunk breadth threshold
    pub lowbiomass: bool,

    /// Keep raw (unfiltered) profile output
    pub keep_raw: bool,

    /// Write optional per-read taxonomic assignments.
    pub save_classifications: bool,

    /// Force mapped-read basis for the chunk-breadth estimate.
    pub by_aligned: bool,

    /// Minimum identity accepted for the genus-level fallback path.
    /// Comparison is inclusive (`identity >= genus_identity`). Default 0.80.
    pub genus_identity: f64,

    /// Mapping-rate cut-off for the low-map path. A run is considered
    /// low-map (genus-fallback gate active, auto `--by-mapped` forced on)
    /// when `mapped_reads / total_reads < low_map_rate_threshold`. The
    /// comparison is strict. Default 0.30.
    pub low_map_rate_threshold: f64,

    /// Enable genus-level assignment for mappings that miss species cutoffs.
    pub genus_fallback: bool,

    /// Require both p-value columns to meet the cutoff when both are available.
    pub strict: bool,

    /// Path to pathogen-host mapping TSV (optional)
    pub pathogen_host: Option<String>,

    /// Host taxonomy ID for pathogen filtering (optional)
    pub host: Option<String>,

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

    /// Enable strain-level profiling
    pub strain: bool,
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
    /// the chunk-breadth basis / genus-fallback gate.
    ///
    /// Requires all of:
    /// 1. sequencer is a read platform — Illumina, Nanopore, or PacBio
    ///    (case-insensitive). Assembly / unknown input is excluded because a
    ///    contig "mapping rate" is not comparable to a per-read rate.
    /// 2. pathogen-detection mode is off (`host` is `None`),
    /// 3. `mapping_rate < low_map_rate_threshold` (strict `<`).
    ///
    /// Non-finite `mapping_rate` values always return `false`.
    pub fn is_low_map_run(&self, mapping_rate: f64) -> bool {
        let read_platform = self.sequencer.is_read_platform();
        mapping_rate.is_finite()
            && read_platform
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

    /// Resolve the effective mapped-vs-total basis for chunk-breadth
    /// estimation.
    ///
    /// - If the user passed `--by-mapped`, always use mapped basis.
    /// - Otherwise auto-force mapped basis for low-map read-platform
    ///   (Illumina/Nanopore/PacBio) non-pathogen runs (see
    ///   [`Self::is_low_map_run`]).
    /// - Fall back to total-read basis in all other cases.
    pub fn resolve_by_aligned(&self, mapping_rate: f64) -> bool {
        self.by_aligned || self.is_low_map_run(mapping_rate)
    }

    /// Read count passed into community chunk-breadth estimation when
    /// `--min-chunk-breadth` is unset.
    ///
    /// Uses mapped vs total reads according to [`Self::resolve_by_aligned`].
    /// Community-specific fallback/scaling rules that depend on reference DB
    /// layout are applied by the community profiler, not by this generic
    /// config helper.
    pub fn chunk_breadth_basis_reads(
        &self,
        mapping_rate: f64,
        num_reads: usize,
        aligned_reads: usize,
    ) -> usize {
        if self.resolve_by_aligned(mapping_rate) {
            aligned_reads
        } else {
            num_reads
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
         \texact_ties                  {}\n\
         \tstrain                       {}\n\
         \tthreads                      {}\n\
         \tbatch_size                   {}\n\
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
         \tby_mapped                    {}\n\
         \tgenus_fallback               {}\n\
         \tgenus_identity               {:.3}\n\
         \tlow_map_rate_threshold       {:.3}\n\
         \tlowbiomass                   {}\n\
         \tstrict                       {}\n\
         \tkeep_raw                     {}\n\
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
        cfg.exact_ties,
        cfg.strain,
        cfg.threads,
        cfg.batch_size,
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
        cfg.by_aligned,
        cfg.genus_fallback,
        cfg.genus_identity,
        cfg.low_map_rate_threshold,
        cfg.lowbiomass,
        cfg.strict,
        cfg.keep_raw,
        cfg.save_classifications,
        cfg.lenient_mapping_input,
        cfg.suppress_mapping_load_progress,
        resolved.total_reads,
        resolved.mapped_reads,
        resolved.mapping_rate,
        display_option(cfg.pathogen_host.as_deref()),
        display_option(cfg.host.as_deref()),
        display_option(cfg.dmp_dir.as_deref()),
        cfg.outprefix,
    );
}

fn display_option<T: ToString>(value: Option<T>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "None".to_string())
}

/// Main profiler that dispatches to the appropriate profiling strategy.
///
/// This is a facade that delegates to either `CommunityProfiler` or
/// `StrainProfiler` based on the `strain` configuration option.
///
/// # Usage
///
/// ```ignore
/// let profiler = Profiler::new(config, taxonomy)?;
/// profiler.run()?;
/// ```
pub struct Profiler {
    /// Configuration for the profiling run
    cfg: ProfilerConfig,

    /// Taxonomy database (shared via Arc for thread safety)
    taxonomy: Arc<Taxonomy>,
}

impl Profiler {
    /// Create a new profiler with the given configuration and taxonomy.
    ///
    /// The taxonomy is wrapped in Arc for efficient sharing across threads.
    ///
    /// # Arguments
    ///
    /// * `cfg` - Profiler configuration
    /// * `taxonomy` - NCBI taxonomy database
    ///
    /// # Returns
    ///
    /// A new Profiler instance ready to run.
    pub fn new(cfg: ProfilerConfig, taxonomy: Taxonomy) -> Result<Self> {
        Ok(Self {
            cfg,
            taxonomy: Arc::new(taxonomy),
        })
    }

    /// Run the profiling pipeline.
    ///
    /// Dispatches to either StrainProfiler or CommunityProfiler based on
    /// the `strain` configuration option.
    ///
    /// # Returns
    ///
    /// Ok(()) on success, or an error if profiling fails.
    pub fn run(&self) -> Result<()> {
        if self.cfg.strain {
            // Strain-level profiling for higher resolution
            strain::StrainProfiler::new(self.cfg.clone(), self.taxonomy.clone())?.run()
        } else {
            // Species-level community profiling (default)
            community::CommunityProfiler::new(self.cfg.clone(), self.taxonomy.clone())?.run()
        }
    }

    /// Run profiling from structured in-process Metam mappings, avoiding
    /// mapping TSV serialization and parsing.
    pub fn run_with_metam_mappings(
        &self,
        receiver: MetamMappingReceiver,
        source: &str,
    ) -> Result<()> {
        if self.cfg.strain {
            strain::StrainProfiler::new(self.cfg.clone(), self.taxonomy.clone())?
                .run_with_metam_mappings(receiver, source)
        } else {
            community::CommunityProfiler::new(self.cfg.clone(), self.taxonomy.clone())?
                .run_with_metam_mappings(receiver, source)
        }
    }
}

#[cfg(test)]
mod reference_taxonomy_tests {
    use super::{
        parse_reference_taxonomy, taxonomy_adjustment_description, ProfileRank, ReferenceTaxonomy,
        ReferenceTaxonomyResolver,
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
        equivalent_read_count_passes_minimum, push_string_keyed_value,
        record_genome_sampling_scale, sampling_scale_from_name_parts, sampling_scale_from_sizes,
        strip_pair_suffix, ReadSupport,
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
    fn scaling_from_name_parts_preserves_legacy_sam_layout() {
        let parts = ["asm", "101", "100", "contig", "1000", "250"];
        assert_close(sampling_scale_from_name_parts(&parts), 4.0);
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
    fn mixed_full_and_fractional_support_is_scaled_before_taxon_aggregation() {
        let mut support = ReadSupport::default();
        support.add_weighted(1.0, 1.0, 1.0);
        support.add_weighted(1.0, 10.0, 1.0);

        assert_close(support.observed, 2.0);
        assert_close(support.equivalent, 11.0);
        assert!(equivalent_read_count_passes_minimum(
            support.equivalent,
            Some(10)
        ));
        assert!(!equivalent_read_count_passes_minimum(
            support.equivalent,
            Some(12)
        ));
    }

    #[test]
    fn em_weight_is_applied_to_observed_and_equivalent_support() {
        let mut support = ReadSupport::default();
        // Two ambiguous reads: one from a full genome (1x) and one from a
        // 10%-indexed genome (10x), assigned to this taxon at 25%.
        support.add_weighted(2.0, 11.0, 0.25);

        assert_close(support.observed, 0.5);
        assert_close(support.equivalent, 2.75);
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
    use crate::progress::ProgressOptions;
    use crate::sequencer::Sequencer;

    use super::{scale_automatic_chunk_breadth, ExactTiePolicy, ProfileMode, ProfilerConfig};

    fn minimal_cfg() -> ProfilerConfig {
        ProfilerConfig {
            sam: String::new(),
            mapper_mappings: None,
            mapping_profile_memory: None,
            db_metadata: None,
            sequencer: Sequencer::Illumina,
            batch_size: 1,
            max_em_iterations: 200,
            lenient_mapping_input: false,
            is_paired: false,
            exact_ties: ExactTiePolicy::GuardedMap,
            identity: 0.97,
            mapped_len: 50,
            breadth: None,
            chunk_breadth: None,
            chunk_breadth_scale: 1.0,
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            fraction: 0.7,
            lowbiomass: false,
            keep_raw: false,
            save_classifications: false,
            by_aligned: false,
            genus_identity: 0.80,
            low_map_rate_threshold: 0.30,
            genus_fallback: false,
            strict: false,
            pathogen_host: None,
            host: None,
            threads: 1,
            verbose: false,
            very_verbose: false,
            progress: ProgressOptions::default(),
            suppress_mapping_load_progress: false,
            outprefix: "t".to_string(),
            mode: ProfileMode::Default,
            dmp_dir: None,
            strain: false,
        }
    }

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
    fn chunk_basis_genus_fallback_non_low_map_uses_total_reads_without_by_aligned() {
        let mut cfg = minimal_cfg();
        cfg.genus_fallback = true;
        // High mapping rate => not low-map; no --by-mapped flag.
        let basis = cfg.chunk_breadth_basis_reads(0.99, 1_000_000, 500_000);
        assert_eq!(
            basis, 1_000_000,
            "explicit genus-fallback alone must not force mapped-read basis on non-low-map runs"
        );
    }

    #[test]
    fn chunk_basis_genus_fallback_non_low_map_honours_by_aligned() {
        let mut cfg = minimal_cfg();
        cfg.genus_fallback = true;
        cfg.by_aligned = true;
        let basis = cfg.chunk_breadth_basis_reads(0.99, 1_000_000, 500_000);
        assert_eq!(
            basis, 500_000,
            "--genus-fallback must not override an explicit --by-mapped basis"
        );
    }

    #[test]
    fn chunk_basis_low_map_unaffected_by_genus_fallback_flag_for_basis_rule() {
        let mut cfg = minimal_cfg();
        cfg.genus_fallback = true;
        // Low-map Illumina: auto by-mapped => mapped basis unless user overrode.
        let basis = cfg.chunk_breadth_basis_reads(0.10, 1_000_000, 100_000);
        assert_eq!(basis, 100_000);
    }

    #[test]
    fn low_map_run_covers_read_platforms_but_not_assembly() {
        let mut cfg = minimal_cfg();
        // Every read platform honours the low-map auto-trigger, so genus
        // fallback and the mapped-read basis reach long reads without an
        // explicit flag. Assembly input is excluded.
        for seq in [Sequencer::Illumina, Sequencer::Nanopore, Sequencer::PacBio] {
            cfg.sequencer = seq;
            assert!(cfg.is_low_map_run(0.10), "{seq}: low rate must be low-map");
            assert!(!cfg.is_low_map_run(0.99), "{seq}: high rate is not low-map");
            // Low-map runs flip the chunk-breadth basis to aligned reads.
            assert_eq!(
                cfg.chunk_breadth_basis_reads(0.10, 1_000_000, 100_000),
                100_000,
                "{seq}"
            );
        }
        // Assembly input is not a read platform -> never auto low-map.
        cfg.sequencer = Sequencer::Assembly;
        assert!(!cfg.is_low_map_run(0.10));
        // Pathogen-detection mode (host set) suppresses it on any platform.
        cfg.sequencer = Sequencer::Nanopore;
        cfg.host = Some("human".to_string());
        assert!(!cfg.is_low_map_run(0.10));
    }
}

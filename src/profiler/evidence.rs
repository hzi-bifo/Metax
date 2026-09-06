//! Deterministic coordinate-evidence resolution for native Metam mappings.
//!
//! Presence and abundance deliberately use different views of a fragment:
//! every valid terminal-taxon candidate remains available to abundance EM,
//! while only a clearly separated best endpoint, or an exact short-read tie
//! that clears the guarded-MAP posterior threshold, may contribute coordinate
//! coverage. This keeps unresolved conserved mappings ambiguous instead of
//! allowing them to manufacture breadth.

use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::sync::Arc;

use metam::format::canonical_mapping_identity;
use metam::mapper::CandidateBoundary;

use super::statistics::{calc_breadth_pvalue, Tail};
use super::{ExactTiePolicy, ProfileRank, ProfilerConfig};

const SCORE_EPSILON: f64 = 1e-9;
const FRACTION_F32_EPSILON: f64 = f32::EPSILON as f64;
pub(crate) const PAIRED_ERROR_MARGIN: f64 = 1.0;
pub(crate) const SINGLE_ERROR_MARGIN: f64 = 0.5;
pub(crate) const LONG_RELATIVE_SCORE_MARGIN: f64 = 0.01;
pub(crate) const GUARDED_MAP_MIN_POSTERIOR: f64 = 0.95;

/// Role carried explicitly by one native mapping row.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MappingRole {
    /// A single-end or long-read mapping.
    Unpaired,
    /// R1 in an explicitly identified paired placement.
    Mate1,
    /// R2 in an explicitly identified paired placement.
    Mate2,
}

/// Compact mapping representation shared by legacy and native community paths.
#[derive(Clone, Debug)]
pub(crate) struct EvidenceAlignment {
    pub(crate) asm: Arc<str>,
    pub(crate) acc: Arc<str>,
    pub(crate) profile_taxid: u32,
    /// Most-specific rank supported by this reference's taxonomy metadata.
    pub(crate) profile_rank: ProfileRank,
    pub(crate) reference_id: u32,
    pub(crate) assembly_id: u32,
    pub(crate) role: MappingRole,
    /// Read-local mapper-selected pair identifier; zero for unpaired rows.
    pub(crate) pair_id: u32,
    /// Taxonomic candidate levels that may be incomplete for this read/pair.
    pub(crate) candidate_boundary: CandidateBoundary,
    pub(crate) sampling_scale: f32,
    pub(crate) asm_len: f64,
    pub(crate) identity: f32,
    pub(crate) fraction: f32,
    pub(crate) mapped_len: u32,
    pub(crate) span: (u32, u32),
    pub(crate) fixed_chunk: (u32, u32),
    pub(crate) flex_chunk: (u32, u32),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum EvidenceMode {
    Paired,
    Single,
    Long,
}

impl EvidenceMode {
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Paired => "paired margin-1",
            Self::Single => "single-end margin-0.5",
            Self::Long => "long-read relative-margin-0.01",
        }
    }

    #[inline]
    pub(crate) fn supports_guarded_map(self) -> bool {
        matches!(self, Self::Paired | Self::Single)
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct EvidencePolicy {
    pub(crate) mode: EvidenceMode,
    pub(crate) identity: f64,
    pub(crate) mapped_len: usize,
    pub(crate) fraction: f64,
    pub(crate) collect_exact_ties: bool,
}

impl EvidencePolicy {
    pub(crate) fn from_config(cfg: &ProfilerConfig) -> Option<Self> {
        // Strain profiling has its own evidence implementation; never
        // activate the community policy for a strain configuration.
        if cfg.strain {
            return None;
        }
        let mode = if cfg.is_paired && cfg.sequencer.is_illumina() {
            EvidenceMode::Paired
        } else if !cfg.is_paired && cfg.sequencer.is_illumina() {
            EvidenceMode::Single
        } else if !cfg.is_paired && cfg.sequencer.is_long_read() {
            EvidenceMode::Long
        } else {
            return None;
        };
        Some(Self {
            mode,
            identity: cfg.identity,
            mapped_len: cfg.mapped_len,
            fraction: cfg.fraction,
            collect_exact_ties: mode.supports_guarded_map()
                && cfg.exact_ties == ExactTiePolicy::GuardedMap,
        })
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct CandidatePlacement {
    pub(crate) taxid: u32,
    pub(crate) profile_rank: ProfileRank,
    pub(crate) reference_id: u32,
    pub(crate) assembly_id: u32,
    pub(crate) sampling_scale: f64,
    pub(crate) asm_len: f64,
    pub(crate) r1_mapped: u32,
    pub(crate) both_mapped: u32,
    pub(crate) pair_errors: f64,
    pub(crate) long_score: f64,
    pub(crate) span: (u32, u32),
    pub(crate) fixed_chunk: (u32, u32),
    pub(crate) flex_chunk: (u32, u32),
    /// Stable order of the representative mapping in the mapper output.
    /// This is biological tie-preserving metadata, unlike a numeric reference ID.
    pub(crate) input_order: usize,
}

/// Coordinate-only payload retained for an exact short-read tie. Abundance already
/// lives in `ReadAssignment`, so keeping score/depth/string fields here would
/// multiply memory without adding information to the guarded MAP pass.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ExactTiePlacement {
    pub(crate) asm_len: f64,
    pub(crate) reference_id: u32,
    pub(crate) assembly_id: u32,
    pub(crate) span: (u32, u32),
    pub(crate) fixed_chunk: (u32, u32),
    pub(crate) flex_chunk: (u32, u32),
}

impl From<&CandidatePlacement> for ExactTiePlacement {
    fn from(placement: &CandidatePlacement) -> Self {
        Self {
            asm_len: placement.asm_len,
            reference_id: placement.reference_id,
            assembly_id: placement.assembly_id,
            span: placement.span,
            fixed_chunk: placement.fixed_chunk,
            flex_chunk: placement.flex_chunk,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ExactTieCandidate {
    pub(crate) taxid: u32,
    placement_start: usize,
    placement_len: usize,
}

/// One recorded exact-tie fragment.
///
/// [`Self::tied`] is the exact top-score prefix that guarded MAP may choose
/// between. [`Self::all`] additionally contains the strictly worse candidates
/// of the same fragment: those can never win, but once MAP resolves the
/// comparison they were still compared against the winner, so they belong in
/// every candidate's resolved-comparison denominator. Crediting only the tied
/// prefix understates that denominator for a taxon which repeatedly loses
/// ties, inflating the best-fragment fraction that the rescue tiers read.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ExactTieFragment<'a> {
    candidates: &'a [ExactTieCandidate],
    tied: usize,
}

impl<'a> ExactTieFragment<'a> {
    /// Candidates eligible to win the tie. Always at least two entries.
    #[inline]
    pub(crate) fn tied(self) -> &'a [ExactTieCandidate] {
        &self.candidates[..self.tied]
    }

    /// Every candidate compared in this fragment, tied prefix first.
    #[inline]
    pub(crate) fn all(self) -> &'a [ExactTieCandidate] {
        self.candidates
    }
}

/// Flat, append-friendly storage for exact top-score short-read ties.
///
/// Flat vectors avoid one or more heap allocations per tied fragment, which
/// matters when conserved regions generate hundreds of thousands of ties. This
/// structure remains empty whenever the resolved policy is `ambiguous`.
#[derive(Debug, Default)]
pub(crate) struct ExactTies {
    fragment_offsets: Vec<usize>,
    /// Exact top-score prefix length of each recorded fragment, parallel to
    /// the windows of `fragment_offsets`. Never exceeds the fragment's
    /// candidate count, because `record` derives it from the same slice it
    /// then stores in full.
    fragment_tied: Vec<u32>,
    candidates: Vec<ExactTieCandidate>,
    placements: Vec<ExactTiePlacement>,
}

impl ExactTies {
    pub(crate) fn record(&mut self, candidates: &[TerminalCandidate]) {
        if candidates
            .iter()
            .any(|candidate| candidate.candidate_boundary.profile_incomplete())
        {
            return;
        }
        let tied = exact_short_read_tie_prefix_len(candidates);
        if tied < 2 {
            return;
        }
        if self.fragment_offsets.is_empty() {
            self.fragment_offsets.push(0);
        }
        // Only a tied candidate can win, so only the tied prefix needs its
        // coordinate placements materialized. The remaining candidates are
        // stored for denominator accounting alone and carry no placements.
        for (index, candidate) in candidates.iter().enumerate() {
            let placement_start = self.placements.len();
            let placement_len = if index < tied {
                self.placements
                    .extend(candidate.placements.iter().map(ExactTiePlacement::from));
                candidate.placements.len()
            } else {
                0
            };
            self.candidates.push(ExactTieCandidate {
                taxid: candidate.taxid,
                placement_start,
                placement_len,
            });
        }
        self.fragment_tied.push(tied as u32);
        self.fragment_offsets.push(self.candidates.len());
    }

    pub(crate) fn append(&mut self, mut other: Self) {
        if other.fragment_count() == 0 {
            return;
        }
        if self.fragment_offsets.is_empty() {
            self.fragment_offsets.push(0);
        }
        let candidate_base = self.candidates.len();
        let placement_base = self.placements.len();
        for candidate in &mut other.candidates {
            candidate.placement_start += placement_base;
        }
        self.fragment_offsets.extend(
            other
                .fragment_offsets
                .into_iter()
                .skip(1)
                .map(|offset| candidate_base + offset),
        );
        self.candidates.append(&mut other.candidates);
        self.placements.append(&mut other.placements);
        self.fragment_tied.append(&mut other.fragment_tied);
    }

    #[inline]
    pub(crate) fn fragment_count(&self) -> usize {
        self.fragment_offsets.len().saturating_sub(1)
    }

    pub(crate) fn fragments(&self) -> impl Iterator<Item = ExactTieFragment<'_>> {
        self.fragment_offsets
            .windows(2)
            .zip(&self.fragment_tied)
            .map(|(bounds, &tied)| ExactTieFragment {
                candidates: &self.candidates[bounds[0]..bounds[1]],
                tied: tied as usize,
            })
    }

    pub(crate) fn placements(&self, candidate: ExactTieCandidate) -> &[ExactTiePlacement] {
        let end = candidate.placement_start + candidate.placement_len;
        &self.placements[candidate.placement_start..end]
    }
}

fn exact_short_read_tie_prefix_len(candidates: &[TerminalCandidate]) -> usize {
    let Some(best) = candidates.first() else {
        return 0;
    };
    candidates
        .iter()
        .take_while(|candidate| {
            (candidate.best.pair_errors - best.best.pair_errors).abs() <= SCORE_EPSILON
        })
        .count()
}

/// Resolve an already verified exact top-score tie using fixed, validated
/// safeguards. Zero/invalid priors are ineligible; one posterior must be
/// uniquely largest and account for at least 95% of eligible prior mass.
pub(crate) fn guarded_map_winner(
    candidates: &[ExactTieCandidate],
    priors: &HashMap<u32, f64>,
) -> Option<(ExactTieCandidate, f64)> {
    if candidates.len() < 2 {
        return None;
    }
    let mut denominator = 0.0;
    let mut winner: Option<(ExactTieCandidate, f64)> = None;
    let mut second_prior = 0.0f64;
    for &candidate in candidates {
        let prior = priors.get(&candidate.taxid).copied().unwrap_or(0.0);
        if !prior.is_finite() || prior <= 0.0 {
            continue;
        }
        denominator += prior;
        match winner {
            None => winner = Some((candidate, prior)),
            Some((_, winner_prior)) if prior > winner_prior => {
                second_prior = second_prior.max(winner_prior);
                winner = Some((candidate, prior));
            }
            Some(_) => second_prior = second_prior.max(prior),
        }
    }
    let (candidate, winner_prior) = winner?;
    let posterior = winner_prior / denominator;
    if !posterior.is_finite()
        || winner_prior <= second_prior * (1.0 + 1e-12)
        || posterior < GUARDED_MAP_MIN_POSTERIOR
    {
        return None;
    }
    Some((candidate, posterior))
}

#[derive(Clone, Debug)]
pub(crate) struct TerminalCandidate {
    pub(crate) taxid: u32,
    pub(crate) profile_rank: ProfileRank,
    pub(crate) best: CandidatePlacement,
    /// Taxonomic candidate levels that may be incomplete.
    pub(crate) candidate_boundary: CandidateBoundary,
    /// One independently evaluated placement per represented assembly.
    pub(crate) placements: Vec<CandidatePlacement>,
}

/// One independent terminal-taxon/assembly coordinate system. Reference
/// identifiers distinguish contigs within the assembly, so their covered
/// bases may be summed after intervals have been merged per contig.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct PresenceKey {
    pub(crate) taxid: u32,
    pub(crate) profile_rank: ProfileRank,
    pub(crate) assembly_id: u32,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct PresenceInterval {
    pub(crate) reference_id: u32,
    pub(crate) span: (u32, u32),
    pub(crate) fixed_chunk: (u32, u32),
    pub(crate) flex_chunk: (u32, u32),
}

#[derive(Clone, Debug, Default)]
pub(crate) struct PresenceAssembly {
    /// Effective indexed assembly size. This is the coordinate space that can
    /// actually receive mappings in a fractional index.
    pub(crate) asm_len: f64,
    pub(crate) count: u64,
    pub(crate) reference_span_bases: u64,
    pub(crate) intervals: Vec<PresenceInterval>,
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct TaxonEvidence {
    pub(crate) candidate_fragments: u64,
    /// Candidate comparisons for which the policy selected a direct or
    /// guarded-MAP winner. Exact ties left unresolved are excluded.
    pub(crate) resolved_candidate_fragments: u64,
    /// Fragments resolved to this taxon either directly (unique or
    /// margin-separated) or by a successful guarded-MAP exact-tie decision.
    /// Unresolved ties are deliberately credited to no taxon.
    pub(crate) best_fragments: u64,
    /// Fragments whose winning candidate was unique or margin-resolved — never
    /// taxid-tiebroken on an exact score tie. High-specificity low-amount
    /// rescue uses this (not `best_fragments`) so a conserved region's
    /// ambiguous ties cannot masquerade as own-species support.
    pub(crate) clear_best_fragments: u64,
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct EvidenceCoverage {
    pub(crate) count: usize,
    pub(crate) depth: f64,
    pub(crate) breadth: f64,
    pub(crate) expected_breadth: f64,
    pub(crate) fixed_chunk: f64,
    pub(crate) flex_chunk: f64,
    pub(crate) expected_flex_chunk: f64,
    pub(crate) oebr: f64,
    pub(crate) coebr: f64,
    pub(crate) breadth_pvalue: Option<f64>,
    pub(crate) chunk_pvalue: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct EvidenceThresholds {
    pub(crate) mode: EvidenceMode,
    pub(crate) min_breadth: f64,
    /// Chunk-breadth cutoff used by the ordinary standard evidence tier.
    pub(crate) standard_min_chunk_breadth: f64,
    /// Original, unscaled cutoff retained by rescue and guarded-MAP policy.
    pub(crate) safeguard_min_chunk_breadth: f64,
    pub(crate) chunk_filter_enabled: bool,
    pub(crate) min_oebr: Option<f64>,
    pub(crate) min_coebr: Option<f64>,
    pub(crate) strict: bool,
}

/// Independent, explainable reasons why a non-standard species witness may be
/// retained. More than one reason can apply; keeping the signals separate lets
/// guarded MAP use the appropriately strict pre-MAP evidence requirement.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct RescueKinds {
    pub(crate) platform: bool,
    pub(crate) high_specificity_low_amount: bool,
    pub(crate) low_abundance_pe: bool,
}

impl RescueKinds {
    #[inline]
    pub(crate) fn any(self) -> bool {
        self.platform || self.high_specificity_low_amount || self.low_abundance_pe
    }
}

/// Independent pre-MAP fixed-chunk evidence required before a post-MAP
/// witness may become a new species call. Encoding the four calibrated guards
/// as an enum keeps the policy attached to the decision without storing an
/// additional `f64` in every witness.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum PremapChunkGuard {
    #[default]
    None,
    FivePercent,
    Half,
    Full,
}

impl PremapChunkGuard {
    #[inline]
    pub(crate) fn fraction(self) -> Option<f64> {
        match self {
            Self::None => None,
            Self::FivePercent => Some(0.05),
            Self::Half => Some(0.50),
            Self::Full => Some(1.00),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct WitnessDecision {
    pub(crate) metrics: EvidenceCoverage,
    pub(crate) standard: bool,
    pub(crate) rescues: RescueKinds,
    pub(crate) required_premap_guard: PremapChunkGuard,
}

impl WitnessDecision {
    #[inline]
    pub(crate) fn accepted(self) -> bool {
        self.standard || self.rescues.any()
    }

    #[inline]
    pub(crate) fn rescued(self) -> bool {
        !self.standard && self.rescues.any()
    }
}

#[derive(Debug)]
struct ChunkAmount {
    observed: f64,
    standard_cutoff: f64,
    safeguard_cutoff: f64,
    filter_enabled: bool,
}

impl ChunkAmount {
    #[inline]
    fn passes_standard(&self) -> bool {
        !self.filter_enabled || self.observed > self.standard_cutoff
    }

    #[inline]
    fn exceeds_safeguard_fraction(&self, fraction: f64) -> bool {
        // Keep the established multiplication form rather than dividing to an
        // explicit ratio. Besides avoiding work in the hot path, this preserves
        // exact floating-point behavior at strict policy boundaries.
        self.filter_enabled && self.observed > self.safeguard_cutoff * fraction
    }

    #[inline]
    fn has_positive_safeguard(&self) -> bool {
        self.filter_enabled && self.safeguard_cutoff > 0.0
    }
}

/// Shared evidence features consumed by every presence tier. Constructing
/// these once prevents subtly different denominator or p-value definitions
/// from drifting between standard, rescue, and post-MAP evaluation.
#[derive(Debug)]
struct WitnessFeatures {
    amount: ChunkAmount,
    count: usize,
    depth: f64,
    breadth_ok: bool,
    oebr: f64,
    coebr: f64,
    resolved_best_fraction: f64,
    all_candidate_best_fraction: f64,
    clear_best_fraction: f64,
    clear_best_fragments: u64,
    standard_pvalue: bool,
    rescue_pvalue: bool,
}

impl WitnessFeatures {
    #[inline]
    fn new(
        metrics: &EvidenceCoverage,
        evidence: &TaxonEvidence,
        thresholds: EvidenceThresholds,
    ) -> Self {
        // Short-read specificity is measured among comparisons for which
        // direct evidence or guarded MAP selected a winner. Exact unresolved
        // ties remain only in the independent all-candidate denominator. Long
        // reads use the all-candidate denominator.
        let resolved_denominator = match thresholds.mode {
            EvidenceMode::Paired | EvidenceMode::Single => evidence.resolved_candidate_fragments,
            EvidenceMode::Long => evidence.candidate_fragments,
        };
        Self {
            amount: ChunkAmount {
                observed: metrics.fixed_chunk,
                standard_cutoff: thresholds.standard_min_chunk_breadth,
                safeguard_cutoff: thresholds.safeguard_min_chunk_breadth,
                filter_enabled: thresholds.chunk_filter_enabled,
            },
            count: metrics.count,
            depth: metrics.depth,
            breadth_ok: metrics.breadth >= thresholds.min_breadth,
            oebr: metrics.oebr,
            coebr: metrics.coebr,
            resolved_best_fraction: fraction(evidence.best_fragments, resolved_denominator),
            all_candidate_best_fraction: fraction(
                evidence.best_fragments,
                evidence.candidate_fragments,
            ),
            clear_best_fraction: fraction(
                evidence.clear_best_fragments,
                evidence.candidate_fragments,
            ),
            clear_best_fragments: evidence.clear_best_fragments,
            standard_pvalue: pvalue_passes(
                metrics.breadth_pvalue,
                metrics.chunk_pvalue,
                thresholds.strict,
            ),
            // Rescue deliberately omits the default p-value gate. An explicit
            // --strict remains an explicit user request and is honored.
            rescue_pvalue: !thresholds.strict
                || pvalue_passes(metrics.breadth_pvalue, metrics.chunk_pvalue, true),
        }
    }
}

#[derive(Debug)]
struct RatioThresholds {
    default_lower: f64,
    default_coebr_lower: f64,
    standard_oebr_min: f64,
    standard_coebr_min: f64,
}

/// Species-level OEBR floor. Long-read COEBR uses the same floor; short reads
/// use a more conservative COEBR floor.
const SPECIES_RATIO_LOWER: f64 = 0.70;
const SHORT_READ_SPECIES_COEBR_LOWER: f64 = 0.75;
/// OEBR/COEBR floor for genus and every higher rank (0.65). Higher ranks
/// aggregate reads across several genomes, so their coverage is inherently more
/// diffuse and a lower shape floor is appropriate; this is the established genus
/// fallback bound, now shared by all ranks above species.
const HIGHER_RANK_RATIO_LOWER: f64 = 0.65;
pub(crate) const STANDARD_RATIO_UPPER: f64 = 1.50;
const PVALUE_CUTOFF: f64 = 1e-5;

const SHORT_READ_ALL_CANDIDATE_BEST_MIN: f64 = 0.20;
const PAIRED_RESOLVED_BEST_MIN: f64 = 0.75;
const SINGLE_RESOLVED_BEST_MIN: f64 = 0.85;

const LOW_ABUNDANCE_MIN_COUNT: usize = 20;
const LOW_ABUNDANCE_MIN_DEPTH: f64 = 0.001;
const LOW_ABUNDANCE_MIN_CLEAR_FRAGMENTS: u64 = 10;
const LOW_ABUNDANCE_MIN_CLEAR_FRACTION: f64 = 0.25;
const LOW_ABUNDANCE_MIN_OEBR: f64 = 0.95;
const LOW_ABUNDANCE_MIN_COEBR: f64 = 0.60;

/// Default lower OEBR/COEBR bounds used by the standard evidence tier.
///
/// This is the single source of truth for direct terminal witnesses and
/// majority-LCA output filtering. LCA rows use only these standard bounds;
/// they do not inherit any species rescue policy.
#[inline]
pub(crate) fn standard_ratio_lower_bounds(rank: ProfileRank, mode: EvidenceMode) -> (f64, f64) {
    let oebr = if rank.is_species() {
        SPECIES_RATIO_LOWER
    } else {
        HIGHER_RANK_RATIO_LOWER
    };
    let coebr = if rank.is_species() && matches!(mode, EvidenceMode::Paired | EvidenceMode::Single)
    {
        SHORT_READ_SPECIES_COEBR_LOWER
    } else {
        oebr
    };
    (oebr, coebr)
}

/// Reduce one assembly in-place. Reusing the same interval vector for the
/// three sorts avoids three full-size temporary vectors and keeps peak memory
/// proportional to the stored evidence itself.
pub(crate) fn calculate_coverage(
    assembly: &mut PresenceAssembly,
    multiple_testing_taxa: usize,
) -> EvidenceCoverage {
    let genome_size = assembly.asm_len.max(1.0);
    let span_bases = merged_interval_bases(&mut assembly.intervals, |record| record.span);
    let fixed_bases = merged_interval_bases(&mut assembly.intervals, |record| record.fixed_chunk);
    let flex_bases = merged_interval_bases(&mut assembly.intervals, |record| record.flex_chunk);
    let count = usize::try_from(assembly.count).unwrap_or(usize::MAX);
    let depth = assembly.reference_span_bases as f64 / genome_size;
    let breadth = (span_bases as f64 / genome_size).min(1.0);
    let fixed_chunk = (fixed_bases as f64 / genome_size).min(1.0);
    let flex_chunk = (flex_bases as f64 / genome_size).min(1.0);
    let expected_breadth = expected_breadth(depth, assembly.count as f64);
    // Deliberately an assembly-level expectation: it assumes a full window per
    // mapping, while observed chunks are clamped to the contig each mapping
    // landed on. The mismatch is load-bearing rather than an approximation to
    // tidy up -- coverage confined to short contigs then yields a small
    // observed chunk against a genome-wide expectation, and that collapsing
    // ratio is how contamination and cross-mapping are rejected. Deriving the
    // bin count per contig instead would raise COEBR for exactly those cases,
    // so it needs joint recalibration of the ratio floors, not a local fix.
    let flex_bins = (genome_size.sqrt() / 2.0).floor().max(1.0);
    let expected_flex_chunk = 1.0 - (1.0 - 1.0 / flex_bins).powf(assembly.count as f64);
    let oebr = ratio_or_nan(breadth, expected_breadth);
    let coebr = ratio_or_nan(flex_chunk, expected_flex_chunk);
    let tests = multiple_testing_taxa.max(1) as f64;
    let breadth_pvalue =
        adjusted_deviation_pvalue(count, Some(depth), breadth, expected_breadth, oebr, tests);
    let chunk_pvalue =
        adjusted_deviation_pvalue(count, None, flex_chunk, expected_flex_chunk, coebr, tests);
    EvidenceCoverage {
        count,
        depth,
        breadth,
        expected_breadth,
        fixed_chunk,
        flex_chunk,
        expected_flex_chunk,
        oebr,
        coebr,
        breadth_pvalue,
        chunk_pvalue,
    }
}

fn merged_interval_bases<F>(intervals: &mut [PresenceInterval], select: F) -> u64
where
    F: Fn(&PresenceInterval) -> (u32, u32) + Copy,
{
    if intervals.is_empty() {
        return 0;
    }
    intervals.sort_unstable_by_key(|record| {
        let interval = select(record);
        (record.reference_id, interval.0, interval.1)
    });
    let mut total = 0u64;
    let mut reference = intervals[0].reference_id;
    let mut current = select(&intervals[0]);
    for record in &intervals[1..] {
        let interval = select(record);
        if record.reference_id == reference && interval.0 <= current.1 {
            current.1 = current.1.max(interval.1);
        } else {
            total = total.saturating_add(u64::from(current.1.saturating_sub(current.0)));
            reference = record.reference_id;
            current = interval;
        }
    }
    total.saturating_add(u64::from(current.1.saturating_sub(current.0)))
}

/// Expected fraction of a reference covered by `count` reads under the
/// Lander-Waterman uniform-placement occupancy model, given mean coverage
/// `depth` (mapped bases / reference length). Shared by the community and
/// strain profilers so the occupancy formula lives in exactly one place.
#[inline]
pub(crate) fn expected_breadth(depth: f64, count: f64) -> f64 {
    if count <= f64::EPSILON {
        0.0
    } else {
        let mean_span_fraction = (depth / count).clamp(0.0, 1.0);
        1.0 - (1.0 - mean_span_fraction).powf(count)
    }
}

#[inline]
fn ratio_or_nan(observed: f64, expected: f64) -> f64 {
    if expected > 0.0 {
        observed / expected
    } else {
        f64::NAN
    }
}

fn adjusted_deviation_pvalue(
    count: usize,
    depth: Option<f64>,
    observed: f64,
    expected: f64,
    ratio: f64,
    tests: f64,
) -> Option<f64> {
    if !ratio.is_finite() {
        return None;
    }
    if (0.90..=1.20).contains(&ratio) {
        return Some(1.0);
    }
    let side = if observed < expected {
        Tail::Lower
    } else {
        Tail::Upper
    };
    calc_breadth_pvalue(count, depth, observed, Some(expected), side)
        .map(|value| (value * tests).min(1.0))
}

#[cfg(test)]
pub(crate) fn evaluate_witness(
    metrics: EvidenceCoverage,
    evidence: TaxonEvidence,
    thresholds: EvidenceThresholds,
) -> WitnessDecision {
    evaluate_ranked_witness(metrics, evidence, thresholds, ProfileRank::Species)
}

/// Evaluate one terminal taxon's assembly-coherent witness. Species may use
/// the rescue tiers; higher-rank terminal references use only the standard
/// evidence rule. Genus uses the 0.65 default lower ratio bound.
pub(crate) fn evaluate_ranked_witness(
    metrics: EvidenceCoverage,
    evidence: TaxonEvidence,
    thresholds: EvidenceThresholds,
    rank: ProfileRank,
) -> WitnessDecision {
    let ratios = ratio_thresholds(rank, thresholds);
    let features = WitnessFeatures::new(&metrics, &evidence, thresholds);
    let standard = standard_tier(&features, &ratios);

    // A standard witness does not need a rescue label. Keeping the rescue
    // signals only for non-standard witnesses prevents a relaxed MAP guard
    // from weakening the ordinary standard-evidence path. Evaluate every
    // rescue independently because multiple reasons can match and the PE
    // low-abundance reason carries its own pre-MAP guard.
    let rescues = if rank.is_species() && features.breadth_ok && !standard {
        RescueKinds {
            platform: ordinary_platform_rescue(&features, thresholds, &ratios),
            high_specificity_low_amount: high_specificity_low_amount_rescue(&features, &ratios),
            low_abundance_pe: low_abundance_pe_rescue(&features, thresholds),
        }
    } else {
        RescueKinds::default()
    };
    let required_premap_guard = required_premap_guard(thresholds.mode, rank, standard, rescues);
    WitnessDecision {
        metrics,
        standard,
        rescues,
        required_premap_guard,
    }
}

#[inline]
fn ratio_thresholds(rank: ProfileRank, thresholds: EvidenceThresholds) -> RatioThresholds {
    // Species use the 0.70 OEBR floor. Short-read species use the 0.75 COEBR
    // floor; long reads use 0.70. Genus and every higher rank use 0.65 because
    // their coverage is inherently more diffuse after aggregation.
    let (default_lower, default_coebr_lower) = standard_ratio_lower_bounds(rank, thresholds.mode);
    RatioThresholds {
        default_lower,
        default_coebr_lower,
        standard_oebr_min: thresholds.min_oebr.unwrap_or(default_lower),
        standard_coebr_min: thresholds.min_coebr.unwrap_or(default_coebr_lower),
    }
}

#[inline]
fn standard_tier(features: &WitnessFeatures, ratios: &RatioThresholds) -> bool {
    features.breadth_ok
        && features.amount.passes_standard()
        && finite_between(
            features.oebr,
            ratios.standard_oebr_min,
            STANDARD_RATIO_UPPER,
        )
        && finite_between(
            features.coebr,
            ratios.standard_coebr_min,
            STANDARD_RATIO_UPPER,
        )
        && features.standard_pvalue
}

/// Ordinary platform rescue retains a coherent species below the full chunk
/// cutoff. Short reads additionally require at least 20% support among all
/// candidate opportunities so a small resolvable subset cannot dominate a
/// much larger unresolved conserved background.
#[inline]
fn ordinary_platform_rescue(
    features: &WitnessFeatures,
    thresholds: EvidenceThresholds,
    ratios: &RatioThresholds,
) -> bool {
    let (amount_fraction, best_fraction) = match thresholds.mode {
        EvidenceMode::Paired => (0.50, PAIRED_RESOLVED_BEST_MIN),
        EvidenceMode::Single => (0.40, SINGLE_RESOLVED_BEST_MIN),
        EvidenceMode::Long => (0.50, 0.50),
    };
    let ratios_pass = match thresholds.mode {
        EvidenceMode::Paired | EvidenceMode::Single => {
            let oebr_min = thresholds.min_oebr.unwrap_or(0.0);
            let coebr_min = thresholds.min_coebr.unwrap_or(ratios.default_coebr_lower);
            features.oebr.is_finite()
                && features.oebr >= oebr_min
                && finite_between(features.coebr, coebr_min, STANDARD_RATIO_UPPER)
        }
        EvidenceMode::Long => {
            let oebr_min = thresholds.min_oebr.unwrap_or(ratios.default_lower);
            let coebr_min = thresholds.min_coebr.unwrap_or(ratios.default_lower);
            features.oebr.is_finite()
                && features.coebr.is_finite()
                && features.oebr >= oebr_min
                && features.coebr >= coebr_min
        }
    };
    features.amount.exceeds_safeguard_fraction(amount_fraction)
        && features.resolved_best_fraction >= best_fraction
        && (thresholds.mode == EvidenceMode::Long
            || features.all_candidate_best_fraction >= SHORT_READ_ALL_CANDIDATE_BEST_MIN)
        && ratios_pass
        && features.rescue_pvalue
}

/// High-specificity low-amount rescue uses direct, unambiguous winners only.
/// Guarded-MAP wins deliberately do not contribute to clear-best support.
#[inline]
fn high_specificity_low_amount_rescue(
    features: &WitnessFeatures,
    ratios: &RatioThresholds,
) -> bool {
    features.clear_best_fraction >= 0.9
        && features.depth >= 0.01
        && features.count >= 3
        && finite_between(
            features.oebr,
            ratios.standard_oebr_min,
            STANDARD_RATIO_UPPER,
        )
        && finite_between(
            features.coebr,
            ratios.standard_coebr_min,
            STANDARD_RATIO_UPPER,
        )
        && features.standard_pvalue
}

/// Low-abundance coherent-evidence rescue is PE/species-only. Its smoothly
/// increasing amount requirement approaches the ordinary half-cutoff as C
/// grows, while independent clear support and near-expected OEBR prevent
/// conserved cross-mapping from manufacturing a low-abundance call.
#[inline]
fn low_abundance_pe_rescue(features: &WitnessFeatures, thresholds: EvidenceThresholds) -> bool {
    let oebr_min = thresholds.min_oebr.unwrap_or(LOW_ABUNDANCE_MIN_OEBR);
    let coebr_min = thresholds.min_coebr.unwrap_or(LOW_ABUNDANCE_MIN_COEBR);
    thresholds.mode == EvidenceMode::Paired
        && features.amount.has_positive_safeguard()
        && features
            .amount
            .exceeds_safeguard_fraction(low_abundance_amount_fraction(
                features.amount.safeguard_cutoff,
            ))
        && features.count >= LOW_ABUNDANCE_MIN_COUNT
        && features.depth >= LOW_ABUNDANCE_MIN_DEPTH
        && features.clear_best_fragments >= LOW_ABUNDANCE_MIN_CLEAR_FRAGMENTS
        && features.clear_best_fraction >= LOW_ABUNDANCE_MIN_CLEAR_FRACTION
        && features.resolved_best_fraction >= PAIRED_RESOLVED_BEST_MIN
        && features.all_candidate_best_fraction >= SHORT_READ_ALL_CANDIDATE_BEST_MIN
        && finite_between(features.oebr, oebr_min, STANDARD_RATIO_UPPER)
        && finite_between(features.coebr, coebr_min, STANDARD_RATIO_UPPER)
        && features.rescue_pvalue
}

#[inline]
fn required_premap_guard(
    mode: EvidenceMode,
    rank: ProfileRank,
    standard: bool,
    rescues: RescueKinds,
) -> PremapChunkGuard {
    if !rank.is_species() {
        return PremapChunkGuard::None;
    }
    match mode {
        EvidenceMode::Long => PremapChunkGuard::None,
        EvidenceMode::Single if standard => PremapChunkGuard::Full,
        EvidenceMode::Single => PremapChunkGuard::None,
        EvidenceMode::Paired if standard => PremapChunkGuard::Half,
        EvidenceMode::Paired if rescues.low_abundance_pe => PremapChunkGuard::FivePercent,
        EvidenceMode::Paired if rescues.any() => PremapChunkGuard::Half,
        EvidenceMode::Paired => PremapChunkGuard::None,
    }
}

#[inline]
fn fraction(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

#[inline]
fn low_abundance_amount_fraction(min_chunk_breadth: f64) -> f64 {
    (0.20 + 0.45 * min_chunk_breadth.clamp(0.0, 1.0)).min(0.50)
}

#[inline]
fn finite_between(value: f64, lower: f64, upper: f64) -> bool {
    value.is_finite() && lower.is_finite() && upper.is_finite() && value >= lower && value <= upper
}

#[inline]
fn pvalue_passes(breadth: Option<f64>, chunk: Option<f64>, strict: bool) -> bool {
    match (breadth, chunk) {
        (Some(breadth), Some(chunk)) if strict => {
            breadth >= PVALUE_CUTOFF && chunk >= PVALUE_CUTOFF
        }
        (Some(breadth), Some(chunk)) => !(breadth < PVALUE_CUTOFF && chunk < PVALUE_CUTOFF),
        _ => !strict,
    }
}

/// Prefer a passing coherent assembly over a failing one: standard first,
/// then rescue, then the largest fixed-chunk breadth. Testing breadth before
/// rescue status could let a larger but incoherent assembly hide a smaller
/// assembly that independently satisfies the rescue rule.
pub(crate) fn witness_is_better(left: WitnessDecision, right: WitnessDecision) -> bool {
    if left.standard != right.standard {
        return left.standard;
    }
    if left.rescued() != right.rescued() {
        return left.rescued();
    }
    if left.metrics.fixed_chunk != right.metrics.fixed_chunk {
        return left.metrics.fixed_chunk > right.metrics.fixed_chunk;
    }
    if left.metrics.breadth != right.metrics.breadth {
        return left.metrics.breadth > right.metrics.breadth;
    }
    if left.metrics.depth != right.metrics.depth {
        return left.metrics.depth > right.metrics.depth;
    }
    left.metrics.count > right.metrics.count
}

/// Reusable per-worker candidate storage. Hash-table capacities survive from
/// one fragment to the next, avoiding allocator traffic in the hot loop.
#[derive(Default)]
pub(crate) struct CandidateScratch {
    by_taxon_assembly: HashMap<(u32, ProfileRank, u32), CandidatePlacement>,
    by_taxon: HashMap<(u32, ProfileRank), TerminalCandidate>,
    paired_rows: HashMap<u32, ExplicitPairRows>,
    pair_order: Vec<u32>,
    output: Vec<TerminalCandidate>,
}

#[derive(Clone, Copy, Debug, Default)]
struct ExplicitPairRows {
    mate1: Option<usize>,
    mate2: Option<usize>,
    malformed: bool,
}

impl CandidateScratch {
    pub(crate) fn build<'a>(
        &'a mut self,
        rows: &[EvidenceAlignment],
        policy: EvidencePolicy,
    ) -> &'a [TerminalCandidate] {
        let candidate_boundary = rows
            .iter()
            .fold(CandidateBoundary::COMPLETE, |boundary, mapping| {
                boundary.union(mapping.candidate_boundary)
            });
        self.by_taxon_assembly.clear();
        self.by_taxon.clear();
        self.paired_rows.clear();
        self.pair_order.clear();
        self.output.clear();

        match policy.mode {
            EvidenceMode::Paired => self.build_paired(rows, policy),
            EvidenceMode::Single => self.build_single(rows, policy),
            EvidenceMode::Long => self.build_long(rows, policy),
        }

        for (_, placement) in self.by_taxon_assembly.drain() {
            match self
                .by_taxon
                .entry((placement.taxid, placement.profile_rank))
            {
                Entry::Vacant(entry) => {
                    entry.insert(TerminalCandidate {
                        taxid: placement.taxid,
                        profile_rank: placement.profile_rank,
                        best: placement,
                        candidate_boundary,
                        placements: vec![placement],
                    });
                }
                Entry::Occupied(mut entry) => {
                    let taxon = entry.get_mut();
                    if placement_is_better(&placement, &taxon.best, policy.mode) {
                        taxon.best = placement;
                    }
                    taxon.placements.push(placement);
                }
            }
        }
        self.output
            .extend(self.by_taxon.drain().map(|(_, value)| value));
        for candidate in &mut self.output {
            candidate
                .placements
                .sort_unstable_by_key(|placement| placement.assembly_id);
        }
        self.output
            .sort_by(|left, right| candidate_order(left, right, policy.mode));
        &self.output
    }

    fn build_paired(&mut self, rows: &[EvidenceAlignment], policy: EvidencePolicy) {
        for (index, row) in rows.iter().enumerate() {
            if row.pair_id == 0 {
                continue;
            }
            let entry = match self.paired_rows.entry(row.pair_id) {
                Entry::Vacant(entry) => {
                    self.pair_order.push(row.pair_id);
                    entry.insert(ExplicitPairRows::default())
                }
                Entry::Occupied(entry) => entry.into_mut(),
            };
            let slot = match row.role {
                MappingRole::Mate1 => &mut entry.mate1,
                MappingRole::Mate2 => &mut entry.mate2,
                MappingRole::Unpaired => {
                    entry.malformed = true;
                    continue;
                }
            };
            if slot.replace(index).is_some() {
                entry.malformed = true;
            }
        }

        for order_index in 0..self.pair_order.len() {
            let pair_id = self.pair_order[order_index];
            let Some(pair) = self.paired_rows.get(&pair_id).copied() else {
                continue;
            };
            let (Some(r1_index), Some(r2_index)) = (pair.mate1, pair.mate2) else {
                continue;
            };
            if pair.malformed {
                continue;
            }
            let r1 = &rows[r1_index];
            let r2 = &rows[r2_index];
            if r1.reference_id != r2.reference_id
                || r1.profile_taxid != r2.profile_taxid
                || r1.profile_rank != r2.profile_rank
                || r1.assembly_id != r2.assembly_id
                || !paired_mate1_qualifies(r1, policy)
                || !paired_mate2_qualifies(r2, policy)
            {
                continue;
            }
            let placement = CandidatePlacement {
                taxid: r1.profile_taxid,
                profile_rank: r1.profile_rank,
                reference_id: r1.reference_id,
                assembly_id: r1.assembly_id,
                sampling_scale: f64::from(r1.sampling_scale),
                asm_len: r1.asm_len,
                r1_mapped: r1.mapped_len,
                both_mapped: r1.mapped_len.saturating_add(r2.mapped_len),
                pair_errors: estimated_errors(r1) + estimated_errors(r2),
                long_score: 0.0,
                span: r1.span,
                fixed_chunk: r1.fixed_chunk,
                flex_chunk: r1.flex_chunk,
                input_order: r1_index.min(r2_index),
            };
            self.insert_placement(placement, policy.mode);
        }
    }

    fn build_single(&mut self, rows: &[EvidenceAlignment], policy: EvidencePolicy) {
        for (input_order, mapping) in rows.iter().enumerate() {
            if !single_mapping_qualifies(mapping, policy) {
                continue;
            }
            let placement = CandidatePlacement {
                taxid: mapping.profile_taxid,
                profile_rank: mapping.profile_rank,
                reference_id: mapping.reference_id,
                assembly_id: mapping.assembly_id,
                sampling_scale: f64::from(mapping.sampling_scale),
                asm_len: mapping.asm_len,
                r1_mapped: mapping.mapped_len,
                both_mapped: mapping.mapped_len,
                pair_errors: estimated_errors(mapping),
                long_score: 0.0,
                span: mapping.span,
                fixed_chunk: mapping.fixed_chunk,
                flex_chunk: mapping.flex_chunk,
                input_order,
            };
            self.insert_placement(placement, policy.mode);
        }
    }

    fn build_long(&mut self, rows: &[EvidenceAlignment], policy: EvidencePolicy) {
        for (input_order, mapping) in rows.iter().enumerate() {
            if !long_mapping_qualifies(mapping, policy) {
                continue;
            }
            let placement = CandidatePlacement {
                taxid: mapping.profile_taxid,
                profile_rank: mapping.profile_rank,
                reference_id: mapping.reference_id,
                assembly_id: mapping.assembly_id,
                sampling_scale: f64::from(mapping.sampling_scale),
                asm_len: mapping.asm_len,
                r1_mapped: mapping.mapped_len,
                both_mapped: mapping.mapped_len,
                pair_errors: estimated_errors(mapping),
                long_score: long_score(mapping),
                span: mapping.span,
                fixed_chunk: mapping.fixed_chunk,
                flex_chunk: mapping.flex_chunk,
                input_order,
            };
            self.insert_placement(placement, policy.mode);
        }
    }

    fn insert_placement(&mut self, placement: CandidatePlacement, mode: EvidenceMode) {
        let key = (
            placement.taxid,
            placement.profile_rank,
            placement.assembly_id,
        );
        match self.by_taxon_assembly.entry(key) {
            Entry::Vacant(entry) => {
                entry.insert(placement);
            }
            Entry::Occupied(mut entry) => {
                if placement_is_better(&placement, entry.get(), mode) {
                    entry.insert(placement);
                }
            }
        }
    }
}

#[inline]
fn paired_mate1_qualifies(mapping: &EvidenceAlignment, policy: EvidencePolicy) -> bool {
    evidence_identity(mapping) >= policy.identity
        && mapping.mapped_len as usize >= policy.mapped_len
        && fraction_at_least(mapping.fraction, policy.fraction)
        && mapping.span.1 > mapping.span.0
}

#[inline]
fn paired_mate2_qualifies(mapping: &EvidenceAlignment, policy: EvidencePolicy) -> bool {
    mapping.mapped_len as usize >= policy.mapped_len
        && fraction_at_least(mapping.fraction, policy.fraction)
        && mapping.span.1 > mapping.span.0
}

#[inline]
fn single_mapping_qualifies(mapping: &EvidenceAlignment, policy: EvidencePolicy) -> bool {
    evidence_identity(mapping) >= policy.identity
        && mapping.mapped_len as usize >= policy.mapped_len
        && fraction_at_least(mapping.fraction, policy.fraction)
        && mapping.span.1 > mapping.span.0
}

#[inline]
fn long_mapping_qualifies(mapping: &EvidenceAlignment, policy: EvidencePolicy) -> bool {
    evidence_identity(mapping) >= policy.identity
        && (mapping.mapped_len as usize >= policy.mapped_len
            || fraction_at_least(mapping.fraction, policy.fraction))
        && mapping.span.1 > mapping.span.0
}

#[inline]
fn estimated_errors(mapping: &EvidenceAlignment) -> f64 {
    (1.0 - evidence_identity(mapping)).max(0.0) * f64::from(mapping.mapped_len)
}

#[inline]
fn long_score(mapping: &EvidenceAlignment) -> f64 {
    f64::from(mapping.mapped_len) * (2.0 * evidence_identity(mapping) - 1.0)
}

#[inline]
fn evidence_identity(mapping: &EvidenceAlignment) -> f64 {
    // Mapping identities are canonical five-decimal values at the mapper /
    // profiler boundary. The compact in-memory record stores them as f32;
    // canonicalizing after widening recovers the exact boundary value and
    // prevents f32 representation noise from splitting an exact score tie.
    canonical_mapping_identity(f64::from(mapping.identity))
}

#[inline]
fn fraction_at_least(fraction: f32, threshold: f64) -> bool {
    // An exact threshold such as 105/150 = 0.70 may round just below 0.70 in
    // f32. One f32 epsilon restores inclusive-threshold semantics without
    // materially relaxing the biological cutoff.
    f64::from(fraction) + FRACTION_F32_EPSILON >= threshold
}

#[inline]
fn placement_is_better(
    left: &CandidatePlacement,
    right: &CandidatePlacement,
    mode: EvidenceMode,
) -> bool {
    match mode {
        EvidenceMode::Paired | EvidenceMode::Single => {
            // Short-read placements are independently ranked; use the stable
            // reference ID only after all biological evidence ties.
            left.pair_errors < right.pair_errors - SCORE_EPSILON
                || ((left.pair_errors - right.pair_errors).abs() <= SCORE_EPSILON
                    && (left.both_mapped > right.both_mapped
                        || (left.both_mapped == right.both_mapped
                            && left.reference_id < right.reference_id)))
        }
        EvidenceMode::Long => {
            // Preserve the mapper's first-ranked long-read placement on an
            // exact semantic tie. A numeric reference ID is metadata, not
            // biological evidence; using it here can move an otherwise
            // identical candidate to a different contig and manufacture a
            // different coverage pattern.
            left.long_score > right.long_score + SCORE_EPSILON
                || ((left.long_score - right.long_score).abs() <= SCORE_EPSILON
                    && (left.pair_errors < right.pair_errors - SCORE_EPSILON
                        || ((left.pair_errors - right.pair_errors).abs() <= SCORE_EPSILON
                            && (left.both_mapped > right.both_mapped
                                || (left.both_mapped == right.both_mapped
                                    && left.input_order < right.input_order)))))
        }
    }
}

fn candidate_order(
    left: &TerminalCandidate,
    right: &TerminalCandidate,
    mode: EvidenceMode,
) -> Ordering {
    match mode {
        EvidenceMode::Paired | EvidenceMode::Single => left
            .best
            .pair_errors
            .partial_cmp(&right.best.pair_errors)
            .unwrap_or(Ordering::Equal)
            .then_with(|| left.taxid.cmp(&right.taxid))
            .then_with(|| left.profile_rank.cmp(&right.profile_rank)),
        EvidenceMode::Long => right
            .best
            .long_score
            .partial_cmp(&left.best.long_score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| left.best.input_order.cmp(&right.best.input_order))
            .then_with(|| left.taxid.cmp(&right.taxid))
            .then_with(|| left.profile_rank.cmp(&right.profile_rank)),
    }
}

pub(crate) fn resolved_candidate_index(
    candidates: &[TerminalCandidate],
    mode: EvidenceMode,
) -> Option<usize> {
    if candidates.is_empty()
        || candidates
            .iter()
            .any(|candidate| candidate.candidate_boundary.profile_incomplete())
    {
        return None;
    }
    if candidates.len() == 1 {
        return Some(0);
    }
    let resolved = match mode {
        EvidenceMode::Paired => {
            candidates[1].best.pair_errors - candidates[0].best.pair_errors >= PAIRED_ERROR_MARGIN
        }
        EvidenceMode::Single => {
            candidates[1].best.pair_errors - candidates[0].best.pair_errors >= SINGLE_ERROR_MARGIN
        }
        EvidenceMode::Long => {
            let best = candidates[0].best.long_score;
            best > f64::EPSILON
                && (best - candidates[1].best.long_score) / best >= LONG_RELATIVE_SCORE_MARGIN
        }
    };
    resolved.then_some(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mapping(
        species: u32,
        reference: u32,
        assembly: u32,
        identity: f32,
        mapped: u32,
        read_len: u32,
        role: MappingRole,
    ) -> EvidenceAlignment {
        mapping_at_rank(
            species,
            ProfileRank::Species,
            reference,
            assembly,
            identity,
            mapped,
            read_len,
            role,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn mapping_at_rank(
        species: u32,
        profile_rank: ProfileRank,
        reference: u32,
        assembly: u32,
        identity: f32,
        mapped: u32,
        read_len: u32,
        role: MappingRole,
    ) -> EvidenceAlignment {
        EvidenceAlignment {
            asm: Arc::<str>::from(format!("asm{assembly}")),
            acc: Arc::<str>::from(format!("ref{reference}")),
            profile_taxid: species,
            profile_rank,
            reference_id: reference,
            assembly_id: assembly,
            role,
            pair_id: if matches!(role, MappingRole::Mate1 | MappingRole::Mate2) {
                reference.max(1)
            } else {
                0
            },
            candidate_boundary: CandidateBoundary::COMPLETE,
            sampling_scale: 1.0,
            asm_len: 1_000_000.0,
            identity,
            fraction: mapped as f32 / read_len as f32,
            mapped_len: mapped,
            span: (100, 100 + mapped),
            fixed_chunk: (1, 1_100 + mapped),
            flex_chunk: (1, 2_100 + mapped),
        }
    }

    fn policy(mode: EvidenceMode) -> EvidencePolicy {
        EvidencePolicy {
            mode,
            identity: 0.97,
            mapped_len: 50,
            fraction: 0.70,
            collect_exact_ties: false,
        }
    }

    fn thresholds(mode: EvidenceMode) -> EvidenceThresholds {
        EvidenceThresholds {
            mode,
            min_breadth: 0.0,
            standard_min_chunk_breadth: 0.60,
            safeguard_min_chunk_breadth: 0.60,
            chunk_filter_enabled: true,
            min_oebr: None,
            min_coebr: None,
            strict: false,
        }
    }

    fn metrics(oebr: f64, coebr: f64, fixed_chunk: f64) -> EvidenceCoverage {
        EvidenceCoverage {
            count: 100,
            depth: 0.1,
            breadth: 0.1,
            expected_breadth: 0.1,
            fixed_chunk,
            flex_chunk: 0.1,
            expected_flex_chunk: 0.1,
            oebr,
            coebr,
            breadth_pvalue: Some(1.0),
            chunk_pvalue: Some(1.0),
        }
    }

    fn taxon_evidence(
        candidates: u64,
        resolved_candidates: u64,
        best: u64,
        clear_best: u64,
    ) -> TaxonEvidence {
        TaxonEvidence {
            candidate_fragments: candidates,
            resolved_candidate_fragments: resolved_candidates,
            best_fragments: best,
            clear_best_fragments: clear_best,
        }
    }

    #[test]
    fn paired_candidates_use_explicit_pair_ids_not_row_adjacency() {
        let mut rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        rows[0].pair_id = 101;
        rows[2].pair_id = 101;
        rows[1].pair_id = 202;
        rows[3].pair_id = 202;

        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 2);
        assert_eq!(candidates[0].best.both_mapped, 200);
        assert_eq!(candidates[1].best.both_mapped, 200);
    }

    #[test]
    fn paired_candidate_requires_identity_only_on_mate1() {
        let rows = vec![
            mapping(1, 10, 1, 0.98, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.80, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].best.both_mapped, 200);
    }

    #[test]
    fn compact_identity_and_fraction_preserve_canonical_inclusive_thresholds() {
        let rows = vec![
            mapping(1, 10, 1, 0.97, 70, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.80, 70, 100, MappingRole::Mate2),
        ];
        // 0.70 is represented slightly below 0.70 as f32 on common targets;
        // the inclusive profiler cutoff must nevertheless accept it.
        assert!(f64::from(rows[0].fraction) < 0.70);
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(evidence_identity(&rows[0]), 0.97);
    }

    #[test]
    fn exact_cross_species_pair_tie_remains_unresolved() {
        let rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(
            resolved_candidate_index(candidates, EvidenceMode::Paired),
            None
        );
    }

    #[test]
    fn guarded_map_resolves_only_a_unique_95_percent_exact_tie_winner() {
        let rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        let mut ties = ExactTies::default();
        ties.record(candidates);
        assert_eq!(ties.fragment_count(), 1);
        let fragment = ties.fragments().next().unwrap().tied();

        let priors = HashMap::from([(1, 95.0), (2, 5.0)]);
        let (winner, posterior) = guarded_map_winner(fragment, &priors).unwrap();
        assert_eq!(winner.taxid, 1);
        assert!((posterior - 0.95).abs() < 1e-12);
        assert_eq!(ties.placements(winner).len(), 1);

        let below_threshold = HashMap::from([(1, 94.0), (2, 6.0)]);
        assert!(guarded_map_winner(fragment, &below_threshold).is_none());
        let equal = HashMap::from([(1, 1.0), (2, 1.0)]);
        assert!(guarded_map_winner(fragment, &equal).is_none());
    }

    #[test]
    fn exact_tie_records_losing_candidates_for_denominator_accounting() {
        // Species 1 and 2 tie exactly; species 3 is strictly worse but still
        // clears the read-level gates, so it was a candidate in the same
        // comparison. It must reach denominator accounting without becoming
        // eligible to win the tie.
        let rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
            mapping(3, 30, 3, 0.98, 100, 100, MappingRole::Mate1),
            mapping(3, 30, 3, 0.98, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 3);
        let mut ties = ExactTies::default();
        ties.record(candidates);
        assert_eq!(ties.fragment_count(), 1);

        let fragment = ties.fragments().next().unwrap();
        assert_eq!(
            fragment
                .tied()
                .iter()
                .map(|candidate| candidate.taxid)
                .collect::<Vec<_>>(),
            vec![1, 2],
            "only the exact top-score prefix may win the tie"
        );
        assert_eq!(
            fragment
                .all()
                .iter()
                .map(|candidate| candidate.taxid)
                .collect::<Vec<_>>(),
            vec![1, 2, 3],
            "every compared candidate must reach denominator accounting"
        );
        // A losing candidate carries no coordinates: it can never win, so
        // materializing its placements would cost memory for nothing.
        assert!(ties.placements(fragment.all()[2]).is_empty());
        assert_eq!(ties.placements(fragment.tied()[0]).len(), 1);

        // The worse candidate must not take the tie even with an overwhelming
        // prior, because it is not part of the tied prefix.
        let priors = HashMap::from([(3, 1_000.0), (1, 95.0), (2, 5.0)]);
        let (winner, _) = guarded_map_winner(fragment.tied(), &priors).unwrap();
        assert_eq!(winner.taxid, 1);
    }

    #[test]
    fn nonzero_pair_score_gap_is_not_stored_as_an_exact_tie() {
        let rows = vec![
            mapping(1, 10, 1, 1.00, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 1.00, 100, 100, MappingRole::Mate2),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        let mut ties = ExactTies::default();
        ties.record(candidates);
        assert_eq!(ties.fragment_count(), 0);
    }

    #[test]
    fn one_base_equivalent_pair_margin_resolves_best_species() {
        let rows = vec![
            mapping(1, 10, 1, 1.00, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 1.00, 100, 100, MappingRole::Mate2),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate1),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(
            resolved_candidate_index(candidates, EvidenceMode::Paired),
            Some(0)
        );
        assert_eq!(candidates[0].taxid, 1);
    }

    #[test]
    fn single_end_uses_the_validated_half_error_margin() {
        let mut scratch = CandidateScratch::default();
        let resolved_rows = vec![
            mapping(1, 10, 1, 1.000, 100, 100, MappingRole::Unpaired),
            mapping(2, 20, 2, 0.995, 100, 100, MappingRole::Unpaired),
        ];
        let candidates = scratch.build(&resolved_rows, policy(EvidenceMode::Single));
        assert_eq!(candidates.len(), 2);
        assert_eq!(
            resolved_candidate_index(candidates, EvidenceMode::Single),
            Some(0)
        );

        let unresolved_rows = vec![
            mapping(1, 10, 1, 1.000, 100, 100, MappingRole::Unpaired),
            mapping(2, 20, 2, 0.996, 100, 100, MappingRole::Unpaired),
        ];
        let candidates = scratch.build(&unresolved_rows, policy(EvidenceMode::Single));
        assert_eq!(
            resolved_candidate_index(candidates, EvidenceMode::Single),
            None
        );
    }

    #[test]
    fn candidate_boundary_keeps_candidate_but_blocks_coordinate_resolution() {
        let mut row = mapping(101, 1, 1, 0.99, 100, 100, MappingRole::Unpaired);
        row.candidate_boundary = CandidateBoundary::from_bits(1).unwrap();
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&[row], policy(EvidenceMode::Single));

        assert_eq!(candidates.len(), 1);
        assert!(candidates[0].candidate_boundary.profile_incomplete());
        assert_eq!(
            resolved_candidate_index(candidates, EvidenceMode::Single),
            None
        );
    }

    #[test]
    fn terminal_only_boundary_does_not_block_community_resolution() {
        let mut row = mapping(101, 1, 1, 0.99, 100, 100, MappingRole::Unpaired);
        row.candidate_boundary = CandidateBoundary::from_bits(2).unwrap();
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&[row], policy(EvidenceMode::Single));

        assert_eq!(candidates.len(), 1);
        assert!(!candidates[0].candidate_boundary.profile_incomplete());
        assert!(candidates[0].candidate_boundary.terminal_incomplete());
        assert_eq!(
            resolved_candidate_index(candidates, EvidenceMode::Single),
            Some(0)
        );
    }

    #[test]
    fn single_end_exact_ties_are_available_to_guarded_map() {
        let rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Unpaired),
            mapping(2, 20, 2, 0.99, 100, 100, MappingRole::Unpaired),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Single));
        assert_eq!(
            resolved_candidate_index(candidates, EvidenceMode::Single),
            None
        );
        let mut ties = ExactTies::default();
        ties.record(candidates);
        assert_eq!(ties.fragment_count(), 1);
        let winner = guarded_map_winner(
            ties.fragments().next().unwrap().tied(),
            &HashMap::from([(1, 95.0), (2, 5.0)]),
        )
        .unwrap();
        assert_eq!(winner.0.taxid, 1);
    }

    #[test]
    fn long_score_prefers_substantially_longer_near_exact_mapping() {
        let mut long_policy = policy(EvidenceMode::Long);
        long_policy.identity = 0.92;
        long_policy.mapped_len = 500;
        long_policy.fraction = 0.20;
        let rows = vec![
            mapping(1, 10, 1, 1.00, 500, 2_000, MappingRole::Unpaired),
            mapping(2, 20, 2, 0.98, 1_000, 2_000, MappingRole::Unpaired),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, long_policy);
        assert_eq!(candidates[0].taxid, 2);
        assert_eq!(
            resolved_candidate_index(candidates, EvidenceMode::Long),
            Some(0)
        );
    }

    #[test]
    fn long_exact_placement_tie_preserves_mapper_order() {
        let mut long_policy = policy(EvidenceMode::Long);
        long_policy.identity = 0.92;
        long_policy.mapped_len = 500;
        long_policy.fraction = 0.20;
        let rows = vec![
            mapping(1, 20, 1, 0.99, 1_000, 2_000, MappingRole::Unpaired),
            mapping(1, 10, 1, 0.99, 1_000, 2_000, MappingRole::Unpaired),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, long_policy);
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].placements.len(), 1);
        assert_eq!(candidates[0].placements[0].reference_id, 20);
    }

    #[test]
    fn long_exact_tie_across_assemblies_preserves_mapper_order_deterministically() {
        let mut long_policy = policy(EvidenceMode::Long);
        long_policy.identity = 0.92;
        long_policy.mapped_len = 500;
        long_policy.fraction = 0.20;
        let rows = vec![
            mapping(1, 20, 2, 0.99, 1_000, 2_000, MappingRole::Unpaired),
            mapping(1, 10, 1, 0.99, 1_000, 2_000, MappingRole::Unpaired),
        ];
        let mut scratch = CandidateScratch::default();
        for _ in 0..100 {
            let candidates = scratch.build(&rows, long_policy);
            assert_eq!(candidates.len(), 1);
            assert_eq!(candidates[0].placements.len(), 2);
            assert_eq!(candidates[0].best.reference_id, 20);
            assert_eq!(candidates[0].best.assembly_id, 2);
        }
    }

    #[test]
    fn mixed_terminal_ranks_remain_distinct_abundance_candidates() {
        let rows = vec![
            mapping_at_rank(
                111,
                ProfileRank::Genus,
                10,
                1,
                0.99,
                100,
                100,
                MappingRole::Unpaired,
            ),
            mapping_at_rank(
                1111,
                ProfileRank::Species,
                20,
                2,
                0.99,
                100,
                100,
                MappingRole::Unpaired,
            ),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Single));
        assert_eq!(candidates.len(), 2);
        assert!(candidates.iter().any(
            |candidate| candidate.taxid == 111 && candidate.profile_rank == ProfileRank::Genus
        ));
        assert!(candidates
            .iter()
            .any(|candidate| candidate.taxid == 1111
                && candidate.profile_rank == ProfileRank::Species));
        assert_eq!(
            resolved_candidate_index(candidates, EvidenceMode::Single),
            None
        );
    }

    #[test]
    fn long_read_mixed_rank_tie_remains_ambiguous_and_ordered() {
        let mut long_policy = policy(EvidenceMode::Long);
        long_policy.identity = 0.92;
        long_policy.mapped_len = 500;
        long_policy.fraction = 0.20;
        let rows = vec![
            mapping_at_rank(
                111,
                ProfileRank::Genus,
                10,
                1,
                0.98,
                1_000,
                2_000,
                MappingRole::Unpaired,
            ),
            mapping_at_rank(
                1111,
                ProfileRank::Species,
                20,
                2,
                0.98,
                1_000,
                2_000,
                MappingRole::Unpaired,
            ),
        ];
        let mut scratch = CandidateScratch::default();
        for _ in 0..100 {
            let candidates = scratch.build(&rows, long_policy);
            assert_eq!(candidates.len(), 2);
            assert_eq!(candidates[0].taxid, 111);
            assert_eq!(candidates[0].profile_rank, ProfileRank::Genus);
            assert_eq!(
                resolved_candidate_index(candidates, EvidenceMode::Long),
                None
            );
        }
    }

    #[test]
    fn paired_exact_placement_tie_uses_stable_reference_id() {
        let rows = vec![
            mapping(1, 20, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 20, 1, 0.99, 100, 100, MappingRole::Mate2),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].placements.len(), 1);
        assert_eq!(candidates[0].placements[0].reference_id, 10);
    }

    #[test]
    fn placements_from_two_assemblies_remain_separate_witnesses() {
        let rows = vec![
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate1),
            mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Mate2),
            mapping(1, 20, 2, 0.98, 100, 100, MappingRole::Mate1),
            mapping(1, 20, 2, 0.98, 100, 100, MappingRole::Mate2),
        ];
        let mut scratch = CandidateScratch::default();
        let candidates = scratch.build(&rows, policy(EvidenceMode::Paired));
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].placements.len(), 2);
        assert_ne!(
            candidates[0].placements[0].assembly_id,
            candidates[0].placements[1].assembly_id
        );
    }

    #[test]
    fn coverage_merges_within_contigs_but_not_across_contigs() {
        let mut assembly = PresenceAssembly {
            asm_len: 1_000.0,
            count: 3,
            reference_span_bases: 250,
            intervals: vec![
                PresenceInterval {
                    reference_id: 1,
                    span: (0, 100),
                    fixed_chunk: (0, 200),
                    flex_chunk: (0, 300),
                },
                PresenceInterval {
                    reference_id: 1,
                    span: (50, 150),
                    fixed_chunk: (100, 300),
                    flex_chunk: (200, 400),
                },
                PresenceInterval {
                    reference_id: 2,
                    span: (0, 50),
                    fixed_chunk: (0, 100),
                    flex_chunk: (0, 100),
                },
            ],
        };
        let coverage = calculate_coverage(&mut assembly, 1);
        assert!((coverage.breadth - 0.20).abs() < 1e-12);
        assert!((coverage.fixed_chunk - 0.40).abs() < 1e-12);
        assert!((coverage.flex_chunk - 0.50).abs() < 1e-12);
        assert!((coverage.depth - 0.25).abs() < 1e-12);
    }

    #[test]
    fn paired_rescue_keeps_low_oebr_when_chunk_distribution_is_coherent() {
        let decision = evaluate_witness(
            metrics(0.10, 1.0, 0.40),
            taxon_evidence(10, 10, 8, 8),
            thresholds(EvidenceMode::Paired),
        );
        assert!(!decision.standard);
        assert!(decision.rescues.platform);
    }

    #[test]
    fn platform_rescue_does_not_credit_unresolved_ties_to_taxid_order() {
        let decision = evaluate_witness(
            metrics(0.10, 1.0, 0.40),
            taxon_evidence(10, 0, 0, 0),
            thresholds(EvidenceMode::Paired),
        );
        assert!(!decision.standard);
        assert!(
            !decision.rescued(),
            "unresolved ties must not satisfy the platform best-support gate"
        );
    }

    #[test]
    fn direct_genus_witness_uses_genus_standard_without_species_rescue() {
        let evidence = taxon_evidence(10, 10, 10, 10);
        // 0.68 sits between the genus floor (0.65) and the species floor (0.70):
        // genus accepts it as a standard witness, species does not.
        let mut coverage = metrics(0.68, 0.68, 0.70);
        coverage.breadth = 0.10;
        coverage.fixed_chunk = 0.80;
        coverage.breadth_pvalue = Some(1.0);
        coverage.chunk_pvalue = Some(1.0);

        let genus = evaluate_ranked_witness(
            coverage,
            evidence,
            thresholds(EvidenceMode::Single),
            ProfileRank::Genus,
        );
        assert!(genus.standard);
        assert!(!genus.rescued());
        assert_eq!(genus.required_premap_guard, PremapChunkGuard::None);

        let species = evaluate_ranked_witness(
            coverage,
            evidence,
            thresholds(EvidenceMode::Single),
            ProfileRank::Species,
        );
        assert!(!species.standard);
    }

    #[test]
    fn standard_ratio_defaults_are_rank_and_platform_aware() {
        assert_eq!(
            standard_ratio_lower_bounds(ProfileRank::Species, EvidenceMode::Paired),
            (0.70, 0.75)
        );
        assert_eq!(
            standard_ratio_lower_bounds(ProfileRank::Species, EvidenceMode::Single),
            (0.70, 0.75)
        );
        assert_eq!(
            standard_ratio_lower_bounds(ProfileRank::Species, EvidenceMode::Long),
            (0.70, 0.70)
        );
        for rank in [
            ProfileRank::Genus,
            ProfileRank::Family,
            ProfileRank::Order,
            ProfileRank::Class,
            ProfileRank::Phylum,
            ProfileRank::Kingdom,
        ] {
            for mode in [
                EvidenceMode::Paired,
                EvidenceMode::Single,
                EvidenceMode::Long,
            ] {
                assert_eq!(standard_ratio_lower_bounds(rank, mode), (0.65, 0.65));
            }
        }
    }

    #[test]
    fn short_read_platform_rescue_uses_resolved_and_all_candidate_guards() {
        let policy = thresholds(EvidenceMode::Paired);
        let coverage = metrics(0.10, 1.0, 0.40);

        let too_small_resolvable_subset =
            evaluate_witness(coverage, taxon_evidence(100, 20, 15, 0), policy);
        assert!(
            !too_small_resolvable_subset.rescued(),
            "75% of resolved comparisons is insufficient when fewer than 20% of all comparisons support the taxon"
        );

        let coherent = evaluate_witness(coverage, taxon_evidence(100, 20, 20, 0), policy);
        assert!(coherent.rescues.platform);
    }

    #[test]
    fn single_end_rescue_requires_two_fifths_chunk_and_85_percent_resolved_support() {
        let policy = thresholds(EvidenceMode::Single);
        let at_amount_boundary = evaluate_witness(
            metrics(0.10, 1.0, 0.24),
            taxon_evidence(20, 20, 20, 20),
            policy,
        );
        assert!(
            !at_amount_boundary.rescued(),
            "the 40 percent amount comparison is strict"
        );

        let insufficient_support = evaluate_witness(
            metrics(0.10, 1.0, 0.25),
            taxon_evidence(20, 20, 16, 16),
            policy,
        );
        assert!(!insufficient_support.rescued());

        let validated = evaluate_witness(
            metrics(0.10, 1.0, 0.25),
            taxon_evidence(20, 20, 17, 17),
            policy,
        );
        assert!(!validated.standard);
        assert!(validated.rescues.platform);
    }

    #[test]
    fn long_rescue_retains_oebr_floor_but_has_no_upper_ratio_bound() {
        let evidence = taxon_evidence(10, 10, 5, 5);
        // 0.68 is below the 0.70 species floor the long rescue enforces on OEBR.
        let low = evaluate_witness(
            metrics(0.68, 1.0, 0.40),
            evidence,
            thresholds(EvidenceMode::Long),
        );
        assert!(!low.rescued());
        let high = evaluate_witness(
            metrics(2.0, 3.0, 0.40),
            evidence,
            thresholds(EvidenceMode::Long),
        );
        assert!(high.rescues.platform);
    }

    #[test]
    fn paired_platform_rescue_requires_half_chunk_and_75_percent_resolved_support() {
        let policy = thresholds(EvidenceMode::Paired);
        // best_fragments 8/10 clears the platform 0.75 gate but remains below
        // the high-specificity tier's 0.9 threshold, isolating the half-cutoff
        // amount from that independent rescue.
        let too_little = evaluate_witness(
            metrics(1.0, 1.0, 0.30),
            taxon_evidence(10, 10, 8, 8),
            policy,
        );
        assert!(!too_little.rescued(), "the cutoff comparison is strict");
        let minority = evaluate_witness(
            metrics(1.0, 1.0, 0.40),
            taxon_evidence(10, 10, 4, 4),
            policy,
        );
        assert!(!minority.rescued());
    }

    #[test]
    fn paired_low_abundance_rescue_accepts_only_coherent_supported_evidence() {
        let policy = thresholds(EvidenceMode::Paired);
        let evidence = taxon_evidence(40, 20, 15, 10);
        let mut coverage = metrics(0.95, 0.60, 0.29);
        coverage.count = 20;
        coverage.depth = 0.001;

        let decision = evaluate_witness(coverage, evidence, policy);
        assert!(!decision.standard);
        assert!(!decision.rescues.platform);
        assert!(!decision.rescues.high_specificity_low_amount);
        assert!(decision.rescues.low_abundance_pe);

        let amount_boundary = policy.safeguard_min_chunk_breadth
            * low_abundance_amount_fraction(policy.safeguard_min_chunk_breadth);
        coverage.fixed_chunk = amount_boundary;
        assert!(
            !evaluate_witness(coverage, evidence, policy)
                .rescues
                .low_abundance_pe,
            "the amount threshold is intentionally strict"
        );

        coverage.fixed_chunk = 0.29;
        coverage.count = 19;
        assert!(!evaluate_witness(coverage, evidence, policy).rescued());
        coverage.count = 20;
        coverage.depth = 0.000_999;
        assert!(!evaluate_witness(coverage, evidence, policy).rescued());
        coverage.depth = 0.001;
        assert!(
            !evaluate_witness(coverage, taxon_evidence(40, 20, 15, 9), policy).rescued(),
            "both the absolute and fractional clear-support gates must hold"
        );
        assert!(
            !evaluate_witness(coverage, taxon_evidence(40, 20, 14, 10), policy).rescued(),
            "resolved-best support below 75% must fail"
        );
    }

    #[test]
    fn scaled_standard_cutoff_does_not_relax_rescue_amount_safeguards() {
        let mut policy = thresholds(EvidenceMode::Paired);
        policy.standard_min_chunk_breadth = 0.42;
        policy.safeguard_min_chunk_breadth = 0.60;

        let standard = evaluate_witness(
            metrics(1.0, 1.0, 0.43),
            taxon_evidence(10, 10, 5, 5),
            policy,
        );
        assert!(standard.standard, "the scaled standard cutoff must apply");

        // 0.25 clears half of the scaled 0.42 cutoff, but not half of the
        // original 0.60 safeguard. It must therefore not enter the ordinary
        // paired rescue tier.
        let below_original_half = evaluate_witness(
            metrics(1.0, 1.0, 0.25),
            taxon_evidence(10, 10, 8, 8),
            policy,
        );
        assert!(!below_original_half.standard);
        assert!(!below_original_half.rescued());

        let above_original_half = evaluate_witness(
            metrics(1.0, 1.0, 0.31),
            taxon_evidence(10, 10, 8, 8),
            policy,
        );
        assert!(above_original_half.rescues.platform);

        // The low-abundance amount curve must likewise be evaluated from the
        // original 0.60 cutoff, not the scaled standard cutoff.
        let evidence = taxon_evidence(40, 20, 15, 10);
        let mut low_abundance = metrics(0.95, 0.60, 0.25);
        low_abundance.count = 20;
        low_abundance.depth = 0.001;
        assert!(!evaluate_witness(low_abundance, evidence, policy).rescued());
        low_abundance.fixed_chunk = 0.29;
        assert!(
            evaluate_witness(low_abundance, evidence, policy)
                .rescues
                .low_abundance_pe
        );
    }

    #[test]
    fn decisions_attach_the_mode_and_tier_specific_premap_guard() {
        let paired_standard = evaluate_witness(
            metrics(1.0, 1.0, 0.70),
            taxon_evidence(10, 10, 5, 5),
            thresholds(EvidenceMode::Paired),
        );
        assert!(paired_standard.accepted());
        assert!(paired_standard.standard);
        assert_eq!(
            paired_standard.required_premap_guard,
            PremapChunkGuard::Half
        );

        let paired_platform = evaluate_witness(
            metrics(0.10, 1.0, 0.40),
            taxon_evidence(10, 10, 8, 8),
            thresholds(EvidenceMode::Paired),
        );
        assert!(paired_platform.rescues.platform);
        assert_eq!(
            paired_platform.required_premap_guard,
            PremapChunkGuard::Half
        );

        let single_standard = evaluate_witness(
            metrics(1.0, 1.0, 0.70),
            taxon_evidence(10, 10, 5, 5),
            thresholds(EvidenceMode::Single),
        );
        assert_eq!(
            single_standard.required_premap_guard,
            PremapChunkGuard::Full
        );

        let single_rescue = evaluate_witness(
            metrics(0.10, 1.0, 0.25),
            taxon_evidence(20, 20, 17, 17),
            thresholds(EvidenceMode::Single),
        );
        assert!(single_rescue.rescues.platform);
        assert_eq!(
            single_rescue.required_premap_guard,
            PremapChunkGuard::None,
            "SE MAP-only additions remain standard-tier only"
        );

        let long_rescue = evaluate_witness(
            metrics(1.0, 1.0, 0.40),
            taxon_evidence(10, 10, 5, 5),
            thresholds(EvidenceMode::Long),
        );
        assert!(long_rescue.rescues.platform);
        assert_eq!(
            long_rescue.required_premap_guard,
            PremapChunkGuard::None,
            "guarded MAP is disabled for long reads"
        );
    }

    #[test]
    fn low_abundance_reason_wins_when_multiple_rescues_match() {
        let mut coverage = metrics(0.95, 1.0, 0.40);
        coverage.count = 20;
        coverage.depth = 0.01;
        let decision = evaluate_witness(
            coverage,
            taxon_evidence(40, 20, 15, 10),
            thresholds(EvidenceMode::Paired),
        );
        assert!(!decision.standard);
        assert!(decision.rescues.platform);
        assert!(decision.rescues.low_abundance_pe);
        assert_eq!(
            decision.required_premap_guard,
            PremapChunkGuard::FivePercent,
            "the independent low-abundance reason retains its calibrated guard"
        );
    }

    #[test]
    fn paired_low_abundance_rescue_honors_ratio_overrides_and_mode_scope() {
        let mut coverage = metrics(0.95, 0.60, 0.29);
        coverage.count = 20;
        coverage.depth = 0.001;
        let evidence = taxon_evidence(40, 20, 15, 10);
        let paired = thresholds(EvidenceMode::Paired);
        assert!(
            evaluate_witness(coverage, evidence, paired)
                .rescues
                .low_abundance_pe
        );

        coverage.oebr = 0.949;
        assert!(!evaluate_witness(coverage, evidence, paired).rescued());
        coverage.oebr = 0.95;
        coverage.coebr = 0.599;
        assert!(!evaluate_witness(coverage, evidence, paired).rescued());
        coverage.coebr = 0.60;

        let mut explicit_oebr = paired;
        explicit_oebr.min_oebr = Some(0.96);
        assert!(!evaluate_witness(coverage, evidence, explicit_oebr).rescued());
        let mut explicit_coebr = paired;
        explicit_coebr.min_coebr = Some(0.61);
        assert!(!evaluate_witness(coverage, evidence, explicit_coebr).rescued());

        assert!(
            !evaluate_witness(coverage, evidence, thresholds(EvidenceMode::Single)).rescued(),
            "the low-abundance extension is not calibrated for single-end reads"
        );
        assert!(
            !evaluate_witness(coverage, evidence, thresholds(EvidenceMode::Long)).rescued(),
            "the low-abundance extension is not calibrated for long reads"
        );
        assert!(
            !evaluate_ranked_witness(coverage, evidence, paired, ProfileRank::Genus,).rescued(),
            "higher-rank terminal references use only the standard rule"
        );
    }

    #[test]
    fn paired_low_abundance_rescue_honors_explicit_strict_mode() {
        let mut coverage = metrics(0.95, 0.60, 0.29);
        coverage.count = 20;
        coverage.depth = 0.001;
        coverage.breadth_pvalue = Some(1e-10);
        coverage.chunk_pvalue = Some(1e-10);
        let evidence = taxon_evidence(40, 20, 15, 10);

        assert!(
            evaluate_witness(coverage, evidence, thresholds(EvidenceMode::Paired))
                .rescues
                .low_abundance_pe
        );
        let mut strict = thresholds(EvidenceMode::Paired);
        strict.strict = true;
        assert!(!evaluate_witness(coverage, evidence, strict).rescued());
    }

    #[test]
    fn low_abundance_amount_curve_is_bounded_and_monotonic() {
        let values = [0.0, 0.10, 0.60, 1.0].map(low_abundance_amount_fraction);
        assert_eq!(values[0], 0.20);
        assert!((values[1] - 0.245).abs() < 1e-12);
        assert!((values[2] - 0.47).abs() < 1e-12);
        assert_eq!(values[3], 0.50);
        assert!(values.windows(2).all(|pair| pair[0] <= pair[1]));
    }

    #[test]
    fn short_read_species_uses_070_oebr_and_075_coebr_defaults() {
        // OEBR 0.72 clears the species OEBR floor while COEBR 0.76 clears the
        // independently calibrated short-read COEBR floor.
        let decision = evaluate_witness(
            metrics(0.72, 0.76, 0.80),
            taxon_evidence(10, 10, 5, 5),
            thresholds(EvidenceMode::Paired),
        );
        assert!(decision.standard);

        let low_coebr = evaluate_witness(
            metrics(0.72, 0.72, 0.80),
            taxon_evidence(10, 10, 5, 5),
            thresholds(EvidenceMode::Paired),
        );
        assert!(!low_coebr.standard);
    }

    #[test]
    fn coverage_rescue_admits_below_cutoff_high_specificity_species() {
        // Tier-C: fixed_chunk 0.20 is below the full cutoff (0.60) and below the
        // platform half-cutoff (0.30), so neither standard nor platform rescue
        // applies — but best_fraction 0.95 with real depth/count and in-band
        // ratios qualifies it for the high-specificity low-amount rescue.
        let decision = evaluate_witness(
            metrics(1.0, 1.0, 0.20),
            taxon_evidence(20, 20, 19, 19),
            thresholds(EvidenceMode::Paired),
        );
        assert!(!decision.standard);
        assert!(decision.rescues.high_specificity_low_amount);
    }

    #[test]
    fn coverage_rescue_requires_high_specificity_depth_and_count() {
        // fixed_chunk 0.20 lies below both cutoffs, so any acceptance here comes
        // only from tier-C. Each variant removes one of its gates.
        let base = thresholds(EvidenceMode::Paired);
        let ev = |cand, best| taxon_evidence(cand, cand, best, best);
        // best_fraction 0.85 < 0.9 -> no coverage rescue.
        assert!(!evaluate_witness(metrics(1.0, 1.0, 0.20), ev(20, 17), base).rescued());
        // coverage depth below 0.01 -> no coverage rescue.
        let mut shallow = metrics(1.0, 1.0, 0.20);
        shallow.depth = 0.005;
        assert!(!evaluate_witness(shallow, ev(20, 19), base).rescued());
        // fewer than three fragments -> no coverage rescue.
        let mut sparse = metrics(1.0, 1.0, 0.20);
        sparse.count = 2;
        assert!(!evaluate_witness(sparse, ev(20, 19), base).rescued());
    }

    #[test]
    fn coverage_rescue_honors_explicit_oebr_override() {
        // OEBR/COEBR 0.72 clear the default 0.70 species band, and fixed_chunk
        // 0.20 sits below both the full (0.60) and half (0.30) cutoffs, so the
        // only possible acceptance is tier-C.
        let cov = metrics(0.72, 0.72, 0.20);
        let ev = taxon_evidence(20, 20, 19, 19);
        // Explicitly lower COEBR to isolate the OEBR override from the 0.75
        // short-read COEBR default.
        let mut default = thresholds(EvidenceMode::Paired);
        default.min_coebr = Some(0.70);
        // Default OEBR band (no explicit OEBR override): tier-C admits it.
        let admitted = evaluate_witness(cov, ev, default);
        assert!(!admitted.standard);
        assert!(
            admitted.rescues.high_specificity_low_amount,
            "0.72 clears the default 0.70 tier-C OEBR band"
        );
        // Explicit --min-oebr 0.85 must be honored by tier-C, not just standard.
        let mut strict = default;
        strict.min_oebr = Some(0.85);
        let rejected = evaluate_witness(cov, ev, strict);
        assert!(!rejected.standard);
        assert!(
            !rejected.rescued(),
            "explicit --min-oebr 0.85 must reject the tier-C rescue at 0.72"
        );
    }

    #[test]
    fn coverage_rescue_uses_direct_clear_winners_not_map_only_support() {
        // Defensive distinction: 19/20 fragments have resolved-best support
        // (which may include guarded-MAP exact ties), but only 4 are direct
        // unique or margin-resolved winners. fixed_chunk 0.20 is below both
        // cutoffs, so tier-C is the only possible path.
        let ambiguous = evaluate_witness(
            metrics(1.0, 1.0, 0.20),
            taxon_evidence(20, 20, 19, 4), // clear_best_fraction 0.20 must not
            thresholds(EvidenceMode::Paired),
        );
        assert!(!ambiguous.standard);
        assert!(
            !ambiguous.rescued(),
            "tier-C must key off direct clear winners, not MAP-only support"
        );
        // The same taxon with genuine unambiguous wins is still rescued.
        let resolved = evaluate_witness(
            metrics(1.0, 1.0, 0.20),
            taxon_evidence(20, 20, 19, 19),
            thresholds(EvidenceMode::Paired),
        );
        assert!(
            resolved.rescues.high_specificity_low_amount,
            "unambiguous wins restore the tier-C rescue"
        );
    }

    #[test]
    fn standard_witness_is_preferred_before_fixed_chunk_size() {
        let standard = WitnessDecision {
            metrics: metrics(1.0, 1.0, 0.61),
            standard: true,
            rescues: RescueKinds {
                platform: true,
                ..RescueKinds::default()
            },
            required_premap_guard: PremapChunkGuard::None,
        };
        let larger_rescue = WitnessDecision {
            metrics: metrics(0.1, 1.0, 0.90),
            standard: false,
            rescues: RescueKinds {
                platform: true,
                ..RescueKinds::default()
            },
            required_premap_guard: PremapChunkGuard::None,
        };
        assert!(witness_is_better(standard, larger_rescue));
        assert!(!witness_is_better(larger_rescue, standard));
    }

    #[test]
    fn rescue_witness_is_preferred_over_larger_rejected_witness() {
        let rescue = WitnessDecision {
            metrics: metrics(1.0, 1.0, 0.40),
            standard: false,
            rescues: RescueKinds {
                platform: true,
                ..RescueKinds::default()
            },
            required_premap_guard: PremapChunkGuard::None,
        };
        let larger_rejected = WitnessDecision {
            metrics: metrics(0.5, 0.5, 0.90),
            standard: false,
            rescues: RescueKinds::default(),
            required_premap_guard: PremapChunkGuard::None,
        };
        assert!(witness_is_better(rescue, larger_rejected));
        assert!(!witness_is_better(larger_rejected, rescue));
    }

    #[test]
    fn explicit_strict_mode_also_applies_to_rescue() {
        let mut metric = metrics(0.1, 1.0, 0.40);
        metric.breadth_pvalue = Some(1e-8);
        let mut policy = thresholds(EvidenceMode::Paired);
        policy.strict = true;
        let decision = evaluate_witness(metric, taxon_evidence(10, 10, 10, 10), policy);
        assert!(!decision.rescued());
    }
}

use std::{
    collections::{HashMap, HashSet},
    fmt,
    hash::{BuildHasherDefault, Hasher},
    sync::Arc,
};

use anyhow::{bail, Result};
use serde::Serialize;

use crate::{
    index::{Index, LoadedIndex, ReferenceView, SeedLookupMode, SeedOccurrence},
    seeds::seeds_with_strobe_start_distance,
    sequence::{base_code, projected_overlap, revcomp, revcomp_into},
};

pub struct RuntimeIndex {
    index: LoadedIndex,
    ref_lookup: HashMap<String, usize>,
    lookup_mode: SeedLookupMode,
    reference_candidate_groups: Option<Box<[ReferenceCandidateGroups]>>,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct RuntimeIndexOptions {
    pub allow_single_strobe_matches: bool,
}

impl RuntimeIndex {
    pub fn new(index: Index) -> Self {
        Self::from_loaded(LoadedIndex::Owned(index))
    }

    pub fn from_loaded(index: LoadedIndex) -> Self {
        Self::from_loaded_with_options(index, RuntimeIndexOptions::default())
    }

    pub fn from_loaded_with_options(index: LoadedIndex, options: RuntimeIndexOptions) -> Self {
        let ref_lookup = (0..index.reference_count())
            .map(|i| (index.reference_id(i as u32).to_string(), i))
            .collect();
        let lookup_mode = if options.allow_single_strobe_matches {
            SeedLookupMode::IgnoreRandstrobeFingerprint
        } else {
            SeedLookupMode::IndexDefault
        };
        Self {
            index,
            ref_lookup,
            lookup_mode,
            reference_candidate_groups: None,
        }
    }

    /// Install caller-provided taxonomic grouping for bounded candidate
    /// retention. The table is indexed by the mapping index's numeric
    /// reference ID and is therefore validated once, outside the read hot
    /// path. Without a table, every reference is treated as its own group.
    pub fn set_reference_candidate_groups(
        &mut self,
        groups: Vec<ReferenceCandidateGroups>,
    ) -> Result<()> {
        if groups.len() != self.index.reference_count() {
            bail!(
                "reference candidate-group count {} does not match index reference count {}",
                groups.len(),
                self.index.reference_count()
            );
        }
        self.reference_candidate_groups = Some(groups.into_boxed_slice());
        Ok(())
    }

    pub fn lookup(&self, hash: u64) -> crate::index::SeedOccurrences<'_> {
        self.index.lookup_with_mode(hash, self.lookup_mode)
    }

    fn header(&self) -> &crate::index::IndexHeader {
        self.index.header()
    }

    fn reference(&self, ref_id: u32) -> ReferenceView<'_> {
        self.index.reference(ref_id)
    }

    fn reference_id(&self, ref_id: u32) -> &str {
        self.index.reference_id(ref_id)
    }

    fn reference_index_by_id(&self, ref_id: &str) -> Option<u32> {
        self.ref_lookup.get(ref_id).map(|i| *i as u32)
    }

    #[inline]
    fn reference_candidate_groups(&self, ref_id: u32) -> ReferenceCandidateGroups {
        self.reference_candidate_groups
            .as_deref()
            .map(|groups| groups[ref_id as usize])
            .unwrap_or_else(|| ReferenceCandidateGroups::for_reference(ref_id))
    }
}

/// Taxonomic endpoints used only to diversify the bounded mapper output.
///
/// `profile_taxon` is the normal community endpoint (species when known,
/// otherwise the lowest known higher rank). `terminal_taxon` is the most
/// specific declared endpoint (strain/subspecies when known, otherwise the
/// profile endpoint). Keeping both makes saved mappings suitable for either
/// community or strain profiling without increasing the number of rows.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct ReferenceCandidateGroups {
    pub profile_taxon: u32,
    pub terminal_taxon: u32,
}

impl ReferenceCandidateGroups {
    pub const fn new(profile_taxon: u32, terminal_taxon: u32) -> Self {
        Self {
            profile_taxon,
            terminal_taxon,
        }
    }

    const fn for_reference(ref_id: u32) -> Self {
        Self::new(ref_id, ref_id)
    }
}

/// Which taxonomic candidate sets may be incomplete at a configured mapper
/// boundary. The two bits share one serialized byte in saved mappings.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct CandidateBoundary(u8);

impl CandidateBoundary {
    const PROFILE_BIT: u8 = 1;
    const TERMINAL_BIT: u8 = 2;
    const VALID_BITS: u8 = Self::PROFILE_BIT | Self::TERMINAL_BIT;

    pub const COMPLETE: Self = Self(0);

    pub const fn from_bits(bits: u8) -> Option<Self> {
        if bits & !Self::VALID_BITS == 0 {
            Some(Self(bits))
        } else {
            None
        }
    }

    pub const fn bits(self) -> u8 {
        self.0
    }

    pub const fn profile_incomplete(self) -> bool {
        self.0 & Self::PROFILE_BIT != 0
    }

    pub const fn terminal_incomplete(self) -> bool {
        self.0 & Self::TERMINAL_BIT != 0
    }

    pub const fn is_complete(self) -> bool {
        self.0 == 0
    }

    const fn with_profile_incomplete(self) -> Self {
        Self(self.0 | Self::PROFILE_BIT)
    }

    const fn with_terminal_incomplete(self) -> Self {
        Self(self.0 | Self::TERMINAL_BIT)
    }

    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }
}

#[derive(Clone, Debug)]
pub struct Mapping {
    pub read_id: String,
    pub ref_id: String,
    pub ref_start: u32,
    pub ref_end: u32,
    pub identity: f64,
    pub read_len: usize,
    pub mapped_bases: usize,
    pub mate_ref_status: MateRefStatus,
    pub strand: i8,
    pub seed_count: usize,
    pub score: i64,
    pub second_score: i64,
    pub rank: usize,
    pub num_mappings: usize,
    pub matched_bases: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub pair_score: i64,
    pub adjusted_score: i64,
    pub template_length: i64,
    pub pair_status: PairStatus,
    /// Mate number in paired-end output: 0 for unpaired reads, 1 for R1, 2 for R2.
    pub mate_number: u8,
    /// Read-local identifier of an explicitly selected paired placement.
    /// Zero means that this mapping is not part of a selected pair.
    pub pair_id: u32,
    /// Internal grouping carried through mate rescue and final output
    /// limiting. It is not serialized in mapping TSV rows.
    pub candidate_groups: ReferenceCandidateGroups,
    /// Taxonomic candidate levels that may be incomplete after configured
    /// evaluation/output limits.
    pub candidate_boundary: CandidateBoundary,
}

#[derive(Clone, Copy, Debug)]
pub enum IdentityMode {
    Gapped,
    Ungapped,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LongReadExtension {
    None,
    Candidate,
    Retained,
}

#[derive(Clone, Copy, Debug)]
pub struct MapParams {
    pub top_n: usize,
    pub score_ratio: f64,
    pub min_seeds: usize,
    pub min_identity: f64,
    pub min_mapped_frac: f64,
    /// Absolute mapped-length floor (bases) for the long-read length-OR-fraction
    /// coverage gate. `Some(n)` (Ont preset) admits a hit that clears `n` mapped
    /// bases even when it falls below `min_mapped_frac`; `None` (Illumina) keeps
    /// the fraction-only gate.
    pub min_mapped_len: Option<usize>,
    pub no_identity_filter: bool,
    pub ignore_n: bool,
    pub identity_mode: IdentityMode,
    pub diag_bin_size: i64,
    pub extension_max_mismatches: usize,
    pub full_read_interval: bool,
    pub fast_long_read_identity: bool,
    pub long_read_chaining: bool,
    pub long_read_extension: LongReadExtension,
    pub query_seed_stride: usize,
    pub long_read_rescue: bool,
    pub long_read_rescue_limit: Option<usize>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct MapProfile {
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub emitted_records: u64,
    pub total_query_seeds: u64,
    /// Occurrences returned by index lookup after any full-randstrobe verification.
    pub total_seed_hits: u64,
    pub candidate_groups: u64,
    pub candidates_before_filtering: u64,
    pub candidates_after_filtering: u64,
    pub identity_comparisons: u64,
    pub elapsed_seconds: f64,
    pub reads_per_second: f64,
}

/// Reusable per-worker scratch for mapping a sequence.
///
/// The public `map_read` convenience wrapper creates this scratch on demand.
/// High-throughput callers should use `map_read_with_workspace` and retain one
/// workspace per worker so long-read edit-distance buffers are not allocated
/// again for every candidate mapping.
#[derive(Default)]
pub struct MappingWorkspace {
    edit_distance: EditDistanceWorkspace,
    candidate_selection: CandidateSelectionWorkspace,
    mate_rescue: MateRescueWorkspace,
}

type FastU32Set = HashSet<u32, BuildHasherDefault<FastHasher>>;

#[derive(Default)]
struct CandidateSelectionWorkspace {
    selected: Vec<bool>,
    profile_groups: FastU32Set,
    terminal_groups: FastU32Set,
}

#[derive(Clone, Copy, Debug)]
struct RescueSeedHit {
    query_pos: u32,
    occurrence: SeedOccurrence,
}

/// Reusable seed and grouping storage for paired mate rescue.
///
/// A mate can have many candidate placements. Rescanning all query seeds and
/// repeating every index lookup once per placement is especially expensive
/// with `--min-seeds 1`. The hit lists below are built once for the rescued
/// read, then filtered into each candidate mate window without changing the
/// original window-specific grouping or scoring.
#[derive(Default)]
struct MateRescueWorkspace {
    reverse_complement: Vec<u8>,
    forward_hits: Vec<RescueSeedHit>,
    reverse_hits: Vec<RescueSeedHit>,
}

#[derive(Default)]
struct EditDistanceWorkspace {
    /// Query-position masks for A, C, G, T, plus an intentionally empty
    /// ambiguous-reference mask.
    peq: [Vec<u64>; 5],
    positive: Vec<u64>,
    negative: Vec<u64>,
}

impl MapProfile {
    pub fn add(&mut self, other: &Self) {
        self.total_reads += other.total_reads;
        self.mapped_reads += other.mapped_reads;
        self.emitted_records += other.emitted_records;
        self.total_query_seeds += other.total_query_seeds;
        self.total_seed_hits += other.total_seed_hits;
        self.candidate_groups += other.candidate_groups;
        self.candidates_before_filtering += other.candidates_before_filtering;
        self.candidates_after_filtering += other.candidates_after_filtering;
        self.identity_comparisons += other.identity_comparisons;
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct CandidateKey {
    ref_id: u32,
    strand: i8,
    diag_bin: i64,
}

type CandidateMap = HashMap<CandidateKey, CandidateAccumulator, BuildHasherDefault<FastHasher>>;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct ChainKey {
    ref_id: u32,
    strand: i8,
}

type ChainMap = HashMap<ChainKey, Vec<Anchor>, BuildHasherDefault<FastHasher>>;
type Interval = (usize, usize, usize, usize);
type FastExtensionChoice = (Interval, Identity, f64, i64);

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Anchor {
    q: u32,
    r: u32,
}

#[derive(Default)]
struct FastHasher(u64);

impl FastHasher {
    #[inline]
    fn add(&mut self, value: u64) {
        self.0 ^= value.wrapping_add(0x9E37_79B9_7F4A_7C15);
        self.0 = self.0.rotate_left(27).wrapping_mul(0x94D0_49BB_1331_11EB);
    }
}

impl Hasher for FastHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0 ^ (self.0 >> 32)
    }

    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for chunk in &mut chunks {
            self.add(u64::from_ne_bytes(
                chunk.try_into().expect("8-byte hash chunk"),
            ));
        }
        let mut tail = 0u64;
        for (shift, byte) in chunks.remainder().iter().enumerate() {
            tail |= u64::from(*byte) << (shift * 8);
        }
        self.add(tail);
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(u64::from(i));
    }

    #[inline]
    fn write_i8(&mut self, i: i8) {
        self.add(i as i64 as u64);
    }

    #[inline]
    fn write_i64(&mut self, i: i64) {
        self.add(i as u64);
    }
}

#[derive(Clone, Debug)]
struct CandidateAccumulator {
    seed_count: usize,
    query_seed_start: usize,
    query_seed_end: usize,
    ref_seed_start: usize,
    ref_seed_end: usize,
    diag_sum: i128,
    diag_min: i64,
    diag_max: i64,
    last_q: Option<u32>,
}

#[derive(Clone, Debug)]
struct Candidate {
    ref_id: u32,
    strand: i8,
    diag: i64,
    seed_count: usize,
    score: i64,
    query_seed_start: usize,
    query_seed_end: usize,
    ref_seed_start: usize,
    ref_seed_end: usize,
}

#[derive(Clone, Copy, Debug)]
struct Identity {
    matches: usize,
    // Query/read bases compared or mapped. This is separate from the identity
    // denominator, which may include indel bases for gapped estimates.
    compared: usize,
    denom: usize,
    value: f64,
}

const LONG_READ_PRIMARY_EVAL_MIN: usize = 32;
const LONG_READ_PRIMARY_EVAL_MAX: usize = 128;
const LONG_READ_RESCUE_EVAL_MIN: usize = 128;
const LONG_READ_RESCUE_EVAL_MAX: usize = 512;
const LONG_READ_RESCUE_IDENTITY: f64 = 0.92;
const LONG_READ_EXTEND_XDROP: i32 = 80;
const LONG_READ_EXTEND_MAX_SIDE: usize = 10_000;
const FAST_EXTEND_MIN_FLANK_DENOM_FOR_FILTER: usize = 64;
// Minimum identity required of the *added* flank for a fast extension to be
// accepted. Anchored at the banded x-drop break-even point: with match +2 /
// mismatch -3 the alignment drifts positive above ~0.6 identity, so --extend
// keeps extending through noisy-but-homologous read ends down to ~60% identity
// and stops at non-homologous sequence. Fast-extend uses the same floor so it
// completes the same regions --extend would, while still rejecting
// non-homologous tails (random DNA is ~0.25 identity). Crucially this is
// decoupled from the core identity, so a noisy flank on a high-identity core is
// no longer rejected just for being noisier than the core.
const FAST_EXTEND_MIN_FLANK_IDENTITY: f64 = 0.60;

#[derive(Clone, Copy, Debug, Default)]
struct DpCell {
    score: i32,
    matches: usize,
    compared: usize,
    denom: usize,
}

impl CandidateAccumulator {
    fn new(q: u32, diag: i64, k: usize) -> Self {
        Self {
            seed_count: 1,
            query_seed_start: q as usize,
            query_seed_end: q as usize + k,
            ref_seed_start: (diag + i64::from(q)).max(0) as usize,
            ref_seed_end: (diag + i64::from(q)).max(0) as usize + k,
            diag_sum: diag as i128,
            diag_min: diag,
            diag_max: diag,
            last_q: Some(q),
        }
    }

    fn add(&mut self, q: u32, diag: i64, k: usize) {
        if self.last_q == Some(q) {
            return;
        }
        self.seed_count += 1;
        self.query_seed_start = self.query_seed_start.min(q as usize);
        self.query_seed_end = self.query_seed_end.max(q as usize + k);
        let r = (diag + i64::from(q)).max(0) as usize;
        self.ref_seed_start = self.ref_seed_start.min(r);
        self.ref_seed_end = self.ref_seed_end.max(r + k);
        self.diag_sum += diag as i128;
        self.diag_min = self.diag_min.min(diag);
        self.diag_max = self.diag_max.max(diag);
        self.last_q = Some(q);
    }

    fn into_candidate(self, key: CandidateKey, k: usize, min_seeds: usize) -> Option<Candidate> {
        if self.seed_count < min_seeds {
            return None;
        }
        let diag = (self.diag_sum / self.seed_count as i128) as i64;
        let spread = self.diag_max - self.diag_min;
        let span = self.query_seed_end.saturating_sub(self.query_seed_start);
        let score = (self.seed_count * k + span) as i64 - spread.abs();
        Some(Candidate {
            ref_id: key.ref_id,
            strand: key.strand,
            diag,
            seed_count: self.seed_count,
            score,
            query_seed_start: self.query_seed_start,
            query_seed_end: self.query_seed_end,
            ref_seed_start: self.ref_seed_start,
            ref_seed_end: self.ref_seed_end,
        })
    }
}

pub fn map_read(
    read_id: &str,
    read: &[u8],
    index: &RuntimeIndex,
    params: MapParams,
    profile: &mut MapProfile,
) -> Vec<Mapping> {
    let mut workspace = MappingWorkspace::default();
    map_read_with_workspace(read_id, read, index, params, profile, &mut workspace)
}

pub fn map_read_with_workspace(
    read_id: &str,
    read: &[u8],
    index: &RuntimeIndex,
    params: MapParams,
    profile: &mut MapProfile,
    workspace: &mut MappingWorkspace,
) -> Vec<Mapping> {
    if read.is_empty() || params.top_n == 0 {
        return Vec::new();
    }

    let (rc, mut candidates) = collect_candidates(read, index, params, profile, |_| true);

    profile.candidates_before_filtering += candidates.len() as u64;
    if candidates.is_empty() {
        return Vec::new();
    }
    candidates.sort_unstable_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| a.ref_id.cmp(&b.ref_id))
            .then_with(|| a.diag.cmp(&b.diag))
    });
    let best_score = candidates[0].score.max(1) as f64;
    let second_score = candidates.get(1).map(|c| c.score).unwrap_or(0);

    let mut mappings = Vec::new();
    let mut evaluated = vec![false; candidates.len()];
    let primary_limit = candidate_eval_limit(params);
    let primary_indices =
        select_score_ranked_candidates(&candidates, best_score, params.score_ratio, primary_limit);
    let mut candidate_evaluation_limited = score_ranked_candidate_limit_reached(
        &candidates,
        best_score,
        params.score_ratio,
        primary_indices.len(),
    );
    let mut rescue_evaluates_beyond_score_ratio = false;
    for index_in_candidates in primary_indices {
        evaluated[index_in_candidates] = true;
        if let Some(mapping) = evaluate_candidate(
            read_id,
            read,
            rc.as_slice(),
            &candidates[index_in_candidates],
            index,
            params,
            second_score,
            profile,
            workspace,
        ) {
            mappings.push(mapping);
        }
    }

    if should_rescue_long_read(&mappings, params) {
        rescue_evaluates_beyond_score_ratio = true;
        let rescue_limit = long_read_rescue_eval_limit(params).min(candidates.len());
        let mut evaluated_count = evaluated.iter().filter(|seen| **seen).count();
        for (index_in_candidates, candidate) in candidates.iter().enumerate() {
            if evaluated_count >= rescue_limit {
                break;
            }
            if evaluated[index_in_candidates] {
                continue;
            }
            evaluated[index_in_candidates] = true;
            evaluated_count += 1;
            if let Some(mapping) = evaluate_candidate(
                read_id,
                read,
                rc.as_slice(),
                candidate,
                index,
                params,
                second_score,
                profile,
                workspace,
            ) {
                mappings.push(mapping);
            }
        }
        // Rescue continues beyond the primary cap. If it exhausts all
        // candidates, the earlier primary boundary is no longer incomplete.
        candidate_evaluation_limited =
            evaluated_count >= rescue_limit && evaluated.iter().any(|seen| !seen);
    }
    mappings.sort_unstable_by(mapping_order);
    mappings.dedup_by(|a, b| {
        a.ref_id == b.ref_id
            && a.ref_start == b.ref_start
            && a.ref_end == b.ref_end
            && a.strand == b.strand
    });
    let mut candidate_boundary = if candidate_evaluation_limited {
        unevaluated_candidate_boundary(
            &candidates,
            &evaluated,
            &mappings,
            index,
            (!rescue_evaluates_beyond_score_ratio).then_some(best_score * params.score_ratio),
            &mut workspace.candidate_selection,
        )
    } else {
        CandidateBoundary::COMPLETE
    };
    candidate_boundary = candidate_boundary.union(retain_taxon_diverse_mappings(
        &mut mappings,
        params.top_n,
        &mut workspace.candidate_selection,
    ));
    if params.long_read_extension == LongReadExtension::Retained {
        extend_retained_long_read_mappings(
            read,
            rc.as_slice(),
            index,
            params,
            profile,
            &mut mappings,
            workspace,
        );
        dedup_mappings_preserve_order(&mut mappings);
        // Extension re-scores each mapping in place, so the pre-extension
        // ordering may no longer be best-first. Re-sort so the reported
        // primary (rank 1) is the highest-scoring reference after extension.
        mappings.sort_unstable_by(mapping_order);
    }
    if !candidate_boundary.is_complete() {
        for mapping in &mut mappings {
            mapping.candidate_boundary = mapping.candidate_boundary.union(candidate_boundary);
        }
    }
    profile.candidates_after_filtering += mappings.len() as u64;
    set_ranks(&mut mappings);
    mappings
}

fn dedup_mappings_preserve_order(mappings: &mut Vec<Mapping>) {
    let mut seen = HashMap::with_capacity(mappings.len());
    let mut retained = Vec::with_capacity(mappings.len());
    for mapping in mappings.drain(..) {
        let key = (
            mapping.ref_id.clone(),
            mapping.ref_start,
            mapping.ref_end,
            mapping.strand,
        );
        if let Some(&index) = seen.get(&key) {
            if is_better_duplicate(&mapping, &retained[index]) {
                retained[index] = mapping;
            }
        } else {
            seen.insert(key, retained.len());
            retained.push(mapping);
        }
    }
    *mappings = retained;
}

fn is_better_duplicate(candidate: &Mapping, current: &Mapping) -> bool {
    candidate
        .mapped_bases
        .cmp(&current.mapped_bases)
        .then_with(|| interval_span_mapping(candidate).cmp(&interval_span_mapping(current)))
        .then_with(|| candidate.identity.total_cmp(&current.identity))
        .then_with(|| candidate.score.cmp(&current.score))
        .is_gt()
}

fn interval_span_mapping(mapping: &Mapping) -> usize {
    (mapping.ref_end as usize).saturating_sub(mapping.ref_start as usize)
}

fn interval_ref_span(interval: Interval) -> usize {
    interval.3.saturating_sub(interval.2)
}

fn interval_query_span(interval: Interval) -> usize {
    interval.1.saturating_sub(interval.0)
}

#[allow(clippy::too_many_arguments)]
fn evaluate_candidate(
    read_id: &str,
    read: &[u8],
    rc: &[u8],
    candidate: &Candidate,
    index: &RuntimeIndex,
    params: MapParams,
    second_score: i64,
    profile: &mut MapProfile,
    workspace: &mut MappingWorkspace,
) -> Option<Mapping> {
    let reference = index.reference(candidate.ref_id);
    let query = if candidate.strand == 1 { read } else { rc };
    let interval = if params.long_read_chaining {
        chain_interval(query.len(), reference.len(), candidate, false)
    } else {
        let (seed_start, seed_end) = if params.full_read_interval {
            (0, query.len())
        } else {
            (candidate.query_seed_start, candidate.query_seed_end)
        };
        refine_interval(
            query,
            &reference,
            candidate.diag,
            seed_start,
            seed_end,
            params.extension_max_mismatches,
        )
    };
    interval.and_then(|interval| {
        evaluate_candidate_interval(
            read_id,
            read,
            query,
            &reference,
            candidate,
            params,
            second_score,
            interval,
            index,
            profile,
            workspace,
        )
    })
}

#[allow(clippy::too_many_arguments)]
fn evaluate_candidate_interval(
    read_id: &str,
    read: &[u8],
    query: &[u8],
    reference: &ReferenceView<'_>,
    candidate: &Candidate,
    params: MapParams,
    second_score: i64,
    interval: (usize, usize, usize, usize),
    index: &RuntimeIndex,
    profile: &mut MapProfile,
    workspace: &mut MappingWorkspace,
) -> Option<Mapping> {
    let (mut interval, mut ident, mut mapped_frac, mut passes) =
        evaluate_interval_identity(read, query, reference, interval, params, profile, workspace);
    if params.long_read_chaining
        && params.long_read_extension == LongReadExtension::Candidate
        && should_attempt_long_read_extension(ident, passes, params)
    {
        let extended = extend_long_read_interval(query, reference, interval, params.ignore_n);
        if extended != interval {
            let (extended_interval, extended_ident, extended_mapped_frac, extended_passes) =
                evaluate_interval_identity(
                    read, query, reference, extended, params, profile, workspace,
                );
            if should_use_extended_interval(
                ident,
                passes,
                extended_ident,
                extended_passes,
                extended_mapped_frac,
            ) {
                interval = extended_interval;
                ident = extended_ident;
                mapped_frac = extended_mapped_frac;
                passes = extended_passes;
            }
        }
    }
    if !passes && !length_rescues(&ident, params) {
        return None;
    }
    let score = mapping_score(
        candidate.score,
        ident,
        mapped_frac,
        params.long_read_chaining,
    );
    let (qs, qe, rs, re) = interval;
    Some(Mapping {
        read_id: read_id.to_string(),
        ref_id: index.reference_id(candidate.ref_id).to_string(),
        ref_start: rs as u32,
        ref_end: re as u32,
        identity: ident.value,
        read_len: read.len(),
        mapped_bases: ident.compared,
        mate_ref_status: MateRefStatus::NotApplicable,
        strand: candidate.strand,
        seed_count: candidate.seed_count,
        score,
        second_score,
        rank: 0,
        num_mappings: 0,
        matched_bases: ident.matches,
        query_start: qs,
        query_end: qe,
        pair_score: 0,
        adjusted_score: score,
        template_length: 0,
        pair_status: PairStatus::Single,
        mate_number: 0,
        pair_id: 0,
        candidate_groups: index.reference_candidate_groups(candidate.ref_id),
        candidate_boundary: CandidateBoundary::COMPLETE,
    })
}

fn extend_retained_long_read_mappings(
    read: &[u8],
    rc: &[u8],
    index: &RuntimeIndex,
    params: MapParams,
    profile: &mut MapProfile,
    mappings: &mut [Mapping],
    workspace: &mut MappingWorkspace,
) {
    if !params.long_read_chaining {
        return;
    }
    for mapping in mappings {
        let Some(ref_id) = index.reference_index_by_id(&mapping.ref_id) else {
            continue;
        };
        let reference = index.reference(ref_id);
        let query = if mapping.strand == 1 { read } else { rc };
        let interval = (
            mapping.query_start,
            mapping.query_end,
            mapping.ref_start as usize,
            mapping.ref_end as usize,
        );
        if interval.0 >= interval.1 || interval.2 >= interval.3 {
            continue;
        }
        let projections = fast_extension_projections(query.len(), reference.len(), interval);
        if projections.is_empty() {
            continue;
        }
        let base = Identity {
            matches: mapping.matched_bases,
            compared: mapping.mapped_bases,
            denom: identity_denom(mapping.identity, mapping.matched_bases),
            value: mapping.identity,
        };
        let mut best_extension: Option<FastExtensionChoice> = None;
        for extended in projections {
            if interval_ref_span(extended) < interval_ref_span(interval)
                || interval_query_span(extended) < interval_query_span(interval)
            {
                continue;
            }
            let ident = projected_flank_identity(
                query, &reference, interval, extended, base, params, profile, workspace,
            );
            let mapped_frac = ident.compared as f64 / read.len() as f64;
            let passes = params.no_identity_filter
                || (ident.compared > 0
                    && ident.value >= params.min_identity
                    && mapped_frac >= params.min_mapped_frac);
            if !should_use_fast_extended_interval(base, ident, passes, mapped_frac, params) {
                continue;
            }
            let score =
                long_read_mapping_score(retained_support_score(mapping), ident, mapped_frac);
            if best_extension
                .as_ref()
                .map(|(_, prev_ident, prev_mapped_frac, prev_score)| {
                    score > *prev_score
                        || (score == *prev_score
                            && (ident.value, mapped_frac) > (prev_ident.value, *prev_mapped_frac))
                })
                .unwrap_or(true)
            {
                best_extension = Some((extended, ident, mapped_frac, score));
            }
        }
        let Some((extended, ident, _mapped_frac, score)) = best_extension else {
            continue;
        };
        let (qs, qe, rs, re) = extended;
        mapping.query_start = qs;
        mapping.query_end = qe;
        mapping.ref_start = rs as u32;
        mapping.ref_end = re as u32;
        mapping.identity = ident.value;
        mapping.mapped_bases = ident.compared;
        mapping.matched_bases = ident.matches;
        mapping.score = score;
        mapping.adjusted_score = score;
    }
}

fn evaluate_interval_identity(
    read: &[u8],
    query: &[u8],
    reference: &ReferenceView<'_>,
    interval: (usize, usize, usize, usize),
    params: MapParams,
    profile: &mut MapProfile,
    workspace: &mut MappingWorkspace,
) -> ((usize, usize, usize, usize), Identity, f64, bool) {
    let (qs, qe, rs, re) = interval;
    let ident = identity(
        &query[qs..qe],
        reference,
        rs,
        re,
        params.ignore_n,
        params.identity_mode,
        params.fast_long_read_identity,
        &mut workspace.edit_distance,
    );
    profile.identity_comparisons += 1;
    let mapped_frac = ident.compared as f64 / read.len() as f64;
    let passes = params.no_identity_filter
        || (ident.compared > 0
            && ident.value >= params.min_identity
            && mapped_frac >= params.min_mapped_frac);
    (interval, ident, mapped_frac, passes)
}

fn should_attempt_long_read_extension(ident: Identity, passes: bool, params: MapParams) -> bool {
    passes
        || params.no_identity_filter
        || (ident.compared > 0 && ident.value + 0.05 >= params.min_identity)
}

fn should_use_extended_interval(
    base: Identity,
    base_passes: bool,
    extended: Identity,
    extended_passes: bool,
    extended_mapped_frac: f64,
) -> bool {
    if !extended_passes {
        return false;
    }
    if !base_passes {
        return true;
    }
    extended.compared > base.compared && extended_mapped_frac.is_finite()
}

/// Long-read length rescue: a mapping that misses `min_mapped_frac` is still
/// admitted when its absolute mapped length clears `min_mapped_len` (set only
/// for the Ont preset) and it still meets the identity cutoff. This mirrors the
/// profiler's length-OR-fraction coverage gate so that large, low-fraction
/// long-read hits survive to profiling. Returns false when `min_mapped_len` is
/// `None` (Illumina), leaving the fraction-only gate unchanged.
fn length_rescues(ident: &Identity, params: MapParams) -> bool {
    ident.compared > 0
        && ident.value >= params.min_identity
        && params
            .min_mapped_len
            .is_some_and(|min_len| ident.compared >= min_len)
}

fn should_use_fast_extended_interval(
    base: Identity,
    extended: Identity,
    extended_passes: bool,
    extended_mapped_frac: f64,
    params: MapParams,
) -> bool {
    // Two gates, no longer tied to the core identity:
    //  1. should_use_extended_interval — the extension must add mapped bases and
    //     the extended mapping must still satisfy the read-level identity/coverage
    //     filters (this rejects large non-homologous tails outright).
    //  2. fast_extension_flank_passes — the *added* flank must be clearly
    //     homologous, so we never glue on a short non-homologous tail that slips
    //     under the read-level filter.
    should_use_extended_interval(base, true, extended, extended_passes, extended_mapped_frac)
        && fast_extension_flank_passes(base, extended, params.min_identity)
}

fn fast_extension_flank_passes(base: Identity, extended: Identity, min_identity: f64) -> bool {
    let added_denom = extended.denom.saturating_sub(base.denom);
    if added_denom < FAST_EXTEND_MIN_FLANK_DENOM_FOR_FILTER {
        // Tiny flank: trust the read-level filter rather than a noisy local estimate.
        return true;
    }
    let added_matches = extended.matches.saturating_sub(base.matches);
    let flank_identity = added_matches as f64 / added_denom as f64;
    // Demand a little more on high-identity datasets, but never below the x-drop
    // break-even floor (~0.6). Decoupled from base.value so noisy ONT read ends
    // on a clean core are still accepted.
    let min_flank_identity = (min_identity - 0.30).max(FAST_EXTEND_MIN_FLANK_IDENTITY);
    flank_identity >= min_flank_identity
}

fn extend_long_read_interval(
    query: &[u8],
    reference: &ReferenceView<'_>,
    interval: (usize, usize, usize, usize),
    ignore_n: bool,
) -> (usize, usize, usize, usize) {
    let (mut qs, mut qe, mut rs, mut re) = interval;
    let left_q = qs.min(LONG_READ_EXTEND_MAX_SIDE);
    let left_r = rs.min(left_q.saturating_add(extension_slack(left_q)));
    if left_q > 0 && left_r > 0 {
        let query_left = query[qs - left_q..qs]
            .iter()
            .rev()
            .copied()
            .collect::<Vec<_>>();
        let reference_left = reference_segment(reference, rs - left_r, rs, true);
        let (q_ext, r_ext) = banded_xdrop_extend(
            &query_left,
            &reference_left,
            ignore_n,
            LONG_READ_EXTEND_XDROP,
        );
        qs -= q_ext;
        rs -= r_ext;
    }

    let right_q = query
        .len()
        .saturating_sub(qe)
        .min(LONG_READ_EXTEND_MAX_SIDE);
    let right_r = reference
        .len()
        .saturating_sub(re)
        .min(right_q.saturating_add(extension_slack(right_q)));
    if right_q > 0 && right_r > 0 {
        let query_right = &query[qe..qe + right_q];
        let reference_right = reference_segment(reference, re, re + right_r, false);
        let (q_ext, r_ext) = banded_xdrop_extend(
            query_right,
            &reference_right,
            ignore_n,
            LONG_READ_EXTEND_XDROP,
        );
        qe += q_ext;
        re += r_ext;
    }
    (qs, qe, rs, re)
}

fn project_full_query_interval(
    query_len: usize,
    ref_len: usize,
    interval: (usize, usize, usize, usize),
) -> Option<(usize, usize, usize, usize)> {
    let (qs, qe, rs, re) = interval;
    if qs >= qe || rs >= re || query_len == 0 {
        return None;
    }
    let q_span = qe - qs;
    let r_span = re - rs;
    if q_span == 0 || r_span == 0 {
        return None;
    }
    let right_q = query_len.saturating_sub(qe);
    let left_r = scale_query_to_ref(qs, q_span, r_span).min(qs + fast_projection_slack(qs));
    let right_r =
        scale_query_to_ref(right_q, q_span, r_span).min(right_q + fast_projection_slack(right_q));
    let projected = (
        0,
        query_len,
        rs.saturating_sub(left_r),
        re.saturating_add(right_r).min(ref_len),
    );
    if projected.2 >= projected.3 || projected == interval {
        None
    } else {
        Some(projected)
    }
}

fn fast_extension_projections(
    query_len: usize,
    ref_len: usize,
    interval: (usize, usize, usize, usize),
) -> Vec<(usize, usize, usize, usize)> {
    let mut projections = Vec::with_capacity(2);
    if let Some(projected) = project_full_query_interval(query_len, ref_len, interval) {
        projections.push(projected);
    }
    if let Some(projected) = project_equal_query_interval(query_len, ref_len, interval) {
        if !projections.contains(&projected) {
            projections.push(projected);
        }
    }
    projections
}

fn project_equal_query_interval(
    query_len: usize,
    ref_len: usize,
    interval: (usize, usize, usize, usize),
) -> Option<(usize, usize, usize, usize)> {
    let (qs, qe, rs, re) = interval;
    if qs >= qe || rs >= re || query_len == 0 {
        return None;
    }
    let right_q = query_len.saturating_sub(qe);
    let left_r = qs.min(rs);
    let right_r = right_q.min(ref_len.saturating_sub(re));
    let projected = (0, query_len, rs - left_r, re + right_r);
    if projected.2 >= projected.3 || projected == interval {
        None
    } else {
        Some(projected)
    }
}

fn fast_projection_slack(query_bases: usize) -> usize {
    query_bases / 10 + 64
}

/// Identity of a long-read extension flank. Unlike [`identity`]/[`gapped_identity`],
/// this never takes the `n == m` ungapped shortcut: extension flanks are short but
/// indel-rich (noisy ONT read ends), and an ungapped comparison goes out of frame
/// after the first indel and reports a near-random identity — which made
/// fast-extend reject legitimate flanks and leave the mapped region incomplete.
fn flank_identity(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ref_start: usize,
    ref_end: usize,
    ignore_n: bool,
    mode: IdentityMode,
    fast_long_read_identity: bool,
    workspace: &mut EditDistanceWorkspace,
) -> Identity {
    match mode {
        // The user explicitly opted out of gap handling.
        IdentityMode::Ungapped => ungapped_identity(query, reference, ref_start, ref_end, ignore_n),
        IdentityMode::Gapped => {
            if query.is_empty() || ref_end <= ref_start {
                return ungapped_identity(query, reference, ref_start, ref_end, ignore_n);
            }
            if fast_long_read_identity && query.len().max(ref_end - ref_start) > 512 {
                fast_edit_identity(query, reference, ref_start, ref_end, ignore_n, workspace)
            } else {
                banded_gapped_identity(query, reference, ref_start, ref_end, ignore_n)
            }
        }
    }
}

fn projected_flank_identity(
    query: &[u8],
    reference: &ReferenceView<'_>,
    base_interval: (usize, usize, usize, usize),
    projected_interval: (usize, usize, usize, usize),
    base: Identity,
    params: MapParams,
    profile: &mut MapProfile,
    workspace: &mut MappingWorkspace,
) -> Identity {
    let (base_qs, base_qe, base_rs, base_re) = base_interval;
    let (proj_qs, proj_qe, proj_rs, proj_re) = projected_interval;
    let mut matches = base.matches;
    let mut compared = base.compared;
    let mut denom = base.denom;

    if proj_qs < base_qs || proj_rs < base_rs {
        let q_start = proj_qs.min(base_qs);
        let q_end = base_qs.min(query.len());
        let r_start = proj_rs.min(base_rs);
        let r_end = base_rs.min(reference.len());
        if q_start < q_end || r_start < r_end {
            let ident = flank_identity(
                &query[q_start..q_end],
                reference,
                r_start,
                r_end,
                params.ignore_n,
                params.identity_mode,
                params.fast_long_read_identity,
                &mut workspace.edit_distance,
            );
            profile.identity_comparisons += 1;
            matches += ident.matches;
            compared += ident.compared;
            denom += ident.denom;
        }
    }

    if base_qe < proj_qe || base_re < proj_re {
        let q_start = base_qe.min(query.len());
        let q_end = proj_qe.min(query.len());
        let r_start = base_re.min(reference.len());
        let r_end = proj_re.min(reference.len());
        if q_start < q_end || r_start < r_end {
            let ident = flank_identity(
                &query[q_start..q_end],
                reference,
                r_start,
                r_end,
                params.ignore_n,
                params.identity_mode,
                params.fast_long_read_identity,
                &mut workspace.edit_distance,
            );
            profile.identity_comparisons += 1;
            matches += ident.matches;
            compared += ident.compared;
            denom += ident.denom;
        }
    }

    Identity {
        matches,
        compared: compared.min(query.len()),
        denom,
        value: if denom == 0 {
            0.0
        } else {
            matches as f64 / denom as f64
        },
    }
}

fn identity_denom(identity: f64, matches: usize) -> usize {
    if identity > 0.0 && identity.is_finite() {
        ((matches as f64) / identity).round().max(matches as f64) as usize
    } else {
        matches
    }
}

fn extension_slack(query_bases: usize) -> usize {
    query_bases / 5 + 128
}

fn reference_segment(
    reference: &ReferenceView<'_>,
    start: usize,
    end: usize,
    reverse: bool,
) -> Vec<u8> {
    if reverse {
        (start..end)
            .rev()
            .map(|pos| reference.base_at(pos))
            .collect()
    } else {
        (start..end).map(|pos| reference.base_at(pos)).collect()
    }
}

fn banded_xdrop_extend(
    query: &[u8],
    reference: &[u8],
    ignore_n: bool,
    xdrop: i32,
) -> (usize, usize) {
    let n = query.len();
    let m = reference.len();
    if n == 0 || m == 0 {
        return (0, 0);
    }
    let band = n.abs_diff(m).saturating_add(96).max(n.max(m) / 20);
    let band = band.min(n.max(m)).max(16);
    let impossible = i32::MIN / 4;
    let mut prev_start = 0usize;
    let mut prev_end = m.min(band);
    let mut prev = (prev_start..=prev_end)
        .map(|j| -(j as i32) * 2)
        .collect::<Vec<_>>();
    let mut best_score = 0i32;
    let mut best = (0usize, 0usize);

    for i in 1..=n {
        let start = i.saturating_sub(band);
        let end = i.saturating_add(band).min(m);
        let mut curr = vec![impossible; end - start + 1];
        let mut row_best = impossible;
        for j in start..=end {
            let idx = j - start;
            if j == 0 {
                curr[idx] = -(i as i32) * 2;
            } else {
                let diag = if (prev_start..=prev_end).contains(&(j - 1)) {
                    prev[j - 1 - prev_start]
                        + extension_subst_score(query[i - 1], reference[j - 1], ignore_n)
                } else {
                    impossible
                };
                let deletion = if (prev_start..=prev_end).contains(&j) {
                    prev[j - prev_start] - 2
                } else {
                    impossible
                };
                let insertion = if j > start {
                    curr[idx - 1] - 2
                } else {
                    impossible
                };
                curr[idx] = diag.max(deletion).max(insertion);
            }
            row_best = row_best.max(curr[idx]);
            if curr[idx] > best_score {
                best_score = curr[idx];
                best = (i, j);
            }
        }
        if row_best < best_score - xdrop {
            break;
        }
        prev = curr;
        prev_start = start;
        prev_end = end;
    }
    best
}

fn extension_subst_score(q: u8, r: u8, ignore_n: bool) -> i32 {
    let qb = base_code(q);
    let rb = base_code(r);
    if ignore_n && (qb.is_none() || rb.is_none()) {
        0
    } else if matches!((qb, rb), (Some(a), Some(b)) if a == b) {
        2
    } else {
        -3
    }
}

fn mapping_score(
    candidate_score: i64,
    ident: Identity,
    mapped_frac: f64,
    long_read_chaining: bool,
) -> i64 {
    if !long_read_chaining {
        return candidate_score;
    }
    long_read_mapping_score(candidate_score.clamp(0, 100_000), ident, mapped_frac)
}

fn long_read_mapping_score(support_score: i64, ident: Identity, mapped_frac: f64) -> i64 {
    let identity_score = (ident.value * 10_000_000.0).round() as i64;
    let coverage_score = (mapped_frac * 1_000_000.0).round() as i64;
    identity_score + coverage_score + ident.matches as i64 + support_score
}

fn retained_support_score(mapping: &Mapping) -> i64 {
    let identity_score = (mapping.identity * 10_000_000.0).round() as i64;
    let coverage_score =
        (mapping.mapped_bases as f64 / mapping.read_len as f64 * 1_000_000.0).round() as i64;
    (mapping.score - identity_score - coverage_score - mapping.matched_bases as i64)
        .clamp(0, 100_000)
}

fn mapping_order(a: &Mapping, b: &Mapping) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| b.mapped_bases.cmp(&a.mapped_bases))
        .then_with(|| b.identity.total_cmp(&a.identity))
        .then_with(|| a.ref_id.cmp(&b.ref_id))
}

fn candidate_eval_limit(params: MapParams) -> usize {
    if params.long_read_chaining {
        params
            .top_n
            .saturating_mul(4)
            .clamp(LONG_READ_PRIMARY_EVAL_MIN, LONG_READ_PRIMARY_EVAL_MAX)
    } else {
        params.top_n.saturating_mul(4).max(16)
    }
}

fn long_read_rescue_eval_limit(params: MapParams) -> usize {
    let default = params
        .top_n
        .saturating_mul(24)
        .clamp(LONG_READ_RESCUE_EVAL_MIN, LONG_READ_RESCUE_EVAL_MAX);
    params
        .long_read_rescue_limit
        .map(|limit| default.min(limit.max(1)))
        .unwrap_or(default)
}

fn select_score_ranked_candidates(
    candidates: &[Candidate],
    best_score: f64,
    score_ratio: f64,
    limit: usize,
) -> Vec<usize> {
    let min_score = best_score * score_ratio;
    let mut selected = Vec::with_capacity(limit.min(candidates.len()));
    for (i, candidate) in candidates.iter().enumerate() {
        if (candidate.score as f64) < min_score {
            break;
        }
        selected.push(i);
        if selected.len() >= limit {
            break;
        }
    }
    selected
}

#[inline]
fn score_ranked_candidate_limit_reached(
    candidates: &[Candidate],
    best_score: f64,
    score_ratio: f64,
    selected: usize,
) -> bool {
    let min_score = best_score * score_ratio;
    candidates
        .get(selected)
        .is_some_and(|candidate| candidate.score as f64 >= min_score)
}

fn unevaluated_candidate_boundary(
    candidates: &[Candidate],
    evaluated: &[bool],
    mappings: &[Mapping],
    index: &RuntimeIndex,
    min_candidate_score: Option<f64>,
    selection: &mut CandidateSelectionWorkspace,
) -> CandidateBoundary {
    selection.profile_groups.clear();
    selection.terminal_groups.clear();
    for mapping in mappings {
        selection
            .profile_groups
            .insert(mapping.candidate_groups.profile_taxon);
        selection
            .terminal_groups
            .insert(mapping.candidate_groups.terminal_taxon);
    }

    let mut boundary = CandidateBoundary::COMPLETE;
    for (candidate, seen) in candidates.iter().zip(evaluated) {
        if *seen || min_candidate_score.is_some_and(|minimum| (candidate.score as f64) < minimum) {
            continue;
        }
        let groups = index.reference_candidate_groups(candidate.ref_id);
        if !selection.profile_groups.contains(&groups.profile_taxon) {
            boundary = boundary.with_profile_incomplete();
        }
        if !selection.terminal_groups.contains(&groups.terminal_taxon) {
            boundary = boundary.with_terminal_incomplete();
        }
        if boundary.bits() == CandidateBoundary::VALID_BITS {
            break;
        }
    }
    boundary
}

/// Keep a bounded set that represents taxonomic diversity before spending any
/// remaining slots on repeated positional placements. Input and output remain
/// in their original best-first order.
fn retain_taxon_diverse_mappings(
    mappings: &mut Vec<Mapping>,
    top_n: usize,
    selection: &mut CandidateSelectionWorkspace,
) -> CandidateBoundary {
    let mut boundary = mappings
        .iter()
        .fold(CandidateBoundary::COMPLETE, |boundary, mapping| {
            boundary.union(mapping.candidate_boundary)
        });
    if mappings.len() <= top_n {
        return boundary;
    }
    if top_n == 0 {
        mappings.clear();
        return boundary
            .with_profile_incomplete()
            .with_terminal_incomplete();
    }

    selection.selected.clear();
    selection.selected.resize(mappings.len(), false);
    selection.profile_groups.clear();
    selection.terminal_groups.clear();
    let mut selected_count = 0usize;

    // First preserve the best mapping for every community-profile endpoint.
    for (index, mapping) in mappings.iter().enumerate() {
        let group = mapping.candidate_groups.profile_taxon;
        if selection.profile_groups.contains(&group) {
            continue;
        }
        if selected_count == top_n {
            boundary = boundary.with_profile_incomplete();
            continue;
        }
        selection.profile_groups.insert(group);
        selection.selected[index] = true;
        selected_count += 1;
    }

    // Then preserve most-specific endpoints (notably strains) in the slots
    // left after broader profile diversity has been represented.
    for (index, mapping) in mappings.iter().enumerate() {
        if selection.selected[index] {
            selection
                .terminal_groups
                .insert(mapping.candidate_groups.terminal_taxon);
        }
    }
    for (index, mapping) in mappings.iter().enumerate() {
        let group = mapping.candidate_groups.terminal_taxon;
        if selection.terminal_groups.contains(&group) {
            continue;
        }
        if selected_count == top_n {
            boundary = boundary.with_terminal_incomplete();
            continue;
        }
        selection.terminal_groups.insert(group);
        selection.selected[index] = true;
        selected_count += 1;
    }

    // Preserve the historical best-position ordering for any capacity that
    // taxonomic representatives did not consume.
    if selected_count < top_n {
        for selected in &mut selection.selected {
            if !*selected {
                *selected = true;
                selected_count += 1;
                if selected_count == top_n {
                    break;
                }
            }
        }
    }

    let mut index = 0usize;
    mappings.retain(|_| {
        let keep = selection.selected[index];
        index += 1;
        keep
    });
    boundary
}

fn should_rescue_long_read(mappings: &[Mapping], params: MapParams) -> bool {
    if !params.long_read_chaining || !params.long_read_rescue {
        return false;
    }
    let rescue_threshold = params.min_identity.max(LONG_READ_RESCUE_IDENTITY);
    mappings
        .iter()
        .map(|mapping| mapping.identity)
        .max_by(f64::total_cmp)
        .is_none_or(|best_identity| best_identity < rescue_threshold)
}

fn collect_candidates<F>(
    read: &[u8],
    index: &RuntimeIndex,
    params: MapParams,
    profile: &mut MapProfile,
    mut accept: F,
) -> (Vec<u8>, Vec<Candidate>)
where
    F: FnMut(SeedOccurrence) -> bool,
{
    if params.long_read_chaining {
        return collect_chained_candidates(read, index, params, profile, accept);
    }

    let rc = revcomp(read);
    let k = index.header().k;
    let mut groups: CandidateMap =
        HashMap::with_capacity_and_hasher(512, BuildHasherDefault::default());
    for (query, strand) in [(read, 1i8), (rc.as_slice(), -1i8)] {
        let query_seeds = seeds_with_strobe_start_distance(
            query,
            index.header().k,
            index.header().w,
            index.header().seed_mode,
            index.header().strobe_start_dist,
        );
        let stride = params.query_seed_stride.max(1);
        let mut used_query_seeds = 0u64;
        for (seed_idx, seed) in query_seeds.into_iter().enumerate() {
            if seed_idx % stride != 0 {
                continue;
            }
            used_query_seeds += 1;
            let occs = index.lookup(seed.hash);
            for occ in occs.iter() {
                profile.total_seed_hits += 1;
                if !accept(occ) {
                    continue;
                }
                let diag = i64::from(occ.pos) - i64::from(seed.pos);
                let key = CandidateKey {
                    ref_id: occ.ref_id,
                    strand,
                    diag_bin: diag_bin(diag, params.diag_bin_size),
                };
                groups
                    .entry(key)
                    .and_modify(|acc| acc.add(seed.pos, diag, k))
                    .or_insert_with(|| CandidateAccumulator::new(seed.pos, diag, k));
            }
        }
        profile.total_query_seeds += used_query_seeds;
    }
    profile.candidate_groups += groups.len() as u64;
    let candidates = groups
        .into_iter()
        .filter_map(|(key, acc)| acc.into_candidate(key, k, params.min_seeds))
        .collect();
    (rc, candidates)
}

fn collect_chained_candidates<F>(
    read: &[u8],
    index: &RuntimeIndex,
    params: MapParams,
    profile: &mut MapProfile,
    mut accept: F,
) -> (Vec<u8>, Vec<Candidate>)
where
    F: FnMut(SeedOccurrence) -> bool,
{
    let rc = revcomp(read);
    let mut groups: ChainMap = HashMap::with_capacity_and_hasher(64, BuildHasherDefault::default());
    for (query, strand) in [(read, 1i8), (rc.as_slice(), -1i8)] {
        let query_seeds = seeds_with_strobe_start_distance(
            query,
            index.header().k,
            index.header().w,
            index.header().seed_mode,
            index.header().strobe_start_dist,
        );
        let stride = params.query_seed_stride.max(1);
        let mut used_query_seeds = 0u64;
        for (seed_idx, seed) in query_seeds.into_iter().enumerate() {
            if seed_idx % stride != 0 {
                continue;
            }
            used_query_seeds += 1;
            let occs = index.lookup(seed.hash);
            for occ in occs.iter() {
                profile.total_seed_hits += 1;
                if !accept(occ) {
                    continue;
                }
                groups
                    .entry(ChainKey {
                        ref_id: occ.ref_id,
                        strand,
                    })
                    .or_default()
                    .push(Anchor {
                        q: seed.pos,
                        r: occ.pos,
                    });
            }
        }
        profile.total_query_seeds += used_query_seeds;
    }
    profile.candidate_groups += groups.len() as u64;

    let k = index.header().k;
    let read_len = read.len();
    let mut candidates = Vec::with_capacity(groups.len());
    for (key, mut anchors) in groups {
        if let Some(candidate) = chain_anchors(key, &mut anchors, k, read_len, params.min_seeds) {
            candidates.push(candidate);
        }
    }
    (rc, candidates)
}

fn prepare_mate_rescue_hits(
    read: &[u8],
    index: &RuntimeIndex,
    params: MapParams,
    profile: &mut MapProfile,
    workspace: &mut MateRescueWorkspace,
) {
    workspace.forward_hits.clear();
    workspace.reverse_hits.clear();
    revcomp_into(read, &mut workspace.reverse_complement);

    collect_rescue_seed_hits(read, index, params, profile, &mut workspace.forward_hits);
    collect_rescue_seed_hits(
        &workspace.reverse_complement,
        index,
        params,
        profile,
        &mut workspace.reverse_hits,
    );
}

fn collect_rescue_seed_hits(
    query: &[u8],
    index: &RuntimeIndex,
    params: MapParams,
    profile: &mut MapProfile,
    output: &mut Vec<RescueSeedHit>,
) {
    let query_seeds = seeds_with_strobe_start_distance(
        query,
        index.header().k,
        index.header().w,
        index.header().seed_mode,
        index.header().strobe_start_dist,
    );
    let stride = params.query_seed_stride.max(1);
    let mut used_query_seeds = 0u64;
    for (seed_index, seed) in query_seeds.into_iter().enumerate() {
        if seed_index % stride != 0 {
            continue;
        }
        used_query_seeds += 1;
        for occurrence in index.lookup(seed.hash).iter() {
            profile.total_seed_hits += 1;
            output.push(RescueSeedHit {
                query_pos: seed.pos,
                occurrence,
            });
        }
    }
    profile.total_query_seeds += used_query_seeds;
}

fn collect_cached_rescue_candidates(
    read_len: usize,
    ref_id: u32,
    window_start: usize,
    window_end: usize,
    index: &RuntimeIndex,
    params: MapParams,
    profile: &mut MapProfile,
    workspace: &mut MateRescueWorkspace,
) -> Vec<Candidate> {
    let k = index.header().k;
    if params.long_read_chaining {
        let mut groups: ChainMap =
            HashMap::with_capacity_and_hasher(64, BuildHasherDefault::default());
        for (hits, strand) in [
            (workspace.forward_hits.as_slice(), 1i8),
            (workspace.reverse_hits.as_slice(), -1i8),
        ] {
            for hit in hits {
                let occurrence = hit.occurrence;
                if occurrence.ref_id != ref_id
                    || !(window_start..window_end).contains(&(occurrence.pos as usize))
                {
                    continue;
                }
                groups
                    .entry(ChainKey { ref_id, strand })
                    .or_default()
                    .push(Anchor {
                        q: hit.query_pos,
                        r: occurrence.pos,
                    });
            }
        }
        profile.candidate_groups += groups.len() as u64;
        return groups
            .into_iter()
            .filter_map(|(key, mut anchors)| {
                chain_anchors(key, &mut anchors, k, read_len, params.min_seeds)
            })
            .collect();
    }

    let mut groups: CandidateMap =
        HashMap::with_capacity_and_hasher(512, BuildHasherDefault::default());
    for (hits, strand) in [
        (workspace.forward_hits.as_slice(), 1i8),
        (workspace.reverse_hits.as_slice(), -1i8),
    ] {
        for hit in hits {
            let occurrence = hit.occurrence;
            if occurrence.ref_id != ref_id
                || !(window_start..window_end).contains(&(occurrence.pos as usize))
            {
                continue;
            }
            let diag = i64::from(occurrence.pos) - i64::from(hit.query_pos);
            let key = CandidateKey {
                ref_id,
                strand,
                diag_bin: diag_bin(diag, params.diag_bin_size),
            };
            groups
                .entry(key)
                .and_modify(|acc| acc.add(hit.query_pos, diag, k))
                .or_insert_with(|| CandidateAccumulator::new(hit.query_pos, diag, k));
        }
    }
    profile.candidate_groups += groups.len() as u64;
    groups
        .into_iter()
        .filter_map(|(key, acc)| acc.into_candidate(key, k, params.min_seeds))
        .collect()
}

fn chain_anchors(
    key: ChainKey,
    anchors: &mut Vec<Anchor>,
    k: usize,
    read_len: usize,
    min_seeds: usize,
) -> Option<Candidate> {
    if anchors.len() < min_seeds {
        return None;
    }
    anchors.sort_unstable_by_key(|anchor| (anchor.q, anchor.r));
    anchors.dedup();

    let n = anchors.len();
    if n < min_seeds {
        return None;
    }
    let max_prev = 128usize;
    let max_gap = read_len.saturating_mul(2).clamp(2_000, 50_000) as u32;
    let mut scores = vec![k as i64; n];
    let mut lengths = vec![1usize; n];
    let mut predecessors = vec![None; n];
    let mut best = 0usize;

    for i in 0..n {
        let anchor_i = anchors[i];
        let lower = i.saturating_sub(max_prev);
        for j in (lower..i).rev() {
            let anchor_j = anchors[j];
            if anchor_j.q >= anchor_i.q || anchor_j.r >= anchor_i.r {
                continue;
            }
            let dq = anchor_i.q - anchor_j.q;
            if dq > max_gap {
                break;
            }
            let dr = anchor_i.r - anchor_j.r;
            if dr > max_gap {
                continue;
            }
            let gap_diff = dq.abs_diff(dr) as i64;
            let gap = dq.max(dr) as i64;
            let penalty = gap_diff / 4 + gap / 5000;
            let step = k as i64 - penalty;
            if step <= 0 {
                continue;
            }
            let score = scores[j] + step;
            if score > scores[i] {
                scores[i] = score;
                lengths[i] = lengths[j] + 1;
                predecessors[i] = Some(j);
            }
        }
        if scores[i] > scores[best] || (scores[i] == scores[best] && lengths[i] > lengths[best]) {
            best = i;
        }
    }

    if lengths[best] < min_seeds {
        return None;
    }

    let mut chain = Vec::with_capacity(lengths[best]);
    let mut cursor = Some(best);
    while let Some(i) = cursor {
        chain.push(anchors[i]);
        cursor = predecessors[i];
    }
    chain.reverse();

    let first = chain.first().copied()?;
    let last = chain.last().copied()?;
    let query_seed_start = first.q as usize;
    let query_seed_end = (last.q as usize + k).min(read_len);
    let ref_seed_start = first.r as usize;
    let ref_seed_end = last.r as usize + k;
    let diag_sum = chain
        .iter()
        .map(|anchor| i128::from(anchor.r) - i128::from(anchor.q))
        .sum::<i128>();
    let diag = (diag_sum / chain.len() as i128) as i64;
    let q_span = query_seed_end.saturating_sub(query_seed_start) as i64;
    let r_span = ref_seed_end.saturating_sub(ref_seed_start) as i64;
    let span_bonus = q_span.min(r_span).min(read_len as i64);
    let drift = (q_span - r_span).abs();
    let score = scores[best] + span_bonus / 8 - drift / 4;

    Some(Candidate {
        ref_id: key.ref_id,
        strand: key.strand,
        diag,
        seed_count: chain.len(),
        score,
        query_seed_start,
        query_seed_end,
        ref_seed_start,
        ref_seed_end,
    })
}

fn chain_interval(
    query_len: usize,
    ref_len: usize,
    candidate: &Candidate,
    full_query_interval: bool,
) -> Option<(usize, usize, usize, usize)> {
    if candidate.query_seed_start >= candidate.query_seed_end
        || candidate.ref_seed_start >= candidate.ref_seed_end
    {
        return None;
    }

    let (q_start, q_end) = if full_query_interval {
        (0, query_len)
    } else {
        (
            candidate.query_seed_start.min(query_len),
            candidate.query_seed_end.min(query_len),
        )
    };
    if q_start >= q_end {
        return None;
    }

    let chain_q = candidate
        .query_seed_end
        .saturating_sub(candidate.query_seed_start)
        .max(1);
    let chain_r = candidate
        .ref_seed_end
        .saturating_sub(candidate.ref_seed_start)
        .max(1);
    let left_q = candidate.query_seed_start.saturating_sub(q_start);
    let right_q = q_end.saturating_sub(candidate.query_seed_end);
    let left_r = scale_query_to_ref(left_q, chain_q, chain_r);
    let right_r = scale_query_to_ref(right_q, chain_q, chain_r);
    let ref_start = candidate.ref_seed_start.saturating_sub(left_r);
    let ref_end = candidate.ref_seed_end.saturating_add(right_r).min(ref_len);

    if ref_start >= ref_end {
        None
    } else {
        Some((q_start, q_end, ref_start, ref_end))
    }
}

fn scale_query_to_ref(
    query_bases: usize,
    chain_query_bases: usize,
    chain_ref_bases: usize,
) -> usize {
    let scaled = (query_bases as u128 * chain_ref_bases as u128) / chain_query_bases as u128;
    scaled.min(usize::MAX as u128) as usize
}

fn diag_bin(diag: i64, bin_size: i64) -> i64 {
    diag.div_euclid(bin_size.max(1))
}

fn refine_interval(
    query: &[u8],
    reference: &ReferenceView<'_>,
    diag: i64,
    seed_q_start: usize,
    seed_q_end: usize,
    extension_max_mismatches: usize,
) -> Option<(usize, usize, usize, usize)> {
    let (overlap_ql, overlap_qr, overlap_rl, overlap_rr) =
        projected_overlap(query.len(), reference.len(), diag)?;
    let mut ql = overlap_ql.max(seed_q_start.min(query.len()));
    let mut qr = overlap_qr.min(seed_q_end.min(query.len()));
    let mut rl;
    let mut rr;
    if ql >= qr {
        ql = overlap_ql;
        qr = overlap_qr;
        rl = overlap_rl;
        rr = overlap_rr;
    } else {
        rl = (diag + ql as i64).max(0) as usize;
        let len = (qr - ql).min(reference.len().saturating_sub(rl));
        qr = ql + len;
        rr = rl + len;
    }
    if ql >= qr || rl >= rr {
        return None;
    }

    let mut mismatch_run = 0usize;
    while ql > 0 && rl > 0 {
        let next_q = ql - 1;
        let next_r = rl - 1;
        if matches_reference(query[next_q], reference.code_at(next_r)) {
            mismatch_run = 0;
        } else {
            mismatch_run += 1;
            if mismatch_run > extension_max_mismatches {
                break;
            }
        }
        ql = next_q;
        rl = next_r;
    }
    trim_left_run(query, reference, &mut ql, &mut rl, qr, rr);

    mismatch_run = 0;
    while qr < query.len() && rr < reference.len() {
        if matches_reference(query[qr], reference.code_at(rr)) {
            mismatch_run = 0;
        } else {
            mismatch_run += 1;
            if mismatch_run > extension_max_mismatches {
                break;
            }
        }
        qr += 1;
        rr += 1;
    }
    trim_right_run(query, reference, ql, rl, &mut qr, &mut rr);

    if ql < qr && rl < rr {
        Some((ql, qr, rl, rr))
    } else {
        None
    }
}

fn trim_left_run(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ql: &mut usize,
    rl: &mut usize,
    qr: usize,
    rr: usize,
) {
    while *ql < qr && *rl < rr && !matches_reference(query[*ql], reference.code_at(*rl)) {
        *ql += 1;
        *rl += 1;
    }
}

fn trim_right_run(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ql: usize,
    rl: usize,
    qr: &mut usize,
    rr: &mut usize,
) {
    while *qr > ql && *rr > rl && !matches_reference(query[*qr - 1], reference.code_at(*rr - 1)) {
        *qr -= 1;
        *rr -= 1;
    }
}

/// Whether a query base matches a reference position.
///
/// Takes the reference side as an already-decoded code so the caller can read
/// it straight out of the packed reference. Ambiguity on either side is never
/// a match.
#[inline]
fn matches_reference(q: u8, r: Option<u8>) -> bool {
    matches!((base_code(q), r), (Some(a), Some(b)) if a == b)
}

fn identity(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ref_start: usize,
    ref_end: usize,
    ignore_n: bool,
    mode: IdentityMode,
    fast_long_read_identity: bool,
    workspace: &mut EditDistanceWorkspace,
) -> Identity {
    match mode {
        IdentityMode::Gapped => {
            if fast_long_read_identity && query.len().max(ref_end - ref_start) > 512 {
                fast_edit_identity(query, reference, ref_start, ref_end, ignore_n, workspace)
            } else {
                gapped_identity(query, reference, ref_start, ref_end, ignore_n)
            }
        }
        IdentityMode::Ungapped => ungapped_identity(query, reference, ref_start, ref_end, ignore_n),
    }
}

fn ungapped_identity(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ref_start: usize,
    ref_end: usize,
    ignore_n: bool,
) -> Identity {
    let mut matches = 0usize;
    let mut compared = 0usize;
    let mut denom = 0usize;
    for (i, &q) in query.iter().take(ref_end - ref_start).enumerate() {
        let qb = base_code(q);
        let rb = reference.code_at(ref_start + i);
        if ignore_n && (qb.is_none() || rb.is_none()) {
            continue;
        }
        compared += 1;
        denom += 1;
        if matches!((qb, rb), (Some(a), Some(b)) if a == b) {
            matches += 1;
        }
    }
    Identity {
        matches,
        compared,
        denom,
        value: if denom == 0 {
            0.0
        } else {
            matches as f64 / denom as f64
        },
    }
}

fn gapped_identity(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ref_start: usize,
    ref_end: usize,
    ignore_n: bool,
) -> Identity {
    let n = query.len();
    let m = ref_end - ref_start;
    if n == m && n <= 512 {
        return ungapped_identity(query, reference, ref_start, ref_end, ignore_n);
    }
    if n == 0 || m == 0 {
        return Identity {
            matches: 0,
            compared: n,
            denom: n.max(m),
            value: 0.0,
        };
    }
    if n.saturating_mul(m) > 250_000 {
        return banded_gapped_identity(query, reference, ref_start, ref_end, ignore_n);
    }
    let mut prev = vec![DpCell::default(); m + 1];
    for j in 1..=m {
        let rb = reference.code_at(ref_start + j - 1);
        prev[j] = ref_gap(prev[j - 1], rb, ignore_n);
    }
    let mut curr = vec![DpCell::default(); m + 1];
    for i in 1..=n {
        let qb = base_code(query[i - 1]);
        curr[0] = query_gap(prev[0], qb, ignore_n);
        for j in 1..=m {
            let rb = reference.code_at(ref_start + j - 1);
            let diag = diag_cell(prev[j - 1], qb, rb, ignore_n);
            let up = query_gap(prev[j], qb, ignore_n);
            let left = ref_gap(curr[j - 1], rb, ignore_n);
            curr[j] = best_cell(diag, up, left);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let best = prev[m];

    Identity {
        matches: best.matches,
        compared: best.compared,
        denom: best.denom,
        value: if best.denom == 0 {
            0.0
        } else {
            best.matches as f64 / best.denom as f64
        },
    }
}

fn banded_gapped_identity(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ref_start: usize,
    ref_end: usize,
    ignore_n: bool,
) -> Identity {
    let n = query.len();
    let m = ref_end - ref_start;
    let band = n.abs_diff(m).max(n.max(m) / 20).max(64);
    let impossible = DpCell {
        score: i32::MIN / 4,
        matches: 0,
        compared: 0,
        denom: 0,
    };
    let mut prev = vec![impossible; m + 1];
    prev[0] = DpCell::default();
    for j in 1..=m.min(band) {
        let rb = reference.code_at(ref_start + j - 1);
        prev[j] = ref_gap(prev[j - 1], rb, ignore_n);
    }

    let mut curr = vec![impossible; m + 1];
    for i in 1..=n {
        let qb = base_code(query[i - 1]);
        let j_start = i.saturating_sub(band);
        let j_end = (i + band).min(m);
        // Clearing all m+1 cells made this banded DP cost O(n*m), defeating the
        // band. Only the band plus one cell either side is ever read: this row
        // reads curr[j_start-1] as its left neighbour, and the next row reads
        // prev[j_end+1], which no row wrote. Cells outside that span are never
        // read before being cleared again, so leaving them stale is safe.
        curr[j_start.saturating_sub(1)..=(j_end + 1).min(m)].fill(impossible);
        if j_start == 0 && prev[0].score > impossible.score {
            curr[0] = query_gap(prev[0], qb, ignore_n);
        }
        for j in j_start.max(1)..=j_end {
            let rb = reference.code_at(ref_start + j - 1);
            let mut best = impossible;
            if prev[j - 1].score > impossible.score {
                best = better_cell(best, diag_cell(prev[j - 1], qb, rb, ignore_n));
            }
            if prev[j].score > impossible.score {
                best = better_cell(best, query_gap(prev[j], qb, ignore_n));
            }
            if curr[j - 1].score > impossible.score {
                best = better_cell(best, ref_gap(curr[j - 1], rb, ignore_n));
            }
            curr[j] = best;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let best = prev[m];
    if best.score <= impossible.score {
        return ungapped_identity(query, reference, ref_start, ref_end, ignore_n);
    }
    Identity {
        matches: best.matches,
        compared: best.compared,
        denom: best.denom,
        value: if best.denom == 0 {
            0.0
        } else {
            best.matches as f64 / best.denom as f64
        },
    }
}

/// Oracle for [`banded_gapped_identity`]: identical except that it clears the
/// whole DP row every iteration, as that function did before the clear was
/// restricted to the band. Used only to assert the two agree exactly.
#[cfg(test)]
fn banded_gapped_identity_full_row_clear(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ref_start: usize,
    ref_end: usize,
    ignore_n: bool,
) -> Identity {
    let n = query.len();
    let m = ref_end - ref_start;
    let band = n.abs_diff(m).max(n.max(m) / 20).max(64);
    let impossible = DpCell {
        score: i32::MIN / 4,
        matches: 0,
        compared: 0,
        denom: 0,
    };
    let mut prev = vec![impossible; m + 1];
    prev[0] = DpCell::default();
    for j in 1..=m.min(band) {
        let rb = reference.code_at(ref_start + j - 1);
        prev[j] = ref_gap(prev[j - 1], rb, ignore_n);
    }

    let mut curr = vec![impossible; m + 1];
    for i in 1..=n {
        curr.fill(impossible);
        let qb = base_code(query[i - 1]);
        let j_start = i.saturating_sub(band);
        let j_end = (i + band).min(m);
        if j_start == 0 && prev[0].score > impossible.score {
            curr[0] = query_gap(prev[0], qb, ignore_n);
        }
        for j in j_start.max(1)..=j_end {
            let rb = reference.code_at(ref_start + j - 1);
            let mut best = impossible;
            if prev[j - 1].score > impossible.score {
                best = better_cell(best, diag_cell(prev[j - 1], qb, rb, ignore_n));
            }
            if prev[j].score > impossible.score {
                best = better_cell(best, query_gap(prev[j], qb, ignore_n));
            }
            if curr[j - 1].score > impossible.score {
                best = better_cell(best, ref_gap(curr[j - 1], rb, ignore_n));
            }
            curr[j] = best;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let best = prev[m];
    if best.score <= impossible.score {
        return ungapped_identity(query, reference, ref_start, ref_end, ignore_n);
    }
    Identity {
        matches: best.matches,
        compared: best.compared,
        denom: best.denom,
        value: if best.denom == 0 {
            0.0
        } else {
            best.matches as f64 / best.denom as f64
        },
    }
}

/// Fast exact Levenshtein identity for ordinary DNA symbols.
///
/// The previous long-read path used a scalar diagonal band. That path was
/// quick relative to a full dynamic-programming matrix, but still evaluated
/// hundreds of scalar cells per query base and could overestimate distance
/// when the optimal path temporarily left its heuristic band. The multiword
/// Myers recurrence below computes the exact global edit distance using 64
/// query positions per machine word, so it cannot report a lower identity
/// than a band-limited search would.
///
/// `--ignore-n` has deliberately different denominator semantics, so it keeps
/// using the existing gapped implementation rather than treating N as a fifth
/// ordinary symbol.
fn fast_edit_identity(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ref_start: usize,
    ref_end: usize,
    ignore_n: bool,
    workspace: &mut EditDistanceWorkspace,
) -> Identity {
    if ignore_n {
        return banded_gapped_identity(query, reference, ref_start, ref_end, true);
    }

    let query_len = query.len();
    let reference_len = ref_end - ref_start;
    if query_len == 0 || reference_len == 0 {
        let denom = query_len.max(reference_len);
        return Identity {
            matches: 0,
            compared: query_len,
            denom,
            value: 0.0,
        };
    }

    let edit_distance = exact_bitvector_edit_distance(
        query,
        (ref_start..ref_end).map(|pos| reference.code_at(pos).unwrap_or(4)),
        workspace,
    );
    edit_identity(query_len, reference_len, edit_distance)
}

fn edit_identity(query_len: usize, reference_len: usize, edit_distance: usize) -> Identity {
    let denom = query_len.max(reference_len);
    let matches = denom
        .saturating_sub(edit_distance)
        .min(query_len)
        .min(reference_len);
    Identity {
        matches,
        compared: query_len,
        denom,
        value: if denom == 0 {
            0.0
        } else {
            matches as f64 / denom as f64
        },
    }
}

/// Exact global Levenshtein distance using the arbitrary-width form of the
/// Myers bit-vector recurrence. `reference_codes` must yield A/C/G/T as 0..=3
/// and ambiguous bases as 4. Ambiguous symbols deliberately have no equality
/// mask, so N-vs-N is a mismatch just as it is in the scalar gapped identity
/// path when `--ignore-n` is not active.
fn exact_bitvector_edit_distance<I>(
    query: &[u8],
    reference_codes: I,
    workspace: &mut EditDistanceWorkspace,
) -> usize
where
    I: IntoIterator<Item = u8>,
{
    debug_assert!(!query.is_empty());
    const WORD_BITS: usize = u64::BITS as usize;

    let word_count = query.len().div_ceil(WORD_BITS);
    for masks in &mut workspace.peq {
        masks.resize(word_count, 0);
        masks.fill(0);
    }
    for (position, &base) in query.iter().enumerate() {
        if let Some(code) = base_code(base) {
            workspace.peq[usize::from(code)][position / WORD_BITS] |=
                1u64 << (position % WORD_BITS);
        }
    }

    let tail_bits = query.len() % WORD_BITS;
    let tail_mask = if tail_bits == 0 {
        u64::MAX
    } else {
        (1u64 << tail_bits) - 1
    };
    let last_word = word_count - 1;
    let score_bit = 1u64 << ((query.len() - 1) % WORD_BITS);

    workspace.positive.resize(word_count, u64::MAX);
    workspace.positive.fill(u64::MAX);
    workspace.positive[last_word] = tail_mask;
    workspace.negative.resize(word_count, 0);
    workspace.negative.fill(0);

    let mut score = query.len();
    for reference_code in reference_codes {
        let eq = &workspace.peq[usize::from(reference_code)];
        let mut add_carry = false;
        let mut positive_shift_carry = 1u64;
        let mut negative_shift_carry = 0u64;

        for word in 0..word_count {
            let positive = workspace.positive[word];
            let negative = workspace.negative[word];
            let equality = eq[word];
            let vertical = equality | negative;

            // The Myers addition is one arbitrary-width integer operation.
            // Carry must therefore propagate between the 64-bit limbs.
            let (sum, carry1) = (equality & positive).overflowing_add(positive);
            let (sum, carry2) = sum.overflowing_add(u64::from(add_carry));
            add_carry = carry1 || carry2;

            let horizontal = (sum ^ positive) | equality;
            let positive_horizontal = negative | !(horizontal | positive);
            let negative_horizontal = positive & horizontal;

            if word == last_word {
                if positive_horizontal & score_bit != 0 {
                    score += 1;
                } else if negative_horizontal & score_bit != 0 {
                    score -= 1;
                }
            }

            let next_positive_carry = positive_horizontal >> (WORD_BITS - 1);
            let next_negative_carry = negative_horizontal >> (WORD_BITS - 1);
            let shifted_positive = (positive_horizontal << 1) | positive_shift_carry;
            let shifted_negative = (negative_horizontal << 1) | negative_shift_carry;
            positive_shift_carry = next_positive_carry;
            negative_shift_carry = next_negative_carry;

            let mask = if word == last_word {
                tail_mask
            } else {
                u64::MAX
            };
            workspace.positive[word] = (shifted_negative | !(vertical | shifted_positive)) & mask;
            workspace.negative[word] = (shifted_positive & vertical) & mask;
        }
    }
    score
}

#[cfg(test)]
fn scalar_banded_edit_identity(
    query: &[u8],
    reference: &ReferenceView<'_>,
    ref_start: usize,
    ref_end: usize,
    ignore_n: bool,
) -> Identity {
    if ignore_n {
        return banded_gapped_identity(query, reference, ref_start, ref_end, ignore_n);
    }
    let n = query.len();
    let m = ref_end - ref_start;
    if n == 0 || m == 0 {
        let denom = n.max(m);
        return Identity {
            matches: 0,
            compared: n,
            denom,
            value: 0.0,
        };
    }

    let ref_codes = (ref_start..ref_end)
        .map(|pos| reference.code_at(pos))
        .collect::<Vec<_>>();
    let band = n.abs_diff(m).saturating_add(96).max(n.max(m) / 100);
    let band = band.min(n.max(m));
    let inf = i32::MAX / 4;

    let mut prev_start = 0usize;
    let mut prev_end = m.min(band);
    let mut prev = (prev_start..=prev_end)
        .map(|j| j as i32)
        .collect::<Vec<_>>();
    let mut curr = Vec::with_capacity(prev.len());

    for i in 1..=n {
        let start = i.saturating_sub(band);
        let end = i.saturating_add(band).min(m);
        let query_code = base_code(query[i - 1]);
        curr.clear();
        curr.resize(end - start + 1, inf);
        for j in start..=end {
            let idx = j - start;
            if j == 0 {
                curr[idx] = i as i32;
                continue;
            }
            let subst_cost =
                i32::from(!matches!((query_code, ref_codes[j - 1]), (Some(a), Some(b)) if a == b));
            let diag = if (prev_start..=prev_end).contains(&(j - 1)) {
                prev[j - 1 - prev_start].saturating_add(subst_cost)
            } else {
                inf
            };
            let deletion = if (prev_start..=prev_end).contains(&j) {
                prev[j - prev_start].saturating_add(1)
            } else {
                inf
            };
            let insertion = if j > start {
                curr[idx - 1].saturating_add(1)
            } else {
                inf
            };
            curr[idx] = diag.min(deletion).min(insertion);
        }
        std::mem::swap(&mut prev, &mut curr);
        prev_start = start;
        prev_end = end;
    }

    let edit_distance = if (prev_start..=prev_end).contains(&m) {
        prev[m - prev_start].max(0) as usize
    } else {
        return banded_gapped_identity(query, reference, ref_start, ref_end, ignore_n);
    };
    edit_identity(n, m, edit_distance)
}

fn diag_cell(mut cell: DpCell, qb: Option<u8>, rb: Option<u8>, ignore_n: bool) -> DpCell {
    if matches!((qb, rb), (Some(a), Some(b)) if a == b) {
        cell.score += 2;
        cell.matches += 1;
        cell.compared += 1;
        cell.denom += 1;
    } else if ignore_n && (qb.is_none() || rb.is_none()) {
        // Keep score and denominator unchanged for ignored ambiguous bases.
    } else {
        cell.score -= 1;
        cell.compared += 1;
        cell.denom += 1;
    }
    cell
}

fn query_gap(mut cell: DpCell, qb: Option<u8>, ignore_n: bool) -> DpCell {
    cell.score -= 1;
    if !(ignore_n && qb.is_none()) {
        cell.compared += 1;
        cell.denom += 1;
    }
    cell
}

fn ref_gap(mut cell: DpCell, rb: Option<u8>, ignore_n: bool) -> DpCell {
    cell.score -= 1;
    if !(ignore_n && rb.is_none()) {
        cell.denom += 1;
    }
    cell
}

fn best_cell(diag: DpCell, up: DpCell, left: DpCell) -> DpCell {
    if diag.score >= up.score && diag.score >= left.score {
        diag
    } else if up.score >= left.score {
        up
    } else {
        left
    }
}

fn better_cell(a: DpCell, b: DpCell) -> DpCell {
    if a.score >= b.score {
        a
    } else {
        b
    }
}

pub fn apply_pair_awareness(
    r1: &mut [Mapping],
    r2: &mut [Mapping],
    min_insert: usize,
    max_insert: usize,
    enabled: bool,
) {
    score_against_mates(r1, r2, min_insert, max_insert, enabled);
    score_against_mates(r2, r1, min_insert, max_insert, enabled);
    if enabled {
        r1.sort_unstable_by(|a, b| b.adjusted_score.cmp(&a.adjusted_score));
        r2.sort_unstable_by(|a, b| b.adjusted_score.cmp(&a.adjusted_score));
        set_ranks(r1);
        set_ranks(r2);
    }
}

pub fn has_proper_pair(
    r1: &[Mapping],
    r2: &[Mapping],
    min_insert: usize,
    max_insert: usize,
) -> bool {
    r1.iter().any(|a| {
        r2.iter()
            .any(|b| pair_compatibility(a, b, min_insert, max_insert).relation == PairRelation::ProperPair)
    })
}

pub fn limit_mappings(mappings: &mut Vec<Mapping>, top_n: usize) {
    let mut selection = CandidateSelectionWorkspace::default();
    limit_mappings_with_selection(mappings, top_n, &mut selection);
}

/// Allocation-reusing form used by the high-throughput paired mapper.
pub fn limit_mappings_with_workspace(
    mappings: &mut Vec<Mapping>,
    top_n: usize,
    workspace: &mut MappingWorkspace,
) {
    limit_mappings_with_selection(mappings, top_n, &mut workspace.candidate_selection);
}

fn limit_mappings_with_selection(
    mappings: &mut Vec<Mapping>,
    top_n: usize,
    selection: &mut CandidateSelectionWorkspace,
) {
    let boundary = retain_taxon_diverse_mappings(mappings, top_n, selection);
    for mapping in mappings.iter_mut() {
        mapping.candidate_boundary = mapping.candidate_boundary.union(boundary);
    }
    set_ranks(mappings);
}

fn score_against_mates(
    reads: &mut [Mapping],
    mates: &[Mapping],
    min_insert: usize,
    max_insert: usize,
    enabled: bool,
) {
    for read in reads {
        if mates.is_empty() {
            read.pair_status = PairStatus::MateUnmapped;
            read.template_length = 0;
            read.pair_score = 0;
            read.adjusted_score = read.score;
            continue;
        }
        let mut best = PairCompatibility {
            score: i64::MIN,
            template_length: 0,
            relation: PairRelation::DiscordantRef,
        };
        for mate in mates {
            let candidate = pair_compatibility(read, mate, min_insert, max_insert);
            if candidate.score > best.score {
                best = candidate;
            }
        }
        read.pair_score = best.score.max(0);
        read.template_length = best.template_length;
        read.pair_status = PairStatus::Paired(best.relation);
        read.adjusted_score = if enabled {
            read.score + read.pair_score
        } else {
            read.score
        };
    }
}

/// Which reference a mapping's mate landed on, as reported beside it.
///
/// Only [`Self::Same`] carries meaning downstream; the other variants are
/// diagnostic. Modelling it as an enum removes a `String` allocated once per
/// evaluated candidate and makes the "same reference" test typed rather than a
/// comparison against `"="`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum MateRefStatus {
    /// No mate evaluation applies to this mapping.
    NotApplicable,
    /// The mate mapped to this same reference.
    Same,
    /// The mate produced no mapping at all.
    Unmapped,
    /// The mate mapped to a different reference, named here.
    Other(Arc<str>),
}

impl MateRefStatus {
    #[inline]
    pub fn is_same(&self) -> bool {
        matches!(self, Self::Same)
    }
}

impl fmt::Display for MateRefStatus {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotApplicable => formatter.write_str("NA"),
            Self::Same => formatter.write_str("="),
            Self::Unmapped => formatter.write_str("*"),
            Self::Other(name) => formatter.write_str(name),
        }
    }
}

/// How two mate placements relate on the reference.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairRelation {
    ProperPair,
    SameRefBadInsert,
    SameRefBadOrientation,
    DiscordantRef,
}

impl PairRelation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ProperPair => "proper_pair",
            Self::SameRefBadInsert => "same_ref_bad_insert",
            Self::SameRefBadOrientation => "same_ref_bad_orientation",
            Self::DiscordantRef => "discordant_ref",
        }
    }
}

/// Pairing state of one emitted mapping. Every value is a fixed label, so an
/// enum replaces what was a `String` allocated once per evaluated candidate.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairStatus {
    /// Not evaluated against a mate.
    Single,
    /// Evaluated, but the mate produced no mapping.
    MateUnmapped,
    /// Ordinary pair-aware relation.
    Paired(PairRelation),
    /// Placement recovered by mate rescue, with the relation it achieved.
    Rescued(PairRelation),
}

impl PairStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Single => "single",
            Self::MateUnmapped => "mate_unmapped",
            Self::Paired(relation) => relation.as_str(),
            Self::Rescued(PairRelation::ProperPair) => "rescued_proper_pair",
            Self::Rescued(PairRelation::SameRefBadInsert) => "rescued_same_ref_bad_insert",
            Self::Rescued(PairRelation::SameRefBadOrientation) => {
                "rescued_same_ref_bad_orientation"
            }
            Self::Rescued(PairRelation::DiscordantRef) => "rescued_discordant_ref",
        }
    }
}

impl fmt::Display for PairStatus {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PairCompatibility {
    pub score: i64,
    pub template_length: i64,
    pub relation: PairRelation,
}

pub fn pair_compatibility(
    a: &Mapping,
    b: &Mapping,
    min_insert: usize,
    max_insert: usize,
) -> PairCompatibility {
    if a.ref_id != b.ref_id {
        return PairCompatibility {
            score: -50,
            template_length: 0,
            relation: PairRelation::DiscordantRef,
        };
    }
    let start = a.ref_start.min(b.ref_start) as i64;
    let end = a.ref_end.max(b.ref_end) as i64;
    let template = end - start;
    let inward_facing = inward_facing_pair(a, b);
    if inward_facing && (min_insert as i64..=max_insert as i64).contains(&template) {
        PairCompatibility {
            score: 80,
            template_length: template,
            relation: PairRelation::ProperPair,
        }
    } else if inward_facing {
        PairCompatibility {
            score: 45,
            template_length: template,
            relation: PairRelation::SameRefBadInsert,
        }
    } else {
        PairCompatibility {
            score: 25,
            template_length: template,
            relation: PairRelation::SameRefBadOrientation,
        }
    }
}

/// Return whether two opposite-strand mappings face inward on the reference.
///
/// Mate labels do not determine which read is leftmost: a fragment sampled
/// from the reverse genomic strand swaps the R1/R2 positions. Comparing
/// mapping midpoints makes this relation symmetric and also handles
/// overlapping inserts. Exact midpoint ties are accepted when the strands
/// are opposite because neither outward geometry can be inferred.
#[inline]
fn inward_facing_pair(a: &Mapping, b: &Mapping) -> bool {
    let a_midpoint_twice = u64::from(a.ref_start) + u64::from(a.ref_end);
    let b_midpoint_twice = u64::from(b.ref_start) + u64::from(b.ref_end);
    match (a.strand, b.strand) {
        (1, -1) => a_midpoint_twice <= b_midpoint_twice,
        (-1, 1) => b_midpoint_twice <= a_midpoint_twice,
        _ => false,
    }
}

pub fn set_mate_status(r1: &mut [Mapping], r2: &mut [Mapping]) {
    // The fallback reference name is shared by every mapping that needs it, so
    // it is interned once per side rather than cloned per mapping.
    let best_r1: Option<Arc<str>> = r1.first().map(|m| Arc::from(m.ref_id.as_str()));
    let best_r2: Option<Arc<str>> = r2.first().map(|m| Arc::from(m.ref_id.as_str()));
    assign_mate_status(r1, r2, best_r2.as_ref());
    assign_mate_status(r2, r1, best_r1.as_ref());
}

fn assign_mate_status(reads: &mut [Mapping], mates: &[Mapping], fallback: Option<&Arc<str>>) {
    for read in reads.iter_mut() {
        read.mate_ref_status = if mates.iter().any(|mate| mate.ref_id == read.ref_id) {
            MateRefStatus::Same
        } else {
            match fallback {
                Some(name) => MateRefStatus::Other(Arc::clone(name)),
                None => MateRefStatus::Unmapped,
            }
        };
    }
}

#[allow(clippy::too_many_arguments)]
pub fn rescue_mate_with_workspace(
    read_id: &str,
    read: &[u8],
    index: &RuntimeIndex,
    mate_mappings: &[Mapping],
    params: MapParams,
    min_insert: usize,
    max_insert: usize,
    profile: &mut MapProfile,
    workspace: &mut MappingWorkspace,
) -> Vec<Mapping> {
    if read.is_empty() || mate_mappings.is_empty() {
        return Vec::new();
    }
    prepare_mate_rescue_hits(read, index, params, profile, &mut workspace.mate_rescue);
    let mut best: Option<Mapping> = None;
    for mate in mate_mappings {
        let Some(ref_id) = index.reference_index_by_id(&mate.ref_id) else {
            continue;
        };
        let reference = index.reference(ref_id);
        if reference.is_empty() {
            continue;
        }
        let window_start = (mate.ref_start as usize).saturating_sub(max_insert);
        let window_end = (mate.ref_end as usize)
            .saturating_add(max_insert)
            .min(reference.len());

        let mut candidates = collect_cached_rescue_candidates(
            read.len(),
            ref_id,
            window_start,
            window_end,
            index,
            params,
            profile,
            &mut workspace.mate_rescue,
        );
        candidates.sort_unstable_by(|a, b| b.score.cmp(&a.score));
        for candidate in candidates.into_iter().take(params.top_n.max(4)) {
            let query = if candidate.strand == 1 {
                read
            } else {
                workspace.mate_rescue.reverse_complement.as_slice()
            };
            let (seed_start, seed_end) = if params.full_read_interval {
                (0, query.len())
            } else {
                (candidate.query_seed_start, candidate.query_seed_end)
            };
            let Some((qs, qe, rs, re)) = refine_interval(
                query,
                &reference,
                candidate.diag,
                seed_start,
                seed_end,
                params.extension_max_mismatches,
            ) else {
                continue;
            };
            if re <= window_start || rs >= window_end {
                continue;
            }
            let ident = identity(
                query,
                &reference,
                rs,
                re,
                params.ignore_n,
                params.identity_mode,
                params.fast_long_read_identity,
                &mut workspace.edit_distance,
            );
            profile.identity_comparisons += 1;
            let mapped_frac = ident.compared as f64 / read.len() as f64;
            let passes = params.no_identity_filter
                || (ident.compared > 0
                    && ident.value >= params.min_identity
                    && mapped_frac >= params.min_mapped_frac);
            if !passes {
                continue;
            }
            let base_score = candidate.score + (ident.value * 1000.0) as i64;
            let mut mapping = Mapping {
                read_id: read_id.to_string(),
                ref_id: index.reference_id(ref_id).to_string(),
                ref_start: rs as u32,
                ref_end: re as u32,
                identity: ident.value,
                read_len: read.len(),
                mapped_bases: ident.compared,
                mate_ref_status: MateRefStatus::NotApplicable,
                strand: candidate.strand,
                seed_count: candidate.seed_count,
                score: base_score,
                second_score: 0,
                rank: 1,
                num_mappings: 1,
                matched_bases: ident.matches,
                query_start: qs,
                query_end: qe,
                pair_score: 0,
                adjusted_score: base_score,
                template_length: 0,
                pair_status: PairStatus::Single,
                mate_number: 0,
                pair_id: 0,
                candidate_groups: index.reference_candidate_groups(ref_id),
                candidate_boundary: CandidateBoundary::COMPLETE,
            };
            let pair = pair_compatibility(&mapping, mate, min_insert, max_insert);
            mapping.pair_score = pair.score.max(0);
            mapping.template_length = pair.template_length;
            mapping.pair_status = PairStatus::Rescued(pair.relation);
            mapping.score = base_score + mapping.pair_score;
            mapping.adjusted_score = mapping.score;
            if best
                .as_ref()
                .map(|prev| {
                    mapping.score > prev.score
                        || (mapping.score == prev.score && mapping.identity > prev.identity)
                })
                .unwrap_or(true)
            {
                best = Some(mapping);
            }
        }

        if best.is_none() && window_end.saturating_sub(window_start) <= 256 {
            best = brute_force_rescue(
                read_id,
                read,
                &workspace.mate_rescue.reverse_complement,
                &reference,
                ref_id,
                mate,
                params,
                min_insert,
                max_insert,
                window_start,
                window_end,
                index,
                profile,
                &mut workspace.edit_distance,
            );
        }
    }
    best.into_iter().collect()
}

#[allow(clippy::too_many_arguments)]
fn brute_force_rescue(
    read_id: &str,
    read: &[u8],
    rc: &[u8],
    reference: &ReferenceView<'_>,
    ref_id: u32,
    mate: &Mapping,
    params: MapParams,
    min_insert: usize,
    max_insert: usize,
    window_start: usize,
    window_end: usize,
    index: &RuntimeIndex,
    profile: &mut MapProfile,
    edit_distance: &mut EditDistanceWorkspace,
) -> Option<Mapping> {
    let mut best = None;
    for (query, strand) in [(read, 1i8), (rc, -1i8)] {
        let min_len = query.len().saturating_sub(8).max(1);
        let max_len = query.len().saturating_add(8).min(reference.len());
        for ref_len in min_len..=max_len {
            if window_end.saturating_sub(window_start) < ref_len {
                continue;
            }
            for pos in window_start..=window_end - ref_len {
                let ident = identity(
                    query,
                    reference,
                    pos,
                    pos + ref_len,
                    params.ignore_n,
                    params.identity_mode,
                    params.fast_long_read_identity,
                    edit_distance,
                );
                profile.identity_comparisons += 1;
                let mapped_frac = ident.compared as f64 / read.len() as f64;
                let passes = params.no_identity_filter
                    || (ident.compared > 0
                        && ident.value >= params.min_identity
                        && mapped_frac >= params.min_mapped_frac);
                if !passes {
                    continue;
                }
                let base_score = (ident.value * 1000.0) as i64 + ident.compared as i64;
                let mut mapping = Mapping {
                    read_id: read_id.to_string(),
                    ref_id: index.reference_id(ref_id).to_string(),
                    ref_start: pos as u32,
                    ref_end: (pos + ref_len) as u32,
                    identity: ident.value,
                    read_len: read.len(),
                    mapped_bases: ident.compared,
                    mate_ref_status: MateRefStatus::NotApplicable,
                    strand,
                    seed_count: 0,
                    score: base_score,
                    second_score: 0,
                    rank: 1,
                    num_mappings: 1,
                    matched_bases: ident.matches,
                    query_start: 0,
                    query_end: read.len(),
                    pair_score: 0,
                    adjusted_score: base_score,
                    template_length: 0,
                    pair_status: PairStatus::Single,
                    mate_number: 0,
                    pair_id: 0,
                    candidate_groups: index.reference_candidate_groups(ref_id),
                    candidate_boundary: CandidateBoundary::COMPLETE,
                };
                let pair = pair_compatibility(&mapping, mate, min_insert, max_insert);
                mapping.pair_score = pair.score.max(0);
                mapping.template_length = pair.template_length;
                mapping.pair_status = PairStatus::Rescued(pair.relation);
                mapping.score = base_score + mapping.pair_score;
                mapping.adjusted_score = mapping.score;
                if best
                    .as_ref()
                    .map(|prev: &Mapping| {
                        mapping.score > prev.score
                            || (mapping.score == prev.score && mapping.identity > prev.identity)
                    })
                    .unwrap_or(true)
                {
                    best = Some(mapping);
                }
            }
        }
    }
    best
}

fn set_ranks(mappings: &mut [Mapping]) {
    let n = mappings.len();
    for (i, mapping) in mappings.iter_mut().enumerate() {
        mapping.rank = i + 1;
        mapping.num_mappings = n;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeds::seeds;
    use crate::{
        index::{
            Index, IndexHeader, Reference, SeedOccurrence, SeedTable, INDEX_MAGIC, INDEX_VERSION,
        },
        seeds::SeedMode,
        sequence::{pack_dna_2bit, BitMask},
    };

    fn tiny_runtime_index() -> RuntimeIndex {
        let seq = b"ACGTACGTACGTACGTACGT".to_vec();
        let reference = Reference {
            id: "r1".to_string(),
            len: seq.len() as u64,
            packed_bases: pack_dna_2bit(&seq),
            ambiguous: BitMask::from_ambiguous(&seq),
        };
        let mut seed_rows = seeds(&seq, 8, 4, SeedMode::Minimizer)
            .into_iter()
            .map(|s| {
                (
                    s.hash,
                    SeedOccurrence {
                        ref_id: 0,
                        pos: s.pos,
                    },
                )
            })
            .collect::<Vec<_>>();
        seed_rows
            .sort_unstable_by_key(|(hash, occurrence)| (*hash, occurrence.ref_id, occurrence.pos));
        let mut seed_table = SeedTable::default();
        let mut i = 0usize;
        while i < seed_rows.len() {
            let hash = seed_rows[i].0;
            let mut j = i + 1;
            while j < seed_rows.len() && seed_rows[j].0 == hash {
                j += 1;
            }
            let occurrences = seed_rows[i..j]
                .iter()
                .map(|(_, occurrence)| *occurrence)
                .collect::<Vec<_>>();
            seed_table
                .push_occurrence_group(hash, &occurrences)
                .expect("valid tiny seed group");
            i = j;
        }
        let index = Index {
            header: IndexHeader {
                magic: INDEX_MAGIC,
                version: INDEX_VERSION,
                metam_version: "test".to_string(),
                k: 8,
                w: 4,
                seed_mode: SeedMode::Minimizer,
                max_occ: 100,
                hash_bits: 64,
                exact_randstrobe_match: false,
                strobe_start_dist: None,
                reference_count: 1,
                total_reference_bases: seq.len() as u64,
                indexed_seed_hashes: seed_table.group_count() as u64,
                indexed_seed_occurrences: seed_table.occurrence_count() as u64,
                filtered_seed_hashes: 0,
                reference_storage_bytes: 0,
                created_unix_seconds: 0,
            },
            references: vec![reference],
            seed_table,
        };
        RuntimeIndex::new(index)
    }

    fn deterministic_read(len: usize) -> Vec<u8> {
        let mut state = 0x6eed_5eed_1234_5678u64;
        let mut out = Vec::with_capacity(len);
        for i in 0..len {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407 + i as u64);
            out.push(match (state >> 62) & 3 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            });
        }
        out
    }

    #[test]
    fn band_limited_row_clear_matches_full_row_clear() {
        // The band-limited clear is an optimization only: it must agree with the
        // previous full-row clear on every input, including the shapes that
        // exercise the band edges (query longer than reference and vice versa,
        // spans below and above the 64-cell band floor, and ambiguous bases).
        let mut state = 0x51ed_5eed_0bad_c0deu64;
        let mut compared = 0usize;
        for case in 0..240usize {
            let query_len = 1 + (case * 7 % 300);
            let ref_len = 1 + (case * 13 % 400);
            let query = deterministic_bases(&mut state, query_len);
            let reference_sequence = deterministic_bases(&mut state, ref_len);
            let reference = owned_reference("band", &reference_sequence);
            let reference = ReferenceView::Owned(&reference);
            for ignore_n in [false, true] {
                let banded =
                    banded_gapped_identity(&query, &reference, 0, reference.len(), ignore_n);
                let oracle = banded_gapped_identity_full_row_clear(
                    &query,
                    &reference,
                    0,
                    reference.len(),
                    ignore_n,
                );
                assert_eq!(
                    (banded.matches, banded.compared, banded.denom),
                    (oracle.matches, oracle.compared, oracle.denom),
                    "case {case} n={query_len} m={ref_len} ignore_n={ignore_n}"
                );
                assert_eq!(banded.value.to_bits(), oracle.value.to_bits());
                compared += 1;
            }
        }
        assert_eq!(compared, 480);
    }

    fn deterministic_bases(state: &mut u64, len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| {
                *state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                match (*state >> 61) & 7 {
                    0 => b'N',
                    1 | 2 => b'A',
                    3 | 4 => b'C',
                    5 => b'G',
                    _ => b'T',
                }
            })
            .collect()
    }

    fn short_test_params() -> MapParams {
        MapParams {
            top_n: 5,
            score_ratio: 0.8,
            min_seeds: 1,
            min_identity: 0.8,
            min_mapped_frac: 0.5,
            min_mapped_len: None,
            no_identity_filter: false,
            ignore_n: false,
            identity_mode: IdentityMode::Gapped,
            diag_bin_size: 1,
            extension_max_mismatches: 4,
            full_read_interval: false,
            fast_long_read_identity: false,
            long_read_chaining: false,
            long_read_extension: LongReadExtension::None,
            query_seed_stride: 1,
            long_read_rescue: false,
            long_read_rescue_limit: None,
        }
    }

    fn candidate_signature(
        candidate: &Candidate,
    ) -> (u32, i8, i64, usize, i64, usize, usize, usize, usize) {
        (
            candidate.ref_id,
            candidate.strand,
            candidate.diag,
            candidate.seed_count,
            candidate.score,
            candidate.query_seed_start,
            candidate.query_seed_end,
            candidate.ref_seed_start,
            candidate.ref_seed_end,
        )
    }

    #[test]
    fn cached_mate_rescue_hits_preserve_windowed_candidates() {
        let index = tiny_runtime_index();
        let read = b"ACGTACGTACGT";

        for chained in [false, true] {
            let mut params = short_test_params();
            params.long_read_chaining = chained;
            let mut workspace = MateRescueWorkspace::default();
            let mut cached_profile = MapProfile::default();
            prepare_mate_rescue_hits(read, &index, params, &mut cached_profile, &mut workspace);
            let cached_query_seeds = cached_profile.total_query_seeds;
            let cached_seed_hits = cached_profile.total_seed_hits;

            for (window_start, window_end) in [(0usize, 20usize), (4, 16)] {
                let mut expected_profile = MapProfile::default();
                let (_, expected) =
                    collect_candidates(read, &index, params, &mut expected_profile, |occurrence| {
                        occurrence.ref_id == 0
                            && (window_start..window_end).contains(&(occurrence.pos as usize))
                    });
                let actual = collect_cached_rescue_candidates(
                    read.len(),
                    0,
                    window_start,
                    window_end,
                    &index,
                    params,
                    &mut cached_profile,
                    &mut workspace,
                );
                let expected = expected.iter().map(candidate_signature).collect::<Vec<_>>();
                let actual = actual.iter().map(candidate_signature).collect::<Vec<_>>();
                assert_eq!(actual, expected, "chained={chained}");
            }

            assert_eq!(cached_profile.total_query_seeds, cached_query_seeds);
            assert_eq!(cached_profile.total_seed_hits, cached_seed_hits);
        }
    }

    fn owned_reference(id: &str, sequence: &[u8]) -> Reference {
        Reference {
            id: id.to_string(),
            len: sequence.len() as u64,
            packed_bases: pack_dna_2bit(sequence),
            ambiguous: BitMask::from_ambiguous(sequence),
        }
    }

    fn mutate_every(seq: &[u8], stride: usize) -> Vec<u8> {
        seq.iter()
            .enumerate()
            .map(|(i, base)| {
                if i % stride == 0 {
                    match base {
                        b'A' => b'C',
                        b'C' => b'G',
                        b'G' => b'T',
                        b'T' => b'A',
                        _ => b'N',
                    }
                } else {
                    *base
                }
            })
            .collect()
    }

    fn full_edit_distance(query: &[u8], reference: &[u8]) -> usize {
        let mut previous = (0..=reference.len()).collect::<Vec<_>>();
        let mut current = vec![0usize; reference.len() + 1];
        for (i, &query_base) in query.iter().enumerate() {
            current[0] = i + 1;
            let query_code = base_code(query_base);
            for (j, &reference_base) in reference.iter().enumerate() {
                let reference_code = base_code(reference_base);
                let substitution = usize::from(
                    !matches!((query_code, reference_code), (Some(a), Some(b)) if a == b),
                );
                current[j + 1] = previous[j + 1]
                    .saturating_add(1)
                    .min(current[j].saturating_add(1))
                    .min(previous[j].saturating_add(substitution));
            }
            std::mem::swap(&mut previous, &mut current);
        }
        previous[reference.len()]
    }

    #[test]
    fn exact_bitvector_distance_matches_full_dp_across_word_boundaries() {
        let lengths = [1usize, 2, 7, 31, 63, 64, 65, 95, 127, 128, 129, 191];
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut workspace = EditDistanceWorkspace::default();

        for query_len in lengths {
            for reference_len in lengths {
                for _ in 0..3 {
                    let mut sequence = |len: usize| {
                        (0..len)
                            .map(|_| {
                                state = state
                                    .wrapping_mul(6364136223846793005)
                                    .wrapping_add(1442695040888963407);
                                match state % 11 {
                                    0 => b'N',
                                    1..=3 => b'A',
                                    4..=5 => b'C',
                                    6..=7 => b'G',
                                    _ => b'T',
                                }
                            })
                            .collect::<Vec<_>>()
                    };
                    let query = sequence(query_len);
                    let reference = sequence(reference_len);
                    let expected = full_edit_distance(&query, &reference);
                    let observed = exact_bitvector_edit_distance(
                        &query,
                        reference.iter().map(|base| base_code(*base).unwrap_or(4)),
                        &mut workspace,
                    );
                    assert_eq!(
                        observed, expected,
                        "query_len={query_len}, reference_len={reference_len}"
                    );
                }
            }
        }
    }

    #[test]
    fn exact_bitvector_workspace_is_safe_to_reuse_for_smaller_queries() {
        let mut workspace = EditDistanceWorkspace::default();
        let large_query = deterministic_read(513);
        let large_reference = mutate_every(&large_query, 17);
        let large = exact_bitvector_edit_distance(
            &large_query,
            large_reference
                .iter()
                .map(|base| base_code(*base).unwrap_or(4)),
            &mut workspace,
        );
        assert_eq!(large, full_edit_distance(&large_query, &large_reference));

        let query = b"ACGTNAC";
        let reference = b"AGGTNTC";
        let observed = exact_bitvector_edit_distance(
            query,
            reference.iter().map(|base| base_code(*base).unwrap_or(4)),
            &mut workspace,
        );
        assert_eq!(observed, full_edit_distance(query, reference));
    }

    #[test]
    fn exact_bitvector_distance_handles_long_indel_drift_and_ambiguous_bases() {
        let mut query = deterministic_read(1_025);
        query[64] = b'N';
        query[512] = b'u';

        // Exercise both insertion and deletion drift beyond several 64-bit
        // boundaries. Ambiguous symbols use code 4 on the reference side but
        // have no query equality mask, so every ambiguous comparison costs
        // one unless `--ignore-n` selects the scalar path.
        let mut reference = Vec::with_capacity(query.len());
        reference.extend_from_slice(&query[..197]);
        reference.extend_from_slice(b"ACGTNNN");
        reference.extend_from_slice(&query[197..733]);
        reference.extend_from_slice(&query[750..]);
        for position in (41..reference.len()).step_by(97) {
            reference[position] = match reference[position].to_ascii_uppercase() {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }

        let expected = full_edit_distance(&query, &reference);
        let mut workspace = EditDistanceWorkspace::default();
        let observed = exact_bitvector_edit_distance(
            &query,
            reference.iter().map(|base| base_code(*base).unwrap_or(4)),
            &mut workspace,
        );
        assert_eq!(observed, expected);
    }

    #[test]
    fn exact_long_read_identity_never_falls_below_legacy_banded_identity() {
        let mut state = 0xd1b5_4a32_d192_ed03u64;
        let mut workspace = EditDistanceWorkspace::default();

        // Deterministic pseudo-random substitutions and indels exercise long
        // intervals, ambiguous symbols, unequal lengths, and workspace reuse.
        // The legacy dynamic program minimizes over a subset of the paths
        // considered by exact global Levenshtein distance, so its identity may
        // be equal or lower, but never higher.
        for case in 0..24usize {
            let query_len = 513 + (state as usize % 1_536);
            let query = deterministic_bases(&mut state, query_len);
            let mut reference_sequence = query.clone();

            let substitutions = 1 + query_len / (31 + case % 23);
            for _ in 0..substitutions {
                state = state
                    .wrapping_mul(2862933555777941757)
                    .wrapping_add(3037000493);
                let position = state as usize % reference_sequence.len();
                reference_sequence[position] = match reference_sequence[position] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    b'T' => b'A',
                    _ => b'T',
                };
            }

            let inserted_len = 1 + (case * 17 % 79);
            let inserted = deterministic_bases(&mut state, inserted_len);
            let insert_at = 1 + state as usize % (reference_sequence.len() - 1);
            reference_sequence.splice(insert_at..insert_at, inserted);

            if case % 2 == 0 {
                let deleted_len = 1 + (case * 11 % 61);
                if reference_sequence.len() > deleted_len + 2 {
                    let delete_at = 1
                        + (state.rotate_left(17) as usize
                            % (reference_sequence.len() - deleted_len - 1));
                    reference_sequence.drain(delete_at..delete_at + deleted_len);
                }
            }

            let reference = owned_reference("property", &reference_sequence);
            let reference = ReferenceView::Owned(&reference);
            let legacy = scalar_banded_edit_identity(&query, &reference, 0, reference.len(), false);
            let exact = fast_edit_identity(
                &query,
                &reference,
                0,
                reference.len(),
                false,
                &mut workspace,
            );
            assert!(
                exact.matches >= legacy.matches,
                "case={case} query_len={} reference_len={} exact={} legacy={}",
                query.len(),
                reference.len(),
                exact.value,
                legacy.value
            );
        }
    }

    #[test]
    fn exact_identity_can_cross_the_long_read_rescue_boundary() {
        // Keep both sequences the same length, but insert 128 bases near the
        // left and remove 128 much later. The optimal path temporarily drifts
        // beyond the legacy 96-base diagonal band and then returns. Exact
        // identity recovers the homologous middle; the banded approximation
        // cannot, which is precisely the case that motivated the Myers kernel.
        let query = deterministic_read(4_096);
        let mut state = 0xa076_1d64_78bd_642fu64;
        let inserted = deterministic_bases(&mut state, 128);
        let mut reference_sequence = Vec::with_capacity(query.len());
        reference_sequence.extend_from_slice(&query[..512]);
        reference_sequence.extend_from_slice(&inserted);
        reference_sequence.extend_from_slice(&query[512..3_072]);
        reference_sequence.extend_from_slice(&query[3_200..]);
        assert_eq!(reference_sequence.len(), query.len());

        let reference = owned_reference("drift", &reference_sequence);
        let reference = ReferenceView::Owned(&reference);
        let legacy = scalar_banded_edit_identity(&query, &reference, 0, reference.len(), false);
        let mut workspace = EditDistanceWorkspace::default();
        let exact = fast_edit_identity(
            &query,
            &reference,
            0,
            reference.len(),
            false,
            &mut workspace,
        );

        assert!(exact.value >= LONG_READ_RESCUE_IDENTITY, "{exact:?}");
        assert!(legacy.value < LONG_READ_RESCUE_IDENTITY, "{legacy:?}");
        assert!(exact.value > legacy.value);

        let mut params = default_long_read_test_params();
        params.min_identity = 0.80;
        let legacy_mapping = test_mapping("read", "drift", 0, 100, 100, legacy.value, 100);
        let exact_mapping = test_mapping("read", "drift", 0, 100, 100, exact.value, 100);
        assert!(should_rescue_long_read(&[legacy_mapping], params));
        assert!(!should_rescue_long_read(&[exact_mapping], params));
    }

    fn biased_long_read_runtime_index(read: &[u8]) -> RuntimeIndex {
        let wrong = mutate_every(read, 7);
        let correct = read.to_vec();
        let references = vec![
            Reference {
                id: "wrong_species".to_string(),
                len: wrong.len() as u64,
                packed_bases: pack_dna_2bit(&wrong),
                ambiguous: BitMask::from_ambiguous(&wrong),
            },
            Reference {
                id: "correct_species".to_string(),
                len: correct.len() as u64,
                packed_bases: pack_dna_2bit(&correct),
                ambiguous: BitMask::from_ambiguous(&correct),
            },
        ];

        let query_seeds = seeds(read, 8, 4, SeedMode::Minimizer);
        assert!(query_seeds.len() > 20);
        let mut seed_rows = Vec::new();
        for seed in query_seeds.iter().step_by(2).take(80) {
            seed_rows.push((
                seed.hash,
                SeedOccurrence {
                    ref_id: 0,
                    pos: seed.pos,
                },
            ));
        }
        let mut correct_indices = vec![0, query_seeds.len() / 2, query_seeds.len() - 1];
        correct_indices.sort_unstable();
        correct_indices.dedup();
        for i in correct_indices {
            let seed = query_seeds[i];
            seed_rows.push((
                seed.hash,
                SeedOccurrence {
                    ref_id: 1,
                    pos: seed.pos,
                },
            ));
        }
        seed_rows
            .sort_unstable_by_key(|(hash, occurrence)| (*hash, occurrence.ref_id, occurrence.pos));

        let mut seed_table = SeedTable::default();
        let mut i = 0usize;
        while i < seed_rows.len() {
            let hash = seed_rows[i].0;
            let mut j = i + 1;
            while j < seed_rows.len() && seed_rows[j].0 == hash {
                j += 1;
            }
            let occurrences = seed_rows[i..j]
                .iter()
                .map(|(_, occurrence)| *occurrence)
                .collect::<Vec<_>>();
            seed_table
                .push_occurrence_group(hash, &occurrences)
                .expect("valid biased seed group");
            i = j;
        }
        seed_table.rebuild_lookup_index();

        let index = Index {
            header: IndexHeader {
                magic: INDEX_MAGIC,
                version: INDEX_VERSION,
                metam_version: "test".to_string(),
                k: 8,
                w: 4,
                seed_mode: SeedMode::Minimizer,
                max_occ: 100,
                hash_bits: 64,
                exact_randstrobe_match: false,
                strobe_start_dist: None,
                reference_count: references.len(),
                total_reference_bases: references.iter().map(|r| r.len).sum(),
                indexed_seed_hashes: seed_table.group_count() as u64,
                indexed_seed_occurrences: seed_table.occurrence_count() as u64,
                filtered_seed_hashes: 0,
                reference_storage_bytes: 0,
                created_unix_seconds: 0,
            },
            references,
            seed_table,
        };
        RuntimeIndex::new(index)
    }

    #[test]
    fn maps_forward_read() {
        let idx = tiny_runtime_index();
        let mut profile = MapProfile::default();
        let out = map_read(
            "r",
            b"ACGTACGTACGT",
            &idx,
            MapParams {
                top_n: 5,
                score_ratio: 0.8,
                min_seeds: 1,
                min_identity: 0.8,
                min_mapped_frac: 0.5,
                min_mapped_len: None,
                no_identity_filter: false,
                ignore_n: false,
                identity_mode: IdentityMode::Gapped,
                diag_bin_size: 1,
                extension_max_mismatches: 4,
                full_read_interval: false,
                fast_long_read_identity: false,
                long_read_chaining: false,
                long_read_extension: LongReadExtension::None,
                query_seed_stride: 1,
                long_read_rescue: false,
                long_read_rescue_limit: None,
            },
            &mut profile,
        );
        assert!(!out.is_empty());
        assert_eq!(out[0].ref_id, "r1");
        assert!(out[0].identity >= 0.99);
    }

    #[test]
    fn long_read_rescue_reranks_low_seed_true_hit_by_identity() {
        let read = deterministic_read(1200);
        let idx = biased_long_read_runtime_index(&read);
        let mut profile = MapProfile::default();
        let out = map_read(
            "ont",
            &read,
            &idx,
            MapParams {
                top_n: 1,
                score_ratio: 0.8,
                min_seeds: 2,
                min_identity: 0.8,
                min_mapped_frac: 0.3,
                min_mapped_len: None,
                no_identity_filter: false,
                ignore_n: false,
                identity_mode: IdentityMode::Gapped,
                diag_bin_size: 128,
                extension_max_mismatches: 40,
                full_read_interval: true,
                fast_long_read_identity: true,
                long_read_chaining: true,
                long_read_extension: LongReadExtension::None,
                query_seed_stride: 1,
                long_read_rescue: true,
                long_read_rescue_limit: None,
            },
            &mut profile,
        );

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].ref_id, "correct_species");
        assert!(out[0].identity > 0.99, "{out:?}");
        assert!(profile.identity_comparisons > 1);
    }

    #[test]
    fn long_read_extend_option_expands_chain_interval() {
        let read = deterministic_read(1600);
        let reference = Reference {
            id: "r1".to_string(),
            len: read.len() as u64,
            packed_bases: pack_dna_2bit(&read),
            ambiguous: BitMask::from_ambiguous(&read),
        };
        let mut seed_rows = seeds(&read, 8, 4, SeedMode::Minimizer)
            .into_iter()
            .filter(|seed| (450..=1050).contains(&(seed.pos as usize)))
            .map(|seed| {
                (
                    seed.hash,
                    SeedOccurrence {
                        ref_id: 0,
                        pos: seed.pos,
                    },
                )
            })
            .collect::<Vec<_>>();
        assert!(seed_rows.len() > 10);
        seed_rows
            .sort_unstable_by_key(|(hash, occurrence)| (*hash, occurrence.ref_id, occurrence.pos));
        let mut seed_table = SeedTable::default();
        let mut i = 0usize;
        while i < seed_rows.len() {
            let hash = seed_rows[i].0;
            let mut j = i + 1;
            while j < seed_rows.len() && seed_rows[j].0 == hash {
                j += 1;
            }
            let occurrences = seed_rows[i..j]
                .iter()
                .map(|(_, occurrence)| *occurrence)
                .collect::<Vec<_>>();
            seed_table
                .push_occurrence_group(hash, &occurrences)
                .expect("valid extension seed group");
            i = j;
        }
        seed_table.rebuild_lookup_index();
        let idx = RuntimeIndex::new(Index {
            header: IndexHeader {
                magic: INDEX_MAGIC,
                version: INDEX_VERSION,
                metam_version: "test".to_string(),
                k: 8,
                w: 4,
                seed_mode: SeedMode::Minimizer,
                max_occ: 100,
                hash_bits: 64,
                exact_randstrobe_match: false,
                strobe_start_dist: None,
                reference_count: 1,
                total_reference_bases: read.len() as u64,
                indexed_seed_hashes: seed_table.group_count() as u64,
                indexed_seed_occurrences: seed_table.occurrence_count() as u64,
                filtered_seed_hashes: 0,
                reference_storage_bytes: 0,
                created_unix_seconds: 0,
            },
            references: vec![reference],
            seed_table,
        });
        let params = MapParams {
            top_n: 1,
            score_ratio: 0.8,
            min_seeds: 2,
            min_identity: 0.9,
            min_mapped_frac: 0.2,
            min_mapped_len: None,
            no_identity_filter: false,
            ignore_n: false,
            identity_mode: IdentityMode::Gapped,
            diag_bin_size: 128,
            extension_max_mismatches: 40,
            full_read_interval: true,
            fast_long_read_identity: true,
            long_read_chaining: true,
            long_read_extension: LongReadExtension::None,
            query_seed_stride: 1,
            long_read_rescue: true,
            long_read_rescue_limit: None,
        };

        let mut default_profile = MapProfile::default();
        let default = map_read("ont", &read, &idx, params, &mut default_profile);
        assert_eq!(default.len(), 1);
        assert!(default[0].query_start > 0, "{default:?}");
        assert!(default[0].query_end < read.len(), "{default:?}");

        let mut extended_profile = MapProfile::default();
        let extended = map_read(
            "ont",
            &read,
            &idx,
            MapParams {
                long_read_extension: LongReadExtension::Candidate,
                ..params
            },
            &mut extended_profile,
        );
        assert_eq!(extended.len(), 1);
        assert_eq!(extended[0].query_start, 0);
        assert_eq!(extended[0].query_end, read.len());
        assert!(extended[0].mapped_bases > default[0].mapped_bases);
        assert!(extended_profile.identity_comparisons > default_profile.identity_comparisons);

        let mut fast_extended_profile = MapProfile::default();
        let fast_extended = map_read(
            "ont",
            &read,
            &idx,
            MapParams {
                long_read_extension: LongReadExtension::Retained,
                ..params
            },
            &mut fast_extended_profile,
        );
        assert_eq!(fast_extended.len(), 1);
        assert_eq!(fast_extended[0].query_start, 0);
        assert_eq!(fast_extended[0].query_end, read.len());
        assert!(fast_extended[0].mapped_bases > default[0].mapped_bases);
        assert!(fast_extended_profile.identity_comparisons > default_profile.identity_comparisons);
    }

    #[test]
    fn retained_dedup_preserves_longest_duplicate_mapping() {
        let mut mappings = vec![
            test_mapping("r", "ref", 10, 50, 40, 0.99, 100),
            test_mapping("r", "ref", 10, 50, 80, 0.98, 120),
            test_mapping("r", "other", 10, 50, 70, 0.99, 110),
        ];

        dedup_mappings_preserve_order(&mut mappings);

        assert_eq!(mappings.len(), 2);
        assert_eq!(mappings[0].ref_id, "ref");
        assert_eq!(mappings[0].mapped_bases, 80);
        assert_eq!(mappings[1].ref_id, "other");
        assert_eq!(mappings[1].mapped_bases, 70);
    }

    #[test]
    fn fast_extend_rejects_non_homologous_added_flank() {
        // A flank well below the x-drop break-even (~0.6) is non-homologous and
        // must be rejected even though the overall mapping still passes the
        // read-level identity filter.
        let params = default_long_read_test_params();
        let base = Identity {
            matches: 2850,
            compared: 3000,
            denom: 3000,
            value: 0.95,
        };
        let extended = Identity {
            matches: 2890, // added flank: 40 matches / 100 denom = 0.40
            compared: 3100,
            denom: 3100,
            value: 2890.0 / 3100.0,
        };

        assert!(
            extended.value >= params.min_identity,
            "overall still passes"
        );
        assert!(!should_use_fast_extended_interval(
            base, extended, true, 1.0, params
        ));
    }

    #[test]
    fn fast_extend_accepts_noisy_homologous_flank() {
        // Real ONT read ends are noisy (~0.75-0.85) but homologous. They sit
        // below the core identity yet above the x-drop floor, so fast-extend
        // must accept them (this is the incomplete-region regression).
        let params = default_long_read_test_params();
        let base = Identity {
            matches: 2850,
            compared: 3000,
            denom: 3000,
            value: 0.95,
        };
        let extended = Identity {
            matches: 2928, // added flank: 78 matches / 100 denom = 0.78
            compared: 3100,
            denom: 3100,
            value: 2928.0 / 3100.0,
        };

        assert!(should_use_fast_extended_interval(
            base, extended, true, 1.0, params
        ));
    }

    #[test]
    fn fast_extend_accepts_high_identity_added_flank() {
        let params = default_long_read_test_params();
        let base = Identity {
            matches: 2850,
            compared: 3000,
            denom: 3000,
            value: 0.95,
        };
        let extended = Identity {
            matches: 2940,
            compared: 3100,
            denom: 3100,
            value: 2940.0 / 3100.0,
        };

        assert!(should_use_fast_extended_interval(
            base, extended, true, 1.0, params
        ));
    }

    #[test]
    fn length_rescue_admits_long_low_fraction_long_read_hits() {
        let mut params = default_long_read_test_params();
        params.min_identity = 0.80;
        params.min_mapped_len = Some(300);
        // 320 bp mapped at 0.85 identity clears the length floor -> rescued even
        // when the mapped fraction would be low.
        let long_hit = Identity {
            matches: 300,
            compared: 320,
            denom: 320,
            value: 0.85,
        };
        assert!(length_rescues(&long_hit, params));
        // Below the 300 bp floor -> not rescued.
        let short_hit = Identity {
            matches: 250,
            compared: 280,
            denom: 280,
            value: 0.85,
        };
        assert!(!length_rescues(&short_hit, params));
        // Long enough but identity under the cutoff -> not rescued.
        let low_identity = Identity {
            matches: 240,
            compared: 320,
            denom: 340,
            value: 0.75,
        };
        assert!(!length_rescues(&low_identity, params));
        // Illumina (no length floor) never rescues.
        params.min_mapped_len = None;
        assert!(!length_rescues(&long_hit, params));
    }

    #[test]
    fn proper_pair_requires_inward_facing_geometry() {
        let mut left = test_mapping("r/1", "ref", 100, 200, 100, 1.0, 100);
        let mut right = test_mapping("r/2", "ref", 300, 400, 100, 1.0, 100);
        left.strand = 1;
        right.strand = -1;
        assert_eq!(
            pair_compatibility(&left, &right, 0, 1_000).relation,
            PairRelation::ProperPair
        );
        assert_eq!(
            pair_compatibility(&right, &left, 0, 1_000).relation,
            PairRelation::ProperPair
        );

        // R1/R2 labels can swap genomic sides for a reverse-oriented
        // fragment; geometry, not mate label, determines proper pairing.
        left.strand = -1;
        right.strand = 1;
        assert_eq!(
            pair_compatibility(&left, &right, 0, 1_000).relation,
            PairRelation::SameRefBadOrientation
        );
        assert_eq!(
            pair_compatibility(&right, &left, 0, 1_000).relation,
            PairRelation::SameRefBadOrientation
        );

        std::mem::swap(&mut left.ref_start, &mut right.ref_start);
        std::mem::swap(&mut left.ref_end, &mut right.ref_end);
        assert_eq!(
            pair_compatibility(&left, &right, 0, 1_000).relation,
            PairRelation::ProperPair
        );
    }

    #[test]
    fn inward_pair_outside_insert_bounds_is_not_proper() {
        let mut left = test_mapping("r/1", "ref", 100, 200, 100, 1.0, 100);
        let mut right = test_mapping("r/2", "ref", 300, 400, 100, 1.0, 100);
        left.strand = 1;
        right.strand = -1;
        let compatibility = pair_compatibility(&left, &right, 350, 500);
        assert_eq!(compatibility.template_length, 300);
        assert_eq!(compatibility.relation, PairRelation::SameRefBadInsert);
    }

    #[test]
    fn fast_identity_counts_ambiguous_pairs_as_mismatches_unless_ignored() {
        let mut query = vec![b'A'; 513];
        query[256] = b'N';
        let reference = owned_reference("ambiguous", &query);
        let reference = ReferenceView::Owned(&reference);
        let mut workspace = EditDistanceWorkspace::default();
        let counted = fast_edit_identity(&query, &reference, 0, query.len(), false, &mut workspace);
        assert_eq!(counted.matches, 512);
        assert_eq!(counted.denom, 513);
        assert_eq!(counted.value.to_bits(), (512.0f64 / 513.0).to_bits());

        let ignored = fast_edit_identity(&query, &reference, 0, query.len(), true, &mut workspace);
        assert_eq!(ignored.matches, 512);
        assert_eq!(ignored.denom, 512);
        assert_eq!(ignored.value.to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn output_limit_preserves_taxon_diversity_before_repeated_positions() {
        let mut mappings = vec![
            with_groups(
                test_mapping("read", "species-a-best", 0, 100, 100, 1.0, 100),
                10,
                11,
            ),
            with_groups(
                test_mapping("read", "species-a-repeat", 10, 110, 100, 0.99, 99),
                10,
                11,
            ),
            with_groups(
                test_mapping("read", "species-b", 0, 100, 100, 0.98, 98),
                20,
                21,
            ),
        ];
        limit_mappings(&mut mappings, 2);
        assert_eq!(
            mappings
                .iter()
                .map(|mapping| mapping.ref_id.as_str())
                .collect::<Vec<_>>(),
            vec!["species-a-best", "species-b"]
        );
        assert!(mappings
            .iter()
            .all(|mapping| mapping.candidate_boundary.is_complete()));
    }

    #[test]
    fn output_limit_keeps_original_best_positions_when_one_taxon_is_present() {
        let mut mappings = vec![
            with_groups(test_mapping("read", "best", 0, 100, 100, 1.0, 100), 10, 11),
            with_groups(
                test_mapping("read", "second", 10, 110, 100, 0.99, 99),
                10,
                11,
            ),
            with_groups(
                test_mapping("read", "third", 20, 120, 100, 0.98, 98),
                10,
                11,
            ),
        ];
        limit_mappings(&mut mappings, 2);
        assert_eq!(
            mappings
                .iter()
                .map(|mapping| mapping.ref_id.as_str())
                .collect::<Vec<_>>(),
            vec!["best", "second"]
        );
        assert!(mappings
            .iter()
            .all(|mapping| mapping.candidate_boundary.is_complete()));
    }

    #[test]
    fn output_limit_distinguishes_terminal_from_profile_incompleteness() {
        let mut mappings = vec![
            with_groups(
                test_mapping("read", "species-a-strain-1", 0, 100, 100, 1.0, 100),
                10,
                11,
            ),
            with_groups(
                test_mapping("read", "species-a-strain-2", 0, 100, 100, 0.99, 99),
                10,
                12,
            ),
            with_groups(
                test_mapping("read", "species-b", 0, 100, 100, 0.98, 98),
                20,
                21,
            ),
        ];
        limit_mappings(&mut mappings, 2);
        assert_eq!(
            mappings
                .iter()
                .map(|mapping| mapping.ref_id.as_str())
                .collect::<Vec<_>>(),
            vec!["species-a-strain-1", "species-b"]
        );
        assert!(mappings.iter().all(|mapping| {
            !mapping.candidate_boundary.profile_incomplete()
                && mapping.candidate_boundary.terminal_incomplete()
        }));
    }

    #[test]
    fn output_limit_marks_unrepresented_profile_taxa() {
        let mut mappings = vec![
            with_groups(test_mapping("read", "a", 0, 100, 100, 1.0, 100), 10, 11),
            with_groups(test_mapping("read", "b", 0, 100, 100, 0.99, 99), 20, 21),
            with_groups(test_mapping("read", "c", 0, 100, 100, 0.98, 98), 30, 31),
        ];
        limit_mappings(&mut mappings, 2);
        assert_eq!(mappings.len(), 2);
        assert!(mappings.iter().all(|mapping| {
            mapping.candidate_boundary.profile_incomplete()
                && mapping.candidate_boundary.terminal_incomplete()
        }));
    }

    #[test]
    fn output_limit_ignores_repeated_placements_on_the_same_reference() {
        let mut mappings = vec![
            test_mapping("read", "ref", 0, 100, 100, 1.0, 100),
            test_mapping("read", "ref", 10, 110, 100, 0.99, 99),
        ];
        limit_mappings(&mut mappings, 1);
        assert_eq!(mappings.len(), 1);
        assert!(mappings[0].candidate_boundary.is_complete());
    }

    #[test]
    fn prior_candidate_boundary_survives_later_nontruncating_limit() {
        let mut mappings = vec![test_mapping("read", "ref", 0, 100, 100, 1.0, 100)];
        mappings[0].candidate_boundary = CandidateBoundary::from_bits(1).unwrap();
        limit_mappings(&mut mappings, 2);
        assert!(mappings[0].candidate_boundary.profile_incomplete());
    }

    fn default_long_read_test_params() -> MapParams {
        MapParams {
            top_n: 1,
            score_ratio: 0.8,
            min_seeds: 2,
            min_identity: 0.9,
            min_mapped_frac: 0.2,
            min_mapped_len: None,
            no_identity_filter: false,
            ignore_n: false,
            identity_mode: IdentityMode::Gapped,
            diag_bin_size: 128,
            extension_max_mismatches: 40,
            full_read_interval: true,
            fast_long_read_identity: true,
            long_read_chaining: true,
            long_read_extension: LongReadExtension::Retained,
            query_seed_stride: 1,
            long_read_rescue: true,
            long_read_rescue_limit: None,
        }
    }

    fn test_mapping(
        read_id: &str,
        ref_id: &str,
        ref_start: u32,
        ref_end: u32,
        mapped_bases: usize,
        identity: f64,
        score: i64,
    ) -> Mapping {
        Mapping {
            read_id: read_id.to_string(),
            ref_id: ref_id.to_string(),
            ref_start,
            ref_end,
            identity,
            read_len: 100,
            mapped_bases,
            mate_ref_status: MateRefStatus::NotApplicable,
            strand: 1,
            seed_count: 1,
            score,
            second_score: 0,
            rank: 0,
            num_mappings: 0,
            matched_bases: mapped_bases,
            query_start: 0,
            query_end: mapped_bases,
            pair_score: 0,
            adjusted_score: score,
            template_length: 0,
            pair_status: PairStatus::Single,
            mate_number: 0,
            pair_id: 0,
            candidate_groups: ReferenceCandidateGroups::new(1, 1),
            candidate_boundary: CandidateBoundary::COMPLETE,
        }
    }

    fn with_groups(mut mapping: Mapping, profile_taxon: u32, terminal_taxon: u32) -> Mapping {
        mapping.candidate_groups = ReferenceCandidateGroups::new(profile_taxon, terminal_taxon);
        mapping
    }
}

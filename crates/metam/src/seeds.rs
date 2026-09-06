// Strobealign-derived seeding for Metam.
//
// This module adapts the syncmer and randstrobe algorithms from the local
// MIT-licensed Strobealign Rust codebase (`strobealign/src/seeding/*`) for
// Metam's no-CIGAR coordinate mapping index. The implementation is integrated
// with Metam's index format and omits Strobealign-specific profile plumbing.

use std::{cmp::min, collections::VecDeque};

use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};

use crate::sequence::base_code;

const PRIME_1: u64 = 0x9E37_79B1_85EB_CA87;
const PRIME_2: u64 = 0xC2B2_AE3D_27D4_EB4F;
const PRIME_3: u64 = 0x1656_67B1_9E37_79F9;
const PRIME_4: u64 = 0x85EB_CA77_C2B2_AE63;
const PRIME_5: u64 = 0x27D4_EB2F_1656_67C5;
const STROBE2_OFFSET_BITS: u32 = 8;
const DEFAULT_AUX_LEN: u8 = 17;
pub const DEFAULT_STROBE_START_DIST_MIN: u16 = 30;
pub const DEFAULT_STROBE_START_DIST_MAX: u16 = 90;
pub const MAX_STROBE_START_DIST: u16 = 255;

/// Meaningful bits in the final combined randstrobe hash.
///
/// The index relies on this mask being contiguous and anchored at bit 63 so it
/// can split the hash into a high-bit lookup key and the remaining fingerprint.
pub(crate) const RANDSTROBE_HASH_MASK: u64 = 0xFFFF_FFFF_FFFF_FE00;
pub(crate) const RANDSTROBE_HASH_LOW_ZERO_BITS: u8 = RANDSTROBE_HASH_MASK.trailing_zeros() as u8;
pub(crate) const RANDSTROBE_HASH_SPAN_BITS: u8 = u64::BITS as u8 - RANDSTROBE_HASH_LOW_ZERO_BITS;

const _: () = {
    assert!(RANDSTROBE_HASH_MASK != 0);
    assert!(RANDSTROBE_HASH_MASK == u64::MAX << RANDSTROBE_HASH_MASK.trailing_zeros());
};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum SeedMode {
    Minimizer,
    Randstrobe,
}

impl SeedMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Minimizer => "minimizer",
            Self::Randstrobe => "randstrobe",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Seed {
    pub hash: u64,
    pub pos: u32,
    pub span: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct StrobeStartDistance {
    pub min: u16,
    pub max: u16,
}

impl StrobeStartDistance {
    pub fn new(min: u16, max: u16) -> Result<Self> {
        if min == 0 {
            bail!("--strobe-start-dist-min must be greater than 0");
        }
        if max < min {
            bail!(
                "--strobe-start-dist-max must be greater than or equal to --strobe-start-dist-min"
            );
        }
        if max > MAX_STROBE_START_DIST {
            bail!("--strobe-start-dist-max must be at most {MAX_STROBE_START_DIST}");
        }
        Ok(Self { min, max })
    }
}

impl Default for StrobeStartDistance {
    fn default() -> Self {
        Self {
            min: DEFAULT_STROBE_START_DIST_MIN,
            max: DEFAULT_STROBE_START_DIST_MAX,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct SyncmerParameters {
    k: usize,
    s: usize,
    t: usize,
}

#[derive(Clone, Copy, Debug)]
struct RandstrobeParameters {
    window: Strobe2Window,
    q: u64,
    main_hash_mask: u64,
    forward_main_hash_mask: u64,
    partial_orientation_pos: u32,
}

#[derive(Clone, Copy, Debug)]
enum Strobe2Window {
    LegacySyncmerRank {
        w_min: usize,
        w_max: usize,
        max_dist: u16,
    },
    StartDistance(StrobeStartDistance),
}

#[derive(Clone, Copy, Debug)]
struct Syncmer {
    hash: u64,
    is_forward: bool,
    position: usize,
}

#[derive(Clone, Copy, Debug)]
struct Randstrobe {
    hash: u64,
    strobe1_pos: usize,
    strobe2_pos: usize,
}

pub fn validate_seed_params(k: usize, w: usize) -> Result<()> {
    if !(8..=32).contains(&k) {
        bail!("-k must be between 8 and 32 for Metam's Strobealign-derived seeds");
    }
    if w == 0 || w > 255 {
        bail!("-w must be between 1 and 255");
    }
    Ok(())
}

pub fn seeds(seq: &[u8], k: usize, w: usize, mode: SeedMode) -> Vec<Seed> {
    seeds_with_strobe_start_distance(seq, k, w, mode, None)
}

pub fn seeds_with_strobe_start_distance(
    seq: &[u8],
    k: usize,
    w: usize,
    mode: SeedMode,
    strobe_start_dist: Option<StrobeStartDistance>,
) -> Vec<Seed> {
    match mode {
        SeedMode::Minimizer => minimizer_seeds(seq, k, w),
        SeedMode::Randstrobe => {
            randstrobe_seeds_with_strobe_start_distance(seq, k, w, strobe_start_dist)
        }
    }
}

pub fn randstrobe_seeds(seq: &[u8], k: usize, w: usize) -> Vec<Seed> {
    randstrobe_seeds_with_strobe_start_distance(seq, k, w, None)
}

pub fn randstrobe_seeds_with_strobe_start_distance(
    seq: &[u8],
    k: usize,
    w: usize,
    strobe_start_dist: Option<StrobeStartDistance>,
) -> Vec<Seed> {
    if seq.len() < k {
        return Vec::new();
    }
    let syncmer = syncmer_parameters(k, w);
    let randstrobe = randstrobe_parameters(k, w, strobe_start_dist);
    let syncmers = SyncmerIterator::new(seq, syncmer);
    let raw = RandstrobeIterator::new(syncmers, randstrobe)
        .map(|seed| Seed {
            hash: seed.hash,
            pos: seed.strobe1_pos as u32,
            span: (seed.strobe2_pos + k - seed.strobe1_pos) as u32,
        })
        .collect();
    thin_position_windows(raw, sparse_window(w))
}

pub fn minimizer_seeds(seq: &[u8], k: usize, w: usize) -> Vec<Seed> {
    if seq.len() < k || w == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut last_pos = None;
    let mut deque: VecDeque<(usize, u64)> = VecDeque::new();
    let hashes = rolling_canonical_hashes(seq, k);

    for (i, hash) in hashes.iter().enumerate() {
        if let Some(hash) = hash {
            while deque.back().map(|(_, h)| *h > *hash).unwrap_or(false) {
                deque.pop_back();
            }
            deque.push_back((i, *hash));
        }
        while deque.front().map(|(pos, _)| pos + w <= i).unwrap_or(false) {
            deque.pop_front();
        }
        if i + 1 >= w {
            if let Some((pos, hash)) = deque.front().copied() {
                if last_pos != Some(pos) {
                    out.push(Seed {
                        hash,
                        pos: pos as u32,
                        span: k as u32,
                    });
                    last_pos = Some(pos);
                }
            }
        }
    }
    out
}

fn syncmer_parameters(k: usize, w: usize) -> SyncmerParameters {
    let delta = (w / 2 + 2).clamp(2, 9).min(k.saturating_sub(2));
    let delta = if delta.is_multiple_of(2) {
        delta
    } else {
        delta + 1
    }
    .min(k - 1);
    let s = k - delta;
    SyncmerParameters {
        k,
        s,
        t: (k - s) / 2 + 1,
    }
}

fn thin_position_windows(seeds: Vec<Seed>, window: usize) -> Vec<Seed> {
    if seeds.is_empty() || window <= 1 {
        return seeds;
    }
    let mut out = Vec::with_capacity(seeds.len().div_ceil(window).max(1));
    let mut current_bin: Option<usize> = None;
    let mut best: Option<Seed> = None;
    for seed in seeds {
        let bin = seed.pos as usize / window;
        if current_bin != Some(bin) {
            if let Some(seed) = best.take() {
                out.push(seed);
            }
            current_bin = Some(bin);
            best = Some(seed);
        } else if best
            .as_ref()
            .map(|prev| seed.hash < prev.hash || (seed.hash == prev.hash && seed.pos < prev.pos))
            .unwrap_or(true)
        {
            best = Some(seed);
        }
    }
    if let Some(seed) = best {
        out.push(seed);
    }
    out
}

fn sparse_window(w: usize) -> usize {
    w.saturating_mul(2).div_ceil(3).max(1)
}

fn randstrobe_parameters(
    k: usize,
    w: usize,
    strobe_start_dist: Option<StrobeStartDistance>,
) -> RandstrobeParameters {
    let q = (1u64 << 8) - 1;
    let main_hash_mask = !0u64 << (9 + DEFAULT_AUX_LEN);
    let partial_orientation_pos = main_hash_mask.trailing_zeros() - 1;
    let window = if let Some(distance) = strobe_start_dist {
        Strobe2Window::StartDistance(distance)
    } else {
        let w_min = (w / 2).max(1);
        let w_max = w.max(w_min);
        let max_dist = k.saturating_mul(6).max(w_max + k).min(255) as u16;
        Strobe2Window::LegacySyncmerRank {
            w_min,
            w_max,
            max_dist,
        }
    };
    RandstrobeParameters {
        window,
        q,
        main_hash_mask,
        forward_main_hash_mask: main_hash_mask | (1u64 << partial_orientation_pos),
        partial_orientation_pos,
    }
}

struct SyncmerIterator<'a> {
    seq: &'a [u8],
    params: SyncmerParameters,
    kmask: u64,
    smask: u64,
    kshift: usize,
    sshift: usize,
    qs: VecDeque<u64>,
    qs_min_val: u64,
    valid_len: usize,
    xk: [u64; 2],
    xs: [u64; 2],
    i: usize,
    exhausted: bool,
}

impl<'a> SyncmerIterator<'a> {
    fn new(seq: &'a [u8], params: SyncmerParameters) -> Self {
        Self {
            seq,
            params,
            kmask: if params.k == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * params.k)) - 1
            },
            smask: (1u64 << (2 * params.s)) - 1,
            kshift: (params.k - 1) * 2,
            sshift: (params.s - 1) * 2,
            qs: VecDeque::new(),
            qs_min_val: u64::MAX,
            valid_len: 0,
            xk: [0, 0],
            xs: [0, 0],
            i: 0,
            exhausted: false,
        }
    }
}

impl Iterator for SyncmerIterator<'_> {
    type Item = Syncmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }
        for i in self.i..self.seq.len() {
            if let Some(c) = base_code(self.seq[i]) {
                self.xk[0] = ((self.xk[0] << 2) | u64::from(c)) & self.kmask;
                self.xk[1] = (self.xk[1] >> 2) | (u64::from(3 - c) << self.kshift);
                self.xs[0] = ((self.xs[0] << 2) | u64::from(c)) & self.smask;
                self.xs[1] = (self.xs[1] >> 2) | (u64::from(3 - c) << self.sshift);
                self.valid_len += 1;
                if self.valid_len < self.params.s {
                    continue;
                }

                let ys = min(self.xs[0], self.xs[1]);
                let hash_s = xxh64(ys);
                self.qs.push_back(hash_s);
                let window = self.params.k - self.params.s + 1;
                if self.qs.len() < window {
                    continue;
                }
                if self.qs.len() == window {
                    self.qs_min_val = self.qs.iter().copied().min().unwrap_or(u64::MAX);
                } else {
                    let front = self.qs.pop_front().expect("non-empty syncmer queue");
                    if front == self.qs_min_val {
                        self.qs_min_val = self.qs.iter().copied().min().unwrap_or(u64::MAX);
                    } else if hash_s < self.qs_min_val {
                        self.qs_min_val = hash_s;
                    }
                }

                if self.qs[self.params.t - 1] == self.qs_min_val {
                    let yk = min(self.xk[0], self.xk[1]);
                    let syncmer = Syncmer {
                        hash: xxh64(yk),
                        is_forward: self.xk[0] <= self.xk[1],
                        position: i + 1 - self.params.k,
                    };
                    self.i = i + 1;
                    return Some(syncmer);
                }
            } else {
                self.qs_min_val = u64::MAX;
                self.valid_len = 0;
                self.xs = [0, 0];
                self.xk = [0, 0];
                self.qs.clear();
            }
        }
        self.exhausted = true;
        None
    }
}

struct RandstrobeIterator<I: Iterator<Item = Syncmer>> {
    params: RandstrobeParameters,
    syncmers: VecDeque<Syncmer>,
    inner: I,
}

impl<I: Iterator<Item = Syncmer>> RandstrobeIterator<I> {
    fn new(inner: I, params: RandstrobeParameters) -> Self {
        Self {
            params,
            syncmers: VecDeque::new(),
            inner,
        }
    }

    fn fill_lookahead(&mut self) {
        match self.params.window {
            Strobe2Window::LegacySyncmerRank { w_max, .. } => {
                while self.syncmers.len() <= w_max {
                    if let Some(syncmer) = self.inner.next() {
                        self.syncmers.push_back(syncmer);
                    } else {
                        break;
                    }
                }
            }
            Strobe2Window::StartDistance(distance) => {
                if self.syncmers.is_empty() {
                    if let Some(syncmer) = self.inner.next() {
                        self.syncmers.push_back(syncmer);
                    }
                }
                let Some(strobe1) = self.syncmers.front().copied() else {
                    return;
                };
                let max_position = strobe1.position.saturating_add(usize::from(distance.max));
                while self
                    .syncmers
                    .back()
                    .map(|syncmer| syncmer.position <= max_position)
                    .unwrap_or(true)
                {
                    let Some(syncmer) = self.inner.next() else {
                        break;
                    };
                    let past_window = syncmer.position > max_position;
                    self.syncmers.push_back(syncmer);
                    if past_window {
                        break;
                    }
                }
            }
        }
    }
}

impl<I: Iterator<Item = Syncmer>> Iterator for RandstrobeIterator<I> {
    type Item = Randstrobe;

    fn next(&mut self) -> Option<Self::Item> {
        self.fill_lookahead();
        if self.syncmers.is_empty() {
            return None;
        }

        let strobe1 = self.syncmers[0];
        let mut min_val = u64::MAX;
        let mut strobe2 = strobe1;
        let first_candidate = match self.params.window {
            Strobe2Window::LegacySyncmerRank { w_min, .. } => w_min,
            Strobe2Window::StartDistance(_) => 1,
        };

        for i in first_candidate..self.syncmers.len() {
            let candidate = self.syncmers[i];
            let distance = candidate.position.saturating_sub(strobe1.position);
            match self.params.window {
                Strobe2Window::LegacySyncmerRank { max_dist, .. } => {
                    if distance > usize::from(max_dist) {
                        break;
                    }
                }
                Strobe2Window::StartDistance(strobe_start_dist) => {
                    if distance < usize::from(strobe_start_dist.min) {
                        continue;
                    }
                    if distance > usize::from(strobe_start_dist.max) {
                        break;
                    }
                }
            }
            let selector = (strobe1.hash ^ candidate.hash) & self.params.q;
            let ones = selector.count_ones() as u64;
            if ones < min_val {
                min_val = ones;
                strobe2 = candidate;
            }
        }
        self.syncmers.pop_front();

        Some(Randstrobe {
            hash: randstrobe_hash(strobe1, strobe2, &self.params),
            strobe1_pos: strobe1.position,
            strobe2_pos: strobe2.position,
        })
    }
}

fn randstrobe_hash(strobe1: Syncmer, strobe2: Syncmer, params: &RandstrobeParameters) -> u64 {
    let is_forward1 = u64::from(strobe1.is_forward);
    let is_forward2 = u64::from(strobe2.is_forward);
    ((strobe1.hash & params.main_hash_mask)
        | (strobe2.hash & !params.forward_main_hash_mask)
        | (is_forward1 << params.partial_orientation_pos)
        | (is_forward2 << STROBE2_OFFSET_BITS))
        & RANDSTROBE_HASH_MASK
}

fn rolling_canonical_hashes(seq: &[u8], k: usize) -> Vec<Option<u64>> {
    let n = seq.len().saturating_sub(k) + 1;
    let mut out = Vec::with_capacity(n);
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let rc_shift = 2 * (k - 1);
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let mut valid = 0usize;

    for (i, &base) in seq.iter().enumerate() {
        if let Some(code) = base_code(base) {
            fwd = ((fwd << 2) | u64::from(code)) & mask;
            rev = (rev >> 2) | (u64::from(3 - code) << rc_shift);
            valid += 1;
        } else {
            fwd = 0;
            rev = 0;
            valid = 0;
        }
        if i + 1 >= k {
            out.push((valid >= k).then(|| xxh64(fwd.min(rev))));
        }
    }
    out
}

#[inline]
pub fn xxh64(input: u64) -> u64 {
    let mut result = PRIME_5.wrapping_add(8);
    result ^= input
        .wrapping_mul(PRIME_2)
        .rotate_left(31)
        .wrapping_mul(PRIME_1);
    result = result
        .rotate_left(27)
        .wrapping_mul(PRIME_1)
        .wrapping_add(PRIME_4);
    result ^= result >> 33;
    result = result.wrapping_mul(PRIME_2);
    result ^= result >> 29;
    result = result.wrapping_mul(PRIME_3);
    result ^= result >> 32;
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::revcomp;

    #[test]
    fn minimizer_hashes_are_reverse_complement_stable() {
        let a = minimizer_seeds(b"ACGTACGTACGTACGT", 8, 4);
        let b = minimizer_seeds(&revcomp(b"ACGTACGTACGTACGT"), 8, 4);
        assert_eq!(
            a.iter().map(|s| s.hash).collect::<Vec<_>>(),
            b.iter().map(|s| s.hash).collect::<Vec<_>>()
        );
    }

    #[test]
    fn randstrobes_have_positions() {
        let s = randstrobe_seeds(b"ACGTACGTACGTACGTACGTACGTACGT", 10, 5);
        assert!(!s.is_empty());
        assert!(s.iter().all(|x| x.span >= 10));
    }

    #[test]
    fn randstrobe_start_distance_bounds_seed_span() {
        let mut seq = Vec::with_capacity(2_000);
        let mut state = 0x6eed_5eed_1234_5678u64;
        for i in 0..2_000 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407 + i as u64);
            seq.push(match (state >> 62) & 3 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            });
        }
        let k = 15;
        let distance = StrobeStartDistance::new(30, 90).unwrap();
        let seeds = randstrobe_seeds_with_strobe_start_distance(&seq, k, 15, Some(distance));
        let paired = seeds.iter().filter(|seed| seed.span > k as u32).count();

        assert!(paired > 20, "expected paired randstrobes, got {paired}");
        assert!(seeds.iter().all(|seed| {
            seed.span == k as u32
                || ((u32::from(distance.min) + k as u32)..=(u32::from(distance.max) + k as u32))
                    .contains(&seed.span)
        }));
    }
}

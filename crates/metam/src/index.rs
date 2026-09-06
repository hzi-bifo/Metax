use std::{
    collections::HashSet,
    fs::{self, File, OpenOptions},
    io::{copy, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    process,
    time::{SystemTime, UNIX_EPOCH},
};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::{
    io::{open_reader, FastaReader, Record},
    seeds::{
        seeds_with_strobe_start_distance, validate_seed_params, Seed, SeedMode,
        StrobeStartDistance, RANDSTROBE_HASH_LOW_ZERO_BITS, RANDSTROBE_HASH_SPAN_BITS,
    },
    sequence::{code_base, pack_dna_2bit, unpack_dna_2bit, BitMask},
};

pub const INDEX_MAGIC: [u8; 8] = *b"METAMM02";
pub const INDEX_VERSION: u32 = 6;
const MIN_SUPPORTED_INDEX_VERSION: u32 = 4;
const INDEX_CACHE_MAGIC: [u8; 8] = *b"METAMMC1";
const INDEX_CACHE_VERSION: u32 = 5;
const INDEX_CACHE_IO_BUFFER_BYTES: usize = 8 * 1024 * 1024;
const INDEX_SEQUENCE_CHUNK: usize = 1_000_000;
const MAX_OPEN_PARTITION_WRITERS: usize = 4_096;
const OPEN_FILE_RESERVE: usize = 64;
const FALLBACK_OPEN_FILE_LIMIT: usize = 256;
pub const DEFAULT_INDEX_SEED_CHUNK_MIB: usize = 1024;
pub const DEFAULT_INDEX_SEED_CHUNK_BYTES: usize = DEFAULT_INDEX_SEED_CHUNK_MIB * 1024 * 1024;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Reference {
    pub id: String,
    pub len: u64,
    pub packed_bases: Vec<u8>,
    pub ambiguous: BitMask,
}

impl Reference {
    pub fn from_record(record: &Record) -> Self {
        Self {
            id: record.id.clone(),
            len: record.seq.len() as u64,
            packed_bases: pack_dna_2bit(&record.seq),
            ambiguous: BitMask::from_ambiguous(&record.seq),
        }
    }

    pub fn sequence(&self) -> Vec<u8> {
        unpack_dna_2bit(&self.packed_bases, self.len as usize, &self.ambiguous)
    }

    /// Two-bit code of the base at `pos`, or `None` where the reference is
    /// ambiguous.
    ///
    /// Comparison callers want the code, not the letter. Routing them through
    /// [`Self::base_at`] would spell the code out as ASCII only for
    /// `base_code` to read it straight back, once for every base of every
    /// candidate considered.
    pub fn code_at(&self, pos: usize) -> Option<u8> {
        if self.ambiguous.get(pos) {
            None
        } else {
            Some((self.packed_bases[pos / 4] >> (6 - 2 * (pos % 4))) & 3)
        }
    }

    pub fn base_at(&self, pos: usize) -> u8 {
        match self.code_at(pos) {
            Some(code) => code_base(code),
            None => b'N',
        }
    }

    pub fn len(&self) -> usize {
        self.len as usize
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedOccurrence {
    pub ref_id: u32,
    pub pos: u32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SeedTable {
    pub hash_bits: u8,
    pub hashes32: Vec<u32>,
    pub hashes64: Vec<u64>,
    pub occurrences: Vec<u64>,
    #[serde(default)]
    fingerprint_bits: u8,
    #[serde(default)]
    fingerprints: Vec<u8>,
    lookup_index: SeedLookupIndex,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
struct SeedLookupIndex {
    prefix_bits: u8,
    offsets: Vec<u64>,
}

fn build_seed_lookup_index<F>(hash_bits: u8, len: usize, mut hash_at: F) -> SeedLookupIndex
where
    F: FnMut(usize) -> u64,
{
    if len == 0 {
        return SeedLookupIndex::default();
    }
    let prefix_bits = hash_bits.min(20);
    let bucket_count = 1usize << prefix_bits;
    let mut offsets = vec![0u64; bucket_count + 1];
    let suffix_bits = hash_bits - prefix_bits;
    let mut row = 0usize;
    for bucket in 0..bucket_count {
        if bucket + 1 == bucket_count {
            row = len;
        } else {
            let upper = ((bucket + 1) as u64) << suffix_bits;
            while row < len && hash_at(row) < upper {
                row += 1;
            }
        }
        offsets[bucket + 1] = row as u64;
    }
    SeedLookupIndex {
        prefix_bits,
        offsets,
    }
}

impl SeedTable {
    pub fn new(hash_bits: u8) -> Self {
        Self::new_with_exact_randstrobe_match(hash_bits, false)
    }

    fn new_with_exact_randstrobe_match(hash_bits: u8, exact_randstrobe_match: bool) -> Self {
        Self {
            hash_bits,
            hashes32: Vec::new(),
            hashes64: Vec::new(),
            occurrences: Vec::new(),
            fingerprint_bits: if exact_randstrobe_match {
                randstrobe_fingerprint_bits(hash_bits)
            } else {
                0
            },
            fingerprints: Vec::new(),
            lookup_index: SeedLookupIndex::default(),
        }
    }

    pub fn push_occurrence_group(
        &mut self,
        hash: u64,
        occurrences: &[SeedOccurrence],
    ) -> Result<()> {
        if occurrences.len() > u32::MAX as usize {
            bail!("seed occurrence group exceeds index format limit");
        }
        self.push_hash_copies(hash, occurrences.len())?;
        self.occurrences.extend(
            occurrences
                .iter()
                .map(|occurrence| pack_occurrence(occurrence.ref_id, occurrence.pos)),
        );
        Ok(())
    }

    fn push_group(&mut self, hash: u64, seeds: &[RawSeed]) -> Result<()> {
        if seeds.len() > u32::MAX as usize {
            bail!("seed occurrence group exceeds index format limit");
        }
        self.push_hash_copies(hash, seeds.len())?;
        if self.fingerprint_bits > 0 {
            for seed in seeds {
                push_packed_fingerprint(
                    &mut self.fingerprints,
                    randstrobe_fingerprint(seed.hash, self.hash_bits),
                    self.fingerprint_bits,
                );
            }
        }
        self.occurrences.extend(
            seeds
                .iter()
                .map(|seed| pack_occurrence(seed.ref_id, seed.pos)),
        );
        Ok(())
    }

    fn push_hash_copies(&mut self, hash: u64, count: usize) -> Result<()> {
        if self.hash_bits <= 32 {
            if hash > u64::from(u32::MAX) {
                bail!("seed hash exceeds 32-bit compact table");
            }
            self.hashes32
                .extend(std::iter::repeat_n(hash as u32, count));
        } else {
            self.hashes64.extend(std::iter::repeat_n(hash, count));
        }
        Ok(())
    }

    pub fn occurrence_count(&self) -> usize {
        self.occurrences.len()
    }

    pub fn group_count(&self) -> usize {
        if self.hash_bits <= 32 {
            count_groups_u32(&self.hashes32)
        } else {
            count_groups_u64(&self.hashes64)
        }
    }

    fn hash_count(&self) -> usize {
        if self.hash_bits <= 32 {
            self.hashes32.len()
        } else {
            self.hashes64.len()
        }
    }

    pub fn rebuild_lookup_index(&mut self) {
        let len = self.occurrence_count();
        let lookup_index = build_seed_lookup_index(self.hash_bits, len, |row| self.hash_at(row));
        self.lookup_index = lookup_index;
    }

    #[cfg(test)]
    fn lookup(&self, key: u64, hash: u64) -> SeedOccurrences<'_> {
        self.lookup_with_mode(key, hash, SeedLookupMode::IndexDefault)
    }

    fn lookup_with_mode(&self, key: u64, hash: u64, mode: SeedLookupMode) -> SeedOccurrences<'_> {
        if self.hash_bits <= 32 {
            let Ok(key) = u32::try_from(key) else {
                return empty_seed_occurrences();
            };
            let (lo, hi) = self.lookup_bounds(u64::from(key));
            let hashes = &self.hashes32[lo..hi];
            let start = lo + hashes.partition_point(|h| *h < key);
            let end = lo + hashes.partition_point(|h| *h <= key);
            self.seed_occurrences(hash, start, end, mode)
        } else {
            let (lo, hi) = self.lookup_bounds(key);
            let hashes = &self.hashes64[lo..hi];
            let start = lo + hashes.partition_point(|h| *h < key);
            let end = lo + hashes.partition_point(|h| *h <= key);
            self.seed_occurrences(hash, start, end, mode)
        }
    }

    fn seed_occurrences(
        &self,
        hash: u64,
        start: usize,
        end: usize,
        mode: SeedLookupMode,
    ) -> SeedOccurrences<'_> {
        let fingerprint_bits = if mode.use_fingerprint_filter() {
            self.fingerprint_bits
        } else {
            0
        };
        let width = fingerprint_bytes(fingerprint_bits);
        let fingerprints = if width == 0 {
            &[]
        } else {
            &self.fingerprints[start * width..end * width]
        };
        SeedOccurrences {
            data: &self.occurrences[start..end],
            fingerprints,
            fingerprint_bits,
            target_fingerprint: if fingerprint_bits == 0 {
                0
            } else {
                randstrobe_fingerprint(hash, self.hash_bits)
            },
        }
    }

    fn lookup_bounds(&self, key: u64) -> (usize, usize) {
        if self.lookup_index.offsets.is_empty() {
            return (0, self.occurrence_count());
        }
        let prefix_bits = self.lookup_index.prefix_bits;
        let suffix_bits = self.hash_bits - prefix_bits;
        let bucket = (key >> suffix_bits) as usize;
        (
            self.lookup_index.offsets[bucket] as usize,
            self.lookup_index.offsets[bucket + 1] as usize,
        )
    }

    fn hash_at(&self, i: usize) -> u64 {
        if self.hash_bits <= 32 {
            u64::from(self.hashes32[i])
        } else {
            self.hashes64[i]
        }
    }

    fn fingerprint_bytes_per_occurrence(&self) -> usize {
        fingerprint_bytes(self.fingerprint_bits)
    }

    fn expected_fingerprint_bytes(&self) -> Option<usize> {
        self.occurrence_count()
            .checked_mul(self.fingerprint_bytes_per_occurrence())
    }
}

impl Default for SeedTable {
    fn default() -> Self {
        Self::new(64)
    }
}

fn count_groups_u32(values: &[u32]) -> usize {
    values
        .iter()
        .enumerate()
        .filter(|(i, value)| *i == 0 || values[*i - 1] != **value)
        .count()
}

fn count_groups_u64(values: &[u64]) -> usize {
    values
        .iter()
        .enumerate()
        .filter(|(i, value)| *i == 0 || values[*i - 1] != **value)
        .count()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CacheSeedKeyLayout {
    /// One key row per occurrence. This is smallest for very sparse tables.
    Repeated = 0,
    /// One key row per distinct key, plus offsets into the occurrence array.
    Grouped = 1,
}

impl CacheSeedKeyLayout {
    fn from_code(code: u8) -> Result<Self> {
        match code {
            0 => Ok(Self::Repeated),
            1 => Ok(Self::Grouped),
            _ => bail!("invalid seed-key layout in index cache"),
        }
    }
}

enum CompactGroupOffsets {
    U32(Vec<u32>),
    U64(Vec<u64>),
}

impl CompactGroupOffsets {
    fn with_capacity(occurrence_count: usize, group_count: usize) -> Result<Self> {
        let capacity = group_count
            .checked_add(1)
            .context("index cache seed-group count overflow")?;
        if u32::try_from(occurrence_count).is_ok() {
            Ok(Self::U32(Vec::with_capacity(capacity)))
        } else {
            Ok(Self::U64(Vec::with_capacity(capacity)))
        }
    }

    fn push(&mut self, offset: usize) -> Result<()> {
        match self {
            Self::U32(offsets) => offsets
                .push(u32::try_from(offset).context("index cache seed-group offset exceeds u32")?),
            Self::U64(offsets) => offsets.push(offset as u64),
        }
        Ok(())
    }

    #[cfg(test)]
    fn width(&self) -> u8 {
        match self {
            Self::U32(_) => 4,
            Self::U64(_) => 8,
        }
    }

    #[cfg(test)]
    fn len(&self) -> usize {
        match self {
            Self::U32(offsets) => offsets.len(),
            Self::U64(offsets) => offsets.len(),
        }
    }
}

#[cfg(test)]
struct CompactSeedKeys {
    hashes32: Vec<u32>,
    hashes64: Vec<u64>,
    group_offsets: CompactGroupOffsets,
    lookup_index: SeedLookupIndex,
}

#[cfg(test)]
impl CompactSeedKeys {
    fn from_seed_table(seed_table: &SeedTable, expected_group_count: usize) -> Result<Self> {
        let occurrence_count = seed_table.occurrence_count();
        let mut group_offsets =
            CompactGroupOffsets::with_capacity(occurrence_count, expected_group_count)?;
        group_offsets.push(0)?;

        let (hashes32, hashes64) = if seed_table.hash_bits <= 32 {
            (
                compact_repeated_keys(
                    &seed_table.hashes32,
                    expected_group_count,
                    &mut group_offsets,
                )?,
                Vec::new(),
            )
        } else {
            (
                Vec::new(),
                compact_repeated_keys(
                    &seed_table.hashes64,
                    expected_group_count,
                    &mut group_offsets,
                )?,
            )
        };
        if occurrence_count > 0 {
            group_offsets.push(occurrence_count)?;
        }

        let key_count = if seed_table.hash_bits <= 32 {
            hashes32.len()
        } else {
            hashes64.len()
        };
        let lookup_index = build_seed_lookup_index(seed_table.hash_bits, key_count, |row| {
            if seed_table.hash_bits <= 32 {
                u64::from(hashes32[row])
            } else {
                hashes64[row]
            }
        });
        Ok(Self {
            hashes32,
            hashes64,
            group_offsets,
            lookup_index,
        })
    }
}

#[cfg(test)]
fn compact_repeated_keys<T: Copy + Ord>(
    values: &[T],
    expected_group_count: usize,
    group_offsets: &mut CompactGroupOffsets,
) -> Result<Vec<T>> {
    let mut compact = Vec::with_capacity(expected_group_count);
    let mut previous = None;
    for (row, &value) in values.iter().enumerate() {
        if let Some(previous) = previous {
            if value < previous {
                bail!("index seed hash rows are not sorted");
            }
            if value == previous {
                continue;
            }
            group_offsets.push(row)?;
        }
        compact.push(value);
        previous = Some(value);
    }
    if compact.len() != expected_group_count {
        bail!("index metadata/hash group count mismatch");
    }
    Ok(compact)
}

#[cfg(test)]
enum CacheSeedKeys {
    Repeated,
    Grouped(CompactSeedKeys),
}

fn choose_cache_seed_key_layout(
    hash_bits: u8,
    occurrence_count: usize,
    group_count: usize,
) -> CacheSeedKeyLayout {
    let key_width = if hash_bits <= 32 {
        std::mem::size_of::<u32>()
    } else {
        std::mem::size_of::<u64>()
    };
    let group_offset_width = if u32::try_from(occurrence_count).is_ok() {
        std::mem::size_of::<u32>()
    } else {
        std::mem::size_of::<u64>()
    };
    let repeated_bytes = occurrence_count as u128 * key_width as u128;
    let grouped_bytes = group_count as u128 * key_width as u128
        + group_count.saturating_add(1) as u128 * group_offset_width as u128;

    if group_count > 0 && grouped_bytes < repeated_bytes {
        CacheSeedKeyLayout::Grouped
    } else {
        CacheSeedKeyLayout::Repeated
    }
}

#[cfg(test)]
impl CacheSeedKeys {
    fn for_index(index: &Index) -> Result<Self> {
        let occurrence_count = index.seed_table.occurrence_count();
        let group_count = usize::try_from(index.header.indexed_seed_hashes)
            .context("index cache distinct seed-key count exceeds address space")?;

        // Grouping is a large win for common/repetitive keys, but can enlarge
        // singleton-heavy tables because every distinct key needs an offset.
        // Choose from exact array sizes so the cache never grows for this
        // optimization alone.
        if choose_cache_seed_key_layout(index.seed_table.hash_bits, occurrence_count, group_count)
            == CacheSeedKeyLayout::Grouped
        {
            Ok(Self::Grouped(CompactSeedKeys::from_seed_table(
                &index.seed_table,
                group_count,
            )?))
        } else {
            Ok(Self::Repeated)
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SeedOccurrences<'a> {
    data: &'a [u64],
    fingerprints: &'a [u8],
    fingerprint_bits: u8,
    target_fingerprint: u32,
}

fn empty_seed_occurrences() -> SeedOccurrences<'static> {
    SeedOccurrences {
        data: &[],
        fingerprints: &[],
        fingerprint_bits: 0,
        target_fingerprint: 0,
    }
}

impl SeedOccurrences<'_> {
    pub fn len(&self) -> usize {
        if self.fingerprint_bits == 0 {
            self.data.len()
        } else {
            (0..self.data.len())
                .filter(|&i| {
                    packed_fingerprint(self.fingerprints, i, self.fingerprint_bits)
                        == self.target_fingerprint
                })
                .count()
        }
    }

    pub fn is_empty(&self) -> bool {
        self.iter().next().is_none()
    }

    pub fn iter(&self) -> SeedOccurrenceIter<'_> {
        SeedOccurrenceIter {
            inner: self.data.iter(),
            fingerprints: self.fingerprints,
            fingerprint_bits: self.fingerprint_bits,
            target_fingerprint: self.target_fingerprint,
            index: 0,
        }
    }
}

pub struct SeedOccurrenceIter<'a> {
    inner: std::slice::Iter<'a, u64>,
    fingerprints: &'a [u8],
    fingerprint_bits: u8,
    target_fingerprint: u32,
    index: usize,
}

impl Iterator for SeedOccurrenceIter<'_> {
    type Item = SeedOccurrence;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let packed = *self.inner.next()?;
            let index = self.index;
            self.index += 1;
            if self.fingerprint_bits == 0
                || packed_fingerprint(self.fingerprints, index, self.fingerprint_bits)
                    == self.target_fingerprint
            {
                return Some(unpack_occurrence(packed));
            }
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Index {
    pub header: IndexHeader,
    pub references: Vec<Reference>,
    pub seed_table: SeedTable,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IndexHeader {
    pub magic: [u8; 8],
    pub version: u32,
    pub metam_version: String,
    pub k: usize,
    pub w: usize,
    pub seed_mode: SeedMode,
    pub max_occ: usize,
    pub hash_bits: u8,
    #[serde(default)]
    pub exact_randstrobe_match: bool,
    #[serde(default)]
    pub strobe_start_dist: Option<StrobeStartDistance>,
    pub reference_count: usize,
    pub total_reference_bases: u64,
    pub indexed_seed_hashes: u64,
    pub indexed_seed_occurrences: u64,
    pub filtered_seed_hashes: u64,
    pub reference_storage_bytes: u64,
    pub created_unix_seconds: u64,
}

#[derive(Clone, Copy, Debug)]
pub enum BuildProgress {
    Phase(&'static str),
    ReferencesRead(usize),
    ReferencesSeeded(usize, usize),
    Writing,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct BuildStats {
    pub raw_seed_occurrences: usize,
    pub filtered_hashes: usize,
    pub max_partition_seed_occurrences: usize,
    pub temp_partitions: usize,
    pub seed_chunk_bases: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct RawSeed {
    hash: u64,
    ref_id: u32,
    pos: u32,
}

struct PendingReferenceRecord {
    ref_id: u32,
    record: Record,
}

struct SeededReferenceRecord {
    ref_id: u32,
    reference: Reference,
    raw_seeds: Vec<RawSeed>,
}

#[derive(Clone, Copy, Debug, Default)]
struct PartitionedReferenceStats {
    reference_count: usize,
    total_reference_bases: u64,
    reference_storage_bytes: u64,
    raw_seed_occurrences: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
struct IndexSourceFingerprint {
    size: u64,
    modified_secs: u64,
    modified_nanos: u32,
}

struct SerializedIndexMetadata {
    magic: [u8; 8],
    version: u32,
    metam_version: String,
    k: usize,
    w: usize,
    seed_mode: SeedMode,
    max_occ: usize,
    hash_bits: u8,
    exact_randstrobe_match: bool,
    strobe_start_dist: Option<StrobeStartDistance>,
    created_unix_seconds: u64,
    filtered_seed_hashes: u64,
}

pub enum LoadedIndex {
    Owned(Index),
    Mapped(MappedIndex),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct IndexCacheBuildProgress {
    pub compressed_bytes_read: u64,
    pub compressed_bytes_total: u64,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SeedLookupMode {
    #[default]
    IndexDefault,
    IgnoreRandstrobeFingerprint,
}

impl SeedLookupMode {
    fn use_fingerprint_filter(self) -> bool {
        matches!(self, Self::IndexDefault)
    }
}

#[derive(Clone, Debug)]
struct MappedReferenceMeta {
    id: String,
    len: u64,
    packed: std::ops::Range<usize>,
    ambiguous: std::ops::Range<usize>,
}

fn push_nonempty_range(ranges: &mut Vec<std::ops::Range<usize>>, range: std::ops::Range<usize>) {
    if !range.is_empty() {
        ranges.push(range);
    }
}

pub struct MappedReference<'a> {
    id: &'a str,
    len: u64,
    packed_bases: &'a [u8],
    ambiguous: &'a [u8],
}

pub struct MappedSeedTable {
    hash_bits: u8,
    key_layout: CacheSeedKeyLayout,
    hashes32: std::ops::Range<usize>,
    hashes32_len: usize,
    hashes64: std::ops::Range<usize>,
    hashes64_len: usize,
    group_offset_width: u8,
    group_offsets: std::ops::Range<usize>,
    group_offset_len: usize,
    occurrences: std::ops::Range<usize>,
    occurrence_len: usize,
    fingerprint_bits: u8,
    fingerprints: std::ops::Range<usize>,
    fingerprint_len: usize,
    prefix_bits: u8,
    offsets: std::ops::Range<usize>,
    offset_len: usize,
}

pub struct MappedIndex {
    mmap: MmapFile,
    header: IndexHeader,
    references: Vec<MappedReferenceMeta>,
    seed_table: MappedSeedTable,
}

#[cfg(unix)]
struct MmapFile {
    ptr: *mut libc::c_void,
    len: usize,
}

#[cfg(unix)]
unsafe impl Send for MmapFile {}

#[cfg(unix)]
unsafe impl Sync for MmapFile {}

#[cfg(unix)]
impl MmapFile {
    fn map(path: &str) -> Result<Self> {
        use std::os::fd::AsRawFd;

        let file =
            File::open(path).with_context(|| format!("failed to open index cache {path}"))?;
        let len = file
            .metadata()
            .with_context(|| format!("failed to stat index cache {path}"))?
            .len() as usize;
        if len == 0 {
            bail!("empty index cache");
        }
        let ptr = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                len,
                libc::PROT_READ,
                libc::MAP_PRIVATE,
                file.as_raw_fd(),
                0,
            )
        };
        if ptr == libc::MAP_FAILED {
            return Err(std::io::Error::last_os_error())
                .with_context(|| format!("failed to memory-map index cache {path}"));
        }
        Ok(Self { ptr, len })
    }

    fn advise(&self, range: std::ops::Range<usize>, advice: libc::c_int) {
        let start = range.start.min(self.len);
        let end = range.end.min(self.len);
        if start >= end {
            return;
        }
        let page_size = unsafe { libc::sysconf(libc::_SC_PAGESIZE) };
        let page_size = usize::try_from(page_size).ok().filter(|size| *size > 0);
        let Some(page_size) = page_size else {
            return;
        };
        let aligned_start = start / page_size * page_size;
        let len = end - aligned_start;
        unsafe {
            let address = self.ptr.cast::<u8>().add(aligned_start).cast();
            // Advice is a best-effort performance hint. Mapping correctness is
            // independent of whether the platform accepts it.
            libc::madvise(address, len, advice);
        }
    }

    fn as_slice(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr.cast::<u8>(), self.len) }
    }
}

#[cfg(unix)]
impl Drop for MmapFile {
    fn drop(&mut self) {
        unsafe {
            libc::munmap(self.ptr, self.len);
        }
    }
}

#[cfg(not(unix))]
struct MmapFile;

#[cfg(not(unix))]
impl MmapFile {
    fn map(_path: &str) -> Result<Self> {
        bail!("memory-mapped index cache is only supported on Unix-like systems")
    }

    fn as_slice(&self) -> &[u8] {
        &[]
    }
}

impl LoadedIndex {
    pub fn is_memory_mapped(&self) -> bool {
        matches!(self, Self::Mapped(_))
    }

    pub fn header(&self) -> &IndexHeader {
        match self {
            Self::Owned(index) => &index.header,
            Self::Mapped(index) => &index.header,
        }
    }

    pub fn reference_count(&self) -> usize {
        match self {
            Self::Owned(index) => index.references.len(),
            Self::Mapped(index) => index.references.len(),
        }
    }

    pub fn reference_id(&self, ref_id: u32) -> &str {
        match self {
            Self::Owned(index) => &index.references[ref_id as usize].id,
            Self::Mapped(index) => &index.references[ref_id as usize].id,
        }
    }

    pub fn reference(&self, ref_id: u32) -> ReferenceView<'_> {
        match self {
            Self::Owned(index) => ReferenceView::Owned(&index.references[ref_id as usize]),
            Self::Mapped(index) => ReferenceView::Mapped(index.reference(ref_id)),
        }
    }

    pub fn lookup(&self, hash: u64) -> SeedOccurrences<'_> {
        self.lookup_with_mode(hash, SeedLookupMode::IndexDefault)
    }

    pub fn lookup_with_mode(&self, hash: u64, mode: SeedLookupMode) -> SeedOccurrences<'_> {
        match self {
            Self::Owned(index) => index.lookup_with_mode(hash, mode),
            Self::Mapped(index) => index.lookup_with_mode(hash, mode),
        }
    }
}

pub enum ReferenceView<'a> {
    Owned(&'a Reference),
    Mapped(MappedReference<'a>),
}

impl ReferenceView<'_> {
    pub fn id(&self) -> &str {
        match self {
            Self::Owned(reference) => &reference.id,
            Self::Mapped(reference) => reference.id,
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Self::Owned(reference) => reference.len(),
            Self::Mapped(reference) => reference.len as usize,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn base_at(&self, pos: usize) -> u8 {
        match self {
            Self::Owned(reference) => reference.base_at(pos),
            Self::Mapped(reference) => reference.base_at(pos),
        }
    }

    /// See [`Reference::code_at`].
    #[inline]
    pub fn code_at(&self, pos: usize) -> Option<u8> {
        match self {
            Self::Owned(reference) => reference.code_at(pos),
            Self::Mapped(reference) => reference.code_at(pos),
        }
    }
}

impl MappedReference<'_> {
    #[inline]
    fn code_at(&self, pos: usize) -> Option<u8> {
        if mapped_ambiguous_get(self.ambiguous, pos) {
            None
        } else {
            Some((self.packed_bases[pos / 4] >> (6 - 2 * (pos % 4))) & 3)
        }
    }

    fn base_at(&self, pos: usize) -> u8 {
        match self.code_at(pos) {
            Some(code) => code_base(code),
            None => b'N',
        }
    }
}

fn mapped_ambiguous_get(bits: &[u8], i: usize) -> bool {
    bits.get(i / 8)
        .map(|byte| (byte & (1 << (i % 8))) != 0)
        .unwrap_or(false)
}

#[allow(clippy::too_many_arguments)]
pub fn build_index_with_progress<F>(
    ref_path: &str,
    k: usize,
    w: usize,
    max_occ: usize,
    hash_bits: u8,
    seed_mode: SeedMode,
    exact_randstrobe_match: bool,
    strobe_start_dist: Option<StrobeStartDistance>,
    threads: usize,
    seed_chunk_bases: usize,
    mut progress: F,
) -> Result<(Index, BuildStats)>
where
    F: FnMut(BuildProgress),
{
    validate_seed_params(k, w)?;
    if max_occ == 0 {
        bail!("--max-occ must be greater than 0");
    }
    if !(24..=64).contains(&hash_bits) {
        bail!("--hash-bits must be between 24 and 64");
    }
    validate_exact_randstrobe_match(seed_mode, exact_randstrobe_match)?;
    let strobe_start_dist = normalize_strobe_start_distance(seed_mode, strobe_start_dist)?;
    let seed_chunk_bases = normalize_seed_chunk_bases(seed_chunk_bases)?;

    progress(BuildProgress::Phase("reading reference FASTA"));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("failed to create indexing thread pool")?;

    let partition_count = temporary_seed_partition_count(threads);
    let temp_dir = make_temp_dir_in(None)?;
    let _temp_dir_cleanup = TempDirCleanup::new(temp_dir.clone());
    let partition_paths = (0..partition_count)
        .map(|i| temp_dir.join(format!("part-{i:04}.bin")))
        .collect::<Vec<_>>();
    let mut writers = partition_paths
        .iter()
        .map(|path| File::create(path).map(BufWriter::new))
        .collect::<std::io::Result<Vec<_>>>()
        .context("failed to create temporary seed partitions")?;

    progress(BuildProgress::Phase("partitioning seed occurrences"));
    let ref_store_path = temp_dir.join("references.bin");
    let mut ref_store = BufWriter::new(
        File::create(&ref_store_path)
            .with_context(|| format!("failed to create {}", ref_store_path.display()))?,
    );
    let partitioned = partition_references_parallel(
        ref_path,
        k,
        w,
        seed_mode,
        strobe_start_dist,
        hash_bits,
        partition_count,
        &pool,
        seed_chunk_bases,
        &mut writers,
        &mut ref_store,
        &mut progress,
    )?;
    let reference_count = partitioned.reference_count;
    let total_reference_bases = partitioned.total_reference_bases;
    let reference_storage_bytes = partitioned.reference_storage_bytes;
    let raw_seed_occurrences = partitioned.raw_seed_occurrences;
    drop(writers);
    ref_store.flush()?;
    drop(ref_store);

    progress(BuildProgress::Phase(
        "sorting and filtering seed partitions",
    ));
    let mut seed_table =
        SeedTable::new_with_exact_randstrobe_match(hash_bits, exact_randstrobe_match);
    let mut indexed_seed_hashes = 0usize;
    let mut filtered_hashes = 0usize;
    let mut max_partition_seed_occurrences = 0usize;
    for path in &partition_paths {
        let mut raw = read_partition(path)?;
        if raw.is_empty() {
            continue;
        }
        max_partition_seed_occurrences = max_partition_seed_occurrences.max(raw.len());
        pool.install(|| sort_raw_seeds(&mut raw, hash_bits));
        let mut i = 0usize;
        while i < raw.len() {
            let hash = seed_key(raw[i].hash, hash_bits);
            let mut j = i + 1;
            while j < raw.len() && seed_key(raw[j].hash, hash_bits) == hash {
                j += 1;
            }
            let count = j - i;
            if count <= max_occ {
                indexed_seed_hashes += 1;
                seed_table.push_group(hash, &raw[i..j])?;
            } else {
                filtered_hashes += 1;
            }
            i = j;
        }
    }
    seed_table.rebuild_lookup_index();

    progress(BuildProgress::Phase("loading compact reference table"));
    let references = read_references(&ref_store_path, reference_count)?;
    std::fs::remove_dir_all(&temp_dir).with_context(|| {
        format!(
            "failed to remove completed index workspace {}",
            temp_dir.display()
        )
    })?;

    let header = IndexHeader {
        magic: INDEX_MAGIC,
        version: INDEX_VERSION,
        metam_version: env!("CARGO_PKG_VERSION").to_string(),
        k,
        w,
        seed_mode,
        max_occ,
        hash_bits,
        exact_randstrobe_match,
        strobe_start_dist,
        reference_count,
        total_reference_bases,
        indexed_seed_hashes: indexed_seed_hashes as u64,
        indexed_seed_occurrences: seed_table.occurrence_count() as u64,
        filtered_seed_hashes: filtered_hashes as u64,
        reference_storage_bytes,
        created_unix_seconds: SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0),
    };
    let stats = BuildStats {
        raw_seed_occurrences,
        filtered_hashes,
        max_partition_seed_occurrences,
        temp_partitions: partition_count,
        seed_chunk_bases,
    };
    Ok((
        Index {
            header,
            references,
            seed_table,
        },
        stats,
    ))
}

#[allow(clippy::too_many_arguments)]
pub fn build_index_file_with_progress<F>(
    ref_path: &str,
    out_path: &str,
    k: usize,
    w: usize,
    max_occ: usize,
    hash_bits: u8,
    seed_mode: SeedMode,
    exact_randstrobe_match: bool,
    strobe_start_dist: Option<StrobeStartDistance>,
    threads: usize,
    seed_chunk_bases: usize,
    temp_parent: Option<&Path>,
    mut progress: F,
) -> Result<(IndexHeader, BuildStats)>
where
    F: FnMut(BuildProgress),
{
    validate_seed_params(k, w)?;
    if max_occ == 0 {
        bail!("--max-occ must be greater than 0");
    }
    if !(24..=64).contains(&hash_bits) {
        bail!("--hash-bits must be between 24 and 64");
    }
    validate_exact_randstrobe_match(seed_mode, exact_randstrobe_match)?;
    let strobe_start_dist = normalize_strobe_start_distance(seed_mode, strobe_start_dist)?;
    let seed_chunk_bases = normalize_seed_chunk_bases(seed_chunk_bases)?;

    progress(BuildProgress::Phase("reading reference FASTA"));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("failed to create indexing thread pool")?;

    let partition_count = temporary_seed_partition_count(threads);
    let temp_dir = make_temp_dir_in(temp_parent)?;
    let _temp_dir_cleanup = TempDirCleanup::new(temp_dir.clone());
    let partition_paths = (0..partition_count)
        .map(|i| temp_dir.join(format!("part-{i:04}.bin")))
        .collect::<Vec<_>>();
    let mut writers = partition_paths
        .iter()
        .map(|path| File::create(path).map(BufWriter::new))
        .collect::<std::io::Result<Vec<_>>>()
        .context("failed to create temporary seed partitions")?;

    progress(BuildProgress::Phase("partitioning seed occurrences"));
    let ref_store_path = temp_dir.join("references.bin");
    let mut ref_store = BufWriter::new(
        File::create(&ref_store_path)
            .with_context(|| format!("failed to create {}", ref_store_path.display()))?,
    );
    let partitioned = partition_references_parallel(
        ref_path,
        k,
        w,
        seed_mode,
        strobe_start_dist,
        hash_bits,
        partition_count,
        &pool,
        seed_chunk_bases,
        &mut writers,
        &mut ref_store,
        &mut progress,
    )?;
    let reference_count = partitioned.reference_count;
    let total_reference_bases = partitioned.total_reference_bases;
    let reference_storage_bytes = partitioned.reference_storage_bytes;
    let raw_seed_occurrences = partitioned.raw_seed_occurrences;
    drop(writers);
    ref_store.flush()?;
    drop(ref_store);

    progress(BuildProgress::Phase(
        "sorting and filtering seed partitions",
    ));
    let seed_store_path = temp_dir.join("seed-groups.bin");
    let mut seed_store = BufWriter::new(
        File::create(&seed_store_path)
            .with_context(|| format!("failed to create {}", seed_store_path.display()))?,
    );
    let fingerprint_store_path = temp_dir.join("seed-fingerprints.bin");
    let mut fingerprint_store = exact_randstrobe_match
        .then(|| File::create(&fingerprint_store_path).map(BufWriter::new))
        .transpose()
        .with_context(|| format!("failed to create {}", fingerprint_store_path.display()))?;
    let mut indexed_seed_hashes = 0usize;
    let mut indexed_seed_occurrences = 0usize;
    let mut filtered_hashes = 0usize;
    let mut max_partition_seed_occurrences = 0usize;
    let mut previous_hash = 0u64;
    let mut fingerprint_scratch = Vec::new();
    for path in &partition_paths {
        let mut raw = read_partition(path)?;
        if raw.is_empty() {
            continue;
        }
        max_partition_seed_occurrences = max_partition_seed_occurrences.max(raw.len());
        pool.install(|| sort_raw_seeds(&mut raw, hash_bits));
        let mut i = 0usize;
        while i < raw.len() {
            let hash = seed_key(raw[i].hash, hash_bits);
            let mut j = i + 1;
            while j < raw.len() && seed_key(raw[j].hash, hash_bits) == hash {
                j += 1;
            }
            let count = j - i;
            if count <= max_occ {
                indexed_seed_hashes += 1;
                indexed_seed_occurrences += count;
                write_seed_group_body(&mut seed_store, &mut previous_hash, hash, &raw[i..j])?;
                if let Some(writer) = fingerprint_store.as_mut() {
                    write_seed_fingerprints(
                        writer,
                        &raw[i..j],
                        hash_bits,
                        &mut fingerprint_scratch,
                    )?;
                }
            } else {
                filtered_hashes += 1;
            }
            i = j;
        }
    }
    seed_store.flush()?;
    drop(seed_store);
    if let Some(mut writer) = fingerprint_store {
        writer.flush()?;
    }

    let header = IndexHeader {
        magic: INDEX_MAGIC,
        version: INDEX_VERSION,
        metam_version: env!("CARGO_PKG_VERSION").to_string(),
        k,
        w,
        seed_mode,
        max_occ,
        hash_bits,
        exact_randstrobe_match,
        strobe_start_dist,
        reference_count,
        total_reference_bases,
        indexed_seed_hashes: indexed_seed_hashes as u64,
        indexed_seed_occurrences: indexed_seed_occurrences as u64,
        filtered_seed_hashes: filtered_hashes as u64,
        reference_storage_bytes,
        created_unix_seconds: SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0),
    };

    progress(BuildProgress::Writing);
    write_streamed_index_payload(
        &header,
        &ref_store_path,
        &seed_store_path,
        exact_randstrobe_match.then_some(fingerprint_store_path.as_path()),
        out_path,
    )?;
    std::fs::remove_dir_all(&temp_dir).with_context(|| {
        format!(
            "failed to remove completed index workspace {}",
            temp_dir.display()
        )
    })?;

    Ok((
        header,
        BuildStats {
            raw_seed_occurrences,
            filtered_hashes,
            max_partition_seed_occurrences,
            temp_partitions: partition_count,
            seed_chunk_bases,
        },
    ))
}

pub fn write_index(index: &Index, path: &str) -> Result<()> {
    write_compressed_index_atomic(Path::new(path), |encoder| {
        write_index_payload(encoder, index)
    })
}

fn make_temp_dir_in(parent: Option<&Path>) -> Result<PathBuf> {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let default_parent;
    let parent = match parent {
        Some(parent) => parent,
        None => {
            default_parent = std::env::temp_dir();
            default_parent.as_path()
        }
    };
    if !parent.is_dir() {
        bail!(
            "index temporary directory '{}' does not exist or is not a directory",
            parent.display()
        );
    }
    for attempt in 0..100u32 {
        let path = parent.join(format!(".metam-index-{}-{nanos}-{attempt}", process::id()));
        match std::fs::create_dir(&path) {
            Ok(()) => return Ok(path),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(error).with_context(|| format!("failed to create {}", path.display()));
            }
        }
    }
    bail!(
        "failed to create a unique index temporary directory in {}",
        parent.display()
    )
}

fn temporary_seed_partition_count(threads: usize) -> usize {
    temporary_seed_partition_count_for_limit(threads, process_open_file_limit())
}

fn temporary_seed_partition_count_for_limit(threads: usize, open_file_limit: usize) -> usize {
    let requested = threads
        .max(1)
        .saturating_mul(64)
        .checked_next_power_of_two()
        .unwrap_or(4_096)
        .clamp(256, 4_096);
    let reserved = OPEN_FILE_RESERVE.min(open_file_limit / 4);
    let writer_budget = open_file_limit
        .saturating_sub(reserved)
        .clamp(1, MAX_OPEN_PARTITION_WRITERS);
    requested.min(writer_budget)
}

#[cfg(unix)]
fn process_open_file_limit() -> usize {
    let mut limit = std::mem::MaybeUninit::<libc::rlimit>::uninit();
    let result = unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, limit.as_mut_ptr()) };
    if result != 0 {
        return FALLBACK_OPEN_FILE_LIMIT;
    }
    let current = unsafe { limit.assume_init() }.rlim_cur;
    usize::try_from(current).unwrap_or(usize::MAX)
}

#[cfg(not(unix))]
fn process_open_file_limit() -> usize {
    FALLBACK_OPEN_FILE_LIMIT
}

struct TempDirCleanup {
    path: PathBuf,
}

struct TempFileCleanup {
    path: PathBuf,
    armed: bool,
}

impl TempFileCleanup {
    fn new(path: PathBuf) -> Self {
        Self { path, armed: true }
    }

    fn disarm(&mut self) {
        self.armed = false;
    }
}

impl Drop for TempFileCleanup {
    fn drop(&mut self) {
        if self.armed {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

impl TempDirCleanup {
    fn new(path: PathBuf) -> Self {
        Self { path }
    }
}

impl Drop for TempDirCleanup {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

fn write_raw_seed<W: Write>(writer: &mut W, raw: RawSeed) -> Result<()> {
    writer.write_all(&raw.hash.to_le_bytes())?;
    writer.write_all(&raw.ref_id.to_le_bytes())?;
    writer.write_all(&raw.pos.to_le_bytes())?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn partition_references_parallel<F>(
    ref_path: &str,
    k: usize,
    w: usize,
    seed_mode: SeedMode,
    strobe_start_dist: Option<StrobeStartDistance>,
    hash_bits: u8,
    partition_count: usize,
    pool: &rayon::ThreadPool,
    seed_chunk_bases: usize,
    writers: &mut [BufWriter<File>],
    ref_store: &mut BufWriter<File>,
    progress: &mut F,
) -> Result<PartitionedReferenceStats>
where
    F: FnMut(BuildProgress),
{
    let mut reader = FastaReader::new(open_reader(ref_path)?);
    let mut seen_ids = HashSet::new();
    let mut stats = PartitionedReferenceStats::default();

    loop {
        let batch = read_reference_seed_batch(
            &mut reader,
            &mut seen_ids,
            stats.reference_count,
            seed_chunk_bases,
        )?;
        if batch.is_empty() {
            break;
        }

        let mut seeded = pool.install(|| {
            batch
                .into_par_iter()
                .map(|pending| seed_reference_record(pending, k, w, seed_mode, strobe_start_dist))
                .collect::<Result<Vec<_>>>()
        })?;

        seeded.sort_by_key(|record| record.ref_id);
        for record in seeded {
            stats.raw_seed_occurrences += record.raw_seeds.len();
            stats.total_reference_bases += record.reference.len;
            stats.reference_storage_bytes += record.reference.packed_bases.len() as u64
                + record.reference.ambiguous.byte_len() as u64;

            for raw in record.raw_seeds {
                let key = seed_key(raw.hash, hash_bits);
                let partition = partition_for_hash(key, hash_bits, partition_count);
                write_raw_seed(&mut writers[partition], raw)?;
            }
            write_reference(ref_store, &record.reference)?;

            stats.reference_count += 1;
            progress(BuildProgress::ReferencesSeeded(stats.reference_count, 0));
        }
    }

    if stats.reference_count == 0 {
        bail!("reference FASTA contains no records");
    }
    progress(BuildProgress::ReferencesRead(stats.reference_count));
    Ok(stats)
}

fn read_reference_seed_batch(
    reader: &mut FastaReader,
    seen_ids: &mut HashSet<String>,
    reference_count: usize,
    target_bases: usize,
) -> Result<Vec<PendingReferenceRecord>> {
    let mut batch = Vec::new();
    let mut batch_bases = 0usize;
    let mut next_ref_id = reference_count;

    while batch.is_empty() || batch_bases < target_bases {
        let Some(record) = reader.next_record()? else {
            break;
        };
        if !seen_ids.insert(record.id.clone()) {
            bail!(
                "reference FASTA contains duplicate sequence ID {}",
                record.id
            );
        }
        if record.seq.len() > u32::MAX as usize {
            bail!("reference {} exceeds 2^32-1 bases", record.id);
        }
        if next_ref_id > u32::MAX as usize {
            bail!("too many reference records for current index format");
        }
        batch_bases = batch_bases.saturating_add(record.seq.len());
        batch.push(PendingReferenceRecord {
            ref_id: next_ref_id as u32,
            record,
        });
        next_ref_id += 1;
    }

    Ok(batch)
}

fn seed_reference_record(
    pending: PendingReferenceRecord,
    k: usize,
    w: usize,
    seed_mode: SeedMode,
    strobe_start_dist: Option<StrobeStartDistance>,
) -> Result<SeededReferenceRecord> {
    let PendingReferenceRecord { ref_id, record } = pending;
    let mut raw_seeds = Vec::with_capacity(estimated_seed_capacity(record.seq.len(), w));
    emit_record_seeds(&record.seq, k, w, seed_mode, strobe_start_dist, |seed| {
        raw_seeds.push(RawSeed {
            hash: seed.hash,
            ref_id,
            pos: seed.pos,
        });
        Ok(())
    })?;
    Ok(SeededReferenceRecord {
        ref_id,
        reference: Reference::from_record(&record),
        raw_seeds,
    })
}

fn estimated_seed_capacity(seq_len: usize, w: usize) -> usize {
    const MAX_PREALLOCATED_SEEDS_PER_RECORD: usize = 4_000_000;
    (seq_len / w.max(1)).min(MAX_PREALLOCATED_SEEDS_PER_RECORD)
}

fn normalize_seed_chunk_bases(seed_chunk_bases: usize) -> Result<usize> {
    if seed_chunk_bases == 0 {
        bail!("index seed chunk size must be greater than 0");
    }
    Ok(seed_chunk_bases)
}

fn emit_record_seeds<F>(
    seq: &[u8],
    k: usize,
    w: usize,
    mode: SeedMode,
    strobe_start_dist: Option<StrobeStartDistance>,
    mut emit: F,
) -> Result<()>
where
    F: FnMut(Seed) -> Result<()>,
{
    if seq.len() <= INDEX_SEQUENCE_CHUNK {
        for seed in seeds_with_strobe_start_distance(seq, k, w, mode, strobe_start_dist) {
            emit(seed)?;
        }
        return Ok(());
    }

    let max_seed_span = strobe_start_dist
        .map(|distance| usize::from(distance.max).saturating_add(k))
        .unwrap_or_else(|| k.saturating_add(263));
    let overlap = k
        .saturating_add(263)
        .max(w.saturating_add(k))
        .max(max_seed_span);
    let mut start = 0usize;
    while start < seq.len() {
        let core_end = start.saturating_add(INDEX_SEQUENCE_CHUNK).min(seq.len());
        let seq_end = core_end.saturating_add(overlap).min(seq.len());
        for mut seed in
            seeds_with_strobe_start_distance(&seq[start..seq_end], k, w, mode, strobe_start_dist)
        {
            let global_pos = start + seed.pos as usize;
            if global_pos < core_end {
                seed.pos = global_pos as u32;
                emit(seed)?;
            }
        }
        start = core_end;
    }
    Ok(())
}

fn read_partition(path: &Path) -> Result<Vec<RawSeed>> {
    let inner = File::open(path)?;
    let bytes = inner.metadata()?.len() as usize;
    let mut file = BufReader::new(inner);
    if bytes == 0 {
        return Ok(Vec::new());
    }
    if !bytes.is_multiple_of(16) {
        bail!(
            "temporary seed partition has invalid length: {}",
            path.display()
        );
    }
    let mut out = Vec::with_capacity(bytes / 16);
    for _ in 0..bytes / 16 {
        let mut chunk = [0u8; 16];
        file.read_exact(&mut chunk)?;
        let hash = u64::from_le_bytes(chunk[0..8].try_into().expect("hash bytes"));
        let ref_id = u32::from_le_bytes(chunk[8..12].try_into().expect("ref id bytes"));
        let pos = u32::from_le_bytes(chunk[12..16].try_into().expect("pos bytes"));
        out.push(RawSeed { hash, ref_id, pos });
    }
    Ok(out)
}

pub fn read_index(path: &str) -> Result<Index> {
    let file = File::open(path).with_context(|| format!("failed to open index {path}"))?;
    let reader = BufReader::new(file);
    let mut decoder = zstd::stream::Decoder::new(reader)
        .with_context(|| format!("invalid or unsupported index compression in {path}"))?;
    let index =
        read_index_payload(&mut decoder).with_context(|| format!("invalid index {path}"))?;
    validate_index(&index)?;
    Ok(index)
}

pub fn read_index_with_cache(path: &str, cache_path: Option<&str>) -> Result<LoadedIndex> {
    read_index_with_cache_progress(path, cache_path, |_| {})
}

pub fn read_index_with_cache_progress<F>(
    path: &str,
    cache_path: Option<&str>,
    mut on_progress: F,
) -> Result<LoadedIndex>
where
    F: FnMut(IndexCacheBuildProgress),
{
    let Some(cache_path) = cache_path else {
        return read_index(path).map(LoadedIndex::Owned);
    };
    let source = index_source_fingerprint(path)?;
    if let Ok(index) = MappedIndex::open(cache_path, source) {
        return Ok(LoadedIndex::Mapped(index));
    }

    // A cache that cannot be opened by this version is a disposable derived
    // artifact. Remove it before rebuilding so a large obsolete cache and its
    // replacement do not have to coexist on disk. Existing mmaps remain valid
    // on Unix after unlinking the directory entry.
    let _ = fs::remove_file(cache_path);

    // Build the disposable mmap cache directly from the compressed index.
    // The streaming transcode holds only the compact key metadata in memory
    // and writes occurrences straight to their final cache range, rather than
    // materializing every reference, seed key, occurrence and fingerprint at
    // once as a non-streaming build would.
    if write_index_cache_from_source_progress(cache_path, source, path, &mut on_progress).is_ok() {
        if let Ok(mapped) = MappedIndex::open(cache_path, source) {
            return Ok(LoadedIndex::Mapped(mapped));
        }
    }
    read_index(path).map(LoadedIndex::Owned)
}

fn index_source_fingerprint(path: &str) -> Result<IndexSourceFingerprint> {
    let metadata = fs::metadata(path).with_context(|| format!("failed to stat index {path}"))?;
    let modified = metadata.modified().unwrap_or(SystemTime::UNIX_EPOCH);
    let duration = modified
        .duration_since(UNIX_EPOCH)
        .unwrap_or_else(|_| std::time::Duration::from_secs(0));
    Ok(IndexSourceFingerprint {
        size: metadata.len(),
        modified_secs: duration.as_secs(),
        modified_nanos: duration.subsec_nanos(),
    })
}

fn write_index_cache_from_source_progress<F>(
    cache_path: &str,
    source: IndexSourceFingerprint,
    index_path: &str,
    mut on_progress: F,
) -> Result<()>
where
    F: FnMut(IndexCacheBuildProgress),
{
    let cache_path = Path::new(cache_path);
    if let Some(parent) = cache_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).with_context(|| {
            format!(
                "failed to create index cache directory {}",
                parent.display()
            )
        })?;
    }

    let tmp_path = cache_tmp_path(cache_path);
    let _tmp_cleanup = CacheTempFileCleanup(tmp_path.clone());
    {
        let source_file =
            File::open(index_path).with_context(|| format!("failed to open index {index_path}"))?;
        on_progress(IndexCacheBuildProgress {
            compressed_bytes_read: 0,
            compressed_bytes_total: source.size,
        });
        let source_reader = IndexCacheProgressReader {
            inner: source_file,
            bytes_read: 0,
            bytes_total: source.size,
            on_progress: &mut on_progress,
        };
        let source_reader = BufReader::with_capacity(INDEX_CACHE_IO_BUFFER_BYTES, source_reader);
        let mut decoder = zstd::stream::Decoder::new(source_reader)
            .with_context(|| format!("invalid or unsupported index compression in {index_path}"))?;
        let cache_file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .read(true)
            .write(true)
            .open(&tmp_path)
            .with_context(|| format!("failed to create index cache {}", tmp_path.display()))?;
        let mut writer = BufWriter::with_capacity(INDEX_CACHE_IO_BUFFER_BYTES, cache_file);
        transcode_index_cache(&mut decoder, &mut writer, &tmp_path, source)
            .with_context(|| format!("failed to build index cache {}", tmp_path.display()))?;
        writer.flush()?;
    }
    on_progress(IndexCacheBuildProgress {
        compressed_bytes_read: source.size,
        compressed_bytes_total: source.size,
    });
    fs::rename(&tmp_path, cache_path).with_context(|| {
        format!(
            "failed to move index cache {} to {}",
            tmp_path.display(),
            cache_path.display()
        )
    })?;
    Ok(())
}

struct IndexCacheProgressReader<R, F> {
    inner: R,
    bytes_read: u64,
    bytes_total: u64,
    on_progress: F,
}

impl<R, F> Read for IndexCacheProgressReader<R, F>
where
    R: Read,
    F: FnMut(IndexCacheBuildProgress),
{
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.bytes_read = self
            .bytes_read
            .saturating_add(read as u64)
            .min(self.bytes_total);
        (self.on_progress)(IndexCacheBuildProgress {
            compressed_bytes_read: self.bytes_read,
            compressed_bytes_total: self.bytes_total,
        });
        Ok(read)
    }
}

#[derive(Clone, Debug)]
struct StreamingCacheRanges {
    lookup_offsets: std::ops::Range<u64>,
    seed_keys: std::ops::Range<u64>,
    group_offsets: std::ops::Range<u64>,
    occurrences: std::ops::Range<u64>,
    fingerprints: std::ops::Range<u64>,
    total_reference_bases_patch: u64,
    indexed_seed_hashes_patch: u64,
    indexed_seed_occurrences_patch: u64,
    reference_storage_bytes_patch: u64,
    reference_metadata_offset_patch: u64,
    reference_metadata_len_patch: u64,
    reference_metadata_offset: u64,
    reference_metadata_len: u64,
    final_len: u64,
}

struct CacheReferenceMetadata {
    id: String,
    len: u64,
    packed_len: u64,
    ambiguous_len: u64,
}

fn transcode_index_cache<R: Read>(
    reader: &mut R,
    writer: &mut BufWriter<File>,
    cache_tmp_path: &Path,
    source: IndexSourceFingerprint,
) -> Result<()> {
    if !cfg!(target_endian = "little") {
        bail!("memory-mapped index cache requires a little-endian target");
    }

    let metadata = read_serialized_index_metadata(reader)?;
    let reference_count_u64 = read_varint(reader)?;
    let reference_count = usize::try_from(reference_count_u64)
        .context("index reference count exceeds address space")?;
    let mut pos = 0u64;

    write_raw(writer, &mut pos, &INDEX_CACHE_MAGIC)?;
    write_u32_le(writer, &mut pos, INDEX_CACHE_VERSION)?;
    write_u64_le(writer, &mut pos, source.size)?;
    write_u64_le(writer, &mut pos, source.modified_secs)?;
    write_u32_le(writer, &mut pos, source.modified_nanos)?;
    write_raw(writer, &mut pos, &metadata.magic)?;
    write_u32_le(writer, &mut pos, metadata.version)?;
    write_cache_string(writer, &mut pos, &metadata.metam_version)?;
    write_u64_le(writer, &mut pos, metadata.k as u64)?;
    write_u64_le(writer, &mut pos, metadata.w as u64)?;
    write_u8(writer, &mut pos, seed_mode_code(metadata.seed_mode) as u8)?;
    write_u64_le(writer, &mut pos, metadata.max_occ as u64)?;
    write_u8(writer, &mut pos, metadata.hash_bits)?;
    write_u8(writer, &mut pos, u8::from(metadata.exact_randstrobe_match))?;
    let (strobe_start_dist_min, strobe_start_dist_max) = metadata
        .strobe_start_dist
        .map(|distance| (distance.min, distance.max))
        .unwrap_or((0, 0));
    write_u64_le(writer, &mut pos, u64::from(strobe_start_dist_min))?;
    write_u64_le(writer, &mut pos, u64::from(strobe_start_dist_max))?;
    write_u64_le(writer, &mut pos, reference_count_u64)?;
    let total_reference_bases_patch = pos;
    write_u64_le(writer, &mut pos, 0)?;
    let indexed_seed_hashes_patch = pos;
    write_u64_le(writer, &mut pos, 0)?;
    let indexed_seed_occurrences_patch = pos;
    write_u64_le(writer, &mut pos, 0)?;
    write_u64_le(writer, &mut pos, metadata.filtered_seed_hashes)?;
    let reference_storage_bytes_patch = pos;
    write_u64_le(writer, &mut pos, 0)?;
    write_u64_le(writer, &mut pos, metadata.created_unix_seconds)?;

    write_u64_le(writer, &mut pos, reference_count_u64)?;
    let reference_metadata_offset_patch = pos;
    write_u64_le(writer, &mut pos, 0)?;
    let reference_metadata_len_patch = pos;
    write_u64_le(writer, &mut pos, 0)?;
    let mut total_reference_bases = 0u64;
    let mut reference_storage_bytes = 0u64;
    let mut reference_metadata = Vec::with_capacity(reference_count);
    let mut copy_buffer = vec![0u8; INDEX_CACHE_IO_BUFFER_BYTES];
    for _ in 0..reference_count {
        let id = read_string(reader)?;
        let len = read_varint(reader)?;
        let packed_len = read_varint(reader)?;
        let expected_packed_len = len.div_ceil(4);
        if packed_len != expected_packed_len {
            bail!("index reference has an invalid packed sequence length");
        }

        copy_exact_cache_bytes(reader, writer, &mut pos, packed_len, &mut copy_buffer)?;

        let ambiguous_len = read_varint(reader)?;
        let expected_ambiguous_len = len.div_ceil(8);
        if ambiguous_len != 0 && ambiguous_len != expected_ambiguous_len {
            bail!("index reference has an invalid ambiguity-mask length");
        }
        copy_exact_cache_bytes(reader, writer, &mut pos, ambiguous_len, &mut copy_buffer)?;

        reference_metadata.push(CacheReferenceMetadata {
            id,
            len,
            packed_len,
            ambiguous_len,
        });

        total_reference_bases = total_reference_bases
            .checked_add(len)
            .context("index total reference length overflow")?;
        reference_storage_bytes = reference_storage_bytes
            .checked_add(packed_len)
            .and_then(|value| value.checked_add(ambiguous_len))
            .context("index reference storage length overflow")?;
    }

    // Keep all IDs and fixed-size reference metadata together. Opening the
    // cache can now read this compact section sequentially instead of faulting
    // one page near each packed reference sequence on cold/network storage.
    let reference_metadata_offset = pos;
    for reference in &reference_metadata {
        write_cache_string(writer, &mut pos, &reference.id)?;
        write_u64_le(writer, &mut pos, reference.len)?;
        write_u64_le(writer, &mut pos, reference.packed_len)?;
        write_u64_le(writer, &mut pos, reference.ambiguous_len)?;
    }
    let reference_metadata_len = pos
        .checked_sub(reference_metadata_offset)
        .context("index cache reference-metadata length underflow")?;
    drop(reference_metadata);

    let occurrence_count_u64 = read_varint(reader)?;
    let group_count_u64 = read_varint(reader)?;
    let occurrence_count = usize::try_from(occurrence_count_u64)
        .context("index seed occurrence count exceeds address space")?;
    let group_count =
        usize::try_from(group_count_u64).context("index seed group count exceeds address space")?;
    if group_count > occurrence_count {
        bail!("index seed group count exceeds occurrence count");
    }

    let key_layout =
        choose_cache_seed_key_layout(metadata.hash_bits, occurrence_count, group_count);
    let key_count = match key_layout {
        CacheSeedKeyLayout::Repeated => occurrence_count,
        CacheSeedKeyLayout::Grouped => group_count,
    };
    let prefix_bits = if key_count == 0 {
        0
    } else {
        metadata.hash_bits.min(20)
    };
    let lookup_offset_len = if key_count == 0 {
        0
    } else {
        (1usize << prefix_bits) + 1
    };
    let key_width = if metadata.hash_bits <= 32 { 4 } else { 8 };
    let group_offset_width = if key_layout == CacheSeedKeyLayout::Grouped {
        if u32::try_from(occurrence_count).is_ok() {
            4
        } else {
            8
        }
    } else {
        0
    };
    let group_offset_len = if key_layout == CacheSeedKeyLayout::Grouped {
        group_count
            .checked_add(1)
            .context("index cache seed-group count overflow")?
    } else {
        0
    };
    let fingerprint_bits = if metadata.exact_randstrobe_match {
        randstrobe_fingerprint_bits(metadata.hash_bits)
    } else {
        0
    };
    let fingerprint_len = occurrence_count
        .checked_mul(fingerprint_bytes(fingerprint_bits))
        .context("index fingerprint array length overflow")?;

    write_u8(writer, &mut pos, metadata.hash_bits)?;
    write_u8(writer, &mut pos, key_layout as u8)?;
    write_u8(writer, &mut pos, prefix_bits)?;
    write_u64_le(writer, &mut pos, lookup_offset_len as u64)?;
    align_writer(writer, &mut pos, std::mem::align_of::<u64>() as u64)?;
    let lookup_offsets = reserve_cache_range(
        writer,
        &mut pos,
        cache_array_bytes(lookup_offset_len, std::mem::size_of::<u64>())?,
    )?;

    let seed_keys = if metadata.hash_bits <= 32 {
        write_u64_le(writer, &mut pos, key_count as u64)?;
        align_writer(writer, &mut pos, std::mem::align_of::<u32>() as u64)?;
        let range = reserve_cache_range(
            writer,
            &mut pos,
            cache_array_bytes(key_count, std::mem::size_of::<u32>())?,
        )?;
        write_u64_le(writer, &mut pos, 0)?;
        align_writer(writer, &mut pos, std::mem::align_of::<u64>() as u64)?;
        range
    } else {
        write_u64_le(writer, &mut pos, 0)?;
        align_writer(writer, &mut pos, std::mem::align_of::<u32>() as u64)?;
        write_u64_le(writer, &mut pos, key_count as u64)?;
        align_writer(writer, &mut pos, std::mem::align_of::<u64>() as u64)?;
        reserve_cache_range(
            writer,
            &mut pos,
            cache_array_bytes(key_count, std::mem::size_of::<u64>())?,
        )?
    };

    write_u8(writer, &mut pos, group_offset_width as u8)?;
    write_u64_le(writer, &mut pos, group_offset_len as u64)?;
    let group_offsets = if group_offset_width == 0 {
        pos..pos
    } else {
        align_writer(writer, &mut pos, group_offset_width as u64)?;
        reserve_cache_range(
            writer,
            &mut pos,
            cache_array_bytes(group_offset_len, group_offset_width)?,
        )?
    };

    write_u64_le(writer, &mut pos, occurrence_count_u64)?;
    align_writer(writer, &mut pos, std::mem::align_of::<u64>() as u64)?;
    let occurrences = reserve_cache_range(
        writer,
        &mut pos,
        cache_array_bytes(occurrence_count, std::mem::size_of::<u64>())?,
    )?;
    write_u8(writer, &mut pos, fingerprint_bits)?;
    write_u64_le(writer, &mut pos, fingerprint_len as u64)?;
    let fingerprints = reserve_cache_range(writer, &mut pos, fingerprint_len as u64)?;
    let ranges = StreamingCacheRanges {
        lookup_offsets,
        seed_keys,
        group_offsets,
        occurrences,
        fingerprints,
        total_reference_bases_patch,
        indexed_seed_hashes_patch,
        indexed_seed_occurrences_patch,
        reference_storage_bytes_patch,
        reference_metadata_offset_patch,
        reference_metadata_len_patch,
        reference_metadata_offset,
        reference_metadata_len,
        final_len: pos,
    };
    writer.flush()?;
    writer.get_ref().set_len(ranges.final_len)?;

    let mut key_file = OpenOptions::new().write(true).open(cache_tmp_path)?;
    key_file.seek(SeekFrom::Start(ranges.seed_keys.start))?;
    let mut key_writer = BufWriter::with_capacity(INDEX_CACHE_IO_BUFFER_BYTES, key_file);
    let mut occurrence_file = OpenOptions::new().write(true).open(cache_tmp_path)?;
    occurrence_file.seek(SeekFrom::Start(ranges.occurrences.start))?;
    let mut occurrence_writer =
        BufWriter::with_capacity(INDEX_CACHE_IO_BUFFER_BYTES, occurrence_file);

    let mut lookup_values = vec![0u64; lookup_offset_len];
    let mut next_lookup_bucket = usize::from(key_count > 0);
    let mut compact_group_offsets = if key_layout == CacheSeedKeyLayout::Grouped {
        let mut offsets = CompactGroupOffsets::with_capacity(occurrence_count, group_count)?;
        offsets.push(0)?;
        Some(offsets)
    } else {
        None
    };
    let mut key_buffer = Vec::with_capacity(INDEX_CACHE_IO_BUFFER_BYTES);
    let mut occurrence_buffer = Vec::with_capacity(INDEX_CACHE_IO_BUFFER_BYTES);
    let mut previous_hash = 0u64;
    let mut decoded_occurrences = 0usize;
    let mut key_bytes_written = 0u64;
    let mut occurrence_bytes_written = 0u64;

    for group in 0..group_count {
        let hash_delta = read_varint(reader)?;
        let hash = previous_hash
            .checked_add(hash_delta)
            .context("seed hash delta overflow in index")?;
        if group > 0 && hash <= previous_hash {
            bail!("index seed groups are not strictly ordered by hash");
        }
        if metadata.hash_bits < 64 && hash >= (1u64 << metadata.hash_bits) {
            bail!("index seed hash exceeds configured hash width");
        }
        previous_hash = hash;

        let group_occurrences_u64 = read_varint(reader)?;
        let group_occurrences = usize::try_from(group_occurrences_u64)
            .context("seed occurrence group exceeds address space")?;
        if group_occurrences == 0 || group_occurrences > u32::MAX as usize {
            bail!("seed occurrence group exceeds index format limit");
        }
        let group_start = decoded_occurrences;
        decoded_occurrences = decoded_occurrences
            .checked_add(group_occurrences)
            .context("index seed occurrence count overflow")?;
        if decoded_occurrences > occurrence_count {
            bail!("index seed groups exceed declared occurrence count");
        }

        let key_row = match key_layout {
            CacheSeedKeyLayout::Repeated => group_start,
            CacheSeedKeyLayout::Grouped => group,
        };
        update_streaming_lookup_offsets(
            &mut lookup_values,
            &mut next_lookup_bucket,
            metadata.hash_bits,
            prefix_bits,
            hash,
            key_row,
        )?;
        let key_repetitions = match key_layout {
            CacheSeedKeyLayout::Repeated => group_occurrences,
            CacheSeedKeyLayout::Grouped => 1,
        };
        // Production indexes contain tens of millions of groups. Keep this
        // tiny value on the stack instead of allocating once per group.
        let key_bytes = hash.to_le_bytes();
        append_repeated_cache_bytes(
            &mut key_writer,
            &mut key_buffer,
            &key_bytes[..key_width],
            key_repetitions,
        )?;
        key_bytes_written = key_bytes_written
            .checked_add(cache_array_bytes(key_repetitions, key_width)?)
            .context("index cache key length overflow")?;

        let mut previous_ref = 0u32;
        let mut previous_pos = 0u32;
        for _ in 0..group_occurrences {
            let ref_delta = read_varint(reader)?;
            let ref_id = u64::from(previous_ref)
                .checked_add(ref_delta)
                .context("reference ID delta overflow in index")?;
            if ref_id > u64::from(u32::MAX) {
                bail!("reference ID exceeds index format limit");
            }
            let pos_delta_or_abs = read_varint(reader)?;
            let reference_pos = if ref_delta == 0 {
                u64::from(previous_pos)
                    .checked_add(pos_delta_or_abs)
                    .context("reference position delta overflow in index")?
            } else {
                pos_delta_or_abs
            };
            if reference_pos > u64::from(u32::MAX) {
                bail!("reference position exceeds index format limit");
            }
            previous_ref = ref_id as u32;
            previous_pos = reference_pos as u32;
            let packed = pack_occurrence(previous_ref, previous_pos).to_le_bytes();
            append_cache_bytes(&mut occurrence_writer, &mut occurrence_buffer, &packed)?;
            occurrence_bytes_written = occurrence_bytes_written
                .checked_add(8)
                .context("index cache occurrence length overflow")?;
        }
        if let Some(offsets) = compact_group_offsets.as_mut() {
            offsets.push(decoded_occurrences)?;
        }
    }

    if decoded_occurrences != occurrence_count {
        bail!("index seed record count mismatch");
    }
    finish_streaming_lookup_offsets(&mut lookup_values, &mut next_lookup_bucket, key_count);
    flush_cache_bytes(&mut key_writer, &mut key_buffer)?;
    flush_cache_bytes(&mut occurrence_writer, &mut occurrence_buffer)?;
    key_writer.flush()?;
    occurrence_writer.flush()?;
    drop(key_writer);
    drop(occurrence_writer);

    if key_bytes_written != ranges.seed_keys.end - ranges.seed_keys.start {
        bail!("index cache key array length mismatch");
    }
    if occurrence_bytes_written != ranges.occurrences.end - ranges.occurrences.start {
        bail!("index cache occurrence array length mismatch");
    }

    write_u64_slice_at(writer, ranges.lookup_offsets.clone(), &lookup_values)?;
    if let Some(offsets) = &compact_group_offsets {
        write_group_offsets_at(writer, ranges.group_offsets.clone(), offsets)?;
    }
    writer.seek(SeekFrom::Start(ranges.fingerprints.start))?;
    let mut fingerprint_pos = ranges.fingerprints.start;
    copy_exact_cache_bytes(
        reader,
        writer,
        &mut fingerprint_pos,
        fingerprint_len as u64,
        &mut copy_buffer,
    )?;
    if fingerprint_pos != ranges.fingerprints.end {
        bail!("index cache fingerprint array length mismatch");
    }
    let mut trailing = [0u8; 1];
    if reader.read(&mut trailing)? != 0 {
        bail!("index contains trailing decompressed data");
    }

    write_u64_at(
        writer,
        ranges.total_reference_bases_patch,
        total_reference_bases,
    )?;
    write_u64_at(writer, ranges.indexed_seed_hashes_patch, group_count_u64)?;
    write_u64_at(
        writer,
        ranges.indexed_seed_occurrences_patch,
        occurrence_count_u64,
    )?;
    write_u64_at(
        writer,
        ranges.reference_storage_bytes_patch,
        reference_storage_bytes,
    )?;
    write_u64_at(
        writer,
        ranges.reference_metadata_offset_patch,
        ranges.reference_metadata_offset,
    )?;
    write_u64_at(
        writer,
        ranges.reference_metadata_len_patch,
        ranges.reference_metadata_len,
    )?;
    writer.flush()?;
    Ok(())
}

fn cache_array_bytes(count: usize, width: usize) -> Result<u64> {
    let bytes = count
        .checked_mul(width)
        .context("index cache array length overflow")?;
    u64::try_from(bytes).context("index cache array length exceeds file format")
}

fn reserve_cache_range<W: Write + Seek>(
    writer: &mut W,
    pos: &mut u64,
    len: u64,
) -> Result<std::ops::Range<u64>> {
    let start = *pos;
    let end = start
        .checked_add(len)
        .context("index cache file length overflow")?;
    writer.seek(SeekFrom::Start(end))?;
    *pos = end;
    Ok(start..end)
}

fn copy_exact_cache_bytes<R: Read, W: Write>(
    reader: &mut R,
    writer: &mut W,
    pos: &mut u64,
    len: u64,
    buffer: &mut [u8],
) -> Result<()> {
    let mut remaining = len;
    while remaining > 0 {
        let chunk = usize::try_from(remaining.min(buffer.len() as u64))
            .context("index cache copy length exceeds address space")?;
        reader.read_exact(&mut buffer[..chunk])?;
        writer.write_all(&buffer[..chunk])?;
        remaining -= chunk as u64;
    }
    *pos = pos
        .checked_add(len)
        .context("index cache file length overflow")?;
    Ok(())
}

fn append_repeated_cache_bytes<W: Write>(
    writer: &mut W,
    buffer: &mut Vec<u8>,
    pattern: &[u8],
    count: usize,
) -> Result<()> {
    for _ in 0..count {
        append_cache_bytes(writer, buffer, pattern)?;
    }
    Ok(())
}

fn append_cache_bytes<W: Write>(writer: &mut W, buffer: &mut Vec<u8>, bytes: &[u8]) -> Result<()> {
    if buffer.len().saturating_add(bytes.len()) > INDEX_CACHE_IO_BUFFER_BYTES {
        flush_cache_bytes(writer, buffer)?;
    }
    buffer.extend_from_slice(bytes);
    Ok(())
}

fn flush_cache_bytes<W: Write>(writer: &mut W, buffer: &mut Vec<u8>) -> Result<()> {
    if !buffer.is_empty() {
        writer.write_all(buffer)?;
        buffer.clear();
    }
    Ok(())
}

fn update_streaming_lookup_offsets(
    offsets: &mut [u64],
    next_bucket: &mut usize,
    hash_bits: u8,
    prefix_bits: u8,
    hash: u64,
    row: usize,
) -> Result<()> {
    if offsets.is_empty() {
        return Ok(());
    }
    let suffix_bits = hash_bits - prefix_bits;
    let bucket = usize::try_from(hash >> suffix_bits)
        .context("index cache lookup bucket exceeds address space")?;
    let bucket_count = offsets.len() - 1;
    if bucket >= bucket_count {
        bail!("index cache lookup bucket exceeds prefix table");
    }
    while *next_bucket <= bucket {
        offsets[*next_bucket] = row as u64;
        *next_bucket += 1;
    }
    Ok(())
}

fn finish_streaming_lookup_offsets(offsets: &mut [u64], next_bucket: &mut usize, rows: usize) {
    while *next_bucket < offsets.len() {
        offsets[*next_bucket] = rows as u64;
        *next_bucket += 1;
    }
}

fn write_u64_slice_at<W: Write + Seek>(
    writer: &mut W,
    range: std::ops::Range<u64>,
    values: &[u64],
) -> Result<()> {
    writer.seek(SeekFrom::Start(range.start))?;
    let mut pos = range.start;
    write_aligned_u64_slice(writer, &mut pos, values)?;
    if pos != range.end {
        bail!("index cache u64 array length mismatch");
    }
    Ok(())
}

fn write_group_offsets_at<W: Write + Seek>(
    writer: &mut W,
    range: std::ops::Range<u64>,
    offsets: &CompactGroupOffsets,
) -> Result<()> {
    writer.seek(SeekFrom::Start(range.start))?;
    let mut pos = range.start;
    match offsets {
        CompactGroupOffsets::U32(values) => write_aligned_u32_slice(writer, &mut pos, values)?,
        CompactGroupOffsets::U64(values) => write_aligned_u64_slice(writer, &mut pos, values)?,
    }
    if pos != range.end {
        bail!("index cache group-offset array length mismatch");
    }
    Ok(())
}

fn write_u64_at<W: Write + Seek>(writer: &mut W, offset: u64, value: u64) -> Result<()> {
    writer.seek(SeekFrom::Start(offset))?;
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

#[cfg(test)]
fn write_index_cache(path: &str, source: IndexSourceFingerprint, index: &Index) -> Result<()> {
    let path = Path::new(path);
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).with_context(|| {
            format!(
                "failed to create index cache directory {}",
                parent.display()
            )
        })?;
    }
    let tmp_path = cache_tmp_path(path);
    let _tmp_cleanup = CacheTempFileCleanup(tmp_path.clone());
    {
        let file = File::create(&tmp_path)
            .with_context(|| format!("failed to create index cache {}", tmp_path.display()))?;
        let mut writer = BufWriter::new(file);
        write_mapped_index_cache(&mut writer, source, index)
            .with_context(|| format!("failed to write index cache {}", tmp_path.display()))?;
        writer.flush()?;
    }
    fs::rename(&tmp_path, path).with_context(|| {
        format!(
            "failed to move index cache {} to {}",
            tmp_path.display(),
            path.display()
        )
    })?;
    Ok(())
}

fn cache_tmp_path(path: &Path) -> PathBuf {
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("metam-index-cache");
    path.with_file_name(format!("{file_name}.tmp.{}", process::id()))
}

struct CacheTempFileCleanup(PathBuf);

impl Drop for CacheTempFileCleanup {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

#[cfg(test)]
fn write_mapped_index_cache<W: Write>(
    writer: &mut W,
    source: IndexSourceFingerprint,
    index: &Index,
) -> Result<()> {
    if !cfg!(target_endian = "little") {
        bail!("memory-mapped index cache requires a little-endian target");
    }
    let mut pos = 0u64;
    write_raw(writer, &mut pos, &INDEX_CACHE_MAGIC)?;
    write_u32_le(writer, &mut pos, INDEX_CACHE_VERSION)?;
    write_u64_le(writer, &mut pos, source.size)?;
    write_u64_le(writer, &mut pos, source.modified_secs)?;
    write_u32_le(writer, &mut pos, source.modified_nanos)?;

    write_raw(writer, &mut pos, &index.header.magic)?;
    write_u32_le(writer, &mut pos, index.header.version)?;
    write_cache_string(writer, &mut pos, &index.header.metam_version)?;
    write_u64_le(writer, &mut pos, index.header.k as u64)?;
    write_u64_le(writer, &mut pos, index.header.w as u64)?;
    write_u8(
        writer,
        &mut pos,
        seed_mode_code(index.header.seed_mode) as u8,
    )?;
    write_u64_le(writer, &mut pos, index.header.max_occ as u64)?;
    write_u8(writer, &mut pos, index.header.hash_bits)?;
    write_u8(
        writer,
        &mut pos,
        u8::from(index.header.exact_randstrobe_match),
    )?;
    let (strobe_start_dist_min, strobe_start_dist_max) = index
        .header
        .strobe_start_dist
        .map(|distance| (distance.min, distance.max))
        .unwrap_or((0, 0));
    write_u64_le(writer, &mut pos, u64::from(strobe_start_dist_min))?;
    write_u64_le(writer, &mut pos, u64::from(strobe_start_dist_max))?;
    write_u64_le(writer, &mut pos, index.header.reference_count as u64)?;
    write_u64_le(writer, &mut pos, index.header.total_reference_bases)?;
    write_u64_le(writer, &mut pos, index.header.indexed_seed_hashes)?;
    write_u64_le(writer, &mut pos, index.header.indexed_seed_occurrences)?;
    write_u64_le(writer, &mut pos, index.header.filtered_seed_hashes)?;
    write_u64_le(writer, &mut pos, index.header.reference_storage_bytes)?;
    write_u64_le(writer, &mut pos, index.header.created_unix_seconds)?;

    write_u64_le(writer, &mut pos, index.references.len() as u64)?;
    let reference_storage_len = index.references.iter().try_fold(0u64, |total, reference| {
        total
            .checked_add(reference.packed_bases.len() as u64)
            .and_then(|value| value.checked_add(reference.ambiguous.byte_len() as u64))
            .context("index cache reference storage length overflow")
    })?;
    let reference_metadata_len = index.references.iter().try_fold(0u64, |total, reference| {
        total
            .checked_add(32)
            .and_then(|value| value.checked_add(reference.id.len() as u64))
            .context("index cache reference-metadata length overflow")
    })?;
    let reference_metadata_offset = pos
        .checked_add(16)
        .and_then(|value| value.checked_add(reference_storage_len))
        .context("index cache reference-metadata offset overflow")?;
    write_u64_le(writer, &mut pos, reference_metadata_offset)?;
    write_u64_le(writer, &mut pos, reference_metadata_len)?;
    for reference in &index.references {
        write_raw(writer, &mut pos, &reference.packed_bases)?;
        write_raw(writer, &mut pos, &reference.ambiguous.clone().into_bytes())?;
    }
    if pos != reference_metadata_offset {
        bail!("index cache reference storage length mismatch");
    }
    for reference in &index.references {
        write_cache_string(writer, &mut pos, &reference.id)?;
        write_u64_le(writer, &mut pos, reference.len)?;
        write_u64_le(writer, &mut pos, reference.packed_bases.len() as u64)?;
        write_u64_le(writer, &mut pos, reference.ambiguous.byte_len() as u64)?;
    }
    let reference_metadata_end = reference_metadata_offset
        .checked_add(reference_metadata_len)
        .context("index cache reference-metadata end overflow")?;
    if pos != reference_metadata_end {
        bail!("index cache reference-metadata length mismatch");
    }

    let cache_seed_keys = CacheSeedKeys::for_index(index)?;
    let (key_layout, lookup_index, hashes32, hashes64) = match &cache_seed_keys {
        CacheSeedKeys::Repeated => (
            CacheSeedKeyLayout::Repeated,
            &index.seed_table.lookup_index,
            index.seed_table.hashes32.as_slice(),
            index.seed_table.hashes64.as_slice(),
        ),
        CacheSeedKeys::Grouped(compact) => (
            CacheSeedKeyLayout::Grouped,
            &compact.lookup_index,
            compact.hashes32.as_slice(),
            compact.hashes64.as_slice(),
        ),
    };

    write_u8(writer, &mut pos, index.seed_table.hash_bits)?;
    write_u8(writer, &mut pos, key_layout as u8)?;
    write_u8(writer, &mut pos, lookup_index.prefix_bits)?;
    write_u64_le(writer, &mut pos, lookup_index.offsets.len() as u64)?;
    write_aligned_u64_slice(writer, &mut pos, &lookup_index.offsets)?;
    if index.seed_table.hash_bits <= 32 {
        write_u64_le(writer, &mut pos, hashes32.len() as u64)?;
        write_aligned_u32_slice(writer, &mut pos, hashes32)?;
        write_u64_le(writer, &mut pos, 0)?;
        align_writer(writer, &mut pos, 8)?;
    } else {
        write_u64_le(writer, &mut pos, 0)?;
        align_writer(writer, &mut pos, 4)?;
        write_u64_le(writer, &mut pos, hashes64.len() as u64)?;
        write_aligned_u64_slice(writer, &mut pos, hashes64)?;
    }
    match &cache_seed_keys {
        CacheSeedKeys::Repeated => {
            write_u8(writer, &mut pos, 0)?;
            write_u64_le(writer, &mut pos, 0)?;
        }
        CacheSeedKeys::Grouped(compact) => {
            write_u8(writer, &mut pos, compact.group_offsets.width())?;
            write_u64_le(writer, &mut pos, compact.group_offsets.len() as u64)?;
            match &compact.group_offsets {
                CompactGroupOffsets::U32(offsets) => {
                    write_aligned_u32_slice(writer, &mut pos, offsets)?
                }
                CompactGroupOffsets::U64(offsets) => {
                    write_aligned_u64_slice(writer, &mut pos, offsets)?
                }
            }
        }
    }
    write_u64_le(writer, &mut pos, index.seed_table.occurrences.len() as u64)?;
    write_aligned_u64_slice(writer, &mut pos, &index.seed_table.occurrences)?;
    write_u8(writer, &mut pos, index.seed_table.fingerprint_bits)?;
    write_u64_le(writer, &mut pos, index.seed_table.fingerprints.len() as u64)?;
    write_raw(writer, &mut pos, &index.seed_table.fingerprints)?;
    Ok(())
}

fn write_raw<W: Write>(writer: &mut W, pos: &mut u64, bytes: &[u8]) -> Result<()> {
    writer.write_all(bytes)?;
    *pos += bytes.len() as u64;
    Ok(())
}

fn write_u8<W: Write>(writer: &mut W, pos: &mut u64, value: u8) -> Result<()> {
    write_raw(writer, pos, &[value])
}

fn write_u32_le<W: Write>(writer: &mut W, pos: &mut u64, value: u32) -> Result<()> {
    write_raw(writer, pos, &value.to_le_bytes())
}

fn write_u64_le<W: Write>(writer: &mut W, pos: &mut u64, value: u64) -> Result<()> {
    write_raw(writer, pos, &value.to_le_bytes())
}

fn write_cache_string<W: Write>(writer: &mut W, pos: &mut u64, value: &str) -> Result<()> {
    write_cache_bytes(writer, pos, value.as_bytes())
}

fn write_cache_bytes<W: Write>(writer: &mut W, pos: &mut u64, bytes: &[u8]) -> Result<()> {
    write_u64_le(writer, pos, bytes.len() as u64)?;
    write_raw(writer, pos, bytes)
}

fn align_writer<W: Write>(writer: &mut W, pos: &mut u64, align: u64) -> Result<()> {
    let pad = ((*pos).next_multiple_of(align) - *pos) as usize;
    if pad > 0 {
        write_raw(writer, pos, &[0; 8][..pad])?;
    }
    Ok(())
}

fn write_aligned_u32_slice<W: Write>(writer: &mut W, pos: &mut u64, values: &[u32]) -> Result<()> {
    align_writer(writer, pos, std::mem::align_of::<u32>() as u64)?;
    let bytes = unsafe {
        std::slice::from_raw_parts(values.as_ptr().cast::<u8>(), std::mem::size_of_val(values))
    };
    write_raw(writer, pos, bytes)
}

fn write_aligned_u64_slice<W: Write>(writer: &mut W, pos: &mut u64, values: &[u64]) -> Result<()> {
    align_writer(writer, pos, std::mem::align_of::<u64>() as u64)?;
    let bytes = unsafe {
        std::slice::from_raw_parts(values.as_ptr().cast::<u8>(), std::mem::size_of_val(values))
    };
    write_raw(writer, pos, bytes)
}

impl MappedIndex {
    fn open(path: &str, source: IndexSourceFingerprint) -> Result<Self> {
        if !cfg!(target_endian = "little") {
            bail!("memory-mapped index cache requires a little-endian target");
        }
        let mmap = MmapFile::map(path)?;
        let mut cursor = CacheCursor::new(mmap.as_slice());
        let magic = cursor.read_array::<8>()?;
        if magic != INDEX_CACHE_MAGIC {
            bail!("invalid index cache magic");
        }
        let version = cursor.read_u32()?;
        if version != INDEX_CACHE_VERSION {
            bail!("unsupported index cache version");
        }
        let found_source = IndexSourceFingerprint {
            size: cursor.read_u64()?,
            modified_secs: cursor.read_u64()?,
            modified_nanos: cursor.read_u32()?,
        };
        if found_source != source {
            bail!("stale index cache");
        }

        let header_magic = cursor.read_array::<8>()?;
        let header_version = cursor.read_u32()?;
        validate_index_identity(header_magic, header_version)?;
        let metam_version = cursor.read_string()?;
        let k = cursor.read_u64()? as usize;
        let w = cursor.read_u64()? as usize;
        let seed_mode = seed_mode_from_code(u64::from(cursor.read_u8()?))?;
        let max_occ = cursor.read_u64()? as usize;
        let hash_bits = cursor.read_u8()?;
        validate_hash_bits(hash_bits)?;
        let exact_randstrobe_match = match cursor.read_u8()? {
            0 => false,
            1 => true,
            _ => bail!("invalid exact randstrobe match flag in index cache"),
        };
        validate_exact_randstrobe_match(seed_mode, exact_randstrobe_match)?;
        let strobe_start_dist_min = cursor.read_u64()?;
        let strobe_start_dist_max = cursor.read_u64()?;
        let strobe_start_dist = if strobe_start_dist_min == 0 && strobe_start_dist_max == 0 {
            None
        } else {
            let min = u16::try_from(strobe_start_dist_min)
                .context("index cache randstrobe start-distance min exceeds u16")?;
            let max = u16::try_from(strobe_start_dist_max)
                .context("index cache randstrobe start-distance max exceeds u16")?;
            Some(StrobeStartDistance::new(min, max)?)
        };
        validate_strobe_start_distance_metadata(header_version, seed_mode, strobe_start_dist)?;
        let reference_count = cursor.read_u64()? as usize;
        let total_reference_bases = cursor.read_u64()?;
        let indexed_seed_hashes = cursor.read_u64()?;
        let indexed_seed_occurrences = cursor.read_u64()?;
        let filtered_seed_hashes = cursor.read_u64()?;
        let reference_storage_bytes = cursor.read_u64()?;
        let created_unix_seconds = cursor.read_u64()?;

        let ref_count = usize::try_from(cursor.read_u64()?)
            .context("index cache reference count exceeds address space")?;
        if ref_count != reference_count {
            bail!("index cache reference count mismatch");
        }
        let reference_metadata_offset = usize::try_from(cursor.read_u64()?)
            .context("index cache reference-metadata offset exceeds address space")?;
        let reference_metadata_len = usize::try_from(cursor.read_u64()?)
            .context("index cache reference-metadata length exceeds address space")?;
        let reference_data_start = cursor.position();
        if reference_metadata_offset < reference_data_start {
            bail!("index cache reference metadata overlaps its header");
        }
        let reference_metadata_end = reference_metadata_offset
            .checked_add(reference_metadata_len)
            .context("index cache reference-metadata range overflow")?;
        if reference_metadata_end > mmap.as_slice().len() {
            bail!("truncated index cache reference metadata");
        }
        let stored_reference_bytes = reference_metadata_offset - reference_data_start;
        if u64::try_from(stored_reference_bytes).ok() != Some(reference_storage_bytes) {
            bail!("index cache reference storage length mismatch");
        }

        let mut reference_data = CacheCursor::at(mmap.as_slice(), reference_data_start)?;
        let mut reference_metadata =
            CacheCursor::new(&mmap.as_slice()[reference_metadata_offset..reference_metadata_end]);
        let mut references = Vec::with_capacity(ref_count);
        let mut found_total_reference_bases = 0u64;
        for _ in 0..ref_count {
            let id = reference_metadata.read_string()?;
            let len = reference_metadata.read_u64()?;
            let packed_len = reference_metadata.read_u64()?;
            if packed_len != len.div_ceil(4) {
                bail!("index cache reference has an invalid packed sequence length");
            }
            let ambiguous_len = reference_metadata.read_u64()?;
            let expected_ambiguous_len = len.div_ceil(8);
            if ambiguous_len != 0 && ambiguous_len != expected_ambiguous_len {
                bail!("index cache reference has an invalid ambiguity-mask length");
            }
            let packed = reference_data.read_range(
                usize::try_from(packed_len)
                    .context("index cache packed sequence exceeds address space")?,
            )?;
            let ambiguous = reference_data.read_range(
                usize::try_from(ambiguous_len)
                    .context("index cache ambiguity mask exceeds address space")?,
            )?;
            found_total_reference_bases = found_total_reference_bases
                .checked_add(len)
                .context("index cache total reference length overflow")?;
            references.push(MappedReferenceMeta {
                id,
                len,
                packed,
                ambiguous,
            });
        }
        if reference_data.position() != reference_metadata_offset {
            bail!("index cache reference data does not end at its metadata section");
        }
        if reference_metadata.position() != reference_metadata_len {
            bail!("index cache reference-metadata length mismatch");
        }
        if found_total_reference_bases != total_reference_bases {
            bail!("index cache total reference length mismatch");
        }
        cursor = CacheCursor::at(mmap.as_slice(), reference_metadata_end)?;
        let seed_hash_bits = cursor.read_u8()?;
        if seed_hash_bits != hash_bits {
            bail!("index cache seed hash width mismatch");
        }
        let key_layout = CacheSeedKeyLayout::from_code(cursor.read_u8()?)?;
        let prefix_bits = cursor.read_u8()?;
        let offset_len = cursor.read_u64()? as usize;
        let offsets = cursor.read_aligned_range::<u64>(offset_len)?;
        let hashes32_len = cursor.read_u64()? as usize;
        let hashes32 = cursor.read_aligned_range::<u32>(hashes32_len)?;
        let hashes64_len = cursor.read_u64()? as usize;
        let hashes64 = cursor.read_aligned_range::<u64>(hashes64_len)?;
        let group_offset_width = cursor.read_u8()?;
        let group_offset_len = cursor.read_u64()? as usize;
        let group_offsets = match group_offset_width {
            0 => cursor.read_range(0)?,
            4 => cursor.read_aligned_range::<u32>(group_offset_len)?,
            8 => cursor.read_aligned_range::<u64>(group_offset_len)?,
            _ => bail!("invalid seed-group offset width in index cache"),
        };
        let occurrence_len = cursor.read_u64()? as usize;
        let occurrences = cursor.read_aligned_range::<u64>(occurrence_len)?;
        let fingerprint_bits = cursor.read_u8()?;
        let fingerprint_len = cursor.read_u64()? as usize;
        let fingerprints = cursor.read_range(fingerprint_len)?;

        let header = IndexHeader {
            magic: header_magic,
            version: header_version,
            metam_version,
            k,
            w,
            seed_mode,
            max_occ,
            hash_bits,
            exact_randstrobe_match,
            strobe_start_dist,
            reference_count,
            total_reference_bases,
            indexed_seed_hashes,
            indexed_seed_occurrences,
            filtered_seed_hashes,
            reference_storage_bytes,
            created_unix_seconds,
        };
        let index = Self {
            mmap,
            header,
            references,
            seed_table: MappedSeedTable {
                hash_bits,
                key_layout,
                hashes32,
                hashes32_len,
                hashes64,
                hashes64_len,
                group_offset_width,
                group_offsets,
                group_offset_len,
                occurrences,
                occurrence_len,
                fingerprint_bits,
                fingerprints,
                fingerprint_len,
                prefix_bits,
                offsets,
                offset_len,
            },
        };
        index.validate()?;
        index.prefetch_seed_lookup_metadata();
        Ok(index)
    }

    fn prefetch_seed_lookup_metadata(&self) {
        #[cfg(unix)]
        {
            // These compact arrays are small compared with occurrences and
            // are consulted for every seed lookup. Prefetching them costs only
            // their actual pages (roughly a few hundred MiB on the production
            // index). The rest of the mmap retains the kernel's normal access
            // policy so cold NFS reads can use mmap read-around and
            // readahead.
            for range in self.seed_lookup_metadata_ranges() {
                self.mmap.advise(range, libc::MADV_WILLNEED);
            }
        }
    }

    #[cfg(unix)]
    fn seed_lookup_metadata_ranges(&self) -> Vec<std::ops::Range<usize>> {
        let mut ranges = Vec::with_capacity(3);
        push_nonempty_range(&mut ranges, self.seed_table.offsets.clone());
        if self.seed_table.hash_bits <= 32 {
            push_nonempty_range(&mut ranges, self.seed_table.hashes32.clone());
        } else {
            push_nonempty_range(&mut ranges, self.seed_table.hashes64.clone());
        }
        push_nonempty_range(&mut ranges, self.seed_table.group_offsets.clone());
        ranges
    }

    fn validate(&self) -> Result<()> {
        validate_index_identity(self.header.magic, self.header.version)?;
        if self.header.reference_count != self.references.len() {
            bail!("index metadata/reference count mismatch");
        }
        validate_hash_bits(self.header.hash_bits)?;
        validate_exact_randstrobe_match(self.header.seed_mode, self.header.exact_randstrobe_match)?;
        let indexed_seed_occurrences = usize::try_from(self.header.indexed_seed_occurrences)
            .context("index cache seed occurrence count exceeds address space")?;
        if indexed_seed_occurrences != self.seed_table.occurrence_len {
            bail!("index metadata/seed occurrence count mismatch");
        }
        let indexed_seed_hashes = usize::try_from(self.header.indexed_seed_hashes)
            .context("index cache distinct seed-key count exceeds address space")?;
        let hash_count = self.seed_table.hash_count();
        let expected_fingerprint_bits = if self.header.exact_randstrobe_match {
            randstrobe_fingerprint_bits(self.header.hash_bits)
        } else {
            0
        };
        if self.seed_table.fingerprint_bits != expected_fingerprint_bits {
            bail!("index cache randstrobe fingerprint width mismatch");
        }
        let expected_fingerprint_len = self
            .seed_table
            .occurrence_len
            .checked_mul(fingerprint_bytes(expected_fingerprint_bits))
            .context("index cache fingerprint length overflow")?;
        if self.seed_table.fingerprint_len != expected_fingerprint_len {
            bail!("index cache randstrobe fingerprint array length mismatch");
        }
        if self.seed_table.hash_bits <= 32 {
            if self.seed_table.hashes64_len != 0 {
                bail!("index cache contains unexpected 64-bit seed hashes");
            }
            if self.seed_table.hashes32_len == 0 && indexed_seed_hashes > 0 {
                bail!("index cache missing 32-bit seed hashes");
            }
        } else {
            if self.seed_table.hashes32_len != 0 {
                bail!("index cache contains unexpected 32-bit seed hashes");
            }
            if self.seed_table.hashes64_len == 0 && indexed_seed_hashes > 0 {
                bail!("index cache missing 64-bit seed hashes");
            }
        }

        match self.seed_table.key_layout {
            CacheSeedKeyLayout::Repeated => {
                if hash_count != self.seed_table.occurrence_len {
                    bail!("repeated index cache seed keys do not match occurrence count");
                }
                if self.seed_table.group_offset_width != 0 || self.seed_table.group_offset_len != 0
                {
                    bail!("repeated index cache unexpectedly contains seed-group offsets");
                }
            }
            CacheSeedKeyLayout::Grouped => {
                if hash_count != indexed_seed_hashes {
                    bail!("grouped index cache key count does not match index metadata");
                }
                let expected_group_offset_len = hash_count
                    .checked_add(1)
                    .context("index cache seed-group count overflow")?;
                if self.seed_table.group_offset_len != expected_group_offset_len {
                    bail!("grouped index cache has an invalid seed-group offset count");
                }
                let expected_width = if u32::try_from(self.seed_table.occurrence_len).is_ok() {
                    4
                } else {
                    8
                };
                if self.seed_table.group_offset_width != expected_width {
                    bail!("grouped index cache has a non-canonical seed-group offset width");
                }
                let first = self
                    .seed_table
                    .group_offset_at(self.mmap.as_slice(), 0)
                    .context("grouped index cache is missing its first seed-group offset")?;
                let last = self
                    .seed_table
                    .group_offset_at(self.mmap.as_slice(), hash_count)
                    .context("grouped index cache is missing its final seed-group offset")?;
                if first != 0 || last != self.seed_table.occurrence_len {
                    bail!("grouped index cache seed-group endpoints are invalid");
                }
            }
        }

        let expected_prefix_bits = if hash_count == 0 {
            0
        } else {
            self.seed_table.hash_bits.min(20)
        };
        if self.seed_table.prefix_bits != expected_prefix_bits {
            bail!("index cache seed lookup prefix width mismatch");
        }
        let expected_offset_len = if hash_count == 0 {
            0
        } else {
            (1usize << expected_prefix_bits) + 1
        };
        if self.seed_table.offset_len != expected_offset_len {
            bail!("index cache seed lookup offset count mismatch");
        }

        let _ = self.seed_table.offsets(self.mmap.as_slice());
        let _ = self.seed_table.occurrences(self.mmap.as_slice());
        let _ = self.seed_table.fingerprints(self.mmap.as_slice());
        if self.seed_table.hash_bits <= 32 {
            let _ = self.seed_table.hashes32(self.mmap.as_slice());
        } else {
            let _ = self.seed_table.hashes64(self.mmap.as_slice());
        }
        if hash_count > 0 {
            let offsets = self.seed_table.offsets(self.mmap.as_slice());
            if offsets.first().copied() != Some(0)
                || offsets.last().copied() != Some(hash_count as u64)
            {
                bail!("index cache seed lookup endpoints are invalid");
            }
        }
        Ok(())
    }

    fn reference(&self, ref_id: u32) -> MappedReference<'_> {
        let meta = &self.references[ref_id as usize];
        let bytes = self.mmap.as_slice();
        MappedReference {
            id: &meta.id,
            len: meta.len,
            packed_bases: &bytes[meta.packed.clone()],
            ambiguous: &bytes[meta.ambiguous.clone()],
        }
    }

    #[cfg(test)]
    fn lookup(&self, hash: u64) -> SeedOccurrences<'_> {
        let key = seed_key(hash, self.header.hash_bits);
        self.seed_table.lookup(self.mmap.as_slice(), key, hash)
    }

    fn lookup_with_mode(&self, hash: u64, mode: SeedLookupMode) -> SeedOccurrences<'_> {
        let key = seed_key(hash, self.header.hash_bits);
        self.seed_table
            .lookup_with_mode(self.mmap.as_slice(), key, hash, mode)
    }
}

impl MappedSeedTable {
    fn hash_count(&self) -> usize {
        if self.hash_bits <= 32 {
            self.hashes32_len
        } else {
            self.hashes64_len
        }
    }

    #[cfg(test)]
    fn lookup<'a>(&self, bytes: &'a [u8], key: u64, hash: u64) -> SeedOccurrences<'a> {
        self.lookup_with_mode(bytes, key, hash, SeedLookupMode::IndexDefault)
    }

    fn lookup_with_mode<'a>(
        &self,
        bytes: &'a [u8],
        key: u64,
        hash: u64,
        mode: SeedLookupMode,
    ) -> SeedOccurrences<'a> {
        if self.hash_bits <= 32 {
            let Ok(key) = u32::try_from(key) else {
                return empty_seed_occurrences();
            };
            let Some((lo, hi)) = self.lookup_bounds(bytes, u64::from(key)) else {
                return empty_seed_occurrences();
            };
            let all_hashes = self.hashes32(bytes);
            let hashes = &all_hashes[lo..hi];
            let start = lo + hashes.partition_point(|h| *h < key);
            let (start, end) = match self.key_layout {
                CacheSeedKeyLayout::Repeated => {
                    let end = lo + hashes.partition_point(|h| *h <= key);
                    (start, end)
                }
                CacheSeedKeyLayout::Grouped => {
                    if start >= hi || all_hashes[start] != key {
                        return empty_seed_occurrences();
                    }
                    let Some(bounds) = self.group_occurrence_bounds(bytes, start) else {
                        return empty_seed_occurrences();
                    };
                    bounds
                }
            };
            self.seed_occurrences(bytes, hash, start, end, mode)
        } else {
            let Some((lo, hi)) = self.lookup_bounds(bytes, key) else {
                return empty_seed_occurrences();
            };
            let all_hashes = self.hashes64(bytes);
            let hashes = &all_hashes[lo..hi];
            let start = lo + hashes.partition_point(|h| *h < key);
            let (start, end) = match self.key_layout {
                CacheSeedKeyLayout::Repeated => {
                    let end = lo + hashes.partition_point(|h| *h <= key);
                    (start, end)
                }
                CacheSeedKeyLayout::Grouped => {
                    if start >= hi || all_hashes[start] != key {
                        return empty_seed_occurrences();
                    }
                    let Some(bounds) = self.group_occurrence_bounds(bytes, start) else {
                        return empty_seed_occurrences();
                    };
                    bounds
                }
            };
            self.seed_occurrences(bytes, hash, start, end, mode)
        }
    }

    fn seed_occurrences<'a>(
        &self,
        bytes: &'a [u8],
        hash: u64,
        start: usize,
        end: usize,
        mode: SeedLookupMode,
    ) -> SeedOccurrences<'a> {
        let fingerprint_bits = if mode.use_fingerprint_filter() {
            self.fingerprint_bits
        } else {
            0
        };
        let width = fingerprint_bytes(fingerprint_bits);
        let fingerprints = if width == 0 {
            &[]
        } else {
            &self.fingerprints(bytes)[start * width..end * width]
        };
        SeedOccurrences {
            data: &self.occurrences(bytes)[start..end],
            fingerprints,
            fingerprint_bits,
            target_fingerprint: if fingerprint_bits == 0 {
                0
            } else {
                randstrobe_fingerprint(hash, self.hash_bits)
            },
        }
    }

    fn lookup_bounds(&self, bytes: &[u8], key: u64) -> Option<(usize, usize)> {
        let hash_count = self.hash_count();
        if hash_count == 0 {
            return Some((0, 0));
        }
        if self.offset_len < 2 || self.prefix_bits > self.hash_bits {
            return None;
        }
        let suffix_bits = self.hash_bits - self.prefix_bits;
        let bucket = (key >> suffix_bits) as usize;
        let offsets = self.offsets(bytes);
        let lo = usize::try_from(*offsets.get(bucket)?).ok()?;
        let hi = usize::try_from(*offsets.get(bucket.checked_add(1)?)?).ok()?;
        (lo <= hi && hi <= hash_count).then_some((lo, hi))
    }

    fn group_occurrence_bounds(&self, bytes: &[u8], group: usize) -> Option<(usize, usize)> {
        let end_group = group.checked_add(1)?;
        let start = self.group_offset_at(bytes, group)?;
        let end = self.group_offset_at(bytes, end_group)?;
        (start <= end && end <= self.occurrence_len).then_some((start, end))
    }

    fn group_offset_at(&self, bytes: &[u8], index: usize) -> Option<usize> {
        if index >= self.group_offset_len {
            return None;
        }
        match self.group_offset_width {
            4 => usize::try_from(self.group_offsets32(bytes)[index]).ok(),
            8 => usize::try_from(self.group_offsets64(bytes)[index]).ok(),
            _ => None,
        }
    }

    fn hashes32<'a>(&self, bytes: &'a [u8]) -> &'a [u32] {
        typed_slice(bytes, self.hashes32.clone(), self.hashes32_len).expect("valid mmap hashes32")
    }

    fn hashes64<'a>(&self, bytes: &'a [u8]) -> &'a [u64] {
        typed_slice(bytes, self.hashes64.clone(), self.hashes64_len).expect("valid mmap hashes64")
    }

    fn occurrences<'a>(&self, bytes: &'a [u8]) -> &'a [u64] {
        typed_slice(bytes, self.occurrences.clone(), self.occurrence_len)
            .expect("valid mmap occurrences")
    }

    fn fingerprints<'a>(&self, bytes: &'a [u8]) -> &'a [u8] {
        &bytes[self.fingerprints.clone()]
    }

    fn offsets<'a>(&self, bytes: &'a [u8]) -> &'a [u64] {
        typed_slice(bytes, self.offsets.clone(), self.offset_len).expect("valid mmap offsets")
    }

    fn group_offsets32<'a>(&self, bytes: &'a [u8]) -> &'a [u32] {
        typed_slice(bytes, self.group_offsets.clone(), self.group_offset_len)
            .expect("valid mmap 32-bit group offsets")
    }

    fn group_offsets64<'a>(&self, bytes: &'a [u8]) -> &'a [u64] {
        typed_slice(bytes, self.group_offsets.clone(), self.group_offset_len)
            .expect("valid mmap 64-bit group offsets")
    }
}

struct CacheCursor<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> CacheCursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    fn at(bytes: &'a [u8], pos: usize) -> Result<Self> {
        if pos > bytes.len() {
            bail!("truncated index cache");
        }
        Ok(Self { bytes, pos })
    }

    fn position(&self) -> usize {
        self.pos
    }

    fn read_array<const N: usize>(&mut self) -> Result<[u8; N]> {
        Ok(self.read_exact(N)?.try_into().expect("fixed cache array"))
    }

    fn read_u8(&mut self) -> Result<u8> {
        Ok(self.read_exact(1)?[0])
    }

    fn read_u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.read_array::<4>()?))
    }

    fn read_u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.read_array::<8>()?))
    }

    fn read_string(&mut self) -> Result<String> {
        Ok(String::from_utf8(self.read_bytes()?.to_vec())?)
    }

    fn read_bytes(&mut self) -> Result<&'a [u8]> {
        let len = self.read_u64()? as usize;
        self.read_exact(len)
    }

    fn read_aligned_range<T>(&mut self, len: usize) -> Result<std::ops::Range<usize>> {
        let align = std::mem::align_of::<T>();
        self.pos = self.pos.next_multiple_of(align);
        let byte_len = len
            .checked_mul(std::mem::size_of::<T>())
            .context("index cache array length overflow")?;
        let range = self.read_range(byte_len)?;
        typed_slice::<T>(self.bytes, range.clone(), len)?;
        Ok(range)
    }

    fn read_exact(&mut self, len: usize) -> Result<&'a [u8]> {
        let range = self.read_range(len)?;
        Ok(&self.bytes[range])
    }

    fn read_range(&mut self, len: usize) -> Result<std::ops::Range<usize>> {
        let end = self
            .pos
            .checked_add(len)
            .context("index cache offset overflow")?;
        if end > self.bytes.len() {
            bail!("truncated index cache");
        }
        let range = self.pos..end;
        self.pos = end;
        Ok(range)
    }
}

fn typed_slice<T>(bytes: &[u8], range: std::ops::Range<usize>, len: usize) -> Result<&[T]> {
    let expected = len
        .checked_mul(std::mem::size_of::<T>())
        .context("index cache typed slice length overflow")?;
    if range.end > bytes.len() || range.len() != expected {
        bail!("invalid index cache typed slice range");
    }
    let ptr = bytes[range.clone()].as_ptr();
    if !(ptr as usize).is_multiple_of(std::mem::align_of::<T>()) {
        bail!("unaligned index cache typed slice");
    }
    Ok(unsafe { std::slice::from_raw_parts(ptr.cast::<T>(), len) })
}

fn validate_index(index: &Index) -> Result<()> {
    validate_index_identity(index.header.magic, index.header.version)?;
    if index.header.reference_count != index.references.len() {
        bail!("index metadata/reference count mismatch");
    }
    validate_hash_bits(index.header.hash_bits)?;
    validate_exact_randstrobe_match(index.header.seed_mode, index.header.exact_randstrobe_match)?;
    validate_strobe_start_distance_metadata(
        index.header.version,
        index.header.seed_mode,
        index.header.strobe_start_dist,
    )?;
    if index.header.indexed_seed_hashes as usize != index.seed_table.group_count() {
        bail!("index metadata/hash group count mismatch");
    }
    if index.header.indexed_seed_occurrences as usize != index.seed_table.occurrence_count() {
        bail!("index metadata/seed occurrence count mismatch");
    }
    if index.seed_table.hash_count() != index.seed_table.occurrence_count() {
        bail!("index seed table arrays have inconsistent lengths");
    }
    let expected_fingerprint_bits = if index.header.exact_randstrobe_match {
        randstrobe_fingerprint_bits(index.header.hash_bits)
    } else {
        0
    };
    if index.seed_table.fingerprint_bits != expected_fingerprint_bits {
        bail!("index randstrobe fingerprint width mismatch");
    }
    let expected_fingerprint_bytes = index
        .seed_table
        .expected_fingerprint_bytes()
        .context("index fingerprint array length overflow")?;
    if index.seed_table.fingerprints.len() != expected_fingerprint_bytes {
        bail!("index randstrobe fingerprint array length mismatch");
    }
    if index.seed_table.hash_bits != index.header.hash_bits {
        bail!("index seed table hash width mismatch");
    }
    if index.seed_table.hash_bits <= 32 {
        if !index.seed_table.hashes32.windows(2).all(|w| w[0] <= w[1]) {
            bail!("index seed hash rows are not sorted");
        }
    } else if !index.seed_table.hashes64.windows(2).all(|w| w[0] <= w[1]) {
        bail!("index seed hash rows are not sorted");
    }
    Ok(())
}

impl Index {
    pub fn lookup(&self, hash: u64) -> SeedOccurrences<'_> {
        self.lookup_with_mode(hash, SeedLookupMode::IndexDefault)
    }

    pub fn lookup_with_mode(&self, hash: u64, mode: SeedLookupMode) -> SeedOccurrences<'_> {
        let key = seed_key(hash, self.header.hash_bits);
        self.seed_table.lookup_with_mode(key, hash, mode)
    }
}

fn seed_key(hash: u64, hash_bits: u8) -> u64 {
    if hash_bits == 64 {
        hash
    } else {
        hash >> (64 - hash_bits)
    }
}

fn validate_exact_randstrobe_match(
    seed_mode: SeedMode,
    exact_randstrobe_match: bool,
) -> Result<()> {
    if exact_randstrobe_match && seed_mode != SeedMode::Randstrobe {
        bail!("--match-both-strobes requires --seed-mode randstrobe");
    }
    Ok(())
}

fn normalize_strobe_start_distance(
    seed_mode: SeedMode,
    strobe_start_dist: Option<StrobeStartDistance>,
) -> Result<Option<StrobeStartDistance>> {
    match seed_mode {
        SeedMode::Randstrobe => {
            let distance = strobe_start_dist.unwrap_or_default();
            StrobeStartDistance::new(distance.min, distance.max).map(Some)
        }
        SeedMode::Minimizer => Ok(None),
    }
}

fn validate_strobe_start_distance_metadata(
    version: u32,
    seed_mode: SeedMode,
    strobe_start_dist: Option<StrobeStartDistance>,
) -> Result<()> {
    match (seed_mode, strobe_start_dist) {
        (SeedMode::Minimizer, Some(_)) => {
            bail!("minimizer indexes must not record randstrobe start-distance settings")
        }
        (SeedMode::Randstrobe, None) if version >= 6 => {
            bail!("index version {version} missing randstrobe start-distance settings")
        }
        (_, Some(distance)) => {
            StrobeStartDistance::new(distance.min, distance.max)?;
            Ok(())
        }
        _ => Ok(()),
    }
}

fn validate_hash_bits(hash_bits: u8) -> Result<()> {
    if !(24..=64).contains(&hash_bits) {
        bail!("invalid hash_bits in index");
    }
    Ok(())
}

fn validate_index_identity(magic: [u8; 8], version: u32) -> Result<()> {
    if magic != INDEX_MAGIC {
        bail!("invalid index magic");
    }
    if !(MIN_SUPPORTED_INDEX_VERSION..=INDEX_VERSION).contains(&version) {
        bail!(
            "unsupported index version {version}; supported versions are {MIN_SUPPORTED_INDEX_VERSION} through {INDEX_VERSION}"
        );
    }
    Ok(())
}

fn sort_raw_seeds(raw: &mut [RawSeed], hash_bits: u8) {
    raw.par_sort_unstable_by_key(|seed| (seed_key(seed.hash, hash_bits), seed.ref_id, seed.pos));
}

fn randstrobe_fingerprint_bits(hash_bits: u8) -> u8 {
    RANDSTROBE_HASH_SPAN_BITS.saturating_sub(hash_bits)
}

fn fingerprint_bytes(fingerprint_bits: u8) -> usize {
    usize::from(fingerprint_bits.div_ceil(8))
}

fn randstrobe_fingerprint(hash: u64, hash_bits: u8) -> u32 {
    let bits = randstrobe_fingerprint_bits(hash_bits);
    if bits == 0 {
        return 0;
    }
    let mask = (1u64 << bits) - 1;
    ((hash >> RANDSTROBE_HASH_LOW_ZERO_BITS) & mask) as u32
}

fn push_packed_fingerprint(out: &mut Vec<u8>, fingerprint: u32, fingerprint_bits: u8) {
    let bytes = fingerprint.to_le_bytes();
    out.extend_from_slice(&bytes[..fingerprint_bytes(fingerprint_bits)]);
}

fn packed_fingerprint(data: &[u8], index: usize, fingerprint_bits: u8) -> u32 {
    let width = fingerprint_bytes(fingerprint_bits);
    if width == 0 {
        return 0;
    }
    let start = index * width;
    let mut bytes = [0u8; 4];
    bytes[..width].copy_from_slice(&data[start..start + width]);
    u32::from_le_bytes(bytes)
}

fn partition_for_hash(hash: u64, hash_bits: u8, partition_count: usize) -> usize {
    let partition = ((u128::from(hash) * partition_count as u128) >> hash_bits) as usize;
    partition.min(partition_count - 1)
}

fn pack_occurrence(ref_id: u32, pos: u32) -> u64 {
    (u64::from(ref_id) << 32) | u64::from(pos)
}

fn unpack_occurrence(packed: u64) -> SeedOccurrence {
    SeedOccurrence {
        ref_id: (packed >> 32) as u32,
        pos: packed as u32,
    }
}

fn write_index_payload<W: Write>(writer: &mut W, index: &Index) -> Result<()> {
    write_index_metadata(writer, &index.header)?;
    write_varint(writer, index.references.len() as u64)?;
    for reference in &index.references {
        write_reference(writer, reference)?;
    }
    write_seed_table(writer, &index.seed_table)?;
    Ok(())
}

fn write_streamed_index_payload(
    header: &IndexHeader,
    ref_store_path: &Path,
    seed_store_path: &Path,
    fingerprint_store_path: Option<&Path>,
    out_path: &str,
) -> Result<()> {
    write_compressed_index_atomic(Path::new(out_path), |encoder| {
        write_index_metadata(encoder, header)?;
        write_varint(encoder, header.reference_count as u64)?;
        let mut ref_store = BufReader::new(File::open(ref_store_path)?);
        copy(&mut ref_store, encoder)?;
        write_varint(encoder, header.indexed_seed_occurrences)?;
        write_varint(encoder, header.indexed_seed_hashes)?;
        let mut seed_store = BufReader::new(File::open(seed_store_path)?);
        copy(&mut seed_store, encoder)?;
        if let Some(path) = fingerprint_store_path {
            let mut fingerprint_store = BufReader::new(File::open(path)?);
            copy(&mut fingerprint_store, encoder)?;
        }
        Ok(())
    })
}

fn write_compressed_index_atomic<F>(path: &Path, write_payload: F) -> Result<()>
where
    F: FnOnce(&mut zstd::stream::write::Encoder<'static, BufWriter<File>>) -> Result<()>,
{
    let (temp_path, file) = create_sibling_temp_file(path)?;
    let mut cleanup = TempFileCleanup::new(temp_path.clone());
    if let Ok(metadata) = fs::metadata(path) {
        if metadata.is_file() {
            fs::set_permissions(&temp_path, metadata.permissions()).with_context(|| {
                format!(
                    "failed to preserve permissions for replacement index {}",
                    path.display()
                )
            })?;
        }
    }

    let writer = BufWriter::new(file);
    let mut encoder = zstd::stream::Encoder::new(writer, 9).with_context(|| {
        format!(
            "failed to initialize zstd index writer for {}",
            path.display()
        )
    })?;
    write_payload(&mut encoder)?;
    let mut writer = encoder
        .finish()
        .with_context(|| format!("failed to finalize index {}", path.display()))?;
    writer.flush()?;
    let file = writer
        .into_inner()
        .map_err(|error| error.into_error())
        .with_context(|| format!("failed to close index {}", path.display()))?;
    drop(file);

    fs::rename(&temp_path, path)
        .with_context(|| format!("failed to atomically publish index {}", path.display()))?;
    cleanup.disarm();
    Ok(())
}

fn create_sibling_temp_file(path: &Path) -> Result<(PathBuf, File)> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = path
        .file_name()
        .with_context(|| format!("index path {} has no file name", path.display()))?;
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or(0);
    for attempt in 0..100u32 {
        let mut temp_name = std::ffi::OsString::from(".");
        temp_name.push(name);
        temp_name.push(format!(".tmp-{}-{nanos}-{attempt}", process::id()));
        let temp_path = parent.join(temp_name);
        match OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&temp_path)
        {
            Ok(file) => return Ok((temp_path, file)),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(error).with_context(|| {
                    format!("failed to create temporary index beside {}", path.display())
                });
            }
        }
    }
    bail!(
        "failed to create a unique temporary index beside {}",
        path.display()
    )
}

fn write_index_metadata<W: Write>(writer: &mut W, header: &IndexHeader) -> Result<()> {
    writer.write_all(&header.magic)?;
    write_varint(writer, u64::from(header.version))?;
    write_string(writer, &header.metam_version)?;
    write_varint(writer, header.k as u64)?;
    write_varint(writer, header.w as u64)?;
    write_varint(writer, seed_mode_code(header.seed_mode))?;
    write_varint(writer, header.max_occ as u64)?;
    writer.write_all(&[header.hash_bits])?;
    if header.version >= 5 {
        writer.write_all(&[u8::from(header.exact_randstrobe_match)])?;
    }
    if header.version >= 6 {
        let (min, max) = header
            .strobe_start_dist
            .map(|distance| (distance.min, distance.max))
            .unwrap_or((0, 0));
        write_varint(writer, u64::from(min))?;
        write_varint(writer, u64::from(max))?;
    }
    write_varint(writer, header.created_unix_seconds)?;
    write_varint(writer, header.filtered_seed_hashes)?;
    Ok(())
}

fn read_serialized_index_metadata<R: Read>(reader: &mut R) -> Result<SerializedIndexMetadata> {
    let mut magic = [0u8; 8];
    reader.read_exact(&mut magic)?;
    let version = read_varint(reader)? as u32;
    validate_index_identity(magic, version)?;
    let metam_version = read_string(reader)?;
    let k = read_varint(reader)? as usize;
    let w = read_varint(reader)? as usize;
    let seed_mode = seed_mode_from_code(read_varint(reader)?)?;
    let max_occ = read_varint(reader)? as usize;
    let mut hash_bits_buf = [0u8; 1];
    reader.read_exact(&mut hash_bits_buf)?;
    let hash_bits = hash_bits_buf[0];
    validate_hash_bits(hash_bits)?;
    let exact_randstrobe_match = if version >= 5 {
        let mut value = [0u8; 1];
        reader.read_exact(&mut value)?;
        match value[0] {
            0 => false,
            1 => true,
            _ => bail!("invalid exact randstrobe match flag in index"),
        }
    } else {
        false
    };
    validate_exact_randstrobe_match(seed_mode, exact_randstrobe_match)?;
    let strobe_start_dist = if version >= 6 {
        let min = read_varint(reader)?;
        let max = read_varint(reader)?;
        if min == 0 && max == 0 {
            None
        } else {
            let min = u16::try_from(min).context("randstrobe start-distance min exceeds u16")?;
            let max = u16::try_from(max).context("randstrobe start-distance max exceeds u16")?;
            Some(StrobeStartDistance::new(min, max)?)
        }
    } else {
        None
    };
    validate_strobe_start_distance_metadata(version, seed_mode, strobe_start_dist)?;
    let created_unix_seconds = read_varint(reader)?;
    let filtered_seed_hashes = read_varint(reader)?;
    Ok(SerializedIndexMetadata {
        magic,
        version,
        metam_version,
        k,
        w,
        seed_mode,
        max_occ,
        hash_bits,
        exact_randstrobe_match,
        strobe_start_dist,
        created_unix_seconds,
        filtered_seed_hashes,
    })
}

fn read_index_payload<R: Read>(reader: &mut R) -> Result<Index> {
    let SerializedIndexMetadata {
        magic,
        version,
        metam_version,
        k,
        w,
        seed_mode,
        max_occ,
        hash_bits,
        exact_randstrobe_match,
        strobe_start_dist,
        created_unix_seconds,
        filtered_seed_hashes,
    } = read_serialized_index_metadata(reader)?;
    let ref_count = read_varint(reader)? as usize;
    let mut references = Vec::with_capacity(ref_count);
    for _ in 0..ref_count {
        let id = read_string(reader)?;
        let len = read_varint(reader)?;
        let packed_bases = read_bytes(reader)?;
        let ambiguous = BitMask::from_bytes(read_bytes(reader)?);
        references.push(Reference {
            id,
            len,
            packed_bases,
            ambiguous,
        });
    }
    let seed_table = read_seed_table(reader, hash_bits, exact_randstrobe_match)?;
    let total_reference_bases = references.iter().map(|r| r.len).sum();
    let reference_storage_bytes = references
        .iter()
        .map(|r| r.packed_bases.len() as u64 + r.ambiguous.byte_len() as u64)
        .sum();
    let indexed_seed_hashes = seed_table.group_count() as u64;
    let indexed_seed_occurrences = seed_table.occurrence_count() as u64;
    Ok(Index {
        header: IndexHeader {
            magic,
            version,
            metam_version,
            k,
            w,
            seed_mode,
            max_occ,
            hash_bits,
            exact_randstrobe_match,
            strobe_start_dist,
            reference_count: references.len(),
            total_reference_bases,
            indexed_seed_hashes,
            indexed_seed_occurrences,
            filtered_seed_hashes,
            reference_storage_bytes,
            created_unix_seconds,
        },
        references,
        seed_table,
    })
}

fn write_string<W: Write>(writer: &mut W, value: &str) -> Result<()> {
    write_bytes(writer, value.as_bytes())
}

fn read_string<R: Read>(reader: &mut R) -> Result<String> {
    Ok(String::from_utf8(read_bytes(reader)?)?)
}

fn write_bytes<W: Write>(writer: &mut W, bytes: &[u8]) -> Result<()> {
    write_varint(writer, bytes.len() as u64)?;
    writer.write_all(bytes)?;
    Ok(())
}

fn read_bytes<R: Read>(reader: &mut R) -> Result<Vec<u8>> {
    let len = read_varint(reader)? as usize;
    let mut bytes = vec![0u8; len];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn write_reference<W: Write>(writer: &mut W, reference: &Reference) -> Result<()> {
    write_string(writer, &reference.id)?;
    write_varint(writer, reference.len)?;
    write_bytes(writer, &reference.packed_bases)?;
    write_bytes(writer, &reference.ambiguous.clone().into_bytes())?;
    Ok(())
}

fn read_reference<R: Read>(reader: &mut R) -> Result<Reference> {
    Ok(Reference {
        id: read_string(reader)?,
        len: read_varint(reader)?,
        packed_bases: read_bytes(reader)?,
        ambiguous: BitMask::from_bytes(read_bytes(reader)?),
    })
}

fn read_references(path: &Path, count: usize) -> Result<Vec<Reference>> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut references = Vec::with_capacity(count);
    for _ in 0..count {
        references.push(read_reference(&mut reader)?);
    }
    Ok(references)
}

fn write_seed_table<W: Write>(writer: &mut W, seed_table: &SeedTable) -> Result<()> {
    write_varint(writer, seed_table.occurrence_count() as u64)?;
    write_varint(writer, seed_table.group_count() as u64)?;

    let mut previous_hash = 0u64;
    let mut start = 0usize;
    while start < seed_table.occurrence_count() {
        let hash = seed_table.hash_at(start);
        let mut end = start + 1;
        while end < seed_table.occurrence_count() && seed_table.hash_at(end) == hash {
            end += 1;
        }
        let count = end - start;
        write_occurrence_group_body(
            writer,
            &mut previous_hash,
            hash,
            count,
            seed_table.occurrences[start..end]
                .iter()
                .copied()
                .map(unpack_occurrence),
        )?;
        start = end;
    }
    writer.write_all(&seed_table.fingerprints)?;
    Ok(())
}

fn write_seed_group_body<W: Write>(
    writer: &mut W,
    previous_hash: &mut u64,
    hash: u64,
    seeds: &[RawSeed],
) -> Result<()> {
    write_occurrence_group_body(
        writer,
        previous_hash,
        hash,
        seeds.len(),
        seeds.iter().map(|seed| SeedOccurrence {
            ref_id: seed.ref_id,
            pos: seed.pos,
        }),
    )
}

fn write_seed_fingerprints<W: Write>(
    writer: &mut W,
    seeds: &[RawSeed],
    hash_bits: u8,
    scratch: &mut Vec<u8>,
) -> Result<()> {
    let fingerprint_bits = randstrobe_fingerprint_bits(hash_bits);
    let width = fingerprint_bytes(fingerprint_bits);
    if width == 0 {
        scratch.clear();
        return Ok(());
    }
    let required = seeds
        .len()
        .checked_mul(width)
        .context("seed fingerprint buffer length overflow")?;
    scratch.clear();
    if scratch.capacity() < required {
        scratch.reserve(required);
    }
    for seed in seeds {
        push_packed_fingerprint(
            scratch,
            randstrobe_fingerprint(seed.hash, hash_bits),
            fingerprint_bits,
        );
    }
    writer.write_all(scratch)?;
    Ok(())
}

fn write_occurrence_group_body<W, I>(
    writer: &mut W,
    previous_hash: &mut u64,
    hash: u64,
    count: usize,
    occurrences: I,
) -> Result<()>
where
    W: Write,
    I: IntoIterator<Item = SeedOccurrence>,
{
    write_varint(
        writer,
        hash.checked_sub(*previous_hash)
            .context("seed records are not sorted by hash")?,
    )?;
    write_varint(writer, count as u64)?;
    *previous_hash = hash;

    let mut previous_ref = 0u32;
    let mut previous_pos = 0u32;
    for occurrence in occurrences {
        let ref_delta = occurrence
            .ref_id
            .checked_sub(previous_ref)
            .context("seed records are not sorted by reference ID")?;
        write_varint(writer, u64::from(ref_delta))?;
        let pos_delta_or_abs = if ref_delta == 0 {
            occurrence
                .pos
                .checked_sub(previous_pos)
                .context("seed records are not sorted by position")?
        } else {
            occurrence.pos
        };
        write_varint(writer, u64::from(pos_delta_or_abs))?;
        previous_ref = occurrence.ref_id;
        previous_pos = occurrence.pos;
    }
    Ok(())
}

fn read_seed_table<R: Read>(
    reader: &mut R,
    hash_bits: u8,
    exact_randstrobe_match: bool,
) -> Result<SeedTable> {
    let seed_count = read_varint(reader)? as usize;
    let group_count = read_varint(reader)? as usize;
    let mut seed_table =
        SeedTable::new_with_exact_randstrobe_match(hash_bits, exact_randstrobe_match);
    if hash_bits <= 32 {
        seed_table.hashes32.reserve(seed_count);
    } else {
        seed_table.hashes64.reserve(seed_count);
    }
    seed_table.occurrences.reserve(seed_count);
    let mut previous_hash = 0u64;
    for _ in 0..group_count {
        let hash_delta = read_varint(reader)?;
        let hash = previous_hash
            .checked_add(hash_delta)
            .context("seed hash delta overflow in index")?;
        previous_hash = hash;
        let occurrence_count = read_varint(reader)? as usize;
        if occurrence_count > u32::MAX as usize {
            bail!("seed occurrence group exceeds index format limit");
        }
        seed_table.push_hash_copies(hash, occurrence_count)?;
        let mut previous_ref = 0u32;
        let mut previous_pos = 0u32;
        for _ in 0..occurrence_count {
            let ref_delta = read_varint(reader)?;
            let ref_id = u64::from(previous_ref)
                .checked_add(ref_delta)
                .context("reference ID delta overflow in index")?;
            if ref_id > u64::from(u32::MAX) {
                bail!("reference ID exceeds index format limit");
            }
            let pos_delta_or_abs = read_varint(reader)?;
            let pos = if ref_delta == 0 {
                u64::from(previous_pos)
                    .checked_add(pos_delta_or_abs)
                    .context("reference position delta overflow in index")?
            } else {
                pos_delta_or_abs
            };
            if pos > u64::from(u32::MAX) {
                bail!("reference position exceeds index format limit");
            }
            previous_ref = ref_id as u32;
            previous_pos = pos as u32;
            seed_table
                .occurrences
                .push(pack_occurrence(previous_ref, previous_pos));
        }
    }
    if seed_table.occurrence_count() != seed_count {
        bail!("index seed record count mismatch");
    }
    let fingerprint_bytes = seed_table
        .expected_fingerprint_bytes()
        .context("index fingerprint array length overflow")?;
    seed_table.fingerprints.resize(fingerprint_bytes, 0);
    reader.read_exact(&mut seed_table.fingerprints)?;
    seed_table.rebuild_lookup_index();
    Ok(seed_table)
}

fn write_varint<W: Write>(writer: &mut W, mut value: u64) -> Result<()> {
    while value >= 0x80 {
        writer.write_all(&[((value as u8) & 0x7f) | 0x80])?;
        value >>= 7;
    }
    writer.write_all(&[value as u8])?;
    Ok(())
}

fn read_varint<R: Read>(reader: &mut R) -> Result<u64> {
    let mut value = 0u64;
    let mut shift = 0u32;
    loop {
        let mut byte = [0u8; 1];
        reader.read_exact(&mut byte)?;
        value |= u64::from(byte[0] & 0x7f) << shift;
        if byte[0] & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
        if shift >= 64 {
            bail!("invalid varint in index");
        }
    }
}

fn seed_mode_code(mode: SeedMode) -> u64 {
    match mode {
        SeedMode::Minimizer => 0,
        SeedMode::Randstrobe => 1,
    }
}

fn seed_mode_from_code(code: u64) -> Result<SeedMode> {
    match code {
        0 => Ok(SeedMode::Minimizer),
        1 => Ok(SeedMode::Randstrobe),
        _ => bail!("invalid seed mode in index"),
    }
}

pub fn human_bytes(bytes: u64) -> String {
    let units = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0usize;
    while value >= 1024.0 && unit + 1 < units.len() {
        value /= 1024.0;
        unit += 1;
    }
    format!("{value:.2} {}", units[unit])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn colliding_randstrobe_hashes() -> (u64, u64) {
        let key = 0x1234_5678u64;
        let a = (key << 32) | (0x12_345u64 << RANDSTROBE_HASH_LOW_ZERO_BITS);
        let b = (key << 32) | (0x54_321u64 << RANDSTROBE_HASH_LOW_ZERO_BITS);
        assert_eq!(seed_key(a, 32), seed_key(b, 32));
        assert_ne!(randstrobe_fingerprint(a, 32), randstrobe_fingerprint(b, 32));
        (a, b)
    }

    fn collision_group() -> (u64, Vec<RawSeed>) {
        let (a, b) = colliding_randstrobe_hashes();
        (
            seed_key(a, 32),
            vec![
                RawSeed {
                    hash: a,
                    ref_id: 0,
                    pos: 10,
                },
                RawSeed {
                    hash: b,
                    ref_id: 1,
                    pos: 20,
                },
            ],
        )
    }

    fn test_index(seed_table: SeedTable, seed_mode: SeedMode) -> Index {
        let exact_randstrobe_match = seed_table.fingerprint_bits > 0;
        Index {
            header: IndexHeader {
                magic: INDEX_MAGIC,
                version: INDEX_VERSION,
                metam_version: "test".to_string(),
                k: 15,
                w: 25,
                seed_mode,
                max_occ: 50,
                hash_bits: seed_table.hash_bits,
                exact_randstrobe_match,
                strobe_start_dist: (seed_mode == SeedMode::Randstrobe)
                    .then_some(StrobeStartDistance::default()),
                reference_count: 0,
                total_reference_bases: 0,
                indexed_seed_hashes: seed_table.group_count() as u64,
                indexed_seed_occurrences: seed_table.occurrence_count() as u64,
                filtered_seed_hashes: 0,
                reference_storage_bytes: 0,
                created_unix_seconds: 0,
            },
            references: Vec::new(),
            seed_table,
        }
    }

    fn with_test_references(mut index: Index, count: usize) -> Index {
        index.references = (0..count)
            .map(|i| {
                let mut sequence = b"ACGT".repeat(32);
                if i % 2 == 1 {
                    sequence[17] = b'N';
                }
                Reference {
                    id: format!("ref-{i}"),
                    len: sequence.len() as u64,
                    packed_bases: pack_dna_2bit(&sequence),
                    ambiguous: BitMask::from_ambiguous(&sequence),
                }
            })
            .collect();
        index.header.reference_count = index.references.len();
        index.header.total_reference_bases = index.references.iter().map(|r| r.len).sum();
        index.header.reference_storage_bytes = index
            .references
            .iter()
            .map(|r| r.packed_bases.len() as u64 + r.ambiguous.byte_len() as u64)
            .sum();
        index
    }

    #[cfg(unix)]
    fn assert_streaming_cache_matches_owned_writer(index: &Index) {
        let dir = tempfile::TempDir::new().unwrap();
        let index_path = dir.path().join("source.mxi");
        let owned_cache_path = dir.path().join("owned.metamcache");
        let streaming_cache_path = dir.path().join("streaming.metamcache");
        write_index(index, index_path.to_str().unwrap()).unwrap();
        let source = index_source_fingerprint(index_path.to_str().unwrap()).unwrap();
        write_index_cache(owned_cache_path.to_str().unwrap(), source, index).unwrap();
        let mut progress_updates = Vec::new();
        write_index_cache_from_source_progress(
            streaming_cache_path.to_str().unwrap(),
            source,
            index_path.to_str().unwrap(),
            |progress| progress_updates.push(progress),
        )
        .unwrap();

        assert_eq!(
            progress_updates.first(),
            Some(&IndexCacheBuildProgress {
                compressed_bytes_read: 0,
                compressed_bytes_total: source.size,
            })
        );
        assert_eq!(
            progress_updates.last(),
            Some(&IndexCacheBuildProgress {
                compressed_bytes_read: source.size,
                compressed_bytes_total: source.size,
            })
        );
        assert!(progress_updates.windows(2).all(|updates| {
            updates[0].compressed_bytes_read <= updates[1].compressed_bytes_read
                && updates[0].compressed_bytes_total == updates[1].compressed_bytes_total
        }));

        assert_eq!(
            std::fs::read(&streaming_cache_path).unwrap(),
            std::fs::read(&owned_cache_path).unwrap()
        );
        let mapped = MappedIndex::open(streaming_cache_path.to_str().unwrap(), source).unwrap();
        assert_eq!(mapped.header.reference_count, index.header.reference_count);
        for (ref_id, expected) in index.references.iter().enumerate() {
            let found = mapped.reference(ref_id as u32);
            assert_eq!(found.id, expected.id);
            assert_eq!(found.len, expected.len);
            for pos in 0..expected.len() {
                assert_eq!(found.base_at(pos), expected.base_at(pos));
                assert_eq!(found.code_at(pos), expected.code_at(pos));
            }
        }
    }

    #[test]
    fn code_at_matches_decoding_base_at() {
        // Ambiguity is covered explicitly: an all-ACGT reference never reaches
        // the branch the comparison callers depend on.
        let record = Record {
            id: "r".to_string(),
            seq: b"ACGTNNACGTNacgtNNNNTTTTGGGGCCCCAAAAN".to_vec(),
        };
        let reference = Reference::from_record(&record);
        assert!(reference.sequence().contains(&b'N'));
        for pos in 0..reference.len() {
            assert_eq!(
                reference.code_at(pos),
                crate::sequence::base_code(reference.base_at(pos)),
                "position {pos}"
            );
        }
        assert_eq!(reference.code_at(0), Some(0));
        assert_eq!(reference.code_at(4), None);
    }

    #[test]
    fn exact_randstrobe_match_filters_same_key_different_pair() {
        let (key, seeds) = collision_group();
        let mut table = SeedTable::new_with_exact_randstrobe_match(32, true);
        table.push_group(key, &seeds).unwrap();
        table.rebuild_lookup_index();

        let hits_a = table.lookup(key, seeds[0].hash).iter().collect::<Vec<_>>();
        let hits_b = table.lookup(key, seeds[1].hash).iter().collect::<Vec<_>>();

        assert_eq!(hits_a, vec![SeedOccurrence { ref_id: 0, pos: 10 }]);
        assert_eq!(hits_b, vec![SeedOccurrence { ref_id: 1, pos: 20 }]);
    }

    #[test]
    fn relaxed_randstrobe_lookup_can_ignore_fingerprint() {
        let (key, seeds) = collision_group();
        let mut table = SeedTable::new_with_exact_randstrobe_match(32, true);
        table.push_group(key, &seeds).unwrap();
        table.rebuild_lookup_index();

        let strict_hits = table.lookup(key, seeds[0].hash).iter().collect::<Vec<_>>();
        let relaxed_hits = table
            .lookup_with_mode(
                key,
                seeds[0].hash,
                SeedLookupMode::IgnoreRandstrobeFingerprint,
            )
            .iter()
            .collect::<Vec<_>>();

        assert_eq!(strict_hits, vec![SeedOccurrence { ref_id: 0, pos: 10 }]);
        assert_eq!(
            relaxed_hits,
            vec![
                SeedOccurrence { ref_id: 0, pos: 10 },
                SeedOccurrence { ref_id: 1, pos: 20 },
            ]
        );
    }

    #[test]
    fn legacy_randstrobe_lookup_keeps_all_same_key_occurrences() {
        let (key, seeds) = collision_group();
        let mut table = SeedTable::new(32);
        table.push_group(key, &seeds).unwrap();
        table.rebuild_lookup_index();

        let hits = table.lookup(key, seeds[0].hash).iter().collect::<Vec<_>>();
        assert_eq!(
            hits,
            vec![
                SeedOccurrence { ref_id: 0, pos: 10 },
                SeedOccurrence { ref_id: 1, pos: 20 },
            ]
        );
    }

    #[test]
    fn full_hash_span_needs_no_secondary_fingerprint() {
        let hash_bits = RANDSTROBE_HASH_SPAN_BITS;
        let a = 0xFEDC_BA98_7654_3200u64 & crate::seeds::RANDSTROBE_HASH_MASK;
        let b = a ^ (1u64 << RANDSTROBE_HASH_LOW_ZERO_BITS);
        let key_a = seed_key(a, hash_bits);
        let key_b = seed_key(b, hash_bits);
        assert_ne!(key_a, key_b);

        let mut rows = vec![
            (
                key_a,
                RawSeed {
                    hash: a,
                    ref_id: 0,
                    pos: 10,
                },
            ),
            (
                key_b,
                RawSeed {
                    hash: b,
                    ref_id: 1,
                    pos: 20,
                },
            ),
        ];
        rows.sort_unstable_by_key(|(key, _)| *key);

        let mut table = SeedTable::new_with_exact_randstrobe_match(hash_bits, true);
        for (key, seed) in rows {
            table.push_group(key, &[seed]).unwrap();
        }
        table.rebuild_lookup_index();

        assert_eq!(table.fingerprint_bits, 0);
        assert!(table.fingerprints.is_empty());
        assert_eq!(table.lookup(key_a, a).len(), 1);
        assert_eq!(table.lookup(key_b, b).len(), 1);
    }

    #[test]
    fn invalid_hash_width_is_rejected_before_seed_table_allocation() {
        let mut encoded = Vec::new();
        encoded.extend_from_slice(&INDEX_MAGIC);
        write_varint(&mut encoded, u64::from(INDEX_VERSION)).unwrap();
        write_string(&mut encoded, "test").unwrap();
        write_varint(&mut encoded, 15).unwrap();
        write_varint(&mut encoded, 25).unwrap();
        write_varint(&mut encoded, seed_mode_code(SeedMode::Randstrobe)).unwrap();
        write_varint(&mut encoded, 50).unwrap();
        encoded.push(65);

        let error = read_index_payload(&mut encoded.as_slice()).unwrap_err();
        assert!(error.to_string().contains("invalid hash_bits"));
    }

    #[test]
    fn parallel_reference_seeding_matches_single_thread_index() {
        use std::io::Write as _;

        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        let seq1 = "ACGTGCAATGCCGTAACCGGTTAACCGGTTACGATCGATGCTAGCTAGGCTA".repeat(40);
        let seq2 = "TTGACCGTACGATGATCCGATCGGATCCGTAACGTTAGCACTGATCGATCGA".repeat(40);
        let seq3 = "GGGATCGTACCATGACCTTAGGCTAATCGGATCGTACCGATGCTAGGATCCA".repeat(40);
        writeln!(fasta, ">ref1\n{seq1}\n>ref2\n{seq2}\n>ref3\n{seq3}").unwrap();
        let path = fasta.path().to_str().unwrap();

        let build = |threads| {
            build_index_with_progress(
                path,
                15,
                25,
                10_000,
                32,
                SeedMode::Randstrobe,
                true,
                Some(StrobeStartDistance::default()),
                threads,
                DEFAULT_INDEX_SEED_CHUNK_BYTES,
                |_| {},
            )
            .unwrap()
        };
        let (single, single_stats) = build(1);
        // Sixteen workers request more temporary hash partitions than the
        // single-threaded build on ordinary file-descriptor limits. The
        // decoded index must remain identical regardless of that partitioning.
        let (parallel, parallel_stats) = build(16);

        assert_eq!(
            single_stats.raw_seed_occurrences,
            parallel_stats.raw_seed_occurrences
        );
        assert_eq!(single_stats.filtered_hashes, parallel_stats.filtered_hashes);
        assert_eq!(
            single_stats.max_partition_seed_occurrences,
            parallel_stats.max_partition_seed_occurrences
        );

        assert_eq!(
            single.header.reference_count,
            parallel.header.reference_count
        );
        assert_eq!(
            single.header.total_reference_bases,
            parallel.header.total_reference_bases
        );
        assert_eq!(
            single.header.indexed_seed_hashes,
            parallel.header.indexed_seed_hashes
        );
        assert_eq!(
            single.header.indexed_seed_occurrences,
            parallel.header.indexed_seed_occurrences
        );
        assert_eq!(
            single.header.filtered_seed_hashes,
            parallel.header.filtered_seed_hashes
        );
        assert_eq!(
            single.header.reference_storage_bytes,
            parallel.header.reference_storage_bytes
        );

        assert_eq!(single.references.len(), parallel.references.len());
        for (a, b) in single.references.iter().zip(parallel.references.iter()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.len, b.len);
            assert_eq!(a.packed_bases, b.packed_bases);
            assert_eq!(a.sequence(), b.sequence());
        }

        assert_eq!(single.seed_table.hash_bits, parallel.seed_table.hash_bits);
        assert_eq!(single.seed_table.hashes32, parallel.seed_table.hashes32);
        assert_eq!(single.seed_table.hashes64, parallel.seed_table.hashes64);
        assert_eq!(
            single.seed_table.occurrences,
            parallel.seed_table.occurrences
        );
        assert_eq!(
            single.seed_table.fingerprint_bits,
            parallel.seed_table.fingerprint_bits
        );
        assert_eq!(
            single.seed_table.fingerprints,
            parallel.seed_table.fingerprints
        );
        assert_eq!(
            single.seed_table.lookup_index.prefix_bits,
            parallel.seed_table.lookup_index.prefix_bits
        );
        assert_eq!(
            single.seed_table.lookup_index.offsets,
            parallel.seed_table.lookup_index.offsets
        );
    }

    #[test]
    fn temporary_partition_writers_respect_process_limits_and_hard_cap() {
        assert_eq!(temporary_seed_partition_count_for_limit(1, usize::MAX), 256);
        assert_eq!(
            temporary_seed_partition_count_for_limit(40, usize::MAX),
            MAX_OPEN_PARTITION_WRITERS
        );
        assert_eq!(
            temporary_seed_partition_count_for_limit(40, 256),
            256 - OPEN_FILE_RESERVE
        );
        assert_eq!(temporary_seed_partition_count_for_limit(40, 32), 24);
    }

    #[test]
    fn index_workspace_is_created_under_requested_parent() {
        let parent = tempfile::tempdir().unwrap();
        let workspace = make_temp_dir_in(Some(parent.path())).unwrap();
        assert_eq!(workspace.parent(), Some(parent.path()));
        assert!(workspace.is_dir());
        drop(TempDirCleanup::new(workspace.clone()));
        assert!(!workspace.exists());
    }

    #[test]
    fn failed_atomic_index_write_preserves_existing_index() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("index.mxi");
        std::fs::write(&path, b"existing-index").unwrap();

        let error = write_compressed_index_atomic(&path, |encoder| {
            encoder.write_all(b"incomplete replacement")?;
            bail!("injected payload failure")
        })
        .unwrap_err();

        assert!(error.to_string().contains("injected payload failure"));
        assert_eq!(std::fs::read(&path).unwrap(), b"existing-index");
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
    }

    #[test]
    fn version_six_roundtrip_preserves_randstrobe_fingerprints() {
        let (key, seeds) = collision_group();
        let mut seed_table = SeedTable::new_with_exact_randstrobe_match(32, true);
        seed_table.push_group(key, &seeds).unwrap();
        seed_table.rebuild_lookup_index();
        let index = Index {
            header: IndexHeader {
                magic: INDEX_MAGIC,
                version: INDEX_VERSION,
                metam_version: "test".to_string(),
                k: 15,
                w: 25,
                seed_mode: SeedMode::Randstrobe,
                max_occ: 50,
                hash_bits: 32,
                exact_randstrobe_match: true,
                strobe_start_dist: Some(StrobeStartDistance::default()),
                reference_count: 0,
                total_reference_bases: 0,
                indexed_seed_hashes: 1,
                indexed_seed_occurrences: 2,
                filtered_seed_hashes: 0,
                reference_storage_bytes: 0,
                created_unix_seconds: 0,
            },
            references: Vec::new(),
            seed_table,
        };

        let mut encoded = Vec::new();
        write_index_payload(&mut encoded, &index).unwrap();
        let decoded = read_index_payload(&mut encoded.as_slice()).unwrap();

        assert!(decoded.header.exact_randstrobe_match);
        assert_eq!(
            decoded.header.strobe_start_dist,
            Some(StrobeStartDistance::default())
        );
        assert_eq!(decoded.lookup(seeds[0].hash).len(), 1);
        assert_eq!(decoded.lookup(seeds[1].hash).len(), 1);
        assert_eq!(
            decoded
                .lookup_with_mode(seeds[0].hash, SeedLookupMode::IgnoreRandstrobeFingerprint)
                .len(),
            2
        );
    }

    #[test]
    fn version_five_index_remains_readable_with_legacy_strobe_distance() {
        let (key, seeds) = collision_group();
        let mut seed_table = SeedTable::new_with_exact_randstrobe_match(32, true);
        seed_table.push_group(key, &seeds).unwrap();
        seed_table.rebuild_lookup_index();
        let index = Index {
            header: IndexHeader {
                magic: INDEX_MAGIC,
                version: 5,
                metam_version: "test".to_string(),
                k: 15,
                w: 25,
                seed_mode: SeedMode::Randstrobe,
                max_occ: 50,
                hash_bits: 32,
                exact_randstrobe_match: true,
                strobe_start_dist: None,
                reference_count: 0,
                total_reference_bases: 0,
                indexed_seed_hashes: 1,
                indexed_seed_occurrences: 2,
                filtered_seed_hashes: 0,
                reference_storage_bytes: 0,
                created_unix_seconds: 0,
            },
            references: Vec::new(),
            seed_table,
        };

        let mut encoded = Vec::new();
        write_index_payload(&mut encoded, &index).unwrap();
        let decoded = read_index_payload(&mut encoded.as_slice()).unwrap();

        assert_eq!(decoded.header.version, 5);
        assert!(decoded.header.exact_randstrobe_match);
        assert_eq!(decoded.header.strobe_start_dist, None);
        assert_eq!(decoded.lookup(seeds[0].hash).len(), 1);
        #[cfg(unix)]
        assert_streaming_cache_matches_owned_writer(&index);
    }

    #[cfg(unix)]
    #[test]
    fn streaming_cache_matches_owned_grouped_randstrobe_cache() {
        let (key, mut seeds) = collision_group();
        let second_hash = seeds[1].hash;
        seeds.push(RawSeed {
            hash: second_hash,
            ref_id: 2,
            pos: 30,
        });
        seeds.push(RawSeed {
            hash: second_hash,
            ref_id: 3,
            pos: 40,
        });
        let mut seed_table = SeedTable::new_with_exact_randstrobe_match(32, true);
        seed_table.push_group(key, &seeds).unwrap();
        seed_table.rebuild_lookup_index();
        let index = with_test_references(test_index(seed_table, SeedMode::Randstrobe), 4);

        assert_eq!(
            choose_cache_seed_key_layout(
                index.header.hash_bits,
                index.seed_table.occurrence_count(),
                index.seed_table.group_count(),
            ),
            CacheSeedKeyLayout::Grouped
        );
        assert_streaming_cache_matches_owned_writer(&index);
    }

    #[cfg(unix)]
    #[test]
    fn streaming_cache_matches_owned_repeated_wide_key_cache() {
        let key = 0x1234_5678_9abc_def0u64;
        let occurrences = [
            SeedOccurrence { ref_id: 0, pos: 10 },
            SeedOccurrence { ref_id: 1, pos: 20 },
        ];
        let mut seed_table = SeedTable::new(64);
        seed_table.push_occurrence_group(key, &occurrences).unwrap();
        seed_table.rebuild_lookup_index();
        let index = with_test_references(test_index(seed_table, SeedMode::Minimizer), 2);

        assert_eq!(
            choose_cache_seed_key_layout(
                index.header.hash_bits,
                index.seed_table.occurrence_count(),
                index.seed_table.group_count(),
            ),
            CacheSeedKeyLayout::Repeated
        );
        assert_streaming_cache_matches_owned_writer(&index);
    }

    #[cfg(unix)]
    #[test]
    fn streaming_cache_matches_owned_empty_seed_cache() {
        let mut seed_table = SeedTable::new(32);
        seed_table.rebuild_lookup_index();
        let index = with_test_references(test_index(seed_table, SeedMode::Minimizer), 1);
        assert_streaming_cache_matches_owned_writer(&index);
    }

    #[cfg(unix)]
    #[test]
    fn mapped_cache_relaxed_lookup_can_ignore_fingerprint() {
        let (key, mut seeds) = collision_group();
        let second_hash = seeds[1].hash;
        seeds.push(RawSeed {
            hash: second_hash,
            ref_id: 2,
            pos: 30,
        });
        seeds.push(RawSeed {
            hash: second_hash,
            ref_id: 3,
            pos: 40,
        });
        let mut seed_table = SeedTable::new_with_exact_randstrobe_match(32, true);
        seed_table.push_group(key, &seeds).unwrap();
        seed_table.rebuild_lookup_index();
        let index = test_index(seed_table, SeedMode::Randstrobe);
        let source = IndexSourceFingerprint {
            size: 123,
            modified_secs: 456,
            modified_nanos: 789,
        };
        let dir = tempfile::TempDir::new().unwrap();
        let cache_path = dir.path().join("index.metamcache");
        write_index_cache(cache_path.to_str().unwrap(), source, &index).unwrap();
        let mapped = MappedIndex::open(cache_path.to_str().unwrap(), source).unwrap();

        assert_eq!(
            mapped.header.strobe_start_dist,
            Some(StrobeStartDistance::default())
        );
        assert_eq!(mapped.seed_table.key_layout, CacheSeedKeyLayout::Grouped);
        assert_eq!(mapped.seed_table.hash_count(), 1);
        assert_eq!(mapped.seed_table.group_offset_width, 4);
        assert_eq!(mapped.seed_table.group_offset_len, 2);
        assert_eq!(mapped.lookup(seeds[0].hash).len(), 1);
        assert_eq!(
            mapped
                .lookup_with_mode(seeds[0].hash, SeedLookupMode::IgnoreRandstrobeFingerprint)
                .len(),
            4
        );
    }

    #[cfg(unix)]
    #[test]
    fn mapped_cache_keeps_repeated_keys_when_grouping_would_be_larger() {
        let key = 0x1234_5678u64;
        let occurrences = [
            SeedOccurrence { ref_id: 0, pos: 10 },
            SeedOccurrence { ref_id: 1, pos: 20 },
        ];
        let mut seed_table = SeedTable::new(32);
        seed_table.push_occurrence_group(key, &occurrences).unwrap();
        seed_table.rebuild_lookup_index();
        let index = test_index(seed_table, SeedMode::Minimizer);
        let source = IndexSourceFingerprint {
            size: 1,
            modified_secs: 2,
            modified_nanos: 3,
        };
        let dir = tempfile::TempDir::new().unwrap();
        let cache_path = dir.path().join("sparse.metamcache");
        write_index_cache(cache_path.to_str().unwrap(), source, &index).unwrap();
        let mapped = MappedIndex::open(cache_path.to_str().unwrap(), source).unwrap();

        assert_eq!(mapped.seed_table.key_layout, CacheSeedKeyLayout::Repeated);
        assert_eq!(mapped.seed_table.hash_count(), occurrences.len());
        assert_eq!(mapped.seed_table.group_offset_width, 0);
        assert_eq!(mapped.seed_table.group_offset_len, 0);
        assert_eq!(
            mapped.lookup(key << 32).iter().collect::<Vec<_>>(),
            occurrences
        );
    }

    #[cfg(unix)]
    #[test]
    fn mapped_cache_compacts_32_bit_keys_across_prefix_buckets() {
        let keys = [0x0000_0001u64, 0x1234_5678, 0xffff_fffe];
        let mut seed_table = SeedTable::new(32);
        let mut expected = Vec::new();
        for (group, key) in keys.into_iter().enumerate() {
            let occurrences = (0..4)
                .map(|i| SeedOccurrence {
                    ref_id: (group * 4 + i) as u32,
                    pos: (100 + group * 10 + i) as u32,
                })
                .collect::<Vec<_>>();
            seed_table.push_occurrence_group(key, &occurrences).unwrap();
            expected.push((key, occurrences));
        }
        seed_table.rebuild_lookup_index();
        let index = test_index(seed_table, SeedMode::Minimizer);
        let source = IndexSourceFingerprint {
            size: 4,
            modified_secs: 5,
            modified_nanos: 6,
        };
        let dir = tempfile::TempDir::new().unwrap();
        let cache_path = dir.path().join("compact-32.metamcache");
        write_index_cache(cache_path.to_str().unwrap(), source, &index).unwrap();
        let mapped = MappedIndex::open(cache_path.to_str().unwrap(), source).unwrap();

        assert_eq!(mapped.seed_table.key_layout, CacheSeedKeyLayout::Grouped);
        assert_eq!(mapped.seed_table.hash_count(), expected.len());
        assert_eq!(mapped.seed_table.group_offset_len, expected.len() + 1);
        for (key, occurrences) in expected {
            assert_eq!(
                mapped.lookup(key << 32).iter().collect::<Vec<_>>(),
                occurrences
            );
        }
        assert!(mapped.lookup(0x1234_5679u64 << 32).is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn mapped_cache_compacts_64_bit_keys_without_changing_occurrences() {
        let first_key = 0x0102_0304_0506_0708u64;
        let second_key = 0x1112_1314_1516_1718u64;
        let first_occurrences = [
            SeedOccurrence { ref_id: 0, pos: 10 },
            SeedOccurrence { ref_id: 1, pos: 20 },
            SeedOccurrence { ref_id: 2, pos: 30 },
            SeedOccurrence { ref_id: 3, pos: 40 },
        ];
        let second_occurrences = [
            SeedOccurrence { ref_id: 4, pos: 50 },
            SeedOccurrence { ref_id: 5, pos: 60 },
            SeedOccurrence { ref_id: 6, pos: 70 },
        ];
        let mut seed_table = SeedTable::new(64);
        seed_table
            .push_occurrence_group(first_key, &first_occurrences)
            .unwrap();
        seed_table
            .push_occurrence_group(second_key, &second_occurrences)
            .unwrap();
        seed_table.rebuild_lookup_index();
        let index = test_index(seed_table, SeedMode::Minimizer);
        let source = IndexSourceFingerprint {
            size: 10,
            modified_secs: 20,
            modified_nanos: 30,
        };
        let dir = tempfile::TempDir::new().unwrap();
        let cache_path = dir.path().join("wide.metamcache");
        write_index_cache(cache_path.to_str().unwrap(), source, &index).unwrap();
        let mapped = MappedIndex::open(cache_path.to_str().unwrap(), source).unwrap();

        assert_eq!(mapped.seed_table.key_layout, CacheSeedKeyLayout::Grouped);
        assert_eq!(mapped.seed_table.hash_count(), 2);
        assert_eq!(mapped.seed_table.group_offset_width, 4);
        assert_eq!(mapped.seed_table.group_offset_len, 3);
        assert_eq!(
            mapped.lookup(first_key).iter().collect::<Vec<_>>(),
            first_occurrences
        );
        assert_eq!(
            mapped.lookup(second_key).iter().collect::<Vec<_>>(),
            second_occurrences
        );
        assert!(mapped.lookup(first_key + 1).is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn mapped_cache_handles_an_empty_seed_table() {
        let mut seed_table = SeedTable::new(32);
        seed_table.rebuild_lookup_index();
        let index = test_index(seed_table, SeedMode::Minimizer);
        let source = IndexSourceFingerprint {
            size: 7,
            modified_secs: 8,
            modified_nanos: 9,
        };
        let dir = tempfile::TempDir::new().unwrap();
        let cache_path = dir.path().join("empty.metamcache");
        write_index_cache(cache_path.to_str().unwrap(), source, &index).unwrap();
        let mapped = MappedIndex::open(cache_path.to_str().unwrap(), source).unwrap();

        assert_eq!(mapped.seed_table.key_layout, CacheSeedKeyLayout::Repeated);
        assert_eq!(mapped.seed_table.hash_count(), 0);
        assert_eq!(mapped.seed_table.offset_len, 0);
        assert_eq!(mapped.seed_table.group_offset_len, 0);
        assert!(mapped.lookup(0).is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn stale_cache_version_is_rebuilt_without_rewriting_the_index() {
        let key = 0x1234_5678u64;
        let occurrences = [
            SeedOccurrence { ref_id: 0, pos: 10 },
            SeedOccurrence { ref_id: 1, pos: 20 },
            SeedOccurrence { ref_id: 2, pos: 30 },
            SeedOccurrence { ref_id: 3, pos: 40 },
        ];
        let mut seed_table = SeedTable::new(32);
        seed_table.push_occurrence_group(key, &occurrences).unwrap();
        seed_table.rebuild_lookup_index();
        let index = test_index(seed_table, SeedMode::Minimizer);
        let dir = tempfile::TempDir::new().unwrap();
        let index_path = dir.path().join("index.mxi");
        let cache_path = dir.path().join("index.mxi.metamcache");
        write_index(&index, index_path.to_str().unwrap()).unwrap();
        let original_index = std::fs::read(&index_path).unwrap();

        let mut old_cache = Vec::new();
        old_cache.extend_from_slice(&INDEX_CACHE_MAGIC);
        old_cache.extend_from_slice(&(INDEX_CACHE_VERSION - 1).to_le_bytes());
        std::fs::write(&cache_path, old_cache).unwrap();

        let loaded = read_index_with_cache(
            index_path.to_str().unwrap(),
            Some(cache_path.to_str().unwrap()),
        )
        .unwrap();
        assert!(matches!(loaded, LoadedIndex::Mapped(_)));
        assert_eq!(std::fs::read(&index_path).unwrap(), original_index);
        let rebuilt_cache = std::fs::read(&cache_path).unwrap();
        assert_eq!(&rebuilt_cache[..8], &INDEX_CACHE_MAGIC);
        assert_eq!(
            u32::from_le_bytes(rebuilt_cache[8..12].try_into().unwrap()),
            INDEX_CACHE_VERSION
        );
    }

    #[test]
    fn version_four_index_remains_readable() {
        let (key, seeds) = collision_group();
        let mut seed_table = SeedTable::new(32);
        seed_table.push_group(key, &seeds).unwrap();
        seed_table.rebuild_lookup_index();
        let index = Index {
            header: IndexHeader {
                magic: INDEX_MAGIC,
                version: 4,
                metam_version: "test".to_string(),
                k: 15,
                w: 25,
                seed_mode: SeedMode::Randstrobe,
                max_occ: 50,
                hash_bits: 32,
                exact_randstrobe_match: false,
                strobe_start_dist: None,
                reference_count: 0,
                total_reference_bases: 0,
                indexed_seed_hashes: 1,
                indexed_seed_occurrences: 2,
                filtered_seed_hashes: 0,
                reference_storage_bytes: 0,
                created_unix_seconds: 0,
            },
            references: Vec::new(),
            seed_table,
        };

        let mut encoded = Vec::new();
        write_index_payload(&mut encoded, &index).unwrap();
        let decoded = read_index_payload(&mut encoded.as_slice()).unwrap();

        assert_eq!(decoded.header.version, 4);
        assert!(!decoded.header.exact_randstrobe_match);
        assert_eq!(decoded.header.strobe_start_dist, None);
        assert_eq!(decoded.lookup(seeds[0].hash).len(), 2);
        #[cfg(unix)]
        assert_streaming_cache_matches_owned_writer(&index);
    }
}

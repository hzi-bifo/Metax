use serde::{Deserialize, Serialize};

pub fn base_code(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

pub fn code_base(code: u8) -> u8 {
    match code & 3 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

pub fn normalize_base(base: u8) -> u8 {
    match base_code(base) {
        Some(code) => code_base(code),
        None => b'N',
    }
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut output = Vec::with_capacity(seq.len());
    revcomp_into(seq, &mut output);
    output
}

/// Write the reverse complement into reusable caller-owned storage.
pub(crate) fn revcomp_into(seq: &[u8], output: &mut Vec<u8>) {
    output.clear();
    output.reserve(seq.len());
    output.extend(seq.iter().rev().map(|b| match base_code(*b) {
        Some(0) => b'T',
        Some(1) => b'G',
        Some(2) => b'C',
        Some(3) => b'A',
        _ => b'N',
    }));
}

pub fn pack_dna_2bit(seq: &[u8]) -> Vec<u8> {
    let mut packed = vec![0u8; seq.len().div_ceil(4)];
    for (i, &base) in seq.iter().enumerate() {
        let code = base_code(base).unwrap_or(0);
        packed[i / 4] |= code << (6 - 2 * (i % 4));
    }
    packed
}

pub fn unpack_dna_2bit(packed: &[u8], len: usize, ambiguous: &BitMask) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        if ambiguous.get(i) {
            out.push(b'N');
        } else {
            let code = (packed[i / 4] >> (6 - 2 * (i % 4))) & 3;
            out.push(code_base(code));
        }
    }
    out
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BitMask {
    bits: Vec<u8>,
}

impl BitMask {
    pub fn from_ambiguous(seq: &[u8]) -> Self {
        if seq.iter().all(|base| base_code(*base).is_some()) {
            return Self { bits: Vec::new() };
        }
        let mut mask = Self {
            bits: vec![0; seq.len().div_ceil(8)],
        };
        for (i, base) in seq.iter().enumerate() {
            if base_code(*base).is_none() {
                mask.set(i);
            }
        }
        mask
    }

    pub fn set(&mut self, i: usize) {
        self.bits[i / 8] |= 1 << (i % 8);
    }

    pub fn get(&self, i: usize) -> bool {
        self.bits
            .get(i / 8)
            .map(|byte| (byte & (1 << (i % 8))) != 0)
            .unwrap_or(false)
    }

    pub fn byte_len(&self) -> usize {
        self.bits.len()
    }

    pub fn from_bytes(bits: Vec<u8>) -> Self {
        Self { bits }
    }

    pub fn into_bytes(self) -> Vec<u8> {
        self.bits
    }
}

pub fn projected_overlap(
    query_len: usize,
    ref_len: usize,
    diag: i64,
) -> Option<(usize, usize, usize, usize)> {
    let q_start = if diag < 0 { (-diag) as usize } else { 0 };
    let r_start = if diag > 0 { diag as usize } else { 0 };
    if q_start >= query_len || r_start >= ref_len {
        return None;
    }
    let len = (query_len - q_start).min(ref_len - r_start);
    if len == 0 {
        None
    } else {
        Some((q_start, q_start + len, r_start, r_start + len))
    }
}

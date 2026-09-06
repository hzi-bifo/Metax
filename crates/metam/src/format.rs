use crate::mapper::Mapping;
use std::fmt;

/// Profiler-facing mapping identities are intentionally canonicalized to five
/// decimal places. The mapper continues to use its full-precision identity for
/// candidate filtering and ranking; this boundary representation keeps saved
/// and structured mappings semantically identical while remaining compact.
pub const MAPPING_IDENTITY_SCALE: u32 = 100_000;

#[derive(Clone, Copy, Debug)]
pub struct MappingIdentity(f64);

pub fn mapping_identity(identity: f64) -> MappingIdentity {
    MappingIdentity(identity)
}

pub fn canonical_mapping_identity(identity: f64) -> f64 {
    if !identity.is_finite() {
        return identity;
    }
    let units = identity_units(identity);
    f64::from(units) / f64::from(MAPPING_IDENTITY_SCALE)
}

fn identity_units(identity: f64) -> u32 {
    (identity.clamp(0.0, 1.0) * f64::from(MAPPING_IDENTITY_SCALE)).round() as u32
}

impl fmt::Display for MappingIdentity {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if !self.0.is_finite() {
            return self.0.fmt(formatter);
        }
        let units = identity_units(self.0);
        let whole = units / MAPPING_IDENTITY_SCALE;
        let mut fractional = units % MAPPING_IDENTITY_SCALE;
        if fractional == 0 {
            return write!(formatter, "{whole}");
        }

        let mut width = 5usize;
        while fractional.is_multiple_of(10) {
            fractional /= 10;
            width -= 1;
        }
        write!(formatter, "{whole}.{fractional:0width$}")
    }
}

pub fn format_mapping(mapping: &Mapping, with_strand: bool, debug_output: bool) -> String {
    let mut out = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        mapping.read_id,
        mapping.ref_id,
        mapping.ref_start,
        mapping.ref_end,
        mapping_identity(mapping.identity),
        mapping.read_len,
        mapping.mapped_bases,
        mapping.mate_ref_status,
        mapping.mate_number,
        mapping.pair_id,
        mapping.candidate_boundary.bits()
    );
    if with_strand {
        out.push('\t');
        out.push(if mapping.strand == 1 { '+' } else { '-' });
    }
    if debug_output {
        out.push_str(&format!(
            "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            mapping.seed_count,
            mapping.score,
            mapping.second_score,
            mapping.rank,
            mapping.num_mappings,
            mapping.matched_bases,
            mapping.query_start,
            mapping.query_end,
            mapping.pair_score,
            mapping.adjusted_score,
            mapping.template_length,
            mapping.pair_status.as_str()
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_identity_is_compact_and_five_decimal() {
        assert_eq!(mapping_identity(1.0).to_string(), "1");
        assert_eq!(mapping_identity(0.98).to_string(), "0.98");
        assert_eq!(mapping_identity(0.993_333_333_333).to_string(), "0.99333");
        assert_eq!(mapping_identity(0.950_004).to_string(), "0.95");
        assert_eq!(mapping_identity(0.950_006).to_string(), "0.95001");
    }

    #[test]
    fn canonical_identity_matches_its_serialized_value() {
        for identity in [0.0, 0.8, 0.949_994, 0.949_996, 0.993_333_333, 1.0] {
            let serialized = mapping_identity(identity).to_string();
            let parsed: f64 = serialized.parse().unwrap();
            assert_eq!(
                parsed.to_bits(),
                canonical_mapping_identity(identity).to_bits()
            );
        }
    }
}

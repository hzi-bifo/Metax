//! Shared alignment-format helpers used by both profiler implementations.

/// Compute (identity, matched_len, fraction) from a SAM CIGAR string.
///
/// Primary implementation with an explicit switch.
///
/// - If `gap_compressed == true`:
///   - extended CIGAR (=,X present): identity = matches / (matches + mismatches + gap_events)
///   - non-extended (only M):        identity proxy = M / (M + gap_events)
///
///   where gap_events counts contiguous runs of I and D as 1 each.
///
/// - If `gap_compressed == false`:
///   - extended CIGAR (=,X present): identity = matches / (matches + mismatches + I_bases + D_bases)
///   - non-extended (only M):        identity proxy = M / (M + I_bases + D_bases)
///
/// `fraction` (query coverage) is always:
///     fraction = (M + = + X + I) / (M + = + X + I + S)
///
/// Returned `matched_len`:
///     * extended CIGAR: '='
///     * otherwise:      'M' (proxy)
pub(crate) fn cigar_to_identity_with_opts(cigar: &str, gap_compressed: bool) -> (f64, usize, f64) {
    const COUNT_N_AS_DELETION: bool = false;

    // Length counters
    let mut eq: usize = 0;
    let mut x: usize = 0;
    let mut i_: usize = 0;
    let mut d: usize = 0;
    let mut s: usize = 0;
    let mut m: usize = 0;

    // Gap-event counters (contiguous runs)
    let mut ins_events: usize = 0;
    let mut del_events: usize = 0;

    // Parsing state
    let mut count_buffer: usize = 0;
    let mut saw_any_op: bool = false;

    // Track whether we're inside a contiguous insertion/deletion run
    let mut in_ins_run: bool = false;
    let mut in_del_run: bool = false;

    for &byte in cigar.as_bytes() {
        if byte.is_ascii_digit() {
            count_buffer = count_buffer
                .saturating_mul(10)
                .saturating_add((byte - b'0') as usize);
            continue;
        }

        if count_buffer == 0 {
            // malformed token (e.g. "M") or "0M"; skip
            continue;
        }

        saw_any_op = true;

        match byte {
            b'=' => {
                eq += count_buffer;
                in_ins_run = false;
                in_del_run = false;
            }
            b'X' => {
                x += count_buffer;
                in_ins_run = false;
                in_del_run = false;
            }
            b'M' => {
                m += count_buffer;
                in_ins_run = false;
                in_del_run = false;
            }
            b'I' => {
                i_ += count_buffer;
                if !in_ins_run {
                    ins_events += 1;
                    in_ins_run = true;
                }
                in_del_run = false;
            }
            b'D' => {
                d += count_buffer;
                if !in_del_run {
                    del_events += 1;
                    in_del_run = true;
                }
                in_ins_run = false;
            }
            b'S' => {
                s += count_buffer;
                in_ins_run = false;
                in_del_run = false;
            }
            b'N' if COUNT_N_AS_DELETION => {
                d += count_buffer;
                if !in_del_run {
                    del_events += 1;
                    in_del_run = true;
                }
                in_ins_run = false;
            }
            _ => {
                // Ignore H, P, N (if not counted), etc.
                // Break gap runs because these terminate/interrupt the aligned segment.
                in_ins_run = false;
                in_del_run = false;
            }
        }

        count_buffer = 0;
    }

    if !saw_any_op {
        return (0.0, 0, 0.0);
    }

    // Fraction (query coverage)
    let aligned_query = m + eq + x + i_;
    let query_len = aligned_query + s;
    let fraction = if query_len == 0 {
        0.0
    } else {
        aligned_query as f64 / query_len as f64
    };

    // Identity
    let has_extended = (eq + x) > 0;
    let matched_len: usize = if has_extended { eq } else { m };

    let identity_den: usize = if gap_compressed {
        let gap_events = ins_events + del_events;
        if has_extended {
            eq + x + gap_events
        } else {
            m + gap_events
        }
    } else {
        if has_extended {
            eq + x + i_ + d
        } else {
            m + i_ + d
        }
    };

    let identity = if identity_den == 0 {
        0.0
    } else {
        matched_len as f64 / identity_den as f64
    };

    (identity, matched_len, fraction)
}

/// Compute the gap-compressed identity used by legacy SAM/BAM profiling.
pub(crate) fn cigar_to_identity(cigar: &str) -> (f64, usize, f64) {
    cigar_to_identity_with_opts(cigar, true)
}

#[cfg(test)]
mod tests {
    use super::cigar_to_identity;

    fn assert_close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 1e-9, "{actual} != {expected}");
    }

    #[test]
    fn parses_gap_compressed_identity_and_query_coverage() {
        let (identity, matched_len, fraction) = cigar_to_identity("80=5I5=10X5H");
        assert_close(identity, 85.0 / 96.0);
        assert_eq!(matched_len, 85);
        assert_close(fraction, 1.0);

        let (identity, matched_len, fraction) = cigar_to_identity("50=50X");
        assert_close(identity, 0.5);
        assert_eq!(matched_len, 50);
        assert_close(fraction, 1.0);

        let (identity, matched_len, fraction) = cigar_to_identity("100M");
        assert_close(identity, 1.0);
        assert_eq!(matched_len, 100);
        assert_close(fraction, 1.0);
    }
}

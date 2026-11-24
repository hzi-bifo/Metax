//! Coverage-distribution statistics of the species tests.

use statrs::function::beta::beta_reg;
use statrs::function::erf::erf;
use statrs::function::gamma::ln_gamma;

fn clamp01(x: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    }
    if x > 1.0 {
        return 1.0;
    }
    x
}

/// Check if a value is a finite number (not NaN or infinity).
fn is_finite_number(x: f64) -> bool {
    x.is_finite()
}

/// Infer genome size (n) from read count and depth.
///
/// n = read_count / depth
///
/// Used when genome size is not directly available but can be
/// estimated from coverage statistics.
fn infer_n_from_depth(read_count: usize, depth: f64) -> Option<usize> {
    if read_count == 0 || depth <= 0.0 || !is_finite_number(depth) {
        return None;
    }
    let n = (read_count as f64 / depth).round() as isize;
    Some(n.max(1) as usize)
}

/// Infer genome size (n) from read count and expected breadth.
///
/// Uses the relationship: exp_breadth = 1 - (1 - 1/n)^read_count
/// Solving for n: n = -read_count / ln(1 - exp_breadth)
fn infer_n_from_exp_breadth(read_count: usize, exp_breadth: f64) -> Option<usize> {
    if read_count == 0 {
        return None;
    }
    let p = clamp01(exp_breadth);
    if p <= 0.0 || p >= 1.0 {
        return Some(read_count.max(1));
    }
    let denom = (1.0 - p).ln();
    if !denom.is_finite() || denom == 0.0 {
        return None;
    }
    let n = (-(read_count as f64) / denom).round() as isize;
    Some(n.max(1) as usize)
}

/// Convert breadth (fraction) to discrete count of covered positions.
///
/// k = round(breadth * n), clamped to [0, n]
fn k_from_breadth(breadth: f64, n: usize) -> Option<usize> {
    if !is_finite_number(breadth) || !(0.0..=1.0).contains(&breadth) {
        return None;
    }
    let mut k = (breadth * n as f64).round() as isize;
    if k < 0 {
        k = 0;
    }
    if k > n as isize {
        k = n as isize;
    }
    Some(k as usize)
}

/// Compute Stirling numbers of the second kind S(m, k) for k in [0, kmax].
///
/// S(m, k) counts the number of ways to partition m elements into k
/// non-empty subsets. Used in the exact p-value calculation.
///
/// Uses the recurrence: S(m, k) = k * S(m-1, k) + S(m-1, k-1)
fn stirling2_row(m: usize, kmax: usize) -> Vec<f64> {
    let mut s_prev = vec![0.0_f64; kmax + 1];
    s_prev[0] = 1.0;

    for i in 1..=m {
        let mut s_cur = vec![0.0_f64; kmax + 1];
        let upper = i.min(kmax);
        for k in 1..=upper {
            s_cur[k] = (s_prev[k] * k as f64) + s_prev[k - 1];
        }
        s_prev = s_cur;
    }

    s_prev
}

/// Compute log(k!) using the gamma function.
/// ln(k!) = ln(Gamma(k+1))
fn log_factorial(k: usize) -> f64 {
    ln_gamma(k as f64 + 1.0)
}

/// Compute log of binomial coefficient C(n, k) = n! / (k! * (n-k)!)
/// Uses log-space arithmetic to avoid overflow.
fn log_combination(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

/// Compute log(sum(exp(values))) in a numerically stable way.
///
/// Uses the log-sum-exp trick: log(sum(exp(x))) = m + log(sum(exp(x - m)))
/// where m = max(x). This avoids overflow/underflow.
fn logsumexp(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }
    let m = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !m.is_finite() {
        return m;
    }
    let mut acc = 0.0;
    for v in values {
        acc += (v - m).exp();
    }
    m + acc.ln()
}

/// Which tail of the distribution to test.
#[derive(Clone, Copy)]
pub(crate) enum Tail {
    /// Test if observed is significantly lower than expected
    Lower,
    /// Test if observed is significantly higher than expected
    Upper,
}

/// The expected breadth `p` (clamped) and the number of genome positions `n`
/// that a breadth test runs on; `None` without reads, with an expected breadth
/// of 0 or 1, or when `n` cannot be inferred.
fn breadth_test_inputs(
    read_count: usize,
    depth: Option<f64>,
    exp_breadth: f64,
) -> Option<(f64, usize)> {
    if read_count == 0 {
        return None;
    }
    let p = clamp01(exp_breadth);
    if p <= 0.0 || p >= 1.0 {
        return None;
    }
    let n = depth
        .filter(|d| *d > 0.0)
        .and_then(|d| infer_n_from_depth(read_count, d))
        .or_else(|| infer_n_from_exp_breadth(read_count, p))?;
    Some((p, n))
}

/// Compute exact p-value using Stirling numbers.
///
/// For small sample sizes (n ≤ 3000, m ≤ 100), computes the exact
/// probability using Stirling numbers of the second kind.
///
/// The distribution is based on the "balls into bins" model where
/// m reads are randomly placed into n genome positions.
fn pvalue_exact_stirling(
    read_count: usize,
    depth: Option<f64>,
    breadth: f64,
    exp_breadth: f64,
    side: Tail,
    max_n: usize,
    max_m: usize,
) -> Option<f64> {
    let (_, n) = breadth_test_inputs(read_count, depth, exp_breadth)?;

    if n > max_n || read_count > max_m {
        return None;
    }

    let mut k_obs = k_from_breadth(breadth, n)?;
    let k_max_possible = n.min(read_count);
    if k_obs > k_max_possible {
        k_obs = k_max_possible;
    }

    if k_obs == 0 {
        return match side {
            Tail::Lower => Some(0.0),
            Tail::Upper => Some(1.0),
        };
    }

    let s_values = stirling2_row(read_count, k_max_possible);
    let log_den = (read_count as f64) * (n as f64).ln();

    let pmf_log = |k: usize| -> f64 {
        if k == 0 || k > k_max_possible {
            return f64::NEG_INFINITY;
        }
        let log_comb = log_combination(n, k);
        let log_stir = s_values[k].ln();
        let log_fact = log_factorial(k);
        log_comb + log_stir + log_fact - log_den
    };

    let cdf_log = |k: usize| -> f64 {
        if k == 0 {
            return f64::NEG_INFINITY;
        }
        if k >= k_max_possible {
            return 0.0;
        }
        let mut logs = Vec::with_capacity(k);
        for t in 1..=k {
            logs.push(pmf_log(t));
        }
        logsumexp(&logs)
    };

    let value_log = match side {
        Tail::Lower => cdf_log(k_obs),
        Tail::Upper => {
            // P(X >= k_obs) = 1 - P(X <= k_obs-1)
            let lower_log = if k_obs > 1 {
                cdf_log(k_obs - 1)
            } else {
                f64::NEG_INFINITY
            };
            let lower_prob = if lower_log.is_finite() {
                lower_log.exp()
            } else {
                0.0
            };
            clamp01(1.0 - lower_prob).ln()
        }
    };

    Some(clamp01(value_log.exp()))
}

/// Compute p-value using binomial approximation.
///
/// Uses the regularized incomplete beta function to compute the
/// CDF of a binomial distribution with parameters n and p.
///
/// More efficient than exact calculation for moderate sample sizes.
fn pvalue_binomial_approx(
    read_count: usize,
    depth: Option<f64>,
    breadth: f64,
    exp_breadth: f64,
    side: Tail,
) -> Option<f64> {
    let (p, n) = breadth_test_inputs(read_count, depth, exp_breadth)?;

    let k_obs = k_from_breadth(breadth, n)?;

    let cdf = |k: isize| -> f64 {
        if k < 0 {
            return 0.0;
        }
        if k as usize >= n {
            return 1.0;
        }
        beta_reg((n - k as usize) as f64, (k as usize + 1) as f64, 1.0 - p)
    };

    let val = match side {
        Tail::Lower => cdf(k_obs as isize),
        Tail::Upper => 1.0 - cdf(k_obs as isize - 1),
    };

    Some(clamp01(val))
}

/// Compute p-value using normal approximation to binomial.
///
/// For large samples where np(1-p) ≥ 20 and 0.05 ≤ p ≤ 0.95,
/// the binomial distribution is well-approximated by a normal.
///
/// Uses continuity correction (+/- 0.5) for better accuracy.
fn pvalue_normal_approx(
    read_count: usize,
    depth: Option<f64>,
    breadth: f64,
    exp_breadth: f64,
    side: Tail,
) -> Option<f64> {
    let (p, n) = breadth_test_inputs(read_count, depth, exp_breadth)?;

    let k_obs = k_from_breadth(breadth, n)? as f64;

    let mu = n as f64 * p;
    let var = n as f64 * p * (1.0 - p);
    if var <= 0.0 {
        return None;
    }
    let sd = var.sqrt();

    let norm_cdf = |x: f64| 0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2));

    let val = match side {
        Tail::Lower => {
            let z = (k_obs + 0.5 - mu) / sd;
            norm_cdf(z)
        }
        Tail::Upper => {
            let z = (k_obs - 0.5 - mu) / sd;
            1.0 - norm_cdf(z)
        }
    };

    Some(clamp01(val))
}

/// Calculate p-value for observed coverage breadth.
///
/// Selects the appropriate method based on sample size:
/// 1. Exact Stirling (n ≤ 3000, m ≤ 100) - most accurate
/// 2. Normal approximation (eff ≥ 20, 0.05 ≤ p ≤ 0.95) - fast
/// 3. Binomial approximation - fallback
///
/// # Arguments
/// * `read_count` - Number of reads mapping to this taxon
/// * `depth` - Coverage depth (optional, for inferring genome size)
/// * `breadth` - Observed breadth of coverage
/// * `exp_breadth` - Expected breadth under null hypothesis
/// * `side` - Which tail to test
///
/// # Returns
///
/// P-value in range [0, 1], or None if calculation fails
pub(crate) fn calc_breadth_pvalue(
    read_count: usize,
    depth: Option<f64>,
    breadth: f64,
    exp_breadth: Option<f64>,
    side: Tail,
) -> Option<f64> {
    if read_count == 0 || !is_finite_number(breadth) || !(0.0..=1.0).contains(&breadth) {
        return None;
    }

    let mut exp_breadth = exp_breadth;
    if exp_breadth.is_none() {
        exp_breadth = depth.filter(|d| *d > 0.0).map(|d| 1.0 - (-d).exp());
    }
    let exp_breadth = exp_breadth?;
    let p = clamp01(if exp_breadth >= 1.0 {
        1.0 - 1e-12
    } else {
        exp_breadth
    });
    if p <= 0.0 || p >= 1.0 {
        return None;
    }

    let mut n = depth
        .filter(|d| *d > 0.0)
        .and_then(|d| infer_n_from_depth(read_count, d));
    if n.is_none() {
        n = infer_n_from_exp_breadth(read_count, p);
    }
    let n = n?;
    let eff = (n as f64 * p).min(n as f64 * (1.0 - p));

    if n <= 3000 && read_count <= 100 {
        if let Some(pv) = pvalue_exact_stirling(read_count, depth, breadth, p, side, 3000, 100) {
            return Some(pv);
        }
    }

    if eff >= 20.0 && (0.05..=0.95).contains(&p) {
        if let Some(pv) = pvalue_normal_approx(read_count, depth, breadth, p, side) {
            return Some(pv);
        }
    }

    pvalue_binomial_approx(read_count, depth, breadth, p, side)
}

/// Observed-to-expected coverage ratios that agree with the model. A coverage
/// p-value inside this window is 1.0 without testing: with many reads any small
/// deviation would test as significant, so a p-value would only measure
/// sampling depth.
pub(crate) const NEUTRAL_RATIO: std::ops::RangeInclusive<f64> = 0.9..=1.2;

#[cfg(test)]
mod tests {
    use super::{calc_breadth_pvalue, Tail};

    #[test]
    fn handles_expected_breadth_near_one() {
        let value = calc_breadth_pvalue(5, Some(10.0), 0.8, Some(1.0), Tail::Lower)
            .expect("expected a numeric p-value");
        assert!(value.is_finite());
        assert!((0.0..=1.0).contains(&value));
    }
}

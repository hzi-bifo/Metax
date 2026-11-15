//! Expectation–maximization over candidate units weighted by how well each
//! read matches each unit.
//!
//! A read r that aligns to several units (genomes, species or genera) is
//! shared among them in proportion to `share(u) · L(r | u)`, where
//! `L(r | u) = exp(−λ · errors(r, u))` and `errors = (1 − identity) · mapped
//! bases` of the read's best alignment to the unit. With a per-base error rate
//! `e`, a base that matches has probability `1 − e` and a given wrong base
//! `e / 3`, so each error lowers the likelihood by `λ = ln((1 − e) / (e / 3))`.
//! Reads that match two units equally leave the shares where the other reads
//! put them; reads that match one unit better move share towards it, so a
//! unit supported only by reads that another unit explains better ends near
//! zero.
//!
//! Reads with the same candidates and the same rounded error counts form one
//! group, so an iteration costs the number of distinct groups, not of reads.

// Used by the level-wise engine, which follows in the next commits.
#![allow(dead_code)]

use std::collections::HashMap;

/// Error counts are rounded to this step when reads are grouped.
const ERROR_STEP: f64 = 0.5;
/// Lowest and highest per-base error rate used for λ.
const MIN_ERROR_RATE: f64 = 1e-4;
const MAX_ERROR_RATE: f64 = 0.3;

/// Penalty per error for a per-base error rate `error_rate`.
pub(crate) fn error_penalty(error_rate: f64) -> f64 {
    let e = if error_rate.is_finite() {
        error_rate.clamp(MIN_ERROR_RATE, MAX_ERROR_RATE)
    } else {
        MAX_ERROR_RATE
    };
    ((1.0 - e) / (e / 3.0)).ln()
}

/// One read's candidates: `(unit, errors)` with `errors` the estimated number
/// of errors of the read's best alignment to the unit.
pub(crate) type Candidates = Vec<(u32, f64)>;

/// Reads grouped by identical candidates and rounded error counts.
#[derive(Default)]
pub(crate) struct ReadGroups {
    index: HashMap<Vec<(u32, u32)>, usize>,
    groups: Vec<Group>,
}

struct Group {
    units: Vec<u32>,
    /// Errors above the group's best candidate, in units of [`ERROR_STEP`].
    excess_steps: Vec<u32>,
    weight: f64,
}

impl ReadGroups {
    /// Add a read with weight `weight` (1 for a read or fragment; larger for
    /// scaled reads). Candidates of one unit keep their fewest errors.
    pub(crate) fn add(&mut self, mut candidates: Candidates, weight: f64) {
        if candidates.is_empty() || weight.is_nan() || weight <= 0.0 {
            return;
        }
        candidates.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
        candidates.dedup_by_key(|c| c.0);
        let best = candidates.iter().map(|c| c.1).fold(f64::INFINITY, f64::min);
        let key = candidates
            .iter()
            .map(|&(unit, errors)| (unit, ((errors - best).max(0.0) / ERROR_STEP).round() as u32))
            .collect::<Vec<_>>();
        if let Some(&slot) = self.index.get(&key) {
            self.groups[slot].weight += weight;
            return;
        }
        self.groups.push(Group {
            units: key.iter().map(|k| k.0).collect(),
            excess_steps: key.iter().map(|k| k.1).collect(),
            weight,
        });
        self.index.insert(key, self.groups.len() - 1);
    }

    pub(crate) fn len(&self) -> usize {
        self.groups.len()
    }
}

/// Result of [`run`]: each unit's share of the grouped reads' weight and the
/// posterior of every candidate of every group, in the order the groups were
/// added.
pub(crate) struct EmResult {
    pub(crate) weight: HashMap<u32, f64>,
    pub(crate) posteriors: Vec<Vec<(u32, f64)>>,
    pub(crate) iterations: usize,
}

/// Share the grouped reads among their candidates. `prior` gives the starting
/// weight of units (for example reads that units have alone); units without
/// one start equal. Stops when no unit's share changes by more than `tolerance`
/// (relative to the total) or after `max_iterations` rounds.
pub(crate) fn run(
    groups: &ReadGroups,
    prior: &HashMap<u32, f64>,
    penalty: f64,
    max_iterations: usize,
    tolerance: f64,
) -> EmResult {
    let total: f64 = groups.groups.iter().map(|g| g.weight).sum();
    let mut share: HashMap<u32, f64> = HashMap::new();
    for group in &groups.groups {
        for &unit in &group.units {
            share.entry(unit).or_insert(1.0);
        }
    }
    for (unit, value) in share.iter_mut() {
        if let Some(&p) = prior.get(unit) {
            *value = p.max(f64::MIN_POSITIVE);
        }
    }
    // Likelihood of each candidate relative to its group's best candidate.
    let relative = groups
        .groups
        .iter()
        .map(|g| {
            g.excess_steps
                .iter()
                .map(|&s| (-penalty * ERROR_STEP * f64::from(s)).exp())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut iterations = 0;
    let mut next: HashMap<u32, f64> = HashMap::with_capacity(share.len());
    for iteration in 0..max_iterations.max(1) {
        iterations = iteration + 1;
        next.clear();
        for (group, rel) in groups.groups.iter().zip(&relative) {
            let norm: f64 = group.units.iter().zip(rel).map(|(u, r)| share[u] * r).sum();
            if norm <= 0.0 {
                continue;
            }
            for (unit, r) in group.units.iter().zip(rel) {
                *next.entry(*unit).or_insert(0.0) += group.weight * share[unit] * r / norm;
            }
        }
        let change = share
            .iter()
            .map(|(u, old)| (next.get(u).copied().unwrap_or(0.0) - old).abs())
            .fold(0.0, f64::max);
        for (unit, value) in share.iter_mut() {
            *value = next.get(unit).copied().unwrap_or(0.0);
        }
        if change <= tolerance * total.max(1e-12) {
            break;
        }
    }
    let posteriors = groups
        .groups
        .iter()
        .zip(&relative)
        .map(|(group, rel)| {
            let norm: f64 = group.units.iter().zip(rel).map(|(u, r)| share[u] * r).sum();
            group
                .units
                .iter()
                .zip(rel)
                .map(|(u, r)| (*u, if norm > 0.0 { share[u] * r / norm } else { 0.0 }))
                .collect()
        })
        .collect();
    EmResult {
        weight: share,
        posteriors,
        iterations,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_default(groups: &ReadGroups) -> EmResult {
        run(groups, &HashMap::new(), error_penalty(0.003), 1_000, 1e-9)
    }

    #[test]
    fn penalty_grows_as_the_error_rate_falls() {
        assert!(error_penalty(0.001) > error_penalty(0.01));
        assert!(error_penalty(0.01) > error_penalty(0.1));
        assert!((error_penalty(0.03) - (0.97f64 / 0.01).ln()).abs() < 1e-12);
        assert_eq!(error_penalty(0.0), error_penalty(MIN_ERROR_RATE));
        assert_eq!(error_penalty(f64::NAN), error_penalty(MAX_ERROR_RATE));
    }

    #[test]
    fn reads_matching_one_unit_better_pull_the_tied_reads_to_it() {
        // 900 reads match units 1 and 2 equally; 100 reads match unit 1 with
        // one error fewer. Unit 2 has no read of its own and ends near zero.
        let mut groups = ReadGroups::default();
        for _ in 0..900 {
            groups.add(vec![(1, 0.5), (2, 0.5)], 1.0);
        }
        for _ in 0..100 {
            groups.add(vec![(1, 0.0), (2, 1.0)], 1.0);
        }
        assert_eq!(groups.len(), 2);
        let result = run_default(&groups);
        assert!(result.weight[&1] > 990.0, "unit 1 {}", result.weight[&1]);
        assert!(result.weight[&2] < 10.0, "unit 2 {}", result.weight[&2]);
    }

    #[test]
    fn two_present_units_keep_their_shares() {
        // Unit 1 has 300 better-matching reads, unit 2 has 100; 600 are tied.
        let mut groups = ReadGroups::default();
        for _ in 0..600 {
            groups.add(vec![(1, 0.0), (2, 0.0)], 1.0);
        }
        for _ in 0..300 {
            groups.add(vec![(1, 0.0), (2, 2.0)], 1.0);
        }
        for _ in 0..100 {
            groups.add(vec![(1, 2.0), (2, 0.0)], 1.0);
        }
        let result = run_default(&groups);
        let ratio = result.weight[&1] / result.weight[&2];
        assert!((ratio - 3.0).abs() < 0.1, "ratio {ratio}");
        let total = result.weight[&1] + result.weight[&2];
        assert!((total - 1_000.0).abs() < 1e-6);
    }

    #[test]
    fn posteriors_sum_to_one_and_follow_the_order_of_groups() {
        let mut groups = ReadGroups::default();
        groups.add(vec![(7, 0.0)], 1.0);
        groups.add(vec![(7, 0.0), (9, 3.0)], 2.0);
        let result = run_default(&groups);
        assert_eq!(result.posteriors.len(), 2);
        assert_eq!(result.posteriors[0], vec![(7, 1.0)]);
        let sum: f64 = result.posteriors[1].iter().map(|p| p.1).sum();
        assert!((sum - 1.0).abs() < 1e-12);
        assert!(result.posteriors[1][0].1 > 0.99);
    }

    #[test]
    fn a_unit_keeps_its_fewest_errors_and_reads_are_grouped() {
        let mut groups = ReadGroups::default();
        groups.add(vec![(1, 2.0), (1, 0.0), (2, 0.4)], 1.0);
        groups.add(vec![(2, 0.6), (1, 0.1)], 1.0);
        assert_eq!(
            groups.len(),
            1,
            "same units and rounded excess errors form one group"
        );
        groups.add(vec![], 1.0);
        groups.add(vec![(1, 0.0)], 0.0);
        assert_eq!(groups.len(), 1, "empty reads and zero weights are ignored");
    }

    #[test]
    fn a_prior_sets_the_split_of_reads_that_carry_no_signal() {
        let mut groups = ReadGroups::default();
        for _ in 0..100 {
            groups.add(vec![(1, 0.0), (2, 0.0)], 1.0);
        }
        let prior = HashMap::from([(1, 3.0), (2, 1.0)]);
        let result = run(&groups, &prior, error_penalty(0.003), 1_000, 1e-9);
        assert!((result.weight[&1] - 75.0).abs() < 1e-6);
        assert!((result.weight[&2] - 25.0).abs() < 1e-6);
    }
}

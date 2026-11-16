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
    /// scaled reads) and return the index of its group. Candidates of one unit
    /// keep their fewest errors. Reads without candidates or weight are not
    /// added.
    pub(crate) fn add(&mut self, mut candidates: Candidates, weight: f64) -> Option<usize> {
        if candidates.is_empty() || weight.is_nan() || weight <= 0.0 {
            return None;
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
            return Some(slot);
        }
        self.groups.push(Group {
            units: key.iter().map(|k| k.0).collect(),
            excess_steps: key.iter().map(|k| k.1).collect(),
            weight,
        });
        let slot = self.groups.len() - 1;
        self.index.insert(key, slot);
        Some(slot)
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
/// one start equal. The EM is accelerated by SQUAREM (two EM steps, an
/// extrapolation along their direction, one stabilising step): when nearly
/// all reads tie between near-identical units, a plain EM moves the shares by
/// a tiny fraction per round and would need thousands of rounds. Stops when no
/// unit's share changes by more than `tolerance` (relative to the total) or
/// after `max_iterations` EM steps.
pub(crate) fn run(
    groups: &ReadGroups,
    prior: &HashMap<u32, f64>,
    penalty: f64,
    max_iterations: usize,
    tolerance: f64,
) -> EmResult {
    let mut units: Vec<u32> = groups
        .groups
        .iter()
        .flat_map(|g| g.units.iter().copied())
        .collect();
    units.sort_unstable();
    units.dedup();
    let slot: HashMap<u32, usize> = units.iter().enumerate().map(|(i, &u)| (u, i)).collect();
    // Per group: unit slots and each candidate's likelihood relative to the
    // group's best candidate.
    let indexed = groups
        .groups
        .iter()
        .map(|g| {
            let slots = g.units.iter().map(|u| slot[u]).collect::<Vec<_>>();
            let relative = g
                .excess_steps
                .iter()
                .map(|&s| (-penalty * ERROR_STEP * f64::from(s)).exp())
                .collect::<Vec<_>>();
            (slots, relative, g.weight)
        })
        .collect::<Vec<_>>();
    let total: f64 = indexed.iter().map(|g| g.2).sum();
    let mut share = units
        .iter()
        .map(|u| prior.get(u).map_or(1.0, |&p| p.max(f64::MIN_POSITIVE)))
        .collect::<Vec<_>>();
    normalise(&mut share, total);

    let step = |from: &[f64], into: &mut Vec<f64>| {
        into.clear();
        into.resize(from.len(), 0.0);
        for (slots, relative, weight) in &indexed {
            let norm: f64 = slots.iter().zip(relative).map(|(&u, r)| from[u] * r).sum();
            if norm > 0.0 {
                for (&u, r) in slots.iter().zip(relative) {
                    into[u] += weight * from[u] * r / norm;
                }
            }
        }
    };
    let (mut first, mut second, mut jumped) = (Vec::new(), Vec::new(), Vec::new());
    let mut iterations = 0;
    while iterations < max_iterations.max(1) {
        step(&share, &mut first);
        step(&first, &mut second);
        iterations += 2;
        let r = first
            .iter()
            .zip(&share)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        let v = second
            .iter()
            .zip(&first)
            .zip(&r)
            .map(|((c, b), r)| c - b - r)
            .collect::<Vec<_>>();
        let r_norm = r.iter().map(|x| x * x).sum::<f64>().sqrt();
        let v_norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        let next = if v_norm > 0.0 && r_norm > 0.0 {
            let alpha = (-r_norm / v_norm).min(-1.0);
            jumped.clear();
            jumped.extend(
                share
                    .iter()
                    .zip(&r)
                    .zip(&v)
                    .map(|((s, r), v)| (s - 2.0 * alpha * r + alpha * alpha * v).max(0.0)),
            );
            if jumped.iter().all(|x| x.is_finite()) && jumped.iter().sum::<f64>() > 0.0 {
                normalise(&mut jumped, total);
                step(&jumped, &mut first);
                iterations += 1;
                &first
            } else {
                &second
            }
        } else {
            &second
        };
        let change = next
            .iter()
            .zip(&share)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        share.clone_from(next);
        if change <= tolerance * total.max(1e-12) {
            break;
        }
    }
    let posteriors = indexed
        .iter()
        .map(|(slots, relative, _)| {
            let norm: f64 = slots.iter().zip(relative).map(|(&u, r)| share[u] * r).sum();
            slots
                .iter()
                .zip(relative)
                .map(|(&u, r)| (units[u], if norm > 0.0 { share[u] * r / norm } else { 0.0 }))
                .collect()
        })
        .collect();
    EmResult {
        weight: units.iter().copied().zip(share).collect(),
        posteriors,
        iterations,
    }
}

/// Scale `values` to sum to `total`.
fn normalise(values: &mut [f64], total: f64) {
    let sum: f64 = values.iter().sum();
    if sum > 0.0 {
        for value in values.iter_mut() {
            *value *= total / sum;
        }
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
    fn near_identical_units_are_resolved_within_the_iteration_limit() {
        // 9,990 reads tie between two units and 10 match unit 1 better: a
        // plain EM would leave unit 2 near half after 200 rounds.
        let mut groups = ReadGroups::default();
        groups.add(vec![(1, 0.0), (2, 0.0)], 9_990.0);
        groups.add(vec![(1, 0.0), (2, 1.0)], 10.0);
        let result = run(&groups, &HashMap::new(), error_penalty(0.003), 200, 1e-9);
        assert!(result.weight[&2] < 100.0, "unit 2 {}", result.weight[&2]);
        assert!(
            result.posteriors[0][0].1 > 0.99,
            "tied reads {:?}",
            result.posteriors[0]
        );
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
        assert_eq!(groups.add(vec![(1, 2.0), (1, 0.0), (2, 0.4)], 1.0), Some(0));
        assert_eq!(groups.add(vec![(2, 0.6), (1, 0.1)], 1.0), Some(0));
        assert_eq!(
            groups.len(),
            1,
            "same units and rounded excess errors form one group"
        );
        assert_eq!(groups.add(vec![], 1.0), None);
        assert_eq!(groups.add(vec![(1, 0.0)], 0.0), None);
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

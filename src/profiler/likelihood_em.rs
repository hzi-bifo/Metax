//! Expectation–maximization over candidate units weighted by how well each
//! read matches each unit.
//!
//! A read r that aligns to several units (species, genera or genomes) is
//! shared among them in proportion to `share(u) · L(r | u)`, where
//! `L(r | u) = exp(−λ · errors(r, u))`. The caller gives `errors(r, u)`: the
//! estimated errors of the read's best alignment to the unit (both mates of a
//! pair), plus a penalty for each base that this alignment covers less than
//! the read's longest one. With a per-base error rate `e`, a base that matches
//! has probability `1 − e` and a given wrong base `e / 3`, so each error lowers
//! the likelihood by `λ = ln((1 − e) / (e / 3))`. Reads that match two units
//! equally leave the shares where the other reads put them; reads that match
//! one unit better move share towards it, so a unit supported only by reads
//! that another unit explains better ends near zero.
//!
//! Reads with the same candidates and the same rounded errors above their best
//! candidate form one group, so an iteration costs the number of distinct
//! groups, not of reads.

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

/// One read's candidates: `(unit, errors)` with `errors` the read's errors on
/// the unit as the caller counts them (see the module documentation).
pub(crate) type Candidates = Vec<(u32, f64)>;

/// Reads grouped by identical candidates and rounded error counts.
#[derive(Default)]
pub(crate) struct ReadGroups {
    /// Each group's key, its candidates by unit, each with its errors above
    /// the group's best candidate in units of [`ERROR_STEP`], and the group's
    /// index.
    index: HashMap<Box<[(u32, u32)]>, usize>,
    /// Reads in each group.
    weights: Vec<f64>,
    /// The key of the read being added, kept to look groups up without
    /// allocating.
    key: Vec<(u32, u32)>,
}

impl ReadGroups {
    /// Add a read with its `candidates` and return the index of its group.
    /// Candidates of one unit keep their fewest errors. A read without
    /// candidates is not added. `candidates` is left sorted by unit, one per
    /// unit.
    pub(crate) fn add(&mut self, candidates: &mut Candidates) -> Option<usize> {
        self.add_reads(candidates, 1.0)
    }

    /// Add `reads` reads with the same `candidates`, as [`Self::add`] adds
    /// one.
    pub(crate) fn add_reads(&mut self, candidates: &mut Candidates, reads: f64) -> Option<usize> {
        if candidates.is_empty() {
            return None;
        }
        candidates.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
        candidates.dedup_by_key(|c| c.0);
        let best = candidates.iter().map(|c| c.1).fold(f64::INFINITY, f64::min);
        self.key.clear();
        self.key.extend(
            candidates.iter().map(|&(unit, errors)| {
                (unit, ((errors - best).max(0.0) / ERROR_STEP).round() as u32)
            }),
        );
        if let Some(&slot) = self.index.get(self.key.as_slice()) {
            self.weights[slot] += reads;
            return Some(slot);
        }
        let slot = self.weights.len();
        self.weights.push(reads);
        self.index.insert(self.key.as_slice().into(), slot);
        Some(slot)
    }

    pub(crate) fn len(&self) -> usize {
        self.weights.len()
    }
}

/// Result of [`run`]: each unit's share of the grouped reads and the
/// posterior of every candidate of every group, in the order the groups were
/// added.
pub(crate) struct EmResult {
    pub(crate) weight: HashMap<u32, f64>,
    pub(crate) posteriors: Vec<Vec<(u32, f64)>>,
    pub(crate) iterations: usize,
}

impl EmResult {
    /// The units' shares summed in unit order, so that the total does not
    /// depend on the order of the map.
    pub(crate) fn total_weight(&self) -> f64 {
        let mut shares = self.weight.iter().collect::<Vec<_>>();
        shares.sort_unstable_by_key(|share| share.0);
        shares.into_iter().map(|share| share.1).sum()
    }
}

/// Share the grouped reads among their candidates, every unit starting with
/// an equal share. The EM is accelerated by SQUAREM (two EM steps, an
/// extrapolation along their direction, one stabilising step): when nearly
/// all reads tie between near-identical units, a plain EM moves the shares by
/// a tiny fraction per round and would need thousands of rounds. Stops when no
/// unit's share changes by more than `tolerance` (relative to the total) or
/// after `max_iterations` EM steps.
pub(crate) fn run(
    groups: &ReadGroups,
    penalty: f64,
    max_iterations: usize,
    tolerance: f64,
) -> EmResult {
    let mut keys: Vec<&[(u32, u32)]> = vec![&[]; groups.len()];
    for (key, &group) in &groups.index {
        keys[group] = key;
    }
    let mut units: Vec<u32> = keys
        .iter()
        .flat_map(|key| key.iter().map(|c| c.0))
        .collect();
    units.sort_unstable();
    units.dedup();
    let slot: HashMap<u32, usize> = units.iter().enumerate().map(|(i, &u)| (u, i)).collect();
    // The groups' candidates, one after the other (group `g` from
    // `starts[g]` to `starts[g + 1]`): each candidate's unit slot and its
    // likelihood relative to the group's best candidate.
    let mut starts = Vec::with_capacity(keys.len() + 1);
    let (mut slots, mut relative) = (Vec::new(), Vec::new());
    starts.push(0);
    for key in &keys {
        for &(unit, steps) in key.iter() {
            slots.push(slot[&unit]);
            relative.push((-penalty * ERROR_STEP * f64::from(steps)).exp());
        }
        starts.push(slots.len());
    }
    let weights = &groups.weights;
    let total: f64 = weights.iter().sum();
    let mut share = vec![1.0; units.len()];
    normalise(&mut share, total);

    let step = |from: &[f64], into: &mut Vec<f64>| {
        into.clear();
        into.resize(from.len(), 0.0);
        for (group, weight) in weights.iter().enumerate() {
            let members = starts[group]..starts[group + 1];
            let (slots, relative) = (&slots[members.clone()], &relative[members]);
            let norm: f64 = slots.iter().zip(relative).map(|(&u, r)| from[u] * r).sum();
            if norm > 0.0 {
                for (&u, r) in slots.iter().zip(relative) {
                    into[u] += weight * from[u] * r / norm;
                }
            }
        }
    };
    let (mut first, mut second, mut jumped) = (Vec::new(), Vec::new(), Vec::new());
    let (mut r, mut v) = (Vec::new(), Vec::new());
    let mut iterations = 0;
    while iterations < max_iterations.max(1) {
        step(&share, &mut first);
        step(&first, &mut second);
        iterations += 2;
        r.clear();
        r.extend(first.iter().zip(&share).map(|(a, b)| a - b));
        v.clear();
        v.extend(
            second
                .iter()
                .zip(&first)
                .zip(&r)
                .map(|((c, b), r)| c - b - r),
        );
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
    let posteriors = (0..keys.len())
        .map(|group| {
            let members = starts[group]..starts[group + 1];
            let (slots, relative) = (&slots[members.clone()], &relative[members]);
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
        run(groups, error_penalty(0.003), 1_000, 1e-9)
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
            groups.add(&mut vec![(1, 0.5), (2, 0.5)]);
        }
        for _ in 0..100 {
            groups.add(&mut vec![(1, 0.0), (2, 1.0)]);
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
        for _ in 0..9_990 {
            groups.add(&mut vec![(1, 0.0), (2, 0.0)]);
        }
        for _ in 0..10 {
            groups.add(&mut vec![(1, 0.0), (2, 1.0)]);
        }
        let result = run(&groups, error_penalty(0.003), 200, 1e-9);
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
            groups.add(&mut vec![(1, 0.0), (2, 0.0)]);
        }
        for _ in 0..300 {
            groups.add(&mut vec![(1, 0.0), (2, 2.0)]);
        }
        for _ in 0..100 {
            groups.add(&mut vec![(1, 2.0), (2, 0.0)]);
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
        groups.add(&mut vec![(7, 0.0)]);
        for _ in 0..2 {
            groups.add(&mut vec![(7, 0.0), (9, 3.0)]);
        }
        let result = run_default(&groups);
        assert_eq!(result.posteriors.len(), 2);
        assert_eq!(result.posteriors[0], vec![(7, 1.0)]);
        let sum: f64 = result.posteriors[1].iter().map(|p| p.1).sum();
        assert!((sum - 1.0).abs() < 1e-12);
        assert!(result.posteriors[1][0].1 > 0.99);
    }

    #[test]
    fn the_total_weight_is_summed_in_unit_order() {
        // In unit order the two single reads are lost beside the large share;
        // summed in another order they would add 2.
        let result = EmResult {
            weight: HashMap::from([(3, 1.0), (2, 1e16), (1, 1.0)]),
            posteriors: Vec::new(),
            iterations: 0,
        };
        assert_eq!(result.total_weight(), 1e16);
    }

    #[test]
    fn a_unit_keeps_its_fewest_errors_and_reads_are_grouped() {
        let mut groups = ReadGroups::default();
        assert_eq!(groups.add(&mut vec![(1, 2.0), (1, 0.0), (2, 0.4)]), Some(0));
        assert_eq!(groups.add(&mut vec![(2, 0.6), (1, 0.1)]), Some(0));
        assert_eq!(
            groups.len(),
            1,
            "same units and rounded excess errors form one group"
        );
        assert_eq!(groups.add(&mut vec![(1, 0.0)]), Some(1));
        assert_eq!(groups.add(&mut vec![]), None);
        assert_eq!(groups.len(), 2, "reads without candidates are ignored");
    }
}

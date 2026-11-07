//! Species whose reads are clearly less identical than those of the sample's
//! typical species, reported at their genus.
//!
//! The reads of an organism that the database lacks land on its closest
//! relatives, at an identity well below what reads of an organism in the
//! database reach. Sequencing errors, read length and strain diversity lower
//! the identity of every species in a sample alike, so the measure is taken
//! against the sample itself: the median identity of its reported species.
//! A species far below it is most likely a relative of the organism in the
//! sample, and its reads are reported at the genus they share.

use std::collections::HashMap;

use metam::format::MAPPING_IDENTITY_SCALE;

use super::evidence::{EvidenceAlignment, MappingRole};
use super::output::{ProfileAnnotations, ProfileRow};

/// Fewest species with best placements for a baseline to be taken.
const MIN_JUDGED_SPECIES: usize = 3;

/// Fragments whose best placement lies on a taxon, and the sum of their
/// identities.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct BestPlacements {
    pub(crate) fragments: u64,
    pub(crate) identity_sum: f64,
}

impl BestPlacements {
    pub(crate) fn add(&mut self, other: Self) {
        self.fragments += other.fragments;
        self.identity_sum += other.identity_sum;
    }

    fn mean_identity(self) -> Option<f64> {
        (self.fragments > 0).then(|| self.identity_sum / self.fragments as f64)
    }
}

/// A placement of one read, or of both mates, in exact integer units:
/// identities are five-decimal values, so scores compare without rounding.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Placement {
    /// Mapped bases x (2 x identity - 1), in identity units.
    score: i64,
    bases: u64,
    /// Mapped bases x identity, in identity units.
    matches: u64,
}

impl Placement {
    fn of(row: &EvidenceAlignment) -> Self {
        let scale = i64::from(MAPPING_IDENTITY_SCALE);
        let identity = (f64::from(row.identity).clamp(0.0, 1.0) * scale as f64).round() as i64;
        let bases = i64::from(row.mapped_len);
        Self {
            score: bases * (2 * identity - scale),
            bases: bases as u64,
            matches: (bases * identity) as u64,
        }
    }

    fn join(self, other: Self) -> Self {
        Self {
            score: self.score + other.score,
            bases: self.bases + other.bases,
            matches: self.matches + other.matches,
        }
    }

    fn identity(self) -> f64 {
        self.matches as f64 / (self.bases as f64 * f64::from(MAPPING_IDENTITY_SCALE))
    }
}

/// Keep the first of equally scored placements.
fn keep_best(slot: &mut Option<Placement>, placement: Placement) {
    if slot.is_none_or(|best| placement.score > best.score) {
        *slot = Some(placement);
    }
}

#[derive(Debug)]
struct TaxonPlacements {
    taxid: u32,
    /// Best single placement of each mate (or of the read).
    mates: [Option<Placement>; 2],
    /// Best placement of both mates as one pair on one reference.
    pair: Option<Placement>,
}

impl TaxonPlacements {
    /// A pair is preferred; without one, each mate's best single placement.
    fn best(&self) -> Option<Placement> {
        self.pair.or_else(|| match self.mates {
            [Some(first), Some(second)] => Some(first.join(second)),
            [first, second] => first.or(second),
        })
    }
}

#[derive(Debug)]
struct PairPlacements {
    taxon: usize,
    reference_id: u32,
    pair_id: u32,
    mates: [Option<Placement>; 2],
}

/// Reusable buffers for [`BestPlacementScratch::credit`].
#[derive(Debug, Default)]
pub(crate) struct BestPlacementScratch {
    taxa: Vec<TaxonPlacements>,
    pairs: Vec<PairPlacements>,
}

impl BestPlacementScratch {
    /// Credit one fragment, given all its mappings, to every taxon whose best
    /// placement reaches the fragment's top score.
    pub(crate) fn credit(
        &mut self,
        alignments: &[EvidenceAlignment],
        tally: &mut HashMap<u32, BestPlacements>,
    ) {
        self.taxa.clear();
        self.pairs.clear();
        for row in alignments {
            let placement = Placement::of(row);
            let mate = usize::from(row.role == MappingRole::Mate2);
            let taxon = match self.taxa.iter().position(|t| t.taxid == row.profile_taxid) {
                Some(taxon) => taxon,
                None => {
                    self.taxa.push(TaxonPlacements {
                        taxid: row.profile_taxid,
                        mates: [None; 2],
                        pair: None,
                    });
                    self.taxa.len() - 1
                }
            };
            keep_best(&mut self.taxa[taxon].mates[mate], placement);
            if row.pair_id == 0 {
                continue;
            }
            let pair = match self.pairs.iter().position(|pair| {
                pair.reference_id == row.reference_id && pair.pair_id == row.pair_id
            }) {
                Some(pair) => pair,
                None => {
                    self.pairs.push(PairPlacements {
                        taxon,
                        reference_id: row.reference_id,
                        pair_id: row.pair_id,
                        mates: [None; 2],
                    });
                    self.pairs.len() - 1
                }
            };
            keep_best(&mut self.pairs[pair].mates[mate], placement);
        }
        for pair in &self.pairs {
            if let [Some(first), Some(second)] = pair.mates {
                keep_best(&mut self.taxa[pair.taxon].pair, first.join(second));
            }
        }
        let Some(top) = self
            .taxa
            .iter()
            .filter_map(|taxon| taxon.best())
            .map(|best| best.score)
            .max()
        else {
            return;
        };
        for taxon in &self.taxa {
            if let Some(best) = taxon.best().filter(|best| best.score == top) {
                tally.entry(taxon.taxid).or_default().add(BestPlacements {
                    fragments: 1,
                    identity_sum: best.identity(),
                });
            }
        }
    }
}

/// A reported species, judged by the identity of its best placements.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct JudgedSpecies {
    pub(crate) taxid: u32,
    /// Mean identity of the fragments placed best on the species.
    pub(crate) identity: f64,
    pub(crate) fragments: u64,
}

/// The identity of the sample's typical species.
#[derive(Debug, PartialEq)]
pub(crate) struct IdentityBaseline {
    /// Median of the species' mean identities, weighted by their fragments.
    pub(crate) baseline: f64,
    /// The species that entered the median, by taxid.
    pub(crate) species: Vec<JudgedSpecies>,
}

impl IdentityBaseline {
    /// Whether the mean identity of `species` lies more than `margin` below
    /// the baseline.
    pub(crate) fn is_below(&self, species: &JudgedSpecies, margin: f64) -> bool {
        species.identity < self.baseline - margin
    }

    /// The species more than `margin` below the baseline.
    pub(crate) fn below(&self, margin: f64) -> impl Iterator<Item = &JudgedSpecies> {
        self.species
            .iter()
            .filter(move |species| self.is_below(species, margin))
    }
}

/// The baseline of the reported `species` that have best placements; `None`
/// when fewer than three have.
pub(crate) fn identity_baseline(
    species: impl IntoIterator<Item = u32>,
    tally: &HashMap<u32, BestPlacements>,
) -> Option<IdentityBaseline> {
    let mut species: Vec<JudgedSpecies> = species
        .into_iter()
        .filter_map(|taxid| {
            let placements = tally.get(&taxid)?;
            Some(JudgedSpecies {
                taxid,
                identity: placements.mean_identity()?,
                fragments: placements.fragments,
            })
        })
        .collect();
    if species.len() < MIN_JUDGED_SPECIES {
        return None;
    }
    species.sort_by(|a, b| {
        a.identity
            .total_cmp(&b.identity)
            .then(a.fragments.cmp(&b.fragments))
            .then(a.taxid.cmp(&b.taxid))
    });
    let total: u64 = species.iter().map(|species| species.fragments).sum();
    let mut seen = 0;
    let baseline = species
        .iter()
        .find(|species| {
            seen += species.fragments;
            2 * seen >= total
        })?
        .identity;
    species.sort_unstable_by_key(|species| species.taxid);
    Some(IdentityBaseline { baseline, species })
}

/// Name and annotation columns of a genus row that the profile does not hold
/// yet; `None` when the genus cannot be reported.
pub(crate) type GenusLabel = Option<(String, Option<ProfileAnnotations>)>;

/// Move species rows into the row of their genus: `moves` maps each species
/// to its genus. A genus row adds up the fragments, depth and abundance of its
/// parts and reports the coverage of the largest one (by depth). A species
/// whose genus `label` cannot name keeps its row. Rows stay ordered by depth.
/// Returns the moves made.
pub(crate) fn move_rows_to_genus(
    rows: &mut Vec<ProfileRow>,
    moves: &HashMap<u32, u32>,
    mut label: impl FnMut(u32) -> GenusLabel,
) -> HashMap<u32, u32> {
    if moves.is_empty() {
        return HashMap::new();
    }
    let mut genus_rows: HashMap<u32, Vec<ProfileRow>> = HashMap::new();
    let mut labels: HashMap<u32, GenusLabel> = HashMap::new();
    let mut made = HashMap::new();
    let mut kept = Vec::with_capacity(rows.len());
    for row in rows.drain(..) {
        let genus = moves.get(&row.taxid).copied().filter(|&genus| {
            row.rank == "species"
                && labels
                    .entry(genus)
                    .or_insert_with(|| label(genus))
                    .is_some()
        });
        match genus {
            Some(genus) => {
                made.insert(row.taxid, genus);
                genus_rows.entry(genus).or_default().push(row);
            }
            None => kept.push(row),
        }
    }
    // Existing rows of a receiving genus join the moved species.
    let mut kept_rows = Vec::with_capacity(kept.len());
    for row in kept {
        match genus_rows.get_mut(&row.taxid) {
            Some(parts) if row.rank == "genus" => parts.insert(0, row),
            _ => kept_rows.push(row),
        }
    }
    for (genus, parts) in genus_rows {
        let largest = parts
            .iter()
            .max_by(|a, b| a.depth.total_cmp(&b.depth).then(b.taxid.cmp(&a.taxid)))
            .expect("a genus receives at least one species row");
        let mut merged = largest.clone();
        merged.equivalent_mapcount = parts.iter().map(|row| row.equivalent_mapcount).sum();
        merged.depth = parts.iter().map(|row| row.depth).sum();
        merged.abundance = parts.iter().map(|row| row.abundance).sum();
        if parts[0].rank == "genus" {
            merged.taxname.clone_from(&parts[0].taxname);
            merged.annotations.clone_from(&parts[0].annotations);
        } else if let Some(Some((name, annotations))) = labels.remove(&genus) {
            merged.taxname = name;
            merged.annotations = annotations;
        }
        merged.taxid = genus;
        merged.rank = "genus".to_string();
        kept_rows.push(merged);
    }
    if !made.is_empty() {
        kept_rows.sort_by(|a, b| b.depth.total_cmp(&a.depth).then(a.taxid.cmp(&b.taxid)));
    }
    *rows = kept_rows;
    made
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    fn row(
        taxid: u32,
        reference_id: u32,
        role: MappingRole,
        pair_id: u32,
        identity: f32,
        bases: u32,
    ) -> EvidenceAlignment {
        EvidenceAlignment {
            asm: Arc::from("asm"),
            acc: Arc::from("acc"),
            profile_taxid: taxid,
            profile_rank: crate::profiler::ProfileRank::Species,
            reference_id,
            assembly_id: reference_id,
            role,
            pair_id,
            candidate_boundary: metam::mapper::CandidateBoundary::COMPLETE,
            sampling_scale: 1.0,
            asm_len: 1.0e6,
            identity,
            fraction: 1.0,
            mapped_len: bases,
            read_len: bases,
            span: (1, 1 + bases),
            fixed_chunk: (1, 1 + bases),
            flex_chunk: (1, 1 + bases),
            marker_gene_len: 0,
        }
    }

    fn credit(rows: &[EvidenceAlignment]) -> HashMap<u32, BestPlacements> {
        let mut tally = HashMap::new();
        BestPlacementScratch::default().credit(rows, &mut tally);
        tally
    }

    #[test]
    fn a_fragment_counts_for_the_taxon_of_its_best_placement() {
        use MappingRole::{Mate1, Mate2};
        let tally = credit(&[
            row(1, 10, Mate1, 1, 0.99, 150),
            row(1, 10, Mate2, 1, 0.97, 150),
            row(2, 20, Mate1, 1, 0.95, 150),
            row(2, 20, Mate2, 1, 0.95, 150),
        ]);
        assert_eq!(tally.len(), 1);
        assert_eq!(tally[&1].fragments, 1);
        assert!((tally[&1].identity_sum - 0.98).abs() < 1e-12);
    }

    #[test]
    fn exactly_tied_taxa_all_count() {
        use MappingRole::{Mate1, Mate2};
        let tally = credit(&[
            row(1, 10, Mate1, 1, 0.98, 150),
            row(1, 10, Mate2, 1, 0.98, 150),
            row(2, 20, Mate1, 1, 0.98, 150),
            row(2, 20, Mate2, 1, 0.98, 150),
        ]);
        assert_eq!(tally[&1].fragments, 1);
        assert_eq!(tally[&2].fragments, 1);
    }

    #[test]
    fn a_pair_on_one_reference_is_preferred_over_better_single_mates() {
        use MappingRole::{Mate1, Mate2};
        // Taxon 1: a pair on reference 10, and a better single mate 1 on
        // reference 11. Its placement is the pair (score 2 x 150 x 0.9).
        // Taxon 2 has no pair; its mates' best singles score 2 x 150 x 0.92.
        let tally = credit(&[
            row(1, 10, Mate1, 1, 0.95, 150),
            row(1, 10, Mate2, 1, 0.95, 150),
            row(1, 11, Mate1, 0, 1.0, 150),
            row(2, 20, Mate1, 0, 0.96, 150),
            row(2, 21, Mate2, 0, 0.96, 150),
        ]);
        assert_eq!(tally.keys().copied().collect::<Vec<_>>(), vec![2]);
        assert!((tally[&2].identity_sum - 0.96).abs() < 1e-12);
    }

    #[test]
    fn a_lone_mate_scores_alone() {
        let tally = credit(&[row(7, 70, MappingRole::Unpaired, 0, 0.96, 100)]);
        assert_eq!(tally[&7].fragments, 1);
        assert!((tally[&7].identity_sum - 0.96).abs() < 1e-12);
    }

    fn tally(entries: &[(u32, u64, f64)]) -> HashMap<u32, BestPlacements> {
        entries
            .iter()
            .map(|&(taxid, fragments, mean)| {
                (
                    taxid,
                    BestPlacements {
                        fragments,
                        identity_sum: mean * fragments as f64,
                    },
                )
            })
            .collect()
    }

    #[test]
    fn the_baseline_is_the_fragment_weighted_median() {
        let tally = tally(&[(1, 100, 0.99), (2, 10, 0.95), (3, 10, 0.995), (4, 5, 0.96)]);
        let judged = identity_baseline([1, 2, 3, 4, 5], &tally).unwrap();
        assert!((judged.baseline - 0.99).abs() < 1e-12);
        assert_eq!(judged.species.len(), 4, "species 5 has no best placements");
        let below = |margin| {
            judged
                .below(margin)
                .map(|species| species.taxid)
                .collect::<Vec<_>>()
        };
        assert_eq!(below(0.025), vec![2, 4]);
        assert_eq!(below(0.035), vec![2]);
    }

    #[test]
    fn no_baseline_with_fewer_than_three_species() {
        let tally = tally(&[(1, 100, 0.99), (2, 10, 0.90)]);
        assert_eq!(identity_baseline([1, 2, 3], &tally), None);
    }

    fn profile_row(taxid: u32, rank: &str, depth: f64, abundance: f64) -> ProfileRow {
        ProfileRow {
            taxname: format!("taxon {taxid}"),
            taxid,
            rank: rank.to_string(),
            equivalent_mapcount: 10.0 * depth,
            depth,
            abundance,
            breadth: depth / 100.0,
            expected_breadth: 1.0,
            breadth_pvalue: None,
            fixed_chunk: 1.0,
            flex_chunk: 1.0,
            expected_flex_chunk: 1.0,
            chunk_pvalue: None,
            best_fit: None,
            annotations: None,
        }
    }

    #[test]
    fn moved_species_join_the_row_of_their_genus() {
        let mut rows = vec![
            profile_row(1, "species", 30.0, 50.0),
            profile_row(2, "species", 20.0, 30.0),
            profile_row(100, "genus", 5.0, 10.0),
            profile_row(3, "species", 4.0, 10.0),
        ];
        let moves = HashMap::from([(2, 100), (3, 200)]);
        let made = move_rows_to_genus(&mut rows, &moves, |genus| {
            Some((format!("genus {genus}"), None))
        });
        assert_eq!(made, moves);
        let summary: Vec<_> = rows
            .iter()
            .map(|row| (row.taxid, row.rank.as_str(), row.depth, row.abundance))
            .collect();
        assert_eq!(
            summary,
            vec![
                (1, "species", 30.0, 50.0),
                (100, "genus", 25.0, 40.0),
                (200, "genus", 4.0, 10.0)
            ]
        );
        // The existing genus row keeps its name; the coverage comes from the
        // largest part, species 2.
        assert_eq!(rows[1].taxname, "taxon 100");
        assert_eq!(rows[1].breadth, 0.2);
        assert_eq!(rows[2].taxname, "genus 200");
        let total: f64 = rows.iter().map(|row| row.abundance).sum();
        assert!((total - 100.0).abs() < 1e-12);
    }

    #[test]
    fn a_species_whose_genus_cannot_be_reported_stays() {
        let mut rows = vec![profile_row(1, "species", 30.0, 100.0)];
        let made = move_rows_to_genus(&mut rows, &HashMap::from([(1, 100)]), |_| None);
        assert!(made.is_empty());
        assert_eq!(rows[0].taxid, 1);
    }

    #[test]
    fn rows_stay_untouched_when_nothing_moves() {
        // Out of depth order on purpose: no move, no reordering.
        let mut rows = vec![
            profile_row(2, "species", 1.0, 10.0),
            profile_row(1, "species", 30.0, 90.0),
        ];
        let before = rows.clone();
        assert!(move_rows_to_genus(&mut rows, &HashMap::new(), |_| unreachable!()).is_empty());
        assert_eq!(rows, before);
    }
}

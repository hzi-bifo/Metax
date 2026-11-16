//! Level-wise profiling engine (development; selected with `--level-engine`).
//!
//! Each fragment's candidates are the species whose genomes it aligns to
//! within the read-level gates; candidates whose estimated errors exceed the
//! fragment's best by less than a margin count as tied. A first EM, weighted
//! by how well each fragment matches each species ([`likelihood_em`]), gives a
//! fragment to a species when its posterior there is at least 0.9; such a
//! fragment supports every genome of that species it aligns to. A species is
//! kept when one of its genomes shows enough evenly spread coverage and the
//! species wins enough of the fragments it could explain, and is reported on
//! the kept genome where most of its fragments fit best. A second EM shares
//! the fragments among the kept species, and abundance is the resulting depth
//! on each species' best-fit genome. With `--save-raw-profile` the engine also
//! writes `<prefix>.levels_genomes.tsv`, the tests of every measured genome.

use std::collections::HashMap;

use anyhow::Result;

use super::{CommunityProfiler, MappingData};
use crate::profiler::evidence::{
    calculate_coverage, pvalue_passes, standard_ratio_lower_bounds, CandidatePlacement,
    CandidateScratch, EvidenceCoverage, EvidenceMode, EvidencePolicy, PresenceAssembly,
    PresenceInterval, STANDARD_RATIO_UPPER,
};
use crate::profiler::likelihood_em::{self, ReadGroups};
use crate::profiler::output::{BestFitReport, ProfileRow};
use crate::profiler::{parse_reference_taxonomy, ProfileRank};

mod genomes;

/// Each base that a candidate aligns less than the read's longest alignment
/// counts as this many errors in the EM.
const UNALIGNED_ERROR_WEIGHT: f64 = 0.5;
/// Candidates whose likelihood is below this fraction of the read's best
/// candidate are dropped before the EM.
const MIN_RELATIVE_LIKELIHOOD: f64 = 1e-6;
/// Candidates whose errors exceed the read's best by less than this margin
/// count as tied: short reads by 0.5 estimated errors, long reads by this
/// share of the best candidate's aligned bases. Smaller differences are as
/// likely to come from a strain the database does not hold as from the
/// candidate itself.
const SHORT_TIE_ERRORS: f64 = 0.5;
const LONG_TIE_SHARE: f64 = 0.005;
/// A fragment resolves to a species, and adds to the coverage of that
/// species' genomes, when the first EM gives the species at least this
/// posterior.
const RESOLVE_POSTERIOR: f64 = 0.9;
/// A species is kept only when it wins (posterior >= RESOLVE_POSTERIOR) at
/// least this share of the fragments that list it among their candidates: a
/// species that only collects reads which a present relative explains wins
/// few of them, while a present species wins most of its reads even when its
/// strain differs from every genome of it in the database.
const MIN_WON_SHARE: f64 = 0.3;
/// Share of the read groups' weight that may still change when the EM stops.
const EM_TOLERANCE: f64 = 1e-7;
/// The sample's error rate is measured on this share of the fragments that
/// match their best candidate best.
const ERROR_RATE_SHARE: f64 = 0.75;

/// A fragment's candidate species with their errors (estimated errors of the
/// species' best placement plus the unaligned-base penalty), the read it comes
/// from, and its mapped bases (first mate of a pair) for depth.
struct Fragment {
    read: usize,
    candidates: Vec<(u32, f64)>,
    bases: f64,
}

/// Fragments shared among a set of species: the EM's read groups, each
/// fragment's group, and each group's bases.
struct SharedReads {
    groups: ReadGroups,
    group_of: Vec<Option<usize>>,
    group_bases: Vec<f64>,
}

/// Coverage of one genome from the fragments that resolve to its species.
struct GenomeEvidence {
    taxid: u32,
    assembly: PresenceAssembly,
    coverage: Option<EvidenceCoverage>,
    kept: bool,
}

pub(super) fn profile(profiler: &CommunityProfiler, data: MappingData) -> Result<()> {
    let cfg = &profiler.cfg;
    let policy = EvidencePolicy::from_config(cfg).ok_or_else(|| {
        anyhow::anyhow!(
            "the level-wise engine takes Illumina reads or long reads without --strain for now"
        )
    })?;
    let mode = policy.mode;
    let (fragments, error_rate) = collect_fragments(&data, policy);
    let penalty = likelihood_em::error_penalty(error_rate);
    log::info!(
        target: "PROFILER",
        "Level-wise engine: {} fragments with candidates; error rate {:.5}, error penalty {:.3}.",
        fragments.len(),
        error_rate,
        penalty
    );

    // First EM over every candidate species; a fragment that it gives to one
    // species with posterior >= RESOLVE_POSTERIOR adds to the coverage of
    // that species' genomes.
    let all = share_reads(&fragments, |_| true, penalty);
    let first = likelihood_em::run(
        &all.groups,
        &HashMap::new(),
        penalty,
        cfg.max_em_iterations,
        EM_TOLERANCE,
    );
    let resolved_taxid = first
        .posteriors
        .iter()
        .map(|posteriors| {
            posteriors
                .iter()
                .find(|p| p.1 >= RESOLVE_POSTERIOR)
                .map(|p| p.0)
        })
        .collect::<Vec<_>>();
    // Fragments listing each species among their candidates, and those resolved
    // to it.
    let mut candidate_fragments: HashMap<u32, u64> = HashMap::new();
    let mut resolved_candidate_fragments: HashMap<u32, u64> = HashMap::new();
    let mut won_fragments: HashMap<u32, u64> = HashMap::new();
    for &group in all.group_of.iter().flatten() {
        let resolved = resolved_taxid[group];
        for &(taxid, _) in &first.posteriors[group] {
            *candidate_fragments.entry(taxid).or_default() += 1;
            if resolved.is_some() {
                *resolved_candidate_fragments.entry(taxid).or_default() += 1;
            }
        }
        if let Some(taxid) = resolved {
            *won_fragments.entry(taxid).or_default() += 1;
        }
    }
    let mut genomes: HashMap<u32, GenomeEvidence> = HashMap::new();
    let mut scratch = CandidateScratch::default();
    for (fragment, group) in fragments.iter().zip(&all.group_of) {
        let Some(taxid) = group.and_then(|g| resolved_taxid[g]) else {
            continue;
        };
        let candidates = scratch.build(&data.read_alignments[fragment.read].alignments, policy);
        if let Some(candidate) = candidates.iter().find(|c| c.taxid == taxid) {
            for placement in &candidate.placements {
                let fits = placement.assembly_id == candidate.best.assembly_id;
                add_placement(&mut genomes, taxid, placement, fits);
            }
        }
    }

    let tests = genomes.len().max(1);
    let (oebr_default, coebr_default) = standard_ratio_lower_bounds(ProfileRank::Species, mode);
    let min_oebr = cfg.min_oebr.unwrap_or(oebr_default);
    let min_coebr = cfg.min_coebr.unwrap_or(coebr_default);
    for genome in genomes.values_mut() {
        let coverage = calculate_coverage(&mut genome.assembly, tests);
        let won = won_fragments.get(&genome.taxid).copied().unwrap_or(0) as f64;
        let candidates = candidate_fragments.get(&genome.taxid).copied().unwrap_or(0);
        genome.kept = coverage.fixed_chunk > data.min_chunk_breadth
            && within(coverage.oebr, min_oebr)
            && within(coverage.coebr, min_coebr)
            && won >= MIN_WON_SHARE * candidates.max(1) as f64;
        genome.coverage = Some(coverage);
    }
    if cfg.keep_raw {
        write_genome_diagnostics(
            profiler,
            &data,
            &genomes,
            &candidate_fragments,
            &resolved_candidate_fragments,
            &won_fragments,
            (min_oebr, min_coebr),
        )?;
    }
    // Best-fit genome of each kept species: its kept genome where most of its
    // fragments fit best, then the one with most fragments.
    let mut best_fit: HashMap<u32, u32> = HashMap::new();
    for (&assembly_id, genome) in &genomes {
        if !genome.kept {
            continue;
        }
        let key = |id: u32, g: &GenomeEvidence| {
            (g.assembly.fits, g.assembly.count, std::cmp::Reverse(id))
        };
        let better = best_fit
            .get(&genome.taxid)
            .is_none_or(|&current| key(assembly_id, genome) > key(current, &genomes[&current]));
        if better {
            best_fit.insert(genome.taxid, assembly_id);
        }
    }

    // Second EM over the kept species.
    let kept = share_reads(&fragments, |taxid| best_fit.contains_key(&taxid), penalty);
    let em = likelihood_em::run(
        &kept.groups,
        &HashMap::new(),
        penalty,
        cfg.max_em_iterations,
        EM_TOLERANCE,
    );
    let without_kept = kept.group_of.iter().filter(|g| g.is_none()).count();
    log::info!(
        target: "PROFILER",
        "Level-wise engine: {} species kept, {} read groups, EM {} + {} iterations; {} fragments without a kept species.",
        best_fit.len(),
        kept.groups.len(),
        first.iterations,
        em.iterations,
        without_kept
    );
    let group_bases = &kept.group_bases;
    // Each species' share of the fragments' bases.
    let mut bases: HashMap<u32, f64> = HashMap::new();
    for (posteriors, &group_total) in em.posteriors.iter().zip(group_bases) {
        for &(taxid, posterior) in posteriors {
            *bases.entry(taxid).or_insert(0.0) += group_total * posterior;
        }
    }

    // Genomes within each kept species, for the fragments the species-level EM
    // gives to one species.
    let species_of = kept
        .group_of
        .iter()
        .map(|group| {
            group.and_then(|g| {
                em.posteriors[g]
                    .iter()
                    .find(|p| p.1 >= RESOLVE_POSTERIOR)
                    .map(|p| p.0)
            })
        })
        .collect::<Vec<_>>();
    let genome_result = genomes::resolve(
        &data,
        &fragments,
        &species_of,
        &genomes::GenomeLevel {
            policy,
            error_rate,
            penalty,
            report_strains: cfg.genome_level,
            chunk_cutoff: data.min_chunk_breadth,
            ratio_bounds: (min_oebr, min_coebr),
            max_iterations: cfg.max_em_iterations,
            tolerance: EM_TOLERANCE,
        },
    );
    let strain_taxids = if cfg.genome_level {
        strain_taxids(cfg)
    } else {
        HashMap::new()
    };
    let species_genomes = &genome_result.species;
    let baseline_identity = genome_result.baseline_identity;
    log::info!(
        target: "PROFILER",
        "Level-wise engine: baseline identity {:.5} (median identity of the species to their best-fit genome).",
        baseline_identity
    );
    if cfg.keep_raw && cfg.genome_level {
        write_strain_diagnostics(profiler, &data, &genome_result.diagnostics)?;
    }
    // Estimated ANI of the sample's strain to a genome: its fragments' mean
    // identity to the genome relative to the sample's baseline identity.
    let ani = |identity: f64| (1.0 - (baseline_identity - identity)).min(1.0);

    // Species rows, each followed by the genomes reported on their own.
    let mut blocks: Vec<(ProfileRow, Vec<ProfileRow>)> = Vec::new();
    for (&taxid, &fallback) in &best_fit {
        let found = species_genomes.get(&taxid);
        let assembly_id = found
            .map(|g| g.best_fit)
            .filter(|id| genomes.contains_key(id))
            .unwrap_or(fallback);
        let genome = &genomes[&assembly_id];
        let coverage = genome.coverage.as_ref().expect("measured above");
        let depth = bases.get(&taxid).copied().unwrap_or(0.0) / genome.assembly.asm_len.max(1.0);
        let reads = em.weight.get(&taxid).copied().unwrap_or(0.0);
        let species_name = profiler.taxonomy.get_name(taxid).unwrap_or("unknown");
        let identity = found.map(|g| (g.mean_identity, ani(g.mean_identity)));
        let species = profile_row(
            species_name.to_string(),
            taxid,
            "species",
            (reads, depth),
            coverage,
            data.assembly_names[assembly_id as usize].clone(),
            identity,
        );
        let mut strains = Vec::new();
        for strain in found.map_or(&[][..], |g| &g.strains[..]) {
            let strain_ani = ani(strain.mean_identity);
            if strain_ani.is_nan() || strain_ani < genomes::SAME_STRAIN_ANI {
                continue;
            }
            let accession = data.assembly_names[strain.assembly_id as usize].clone();
            let (name, strain_taxid) = match strain_taxids.get(&*accession) {
                Some(&t) => (
                    profiler
                        .taxonomy
                        .get_name(t)
                        .unwrap_or(species_name)
                        .to_string(),
                    t,
                ),
                None => (format!("{species_name} {accession}"), taxid),
            };
            strains.push(profile_row(
                name,
                strain_taxid,
                "strain",
                (strain.reads, strain.bases / strain.asm_len.max(1.0)),
                &strain.coverage,
                accession,
                Some((strain.mean_identity, strain_ani)),
            ));
        }
        blocks.push((species, strains));
    }
    let total_depth: f64 = blocks.iter().map(|b| b.0.depth).sum();
    let abundance = |depth: f64| {
        if total_depth > 0.0 {
            100.0 * depth / total_depth
        } else {
            0.0
        }
    };
    for (species, strains) in &mut blocks {
        species.abundance = abundance(species.depth);
        for strain in strains {
            strain.abundance = abundance(strain.depth);
        }
    }
    blocks.sort_by(|a, b| {
        b.0.abundance
            .total_cmp(&a.0.abundance)
            .then(a.0.taxid.cmp(&b.0.taxid))
    });
    let rows = blocks
        .into_iter()
        .flat_map(|(species, strains)| std::iter::once(species).chain(strains))
        .collect::<Vec<_>>();

    let mut raw_rows = Vec::new();
    for (&assembly_id, genome) in &genomes {
        let coverage = genome.coverage.as_ref().expect("measured above");
        raw_rows.push(profile_row(
            profiler
                .taxonomy
                .get_name(genome.taxid)
                .unwrap_or("unknown")
                .to_string(),
            genome.taxid,
            "species",
            (genome.assembly.count as f64, coverage.depth),
            coverage,
            data.assembly_names[assembly_id as usize].clone(),
            None,
        ));
    }
    raw_rows.sort_by(|a, b| {
        b.equivalent_mapcount
            .total_cmp(&a.equivalent_mapcount)
            .then(a.taxid.cmp(&b.taxid))
    });
    profiler.write_profile_files(&rows, None, raw_rows)
}

/// Every fragment with at least one candidate species, and the sample's error
/// rate.
fn collect_fragments(data: &MappingData, policy: EvidencePolicy) -> (Vec<Fragment>, f64) {
    let mut scratch = CandidateScratch::default();
    let mut fragments = Vec::with_capacity(data.read_alignments.len());
    let mut best_matches: Vec<(f32, f32)> = Vec::with_capacity(data.read_alignments.len());
    for (read, bundle) in data.read_alignments.iter().enumerate() {
        let candidates = scratch.build(&bundle.alignments, policy);
        let longest = candidates
            .iter()
            .map(|c| f64::from(c.best.both_mapped))
            .fold(0.0, f64::max);
        let mut species = candidates
            .iter()
            .filter(|c| c.profile_rank == ProfileRank::Species)
            .map(|c| {
                let errors = c.best.pair_errors
                    + UNALIGNED_ERROR_WEIGHT * (longest - f64::from(c.best.both_mapped)).max(0.0);
                (c.taxid, errors, c.best.both_mapped, c.best.r1_mapped)
            })
            .collect::<Vec<_>>();
        if species.is_empty() {
            continue;
        }
        species.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        let fewest = species[0].1;
        let margin = match policy.mode {
            EvidenceMode::Long => LONG_TIE_SHARE * f64::from(species[0].2),
            EvidenceMode::Paired | EvidenceMode::Single => SHORT_TIE_ERRORS,
        };
        best_matches.push((fewest as f32, species[0].2 as f32));
        fragments.push(Fragment {
            read,
            candidates: species
                .iter()
                .map(|c| (c.0, if c.1 - fewest < margin { fewest } else { c.1 }))
                .collect(),
            bases: f64::from(species[0].3),
        });
    }
    (fragments, error_rate(&mut best_matches))
}

/// Group the fragments for the EM over the species that `allowed` admits,
/// dropping candidates whose likelihood is below [`MIN_RELATIVE_LIKELIHOOD`]
/// of the fragment's best allowed candidate.
fn share_reads(fragments: &[Fragment], allowed: impl Fn(u32) -> bool, penalty: f64) -> SharedReads {
    let max_excess_errors = -MIN_RELATIVE_LIKELIHOOD.ln() / penalty;
    let mut shared = SharedReads {
        groups: ReadGroups::default(),
        group_of: Vec::with_capacity(fragments.len()),
        group_bases: Vec::new(),
    };
    for fragment in fragments {
        let mut candidates = fragment
            .candidates
            .iter()
            .filter(|c| allowed(c.0))
            .copied()
            .collect::<Vec<_>>();
        let fewest = candidates.iter().map(|c| c.1).fold(f64::INFINITY, f64::min);
        candidates.retain(|c| c.1 - fewest <= max_excess_errors);
        let slot = shared.groups.add(candidates, 1.0);
        if let Some(slot) = slot {
            if slot == shared.group_bases.len() {
                shared.group_bases.push(0.0);
            }
            shared.group_bases[slot] += fragment.bases;
        }
        shared.group_of.push(slot);
    }
    shared
}

/// Errors per aligned base of the fragments' best candidates, over the share
/// [`ERROR_RATE_SHARE`] that match best: the sample's sequencing error rate,
/// little affected by reads from relatives of the database genomes.
fn error_rate(best_matches: &mut [(f32, f32)]) -> f64 {
    best_matches.sort_by(|a, b| (a.0 / a.1.max(1.0)).total_cmp(&(b.0 / b.1.max(1.0))));
    let kept = &best_matches[..(best_matches.len() as f64 * ERROR_RATE_SHARE).ceil() as usize];
    let errors: f64 = kept.iter().map(|m| f64::from(m.0)).sum();
    let bases: f64 = kept.iter().map(|m| f64::from(m.1)).sum();
    if bases > 0.0 {
        errors / bases
    } else {
        0.01
    }
}

fn add_placement(
    genomes: &mut HashMap<u32, GenomeEvidence>,
    taxid: u32,
    placement: &CandidatePlacement,
    fits: bool,
) {
    let genome = genomes
        .entry(placement.assembly_id)
        .or_insert_with(|| GenomeEvidence {
            taxid,
            assembly: PresenceAssembly {
                asm_len: placement.asm_len,
                sampling_scale: placement.sampling_scale,
                ..PresenceAssembly::default()
            },
            coverage: None,
            kept: false,
        });
    let assembly = &mut genome.assembly;
    assembly.count += 1;
    assembly.fits += u64::from(fits);
    assembly.reference_span_bases += u64::from(placement.span.1.saturating_sub(placement.span.0));
    assembly.intervals.push(PresenceInterval {
        reference_id: placement.reference_id,
        span: placement.span,
        mate_span: placement.mate_span,
        fixed_chunk: placement.fixed_chunk,
        flex_chunk: placement.flex_chunk,
    });
}

fn within(ratio: f64, lower: f64) -> bool {
    ratio.is_finite() && ratio >= lower && ratio <= STANDARD_RATIO_UPPER
}

/// One profile row; `amount` is `(reads, depth)`.
fn profile_row(
    name: String,
    taxid: u32,
    rank: &str,
    (reads, depth): (f64, f64),
    coverage: &EvidenceCoverage,
    genome: std::sync::Arc<str>,
    identity: Option<(f64, f64)>,
) -> ProfileRow {
    ProfileRow {
        taxname: name,
        taxid,
        rank: rank.to_string(),
        equivalent_mapcount: reads,
        depth,
        abundance: 0.0,
        breadth: coverage.breadth,
        expected_breadth: coverage.expected_breadth,
        breadth_pvalue: coverage.breadth_pvalue,
        fixed_chunk: coverage.fixed_chunk,
        flex_chunk: coverage.flex_chunk,
        expected_flex_chunk: coverage.expected_flex_chunk,
        chunk_pvalue: coverage.chunk_pvalue,
        best_fit: Some(BestFitReport {
            genome,
            evidence_depth: coverage.depth,
            identity,
        }),
        annotations: None,
    }
}

/// Strain taxid of every genome that declares one.
fn strain_taxids(cfg: &crate::profiler::ProfilerConfig) -> HashMap<std::sync::Arc<str>, u32> {
    let mut out = HashMap::new();
    let Some(database) = cfg.db_metadata.as_deref() else {
        return out;
    };
    let mut previous: Option<&str> = None;
    for reference in &database.references {
        if previous == Some(&*reference.genome) {
            continue;
        }
        previous = Some(&*reference.genome);
        if let Some(strain) = parse_reference_taxonomy(&reference.taxonomy_field)
            .ok()
            .and_then(|t| t.primary_taxid)
        {
            out.insert(reference.genome.clone(), strain);
        }
    }
    out
}

/// Development diagnostics: one row per genome with resolved fragments, with
/// its coverage, the species' candidate and won fragments, and each test.
fn write_genome_diagnostics(
    profiler: &CommunityProfiler,
    data: &MappingData,
    genomes: &HashMap<u32, GenomeEvidence>,
    candidate_fragments: &HashMap<u32, u64>,
    resolved_candidate_fragments: &HashMap<u32, u64>,
    won_fragments: &HashMap<u32, u64>,
    (min_oebr, min_coebr): (f64, f64),
) -> Result<()> {
    use std::io::Write;
    let path = format!("{}.levels_genomes.tsv", profiler.cfg.outprefix);
    let mut out = std::io::BufWriter::new(std::fs::File::create(&path)?);
    writeln!(
        out,
        "taxid\tname\tgenome\tfragments\tfits\tspecies_candidate_fragments\tspecies_resolved_candidate_fragments\tspecies_won_fragments\toebr\tcoebr\tfixed_chunk\tchunk_cutoff\tbreadth_pvalue\tchunk_pvalue\tamount_ok\tevenness_ok\tpvalue_ok"
    )?;
    let mut ids = genomes.keys().copied().collect::<Vec<_>>();
    ids.sort_unstable();
    for id in ids {
        let genome = &genomes[&id];
        let coverage = genome.coverage.as_ref().expect("measured before");
        let p = |v: Option<f64>| v.map_or("NA".to_string(), |v| format!("{v:.3e}"));
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}",
            genome.taxid,
            profiler
                .taxonomy
                .get_name(genome.taxid)
                .unwrap_or("unknown"),
            data.assembly_names[id as usize],
            genome.assembly.count,
            genome.assembly.fits,
            candidate_fragments.get(&genome.taxid).copied().unwrap_or(0),
            resolved_candidate_fragments
                .get(&genome.taxid)
                .copied()
                .unwrap_or(0),
            won_fragments.get(&genome.taxid).copied().unwrap_or(0),
            coverage.oebr,
            coverage.coebr,
            coverage.fixed_chunk,
            data.min_chunk_breadth,
            p(coverage.breadth_pvalue),
            p(coverage.chunk_pvalue),
            coverage.fixed_chunk > data.min_chunk_breadth,
            within(coverage.oebr, min_oebr) && within(coverage.coebr, min_coebr),
            pvalue_passes(
                coverage.breadth_pvalue,
                coverage.chunk_pvalue,
                profiler.cfg.strict
            ),
        )?;
    }
    Ok(())
}

/// Development diagnostics of the genome level: one row per genome that
/// fragments resolve to, and per species its best-fit genome.
fn write_strain_diagnostics(
    profiler: &CommunityProfiler,
    data: &MappingData,
    rows: &[genomes::StrainDiagnostic],
) -> Result<()> {
    use std::io::Write;
    let path = format!("{}.levels_strains.tsv", profiler.cfg.outprefix);
    let mut out = std::io::BufWriter::new(std::fs::File::create(&path)?);
    writeln!(
        out,
        "taxid\tname\tgenome\tbest_fit\tem_share\tresolved_fragments\tstrain_identity_fragments\toebr\tcoebr\tfixed_chunk\tchunk_cutoff\tkept"
    )?;
    for row in rows {
        let (oebr, coebr, chunk) = row
            .coverage
            .as_ref()
            .map_or((f64::NAN, f64::NAN, f64::NAN), |c| {
                (c.oebr, c.coebr, c.fixed_chunk)
            });
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}",
            row.taxid,
            profiler.taxonomy.get_name(row.taxid).unwrap_or("unknown"),
            data.assembly_names[row.assembly_id as usize],
            row.best_fit,
            row.em_share,
            row.resolved,
            row.gated,
            oebr,
            coebr,
            chunk,
            data.min_chunk_breadth,
            row.kept,
        )?;
    }
    Ok(())
}

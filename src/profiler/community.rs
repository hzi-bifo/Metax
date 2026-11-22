//! Community profiling: loading the mappings for the level-wise engine.
//!
//! 1. **Reference metadata**: genome, taxonomy, size and reference type of
//!    every reference come from the database metadata. MANUAL.md describes
//!    the reference format.
//!
//! 2. **Mappings**: Metam coordinate mappings, streamed from the mapper or
//!    replayed from a saved mapping file, are grouped by read (or pair), each
//!    alignment with its coverage spans: the read, a fixed chunk around it, and
//!    a chunk scaled to the genome's size.
//!
//! 3. **Run-level settings**: the mapping rate, the automatic chunk-breadth
//!    cutoff of the species amount test, the inputs of the genus roll-up and
//!    of `--estimate-unclassified`, and the marker-gene model.
//!
//! 4. **Profiling**: the level-wise engine ([`levels`]) judges species, the
//!    genomes within them and the levels above, and writes the outputs.

use std::collections::{hash_map::Entry, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::Path;
use std::str;
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};

use crate::taxonomy::Taxonomy;
use csv::WriterBuilder;

use crate::db::{DatabaseMetadata, ReferenceType};
use crate::mapper::tsv::{
    for_each_metam_mapping_record, for_each_structured_mapping_record, load_shared_metam_profile,
    open_metam_mapping_reader, total_input_units, MetamMappingReceiver, MetamMappingRecord,
    MetamParseStats, MetamProfile,
};
use crate::progress;

use super::evidence::{
    witness_is_better, CandidatePlacement, EvidenceAlignment, GeneHit, MappingRole, WitnessDecision,
};
use super::markers::{FragmentGeometry, FragmentLayout, LengthHistogram, MarkerModel};
use super::output::ProfileRow;
use super::rollup::{GenomeLayout, GenusRollup};
use super::unclassified::{self, GenomeReads, MarkerSpecies, ReferenceKind};
use super::{
    coverage_passes, log_reference_taxonomy_adjustments, log_resolved_profiler_parameters,
    parse_reference_taxonomy, push_string_keyed_value, record_genome_sampling_scale,
    record_reference_taxonomy_adjustment, sampling_scale_from_sizes, scale_automatic_chunk_breadth,
    strip_pair_suffix, ProfileMode, ProfileRank, ProfilerConfig, ReferenceTaxonomy,
    ReferenceTaxonomyAdjustments, ReferenceTaxonomyResolver, ResolvedProfilerParameters,
    StableReadMap, UnresolvedReferences,
};

/// Fixed coordinate radius around each accepted read span for fixed-chunk
/// coverage.
const NATIVE_FIXED_FLANK_LEN: usize = 1_000;

fn record_terminal_rank(
    ranks: &mut HashMap<u32, ProfileRank>,
    taxid: u32,
    rank: ProfileRank,
    reference: &str,
) -> Result<()> {
    match ranks.entry(taxid) {
        Entry::Vacant(entry) => {
            entry.insert(rank);
            Ok(())
        }
        Entry::Occupied(entry) if *entry.get() == rank => Ok(()),
        Entry::Occupied(entry) => bail!(
            "conflicting terminal ranks for taxid {}: '{}' and '{}' (reference '{}')",
            taxid,
            entry.get().as_str(),
            rank.as_str(),
            reference,
        ),
    }
}

mod levels;

/// Information about a single alignment of a read to a reference.
///
/// Derived fields are computed once where the alignment enters the profiler --
/// in [`CommunityProfiler::entry_from_metam_mapping`] -- and stored here ready
/// to use, so the profiler touches no `HashMap` or `String::split` per
/// alignment.
type AlignmentEntry = EvidenceAlignment;

/// Narrow a genomic coordinate/length to `u32` for the packed
/// [`AlignmentEntry`] fields. Coordinates are bounded by `ref_len + 1`, and the
/// index rejects references `>= 2^32` bases, so the `.min` only guards the
/// unreachable 4 Gbp-contig edge — every real value is preserved exactly.
#[inline]
fn clamp_u32(value: usize) -> u32 {
    value.min(u32::MAX as usize) as u32
}

/// Bundle of all alignments for a single read.
///
/// A read may align to multiple reference sequences (multi-mapping).
/// This structure groups all alignments for downstream processing.
#[derive(Clone, Debug)]
struct ReadBundle {
    /// Read identifier (pair suffix stripped for paired-end)
    read_name: Arc<str>,

    /// All alignments for this read, ordered by quality
    alignments: Vec<AlignmentEntry>,
}

/// Reference metadata pre-parsed for the Metam coordinate path.
///
/// Metam emits the reference name on every mapping row. Splitting that
/// reference name and re-looking-up assembly lengths in the hot loop creates
/// avoidable allocator and hash-table work. This cache moves that work to DB
/// load time, where each reference is parsed once.
#[derive(Clone, Debug)]
/// A database reference as the mapping loader needs it; values shared by its
/// genome are in [`MetamAssembly`].
struct MetamReferenceEntry {
    ref_len: usize,
    profile_taxid: u32,
    profile_rank: ProfileRank,
    assembly_id: u32,
    is_marker: bool,
}

/// Values of a genome that every mapping record on it needs.
struct MetamAssembly {
    asm_len: f64,
    flex_len: usize,
    sampling_scale: f32,
}

/// The references that mapping records name, found by name and indexed by
/// reference id. A reference whose genome has no length has no entry, and
/// records on it are skipped.
struct MetamReferences<'a> {
    ids: HashMap<&'a str, u32>,
    entries: Vec<Option<MetamReferenceEntry>>,
    assemblies: Vec<MetamAssembly>,
}

impl<'a> MetamReferences<'a> {
    fn get(&self, name: &str) -> Option<(u32, &MetamReferenceEntry, &MetamAssembly)> {
        let id = *self.ids.get(name)?;
        let entry = self.entries[id as usize].as_ref()?;
        Some((id, entry, &self.assemblies[entry.assembly_id as usize]))
    }
}

/// The profiler of a sample: loads its mappings and profiles them with the
/// level-wise engine.
pub struct CommunityProfiler {
    cfg: ProfilerConfig,

    /// Taxonomy database (shared across threads)
    taxonomy: Arc<Taxonomy>,
}

impl CommunityProfiler {
    /// Create a profiler with the given configuration.
    ///
    /// # Arguments
    /// * `cfg` - Profiler configuration
    /// * `taxonomy` - NCBI taxonomy database (wrapped in Arc)
    pub fn new(cfg: ProfilerConfig, taxonomy: Arc<Taxonomy>) -> Result<Self> {
        Ok(Self { cfg, taxonomy })
    }

    /// Load the saved or staged mapping file and profile it.
    pub fn run(mut self) -> Result<()> {
        log::info!(target: "PROFILER", "Initializing the profiler!");

        let mappings = self
            .cfg
            .mapper_mappings
            .clone()
            .ok_or_else(|| anyhow!("missing Metax mapping file"))?;
        let db = Arc::clone(
            self.cfg
                .db_metadata
                .as_ref()
                .ok_or_else(|| anyhow!("missing reference metadata for native mappings"))?,
        );
        let mapping_data = self.load_metam(Path::new(&mappings), &db)?;
        self.profile(mapping_data)
    }

    pub fn run_with_metam_mappings(
        mut self,
        receiver: MetamMappingReceiver,
        source: &str,
    ) -> Result<()> {
        log::info!(target: "PROFILER", "Initializing the profiler!");
        let db = Arc::clone(
            self.cfg
                .db_metadata
                .as_ref()
                .ok_or_else(|| anyhow!("missing Metax database metadata for streamed mappings"))?,
        );
        let mapping_data = self.load_metam_mappings(receiver, source, &db)?;
        self.profile(mapping_data)
    }

    /// Lower an automatic long-read identity cutoff when the reads align worse
    /// than their base qualities claim, as when a basecaller over-states its
    /// accuracy. The upper quartile of the reads' best alignment identities
    /// then takes the claimed accuracy's place: reads from well-represented
    /// genomes align at their sequencing accuracy, while divergent or poorly
    /// represented organisms lower only the lower quantiles.
    fn lower_identity_to_alignments(&mut self, reads: &[ReadBundle]) {
        let cfg = &self.cfg;
        let Some(claimed) = cfg.identity_check_accuracy else {
            return;
        };
        let score = |aln: &&AlignmentEntry| f64::from(aln.identity) * f64::from(aln.mapped_len);
        let mut best: Vec<f32> = reads
            .iter()
            .filter_map(|read| {
                read.alignments
                    .iter()
                    .filter(|aln| {
                        coverage_passes(
                            aln.mapped_len as usize >= cfg.mapped_len,
                            f64::from(aln.fraction) >= cfg.fraction,
                            cfg.coverage_gate_is_or(),
                        )
                    })
                    .max_by(|a, b| score(a).total_cmp(&score(b)))
                    .map(|aln| aln.identity)
            })
            .collect();
        if best.is_empty() {
            return;
        }
        let upper = best.len() * 3 / 4;
        let aligned = f64::from(*best.select_nth_unstable_by(upper, f32::total_cmp).1);
        let lowered = (aligned < claimed)
            .then(|| crate::adaptive_identity_from_accuracy(aligned, cfg.mode))
            .flatten()
            .filter(|&cutoff| cutoff < cfg.identity);
        if let Some(cutoff) = lowered {
            log::info!(
                target: "PROFILER",
                "Long-read identity threshold lowered to {cutoff:.4}: reads align worse than their base qualities claim (upper quartile of read identities {aligned:.4} < claimed accuracy {claimed:.4}).",
            );
            self.cfg.identity = cutoff;
        }
    }

    fn load_metam_profile(&self, embedded: Option<MetamProfile>) -> Result<MetamProfile> {
        if let Some(profile) = embedded {
            return Ok(profile);
        }
        if let Some(shared) = &self.cfg.mapping_profile_memory {
            return load_shared_metam_profile(shared);
        }
        bail!("mapping file is missing embedded profile counters; rerun mapping with this Metax version")
    }

    fn load_metam(&mut self, mappings_path: &Path, db: &DatabaseMetadata) -> Result<MappingData> {
        if !mappings_path.exists() {
            bail!(
                "Metam mapping file {} does not exist",
                mappings_path.display()
            );
        }
        let reader = open_metam_mapping_reader(mappings_path)?;
        self.load_metam_reader(reader, &mappings_path.display().to_string(), db)
    }

    fn load_metam_reader<R: BufRead>(
        &mut self,
        reader: R,
        source: &str,
        db: &DatabaseMetadata,
    ) -> Result<MappingData> {
        let lenient = self.cfg.lenient_mapping_input;
        self.load_metam_records(source, db, |on_mapping| {
            for_each_metam_mapping_record(reader, source, lenient, on_mapping)
        })
    }

    fn load_metam_mappings(
        &mut self,
        receiver: MetamMappingReceiver,
        source: &str,
        db: &DatabaseMetadata,
    ) -> Result<MappingData> {
        self.load_metam_records(source, db, |on_mapping| {
            for_each_structured_mapping_record(receiver, source, on_mapping)
        })
    }

    fn load_metam_records<F>(
        &mut self,
        source: &str,
        db: &DatabaseMetadata,
        feed: F,
    ) -> Result<MappingData>
    where
        F: FnOnce(
            &mut dyn for<'a> FnMut(MetamMappingRecord<'a>) -> Result<()>,
        ) -> Result<MetamParseStats>,
    {
        log::info!(target: "PROFILER", "Loading Metam coordinate mappings ...");
        let fixed_flank_len = NATIVE_FIXED_FLANK_LEN;
        let mut asm_lengths_raw: HashMap<String, f64> = HashMap::new();
        let mut genome_sampling_scales: HashMap<String, f64> = HashMap::new();
        let mut profile_taxids: HashSet<u32> = HashSet::new();
        let mut terminal_ranks: HashMap<u32, ProfileRank> = HashMap::new();
        let mut reference_taxonomies = Vec::with_capacity(db.references.len());
        let mut taxonomy_adjustments = ReferenceTaxonomyAdjustments::new();
        let mut taxonomy_resolver = ReferenceTaxonomyResolver::new(&self.taxonomy);
        let mut unresolved = UnresolvedReferences::default();
        let has_marker_genomes = db.uniform_reference_type() != Some(ReferenceType::Genome);
        let mut taxids_with_genomes: HashSet<u32> = HashSet::new();
        // The references of a genome share its taxonomy field, so a field
        // equal to the previous one is not parsed again.
        let mut previous_taxonomy: Option<(&str, ReferenceTaxonomy)> = None;
        for reference in &db.references {
            let declared_taxonomy = match previous_taxonomy {
                Some((field, taxonomy)) if field == &*reference.taxonomy_field => taxonomy,
                _ => {
                    let taxonomy = parse_reference_taxonomy(&reference.taxonomy_field)
                        .with_context(|| {
                            format!(
                                "invalid taxonomy metadata for reference '{}'",
                                reference.name
                            )
                        })?;
                    previous_taxonomy = Some((&*reference.taxonomy_field, taxonomy));
                    taxonomy
                }
            };
            let taxonomy = match taxonomy_resolver.resolve(declared_taxonomy, &reference.name) {
                Ok(taxonomy) => taxonomy,
                Err(error) => {
                    unresolved.skip(
                        &reference.genome,
                        error,
                        self.cfg.ignore_unresolved_references,
                    )?;
                    reference_taxonomies.push(None);
                    continue;
                }
            };
            record_reference_taxonomy_adjustment(
                &mut taxonomy_adjustments,
                declared_taxonomy,
                taxonomy,
            );
            match asm_lengths_raw.get_mut(&*reference.genome) {
                Some(length) => *length += reference.length as f64,
                None => {
                    asm_lengths_raw.insert(reference.genome.to_string(), reference.length as f64);
                }
            }
            profile_taxids.insert(taxonomy.profile_taxid);
            if reference.reference_type == ReferenceType::Genome {
                taxids_with_genomes.insert(taxonomy.profile_taxid);
            }
            record_terminal_rank(
                &mut terminal_ranks,
                taxonomy.profile_taxid,
                taxonomy.profile_rank,
                &reference.name,
            )?;
            let scale = sampling_scale_from_sizes(
                reference.full_genome_size,
                reference.effective_indexed_size,
            );
            record_genome_sampling_scale(&mut genome_sampling_scales, &reference.genome, scale)?;
            reference_taxonomies.push(Some(taxonomy));
        }
        unresolved.log();
        log_reference_taxonomy_adjustments(
            &taxonomy_adjustments,
            self.cfg.verbose || self.cfg.very_verbose,
        );
        let asm_flex_lengths: HashMap<String, usize> = asm_lengths_raw
            .iter()
            .map(|(asm, len)| {
                let flex = len.sqrt().round() as usize;
                (asm.clone(), flex.max(1))
            })
            .collect();
        let (species_normalize_raw, mut species_to_genus_raw) =
            self.build_taxonomy_lookups(&profile_taxids);
        let normalize = |taxid: u32| species_normalize_raw.get(&taxid).copied().unwrap_or(taxid);
        let mut marker_only_species: HashMap<u32, Option<u32>> = HashMap::new();
        if has_marker_genomes {
            // A marker gene holds too little of a genome for its reads to
            // stand for a genus or family: reads failing every species there
            // come mostly from the homologous genes of relatives. Only species
            // with genome-type references have a genus here.
            let species_with_genomes: HashSet<u32> = taxids_with_genomes
                .iter()
                .map(|&taxid| normalize(taxid))
                .collect();
            if self.cfg.removed_to_genus {
                marker_only_species = profile_taxids
                    .iter()
                    .map(|&taxid| normalize(taxid))
                    .filter(|species| !species_with_genomes.contains(species))
                    .map(|species| (species, species_to_genus_raw.get(&species).copied()))
                    .collect();
            }
            species_to_genus_raw.retain(|species, _| species_with_genomes.contains(species));
        }
        // The roll-up needs every genus, also on a fractional database.
        let rollup_genus_of = self
            .cfg
            .removed_to_genus
            .then(|| species_to_genus_raw.clone());
        let db_is_subsampled = db.is_fractional;
        let species_normalize = Arc::new(species_normalize_raw);
        let genus_of_species = Arc::new(species_to_genus_raw);
        let mut assembly_ids: HashMap<&str, u32> = HashMap::new();
        let mut assembly_names: Vec<Arc<str>> = Vec::new();
        for reference in &db.references {
            if !assembly_ids.contains_key(&*reference.genome) {
                let next = u32::try_from(assembly_ids.len())
                    .context("database contains more than 2^32 assemblies")?;
                assembly_ids.insert(&*reference.genome, next);
                assembly_names.push(Arc::clone(&reference.genome));
            }
        }
        let genus_rollup = rollup_genus_of.map(|genus_of| {
            let layout = GenomeLayout::new(
                db.references
                    .iter()
                    .map(|reference| (assembly_ids[&*reference.genome], reference.length as u64)),
            );
            Arc::new(GenusRollup {
                genus_of,
                marker_only: marker_only_species,
                layout,
                fractional: db_is_subsampled,
            })
        });
        let assemblies = assembly_names
            .iter()
            .map(|asm| MetamAssembly {
                asm_len: asm_lengths_raw.get(&**asm).copied().unwrap_or(0.0),
                flex_len: asm_flex_lengths
                    .get(&**asm)
                    .copied()
                    .unwrap_or(fixed_flank_len),
                sampling_scale: genome_sampling_scales.get(&**asm).copied().unwrap_or(1.0) as f32,
            })
            .collect::<Vec<_>>();
        let mut metam_refs = MetamReferences {
            ids: HashMap::with_capacity(db.references.len()),
            entries: Vec::with_capacity(db.references.len()),
            assemblies,
        };
        for (reference_index, (reference, taxonomy)) in
            db.references.iter().zip(&reference_taxonomies).enumerate()
        {
            let assembly_id = assembly_ids[&*reference.genome];
            let has_length = metam_refs.assemblies[assembly_id as usize].asm_len > 0.0;
            // A reference without a resolvable taxonomy, or of a genome without
            // length, takes no mapping records.
            let Some(taxonomy) = (*taxonomy).filter(|_| has_length) else {
                metam_refs.entries.push(None);
                continue;
            };
            let raw_profile_taxid = taxonomy.profile_taxid;
            let profile_taxid = species_normalize
                .get(&raw_profile_taxid)
                .copied()
                .unwrap_or(raw_profile_taxid);
            let reference_id = u32::try_from(reference_index)
                .context("database contains more than 2^32 references")?;
            match metam_refs.ids.entry(reference.name.as_str()) {
                Entry::Occupied(_) => bail!(
                    "database metadata contains duplicate reference name '{}'",
                    reference.name
                ),
                Entry::Vacant(slot) => {
                    slot.insert(reference_id);
                }
            }
            metam_refs.entries.push(Some(MetamReferenceEntry {
                ref_len: reference.length,
                profile_taxid,
                profile_rank: taxonomy.profile_rank,
                assembly_id,
                is_marker: reference.reference_type == ReferenceType::Markers,
            }));
        }

        let mut read_map: StableReadMap<AlignmentEntry> = StableReadMap::default();
        let mut num_skipped_reference = 0usize;
        let mut num_alignments = 0usize;
        let mut read_lengths = LengthHistogram::default();
        let mut lone_mates = false;
        let load_progress_options = self.cfg.mapping_load_progress_options();
        let load_progress = progress::counter(load_progress_options, "Loading mapping records");
        let parse_result = feed(&mut |mapping| {
            num_alignments += 1;
            progress::tick_count(
                &load_progress,
                "Loading mapping records",
                num_alignments as u64,
                load_progress_options.interval(),
            );
            if self.cfg.is_paired && !mapping.is_paired_evidence() {
                return Ok(());
            }
            if has_marker_genomes {
                read_lengths.add(mapping.read_len);
                lone_mates |= mapping.is_mate_beyond_reference_end();
            }
            match self.entry_from_metam_mapping(mapping, &metam_refs, fixed_flank_len) {
                Some((read_id, entry)) => {
                    push_string_keyed_value(&mut read_map, read_id, entry);
                }
                None => {
                    num_skipped_reference += 1;
                }
            }
            Ok(())
        });
        load_progress.finish_and_clear();
        let parse_stats = parse_result?;
        if parse_stats.skipped_malformed > 0 {
            log::warn!(
                target: "PROFILER",
                "Skipped {} malformed mapping records from {}.",
                parse_stats.skipped_malformed,
                source
            );
        }

        let aligned_reads = read_map.len();
        let mapped_unit = if self.cfg.is_paired {
            "fragments"
        } else {
            "reads"
        };
        let profile = self.load_metam_profile(parse_stats.profile)?;
        let total_reads = total_input_units(&profile, self.cfg.is_paired).max(aligned_reads);
        let read_alignments = read_map
            .into_iter()
            .map(|(read_name, alignments)| ReadBundle {
                read_name: Arc::<str>::from(read_name),
                alignments,
            })
            .collect::<Vec<_>>();
        self.lower_identity_to_alignments(&read_alignments);
        let marker_model = if has_marker_genomes {
            let geometry =
                sample_fragment_geometry(&self.cfg, &read_alignments, &read_lengths, lone_mates);
            log::info!(target: "PROFILER", "Marker genes: {geometry}.");
            MarkerModel::from_database(db, |genome| assembly_ids[genome], geometry).map(Arc::new)
        } else {
            None
        };
        let unclassified_species = self.cfg.estimate_unclassified.then(|| {
            unclassified::marker_species(
                db.references.iter().zip(&reference_taxonomies).filter_map(
                    |(reference, taxonomy)| {
                        let taxid = taxonomy.as_ref()?.profile_taxid;
                        Some(ReferenceKind {
                            taxid: species_normalize.get(&taxid).copied().unwrap_or(taxid),
                            assembly_id: assembly_ids[&*reference.genome],
                            reference_type: reference.reference_type,
                            declared_genome_size: reference.declared_genome_size,
                        })
                    },
                ),
            )
        });

        let mapping_rate = if total_reads == 0 {
            0.0
        } else {
            aligned_reads as f64 / total_reads as f64
        };
        // The low-map gate asks how much of the sample the reference explains, so
        // fragments count as the genome share they represent.
        let gate_rate =
            genome_equivalent_mapping_rate(&read_alignments, total_reads, marker_model.as_deref());
        if gate_rate != mapping_rate {
            log::info!(
                target: "PROFILER",
                "Genome-equivalent mapping rate for the low-map gate: {:.5}.",
                gate_rate
            );
        }
        let gate_active = low_map_gate_active(&self.cfg, gate_rate);
        let chunk_plan = plan_chunk_breadth(
            &self.cfg,
            gate_rate,
            total_reads,
            aligned_reads,
            db_is_subsampled,
            gate_active,
        );
        let min_breadth = self.cfg.breadth.unwrap_or(0.0);
        let chunk_thresholds =
            resolve_community_chunk_breadth(&self.cfg, chunk_plan.estimate_reads);
        let min_chunk_breadth_value = chunk_thresholds.standard;
        let chunk_breadth_source =
            community_chunk_breadth_source(&self.cfg, chunk_plan, chunk_thresholds);

        log::info!(target: "PROFILER", "Total number of reads: {}", total_reads);
        log::info!(target: "PROFILER", "Number of mapping records: {}", num_alignments);
        log::info!(
            target: "PROFILER",
            "Number of mapped {}: {}",
            mapped_unit,
            aligned_reads
        );
        log::info!(target: "PROFILER", "Metam reported mapped read records: {}", profile.mapped_reads);
        log::info!(target: "PROFILER", "Metam emitted mapping records: {}", profile.emitted_records);
        if self.cfg.lowbiomass && self.cfg.chunk_breadth.is_none() {
            log::info!(
                target: "PROFILER",
                "Automatic minimum chunk breadth filter disabled by --low-biomass."
            );
        } else if self.cfg.chunk_breadth.is_none() {
            if chunk_plan.use_aligned_basis {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from mapped {} ({}).",
                    mapped_unit,
                    aligned_reads
                );
            } else {
                log::info!(
                    target: "PROFILER",
                    "Minimum chunk breadth estimated from total reads ({}).",
                    total_reads
                );
            }
        }
        log_chunk_breadth_scaling(chunk_thresholds);
        log_resolved_profiler_parameters(
            &self.cfg,
            ResolvedProfilerParameters {
                min_breadth,
                min_chunk_breadth: min_chunk_breadth_value,
                chunk_breadth_source,
                chunk_breadth_basis_reads: chunk_plan.estimate_reads,
                total_reads,
                mapped_reads: aligned_reads,
                mapping_rate,
            },
        );

        Ok(MappingData {
            read_alignments,
            min_chunk_breadth: min_chunk_breadth_value,
            has_scaled_genomes: genome_sampling_scales.values().any(|scale| *scale > 1.0),
            total_reads,
            genus_of_species,
            terminal_ranks: Arc::new(terminal_ranks),
            marker_model,
            genus_rollup,
            unclassified_species,
            assembly_names,
        })
    }

    fn entry_from_metam_mapping<'a>(
        &self,
        mapping: MetamMappingRecord<'a>,
        metam_refs: &MetamReferences,
        fixed_flank_len: usize,
    ) -> Option<(&'a str, AlignmentEntry)> {
        if self.cfg.is_paired && !mapping.is_paired_evidence() {
            return None;
        }
        let (reference_id, reference, assembly) = metam_refs.get(mapping.ref_id)?;
        if reference.ref_len == 0 || mapping.read_len == 0 || mapping.mapped_bases == 0 {
            return None;
        }
        let start = mapping.start_1based.clamp(1, reference.ref_len.max(1));
        let mut end = mapping.end_1based_exclusive.max(start.saturating_add(1));
        end = end.min(reference.ref_len.saturating_add(1));
        let midpoint = start + ((end.saturating_sub(start)) / 2);
        let fixed_chunk_start = start.saturating_sub(fixed_flank_len).max(1);
        let fixed_chunk_end = end
            .saturating_add(fixed_flank_len)
            .min(reference.ref_len.saturating_add(1));
        let flex_chunk_start = midpoint.saturating_sub(assembly.flex_len).max(1);
        let flex_chunk_end = midpoint
            .saturating_add(assembly.flex_len)
            .min(reference.ref_len.saturating_add(1));
        let fraction = (mapping.mapped_bases as f64 / mapping.read_len as f64).min(1.0);
        let read_name: &'a str = if self.cfg.is_paired {
            strip_pair_suffix(mapping.read_id)
        } else {
            mapping.read_id
        };
        Some((
            read_name,
            AlignmentEntry {
                profile_taxid: reference.profile_taxid,
                profile_rank: reference.profile_rank,
                reference_id,
                assembly_id: reference.assembly_id,
                role: if self.cfg.is_paired {
                    match mapping.mate_number {
                        1 => MappingRole::Mate1,
                        2 => MappingRole::Mate2,
                        _ => return None,
                    }
                } else {
                    MappingRole::Unpaired
                },
                pair_id: mapping.pair_id,
                sampling_scale: assembly.sampling_scale,
                asm_len: assembly.asm_len,
                identity: mapping.identity as f32,
                fraction: fraction as f32,
                mapped_len: clamp_u32(mapping.mapped_bases),
                read_len: clamp_u32(mapping.read_len),
                span: (clamp_u32(start), clamp_u32(end)),
                fixed_chunk: (clamp_u32(fixed_chunk_start), clamp_u32(fixed_chunk_end)),
                flex_chunk: (clamp_u32(flex_chunk_start), clamp_u32(flex_chunk_end)),
                marker_gene_len: if reference.is_marker {
                    clamp_u32(reference.ref_len)
                } else {
                    0
                },
            },
        ))
    }

    /// Precompute taxonomy lookups used during per-read assignment.
    ///
    /// Given the set of raw taxids of the database references, this produces:
    ///
    /// - `species_normalize`: sparse `raw -> species` map. Only contains
    ///   entries where the raw taxid is **below** species rank and a
    ///   species ancestor exists. Entries where raw is already species
    ///   (or has no species ancestor) are omitted; callers treat a
    ///   missing entry as "use the raw taxid unchanged". This is the
    ///   mechanism that merges multiple sub-species / strains belonging
    ///   to the same species into a single species-level assignment.
    ///
    /// - `species_to_genus`: `species -> genus` map, keyed by the
    ///   **post-normalization** species taxid so that the genus is looked
    ///   up once per species regardless of how many sub-taxon flavours
    ///   share that species.
    ///
    /// Taxids without a resolvable genus are simply absent from
    /// `species_to_genus`.
    fn build_taxonomy_lookups(
        &self,
        raw_taxids: &HashSet<u32>,
    ) -> (HashMap<u32, u32>, HashMap<u32, u32>) {
        let tax = self.taxonomy.as_ref();
        let mut species_normalize: HashMap<u32, u32> = HashMap::new();
        let mut species_to_genus: HashMap<u32, u32> = HashMap::with_capacity(raw_taxids.len());
        for &raw in raw_taxids {
            // Determine the species anchor (or fall back to the raw id
            // when no species ancestor exists, to preserve behaviour for
            // references that already point above species rank).
            let species = match resolve_species_taxid(tax, raw) {
                Some(s) => {
                    if s != raw {
                        species_normalize.insert(raw, s);
                    }
                    s
                }
                None => raw,
            };
            if let Entry::Vacant(slot) = species_to_genus.entry(species) {
                if let Some(g) = resolve_genus_taxid(tax, species) {
                    slot.insert(g);
                }
            }
        }
        (species_normalize, species_to_genus)
    }

    /// Profile the loaded mappings with the level-wise engine ([`levels`]).
    fn profile(&self, mapping_data: MappingData) -> Result<()> {
        levels::profile(self, mapping_data)
    }

    /// The unclassified share (`--estimate-unclassified`): `None` when not
    /// requested, `Some(None)` when it cannot be estimated.
    fn estimate_unclassified(
        &self,
        data: &MappingData,
        rows: &[ProfileRow],
        marker_genome_reads: &HashMap<u32, GenomeReads>,
    ) -> Option<Option<unclassified::UnclassifiedEstimate>> {
        data.unclassified_species.as_ref().map(|species| {
            // Long reads differ in length and in the genes each covers, so
            // their scales come from the reads themselves.
            let scales = unclassified::marker_scales(
                species,
                data.marker_model.as_deref(),
                marker_genome_reads,
                self.cfg.coverage_gate_is_or(),
            );
            let estimate = unclassified::estimate(rows, &scales, data.total_reads);
            match estimate {
                Ok(estimate) => log::info!(
                    target: "PROFILER",
                    "Unclassified estimate: {:.2}% of {} input {} are not explained by the profile.",
                    estimate.share * 100.0,
                    data.total_reads,
                    if self.cfg.is_paired { "fragments" } else { "reads" }
                ),
                Err(taxid) => log::warn!(
                    target: "PROFILER",
                    "Unclassified estimate unavailable: species {taxid} is represented by marker genes that declare no genome size, cannot hold a read of this length, or have no resolved long read."
                ),
            }
            estimate.ok()
        })
    }

    /// Path prefix of the profile files: `<prefix>.pathogen` with `--host`.
    fn profile_prefix(&self) -> String {
        if self.cfg.host.is_some() {
            format!("{}.pathogen", self.cfg.outprefix)
        } else {
            self.cfg.outprefix.clone()
        }
    }

    /// Write `<prefix>.profile.txt`: the `header` text, then the rows and, with
    /// `--estimate-unclassified`, the unclassified row.
    fn write_profile(
        &self,
        header: &str,
        rows: &[ProfileRow],
        unclassified: Option<Option<unclassified::UnclassifiedEstimate>>,
    ) -> Result<()> {
        let mut file = File::create(format!("{}.profile.txt", self.profile_prefix()))?;
        file.write_all(header.as_bytes())?;
        let mut profile_writer = WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(file);
        for row in rows {
            row.write(&mut profile_writer)?;
        }
        if let Some(estimate) = unclassified {
            unclassified::write_row(&mut profile_writer, estimate)?;
        }
        profile_writer.flush()?;
        Ok(())
    }

    /// Path of the per-read classifications, when they are requested and can
    /// be written for this database.
    fn classification_path(&self, data: &MappingData) -> Option<String> {
        self.cfg
            .writes_classifications(data.has_scaled_genomes)
            .then(|| {
                let classification_path = format!("{}.classify.txt.gz", self.cfg.outprefix);
                log::info!(
                    target: "PROFILER",
                    "Writing per-read classifications to {}",
                    classification_path
                );
                classification_path
            })
    }

    fn host_taxid(&self) -> Option<u32> {
        self.cfg.host.as_deref().and_then(|h| h.parse::<u32>().ok())
    }
}

/// Return the preferred coherent witness using a total, order-independent
/// comparison. The taxid tiebreak is required because `witness_is_better`
/// intentionally compares evidence strength only.
#[inline]
fn coherent_witness_max(
    left: (u32, WitnessDecision),
    right: (u32, WitnessDecision),
) -> (u32, WitnessDecision) {
    let right_is_better = witness_is_better(right.1, left.1);
    let equivalent = !right_is_better && !witness_is_better(left.1, right.1);
    if right_is_better || (equivalent && right.0 < left.0) {
        right
    } else {
        left
    }
}

/// Mapped fragments as a share of the input fragments, each counting as the
/// share of a genome it represents: the sampling scale of the genome of its
/// best placement, so that an index of sampled genomes is not mistaken for a
/// sample of unknown or host DNA. Fragments placed only on marker genes, whose
/// genome size is unknown, are left out. Equals the plain mapping rate for a
/// database of whole genomes.
fn genome_equivalent_mapping_rate(
    bundles: &[ReadBundle],
    total_reads: usize,
    markers: Option<&MarkerModel>,
) -> f64 {
    if total_reads == 0 {
        return 0.0;
    }
    let mapped: f64 = bundles
        .iter()
        .filter_map(|bundle| {
            bundle
                .alignments
                .iter()
                .find(|row| markers.is_none_or(|model| !model.is_marker_genome(row.assembly_id)))
                .map(|row| f64::from(row.sampling_scale))
        })
        .fold(0.0, |total, scale| total + scale);
    (mapped / total_reads as f64).min(1.0)
}

/// Read length, template length and read layout of this sample, measured
/// on the mappings, for the coverage expected on marker genes.
fn sample_fragment_geometry(
    cfg: &ProfilerConfig,
    bundles: &[ReadBundle],
    read_lengths: &LengthHistogram,
    lone_mates: bool,
) -> FragmentGeometry {
    let read_len = read_lengths.median().unwrap_or(1);
    let mut templates = LengthHistogram::default();
    if cfg.is_paired {
        for rows in bundles.iter().map(|bundle| &bundle.alignments) {
            for mate1 in rows
                .iter()
                .filter(|row| row.role == MappingRole::Mate1 && row.pair_id != 0)
            {
                if let Some(mate2) = rows
                    .iter()
                    .find(|row| row.role == MappingRole::Mate2 && row.pair_id == mate1.pair_id)
                {
                    let start = mate1.span.0.min(mate2.span.0);
                    let end = mate1.span.1.max(mate2.span.1);
                    templates.add(end.saturating_sub(start) as usize);
                }
            }
        }
    }
    let layout = match (cfg.is_paired, lone_mates) {
        (false, _) if cfg.coverage_gate_is_or() => FragmentLayout::LongReads,
        (false, _) => FragmentLayout::Reads,
        (true, false) => FragmentLayout::Pairs,
        (true, true) => FragmentLayout::PairsAndLoneMates,
    };
    let template_len = templates.median().unwrap_or(read_len);
    FragmentGeometry::new(layout, read_len, template_len, cfg.mapped_len, cfg.fraction)
}

/// Depth that one fragment adds to its candidate: its mapped bases over the
/// genome length. A fragment on a marker gene stands for one fragment of the
/// whole genome: it adds its whole read length, also when the read is clipped
/// at a gene end, over the placements that the genome's genes offer a read of
/// its length, so that marker and genome taxa have comparable depths; a long
/// read adds this for every gene of the genome it covers, as the placements
/// count every gene. Marker genes that cannot hold a read keep mapped bases
/// over their summed length.
fn fragment_depth(
    placement: CandidatePlacement,
    mapped: u32,
    markers: Option<&MarkerModel>,
) -> f64 {
    match markers
        .and_then(|model| model.placements_for_read(placement.assembly_id, placement.read_len))
    {
        Some(placements) => {
            f64::from(1 + placement.other_genes) * f64::from(placement.read_len) / placements
        }
        None => f64::from(mapped) / placement.asm_len.max(1.0),
    }
}

/// Record a resolved fragment's `placement` on a marker genome for the
/// unclassified estimate: for long reads (`long_reads`) the whole-genome reads
/// that its hits on the genome's genes, with `other_gene_hits`, stand for; for
/// other reads one fragment with its gene hits. Placements on other genomes
/// are not recorded.
fn add_marker_genome_reads(
    model: &MarkerModel,
    long_reads: bool,
    placement: &CandidatePlacement,
    other_gene_hits: &[GeneHit],
    reads: &mut HashMap<u32, GenomeReads>,
) {
    if long_reads {
        let genes = std::iter::once(placement.reference_id).chain(
            other_gene_hits
                .iter()
                .filter(|hit| hit.assembly_id == placement.assembly_id)
                .map(|hit| hit.interval.reference_id),
        );
        if let Some((hits, per_hit)) =
            model.whole_genome_reads(placement.assembly_id, placement.read_len, genes)
        {
            reads
                .entry(placement.assembly_id)
                .or_default()
                .add(hits, per_hit);
        }
    } else if model.is_marker_genome(placement.assembly_id) {
        reads
            .entry(placement.assembly_id)
            .or_default()
            .add_fragment(1 + placement.other_genes);
    }
}

// ============================================================================
// Family-level gates and the chunk-breadth cutoff
// ============================================================================

/// Minimum mapped length (inclusive) of a short read without a species-level
/// candidate for the family level; long reads use `cfg.mapped_len` (see
/// [`family_level_min_mapped_len`]).
const FAMILY_LEVEL_MIN_MAPPED_LEN: usize = 50;
/// Minimum aligned fraction (inclusive) of a read without a species-level
/// candidate for the family level.
const FAMILY_LEVEL_MIN_FRACTION: f64 = 0.60;

/// Mapped-length floor of the family level: a fixed 50 bases for Illumina
/// reads; for long reads the species cutoff (`cfg.mapped_len`, from defaults
/// or `-m`).
#[inline]
fn family_level_min_mapped_len(cfg: &ProfilerConfig) -> usize {
    if cfg.sequencer.is_illumina() {
        FAMILY_LEVEL_MIN_MAPPED_LEN
    } else {
        cfg.mapped_len
    }
}

/// Whether the low-map gate is active: the chunk-breadth cutoff is then
/// estimated from the mapped reads ([`plan_chunk_breadth`]).
fn low_map_gate_active(cfg: &ProfilerConfig, mapping_rate: f64) -> bool {
    cfg.is_low_map_run(mapping_rate)
}

/// Run-level decision for the read count fed into the automatic chunk-breadth
/// curve.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ChunkBreadthPlan {
    /// `true` when the unscaled basis is mapped reads; `false` for total reads.
    use_aligned_basis: bool,
    /// `true` when aligned basis was selected by the automatic low-map rule.
    by_aligned_auto: bool,
    /// Unscaled read count selected as the chunk-breadth basis.
    basis_reads: usize,
    /// Read count actually passed to the automatic chunk-breadth curve.
    estimate_reads: usize,
    /// Whether `basis_reads` was multiplied by `LOW_MAP_CHUNK_BASIS_SCALE`.
    scaled: bool,
}

/// Choose the read-count basis used to estimate chunk breadth.
#[inline]
fn plan_chunk_breadth(
    cfg: &ProfilerConfig,
    mapping_rate: f64,
    total_reads: usize,
    aligned_reads: usize,
    db_is_subsampled: bool,
    low_map_gate_active: bool,
) -> ChunkBreadthPlan {
    let chunk_breadth_pinned = cfg.chunk_breadth.is_some();
    let by_aligned_effective = cfg.resolve_by_aligned(mapping_rate);
    let basis_reads = if by_aligned_effective {
        aligned_reads
    } else {
        total_reads
    };

    let low_map_full_database = low_map_gate_active && !db_is_subsampled;
    let scaled = !chunk_breadth_pinned && by_aligned_effective && low_map_full_database;
    let estimate_reads = scale_low_map_chunk_basis(basis_reads, scaled);

    ChunkBreadthPlan {
        use_aligned_basis: by_aligned_effective,
        by_aligned_auto: by_aligned_effective && !cfg.by_aligned && low_map_gate_active,
        basis_reads,
        estimate_reads,
        scaled,
    }
}

/// Multiplier of the mapped reads that a low-map run on a full database
/// estimates its chunk-breadth cutoff from.
const LOW_MAP_CHUNK_BASIS_SCALE: f64 = 1.5;

/// Recall-mode scaling starts only once the automatic cutoff is based on more
/// than one million effective reads/fragments. This leaves shallow samples on
/// the default curve, where lowering an already-small amount gate
/// mostly adds noise.
const RECALL_CHUNK_SCALE_MIN_BASIS: usize = 1_000_000;
const SHORT_RECALL_CHUNK_SCALE_FLOOR: f64 = 0.70;
const LONG_RECALL_CHUNK_SCALE: f64 = 0.80;

/// The automatic chunk-breadth cutoff, after a user scale (`safeguard`) and
/// then recall mode (`standard`, the cutoff the species amount test uses).
#[derive(Clone, Copy, Debug, PartialEq)]
struct CommunityChunkBreadthThresholds {
    /// Calibrated read-count curve before any user or mode adjustment.
    automatic: f64,
    /// The cutoff of the species amount test: `safeguard` times the recall
    /// scale.
    standard: f64,
    /// Automatic cutoff after the explicit user scale.
    safeguard: f64,
    user_scale: f64,
    recall_scale: f64,
}

impl CommunityChunkBreadthThresholds {
    #[inline]
    fn is_user_scaled(self) -> bool {
        self.user_scale != 1.0 && self.automatic > 0.0
    }

    #[inline]
    fn is_recall_scaled(self) -> bool {
        self.recall_scale < 1.0 && self.safeguard > 0.0
    }

    #[inline]
    fn is_scaled(self) -> bool {
        self.is_user_scaled() || self.is_recall_scaled()
    }
}

/// Resolve the chunk-breadth cutoffs once.
fn resolve_community_chunk_breadth(
    cfg: &ProfilerConfig,
    basis_reads: usize,
) -> CommunityChunkBreadthThresholds {
    // Explicit user policy is authoritative in every mode, including
    // low-biomass mode and recall mode.
    if let Some(cutoff) = cfg.chunk_breadth {
        return CommunityChunkBreadthThresholds {
            automatic: cutoff,
            standard: cutoff,
            safeguard: cutoff,
            user_scale: 1.0,
            recall_scale: 1.0,
        };
    }

    // Low-biomass mode has no automatic chunk-breadth filtering, and
    // therefore nothing to scale.
    if cfg.lowbiomass {
        return CommunityChunkBreadthThresholds {
            automatic: 0.0,
            standard: 0.0,
            safeguard: 0.0,
            user_scale: 1.0,
            recall_scale: 1.0,
        };
    }

    let automatic = community_chunk_breadth_cutoff(basis_reads, cfg.sequencer.is_illumina());
    let safeguard = scale_automatic_chunk_breadth(automatic, cfg.chunk_breadth_scale);
    let recall_scale = automatic_recall_chunk_breadth_scale(cfg, basis_reads);
    let standard = if recall_scale == 1.0 {
        safeguard
    } else {
        safeguard * recall_scale
    };
    CommunityChunkBreadthThresholds {
        automatic,
        standard,
        safeguard,
        user_scale: cfg.chunk_breadth_scale,
        recall_scale,
    }
}

/// Scale only an automatic cutoff in recall mode. The short-read
/// curve is continuous at one million reads and approaches 0.70 smoothly:
/// `scale = 0.70 + 0.30 / m`, where `m = basis_reads / 1e6`. Long reads use
/// a fixed 0.80 multiplier above the same depth gate.
#[inline]
fn automatic_recall_chunk_breadth_scale(cfg: &ProfilerConfig, basis_reads: usize) -> f64 {
    if cfg.mode != ProfileMode::Recall
        || cfg.lowbiomass
        || cfg.chunk_breadth.is_some()
        || basis_reads <= RECALL_CHUNK_SCALE_MIN_BASIS
    {
        return 1.0;
    }

    if cfg.sequencer.is_illumina() {
        let millions = basis_reads as f64 / 1_000_000.0;
        (SHORT_RECALL_CHUNK_SCALE_FLOOR + (1.0 - SHORT_RECALL_CHUNK_SCALE_FLOOR) / millions)
            .clamp(SHORT_RECALL_CHUNK_SCALE_FLOOR, 1.0)
    } else {
        LONG_RECALL_CHUNK_SCALE
    }
}

fn community_chunk_breadth_source(
    cfg: &ProfilerConfig,
    plan: ChunkBreadthPlan,
    thresholds: CommunityChunkBreadthThresholds,
) -> &'static str {
    if cfg.chunk_breadth.is_some() {
        "user-supplied --min-chunk-breadth"
    } else if cfg.lowbiomass {
        "disabled by --low-biomass"
    } else {
        match (thresholds.is_scaled(), plan.scaled, plan.use_aligned_basis) {
            (true, true, _) => "recall-scaled auto-estimate from scaled mapped reads",
            (true, false, true) => "recall-scaled auto-estimate from mapped reads",
            (true, false, false) => "recall-scaled auto-estimate from total reads",
            (false, true, _) => "auto-estimated from scaled mapped reads",
            (false, false, true) => "auto-estimated from mapped reads",
            (false, false, false) => "auto-estimated from total reads",
        }
    }
}

fn log_chunk_breadth_scaling(thresholds: CommunityChunkBreadthThresholds) {
    if thresholds.is_user_scaled() {
        let uncapped = thresholds.automatic * thresholds.user_scale;
        if uncapped > 1.0 {
            log::info!(
                target: "PROFILER",
                "User chunk-breadth scale adjusted the automatic cutoff: {:.6} x {} = {:.6}, capped at 1.000000.",
                thresholds.automatic,
                thresholds.user_scale,
                thresholds.safeguard,
            );
        } else {
            log::info!(
                target: "PROFILER",
                "User chunk-breadth scale adjusted the automatic cutoff: {:.6} x {} = {:.6}.",
                thresholds.automatic,
                thresholds.user_scale,
                thresholds.safeguard,
            );
        }
    }
    if thresholds.is_recall_scaled() {
        log::info!(
            target: "PROFILER",
            "Recall mode scaled the automatic chunk-breadth cutoff: {:.6} x {:.4} = {:.6}; OEBR/COEBR thresholds are unchanged.",
            thresholds.safeguard,
            thresholds.recall_scale,
            thresholds.standard,
        );
    }
}

/// Scale the read count used as input to the automatic chunk-breadth curve.
#[inline]
fn scale_low_map_chunk_basis(chunk_basis_reads: usize, scale: bool) -> usize {
    if scale {
        ((chunk_basis_reads as f64) * LOW_MAP_CHUNK_BASIS_SCALE).round() as usize
    } else {
        chunk_basis_reads
    }
}

/// Derive the genus taxid for a species taxid using the NCBI taxonomy
/// lineage (no header parsing).
///
/// Returns `None` if the species has no ancestor at rank "genus".
fn resolve_genus_taxid(taxonomy: &Taxonomy, species_taxid: u32) -> Option<u32> {
    taxonomy
        .get_parents(species_taxid)
        .into_iter()
        .find(|(_, _, rank)| rank == "genus")
        .map(|(taxid, _, _)| taxid)
}

/// Resolve a reference taxid to its species ancestor when one is available.
fn resolve_species_taxid(taxonomy: &Taxonomy, taxid: u32) -> Option<u32> {
    if matches!(taxonomy.get_rank(taxid), Some("species")) {
        return Some(taxid);
    }
    taxonomy
        .get_parents(taxid)
        .into_iter()
        .find(|(_, _, rank)| rank == "species")
        .map(|(id, _, _)| id)
}

/// Depth-scaled minimum chunk-breadth cutoff of the species amount test.
///
/// Both platforms use a bounded, saturating Michaelis-Menten form in
/// `m = num_reads / 1e6`:
///
/// - Illumina: `0.66*m/(m+1) + 0.24*m/(m+100)`.
/// - Nanopore/PacBio: `0.85*m/(m+0.78)`.
///
/// The cutoff rises smoothly with sequencing depth and then plateaus. The
/// 0.95 clamp is defensive; both fitted plateaus are lower.
fn community_chunk_breadth_cutoff(num_reads: usize, illumina: bool) -> f64 {
    let m = num_reads as f64 / 1_000_000.0;
    let value = if illumina {
        0.66 * m / (m + 1.0) + 0.24 * m / (m + 100.0)
    } else {
        0.85 * m / (m + 0.78)
    };
    value.min(0.95)
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Aggregated data from mapping loading.
///
/// Contains all information needed for the profiling phase.
struct MappingData {
    /// Per-read alignment bundles
    read_alignments: Vec<ReadBundle>,

    /// Chunk-breadth cutoff of the species amount test. In deep recall-mode
    /// runs this may be lower than the original automatic cutoff.
    min_chunk_breadth: f64,

    /// Whether at least one genome has a full/effective size multiplier > 1.
    has_scaled_genomes: bool,

    /// Total reads in the input.
    total_reads: usize,

    /// The genus of every species with whole-genome references, also on a
    /// fractional database. Species represented only by marker genes are
    /// absent.
    genus_of_species: Arc<HashMap<u32, u32>>,

    /// Most-specific reference-supported rank for every terminal taxid. Other
    /// rows (genus roll-up, family level) take their rank from the taxonomy.
    terminal_ranks: Arc<HashMap<u32, ProfileRank>>,

    /// Marker genes and fragment geometry when the database has marker genomes.
    marker_model: Option<Arc<MarkerModel>>,

    /// Inputs of the genus roll-up; `None` when it is off (`--strict-coverage`).
    genus_rollup: Option<Arc<GenusRollup>>,

    /// Inputs of `--estimate-unclassified`; `None` when it is off.
    unclassified_species: Option<MarkerSpecies>,

    /// Genome accession of each assembly id.
    assembly_names: Vec<Arc<str>>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profiler::evidence::{EvidenceMode, EvidencePolicy};
    use crate::profiler::ProfileMode;
    use crate::sequencer::Sequencer;

    fn approx_equal(lhs: f64, rhs: f64) {
        assert!(
            (lhs - rhs).abs() < 1e-9,
            "left {} not approximately equal to right {}",
            lhs,
            rhs
        );
    }

    fn placement(assembly_id: u32, mapped: u32, read_len: u32) -> CandidatePlacement {
        CandidatePlacement {
            taxid: 1,
            profile_rank: ProfileRank::Species,
            reference_id: 0,
            assembly_id,
            sampling_scale: 1.0,
            asm_len: 2_000.0,
            r1_mapped: mapped,
            both_mapped: mapped,
            read_len,
            other_genes: 0,
            pair_errors: 0.0,
            long_score: 0.0,
            span: (1, 1 + mapped),
            mate_span: (0, 0),
            fixed_chunk: (1, 1 + mapped),
            flex_chunk: (1, 1 + mapped),
            input_order: 0,
        }
    }

    #[test]
    fn a_long_read_adds_depth_for_every_marker_gene_it_covers() {
        // Long reads of 3,000 bases on genes of 1,000 bases: 500 bases inside
        // a gene, 1,000 - 1,000 + 3,001 = 3,001 placements per gene.
        let geometry = FragmentGeometry::new(FragmentLayout::LongReads, 3_000, 3_000, 500, 0.2);
        let markers = MarkerModel::from_gene_lengths(&[(7, &[1_000, 1_000])], geometry);
        let mut read = placement(7, 900, 3_000);
        approx_equal(fragment_depth(read, 900, Some(&markers)), 3_000.0 / 6_002.0);
        read.other_genes = 1;
        approx_equal(fragment_depth(read, 900, Some(&markers)), 6_000.0 / 6_002.0);
        // A read of another length is placed by its own length.
        let short = placement(7, 900, 1_000);
        approx_equal(
            fragment_depth(short, 900, Some(&markers)),
            1_000.0 / (2.0 * (1_000.0 + 1_001.0 - 400.0)),
        );
    }

    #[test]
    fn a_fragment_on_a_marker_gene_adds_its_read_over_the_genes_placements() {
        // 150 bp pairs of 400 bases with lone mates: genes of 1,000 bases
        // offer 1,191 placements each; a gene of 100 bases holds no read.
        let geometry = FragmentGeometry::new(FragmentLayout::PairsAndLoneMates, 150, 400, 50, 0.7);
        let markers =
            MarkerModel::from_gene_lengths(&[(7, &[1_000, 1_000]), (8, &[100])], geometry);
        let depth = |assembly, mapped, markers| {
            fragment_depth(placement(assembly, mapped, 150), mapped, markers)
        };
        approx_equal(depth(7, 150, Some(&markers)), 150.0 / 2_382.0);
        // A read clipped at a gene end stands for its whole read.
        approx_equal(depth(7, 110, Some(&markers)), 150.0 / 2_382.0);
        // Genome references and marker genes without placements keep mapped
        // bases over the summed length.
        approx_equal(depth(9, 140, Some(&markers)), 140.0 / 2_000.0);
        approx_equal(depth(8, 90, Some(&markers)), 90.0 / 2_000.0);
        approx_equal(depth(7, 140, None), 140.0 / 2_000.0);
    }

    #[test]
    fn expected_breadth_is_clamped_for_high_depths() {
        let expected = crate::profiler::evidence::expected_breadth(10.0, 5.0);
        approx_equal(expected, 1.0);
        assert!(expected.is_finite());

        let zero = crate::profiler::evidence::expected_breadth(10.0, 0.0);
        approx_equal(zero, 0.0);
    }

    // ------------------------------------------------------------------
    // Test helpers.
    // ------------------------------------------------------------------

    /// Minimal `ProfilerConfig` for unit tests. Callers mutate only the
    /// fields relevant to the scenario under test.
    fn make_cfg(sequencer: Sequencer, host: Option<&str>) -> ProfilerConfig {
        ProfilerConfig {
            genome_level: false,
            mapper_mappings: None,
            mapping_profile_memory: None,
            db_metadata: None,
            sequencer,
            max_em_iterations: 200,
            lenient_mapping_input: false,
            is_paired: false,
            identity: 0.95,
            identity_check_accuracy: None,
            noisy_long_reads: false,
            mapped_len: 50,
            breadth: Some(0.0),
            chunk_breadth: Some(0.0),
            chunk_breadth_scale: 1.0,
            min_reads: None,
            min_oebr: None,
            min_coebr: None,
            min_marker_genes: 3,
            relative_identity_gap: None,
            removed_to_genus: false,
            estimate_unclassified: false,
            identity_margin: None,
            fraction: 0.6,
            lowbiomass: false,
            save_decisions: false,
            save_classifications: false,
            by_aligned: false,
            family_identity: 0.80,
            low_map_rate_threshold: 0.30,
            pathogen_host: None,
            host: host.map(|s| s.to_string()),
            threads: 1,
            verbose: false,
            very_verbose: false,
            progress: crate::progress::ProgressOptions::default(),
            suppress_mapping_load_progress: false,
            outprefix: "test".to_string(),
            mode: ProfileMode::Default,
            dmp_dir: None,
            ignore_unresolved_references: false,
        }
    }

    fn test_taxonomy(directory: &Path) -> Arc<Taxonomy> {
        let nodes = directory.join("nodes.dmp");
        let names = directory.join("names.dmp");
        let merged = directory.join("merged.dmp");
        std::fs::write(
            &nodes,
            "1\t|\t1\t|\tno rank\t|\n101\t|\t1\t|\tspecies\t|\n102\t|\t1\t|\tspecies\t|\n",
        )
        .unwrap();
        std::fs::write(
            &names,
            "1\t|\troot\t|\t\t|\tscientific name\t|\n101\t|\tSpecies one\t|\t\t|\tscientific name\t|\n102\t|\tSpecies two\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();
        std::fs::write(&merged, "").unwrap();
        Arc::new(Taxonomy::from_files(&nodes, &names, &merged, None).unwrap())
    }

    #[test]
    fn no_classifications_are_written_for_a_fractional_database() {
        let directory = tempfile::tempdir().unwrap();
        let taxonomy = test_taxonomy(directory.path());
        for fractional in [false, true] {
            let prefix = directory.path().join(format!("fractional-{fractional}"));
            let mut cfg = make_cfg(Sequencer::Illumina, None);
            cfg.save_classifications = true;
            cfg.outprefix = prefix.to_string_lossy().into_owned();
            cfg.db_metadata = Some(Arc::new(DatabaseMetadata {
                metam_index: std::path::PathBuf::from("metax_db.mxi"),
                references: Vec::new(),
                is_fractional: fractional,
            }));
            let profiler = CommunityProfiler::new(cfg, Arc::clone(&taxonomy)).unwrap();
            let mapping_data = MappingData {
                read_alignments: Vec::new(),
                min_chunk_breadth: 0.05,
                has_scaled_genomes: false,
                total_reads: 2,
                genus_of_species: Arc::new(HashMap::new()),
                terminal_ranks: Arc::new(HashMap::new()),
                marker_model: None,
                genus_rollup: None,
                unclassified_species: None,
                assembly_names: Vec::new(),
            };
            profiler.profile(mapping_data).unwrap();
            let written =
                Path::new(&format!("{}.classify.txt.gz", prefix.to_string_lossy())).exists();
            assert_eq!(written, !fractional);
        }
    }

    #[test]
    fn conflicting_terminal_rank_for_one_taxid_is_rejected() {
        let mut ranks = HashMap::new();
        record_terminal_rank(&mut ranks, 111, ProfileRank::Genus, "genus-ref").unwrap();
        record_terminal_rank(&mut ranks, 111, ProfileRank::Genus, "second-genus-ref").unwrap();
        let error =
            record_terminal_rank(&mut ranks, 111, ProfileRank::Species, "species-ref").unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("conflicting terminal ranks for taxid 111"));
        assert!(message.contains("genus"));
        assert!(message.contains("species"));
    }

    #[test]
    fn equivalent_mapping_rate_counts_the_genome_share_of_fragments() {
        let bundle = |scale: f32| {
            let mut row = native_mapping(1, 10, 1, 0.99, 100, 100, MappingRole::Unpaired);
            row.sampling_scale = scale;
            ReadBundle {
                read_name: Arc::from("read"),
                alignments: vec![row],
            }
        };
        assert_eq!(
            genome_equivalent_mapping_rate(&[bundle(1.0), bundle(1.0)], 4, None),
            0.5
        );
        assert_eq!(
            genome_equivalent_mapping_rate(&[bundle(10.0)], 20, None),
            0.5
        );
        assert_eq!(
            genome_equivalent_mapping_rate(&[bundle(10.0)], 4, None),
            1.0
        );
        assert_eq!(genome_equivalent_mapping_rate(&[], 0, None), 0.0);
    }

    fn native_mapping(
        species: u32,
        reference_id: u32,
        assembly_id: u32,
        identity: f32,
        mapped_len: u32,
        read_len: u32,
        role: MappingRole,
    ) -> AlignmentEntry {
        native_mapping_at_rank(
            species,
            ProfileRank::Species,
            reference_id,
            assembly_id,
            identity,
            mapped_len,
            read_len,
            role,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn native_mapping_at_rank(
        species: u32,
        profile_rank: ProfileRank,
        reference_id: u32,
        assembly_id: u32,
        identity: f32,
        mapped_len: u32,
        read_len: u32,
        role: MappingRole,
    ) -> AlignmentEntry {
        AlignmentEntry {
            profile_taxid: species,
            profile_rank,
            reference_id,
            assembly_id,
            role,
            pair_id: if matches!(role, MappingRole::Mate1 | MappingRole::Mate2) {
                reference_id.max(1)
            } else {
                0
            },
            sampling_scale: 1.0,
            asm_len: 1_000.0,
            identity,
            fraction: mapped_len as f32 / read_len as f32,
            mapped_len,
            read_len,
            span: (100, 100 + mapped_len),
            fixed_chunk: (1, 300 + mapped_len),
            flex_chunk: (1, 500 + mapped_len),
            marker_gene_len: 0,
        }
    }

    #[test]
    fn evidence_policy_is_limited_to_supported_inputs() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        assert_eq!(
            EvidencePolicy::from_config(&cfg).map(|policy| policy.mode),
            Some(EvidenceMode::Single)
        );

        cfg.is_paired = true;
        assert_eq!(
            EvidencePolicy::from_config(&cfg).map(|policy| policy.mode),
            Some(EvidenceMode::Paired)
        );

        cfg.is_paired = false;
        cfg.sequencer = Sequencer::Nanopore;
        assert_eq!(
            EvidencePolicy::from_config(&cfg).map(|policy| policy.mode),
            Some(EvidenceMode::Long)
        );

        cfg.is_paired = true;
        assert!(EvidencePolicy::from_config(&cfg).is_none());
    }

    #[test]
    fn automatic_long_read_identity_follows_reads_that_align_worse() {
        let unpaired = MappingRole::Unpaired;
        let reads: Vec<ReadBundle> = [0.80_f32, 0.84, 0.88, 0.89]
            .iter()
            .enumerate()
            .map(|(i, &identity)| ReadBundle {
                read_name: format!("r{i}").into(),
                alignments: vec![
                    native_mapping(101, 10, 1, identity, 1_000, 2_000, unpaired),
                    // Higher identity over fewer bases is not the read's best record.
                    native_mapping(102, 20, 2, identity + 0.05, 600, 2_000, unpaired),
                ],
            })
            .collect();
        let resolved = |claimed: Option<f64>, identity: f64| {
            let mut cfg = make_cfg(Sequencer::Nanopore, None);
            cfg.mapped_len = 500;
            cfg.fraction = 0.20;
            cfg.identity = identity;
            cfg.identity_check_accuracy = claimed;
            let taxonomy = Taxonomy::from_parts(
                HashMap::new(),
                HashMap::new(),
                HashMap::new(),
                HashSet::new(),
            );
            let mut profiler = CommunityProfiler::new(cfg, Arc::new(taxonomy)).unwrap();
            profiler.lower_identity_to_alignments(&reads);
            profiler.cfg.identity
        };
        // Upper quartile 0.89 below the claimed 0.94: cutoff 0.89 - 0.03.
        assert!((resolved(Some(0.94), 0.91) - 0.86).abs() < 1e-6);
        // Reads align as claimed, the cutoff would rise, or it is not automatic.
        assert_eq!(resolved(Some(0.88), 0.85), 0.85);
        assert_eq!(resolved(Some(0.99), 0.80), 0.80);
        assert_eq!(resolved(None, 0.91), 0.91);
    }

    // ------------------------------------------------------------------
    // Low-map gate tests
    // ------------------------------------------------------------------

    #[test]
    fn low_map_gate_activates_only_below_mapping_rate_threshold() {
        let cfg = make_cfg(Sequencer::Illumina, None);
        assert!(low_map_gate_active(&cfg, 0.0));
        assert!(low_map_gate_active(&cfg, 0.29));
        assert!(low_map_gate_active(&cfg, 0.2999999));
        // Exactly at the 0.30 threshold is NOT active (strict `<`).
        assert!(!low_map_gate_active(&cfg, 0.30));
        assert!(!low_map_gate_active(&cfg, 0.50));
        assert!(!low_map_gate_active(&cfg, 1.0));
    }

    #[test]
    fn low_map_gate_disabled_when_host_is_set() {
        // Any host string (e.g. human 9606) disables the low-map gate regardless
        // of sequencer or mapping rate.
        let cfg_host = make_cfg(Sequencer::Illumina, Some("9606"));
        assert!(!low_map_gate_active(&cfg_host, 0.0));
        assert!(!low_map_gate_active(&cfg_host, 0.1));
    }

    #[test]
    fn low_map_gate_covers_long_reads_but_not_assembly() {
        // Nanopore and PacBio honour the low-map auto-trigger too.
        for seq in [Sequencer::Nanopore, Sequencer::PacBio] {
            let cfg = make_cfg(seq, None);
            assert!(
                low_map_gate_active(&cfg, 0.0),
                "expected the gate active for low-map {seq}"
            );
            assert!(
                !low_map_gate_active(&cfg, 0.99),
                "expected the gate inactive for high-map {seq}"
            );
        }
    }

    #[test]
    fn low_map_gate_respects_custom_low_map_rate_threshold() {
        // With the default 0.30 threshold, 0.40 is NOT low-map.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        assert!(!low_map_gate_active(&cfg, 0.40));
        // Raise the threshold: now 0.40 IS low-map (strict `<`).
        cfg.low_map_rate_threshold = 0.50;
        assert!(low_map_gate_active(&cfg, 0.40));
        assert!(low_map_gate_active(&cfg, 0.4999));
        assert!(!low_map_gate_active(&cfg, 0.50));
        // Lower the threshold: even 0.05 no longer triggers when = threshold.
        cfg.low_map_rate_threshold = 0.05;
        assert!(!low_map_gate_active(&cfg, 0.05));
        assert!(low_map_gate_active(&cfg, 0.0499));
    }

    #[test]
    fn resolve_by_aligned_honors_explicit_flag() {
        // --by-mapped flag present: mapped basis is forced regardless of
        // sequencer, host, or mapping rate.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.by_aligned = true;
        assert!(cfg.resolve_by_aligned(0.99));
        cfg.sequencer = Sequencer::Nanopore;
        assert!(cfg.resolve_by_aligned(0.99));
        cfg.host = Some("9606".to_string());
        assert!(cfg.resolve_by_aligned(0.99));
    }

    #[test]
    fn resolve_by_aligned_auto_forces_for_low_map_read_platforms_non_pathogen() {
        // Default cfg has --by-mapped unset (false); auto-rule applies.
        let cfg = make_cfg(Sequencer::Illumina, None);
        assert!(!cfg.by_aligned);
        // Below default 0.30 threshold => aligned basis is auto-forced.
        assert!(cfg.resolve_by_aligned(0.10));
        // At or above threshold => total-read basis.
        assert!(!cfg.resolve_by_aligned(0.30));
        assert!(!cfg.resolve_by_aligned(0.50));
        // Pathogen mode on disables auto-force.
        let cfg_host = make_cfg(Sequencer::Illumina, Some("9606"));
        assert!(!cfg_host.resolve_by_aligned(0.10));
        // Long-read platforms also auto-force on low map.
        let cfg_np = make_cfg(Sequencer::Nanopore, None);
        assert!(cfg_np.resolve_by_aligned(0.10));
    }

    #[test]
    fn scale_low_map_chunk_basis_applies_only_when_scale_true() {
        // Pass-through when scaling is off.
        assert_eq!(scale_low_map_chunk_basis(0, false), 0);
        assert_eq!(scale_low_map_chunk_basis(1, false), 1);
        assert_eq!(scale_low_map_chunk_basis(100_000, false), 100_000);

        // When scaling is on, count is multiplied by LOW_MAP_CHUNK_BASIS_SCALE
        // and rounded to the nearest integer.
        assert_eq!(scale_low_map_chunk_basis(0, true), 0);
        assert_eq!(scale_low_map_chunk_basis(1, true), 2); // 1.5 -> 2
        assert_eq!(scale_low_map_chunk_basis(2, true), 3);
        assert_eq!(scale_low_map_chunk_basis(10, true), 15);
        assert_eq!(scale_low_map_chunk_basis(100, true), 150);
        assert_eq!(scale_low_map_chunk_basis(100_000, true), 150_000);
        // Round-half-to-even-ish: 3 * 1.5 = 4.5 => rounds to 5.
        assert_eq!(scale_low_map_chunk_basis(3, true), 5);
    }

    #[test]
    fn low_map_chunk_basis_scale_constant_is_one_and_a_half() {
        assert!((LOW_MAP_CHUNK_BASIS_SCALE - 1.5).abs() < 1e-12);
    }

    #[test]
    fn chunk_breadth_plan_scales_only_low_map_aligned_basis() {
        let total_reads = 100_000usize;
        let aligned_reads = 10_000usize;

        // --by-mapped alone: mapped basis, same equation, no 1.5x.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.by_aligned = true;
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // Automatic low-map gate on a full DB: aligned basis with 1.5x.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        let gate = low_map_gate_active(&cfg, 0.10);
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, false, gate);
        assert!(plan.use_aligned_basis);
        assert!(plan.by_aligned_auto);
        assert!(plan.scaled);
        assert_eq!(plan.basis_reads, aligned_reads);
        assert_eq!(plan.estimate_reads, 15_000);

        // Without the low-map gate and --by-mapped: total basis, no scale.
        let plan = plan_chunk_breadth(&cfg, 0.90, total_reads, aligned_reads, false, false);
        assert!(!plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.basis_reads, total_reads);
        assert_eq!(plan.estimate_reads, total_reads);

        // A low-map subsampled run takes the mapped basis without scaling.
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        let gate = low_map_gate_active(&cfg, 0.10);
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, true, gate);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.estimate_reads, aligned_reads);

        // User-pinned --min-chunk-breadth is never rescaled.
        cfg.chunk_breadth = Some(0.25);
        let plan = plan_chunk_breadth(&cfg, 0.10, total_reads, aligned_reads, false, gate);
        assert!(plan.use_aligned_basis);
        assert!(!plan.scaled);
        assert_eq!(plan.estimate_reads, aligned_reads);
    }

    #[test]
    fn low_map_gate_active_for_illumina_single_and_paired() {
        // Single-end
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.is_paired = false;
        assert!(low_map_gate_active(&cfg, 0.1));
        // Paired-end
        cfg.is_paired = true;
        assert!(low_map_gate_active(&cfg, 0.1));
    }

    #[test]
    fn community_chunk_breadth_curve_matches_spec() {
        // Short (Illumina): 0.66*m/(m+1) + 0.24*m/(m+100), plateau 0.90.
        let short = |n| community_chunk_breadth_cutoff(n, true);
        assert!((short(0) - 0.0).abs() < 1e-12);
        assert!((short(10_000) - 0.00656).abs() < 1e-4);
        assert!((short(100_000) - 0.06024).abs() < 1e-4);
        assert!((short(1_000_000) - 0.33238).abs() < 1e-4); // C = 0.33 at 1M reads
        assert!((short(3_300_000) - 0.51418).abs() < 1e-4);
        assert!((short(16_600_000) - 0.65667).abs() < 1e-4);
        // approaches the 0.90 plateau at extreme depth, always below the clamp
        let short_deep = short(100_000_000_000);
        assert!(short_deep > 0.899 && short_deep < 0.95);

        // Long (Nanopore/PacBio): 0.85*m/(m+0.78), plateau 0.85.
        let long = |n| community_chunk_breadth_cutoff(n, false);
        assert!((long(0) - 0.0).abs() < 1e-12);
        assert!((long(10_000) - 0.01076).abs() < 1e-4);
        assert!((long(100_000) - 0.09659).abs() < 1e-4);
        assert!((long(1_000_000) - 0.47753).abs() < 1e-4);
        assert!((long(3_300_000) - 0.68750).abs() < 1e-4);
        assert!((long(16_600_000) - 0.81185).abs() < 1e-4);
        let long_deep = long(100_000_000_000);
        assert!(long_deep > 0.849 && long_deep < 0.85);

        // Monotonic non-decreasing in read count on both platforms.
        for &(a, b) in &[
            (1_000usize, 10_000usize),
            (100_000, 1_000_000),
            (1_000_000, 16_600_000),
        ] {
            assert!(
                community_chunk_breadth_cutoff(b, true) >= community_chunk_breadth_cutoff(a, true)
            );
            assert!(
                community_chunk_breadth_cutoff(b, false)
                    >= community_chunk_breadth_cutoff(a, false)
            );
        }
    }

    #[test]
    fn recall_scales_only_deep_automatic_community_chunk_cutoffs() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.mode = ProfileMode::Recall;

        let at_gate = resolve_community_chunk_breadth(&cfg, 1_000_000);
        approx_equal(at_gate.recall_scale, 1.0);
        approx_equal(at_gate.standard, at_gate.safeguard);
        assert_eq!(at_gate.standard.to_bits(), at_gate.safeguard.to_bits());

        let two_million = resolve_community_chunk_breadth(&cfg, 2_000_000);
        approx_equal(two_million.recall_scale, 0.85);
        approx_equal(
            two_million.standard,
            two_million.safeguard * two_million.recall_scale,
        );

        let five_million = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(five_million.recall_scale, 0.76);
        let deep = resolve_community_chunk_breadth(&cfg, 16_600_000);
        approx_equal(deep.recall_scale, 0.70 + 0.30 / 16.6);
        assert!(deep.recall_scale > SHORT_RECALL_CHUNK_SCALE_FLOOR);

        cfg.sequencer = Sequencer::Nanopore;
        let long_at_gate = resolve_community_chunk_breadth(&cfg, 1_000_000);
        approx_equal(long_at_gate.recall_scale, 1.0);
        let long_deep = resolve_community_chunk_breadth(&cfg, 1_000_001);
        approx_equal(long_deep.recall_scale, LONG_RECALL_CHUNK_SCALE);
        approx_equal(
            long_deep.standard,
            long_deep.safeguard * LONG_RECALL_CHUNK_SCALE,
        );

        cfg.sequencer = Sequencer::PacBio;
        approx_equal(
            resolve_community_chunk_breadth(&cfg, 5_000_000).recall_scale,
            LONG_RECALL_CHUNK_SCALE,
        );
    }

    #[test]
    fn recall_chunk_scaling_respects_user_and_nonstandard_policy_boundaries() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.mode = ProfileMode::Recall;
        cfg.chunk_breadth = Some(0.37);
        let explicit = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(explicit.standard, 0.37);
        approx_equal(explicit.safeguard, 0.37);
        approx_equal(explicit.recall_scale, 1.0);
        assert_eq!(explicit.standard.to_bits(), explicit.safeguard.to_bits());

        cfg.chunk_breadth = None;
        cfg.lowbiomass = true;
        let low_biomass = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(low_biomass.standard, 0.0);
        approx_equal(low_biomass.safeguard, 0.0);
        approx_equal(low_biomass.recall_scale, 1.0);

        // A read minimum adds to the automatic cutoff and leaves it unchanged.
        cfg.lowbiomass = false;
        let without_min_reads = resolve_community_chunk_breadth(&cfg, 20_000_000);
        cfg.min_reads = Some(3);
        assert_eq!(
            resolve_community_chunk_breadth(&cfg, 20_000_000),
            without_min_reads
        );

        cfg.min_reads = None;
        cfg.mode = ProfileMode::Default;
        let default = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(default.recall_scale, 1.0);
        assert_eq!(default.standard.to_bits(), default.safeguard.to_bits());
        cfg.mode = ProfileMode::Precision;
        let precision = resolve_community_chunk_breadth(&cfg, 20_000_000);
        approx_equal(precision.recall_scale, 1.0);
        assert_eq!(precision.standard.to_bits(), precision.safeguard.to_bits());
    }

    #[test]
    fn user_chunk_scale_adjusts_automatic_standard_and_safeguard_cutoffs() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;

        let unchanged = resolve_community_chunk_breadth(&cfg, 5_000_000);
        assert_eq!(unchanged.automatic.to_bits(), unchanged.safeguard.to_bits());
        assert_eq!(unchanged.automatic.to_bits(), unchanged.standard.to_bits());

        cfg.chunk_breadth_scale = 0.5;
        let relaxed = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(relaxed.safeguard, relaxed.automatic * 0.5);
        assert_eq!(relaxed.standard.to_bits(), relaxed.safeguard.to_bits());

        cfg.chunk_breadth_scale = 10.0;
        let capped = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(capped.safeguard, 1.0);
        approx_equal(capped.standard, 1.0);
    }

    #[test]
    fn recall_scaling_composes_after_user_chunk_scaling() {
        let mut cfg = make_cfg(Sequencer::Illumina, None);
        cfg.chunk_breadth = None;
        cfg.chunk_breadth_scale = 0.5;
        cfg.mode = ProfileMode::Recall;

        let thresholds = resolve_community_chunk_breadth(&cfg, 5_000_000);
        approx_equal(thresholds.user_scale, 0.5);
        approx_equal(thresholds.recall_scale, 0.76);
        approx_equal(thresholds.safeguard, thresholds.automatic * 0.5);
        approx_equal(
            thresholds.standard,
            thresholds.safeguard * thresholds.recall_scale,
        );
    }
}

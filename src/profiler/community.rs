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
//!    cutoff of the species amount test, the inputs of the levels above
//!    species and of `--estimate-unclassified`, and the marker-gene model.
//!
//! 4. **Profiling**: the level-wise engine ([`levels`]) judges species, the
//!    genomes within them and the levels above, and writes the outputs.

use std::collections::{hash_map::Entry, HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;
use std::str;
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};

use crate::db::{DatabaseMetadata, ReferenceType};
use crate::mapper::tsv::{
    for_each_metam_mapping_record, for_each_structured_mapping_record, load_shared_metam_profile,
    open_metam_mapping_reader, total_input_units, MetamMappingReceiver, MetamMappingRecord,
    MetamParseStats, MetamProfile,
};
use crate::progress;
use crate::taxonomy::Taxonomy;

use super::evidence::{EvidenceAlignment, MappingRole};
use super::markers::{FragmentGeometry, FragmentLayout, LengthHistogram, MarkerModel};
use super::rollup::{AboveSpeciesInputs, GenomeLayout};
use super::unclassified::{self, MarkerSpecies, ReferenceKind};
use super::{
    coverage_passes, log_reference_taxonomy_adjustments, log_resolved_profiler_parameters,
    parse_reference_taxonomy, push_string_keyed_value, record_genome_sampling_scale,
    record_reference_taxonomy_adjustment, sampling_scale_from_sizes, scale_automatic_chunk_breadth,
    strip_pair_suffix, ProfileMode, ProfileRank, Profiler, ProfilerConfig, ReferenceTaxonomy,
    ReferenceTaxonomyAdjustments, ReferenceTaxonomyResolver, ResolvedProfilerParameters,
    StableReadMap, UnresolvedReferences,
};

/// Fixed coordinate radius around each read span for fixed-chunk coverage.
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

/// Narrow a genomic coordinate/length to `u32` for the packed
/// [`EvidenceAlignment`] fields. Coordinates are bounded by `ref_len + 1`, and the
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
    read_name: Box<str>,

    /// All alignments of this read, in the order of the mapping records.
    alignments: Vec<EvidenceAlignment>,
}

#[derive(Clone, Debug)]
/// A database reference as the mapping loader needs it; values shared by its
/// genome are in [`MetamAssembly`].
struct MetamReferenceEntry {
    ref_len: usize,
    profile_taxid: u32,
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

/// Loading the mappings of a sample.
impl Profiler {
    /// The saved or staged mapping file's mappings.
    pub(super) fn load_mapping_file(&mut self) -> Result<MappingData> {
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
        self.load_metam(Path::new(&mappings), &db)
    }

    /// The mappings streamed from the mapper.
    pub(super) fn load_mapping_stream(
        &mut self,
        receiver: MetamMappingReceiver,
        source: &str,
    ) -> Result<MappingData> {
        let db = Arc::clone(
            self.cfg
                .db_metadata
                .as_ref()
                .ok_or_else(|| anyhow!("missing Metax database metadata for streamed mappings"))?,
        );
        self.load_metam_mappings(receiver, source, &db)
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
        let score = |aln: &&EvidenceAlignment| f64::from(aln.identity) * f64::from(aln.mapped_len);
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
        let mut strain_taxids = vec![None; assembly_names.len()];
        let mut previous_genome = None;
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
            // The strain that a genome's references declare, read where they
            // start; the dump may not know it.
            let assembly_id = assembly_ids[&*reference.genome];
            if previous_genome != Some(assembly_id) {
                if let Some(strain) = declared_taxonomy.primary_taxid {
                    strain_taxids[assembly_id as usize] = Some(strain);
                }
                previous_genome = Some(assembly_id);
            }
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
        // The species represented only by marker genomes, with those genomes.
        let marker_species = if has_marker_genomes {
            unclassified::marker_species(
                db.references.iter().zip(&reference_taxonomies).filter_map(
                    |(reference, taxonomy)| {
                        Some(ReferenceKind {
                            taxid: taxonomy.as_ref()?.profile_taxid,
                            assembly_id: assembly_ids[&*reference.genome],
                            reference_type: reference.reference_type,
                            declared_genome_size: reference.declared_genome_size,
                        })
                    },
                ),
            )
        } else {
            MarkerSpecies::new()
        };
        let mut genus_of_species = self.genus_of_taxa(&profile_taxids);
        let marker_only = if self.cfg.above_species {
            marker_species
                .keys()
                .map(|&species| (species, genus_of_species.get(&species).copied()))
                .collect()
        } else {
            HashMap::new()
        };
        if has_marker_genomes {
            // A marker gene holds too little of a genome for its reads to
            // stand for a genus or family: reads failing every species there
            // come mostly from the homologous genes of relatives. Only species
            // with genome-type references have a genus here.
            genus_of_species.retain(|species, _| taxids_with_genomes.contains(species));
        }
        let above_species = self.cfg.above_species.then(|| AboveSpeciesInputs {
            marker_only,
            layout: GenomeLayout::new(
                db.references
                    .iter()
                    .map(|reference| (assembly_ids[&*reference.genome], reference.length as u64)),
            ),
        });
        let assemblies = assembly_names
            .iter()
            .map(|asm| MetamAssembly {
                asm_len: asm_lengths_raw.get(&**asm).copied().unwrap_or(0.0),
                flex_len: asm_lengths_raw
                    .get(&**asm)
                    .map_or(NATIVE_FIXED_FLANK_LEN, |len| {
                        (len.sqrt().round() as usize).max(1)
                    }),
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
                profile_taxid: taxonomy.profile_taxid,
                assembly_id,
                is_marker: reference.reference_type == ReferenceType::Markers,
            }));
        }

        let mut read_map: StableReadMap<EvidenceAlignment> = StableReadMap::default();
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
            if let Some((read_id, entry)) = self.entry_from_metam_mapping(mapping, &metam_refs) {
                push_string_keyed_value(&mut read_map, read_id, entry);
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
                read_name: Box::<str>::from(read_name),
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
        let basis = ChunkBasis::of(
            &self.cfg,
            gate_rate,
            (total_reads, aligned_reads),
            db.is_fractional,
        );
        let chunk_cutoff = ChunkCutoff::resolve(&self.cfg, basis.reads);

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
        match chunk_cutoff.source {
            ChunkSource::LowBiomass => log::info!(
                target: "PROFILER",
                "Automatic minimum chunk breadth filter disabled by --low-biomass."
            ),
            ChunkSource::Automatic if basis.mapped => log::info!(
                target: "PROFILER",
                "Minimum chunk breadth estimated from mapped {} ({}).",
                mapped_unit,
                aligned_reads
            ),
            ChunkSource::Automatic => log::info!(
                target: "PROFILER",
                "Minimum chunk breadth estimated from total reads ({}).",
                total_reads
            ),
            ChunkSource::User => {}
        }
        chunk_cutoff.log_scaling();
        log_resolved_profiler_parameters(
            &self.cfg,
            ResolvedProfilerParameters {
                min_breadth: self.cfg.breadth.unwrap_or(0.0),
                min_chunk_breadth: chunk_cutoff.cutoff,
                chunk_breadth_source: chunk_cutoff.source_in_words(&basis),
                chunk_breadth_basis_reads: basis.reads,
                total_reads,
                mapped_reads: aligned_reads,
                mapping_rate,
            },
        );

        Ok(MappingData {
            read_alignments,
            min_chunk_breadth: chunk_cutoff.cutoff,
            has_scaled_genomes: genome_sampling_scales.values().any(|scale| *scale > 1.0),
            total_reads,
            genus_of_species,
            terminal_ranks,
            marker_model,
            fractional: db.is_fractional,
            above_species,
            unclassified_species: self.cfg.estimate_unclassified.then_some(marker_species),
            assembly_names,
            strain_taxids,
        })
    }

    fn entry_from_metam_mapping<'a>(
        &self,
        mapping: MetamMappingRecord<'a>,
        metam_refs: &MetamReferences,
    ) -> Option<(&'a str, EvidenceAlignment)> {
        let (reference_id, reference, assembly) = metam_refs.get(mapping.ref_id)?;
        if reference.ref_len == 0 || mapping.read_len == 0 || mapping.mapped_bases == 0 {
            return None;
        }
        let start = mapping.start_1based.clamp(1, reference.ref_len.max(1));
        let mut end = mapping.end_1based_exclusive.max(start.saturating_add(1));
        end = end.min(reference.ref_len.saturating_add(1));
        let midpoint = start + ((end.saturating_sub(start)) / 2);
        let fixed_chunk_start = start.saturating_sub(NATIVE_FIXED_FLANK_LEN).max(1);
        let fixed_chunk_end = end
            .saturating_add(NATIVE_FIXED_FLANK_LEN)
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
            EvidenceAlignment {
                profile_taxid: reference.profile_taxid,
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

    /// The genus of each of the references' profile taxa (`taxids`, at species
    /// rank or above) that has one: the taxon itself at genus rank, or its
    /// ancestor there.
    fn genus_of_taxa(&self, taxids: &HashSet<u32>) -> HashMap<u32, u32> {
        let taxonomy = &self.taxonomy;
        taxids
            .iter()
            .filter_map(|&taxid| {
                debug_assert!(
                    resolve_species_taxid(taxonomy, taxid).is_none_or(|species| species == taxid),
                    "profile taxid {taxid} lies below species"
                );
                resolve_genus_taxid(taxonomy, taxid).map(|genus| (taxid, genus))
            })
            .collect()
    }
}

/// Mapped fragments as a share of the input fragments, each counting as the
/// share of a genome it represents: the sampling scale of the genome of its
/// best placement, so that a fractional database is not mistaken for a
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

// ============================================================================
// The chunk-breadth cutoff
// ============================================================================

/// The read count that the automatic chunk-breadth curve is read at.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ChunkBasis {
    /// Mapped reads (`--by-mapped`, or a low-map run), otherwise total reads.
    mapped: bool,
    /// Mapped reads of a low-map run on a full database, without a given
    /// `--min-chunk-breadth`: multiplied by [`LOW_MAP_CHUNK_BASIS_SCALE`].
    scaled: bool,
    reads: usize,
}

impl ChunkBasis {
    /// The basis of a run with `mapping_rate` and `(total, mapped)` reads, on
    /// a `fractional` database or not.
    fn of(
        cfg: &ProfilerConfig,
        mapping_rate: f64,
        (total, mapped): (usize, usize),
        fractional: bool,
    ) -> Self {
        let low_map = cfg.is_low_map_run(mapping_rate);
        let by_mapped = cfg.by_aligned || low_map;
        let scaled = cfg.chunk_breadth.is_none() && low_map && !fractional;
        let reads = if by_mapped { mapped } else { total };
        Self {
            mapped: by_mapped,
            scaled,
            reads: scale_low_map_chunk_basis(reads, scaled),
        }
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

/// Where the chunk-breadth cutoff of the species amount test comes from.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ChunkSource {
    /// `--min-chunk-breadth`.
    User,
    /// `--low-biomass` without `--min-chunk-breadth`: no cutoff.
    LowBiomass,
    /// The automatic curve at the run's read count.
    Automatic,
}

/// The chunk-breadth cutoff of the species amount test: the curve's value
/// (or the given cutoff), times the user scale (`--chunk-breadth-scale`,
/// capped at 1), times the recall-mode scale.
#[derive(Clone, Copy, Debug, PartialEq)]
struct ChunkCutoff {
    source: ChunkSource,
    curve: f64,
    user_scale: f64,
    /// The curve's value times the user scale.
    user_scaled: f64,
    recall_scale: f64,
    /// The cutoff of the species amount test: `user_scaled` times the recall
    /// scale.
    cutoff: f64,
}

impl ChunkCutoff {
    /// The cutoff of a run whose automatic curve is read at `basis_reads`. A
    /// given `--min-chunk-breadth` is taken as it is, in every mode.
    fn resolve(cfg: &ProfilerConfig, basis_reads: usize) -> Self {
        let fixed = |source, cutoff| Self {
            source,
            curve: cutoff,
            user_scale: 1.0,
            user_scaled: cutoff,
            recall_scale: 1.0,
            cutoff,
        };
        if let Some(cutoff) = cfg.chunk_breadth {
            return fixed(ChunkSource::User, cutoff);
        }
        if cfg.lowbiomass {
            return fixed(ChunkSource::LowBiomass, 0.0);
        }
        let curve = chunk_breadth_curve(basis_reads, cfg.sequencer.is_illumina());
        let user_scaled = scale_automatic_chunk_breadth(curve, cfg.chunk_breadth_scale);
        let recall_scale = recall_chunk_scale(cfg, basis_reads);
        Self {
            source: ChunkSource::Automatic,
            curve,
            user_scale: cfg.chunk_breadth_scale,
            user_scaled,
            recall_scale,
            cutoff: user_scaled * recall_scale,
        }
    }

    fn is_user_scaled(&self) -> bool {
        self.user_scale != 1.0 && self.curve > 0.0
    }

    fn is_recall_scaled(&self) -> bool {
        self.recall_scale < 1.0 && self.user_scaled > 0.0
    }

    /// Where the cutoff comes from, in words, with the `basis` of the curve.
    fn source_in_words(&self, basis: &ChunkBasis) -> &'static str {
        match self.source {
            ChunkSource::User => "user-supplied --min-chunk-breadth",
            ChunkSource::LowBiomass => "disabled by --low-biomass",
            ChunkSource::Automatic => {
                let scaled = self.is_user_scaled() || self.is_recall_scaled();
                match (scaled, basis.scaled, basis.mapped) {
                    (true, true, _) => "recall-scaled auto-estimate from scaled mapped reads",
                    (true, false, true) => "recall-scaled auto-estimate from mapped reads",
                    (true, false, false) => "recall-scaled auto-estimate from total reads",
                    (false, true, _) => "auto-estimated from scaled mapped reads",
                    (false, false, true) => "auto-estimated from mapped reads",
                    (false, false, false) => "auto-estimated from total reads",
                }
            }
        }
    }

    /// Log the user and recall-mode scaling of the automatic cutoff.
    fn log_scaling(&self) {
        if self.is_user_scaled() {
            let capped = if self.curve * self.user_scale > 1.0 {
                ", capped at 1.000000"
            } else {
                ""
            };
            log::info!(
                target: "PROFILER",
                "User chunk-breadth scale adjusted the automatic cutoff: {:.6} x {} = {:.6}{capped}.",
                self.curve,
                self.user_scale,
                self.user_scaled,
            );
        }
        if self.is_recall_scaled() {
            log::info!(
                target: "PROFILER",
                "Recall mode scaled the automatic chunk-breadth cutoff: {:.6} x {:.4} = {:.6}; OEBR/COEBR thresholds are unchanged.",
                self.user_scaled,
                self.recall_scale,
                self.cutoff,
            );
        }
    }
}

/// Scale only an automatic cutoff in recall mode. The short-read
/// curve is continuous at one million reads and approaches 0.70 smoothly:
/// `scale = 0.70 + 0.30 / m`, where `m = basis_reads / 1e6`. Long reads use
/// a fixed 0.80 multiplier above the same depth gate.
fn recall_chunk_scale(cfg: &ProfilerConfig, basis_reads: usize) -> f64 {
    if cfg.mode != ProfileMode::Recall || basis_reads <= RECALL_CHUNK_SCALE_MIN_BASIS {
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

/// Scale the read count used as input to the automatic chunk-breadth curve.
fn scale_low_map_chunk_basis(chunk_basis_reads: usize, scale: bool) -> usize {
    if scale {
        ((chunk_basis_reads as f64) * LOW_MAP_CHUNK_BASIS_SCALE).round() as usize
    } else {
        chunk_basis_reads
    }
}

/// The genus of `taxid`: the taxon itself at genus rank, or its ancestor
/// there; `None` without one.
fn resolve_genus_taxid(taxonomy: &Taxonomy, taxid: u32) -> Option<u32> {
    taxonomy
        .get_parents(taxid)
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
fn chunk_breadth_curve(num_reads: usize, illumina: bool) -> f64 {
    let m = num_reads as f64 / 1_000_000.0;
    let value = if illumina {
        0.66 * m / (m + 1.0) + 0.24 * m / (m + 100.0)
    } else {
        0.85 * m / (m + 0.78)
    };
    value.min(0.95)
}

/// Aggregated data from mapping loading.
///
/// Contains all information needed for the profiling phase.
pub(super) struct MappingData {
    /// Per-read alignment bundles
    read_alignments: Vec<ReadBundle>,

    /// Chunk-breadth cutoff of the species amount test.
    min_chunk_breadth: f64,

    /// Whether at least one genome is indexed in part (its full size exceeds
    /// its indexed size).
    has_scaled_genomes: bool,

    /// Total reads in the input.
    total_reads: usize,

    /// The genus of every species with whole-genome references, also on a
    /// fractional database. Species represented only by marker genes are
    /// absent.
    genus_of_species: HashMap<u32, u32>,

    /// Most-specific reference-supported rank for every terminal taxid. Rows
    /// above species take their rank from the taxonomy.
    terminal_ranks: HashMap<u32, ProfileRank>,

    /// Marker genes and fragment geometry when the database has marker genomes.
    marker_model: Option<Arc<MarkerModel>>,

    /// A fractional database: reads not placed at a species are not reported
    /// above species, and no species is kept for even coverage.
    fractional: bool,

    /// Inputs of the levels above species and of the even-coverage species;
    /// `None` when they are off (`--strict-coverage`).
    above_species: Option<AboveSpeciesInputs>,

    /// Inputs of `--estimate-unclassified`; `None` when it is off.
    unclassified_species: Option<MarkerSpecies>,

    /// Genome accession of each assembly id.
    assembly_names: Vec<Arc<str>>,

    /// The strain taxid that each genome declares (`t:`), by assembly id.
    strain_taxids: Vec<Option<u32>>,
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

    #[test]
    fn expected_breadth_is_clamped_for_high_depths() {
        let expected = crate::profiler::evidence::expected_breadth(10.0, 5.0);
        approx_equal(expected, 1.0);
        assert!(expected.is_finite());

        let zero = crate::profiler::evidence::expected_breadth(10.0, 0.0);
        approx_equal(zero, 0.0);
    }

    /// The test configuration for reads of `sequencer`.
    fn cfg(sequencer: Sequencer) -> ProfilerConfig {
        ProfilerConfig {
            sequencer,
            ..ProfilerConfig::for_tests()
        }
    }

    fn test_taxonomy(directory: &Path) -> Taxonomy {
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
        Taxonomy::from_files(&nodes, &names, &merged, None).unwrap()
    }

    #[test]
    fn no_classifications_are_written_for_a_fractional_database() {
        let directory = tempfile::tempdir().unwrap();
        for fractional in [false, true] {
            let prefix = directory.path().join(format!("fractional-{fractional}"));
            let mut cfg = cfg(Sequencer::Illumina);
            cfg.save_classifications = true;
            cfg.outprefix = prefix.to_string_lossy().into_owned();
            cfg.db_metadata = Some(Arc::new(DatabaseMetadata {
                metam_index: std::path::PathBuf::from("metax_db.mxi"),
                references: Vec::new(),
                is_fractional: fractional,
            }));
            let profiler = Profiler::new(cfg, test_taxonomy(directory.path()));
            let mapping_data = MappingData {
                read_alignments: Vec::new(),
                min_chunk_breadth: 0.05,
                has_scaled_genomes: false,
                total_reads: 2,
                genus_of_species: HashMap::new(),
                terminal_ranks: HashMap::new(),
                marker_model: None,
                fractional: false,
                above_species: None,
                unclassified_species: None,
                assembly_names: Vec::new(),
                strain_taxids: Vec::new(),
            };
            profiler.profile(mapping_data).unwrap();
            let written =
                Path::new(&format!("{}.classify.txt.gz", prefix.to_string_lossy())).exists();
            assert_eq!(written, !fractional);
        }
    }

    #[test]
    fn reads_not_placed_at_a_species_go_above_species_only_on_a_full_database() {
        // Genus 100 (family 40) holds species 1001 and 1002, each with 125 reads spread along its 1 Mb genome. Both
        // fail the amount test, so their 250 reads, unique to genus 100, make a genus row on a full database and no
        // row on a fractional database.
        let directory = tempfile::tempdir().unwrap();
        let node = |parent: u32, rank: &str| (parent, rank.to_string());
        let taxonomy = || {
            Taxonomy::from_parts(
                HashMap::from([
                    (1, node(1, "no rank")),
                    (40, node(1, "family")),
                    (100, node(40, "genus")),
                    (1001, node(100, "species")),
                    (1002, node(100, "species")),
                ]),
                HashMap::new(),
                HashMap::new(),
                HashSet::new(),
            )
        };
        let alignment = |taxid: u32, assembly: u32, start: u32| EvidenceAlignment {
            profile_taxid: taxid,
            reference_id: assembly,
            assembly_id: assembly,
            role: MappingRole::Unpaired,
            pair_id: 0,
            sampling_scale: 1.0,
            asm_len: 1_000_000.0,
            identity: 0.99,
            fraction: 1.0,
            mapped_len: 100,
            read_len: 100,
            span: (start, start + 100),
            fixed_chunk: (start, start + 100),
            flex_chunk: (start.saturating_sub(950), start + 1_050),
            marker_gene_len: 0,
        };
        for fractional in [false, true] {
            let prefix = directory.path().join(format!("fractional-{fractional}"));
            let mut cfg = cfg(Sequencer::Illumina);
            cfg.above_species = true;
            cfg.outprefix = prefix.to_string_lossy().into_owned();
            let profiler = Profiler::new(cfg, taxonomy());
            let mapping_data = MappingData {
                read_alignments: (0..250u32)
                    .map(|i| ReadBundle {
                        read_name: Box::from(format!("read{i}")),
                        alignments: vec![alignment(1001 + i % 2, i % 2, (i / 2) * 7_900)],
                    })
                    .collect(),
                min_chunk_breadth: 0.9,
                has_scaled_genomes: false,
                total_reads: 250,
                genus_of_species: HashMap::from([(1001, 100), (1002, 100)]),
                terminal_ranks: HashMap::from([
                    (1001, ProfileRank::Species),
                    (1002, ProfileRank::Species),
                ]),
                marker_model: None,
                fractional,
                above_species: Some(AboveSpeciesInputs {
                    marker_only: HashMap::new(),
                    layout: GenomeLayout::new([(0, 1_000_000), (1, 1_000_000)]),
                }),
                unclassified_species: None,
                assembly_names: vec![Arc::from("GCF_1"), Arc::from("GCF_2")],
                strain_taxids: vec![None; 2],
            };
            profiler.profile(mapping_data).unwrap();
            let profile =
                std::fs::read_to_string(format!("{}.profile.txt", prefix.to_string_lossy()))
                    .unwrap();
            let ranks = profile
                .lines()
                .filter(|line| !line.starts_with('#'))
                .map(|line| line.split('\t').nth(2).unwrap())
                .collect::<Vec<_>>();
            let expected: &[&str] = if fractional { &[] } else { &["genus"] };
            assert_eq!(ranks, expected, "{profile}");
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
                read_name: Box::from("read"),
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
    ) -> EvidenceAlignment {
        EvidenceAlignment {
            profile_taxid: species,
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
        let mut cfg = cfg(Sequencer::Illumina);
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
            let mut cfg = cfg(Sequencer::Nanopore);
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
            let mut profiler = Profiler::new(cfg, taxonomy);
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

    #[test]
    fn the_chunk_basis_is_mapped_reads_on_low_map_runs_or_when_asked() {
        let (total, mapped) = (100_000usize, 10_000usize);
        let basis = |cfg: &ProfilerConfig, rate, fractional| {
            ChunkBasis::of(cfg, rate, (total, mapped), fractional)
        };
        let mut cfg = cfg(Sequencer::Illumina);
        // Total reads at a high mapping rate.
        assert_eq!(
            basis(&cfg, 0.90, false),
            ChunkBasis {
                mapped: false,
                scaled: false,
                reads: total
            }
        );
        // A low-map run on a full database: 1.5 times the mapped reads.
        assert_eq!(
            basis(&cfg, 0.10, false),
            ChunkBasis {
                mapped: true,
                scaled: true,
                reads: 15_000
            }
        );
        // On a fractional database the mapped reads as they are.
        assert_eq!(
            basis(&cfg, 0.10, true),
            ChunkBasis {
                mapped: true,
                scaled: false,
                reads: mapped
            }
        );
        // --by-mapped alone: the mapped reads as they are.
        cfg.by_aligned = true;
        assert_eq!(
            basis(&cfg, 0.90, false),
            ChunkBasis {
                mapped: true,
                scaled: false,
                reads: mapped
            }
        );
        // A given --min-chunk-breadth is never scaled.
        cfg.chunk_breadth = Some(0.25);
        assert!(!basis(&cfg, 0.10, false).scaled);
        // A host turns the low-map rule off.
        let mut cfg = cfg_with_host();
        cfg.chunk_breadth = None;
        assert!(!basis(&cfg, 0.10, false).mapped);
    }

    fn cfg_with_host() -> ProfilerConfig {
        ProfilerConfig {
            host: Some(crate::profiler::Host {
                taxid: 9606,
                pathogen_table: "pathogens.tsv".to_string(),
            }),
            ..cfg(Sequencer::Illumina)
        }
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
        // Half rounds away from zero: 3 * 1.5 = 4.5 => 5.
        assert_eq!(scale_low_map_chunk_basis(3, true), 5);
    }

    #[test]
    fn chunk_breadth_curve_matches_spec() {
        // Short (Illumina): 0.66*m/(m+1) + 0.24*m/(m+100), plateau 0.90.
        let short = |n| chunk_breadth_curve(n, true);
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
        let long = |n| chunk_breadth_curve(n, false);
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
            assert!(chunk_breadth_curve(b, true) >= chunk_breadth_curve(a, true));
            assert!(chunk_breadth_curve(b, false) >= chunk_breadth_curve(a, false));
        }
    }

    #[test]
    fn recall_scales_only_deep_automatic_chunk_cutoffs() {
        let mut cfg = cfg(Sequencer::Illumina);
        cfg.mode = ProfileMode::Recall;

        let at_gate = ChunkCutoff::resolve(&cfg, 1_000_000);
        approx_equal(at_gate.recall_scale, 1.0);
        assert_eq!(at_gate.cutoff.to_bits(), at_gate.user_scaled.to_bits());

        let two_million = ChunkCutoff::resolve(&cfg, 2_000_000);
        approx_equal(two_million.recall_scale, 0.85);
        approx_equal(
            two_million.cutoff,
            two_million.user_scaled * two_million.recall_scale,
        );

        let five_million = ChunkCutoff::resolve(&cfg, 5_000_000);
        approx_equal(five_million.recall_scale, 0.76);
        let deep = ChunkCutoff::resolve(&cfg, 16_600_000);
        approx_equal(deep.recall_scale, 0.70 + 0.30 / 16.6);
        assert!(deep.recall_scale > SHORT_RECALL_CHUNK_SCALE_FLOOR);

        cfg.sequencer = Sequencer::Nanopore;
        let long_at_gate = ChunkCutoff::resolve(&cfg, 1_000_000);
        approx_equal(long_at_gate.recall_scale, 1.0);
        let long_deep = ChunkCutoff::resolve(&cfg, 1_000_001);
        approx_equal(long_deep.recall_scale, LONG_RECALL_CHUNK_SCALE);
        approx_equal(
            long_deep.cutoff,
            long_deep.user_scaled * LONG_RECALL_CHUNK_SCALE,
        );

        cfg.sequencer = Sequencer::PacBio;
        approx_equal(
            ChunkCutoff::resolve(&cfg, 5_000_000).recall_scale,
            LONG_RECALL_CHUNK_SCALE,
        );
    }

    #[test]
    fn a_given_cutoff_and_low_biomass_are_not_scaled() {
        let mut cfg = cfg(Sequencer::Illumina);
        cfg.mode = ProfileMode::Recall;
        cfg.chunk_breadth = Some(0.37);
        let given = ChunkCutoff::resolve(&cfg, 20_000_000);
        assert_eq!(given.source, ChunkSource::User);
        assert_eq!(given.cutoff, 0.37);
        approx_equal(given.recall_scale, 1.0);

        cfg.chunk_breadth = None;
        cfg.lowbiomass = true;
        let low_biomass = ChunkCutoff::resolve(&cfg, 20_000_000);
        assert_eq!(low_biomass.source, ChunkSource::LowBiomass);
        assert_eq!(low_biomass.cutoff, 0.0);
        approx_equal(low_biomass.recall_scale, 1.0);

        // A read minimum leaves the automatic cutoff unchanged.
        cfg.lowbiomass = false;
        let without_min_reads = ChunkCutoff::resolve(&cfg, 20_000_000);
        assert_eq!(without_min_reads.source, ChunkSource::Automatic);
        cfg.min_reads = Some(3);
        assert_eq!(ChunkCutoff::resolve(&cfg, 20_000_000), without_min_reads);

        // Only recall mode scales.
        cfg.min_reads = None;
        for mode in [ProfileMode::Default, ProfileMode::Precision] {
            cfg.mode = mode;
            let cutoff = ChunkCutoff::resolve(&cfg, 20_000_000);
            approx_equal(cutoff.recall_scale, 1.0);
            assert_eq!(cutoff.cutoff.to_bits(), cutoff.curve.to_bits());
        }
    }

    #[test]
    fn the_user_scale_adjusts_the_automatic_cutoff_up_to_one() {
        let mut cfg = cfg(Sequencer::Illumina);
        let unchanged = ChunkCutoff::resolve(&cfg, 5_000_000);
        assert_eq!(unchanged.curve.to_bits(), unchanged.user_scaled.to_bits());
        assert_eq!(unchanged.curve.to_bits(), unchanged.cutoff.to_bits());

        cfg.chunk_breadth_scale = 0.5;
        let relaxed = ChunkCutoff::resolve(&cfg, 5_000_000);
        approx_equal(relaxed.user_scaled, relaxed.curve * 0.5);
        assert_eq!(relaxed.cutoff.to_bits(), relaxed.user_scaled.to_bits());

        cfg.chunk_breadth_scale = 10.0;
        let capped = ChunkCutoff::resolve(&cfg, 5_000_000);
        approx_equal(capped.user_scaled, 1.0);
        approx_equal(capped.cutoff, 1.0);
    }

    #[test]
    fn recall_scaling_composes_after_user_chunk_scaling() {
        let mut cfg = cfg(Sequencer::Illumina);
        cfg.chunk_breadth_scale = 0.5;
        cfg.mode = ProfileMode::Recall;

        let cutoff = ChunkCutoff::resolve(&cfg, 5_000_000);
        approx_equal(cutoff.user_scale, 0.5);
        approx_equal(cutoff.recall_scale, 0.76);
        approx_equal(cutoff.user_scaled, cutoff.curve * 0.5);
        approx_equal(cutoff.cutoff, cutoff.user_scaled * cutoff.recall_scale);
        assert_eq!(
            cutoff.source_in_words(&ChunkBasis {
                mapped: false,
                scaled: false,
                reads: 5_000_000
            }),
            "recall-scaled auto-estimate from total reads"
        );
    }
}

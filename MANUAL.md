# Metax Manual

This manual describes how to install Metax, prepare a reference database, run profiling, and interpret the output and the underlying method.

## 1. Installation

### 1.1 Install with Conda

```bash
conda create -n metax -c zldeng metax
conda activate metax
```

Packages are published for `linux-64`, `linux-aarch64`, and `osx-arm64`. The package contains a single self-contained executable.

### 1.2 Build from source

Install Rust 1.89 or newer and a C compiler, then build:

```bash
cargo build --release
```

The executable is written to `target/release/metax`. Add it to your `PATH` if desired:

```bash
export PATH="$PWD/target/release:$PATH"
```

Check the command-line interface:

```bash
metax --help
metax index --help
metax profile --help
metax profile --mapper-help
metax info --help
```

Run the unit tests and the end-to-end smoke tests:

```bash
cargo test
smoke_tests/run_smoke_tests.sh
```

## 2. Taxonomy Files

Metax needs NCBI-style taxonomy dump files.

### 2.1 Managed taxonomy (default)

When `--dmp-dir` is omitted, Metax uses `~/.metax/taxdmp`. On the first run, or when that copy is missing or invalid, Metax downloads the official NCBI `taxdump.tar.gz`, validates it, and keeps `nodes.dmp`, `names.dmp`, `merged.dmp`, and `delnodes.dmp`. Concurrent Metax processes share an installation lock, so a partially written dump is never used.

### 2.2 Explicit taxonomy directory

```bash
mkdir -p taxonomy
curl -L -o taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz -C taxonomy
metax profile --dmp-dir taxonomy ...
```

An explicit directory is never downloaded, repaired, or modified by Metax. Use it for offline installations, pinned taxonomy snapshots, or alternative taxonomies (e.g. GTDB or ICTV) converted to the same dump-file format. The directory must contain at least `nodes.dmp`, `names.dmp`, and `merged.dmp`.

## 3. Reference Databases

### 3.1 Pre-built database

A pre-built database for Metax 2 is available [here](https://s3.bi.denbi.de/swift/v1/cami3/metax/metax2_refseq_202607_k16w37.tar.gz) (38 GiB). It is based on the July 2026 NCBI RefSeq snapshot and covers bacteria, archaea, viruses, fungi, protozoa and the human reference genome (GRCh38), with up to three genomes per species for bacteria and archaea. Unnamed "sp." isolates are included only where GTDB (release 232) places them in a species of their own. In total, it contains 29,112 genomes.

```bash
tar -xzf metax2_refseq_202607_k16w37.tar.gz
metax profile -d metax2_refseq_202607_k16w37 \
  -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
  -o sample_pe -t 32
```

The archive holds `metax_db.mxi`, `metax_db.tsv.gz`, and `md5sums.txt` with the checksums of both. The first run writes an index cache of about 50 GiB into the unpacked folder, which must therefore be writable.

Databases built for Metax 0.9 use a different index format and must be rebuilt with `metax index`.

### 3.2 Reference FASTA header format

Every sequence header in the FASTA used for indexing must follow:

```text
>genome_id|sequence_id|taxonomy[|genome_size]
```

Fields:

- `genome_id`: identifier of the genome or assembly. All contigs of one genome share this value.
- `sequence_id`: identifier of the contig or sequence, unique within the genome.
- `taxonomy`: one or more `prefix:taxid` pairs separated by `;`, naming the lowest confidently known rank(s).
- `genome_size` (optional): total genome length in bases. When omitted it is computed from the summed contig lengths.

Supported taxonomy prefixes:

| prefix | rank |
| --- | --- |
| `t` | strain / subspecies |
| `s` | species |
| `g` | genus |
| `f` | family |
| `o` | order |
| `c` | class |
| `p` | phylum |
| `k` | kingdom |

Examples:

```text
>GCF_000001405.40|NC_000001.11|t:1234;s:9606|3200000000
>GCF_000005845.2|chromosome|s:562|4641652
>MAG_001|contig_42|g:561
```

Rules:

- If the strain or subspecies is known, include both `t:<strain_taxid>` and `s:<species_taxid>`.
- If only the species or a higher rank is known, include only that rank.
- The sequences of one genome do not need to be adjacent in the FASTA, but all of them must carry the same taxonomy and, when given, the same `genome_size`.
- Malformed headers, unsupported prefixes, duplicate contigs within a genome, and inconsistent taxonomies or genome sizes are rejected with an explanatory error.
- If a provided `genome_size` is larger than the summed contig lengths, indexing warns that only part of the genome is represented. If it is smaller, indexing fails.

Metax keeps the most specific supported `(rank, taxid)` of every genome as its terminal reporting label. A genome labelled `t:1234;s:9606` is reported as species 9606 in community profiling and as strain 1234 in strain profiling. A genome labelled `g:561` is reported directly at genus 561; it is never relabelled as a species.

### 3.3 Header rewriting workflow

Prepare a two-column mapping table with one row per source sequence:

```text
<original_header>	<target_header>
```

Rewrite the headers with `seqkit replace` or a similar tool, then validate the first few headers:

```bash
grep '^>' renamed_genomes.fa | head
```

## 4. Build an Index

### 4.1 Full reference index

```bash
metax index renamed_genomes.fa -o db -t 20
```

Gzip-compressed input is accepted. The output directory contains:

- `metax_db.tsv.gz`: compressed reference, genome, taxonomy, and size metadata;
- `metax_db.mxi`: the mapping index;
- `metax_reference.fasta.gz`: only for fractional indexes, the exact sampled FASTA used to build the index (not required for profiling).

`metax index` refuses to overwrite an existing output path. The completed database is published by a same-filesystem rename, so an interrupted build never leaves a partial database in place.

Pass the database directory to `metax profile --db`.

### 4.2 Fractional index

A fractional index samples segments from each large genome before index construction. It reduces index size, peak memory, and runtime:

```bash
metax index renamed_genomes.fa -o db_frac_10 -f 0.10 -s 42 -l 50000 -m 3000 -t 20
```

Options:

- `-f, --fraction`: fraction of each eligible genome to index (0–1).
- `-s, --subsample-seed`: random seed for reproducible sampling. Default `42`.
- `-l, --segment-length`: target length of the uniformly distributed segments. Default `50000`.
- `-m, --min-segment-length`: discard shorter sampled segments for sufficiently large genomes. Default `3000`.

Behavior:

- Genomes smaller than 500 kb are kept whole.
- Larger genomes are sampled reproducibly using genome-specific seeds, with segments distributed across contigs approximately uniformly.
- The full genome size and the effective indexed size are stored separately.

During profiling with a fractional index, coverage and depth use the effective indexed sequence, and assigned-read support is scaled per genome by `full_size / effective_size` before taxonomic aggregation. Fully indexed genomes keep a factor of 1. Per-read classification output is unavailable for fractional databases.

### 4.3 Marker-gene and mixed databases

A genome can be represented by its assembled sequences (reference type `genome`, the default) or by a set of marker genes with one FASTA record per gene (reference type `markers`), such as the single-copy genes that GTDB uses for its trees. Set the type of every genome, or of single genomes with a two-column table:

```bash
metax index marker_genes.fa -o db_markers --reference-type markers -t 20
metax index genomes_and_markers.fa -o db_mixed --reference-type-table types.tsv --dmp-dir taxonomy -t 20
```

`types.tsv` has one `genome_id<TAB>genome|markers` row per genome whose type differs from `--reference-type`; lines starting with `#` are ignored.

Marker genomes:

- are indexed whole, also with `--fraction`;
- have the summed marker length as genome size: read counts (profile column 4) are the reads mapped to the marker genes, not full-genome equivalents, and coverage refers to the marker genes. Depth (profile column 5) is that of the whole genome: a fragment on a marker gene adds the length of its read over the placements that the genome's genes offer (§7.13). Marker genes are largely single-copy, so relative abundances of marker and genome taxa are comparable even though their read counts are not;
- keep a `genome_size` from the FASTA header that is larger than their summed marker length as the size of the whole genome. It is used only for the unclassified estimate (§7.14) and never scales read counts; `metax info` reports how many marker genomes have one;
- are evaluated by the genes their reads hit rather than by chunk breadth (§7.4);
- never take part in genus fallback (§7.9).

A read counts on a marker gene only when enough of it lies inside the gene; a long read also counts on a gene shorter than that when it aligns to more than half of the gene (§7.4).

In a mixed database, reads from the part of a marker-only genome that is not indexed can map to the genomes of related species of type `genome`. With `--dmp-dir`, `metax index` lists the genera that hold species of both types; give the species of one genus the same type where possible.

Databases with marker genomes write metadata version 3, which adds a `reference_type` column and a `#reference_type genome|markers|mixed` line; earlier Metax releases cannot read it. Databases of genomes keep metadata version 2.

### 4.4 Index options

Seed scheme:

- `-k, --seed-length`: length of each syncmer component. Default `15`.
- `-w, --seed-window`: seed density/window parameter; larger values produce fewer seeds and a smaller index. Default `15`.
- `--max-occ`: discard seed keys occurring at more than N reference positions. Default `500` for an index of marker genes only (`--reference-type markers` without `--reference-type-table` or `--fraction`), where conserved genes share many seeds, and `100` otherwise.
- `--seed-mode`: `randstrobe` (default) or `minimizer`.
- `--hash-bits`: high hash bits stored as the lookup key (24–64). Default `32`.
- `--strobe-start-dist-min` / `--strobe-start-dist-max`: start-to-start distance range between the two strobes. Defaults `30` and `90`.
- `--first-strobe-only`: compatibility mode that matches only the first strobe; more sensitive but less specific than the default full-pair verification.

From index format version 7, a randstrobe seed keeps the orientation of both strobes. An index built by an earlier release keeps its own seeds, and reads are seeded to match them, so it maps as before; rebuild it to use the complete seeds.

Resources:

- `-t, --threads`: worker threads for sampling and parallel seed generation.
- `--seed-chunk-mib`: approximate reference bases held per parallel seeding batch. Default `1024`. A single record larger than this is kept intact.
- `--temp-dir`: existing directory for large temporary seed partitions and sort stores. Fast local scratch is recommended. Without it, temporaries are created inside the hidden staging area beside `--output`. They are removed after indexing.

Progress:

- `-q, --quiet` suppresses progress and summary messages; `--progress` forces live progress output; `--no-progress` disables it; `--progress-interval N` sets the record interval for count-based updates.

### 4.5 Inspect a database

```bash
metax info -d db
metax info -i db/metax_db.mxi
```

`info -d` summarizes the database metadata, including its reference type and number of marker genomes, and, when the `.mxi` file is present, the mapping-index header. `info -i` inspects an index directly. Neither command loads the full index into memory.

## 5. Run Profiling

### 5.1 Illumina paired-end reads

```bash
metax profile -d db --dmp-dir taxonomy \
  -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
  -o sample_pe -t 32
```

Supplying `-I` enables paired-end mode. By default only concordant same-reference pairs are emitted by the mapper; see §5.6 for the emission policy.

### 5.2 Illumina single-end reads

```bash
metax profile -d db --dmp-dir taxonomy -i sample.fastq.gz -o sample_se -t 32
```

### 5.3 Nanopore reads

```bash
metax profile -d db --dmp-dir taxonomy -i ont_reads.fastq.gz -s Nanopore -o sample_ont -t 32
```

### 5.4 PacBio reads

```bash
metax profile -d db --dmp-dir taxonomy -i pacbio_reads.fastq.gz -s PacBio -o sample_pacbio -t 32
```

### 5.5 Profiling modes and read-level thresholds

`--mode` selects an evidence preset: `default` (balanced), `recall` (more permissive), or `precision` (stricter). The profiler's read-level thresholds are:

| Input and profile | Recall | Default | Precision | Gate |
| --- | --- | --- | --- | --- |
| Illumina community | identity 0.96; 50 bases; fraction 0.60 | identity 0.97; 50 bases; fraction 0.70 | identity 0.98; 50 bases; fraction 0.80 | length AND fraction |
| Illumina strain | identity 0.96; 50 bases; fraction 0.65 | identity 0.97; 50 bases; fraction 0.80 | identity 0.98; 60 bases; fraction 0.80 | length AND fraction |
| Nanopore/PacBio community | identity 0.90; 300 bases; fraction 0.20 | identity 0.92; 500 bases; fraction 0.20 | identity 0.94; 500 bases; fraction 0.30 | length OR fraction |
| Nanopore/PacBio strain | identity 0.90; 500 bases; fraction 0.30 | identity 0.92; 500 bases; fraction 0.30 | identity 0.94; 500 bases; fraction 0.30 | length OR fraction |
| Assembly | identity 0.92; 400 bases; fraction 0.55 | identity 0.93; 500 bases; fraction 0.60 | identity 0.95; 600 bases; fraction 0.70 | length AND fraction |

Override any of them explicitly:

```bash
metax profile ... --min-identity 0.95 --mapped-len 100 --min-mapped-fraction 0.65
```

These options decide whether a mapping is strong enough for species/strain assignment. They do not discard lower-identity records emitted by the mapper, which remain available to genus fallback (§7.9).

**Noisy long reads.** For Nanopore/PacBio reads, Metax samples the first 10,000 reads (up to 250 million quality bases) and computes the mean base accuracy from the Phred scores. Reads below `0.95` count as noisy, which changes two defaults:

- **Identity.** For community profiling without `--min-identity`, the identity cutoff becomes `accuracy − margin`, where the margin is `0.03` (default), `0.04` (recall), or `0.02` (precision), never below `0.84`. The identity in the table is the high-quality default; strain mode and an explicit `--min-identity` keep the fixed preset.
- **Seed lookup.** The mapper also looks a seed up by its first strobe alone, as it does for Illumina reads, because both strobes of a seed are rarely free of errors in such reads (§5.6). Many more noisy reads are mapped, at some cost in mapping time.

High-quality samples, FASTA input and short reads keep both defaults; on a database with marker genes, long reads always use this seed lookup (§5.6).

**Qualities that over-state accuracy.** Some basecallers report higher qualities than their reads achieve. After loading the mappings, Metax therefore takes each read's best alignment (largest identity × mapped bases among those passing the length / fraction gate) and the upper quartile of their identities: reads from genomes the database holds well align at their sequencing accuracy, while divergent or poorly represented organisms lower only the lower quantiles. If this quartile is below the sampled accuracy, it takes the accuracy's place in the identity cutoff above (also for samples whose qualities claim `0.95` or more). The cutoff is only ever lowered. This applies to the automatic long-read community cutoff, including one saved in a mapping file by an automatic run; strain mode and an explicit `--min-identity` are unchanged.

### 5.6 Mapper options

Options before `--` configure the profiler; options after `--` configure the mapper, i.e. which candidate mappings reach the profiler. Show the complete reference with:

```bash
metax profile --mapper-help
```

Resolved mapper defaults:

| Input | top-n | score ratio | min seeds | mapper identity | mapper fraction |
| --- | ---: | ---: | ---: | ---: | ---: |
| Illumina single-end / assembly | 5 | 0.90 (0.80 in recall mode; 0.50 on a database with marker genes) | 2 | 0.80 | 0.50 |
| Illumina paired-end | 50 | 0.50 | 1 | 0.80 | 0.50 |
| Nanopore / PacBio | 10 | 0.70 | 3 | 0.80 | 0.30, or 300 mapped bases |

Near-identical marker genes of sister species share most seeds of a short read, and seed sampling can give the gene a read came from fewer seed hits than one of them, so on a database with marker genes an unpaired short read keeps candidates down to half of the best seed score, as a pair does, and alignment identity decides. Illumina reads and noisy long reads (§5.5) find a seed by its first strobe as well; accurate long reads require both strobes of the seed. A short gene holds few seeds, so on a database with marker genes long reads also find a seed by its first strobe, and on a marker gene a chain needs at most two seeds instead of `min seeds`. `--allow-single-strobe-matches` and `--match-both-strobes` override the seed lookup.

The paired-end candidate set is deliberately broad so that the profiler can see competing species and keep conserved regions ambiguous; the modes differ in the evidence stage, not in candidate generation. Example overrides:

```bash
metax profile ... -- --top-n 20 --score-ratio 0.80
```

Faster paired-end mapping: `--top-n 10` retains at most ten placements per mate at the default score ratio. It is considerably faster and leaves species and genus profiles practically unchanged where the database covers the community well. Where the database covers a community poorly, as is typical for soil, it loses recall above species level, so keep the default there.

```bash
metax profile ... -- --top-n 10
```

Note that `--min-identity` appears on both sides of `--` with different meanings: before it, the profiler's species/strain threshold; after it, the mapper's emission cutoff.

Paired-end emission: `--paired-emission concordant-only` (default) emits same-reference pairs and, for a fragment without any pair, its mates whose mate may lie beyond the end of their reference (§7.1); `concordant-first` emits concordant pairs when present and otherwise discordant rows; `all` emits every mate mapping. The latter two are mainly useful together with `--save-mapping` for diagnostics. A pair is written once, also when mate rescue finds again a placement that a mate already has.

Mate rescue: a mate that seeding did not place is searched for in the insert window of its partner's placements. `--mate-rescue-taxa N` (default `8`) extends this search to fragments that already have a proper pair: for up to `N` other taxa per mate placed on that mate only, the missing mate is searched for and kept when that pair matches at least as many bases as the best seeded pair (§7.1). `--mate-rescue-taxa 0` restricts mate rescue to fragments without any proper pair, the cases that Metax 2.0.2 rescued.

### 5.7 Saved and reused data

Normal profiling writes the coordinate mappings to a temporary compressed file, releases the mapping index, profiles the records, and removes the file. The following options change that behavior.

**Keep the mappings.**

```bash
metax profile ... --save-mapping
```

Writes `<prefix>.mapping.tsv.gz` (format in §6.4).

**Re-profile saved mappings** with different profiler settings, without the reads and without remapping. The same database must be supplied:

```bash
metax profile --reuse-mapping sample.mapping.tsv.gz -d db --dmp-dir taxonomy -o sample_reprofile -t 20
```

The file header supplies the platform, pairing state, mapper mode, and the original read-level gates, so `--sequencer`, `-I`, and mapper arguments must not be repeated. Profiler options such as `--mode`, `--strain`, evidence thresholds, coverage filters, genus fallback, and output controls remain available and override the recorded values. Changing `--mode` on reuse changes profiler thresholds only; the saved candidate set cannot be expanded.

**Resume** reuses `<prefix>.mapping.tsv.gz` if it exists and otherwise maps from the supplied reads:

```bash
metax profile ... --resume
```

SAM/BAM input from other mappers is no longer supported; `--reuse-sam` stops with an error. Save the mappings of a run with `--save-mapping` and profile them again with `--reuse-mapping`.

**Stream without temporary disk space.** `--stream-mappings` passes mappings directly to the profiler. It can increase peak memory and cannot be combined with `--save-mapping`.

**Diagnostic outputs.** `--save-raw-profile` writes `<prefix>.rprofile.txt`; `--save-classification` writes `<prefix>.classify.txt` (§6).

### 5.8 Pathogen-focused analysis

```bash
metax profile ... --host 9606 --pathogen-table pathogen_hosts.tsv
```

With `--host`, only pathogens associated with the given host taxid are reported, and three annotation columns are appended to the profile. The table is tab-separated with a header row and four columns:

```text
pathogen_taxid	host_taxids	host_names	diseases
```

Multiple host taxids, host names, or diseases are semicolon-separated. Rows with fewer than four columns are skipped. For viral pathogens, such a table can be compiled from the [Virus-Host Database](https://www.genome.jp/ftp/db/virushostdb/).

### 5.9 Strain profiling (experimental)

```bash
metax profile ... --strain
```

Strain-aware profiling reports at the strain/subspecies label when the database carries `t:` taxids. It runs an independent filter chain with its own identity, chunk-breadth, and ratio presets, does not use the community evidence policy described in §7, and uses a coverage-probability cutoff of `1e-10` instead of the community `1e-5`. Treat strain-level calls as provisional and validate them on your own data.

### 5.10 Resources, robustness, and logging

- `-t, --threads N`: worker threads shared by mapping and profiling. Default `20`.
- `--temp-dir DIR`: directory for temporary staged mappings. Without it, Metax prefers suitable local temporary storage and falls back to network temporary storage or the output directory.
- `--batch-size RECORDS`: profiler work-batch size for file-backed input.
- `--max-em-iterations N`: cap on ambiguous-read EM iterations. Default `200`; converged runs stop earlier.
- `--lenient-mapping-input`: skip malformed rows in reused mapping files with warnings instead of stopping.
- `-v, --verbose` and `--very-verbose`: show resolved settings and progress; the latter includes the complete mapper configuration.
- `--progress`, `--no-progress`, `--progress-interval N`: control live progress displays and periodic mapper progress messages.

## 6. Output Files

For an output prefix `sample`, Metax writes:

- `sample.profile.txt`: final taxonomic profile;
- `sample.log`: run log;
- `sample.rprofile.txt`: post-EM profile before taxon-level filters, with `--save-raw-profile`;
- `sample.classify.txt`: per-read assignments, with `--save-classification`;
- `sample.mapping.tsv.gz`: coordinate mappings, with `--save-mapping`.

With `--host`, the profile files use the prefix `sample.pathogen`.

### 6.1 Profile columns

The profile is tab-separated, has no header row, and contains one row per reported taxon. Numeric fields are written with five decimals.

| Column | Meaning |
| ---: | --- |
| 1 | taxon name |
| 2 | taxid |
| 3 | rank |
| 4 | number of reads, in full-index equivalents (scaled for fractional databases; for marker genomes, the reads mapped to the marker genes) |
| 5 | depth of coverage (for marker genomes, of the whole genome; §7.13) |
| 6 | relative abundance (%), normalized over the reported taxa |
| 7 | observed breadth of coverage |
| 8 | expected breadth of coverage |
| 9 | p-value of the breadth deviation (`NA` when not computed) |
| 10 | fixed chunk breadth |
| 11 | flexible chunk breadth |
| 12 | expected flexible chunk breadth |
| 13 | p-value of the flexible chunk breadth deviation (`NA` when not computed) |

With `--host`, three columns are appended:

| Column | Meaning |
| ---: | --- |
| 14 | host names |
| 15 | host taxids |
| 16 | diseases |

With `--estimate-unclassified`, the profile ends with one more row: name and rank `unclassified`, taxid `0`, column 4 the input fragments (single-end and long reads: reads) that no row explains, column 6 their share of the input in percent, and `NA` in the other columns (§7.14). Column 6 of this row is a share of the input, not part of the normalized abundance of the other rows.

For species represented by marker genes, columns 10 and 11 hold the fraction of marker genes hit, and column 12 the fraction expected to be hit by the same reads; column 13 is the p-value of the genes hit (§7.4).

Terminal rows are normalized once across the whole profile, not per rank. A genus-only genome and species-resolved genomes from the same genus are reported as separate rows; Metax does not add the species rows back into the genus row.

### 6.2 Raw profile

`rprofile.txt` has the same columns as the profile. It starts from all mappings that pass the read-level identity, mapped-length, and mapped-fraction gates, applies taxonomy normalization, majority-LCA handling, and EM refinement, and is written before the minimum-read, breadth, chunk-breadth, ratio, p-value, and coordinate-evidence filters. A species rejected from the final profile therefore appears here together with the statistics that caused its rejection. Abundance is normalized across all raw rows; renormalize if you filter the table yourself.

### 6.3 Per-read classification

`classify.txt` lists each classified read with its assigned taxid, the rank, the candidate taxids, and their EM fractions. A row is written only when at least one of the read's assigned taxa survives into the final profile; remaining fractions are renormalized when filtering removes candidates. This output is unavailable for fractional databases.

### 6.4 Mapping file

`<prefix>.mapping.tsv.gz` is tab-separated and contains two comment-prefixed run headers, one row per retained mapping, and one final sample-count line:

```text
#metax_mapping  version=5  sequencer=Illumina  paired=true  mapper_mode=default  strain=false  profiler_identity=0.97  profiler_mapped_len=50  profiler_fraction=0.7  sampled_mean_base_accuracy=NA
#metax_mapping_database  fractional=false  sequences=71442  genomes=33533  fingerprint=0123456789abcdef0123456789abcdef
S0R0  GCF_000001405.40|NC_000001.11|s:9606|3298430636  1872848  1872998  0.99333  150  150  =  1  1  0
#metax_mapping_sample  total_reads=33237808  mapped_reads=23098244  emitted_records=85758065
```

Mapping columns: read id, reference id, start, end, identity, read length, mapped bases, mate reference status (`=` mate on the same reference, `>` mate may lie beyond the end of this reference, `*` mate unmapped, `NA` single read, or the reference of the mate), mate number (`0`, `1`, or `2`), read-local pair id, and a candidate-boundary mask (`0`–`3`) indicating whether additional endpoints may have been omitted by the candidate cap. Identity is canonicalized to five decimal places. Reference, genome, taxonomy, and size tables are never copied into the mapping file; they come from the database, whose fingerprint is validated on reuse. Headerless or obsolete mapping files are rejected.

## 7. Method

### 7.1 Mapping and candidate generation

The built-in mapper indexes references with syncmer-based randstrobe seeds and retains up to `top-n` candidate mappings per read or mate whose seed score is at least `score-ratio` of the best. Candidates are refined to reference intervals with a per-base identity: an ungapped comparison for short reads, and for long reads a gapped identity computed by exact edit distance on intervals longer than 512 bases. Mapper-level acceptance requires identity of at least `0.80` and the platform's mapped-fraction gate (§5.6).

A long read on a marker gene is also accepted when it aligns to more than half of the gene. A long read can run over several genes of a genome, so when any of its candidates lies on a marker gene, candidates on different parts of the read (read intervals overlapping by less than half of the shorter one, directly or through other candidates) do not compete: the score ratio and `top-n` apply within each part. A candidate on a genome reference joins the parts it overlaps, so a read placed along a genome stays one part. Both strobes of a seed must lie inside a gene, so the seeds of a long read stop short of the gene ends; its alignment on a marker gene is therefore extended into the gene with a banded x-drop alignment for as long as the read keeps matching, and its identity, mapped bases and reference interval are those of the extended alignment.

Before applying `top-n`, the mapper reserves the strongest accepted mapping for each represented profile endpoint (normally species) and each most-specific endpoint (strain when known), then fills the remaining slots with the next best mappings. Repeated placements from one taxon therefore cannot hide another eligible taxon. A candidate-boundary mask records whether endpoints may have been omitted; the profiler keeps such candidates in abundance EM but does not force-resolve an incomplete candidate set. A long read whose records all lie on marker genes is the exception: its omitted candidates are copies of the same genes in further relatives, and in a genus with many species the cut falls on the most conserved genes, so the read is resolved as if its candidate set were complete; otherwise those genes would have no coordinate evidence and the breadth checks would read the gap as missing coverage.

Paired-end mapping is pair-aware: mates are scored together within the insert-size window, and only concordant same-reference pairs are emitted by default. A mate without a seeded placement, or a fragment without any proper pair, triggers mate rescue: the read is searched for in the insert window around its partner's placements. A placement found there is evaluated and accepted like a seeded one, so its identity and mapped bases describe the aligned interval. Only on a reference of at most 256 bases, when no such placement is found, is the read aligned as a whole at every position of the window. (Up to Metax 2.0.2 a rescued mate was aligned as a whole read: bases outside the interval lowered its identity, and it was accepted only when at least the identity cutoff of the whole read matched.)

A mate that remains without a mate on its reference after rescue is kept as a lone mate when the insert window on the side where its mate must lie (downstream of a forward read, upstream of a reverse one) runs past the end of the reference: its mate may come from sequence that the reference does not contain. This applies to marker genes and to every reference of a fractional index, whose ends cut the sequence a fragment came from; in a mixed database each reference follows its own type. On genome references of a full index no mate is kept this way, since the end of an assembled sequence mostly does not cut a fragment. Lone mates are emitted only for fragments without any pair and are marked `>` in mapping files. Other lone mates are dropped.

A proper pair on one taxon does not end the search. Seeding can place a second taxon on one mate only, for example a close relative for which the other mate has too few seed hits, and concordant-only emission would then drop that taxon although it may explain the fragment equally well or better. For each such taxon (profile endpoint, up to `--mate-rescue-taxa` per mate) the missing mate is searched for beside the taxon's best placement. The recovered pair is kept only if it matches at least as many bases as the best pair that seeding placed on one reference, so rescue can add an equally good or better explanation of the fragment but never a worse one. Kept mates are appended after the `top-n` cut and therefore never take a `top-n` slot from a seeded placement.

### 7.2 Read-level gates

Each mapping is tested against the profiler's identity, mapped-length, and mapped-fraction thresholds (§5.5). Short reads and assemblies require both the length and the fraction gate; long reads pass with either. Mappings that fail remain available only to genus fallback.

### 7.3 Coordinate evidence and read assignment

Native single-end Illumina, paired-end Illumina, Nanopore, and PacBio community profiling keep presence evidence separate from abundance ambiguity:

- every qualifying terminal-taxon candidate remains available to fractional EM abundance estimation;
- only a clearly resolved best terminal taxon may create coordinate coverage;
- coordinates are evaluated in one assembly coordinate system at a time — contigs of the same assembly may be combined, but coordinates from different assemblies are never unioned;
- the strongest coherent assembly witness determines presence.

This prevents a read from creating breadth for several taxa merely because it maps to a region they share.

**Single-end Illumina.** Mappings are ranked by `estimated_errors = (1 − identity) × mapped_bases`. One best placement is retained per `(terminal taxid, assembly)`. A read contributes coordinates only when it has a single candidate endpoint or the runner-up is worse by at least `0.5` estimated errors.

**Paired-end Illumina.** Mate 1 must pass the identity threshold, both mates must pass the length and fraction gates, and both must belong to a concordant same-reference pair. Pair errors are the sum over both mates. A fragment contributes coordinates only when it has a single candidate endpoint or `runner_up_pair_errors − best_pair_errors ≥ 1.0`. Mate 1 supplies the mapped bases used for depth and the intervals used for coverage; mate 2 contributes to pair discrimination. A fragment without any pair is placed by its lone mates (§7.1), each of which must pass all read-level gates and is treated like a single-end read. Such a fragment is placed in the same way by a mate of a pair on a marker gene whose partner fails its gates and reaches the first or last base of the gene: that partner runs off the gene and is clipped there, as the mate of a lone mate would be.

**Nanopore/PacBio.** Candidates are ranked by `long_score = mapped_bases × (2 × identity − 1)`. Coordinate evidence requires the best terminal score to exceed the runner-up by at least 1 % of the best score. Exact ties keep the mapper's first-ranked placement.

Exact and near-exact cross-species ties stay in the EM candidate set and create no coordinate evidence, unless resolved by guarded MAP (§7.8). Reads that map competitively across terminal ranks may be reported at their lowest common ancestor when that rank has sufficient retained evidence; otherwise they remain unclassified. Assembly input and strain profiling use the legacy evidence rules instead of this policy; strain profiling uses pairs only.

### 7.4 Coverage metrics

For fixed-chunk coverage, every accepted read span is extended by 1,000 bases on both sides and clamped to its contig. For a witness with assembly coordinate size `G`, `n` resolved fragments, and accepted spans `span_i`:

```text
observed_breadth   = merged_span_bases / G
depth              = sum(span_i) / G
mean_span_fraction = clamp(depth / n, 0, 1)
expected_breadth   = 1 − (1 − mean_span_fraction)^n
OEBR               = observed_breadth / expected_breadth
```

Fixed and flexible chunks are merged per contig in the same way; the flexible chunk count scales with the square root of the genome length. `COEBR` is the observed/expected ratio for flexible chunk breadth.

For a fractional index, `G` is the effective indexed size, and support and abundance are scaled per genome as described in §4.2.

**Marker genomes.** Each reference sequence of a marker genome is one gene, and a fragment counts on a gene only when enough of a read lies inside it. Chunks of equal length would expect more coverage than reads can give, most of all for short genes, so the genes themselves serve as chunks. From the sample's mappings Metax takes the median read length `R`, the median template length `T` of the pairs, and the bases `a` that a read needs inside a gene to pass the mapped-length and mapped-fraction gates (both gates for short reads, either of them for long reads). A long read also counts on a gene when it aligns to more than half of the gene, so it needs `m = min(a, ⌊L/2⌋ + 1)` bases inside a gene of length `L`. A gene of length `L` then offers

```text
reads:                 c = L − 2a + R + 1
long reads:            c = L − 2m + R + 1
pairs:                 c = L − T + 1
pairs and lone mates:  c = min(L − 2a + T + 1, 2 × (L − 2a + R + 1))
```

placements (none when `L < a`; every gene offers placements to long reads). A long read can run over several genes of the genome, and the two mates of a pair can lie on two neighbouring genes; such a fragment counts on each gene it covers. With `n` gene hits (fragments, counting every gene a fragment covers) on the genome and `p_i = c_i / Σc`:

```text
genes_hit                 = genes with at least one placed fragment
expected_genes_hit        = Σ_i 1 − (1 − p_i)^n
fixed and flex chunk      = genes_hit / genes
expected flex chunk       = expected_genes_hit / genes
COEBR                     = genes_hit / expected_genes_hit
expected_breadth          = Σ over marker bases of 1 − (1 − p_i × k_b / s_i)^n, divided by the marker length
```

where `k_b` is the number of the `s_i` possible read starts on gene `i` whose read covers base `b`, so bases near gene ends are expected to be covered less often. The p-value of the genes hit is the exact probability of hitting at most `genes_hit` genes (multiplied by the number of tests), and the breadth p-value is computed as for genomes. A marker species must hit at least `--min-marker-genes` genes (default `3`) in every evidence tier.

### 7.5 Automatic chunk-breadth cutoff

`-c, --min-chunk-breadth` sets the fixed chunk-breadth requirement directly. When it is omitted, Metax estimates it from the read count. With `m = reads / 1,000,000`:

| Platform | cutoff |
| --- | --- |
| Illumina | `0.66·m/(m+1) + 0.24·m/(m+100)` |
| Nanopore / PacBio | `0.85·m/(m+0.78)` |

Both curves are capped at `0.95`. The read-count basis is the total input by default; with `--by-mapped`, or automatically on low-map runs, it is the mapped reads/fragments. A run is a low-map run when its genome-equivalent mapping rate is below `--low-map-rate-threshold` (default `0.30`): each mapped fragment counts as the sampling scale of the genome of its best placement (`1` for a genome indexed whole, `full_size / effective_size` for a sampled one), fragments placed only on marker genes are left out, and the sum is divided by the input fragments and capped at `1`. For a database of whole genomes this is the plain mapping rate; for a fractional or marker database it keeps samples that the reference explains well out of low-map mode.

Adjustments:

- `--chunk-breadth-scale FACTOR` multiplies the automatic cutoff (capped at `1.0`). It conflicts with `--min-chunk-breadth` and `--low-biomass`.
- In community **recall** mode, once the basis exceeds one million reads the standard cutoff is additionally scaled by `0.70 + 0.30/m` for Illumina and by `0.80` for long reads. This scaling affects only the standard presence tier; the safeguard used by rescue tiers and guarded MAP is unaffected.
- `-n, --min-reads N` adds a minimum read count to the automatic cutoff. `-l, --low-biomass` disables the cutoff, and `--min-chunk-breadth 0` turns it off.
- Strain profiling uses its own piecewise cutoff and is not affected by recall scaling.

### 7.6 Ratio filters

```text
oebr  = observed_breadth / expected_breadth
coebr = observed_flex_chunk / expected_flex_chunk
```

Default bounds:

| Rank | OEBR lower | COEBR lower | Upper |
| --- | ---: | ---: | ---: |
| species (short reads) | 0.70 | 0.75 | 1.50 |
| species (long reads) | 0.70 | 0.70 | 1.50 |
| genus and higher | 0.65 | 0.65 | 1.50 |

`--min-oebr` and `--min-coebr` override the lower bounds; the upper bound stays `1.50`. `-b, --min-breadth` additionally requires a minimum observed breadth. Values given explicitly for `--min-breadth`, `--min-chunk-breadth`, `--min-oebr`, and `--min-coebr` apply to every reported taxon, including taxa admitted by a rescue tier (§7.8).

### 7.7 P-values

Species rows carry two p-values: one for the observed breadth against its expectation (column 9) and one for the flexible chunk breadth (column 13). Rows above species generally report `NA`.

- If the observed/expected ratio lies within `[0.9, 1.2]`, the p-value is `1.0`: the observation agrees with the model.
- Otherwise the deviation is tested one-sided, towards whichever side it fell, using an exact occupancy model or its normal/binomial approximation depending on sample size.
- P-values are multiplied by the number of taxa tested and capped at `1.0`.

Community filtering removes a row only when both available p-values are below `1e-5`. With `--strict`, both must be at least `1e-5`.

### 7.8 Rescue tiers and guarded MAP

Standard species evidence requires the chunk-breadth amount, OEBR and COEBR within their bounds, and the p-value policy. Species that miss the standard amount may still be admitted through rescue tiers, where `cutoff` is the automatic (or explicit) chunk-breadth cutoff after `--chunk-breadth-scale`:

- **Single-end:** fixed chunk breadth above `0.40 × cutoff`, resolved-best support for at least 85 % of resolved comparisons and 20 % of all candidate reads, non-negative OEBR, and COEBR in `[0.75, 1.50]`.
- **Paired-end:** fixed chunk breadth above `0.50 × cutoff`, resolved-best support for at least 75 % of resolved comparisons and 20 % of all candidate fragments, non-negative OEBR, and COEBR in `[0.75, 1.50]`.
- **Long reads:** fixed chunk breadth above `0.50 × cutoff`, at least 50 % best support among all candidate reads, and OEBR and COEBR each at least `0.70`.
- **Coverage specificity:** admitted below the half cutoff when at least 90 % of candidate fragments are a unique or margin-resolved winner for the taxon, with depth at least `0.01`, at least three fragments, in-band OEBR/COEBR, and the standard p-value policy.
- **Low-abundance paired-end:** fixed chunk breadth above `cutoff × min(0.50, 0.20 + 0.45 × cutoff)`, at least 20 fragments, depth at least `0.001`, at least 10 direct clear winners, direct-clear support of at least 25 % of all candidates, resolved-best support of at least 75 %, best support of at least 20 %, OEBR in `[0.95, 1.50]`, and COEBR in `[0.60, 1.50]`.

Rescue tiers apply to species only; direct higher-rank endpoints use the standard rule. Explicit `--strict`, `--min-oebr`, and `--min-coebr` remain authoritative for every tier, and explicit breadth and chunk-breadth cutoffs are applied to every reported taxon afterwards.

**Guarded MAP.** Native Illumina community profiling resolves exact top-score cross-species ties by default (`--exact-ties guarded-map`). It seeds an EM prior from species with independent unique-candidate support, requires a single largest posterior of at least `0.95`, preserves every species selected by the margin rule, and requires independent pre-MAP fixed-chunk breadth on the same assembly: above `0.50 × cutoff` for paired-end (or `0.05 × cutoff` for a species meeting the low-abundance tier) and above the full cutoff for single-end. MAP affects coordinate evidence only; abundance still uses fractional EM over all retained candidates. Use `--exact-ties ambiguous` to keep every exact tie ambiguous. Guarded MAP is not applied to long reads, strain profiling, or assembly input, and it falls back to `ambiguous` when `--low-biomass` or `--min-chunk-breadth 0` removes the chunk-breadth safeguard.

### 7.9 Genus fallback

Genus fallback assigns a read to a genus when none of its mappings qualifies for species/strain evidence. It is enabled explicitly with `-g, --genus-fallback`, or automatically for Illumina community runs without `--host` on a full (non-fractional) database when the genome-equivalent mapping rate (§7.5) is below `--low-map-rate-threshold` (default `0.30`). Species represented only by marker genes never take part: reads that fail every species on a marker gene come mostly from the homologous genes of relatives.

A fallback mapping must have identity of at least `--genus-identity` (default `0.80`), a mapped length of at least 50 bases for Illumina or the configured species `--mapped-len` otherwise, and pass the platform-aware genus coverage gate. The genus taxid is resolved from the species taxid through the taxonomy; species/strain evidence always takes precedence.

### 7.10 Reads of removed species at genus

The coverage evidence removes a species whose reads cover its reference too little or unevenly. These reads often still belong to the species' genus: the sample organism may be a divergent strain that shares only part of its genome with the reference, or a new species of a genus the database holds. In community profiling, a fragment none of whose candidate species was kept is credited to the candidate of its best placement, the donor, and is reported at the genus when all its candidates lie in that genus. A donor takes part only when its own evidence spreads over its genome:

- chunk ratio (COEBR, §7.6) of at least `0.40`, and
- window spread of at least `0.50`: the genome of the donor's strongest witness is cut into `min(100, n/5)` equal windows for its `n` resolved fragments, and the spread is the share of windows that hold the start of at least one of them. With fewer than 4 windows the donor does not take part.

Before that, a removed species whose reads cover its whole genome evenly is reported at species instead, when no taxon of its genus is kept by the coverage evidence: at least 200 fragments, OEBR and COEBR between `0.9` and `1.1`, a window spread of at least `0.9`, breadth and chunk p-values of at least `1e-5`, and fixed chunk breadth below the cutoff (§7.5). Such a species, for example a large genome at low depth, has too few reads for the amount tests, but its reads land across its genome as an organism's own reads do. It enters abundance estimation like a kept species, and its fragments do not go to the genus. Species represented by marker genes do not take part, and on a fractional database no species is reported this way. The log reports their number; `--verbose` lists them.

A genus row needs at least 200 fragments from its donors. It enters abundance estimation with the other taxa, and its coverage columns show the donor with the most fragments, as a majority-LCA row shows its witness. `--min-reads` and explicitly set coverage cutoffs apply to it. Species represented by marker genes do not take part. On a fractional database a genus that holds a species kept by the coverage evidence gets no row: reads from the unindexed parts of that species' genome land on the sampled segments of relatives and leave them removed, and the species' abundance already includes them through the scaling to its full genome (§4.2). Host genomes in the database take part like any other: reads of a host whose genome they cover too thinly for a species call are reported at species when they cover it evenly, and otherwise at its genus, with an abundance that follows their small depth. The log reports the number of genus rows; `--verbose` lists the donors of each.

Foreign contigs in a reference, a misassembled region, and DNA shared across genera also leave removed species behind, but their reads reach only that part of the genome, which the window spread rejects. A reference that consists mostly of foreign DNA can pass both checks and produce a false genus row. Both the even-coverage species and the genus rows are on by default for paired-end short reads and long reads; `--strict-coverage` turns both off, and the reads of removed species are then dropped. Single-end short reads are always profiled with strict coverage: there both add more false rows than they recover. Strain profiling does not use them.

### 7.11 Relative identity guard

When part of a genome is missing from the index, as in fractional indexes and for marker genomes, reads from the missing part land on related genomes at lower identity. Such a relative can pass every coverage test, because these reads spread over its sequence like its own reads. `--relative-identity-gap F` drops a reported species when a reported species of the same genus has more depth and a mean mapping identity at least `F` higher; the mean identity of a taxon is taken over every fragment mapped to it, using the fragment's best mapping, before the read-level identity gate. Species without a genus in the taxonomy are not checked. The number of removed species is logged; `--verbose` lists each with the relative that explains it.

The guard is off by default, because it cannot tell such a relative from a species that is truly present but much less abundant than a close relative: the relative's reads land on that species too and lower its mean identity, so it is removed as well. Use the guard when close relatives are not expected in the same sample.

### 7.12 Identity relative to the sample

The reads of an organism that the database lacks land on its closest relatives, at an identity below that of reads from organisms the database holds, and such a relative can pass every coverage test. Sequencing errors, read length and strain diversity lower the identity of all species in a sample alike, so Metax compares each reported species with the sample itself (`--identity-margin [F]`):

- In every fragment, each taxon's placement is scored as mapped bases × (2 × identity − 1): the pair of both mates on one reference when the taxon has one, otherwise the best single placement of each mate. The fragment counts for every taxon that reaches its top score, with the identity of that placement (weighted by mapped bases).
- For every reported species the mean identity of its fragments is taken. The baseline is the median of these means over the reported species, weighted by their fragments.
- A species whose mean lies more than `F` below the baseline is reported at its genus: its row joins the genus row, which adds up fragments, depth and abundance and shows the coverage columns of its largest part, and its classified reads follow.

The rule is on by default with `F = 0.025` for paired-end Illumina community profiling; `--no-identity-margin` turns it off. For single-end and long-read profiling it applies only when `--identity-margin` is given (without a value, `F` is `0.025`); strain profiling does not use it. The rule needs at least 3 reported species with such fragments; a species without a genus in the taxonomy stays. With the relative identity guard (§7.11) it runs after the guard. The log reports the baseline and the number of species moved; `--verbose` lists each with its mean identity.

### 7.13 Abundance estimation

Each fragment adds depth to its candidates: the mapped bases of its read (mate 1 of a pair) over the genome size, which for a fractional index is the effective indexed size. A fragment on a marker gene adds the length of its read, also when the read is clipped at the end of the gene, over `Σ placements(L_i)` for the genome's genes of lengths `L_i` (§7.4); for a long read these are the placements of its own length rather than of the median, since long reads differ widely in length, and a long read adds this for every gene of the genome it covers. A genome of `G` bases offers about `G` placements, so this is the depth of the whole genome, comparable with the depth of a genome taxon. When no gene of a marker genome can hold a read of the sample's length, its fragments add mapped bases over the summed gene length.

Reads with several retained candidates are redistributed by an expectation–maximization procedure that weights the candidates by these depths and iterates until the relative change falls below `1e-5` or `--max-em-iterations` is reached. Abundance is depth normalized across the taxa that pass the final filters. Read counts are reported in full-index equivalents.

### 7.14 Unclassified estimate

With `--estimate-unclassified`, Metax estimates how much of the input the profile does not explain. Every row explains the fragments in its column 4, which for a fractional database are already full-index equivalents. A species represented only by marker genomes explains more than its fragments on the marker genes: a genome of declared size `G` offers `placements(G)` fragment placements against `Σ placements(L_i)` for its genes of lengths `L_i`, with the placements per sequence of §7.4. The placements of the genes count a fragment on every gene it covers, so the ratio is multiplied by the gene hits per fragment of the fragments resolved to the genome (above 1 where the two mates of a pair lie on two genes). Its fragments are multiplied by that ratio, averaged over its marker genomes. Long reads differ in length and can cover several genes, so for them the ratio is measured on the reads resolved to the species' marker genomes. Only genes that need all `a` bases from a read (`⌊L/2⌋ + 1 ≥ a`, §7.4) count: a read's hits on them are found reliably, while hits on shorter genes, which need only more than half of the gene, are often missed. Each hit of a read of length `R` on such a gene stands for `placements(G) / Σ placements(L_i)` whole-genome reads at that length, summed over these genes, and the ratio is their sum over the number of reads. A species that also has genome references counts its fragments as they are. Then

```text
unclassified share = 1 − min(1, explained fragments / input fragments)
```

The estimate is taken after the relative identity guard and before species move to their genus (§7.12). Genus rows for the reads of removed species (§7.10) explain their fragments; other fragments of removed species, and fragments that fail the read-level gates, count as unexplained. When a reported marker species has no declared genome size, none of its genes can hold a read of the sample's length, or (long reads) no read is resolved to its marker genomes, the estimate is unavailable: the log says which species, and the row holds `NA`. The share is also written to the log. The option applies to community profiling and cannot be combined with `--host` or strain profiling.

## 8. Common Recipes

### 8.1 High-sensitivity Illumina run

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample_recall -t 32 --mode recall
```

### 8.2 Stricter final filtering

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample_strict -t 32 \
  --min-reads 10 --min-oebr 0.85 --min-coebr 0.85
```

`--min-reads` applies in addition to the automatic chunk-breadth cutoff.

### 8.3 Explicit chunk-breadth cutoff plus minimum reads

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample_filtered -t 32 \
  --min-reads 10 --min-chunk-breadth 0.20
```

### 8.4 Nanopore run with custom thresholds

```bash
metax profile -d db -i ont.fastq.gz -s Nanopore -o ont_custom -t 32 \
  --min-identity 0.92 --mapped-len 200 --min-mapped-fraction 0.55
```

### 8.5 Map once, profile several ways

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample -t 32 --save-mapping
metax profile --reuse-mapping sample.mapping.tsv.gz -d db -o sample_precision --mode precision
metax profile --reuse-mapping sample.mapping.tsv.gz -d db -o sample_min10 --min-reads 10
```

### 8.6 Low-biomass sample

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample_lowbiomass -t 32 --low-biomass
```

### 8.7 Marker-gene database

```bash
metax index marker_genes.fa -o db_markers --reference-type markers -t 20
metax profile -d db_markers -i R1.fastq.gz -I R2.fastq.gz -o sample_markers -t 32
```

### 8.8 Fractional index with the relative identity guard

For samples in which close relatives are not expected together (§7.11):

```bash
metax profile -d db_frac_10 -i R1.fastq.gz -I R2.fastq.gz -o sample_frac -t 32 \
  --relative-identity-gap 0.01
```

## 9. Limitations

- The built-in mapper is a coordinate mapper, not a full base-level aligner; profiling uses its identity and interval output directly, and normal runs do not produce SAM.
- Per-read classification output is unavailable for fractional databases.
- Profiling reads Metax's own coordinate mappings; SAM/BAM files from other mappers are not supported.
- Strain profiling is experimental (§5.9).
- The genus rows for the reads of removed species cannot recognise a reference that consists mostly of foreign DNA (§7.10).
- The identity margin also reports at genus a species present only as a strain that has diverged clearly from every reference of its species (§7.12).
- The unclassified estimate counts reads failing the read-level gates, and reads of filtered species that no genus row takes, as unexplained, so it errs high; for marker databases it depends on the genome sizes declared in the reference FASTA (§7.14).
- Windows is supported only through WSL.

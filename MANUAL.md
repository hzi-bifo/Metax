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
- Malformed headers, unsupported prefixes, duplicate contigs within a genome, and inconsistent genome sizes are rejected with an explanatory error.
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

### 4.3 Index options

Seed scheme:

- `-k, --seed-length`: length of each syncmer component. Default `15`.
- `-w, --seed-window`: seed density/window parameter; larger values produce fewer seeds and a smaller index. Default `15`.
- `--max-occ`: discard seed keys occurring at more than N reference positions. Default `100`.
- `--seed-mode`: `randstrobe` (default) or `minimizer`.
- `--hash-bits`: high hash bits stored as the lookup key (24–64). Default `32`.
- `--strobe-start-dist-min` / `--strobe-start-dist-max`: start-to-start distance range between the two strobes. Defaults `30` and `90`.
- `--first-strobe-only`: compatibility mode that matches only the first strobe; more sensitive but less specific than the default full-pair verification.

Resources:

- `-t, --threads`: worker threads for sampling and parallel seed generation.
- `--seed-chunk-mib`: approximate reference bases held per parallel seeding batch. Default `1024`. A single record larger than this is kept intact.
- `--temp-dir`: existing directory for large temporary seed partitions and sort stores. Fast local scratch is recommended. Without it, temporaries are created inside the hidden staging area beside `--output`. They are removed after indexing.

Progress:

- `-q, --quiet` suppresses progress and summary messages; `--progress` forces live progress output; `--no-progress` disables it; `--progress-interval N` sets the record interval for count-based updates.

### 4.4 Inspect a database

```bash
metax info -d db
metax info -i db/metax_db.mxi
```

`info -d` summarizes the database metadata and, when the `.mxi` file is present, the mapping-index header. `info -i` inspects an index directly. Neither command loads the full index into memory.

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

**Adaptive long-read identity.** For Nanopore/PacBio community profiling, the identity in the table is the high-quality default. When `--min-identity` is omitted and reads are available, Metax samples the first 10,000 reads (up to 250 million quality bases) and computes the mean base accuracy from the Phred scores. If the accuracy is below `0.95`, the identity cutoff becomes `accuracy − margin`, where the margin is `0.03` (default), `0.04` (recall), or `0.02` (precision), never below `0.84`. High-quality samples, short reads, strain mode, and an explicit `--min-identity` keep the fixed preset.

### 5.6 Mapper options

Options before `--` configure the profiler; options after `--` configure the mapper, i.e. which candidate mappings reach the profiler. Show the complete reference with:

```bash
metax profile --mapper-help
```

Resolved mapper defaults:

| Input | top-n | score ratio | min seeds | mapper identity | mapper fraction |
| --- | ---: | ---: | ---: | ---: | ---: |
| Illumina single-end / assembly | 5 | 0.90 (0.80 in recall mode) | 2 | 0.80 | 0.50 |
| Illumina paired-end | 50 | 0.60 | 1 | 0.80 | 0.50 |
| Nanopore / PacBio | 10 | 0.70 | 3 | 0.80 | 0.30, or 300 mapped bases |

The paired-end candidate set is deliberately broad so that the profiler can see competing species and keep conserved regions ambiguous; the modes differ in the evidence stage, not in candidate generation. Example overrides:

```bash
metax profile ... -- --top-n 20 --score-ratio 0.80
```

Note that `--min-identity` appears on both sides of `--` with different meanings: before it, the profiler's species/strain threshold; after it, the mapper's emission cutoff.

Paired-end emission: `--paired-emission concordant-only` (default) emits only concordant same-reference pairs; `concordant-first` emits concordant pairs when present and otherwise discordant rows; `all` emits every mate mapping. The latter two are mainly useful together with `--save-mapping` for diagnostics.

Mate rescue: a mate that seeding did not place is searched for in the insert window of its partner's placements. `--mate-rescue-taxa N` (default `8`) extends this search to fragments that already have a proper pair: for up to `N` other taxa placed on one mate only, the missing mate is searched for and kept when that pair matches at least as many bases as the best seeded pair (§7.1). `--mate-rescue-taxa 0` restricts mate rescue to fragments without any proper pair, the behaviour up to Metax 2.0.2.

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

**Reuse external alignments** in SAM, gzip-compressed SAM, or BAM format:

```bash
metax profile --reuse-sam sample.bam --dmp-dir taxonomy -s Illumina -o sample_reprofile -t 20
```

The alignments must use extended CIGAR operators that distinguish matches (`=`) from mismatches (`X`), so that identity and mapped fraction can be computed. CRAM is not supported. This path uses the legacy SAM/BAM evidence rules rather than the native coordinate-evidence policy (§7.3).

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
| 4 | number of reads, in full-index equivalents (scaled for fractional databases) |
| 5 | depth of coverage |
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

Mapping columns: read id, reference id, start, end, identity, read length, mapped bases, mate reference status, mate number (`0`, `1`, or `2`), read-local pair id, and a candidate-boundary mask (`0`–`3`) indicating whether additional endpoints may have been omitted by the candidate cap. Identity is canonicalized to five decimal places. Reference, genome, taxonomy, and size tables are never copied into the mapping file; they come from the database, whose fingerprint is validated on reuse. Headerless or obsolete mapping files are rejected.

## 7. Method

### 7.1 Mapping and candidate generation

The built-in mapper indexes references with syncmer-based randstrobe seeds and retains up to `top-n` candidate mappings per read or mate whose seed score is at least `score-ratio` of the best. Candidates are refined to reference intervals with a per-base identity: an ungapped comparison for short reads, and for long reads a gapped identity computed by exact edit distance on intervals longer than 512 bases. Mapper-level acceptance requires identity of at least `0.80` and the platform's mapped-fraction gate (§5.6).

Before applying `top-n`, the mapper reserves the strongest accepted mapping for each represented profile endpoint (normally species) and each most-specific endpoint (strain when known), then fills the remaining slots with the next best mappings. Repeated placements from one taxon therefore cannot hide another eligible taxon. A candidate-boundary mask records whether endpoints may have been omitted; the profiler keeps such candidates in abundance EM but does not force-resolve an incomplete candidate set.

Paired-end mapping is pair-aware: mates are scored together within the insert-size window, and only concordant same-reference pairs are emitted by default. A mate without a seeded placement, or a fragment without any proper pair, triggers mate rescue: the read is searched for in the insert window around its partner's placements.

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

**Paired-end Illumina.** Mate 1 must pass the identity threshold, both mates must pass the length and fraction gates, and both must belong to a concordant same-reference pair. Pair errors are the sum over both mates. A fragment contributes coordinates only when it has a single candidate endpoint or `runner_up_pair_errors − best_pair_errors ≥ 1.0`. Mate 1 supplies the mapped bases used for depth and the intervals used for coverage; mate 2 contributes to pair discrimination.

**Nanopore/PacBio.** Candidates are ranked by `long_score = mapped_bases × (2 × identity − 1)`. Coordinate evidence requires the best terminal score to exceed the runner-up by at least 1 % of the best score. Exact ties keep the mapper's first-ranked placement.

Exact and near-exact cross-species ties stay in the EM candidate set and create no coordinate evidence, unless resolved by guarded MAP (§7.8). Reads that map competitively across terminal ranks may be reported at their lowest common ancestor when that rank has sufficient retained evidence; otherwise they remain unclassified. Assembly input, strain profiling, and `--reuse-sam` use the legacy evidence rules instead of this policy.

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

### 7.5 Automatic chunk-breadth cutoff

`-c, --min-chunk-breadth` sets the fixed chunk-breadth requirement directly. When it is omitted, Metax estimates it from the read count. With `m = reads / 1,000,000`:

| Platform | cutoff |
| --- | --- |
| Illumina | `0.66·m/(m+1) + 0.24·m/(m+100)` |
| Nanopore / PacBio | `0.85·m/(m+0.78)` |

Both curves are capped at `0.95`. The read-count basis is the total input by default; with `--by-mapped`, or automatically on low-map runs (mapping rate below `--low-map-rate-threshold`, default `0.30`), it is the mapped reads/fragments.

Adjustments:

- `--chunk-breadth-scale FACTOR` multiplies the automatic cutoff (capped at `1.0`). It conflicts with `--min-chunk-breadth`, `--min-reads`, and `--low-biomass`.
- In community **recall** mode, once the basis exceeds one million reads the standard cutoff is additionally scaled by `0.70 + 0.30/m` for Illumina and by `0.80` for long reads. This scaling affects only the standard presence tier; the safeguard used by rescue tiers and guarded MAP is unaffected.
- `-n, --min-reads N` disables the automatic cutoff unless `--min-chunk-breadth` is also given. `-l, --low-biomass` disables it entirely.
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

`--min-oebr` and `--min-coebr` override the lower bounds; the upper bound stays `1.50`. `-b, --min-breadth` additionally requires a minimum observed breadth.

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

Rescue tiers apply to species only; direct higher-rank endpoints use the standard rule. Explicit `--strict`, `--min-oebr`, and `--min-coebr` remain authoritative for every tier.

**Guarded MAP.** Native Illumina community profiling resolves exact top-score cross-species ties by default (`--exact-ties guarded-map`). It seeds an EM prior from species with independent unique-candidate support, requires a single largest posterior of at least `0.95`, preserves every species selected by the margin rule, and requires independent pre-MAP fixed-chunk breadth on the same assembly: above `0.50 × cutoff` for paired-end (or `0.05 × cutoff` for a species meeting the low-abundance tier) and above the full cutoff for single-end. MAP affects coordinate evidence only; abundance still uses fractional EM over all retained candidates. Use `--exact-ties ambiguous` to keep every exact tie ambiguous. Guarded MAP is not applied to long reads, strain profiling, assembly input, or `--reuse-sam`, and it falls back to `ambiguous` when `--low-biomass` or `--min-reads` disables the chunk-breadth safeguard.

### 7.9 Genus fallback

Genus fallback assigns a read to a genus when none of its mappings qualifies for species/strain evidence. It is enabled explicitly with `-g, --genus-fallback`, or automatically for Illumina community runs without `--host` on a full (non-fractional) database when the pre-assignment mapping rate is below `--low-map-rate-threshold` (default `0.30`).

A fallback mapping must have identity of at least `--genus-identity` (default `0.80`), a mapped length of at least 50 bases for Illumina or the configured species `--mapped-len` otherwise, and pass the platform-aware genus coverage gate. The genus taxid is resolved from the species taxid through the taxonomy; species/strain evidence always takes precedence.

### 7.10 Abundance estimation

Reads with several retained candidates are redistributed by an expectation–maximization procedure that iterates until the relative change falls below `1e-5` or `--max-em-iterations` is reached. Abundance is normalized across the taxa that pass the final filters. Read counts are reported in full-index equivalents.

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

`--min-reads` without `--min-chunk-breadth` disables the automatic chunk-breadth cutoff.

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

## 9. Limitations

- The built-in mapper is a coordinate mapper, not a full base-level aligner; profiling uses its identity and interval output directly, and normal runs do not produce SAM.
- Per-read classification output is unavailable for fractional databases.
- Strain profiling is experimental (§5.9).
- Windows is supported only through WSL.

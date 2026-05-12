# Metax Manual

This manual describes how to install Metax, prepare indexes, run profiling, and interpret the main filtering behavior.

## 1. Installation

### 1.1 Install with Conda

Configure Bioconda and conda-forge:

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

Install Metax:

```bash
conda create -n metax -c zldeng metax
conda activate metax
```

### 1.2 Build from source

Install Rust and Cargo, then build:

```bash
git clone https://github.com/hzi-bifo/Metax.git
cd Metax
cargo build --release
```

The executable is written to:

```text
target/release/metax
```

Add it to your `PATH` if desired:

```bash
export PATH="$PWD/target/release:$PATH"
```

Metax calls `maCMD` for read alignment. Ensure the MA aligner is available:

```bash
conda create -n metax -c bioconda ma=1.1.4
conda activate metax
which maCMD
```

Check the Metax CLI:

```bash
metax --help
metax index --help
metax profile --help
```

## 2. Taxonomy Files

Metax needs NCBI-style taxonomy dump files. Download and unpack the NCBI taxonomy dump:

```bash
mkdir -p taxonomy
curl -L -o taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz -C taxonomy
```

The directory passed to `--dmp-dir` should contain at least:

```text
nodes.dmp
names.dmp
merged.dmp
```

Alternative taxonomy systems can be used if they are converted to the same dump-file format.

## 3. Reference Databases

You can either download the pre-built Metax database or build a custom database.

### 3.1 Download the Pre-built Metax Database

A pre-built database is available here:

```text
https://research.bifo.helmholtz-hzi.de/downloads/metax/metax_db.tar.xz
```

Download and unpack it:

```bash
mkdir -p metax_db
curl -L -o metax_db.tar.xz \
  https://research.bifo.helmholtz-hzi.de/downloads/metax/metax_db.tar.xz
tar -xJf metax_db.tar.xz -C metax_db
find metax_db -name "metax_db.json"
```

Use the discovered `metax_db.json` path as the `--db` value during profiling.

This database comprises 33,143 RefSeq genomes downloaded on 10 August 2022, spanning bacteria, archaea, viruses, fungi, protozoa, and Homo sapiens.

### 3.2 Download Reference Genomes for a Custom Index

Metax does not require a specific genome source, but the FASTA headers must be converted to the Metax format before indexing.

Common reference sources include:

- NCBI Datasets / RefSeq / GenBank
- GTDB-derived genome collections
- ICTV or virus-specific curated genome sets
- A local curated genome collection

Install the NCBI Datasets CLI separately if you use the following example. Example command for a manageable RefSeq download:

```bash
datasets download genome taxon "Escherichia coli" \
  --assembly-source RefSeq \
  --assembly-level complete,chromosome \
  --include genome \
  --filename ecoli_refseq.zip

unzip ecoli_refseq.zip -d ecoli_refseq
find ecoli_refseq/ncbi_dataset/data -name "*_genomic.fna"
```

For large taxonomic groups, first preview or summarize the download with NCBI Datasets, then select the assemblies you want to include. Record the source, date, filters, and assembly accessions for reproducibility.

### 3.3 Metax FASTA Header Format

Every sequence header in the FASTA used for indexing must follow:

```text
>genome_id|txid|species_txid|sequence_id|genome_size
```

Fields:

- `genome_id`: Unique genome or assembly identifier. All contigs from one genome should share this value.
- `txid`: Taxonomy ID for the genome/assembly.
- `species_txid`: Species-level taxonomy ID. If the assembly taxid is below species, use the species ancestor here.
- `sequence_id`: Unique sequence or contig identifier.
- `genome_size` (optional): Total genome length in bases, not just the current contig length.

Example:

```text
>GCF_000005845.2|562|562|NC_000913.3|4641652
```

If a genome has multiple contigs, use the same `genome_id`, `txid`, `species_txid`, and total `genome_size` for each contig, but use a distinct `sequence_id`.

### 3.4 Example Header Rewriting Workflow

Prepare a mapping table with one row per source FASTA:

```text
<original_header>  <target_header>
```

Then rewrite the FASTA headers with:

`seqit rename genomes.fa --map-file id_map.tsv -o renamed_genomes.fa`

 or 

`seqkit replace`

They can be downloaded/installed from here: [seqit](https://github.com/dawnmy/seqit/) and [seqkit](https://github.com/shenwei356/seqkit)

Validate the first few headers:

```bash
grep '^>' renamed_genomes.fa | head
```

## 4. Build an Index

### 4.1 Full Reference Index

Build a full index:

```bash
metax index renamed_genomes.fa -o db_full
```

The output directory contains a `metax_db.json` file. Use that JSON file for profiling:

```bash
metax profile --db db_full/metax_db.json ...
```

### 4.2 Fractional Index

If your computer/server has limited RAM or computational resources, you may consider using a fractional index to reduce memory usage and runtime. A fractional index samples segments from each genome before index construction:

```bash
metax index references.metax.fa \
  -o db_frac_10 \
  --fraction 0.10 \
  --seed 42 \
  --segment-length 50000 \
  --min-length 3000 \
  --threads 20 \
  --compress
```

Equivalent short options:

```bash
metax index references.metax.fa -o db_frac_10 -f 0.10 -s 42 -l 50000 -m 3000 -t 20 -z
```

Fractional-index options:

- `-f, --fraction`: Fraction of each genome to target. Must be greater than 0 and less than 1.
- `-s, --seed`: Random seed. Default `42`.
- `-l, --segment-length`: Segment length. Default `50000`.
- `-m, --min-length`: Minimum retained segment length for large genomes. Default `3000`.
- `-t, --threads`: Threads used during subsampling.
- `-z, --compress`: Compress the generated `metax_subsampled.fasta` after indexing.

Subsampling behavior:

- Genomes smaller than `500000 bp` are not chunk-subsampled; their full non-empty contigs are retained. Therefore, their detection is not affected by subsampling, even though they typically receive fewer supporting reads because of their smaller genome size.
- Genomes `>= 500000 bp` are sampled to at least `ceil(fraction * genome_size)` bases when possible. 
- Random candidate segments have length `--segment-length`.
- Contigs shorter than `--segment-length` are not used in the random segment pool.
- If the random segment pool cannot cover the target length, Metax falls back to longest-contig segments.
- For genomes larger than `10 * --min-length`, selected segments shorter than `--min-length` are filtered unless all selected segments would be removed.
- Subsampled headers include a sixth field containing the sampled length. This is how profiling detects a subsampled database.

During profiling with a fractional index:

- Reported read counts are scaled by `genome_size / sampled_size`.
- `classify.txt` is skipped.
- Automatic low-map genus fallback is disabled for subsampled databases unless `--genus-fallback` is explicitly supplied.
- If explicit `--genus-fallback` is used on a subsampled database and `--min-cbreadth` is not set, the chunk-breadth estimate uses total reads.

## 5. Run Profiling

### 5.1 Illumina Paired-end Reads

```bash
metax profile \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  --in-seq sample_R1.fastq.gz,sample_R2.fastq.gz \
  --is-paired \
  --sequencer Illumina \
  --outprefix sample_illumina \
  --threads 32
```

Short form:

```bash
metax profile -i sample_R1.fastq.gz,sample_R2.fastq.gz -p \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  -o sample_illumina \
  -t 32
```

### 5.2 Illumina Single-end Reads

```bash
metax profile \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  --in-seq sample.fastq.gz \
  --sequencer Illumina \
  --outprefix sample_illumina_single \
  --threads 32
```

### 5.3 ONT / Nanopore Reads

Use `--sequencer Nanopore`:

```bash
metax profile \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  --in-seq ont_reads.fastq.gz \
  --sequencer Nanopore \
  --outprefix sample_ont \
  --threads 32
```

### 5.4 PacBio Reads

```bash
metax profile \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  --in-seq pacbio_reads.fastq.gz \
  --sequencer PacBio \
  --outprefix sample_pacbio \
  --threads 32
```

### 5.5 Reuse an Existing SAM/BAM File

If reads have already been aligned to a Metax-formatted reference:

```bash
metax profile \
  --reuse-sam sample.sam \
  --dmp-dir taxonomy \
  --sequencer Illumina \
  --is-paired \
  --outprefix sample_reprofile \
  --threads 32
```

Metax supports reuse of external alignments in SAM, compressed SAM, and BAM formats; CRAM is not currently supported. External alignment files must contain extended CIGAR strings that distinguish matches (=) from mismatches (X), rather than using M for both. This is required for accurate calculation of alignment identity, aligned fraction, and downstream coverage-based statistics.


### 5.6 Resume and Compression

Resume from an existing `<outprefix>.sam` or `<outprefix>.sam.gz`:

```bash
metax profile \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  -i sample_R1.fastq.gz,sample_R2.fastq.gz \
  -p \
  -o sample_illumina \
  --resume
```

Compress the generated SAM after profiling:

```bash
metax profile ... --compress-sam
```

### 5.7 Profiling Modes

Choose a preset with:

```bash
metax profile ... --mode default
metax profile ... --mode recall
metax profile ... --mode precision
```

Default alignment-filter thresholds:

| Sequencer / mode | identity | mapped length | aligned fraction |
| --- | ---: | ---: | ---: |
| Illumina default | 0.97 | 50 | 0.70 |
| Illumina recall | 0.95 | 50 | 0.60 |
| Illumina precision | 0.98 | 50 | 0.80 |
| Nanopore default | 0.86 | 250 | 0.60 |
| Nanopore recall | 0.85 | 250 | 0.50 |
| Nanopore precision | 0.90 | 250 | 0.60 |
| PacBio default | 0.86 | 250 | 0.60 |
| PacBio recall | 0.85 | 250 | 0.50 |
| PacBio precision | 0.90 | 250 | 0.60 |

Override these with:

```bash
metax profile ... --identity 0.95 --mapped-len 100 --fraction 0.65
```

## 6. Output Files

For an output prefix `sample`, Metax writes:

- `sample.sam`: Alignment output from maCMD, unless `--reuse-sam` was used.
- `sample.profile.txt`: Final taxonomic profile.
- `sample.classify.txt`: Per-read classifications, skipped for subsampled databases.
- `sample.log`: Run log.

If `--host` is used, profile files use the prefix `sample.pathogen`.

### 6.1 Profile Columns

The profile output is tab-separated and does not currently include a header row.

Standard columns:

| Column | Meaning |
| ---: | --- |
| 1 | taxon name |
| 2 | taxid |
| 3 | rank |
| 4 | read count, scaled for subsampled databases |
| 5 | depth |
| 6 | abundance |
| 7 | observed breadth |
| 8 | expected breadth |
| 9 | breadth p-value field (`cov_prob`) |
| 10 | fixed chunk breadth |
| 11 | flexible chunk breadth |
| 12 | expected flexible chunk breadth |
| 13 | flexible chunk p-value field (`chunk_prob`) |

When pathogen-host filtering is enabled, three columns are appended:

| Column | Meaning |
| ---: | --- |
| 14 | host names |
| 15 | host taxids |
| 16 | diseases |

The final `profile.txt` abundance column is renormalized over taxa that pass final filtering.

## 7. Method
### 7.1 Alignment and Read Assignment

Metax runs maCMD unless `--reuse-sam` is supplied. The maCMD preset is selected from `--sequencer`:

- Illumina single-end: `Illumina`
- Illumina paired-end: `Illumina_Paired`
- Nanopore: `Nanopore`
- PacBio: `PacBio`
- Other values: `Default`

The profiler parses alignments, computes identity and aligned fraction from CIGAR, filters alignments by `identity`, `mapped_len`, and `fraction`, then assigns reads to taxa. Multi-mapping reads are resolved with an EM procedure.

### 7.2 Coverage Metrics

Metax reports:

- `breadth`: fraction of genome covered by at least one read.
- `depth`: average coverage depth.
- `fixed_chunk`: breadth of fixed-size chunks around alignments.
- `flex_chunk`: breadth of variable-size chunks. The flexible chunk count is based on approximately `sqrt(genome_length) / 2`.
- `expected_breadth`: expected breadth under the read/depth model.
- `expected_flex_chunk`: expected flexible chunk breadth.

### 7.3 Minimum Breadth

`--min-breadth` filters taxa by observed breadth before abundance estimation.

Default:

```text
--min-breadth 0.0
```

### 7.4 Minimum Chunk Breadth

`--min-cbreadth` sets the fixed chunk-breadth threshold directly. If it is unset, Metax estimates it from read count.

Estimated Illumina threshold (determined based on the CAMI HMP toy dataset):

| Read count basis | threshold |
| ---: | ---: |
| `<= 100000` | `reads / 1000000` |
| `100001` to `1000000` | `0.11 * reads_in_million + 0.09` |
| `> 1000000` | `min(0.2 + 0.382 * log10(reads_in_million), 0.95)` |

Estimated long-read threshold:

| Read count basis | threshold |
| ---: | ---: |
| `<= 5000` | `0.0` |
| `5001` to `1000000` | `0.3` |
| `> 1000000` | `0.5` |

Read-count basis:

- Default: total reads.
- `--by-aligned`: aligned reads.
- Automatic low-map Illumina non-pathogen runs: aligned reads.
- Low-map means `aligned_reads / total_reads < --low-map-rate-threshold`.
- Default `--low-map-rate-threshold`: `0.30`.

If `--min-reads` is set and `--min-cbreadth` is not set, the automatic fixed chunk-breadth filter is skipped. If both are set, the explicit `--min-cbreadth` filter is still applied.

`--lowbiomass` makes the estimated chunk threshold `0.0`.

### 7.5 Minimum Reads

`--min-reads` filters the final profile by the read-count column:

```text
reported_count >= --min-reads
```

For subsampled databases, this count is already scaled by `genome_size / sampled_size`.

### 7.6 Observed/Expected Ratio Filters

Metax applies two ratio filters:

```text
oebr  = observed_breadth / expected_breadth
coebr = observed_flex_chunk / expected_flex_chunk
```

Community mode defaults:

| Rank | default lower bound | upper bound |
| --- | ---: | ---: |
| species and other non-genus ranks | `0.75` | `1.5` |
| genus | `0.65` | `1.5` |

The upper bound is inclusive in community mode:

```text
ratio <= 1.5
```

`--min-oebr` and `--min-coebr` override the lower bounds independently:

```text
oebr  >= --min-oebr
coebr >= --min-coebr
```

The maximum remains `1.5`. If a minimum is set above `1.5`, no row can pass that ratio.

### 7.7 P-value Logic

Species rows can receive two p-value fields:

- `cov_prob`: p-value for observed breadth relative to expected breadth.
- `chunk_prob`: p-value for observed flexible chunk breadth relative to expected flexible chunk breadth.

Rows outside species rank generally have `NA` p-value fields.

P-value calculation:

- If the observed/expected ratio is within `[0.9, 1.2]`, the p-value is set to `1.0`.
- Otherwise, Metax tests the lower or upper tail depending on whether observed breadth is below or above expected breadth.
- P-values are multiplied by the number of taxa considered and capped at `1.0`.
- The implementation chooses among exact Stirling, normal approximation, and binomial approximation according to sample size and expected breadth.

Because exact Stirling-number-based p-value calculation is computationally intensive, version below v0.9.22 calculated this p-value only for taxa with 2–100 supporting reads and OEBR or COEBR below 0.75. For all other taxa, the p-value was reported as NA. And the filter merely relies on OEBRs.

Community final p-value filter:

- Cutoff: `1e-5`.
- Default mode keeps a row when at least one available p-value is `>= 1e-5`.
- If both available p-values are `< 1e-5`, the row is removed.
- With `--strict`, both available p-values must be `>= 1e-5`.
- In strict mode, genus rows with both p-values absent can still proceed to the other final filters.


### 7.8 Genus Fallback

Genus fallback is used by the community profiler.

Automatic activation requires:

- `--sequencer Illumina`
- no `--host`
- `aligned_reads / total_reads < --low-map-rate-threshold`
- full, non-subsampled database

This fallback is intended for cases in which the reference database does not adequately represent the species present in the sample, resulting in a low mapping rate. Under these conditions, higher-rank profiling can be more informative. By contrast, pathogen detection (--host) typically requires more specific, lower-rank profiling.

The default low-map threshold is:

```text
--low-map-rate-threshold 0.30
```

Explicit activation:

```bash
metax profile ... --genus-fallback
```

Candidate alignment requirements:

- No species-qualified hit exists for the read.
- A genus taxid can be resolved from the species taxid.
- `identity >= --genus-identity`.
- Default `--genus-identity`: `0.80`.
- Aligned fraction is at least `0.60`.
- Mapped length is at least `50` for Illumina.
- For non-Illumina, mapped length follows the configured species `mapped_len` threshold.

Species-qualified hits take precedence over genus fallback.

### 7.9 Host / Pathogen Mode

With `--host`, Metax writes `<outprefix>.pathogen.profile.txt` and filters/annotates rows using the pathogen-host table when supplied:

```bash
metax profile ... --host 9606 --pathogen-host pathogen_hosts.tsv
```

The pathogen-host table is tab-separated with a header:

```text
pathogen_taxid	host_taxids	host_names	diseases
```

Multiple host taxids, host names, or diseases can be semicolon-separated.

For viral pathogens, this table can be compiled from host-association records in the [Virus-Host Database](https://www.genome.jp/ftp/db/virushostdb/).

## 8. Common Recipes

### 8.1 High-sensitivity Illumina run

```bash
metax profile \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  -i R1.fastq.gz,R2.fastq.gz \
  -p \
  -o sample_recall \
  -t 32 \
  --mode recall
```

### 8.2 Stricter final profile filtering

```bash
metax profile \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  -i R1.fastq.gz,R2.fastq.gz \
  -p \
  -o sample_filtered \
  --min-reads 10 \
  --min-oebr 0.85 \
  --min-coebr 0.85
```

When `--min-reads` is set without `--min-cbreadth`, the automatic fixed chunk-breadth filter is skipped.

### 8.3 Explicit chunk-breadth filtering plus minimum reads

```bash
metax profile \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  -i R1.fastq.gz,R2.fastq.gz \
  -p \
  -o sample_filtered2 \
  --min-reads 10 \
  --min-cbreadth 0.20
```

In this case, both the minimum read count and the explicit chunk-breadth threshold are applied.

### 8.4 ONT run with custom alignment thresholds

```bash
metax profile \
  --db db_full/metax_db.json \
  --dmp-dir taxonomy \
  -i ont.fastq.gz \
  --sequencer Nanopore \
  -o ont_custom \
  --identity 0.92 \
  --mapped-len 200 \
  --fraction 0.55
```

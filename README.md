# Metax

Metax is a cross-domain metagenomic taxonomic profiler designed to deliver accurate, robust, and interpretable community composition analyses across bacteria, archaea, eukaryotes, and viruses. Unlike existing profilers, Metax integrates probabilistic modeling of genome coverage to distinguish true community members from artifacts caused by reference contamination, local genomic similarity, or reagent-derived DNA fragments.

Through comprehensive benchmarks on more than 600 samples, Metax demonstrated:

🧬 Species-level accuracy across all domains of life

⚡ Robustness to shallow sequencing and low-biomass, host-dominated samples

🔍 Contamination detection, including reagent-borne DNA and reference misassemblies

🦠 Clinical and environmental relevance, e.g. enables identifying cross kingdom interactions and clarifying tumor microbiome signals

By unifying coverage-informed presence probability estimation with EM-based abundance refinement, Metax overcomes challenges of ambiguous read mapping, database contamination and kitome DNA fragments. These properties make it a powerful tool for microbiome research, clinical metagenomics, and identifying genome misassemblies.

## What's new in Metax 2

- **Built-in read mapper.** Metax no longer depends on the external MA aligner (`maCMD`). Indexing, mapping, and profiling run inside a single `metax` executable.
- **Native short- and long-read evidence.** Single-end and paired-end Illumina, Nanopore, and PacBio reads are profiled with platform-specific coordinate-evidence rules; SAM/BAM parsing is retained only for reusing external alignments.
- **Rank-aware reference labels.** References may be labelled at strain, species, genus, or any higher rank, and Metax reports each genome at the lowest rank its label supports.
- **Lower peak memory.** Mapping and profiling run as separate stages, so the mapping index is released before the profiler starts.
- **Reusable mapping output.** Save compressed coordinate mappings once and re-profile them with different thresholds without remapping.
- **Managed taxonomy.** The NCBI taxonomy dump is downloaded and cached automatically when no `--dmp-dir` is given.
- **More platforms.** Conda packages are published for Linux x86-64, Linux ARM64, and macOS Apple Silicon.

Databases built for Metax 0.9 use a different index format and must be rebuilt with `metax index` (see [Reference database](#reference-database)).

## Installation

### Install with conda

```shell
conda create -n metax -c zldeng metax
conda activate metax
```

Packages are available for `linux-64`, `linux-aarch64`, and `osx-arm64`. The package is self-contained: no separate aligner or runtime library is required.

### Compile from source

Metax requires Rust 1.89 or newer and a C compiler (used for the bundled htslib).

```bash
cargo build --release
```

This produces the binary at `target/release/metax`. Run the bundled smoke tests to verify the build:

```bash
smoke_tests/run_smoke_tests.sh
```

## Download databases

### Taxonomy dmp files

Metax uses NCBI-style taxonomy dump files (`nodes.dmp`, `names.dmp`, `merged.dmp`, `delnodes.dmp`).

- If `--dmp-dir` is omitted, Metax downloads the official NCBI [taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) to `~/.metax/taxdmp` on first use and reuses it afterwards.
- Pass `--dmp-dir <dir>` to use an offline installation, a pinned taxonomy snapshot, or an alternative taxonomy (e.g. GTDB or ICTV) converted to the same dump-file format. An explicit directory is never modified by Metax.

### Reference database

Databases built for Metax 0.9 are not compatible with Metax 2 and must be rebuilt.

To build a custom database:

1. Prepare the reference genomes in FASTA format. Each sequence header must follow

   ```text
   >genome_id|sequence_id|taxonomy[|genome_size]
   ```

   where `taxonomy` lists the lowest confidently known rank(s) as `prefix:taxid` pairs separated by `;`, for example `t:1234;s:9606` (strain and species), `s:562` (species only), or `g:561` (genus only). Supported prefixes are `t` (strain/subspecies), `s`, `g`, `f`, `o`, `c`, `p`, and `k`. All contigs of one genome share the same `genome_id`; `genome_size` is optional.

2. Build the database:

   ```shell
   metax index <fasta_file> -o <database_dir> -t 20
   ```

   Gzip-compressed input is accepted. For machines with limited RAM, `-f <fraction>` builds a **fractional index** from an evenly distributed fraction of each large genome; see the [manual](MANUAL.md#42-fractional-index).

## How to run profiling

```shell
metax profile --help
```

Paired-end Illumina reads:

```shell
metax profile -d <database_dir> \
    -i <r1.fq.gz> -I <r2.fq.gz> \
    -o <output_prefix> \
    -t 20
```

Single-end Illumina reads:

```shell
metax profile -d <database_dir> -i <reads.fq.gz> -o <output_prefix> -t 20
```

Nanopore or PacBio reads:

```shell
metax profile -d <database_dir> -i <reads.fq.gz> -s Nanopore -o <output_prefix> -t 20
```

Add `--dmp-dir <dump_dir>` to use a specific taxonomy directory. The first run against a database builds a memory-mapped index cache beside it; subsequent runs start faster.

## Test datasets

- [CAMI II marine](https://frl.publisso.de/data/frl:6425521/marine/short_read/)
- [CAMI II pathogen detection](https://frl.publisso.de/data/frl:6425521/patmgCAMI2.tar.gz)
- [Other benchmark datasets](https://research.bifo.helmholtz-hzi.de/downloads/metax/benchmark_datasets/)

## Output

- Final taxonomy profile: `<output_prefix>.profile.txt`

  ```text
  column 1:  Taxon name
  column 2:  Taxon ID
  column 3:  Taxon rank
  column 4:  Number of reads (full-index equivalents)
  column 5:  Depth of coverage
  column 6:  Relative abundance (%)
  column 7:  Breadth of coverage
  column 8:  Expected breadth of coverage
  column 9:  P-value of the breadth deviation
  column 10: Fixed chunk breadth of coverage
  column 11: Flex chunk breadth of coverage
  column 12: Expected flex chunk breadth of coverage
  column 13: P-value of the flex chunk breadth deviation
  ```

  When `--host` is used, three columns are appended:

  ```text
  column 14: Host names
  column 15: Host taxonomy IDs
  column 16: Associated diseases
  ```

- Run log: `<output_prefix>.log`

Optional outputs:

- `--save-raw-profile`: `<output_prefix>.rprofile.txt`, the post-EM profile before taxon-level filters.
- `--save-classification`: `<output_prefix>.classify.txt`, per-read taxonomic assignments.
- `--save-mapping`: `<output_prefix>.mapping.tsv.gz`, the coordinate mappings for `--reuse-mapping`.

## Manual

For detailed installation, indexing, profiling, output, and method descriptions, see [MANUAL.md](MANUAL.md).

## FAQ

1. **Which platforms are supported?**
   Conda packages are provided for Linux x86-64, Linux ARM64, and macOS Apple Silicon. Other Unix-like systems can build from source. On Windows, use WSL.

2. **Can I use a database built for Metax 0.9?**
   No. Metax 2 uses a different index format. Rebuild the database with `metax index`.

3. **Why is the first profiling run against a database slower?**
   The first run creates a memory-mapped cache (`metax_db.mxi.metamcache`) beside the index. Later runs open the completed cache directly.

4. **Can I reuse alignments produced by another tool?**
   Yes, with `--reuse-sam`. The SAM/BAM file must use extended CIGAR operators (`=`/`X`) so that identity can be computed. For repeated analyses, prefer `--save-mapping` and `--reuse-mapping`, which avoid remapping.

5. **A species appears in `rprofile.txt` but not in `profile.txt`. Why?**
   `rprofile.txt` is written before the coverage-evidence and support filters. The columns of that row show which criterion the taxon failed.

## Acknowledgement

We thank Gary Robertson for IT support, Dr. Mohammad-Hadi Foroughmand-Araabi for advice on statistical formulations, and Hesham Almessady for software testing.

## AI-assisted development

AI coding agents were used to assist with code refinement, comments, and documentation. The design and implementation of the software remain the responsibility of the authors, and all AI-assisted changes were reviewed and tested by the authors.

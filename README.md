# Metax

Metax is a cross-domain metagenomic taxonomic profiler for bacteria, archaea, eukaryotes and viruses. It maps short or long reads to reference genomes and uses how evenly the reads cover each genome to separate true community members from artifacts of reference contamination, local genomic similarity and reagent-derived DNA. Reads that align to several taxa are shared among them by a likelihood-based EM.

## Installation

### Install with conda

```shell
conda create -n metax -c zldeng metax
conda activate metax
```

Packages are available for `linux-64`, `linux-aarch64` and `osx-arm64`. The package is a single self-contained executable.

### Compile from source

Metax requires Rust 1.89 or later and a C compiler (for the bundled zstd library).

```bash
cargo build --release
smoke_tests/run_smoke_tests.sh
```

The binary is `target/release/metax`.

## Databases

### Pre-built databases

| Database | Contents | Archive | Index cache |
| --- | --- | --- | --- |
| [RefSeq](https://s3.bi.denbi.de/swift/v1/cami3/metax/metax2.1_refseq_202607_k16w37.tar.gz) | 29,112 genomes of bacteria, archaea, viruses, fungi, protozoa and the human reference genome | 40 GiB | 42 GiB |
| [GTDB marker genes](https://s3.bi.denbi.de/swift/v1/cami3/metax/metax2.1_gtdb_r232_reps_markers_k16w37.tar.gz) | the marker genes of the 199,923 species representatives of GTDB release 232 | 12 GiB | 15 GiB |
| [GTDB fractional](https://s3.bi.denbi.de/swift/v1/cami3/metax/metax2.1_gtdb_r232_reps_f010_k16w37.tar.gz) | a tenth of each of the same 199,923 genomes, small genomes whole | 32 GiB | 34 GiB |

Download an archive, unpack it and pass the unpacked folder to `-d`:

```shell
curl -LO https://s3.bi.denbi.de/swift/v1/cami3/metax/metax2.1_refseq_202607_k16w37.tar.gz
tar -xzf metax2.1_refseq_202607_k16w37.tar.gz
metax profile -d metax2.1_refseq_202607_k16w37 -i reads_1.fq.gz -I reads_2.fq.gz -o sample -t 20
```

The GTDB databases are used the same way, for example with `-d metax2.1_gtdb_r232_reps_markers_k16w37`. Each archive holds the taxonomy dump its database was built with, NCBI for RefSeq and GTDB for the GTDB databases, which Metax reads in place, so no `--dmp-dir` is needed. The first run writes the index cache into the folder, which must be writable; to place the cache elsewhere, see the [manual](MANUAL.md#31-pre-built-databases).

RefSeq holds up to three genomes per bacterial or archaeal species, so `--strain` can tell the genomes of a species apart. The GTDB databases represent each bacterial and archaeal species by one genome, its GTDB representative. The marker-gene database is the smallest; it reports species and the levels above them. On the fractional database, reads not placed at a species are not reported above species. The [manual](MANUAL.md#31-pre-built-databases) describes each database.

### Custom database

Each sequence header of the reference FASTA follows

```text
>genome_id|sequence_id|taxonomy[|genome_size]
```

where `taxonomy` lists the lowest known ranks as `prefix:taxid` pairs separated by `;`, for example `t:1234;s:9606`, `s:562` or `g:561`. Build the database with:

```shell
metax index references.fa -o database --dmp-dir taxonomy -t 20
```

The manual describes [fractional indexes](MANUAL.md#42-fractional-index), which index a share of each large genome, and [marker-gene and mixed databases](MANUAL.md#43-marker-gene-and-mixed-databases).

## Run profiling

```shell
# Paired-end Illumina reads
metax profile -d database -i r1.fq.gz -I r2.fq.gz -o sample -t 20

# Single-end Illumina reads
metax profile -d database -i reads.fq.gz -o sample -t 20

# Nanopore or PacBio reads
metax profile -d database -i reads.fq.gz -s Nanopore -o sample -t 20
```

Useful options:

- `--strain`: report the genomes within each kept species.
- `--save-decisions` (`-k`): write how every judged taxon was decided.
- `--save-classification`: write per-read assignments.
- `--save-mapping`, then `--reuse-mapping`: map once and profile several ways.
- `--host` with `--pathogen-table`: filter and annotate the profile for one host.
- `--estimate-unclassified`: add the share of the input that the profile does not explain.

Run `metax profile --help` for every option.

## Output

- `<prefix>.profile.txt`: the profile, with two header lines and one row per reported taxon. Its columns are name, taxid, rank, reads, depth, abundance, seven coverage columns, the depth of the coverage evidence, the best-fit genome, its mean identity and the estimated identity of the sample's strain to it. `--host` adds three annotation columns.
- `<prefix>.log`: the run log.
- `<prefix>.decisions.tsv`, `<prefix>.classify.txt.gz` and `<prefix>.mapping.tsv.gz`: optional tables, written with the options above.

## Manual

[MANUAL.md](MANUAL.md) describes installation, databases, every option, the output files and the method.

## FAQ

1. **Which platforms are supported?**
   Conda packages are provided for Linux x86-64, Linux ARM64 and macOS Apple Silicon. Other Unix-like systems can build from source. On Windows, use WSL.

2. **Why is the first profiling run against a database slower?**
   The first run creates a memory-mapped cache (`metax_db.mxi.metamcache`) beside the index. Later runs open it directly.

3. **Can I reuse alignments produced by another tool?**
   No. Metax profiles its own coordinate mappings. Save them with `--save-mapping` and profile them again with `--reuse-mapping`.

4. **A species is not in `profile.txt`. Why?**
   Run with `--save-decisions`: `<prefix>.decisions.tsv` lists every judged taxon with the tests it passed or failed.

## Acknowledgement

We thank Gary Robertson for IT support, Dr. Mohammad-Hadi Foroughmand-Araabi for advice on statistical formulations, and Hesham Almessady for software testing.

## Citations

If you use Metax in your research, you can cite the following references.
- Deng, Zhi-Luo, Safaei, Nasim, and McHardy, Alice Carolyn (**2026**). Metax enables accurate cross-domain taxonomic profiling of metagenomes. Cell, 189(22), 7196–7211.e7. DOI 10.1016/j.cell.2026.08.024.
- Deng, Zhi-Luo (**2023**). Metax: Improved taxonomic profiling of low-ranking taxa using coverage information (Version 0.80) [Software]. Zenodo. DOI 10.5281/zenodo.10085101.

## AI-assisted development

AI coding agents were used to assist with code refinement, comments, and documentation. The design and implementation of the software remain the responsibility of the authors, and all AI-assisted changes were reviewed and tested by the authors.

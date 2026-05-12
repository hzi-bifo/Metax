# Metax

Metax is a cross‑domain metagenomic taxonomic profiler designed to deliver accurate, robust, and interpretable community composition analyses across bacteria, archaea, eukaryotes, and viruses. Unlike existing profilers, Metax integrates probabilistic modeling of genome coverage to distinguish true community members from artifacts caused by reference contamination, local genomic similarity, or reagent‑derived DNA fragments.

Through comprehensive benchmarks on more than 600 samples, Metax demonstrated:

🧬 Species‑level accuracy across all domains of life

⚡ Robustness to shallow sequencing and low‑biomass, host‑dominated samples

🔍 Contamination detection, including reagent‑borne DNA and reference misassemblies

🦠 Clinical and environmental relevance, e.g. enables identifying cross kingdom interactions and clarifying tumor microbiome signals

By unifying coverage‑informed presence probability estimation with EM‑based abundance refinement, Metax overcomes challenges of ambiguous read mapping, database contamination and kitome DNA fragments. These properties make it a powerful tool for microbiome research, clinical metagenomics, and identifying genome misassemblies.



## Installation

### Install the package using conda:

  Ensure that `bioconda` is included in the Conda channel source file (`~/.condarc`):
  ```
  channels:
     - conda-forge
     - bioconda
  channel_priority: strict
  ```
  
   ```shell
   conda create -n metax -c zldeng metax
   conda activate metax
   ```


### Compile from source

You can build the executable directly with Cargo:

```bash
cargo build --release
```

This produces the binary at `target/release/metax`.

- Ensure the MA aligner (v1.1.4) is available in your PATH. You can install it with `Conda`.

```
conda create -n metax -c bioconda ma=1.1.4
conda activate metax
```

## Download databases

### Taxonmy dmp files
  1. Download the NCBI taxonomy dump (`taxdump.tar.gz`) from [NCBI](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz).
  2. _Optional:_ You can also use an alternative taxonomy source (e.g. GTDB or ICTV), by generating your own dmp files.
  
### Reference database

A pre-built reference database is available at [here](https://research.bifo.helmholtz-hzi.de/downloads/metax/metax_db.tar.xz). It is based on the RefSeq snapshot of 10 August 2022 and includes top genomes (selected with [genome_updater](https://github.com/pirovc/genome_updater)) for each species, prioritizing assemblies flagged as “representative” or “reference” and then selecting the highest assembly level (Complete Genome > Chromosome > Scaffold > Contig). In total, it contains 33,143 genomes from bacteria, archaea, viruses, fungi, protozoa, and Homo sapiens. Another pre-built reference [database](https://research.bifo.helmholtz-hzi.de/downloads/cami/metax_db.tar.xz) is available for CAMI II data benchmarks.


Users can also build a customized reference database by following steps:

1. Prepare the genomes in FASTA format, the header of each sequence should be in the format:
    ```
    >genome_id|txid|species_txid|sequence_id|genome_size
    ```
    Each genome must have a unique genome_id, and each sequence a unique sequence_id. The genome_size field is optional if you don't need subsampling. When using the NCBI taxonomy, txid should be the genome’s NCBI Taxonomy ID, and species_txid the species’ NCBI Taxonomy ID. If you choose a different taxonomy source (e.g. GTDB, ICTV), use the corresponding IDs from your taxonomy dump files.

2. Build the database directly with `metax index`:
    
    ```shell
    metax index <fasta_file> -o <database_dir>
    ```
    The command produces a database named `metax_db` in `<database_dir>`. It may take long to complete, please run it in a `tmux`  or `screen` session.

    Metax supports a **fractional index** for users with limited RAM or CPU resources. A fractional index is built from only an evenly distributed fraction of each reference genome, rather than the full genome sequence. This reduces index size, peak memory usage, and runtime.

    Use `-f <fraction>` during indexing, where `0 < fraction < 1`. For example, `-f 0.1` indexes approximately 10% of each genome, while `-f 0.05` indexes approximately 5%.

    Metax selects non-overlapping genome segments of length `-l/--segment-length` (`50 kb` by default) until the requested fraction is reached. Segments are distributed across the genome to retain representative genome-wide coverage information. The fractional FASTA is written alongside the index.

    During profiling, reported counts are scaled according to the reduced effective genome size. Read-level classification output is disabled when a fractional index is used.

    Useful options:

    - `-s/--seed` (`42` by default): reproducible segment selection.

    - `-m/--min-length` (`3 kb` by default): remove short selected segments from sufficiently long genomes.

    - `-t/--threads`: parallelize fractional reference construction.

    - `-z/--compress`: gzip the fractional FASTA after indexing.


## How to run profiling

You can get the help message for profiling by:

```shell
metax profile --help
```

A typical command to run `metax profile`:

```shell
metax profile --dmp-dir <dump_dir> \
    --db <reference_db> \
    -i <r1>[,<r2>] \
    -p \
    -o <output_prefix> \
    [other options ...]
```
`<dump_dir>`: path to the folder containing all the taxonomy dump files.
`<reference_db>`: path to the json file of the database (i.e. `<database_dir>/metax_db.json`)
The first run (sample) takes a bit longer; subsequent runs will be substantially faster by using the cached database.


## Test datasets
   - [CAMI II marine](https://frl.publisso.de/data/frl:6425521/marine/short_read/)
   - [CAMI II pathogen detecton](https://frl.publisso.de/data/frl:6425521/patmgCAMI2.tar.gz)

   - [Other benchmark datasets](https://research.bifo.helmholtz-hzi.de/downloads/metax/benchmark_datasets/)


## Output

- Final taxonomy profile: `*.profile.txt`
```
column 1: Taxon name
column 2: Taxon ID
column 3: Taxon rank
column 4: Number reads
column 5: Depth of coverage
column 6: Abundance
column 7: Breadth of coverage (B)
column 8: Expected breadth of coverage (EB)
column 9: Likelihood of presence based on breadth
column 10: Fixed chunk breadth of coverage
column 11: Flex chunk breadth of coverage
column 12: Expected flex chunk breadth of coverage (ECB)
column 13: Likelihood of presence based on flex chunk breadth
```

If pathogen detection mode is enabled, the output profile will also include 3 extra columns as below:

```
column 14: The host names
column 15: The host taxonomy IDs
column 16: The relevant diseases
```

- Reads taxonomy classification: `*.classify.txt`

```
column 1: Read name
column 2: taxonomy ID of the most likely taxon
column 3: Rank of the taxon
column 4: Taxonomy IDs of all possible taxa
column 5: Likelihood for each of those possible taxa
```

## Manual

For a detailed installation, indexing, profiling, and method details, see [MANUAL.md](MANUAL.md).

## FAQ

1. What platforms and operating systems does Metax support?

Metax currently supports Linux on x86-64 (64-bit Intel/AMD) systems. Other architectures (e.g., ARM/macOS) are not yet officially supported.

2. Why do I get the error: “Processor 6174 is not supported by this build”?

This error indicates that your CPU does not support some modern instruction sets required by Metax.


## Acknowledgement

We thank Gary Robertson for IT support, Dr. Mohammad-Hadi Foroughmand-Araabi for advice on statistical formulations, and Hesham Almessady for software testing.
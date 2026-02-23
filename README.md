# Metax

Metax is a cross‑domain metagenomic taxonomic profiler designed to deliver accurate, robust, and interpretable community composition analyses across bacteria, archaea, eukaryotes, and viruses. Unlike existing profilers, Metax integrates probabilistic modeling of genome coverage to distinguish true community members from artifacts caused by reference contamination, local genomic similarity, or reagent‑derived DNA fragments.

Through comprehensive benchmarks on more than 500 samples, Metax demonstrated:

🧬 Species‑level accuracy across all domains of life

⚡ Robustness to shallow sequencing and low‑biomass, host‑dominated samples

🔍 Contamination detection, including reagent‑borne DNA and reference misassemblies

🦠 Clinical and environmental relevance, e.g. enables identifying cross kingdom interactions and clarifying tumor microbiome signals

By unifying coverage‑informed presence probability estimation with EM‑based abundance refinement, Metax overcomes challenges of ambiguous read mapping, database contamination and kitome DNA fragments. These properties make it a powerful tool for microbiome research, clinical metagenomics, and identifying genome misassemblies.


## Installation
- Install the package using conda:

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
   
  
<!-- - Alternatively, you can install it manually
   
   Download the pre-built binary from the [releases page](https://github.com/dawnmy/Metax/releases), unpack it, and add the directory to your PATH.
   
   Install the dependencies:
    ```shell
    conda install -c bioconda ma=1.1.4
    ``` -->

## Download databases

- Taxonmy dmp files
  1. Create a `metax_dmp` directory.  
  2. Download the NCBI taxonomy dump (`taxdump.tar.gz`) from: [NCBI](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz).
  3. Extract its contents directly into `metax_dmp`.  
  4. _Optional:_ To use an alternative taxonomy source (e.g. GTDB or ICTV), replace the extracted `taxdump` files in `metax_dmp` with your own dmp files.
  
- Reference database

A pre-built reference database is available at [here](https://research.bifo.helmholtz-hzi.de/downloads/metax/metax_db.tar.xz). It is based on the RefSeq snapshot of 10 August 2022 and includes top genomes for each NCBI taxonomic identifier (txid), prioritizing assemblies flagged as “representative” or “reference” and then selecting the highest assembly level (Complete Genome > Chromosome > Scaffold > Contig). In total, it contains 33,143 genomes from bacteria, archaea, viruses, fungi, protozoa, and Homo sapiens (bavfph).

Another pre-built reference [database](https://research.bifo.helmholtz-hzi.de/downloads/cami/metax_db.tar.xz) is available for CAMI II data benchmarks.


A customized reference database can be created by following steps:

1. Prepare the genomes in fasta format, the header of each sequence should be in the format:
    ```
    >genome_id|txid|species_txid|sequence_id[|genome_size]
    ```
    Each genome must have a unique genome_id, and each sequence a unique sequence_id. The genome_size field is optional but necessary for subsampled reference database creation using `metax index`. When using the NCBI taxonomy, `txid` should be the genome’s NCBI taxonomy ID, and `species_txid` the species’ NCBI taxonomy ID. If you choose a different taxonomy source (e.g. GTDB, ICTV), use the corresponding IDs from your taxonomy dump files.

2. Run the following command to build the database:

    It may take long to complete, please run it in a Tmux session or screen.

    ### For metax version >=0.9.12:

    ```shell
    metax index <fasta_file> -o <database_dir>
    ```
    This command produces a database named `metax_db` in `<database_dir>`.

    To subsample each genome before indexing, provide `-f <fraction>` (where `0 < fraction < 1`). This enables the creation of a database that uses less space and memory while also reducing profiling runtime. But it might      be less sensitive for low read count taxa. The CLI will extract
    evenly distributed, non-overlapping segments of length `-l/--segment-length` (50 Kbp by default) across each genome
    until at least the requested fraction is collected, write a subsampled FASTA alongside the index, and skip generating
    read-level classifications while scaling reported counts during profiling to account for the reduced genome size.
    Use `-s/--seed` (default `42`) to make the segment selection reproducible, and `-m/--min-length` (default 3 Kbp) to
    discard subsampled segments shorter than the threshold for genomes more than ten times longer than that value. Supply
    `-t/--threads` to subsample genomes in parallel; omit it (or set it to 1) to run sequentially. Add `-z/--compress` to
    gzip the subsampled FASTA once it finishes building the index (the uncompressed file is removed after compression).

    
    ### For metax version <0.9.12:
    ```shell
    build_db <fasta_file> -o <database_dir>
    ```

- Pathogen host map file for pathogen detection mode

You may optionally provide a custom pathogen host mapping file to enable Metax to prioritize detection of microorganisms relevant to a specific host of interest.
The mapping file should be a tab-delimited table containing the following columns (see the format of `data/pathogen_host_disease.txt`):
txid (microbial taxon ID), host_txids (associated host taxon IDs), host (host name or label) and diseases (asociated diseases). 

The diseases column can be left blank.

For convenience, we also provide a precompiled virus host mapping file (`data/pathogen_host_disease.txt`) generated from the Virus-Host Database.

3. Test data
   - [CAMI II marine](https://frl.publisso.de/data/frl:6425521/marine/short_read/)
   - [CAMI II pathogen detecton](https://frl.publisso.de/data/frl:6425521/patmgCAMI2.tar.gz)


## How to run


### For version >=0.9.12

```shell
Usage: metax profile [OPTIONS] --outprefix <PREFIX> [-- <EXTRA_ARGS>...]

Arguments:
  [EXTRA_ARGS]...  Additional arguments passed directly to maCMD (use after `--`).

Options:
      --db <DB>                   Path to the maCMD reference database (metax_db.json).
      --dmp-dir <DMP_DIR>         Directory containing the NCBI-style taxonomy dump (dmp files).
  -i, --in-seq <READS>            Comma-separated list of input read files (one or two for Illumina paired-end).
  -o, --outprefix <PREFIX>        Prefix for output files.
  -t, --threads <THREADS>         Number of threads to use for alignment and profiling. [default: 20]
  -r, --resume                    Resume profiling by reusing existing alignment output if present.
      --reuse-sam <SAM>           Existing SAM (or compressed SAM) file to reuse instead of running maCMD.
      --sequencer <TYPE>          Sequencer type (e.g. Nanopore, PacBio, Illumina). [default: Illumina]
  -p, --is-paired                 Treat Illumina inputs as paired-end reads (expects two files).
      --strain                    Enable strain-level profiling outputs.
      --mode <MODE>               Alignment mode preset: default, recall, or precision. [default: default] [possible values: recall, precision, default]
      --batch-size <N>            Maximum number of reads to process per batch.
      --identity <FLOAT>          Minimum alignment identity threshold for retaining a read.
  -m, --mapped-len <LEN>          Minimum mapped read length threshold.
  -b, --breadth <FRACTION>        Minimum breadth of coverage required to report a genome.
      --chunk-breadth <FRACTION>  Manually set the minimum chunk breadth (overrides automatic estimate).
  -f, --fraction <FRACTION>       Minimum aligned fraction a read must cover to be considered.
  -l, --lowbiomass                Apply heuristics tuned for low biomass samples.
  -k, --keep-raw                  Retain the unfiltered rprofile.txt output alongside the final profile.
      --by-aligned                Estimate the minimum chunk breadth using the number of aligned reads.
  -z, --compress-sam              Compress the generated SAM file after profiling completes.
      --pathogen-host <TSV>       Optional TSV mapping pathogen taxids to host metadata for annotation.
      --host <TAXID>              NCBI taxid of the host organism (enables pathogen-specific profiling).
      --verbose                   Log the full command line parameters.
  -h, --help                      Print help
  -V, --version                   Print version
```

```shell
metax --dmp-dir <dump_dir> \
    --db <reference_db> \
    -i <r1>[,<r2>] \
    -o <output_prefix> \
    [other options ...]
```
`<dump_dir>`: path to the `metax_dmp` folder where the dump files located.
`<reference_db>`: path to the json file of the database (e.g. `metax_bavfph/metax_db.json`)
The first run (sample) takes a bit longer; subsequent runs will be substantially faster by using the cached database.

  
### For version <0.9.12

```shell
 Usage: metax [OPTIONS] [EXTRA_ARGS]...

 A taxonomy profiler for metagenomic data

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│    --db                 PATH                                 The reference database file.                          │
│    --dmp_dir            PATH                                 The directory of dmp files.                           │
│    --in_seq         -i  TEXT                                 The input read files separated with comma.            │
│ *  --outprefix      -o  TEXT                                 The prefix of output files. [required]                │
│    --threads        -t  INTEGER                              Number of threads to use.                             │
│    --resume         -r                                       Resume from the last run.                             │
│    --reuse_sam          PATH                                 The sam file to reuse for profiling.                  │
│    --sequencer          [Illumina|Nanopore|PacBio|assembly]  Sequencer used to generate the reads. Default:        │
│                                                              Illumina                                              │
│    --is_paired      -p                                       Whether the reads are paired or not?                  │
│    --strain                                                  Whether profile on strain level? (experimental)       │
│    --mode               [recall|precision|default]           The mode of the profiler. recall: ensure high recall, │
│                                                              precision: ensure high precision, default: use the    │
│                                                              default mode.                                         │
│    --batch_size         INTEGER                              Reduce memory consumption with smaller batch size.    │
│                                                              (Default: 5000 for short, 1000 for long reads)        │
│    --identity           FLOAT                                The sequence identity (matched bases/gap compressed   │
│                                                              len) cutoff to consider a valid assignment. (Default: │
│                                                              0.95 for short, 0.86 for long reads)                  │
│    --mapped_len     -m  INTEGER                              The mapped length cutoff to consider a valid          │
│                                                              assignment. (Default: 50 for short, 250 for long      │
│                                                              reads)                                                │
│    --breadth        -b  FLOAT                                The genome breadth coverage cutoff to consider the    │
│                                                              presence of a genome.                                 │
│    --chunk_breadth      FLOAT                                The genome chunk breadth coverage cutoff to consider  │
│                                                              the presence of a genome.                             │
│    --fraction       -f  FLOAT                                The fraction of matched based in a read to consider a │
│                                                              valid alignment. (Default: 0.6)                       │
│    --lowbiomass     -l                                       Is a low biomass sample? (No coverage filter by       │
│                                                              default for low biomass sample)                       │
│    --keep_raw       -k                                       Keep raw profiling file without statistical           │
│                                                              filtering.                                            │
│    --pathogen_host      PATH                                 The pathogen host table file                          │
│    --host               TEXT                                 The host taxid for pathogen detection                 │
│    --version                                                 Show the version and exit.                            │
│    --help           -h                                       Show this message and exit.                           │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

```

```shell
metax --dmp_dir <dump_dir> \
    --db <reference_db> \
    -i <r1>[,<r2>] \
    -o <output_prefix> \
    [other options ...]
```

`<dump_dir>`: path to the `metax_dmp` folder where the dump files located.
`<reference_db>`: path to the json file of the database (e.g. `metax_bavfph/metax_db.json`)


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
column 2: Name of the most likely taxon
column 3: taxonomy ID of the most likely taxon
column 4: Rank of the taxon
column 5: Names of all possible taxa that the reads originated from
column 6: Taxonomy IDs of all possible taxa
column 7: Likelihood for each of those possible taxa
```

column 2 and column 5 are not included in the version >=9.12


## FAQ

1. What platforms and operating systems does Metax support?

Metax currently supports Linux on x86-64 (64-bit Intel/AMD) systems. Other architectures (e.g., ARM/macOS) are not yet officially supported.

2. Why do I get the error: “Processor 6174 is not supported by this build”?

This error indicates that your CPU does not support some modern instruction sets required by Metax.


## Acknowledgement

We thank Gary Robertson for IT support, Dr. Mohammad-Hadi Foroughmand-Araabi for advice on statistical formulations, and Hesham Almessady for software testing.




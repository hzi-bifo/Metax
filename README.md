# Metax

This repo contains the binary release of Metax, a command-line tool for taxonomic profiling of metagenomic sequences.

## Installation
- Install the package using conda:
   ```shell
   conda install -c zldeng metax
   ```
<!-- - Alternatively, you can install it manually
   
   Download the pre-built binary from the [releases page](https://github.com/dawnmy/Metax/releases), unpack it, and add the directory to your PATH.
   
   Install the dependencies:
    ```shell
    conda install -c bioconda ma=1.1.4
    ``` -->

## Download databases:

- Taxonmy dmp files
  1. Create a `Metax/data/` directory.  
  2. Download the NCBI taxonomy dump (`taxdump.tar.gz`) from:  
     [here]([placeholder](https://research.bifo.helmholtz-hzi.de/downloads/cami/taxdump.tar.gz)) or the latest from [NCBI](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz).
  3. Extract its contents directly into `Metax/data/`.  
  4. _Optional:_ To use an alternative taxonomy source (e.g. GTDB or ICTV), replace the extracted `taxdump` files in `Metax/data/` with your own dmp files.
  
- Reference database

A pre-built reference database is available at [here](https://research.bifo.helmholtz-hzi.de/downloads/cami/metax_db.tar.xz). It is based on the RefSeq snapshot of 10 August 2022 and includes top genomes for each NCBI taxonomic identifier (txid), prioritizing assemblies flagged as “representative” or “reference” and then selecting the highest assembly level (Complete Genome > Chromosome > Scaffold > Contig). In total, it contains 33,143 genomes from bacteria, archaea, viruses, fungi, protozoa, and Homo sapiens.

A customized reference database can be created by following steps:

1. Prepare the genomes in fasta format, the header of each sequence should be in the format:
    ```
    >genome_id|txid|species_txid|sequence_id[|genome_size]
    ```
    Each genome must have a unique genome_id, and each sequence a unique sequence_id. The genome_size field is optional. When using the NCBI taxonomy, txid should be the genome’s NCBI Taxonomy ID, and species_txid the species’ NCBI Taxonomy ID. If you choose a different taxonomy source (e.g. GTDB, ICTV), use the corresponding IDs from your taxonomy dump files.

2. Run the following command to build the database:
    
    <!-- Download and use the `build_db` tool provided in the `utils/` directory in this repository: -->
    ```shell
    build_db <fasta_file> -o <database_dir>
    ```
    
    It may take long to complete, please run it in a Tmux session or screen.

## How to run

You can get the help message by:

```shell
metax --help
```

You must specify the dump directory and the reference database directory when running Metax:

```shell
metax --dmp_dir <dump_dir> \
    --db <reference_db_dir> \
    -i <r1>[,<r2>] \
    -o <output_prefix> \
    [other options ...]
```

## Output

- Final taxonomy profile: `*.profile.txt` and raw (unfiltered) taxonomy profile: `*.rprofile.txt`
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
column 3: taxonomy ID of the taxon
column 4: Rank of the taxon
column 5: Names of all possible taxa that the reads originated from
column 6: Taxonomy IDs of all possible taxa
column 7: Likelihood for each of those taxa
```

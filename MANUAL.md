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

Install Rust 1.89 or later and a C compiler, then build:

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

Metax needs NCBI-style taxonomy dump files. `metax profile` takes the first of:

1. the directory given with `--dmp-dir`;
2. `taxdmp/` in the database directory, if present: a folder of dump files, for example the archive below unpacked by hand. If this folder lacks a dump file or cannot be read, Metax stops rather than use another taxonomy;
3. `taxdmp.tar.gz` in the database directory, the dump that `metax index` stores with the database (§4.1). Metax reads it on every run without unpacking it, so nothing is written for the taxonomy. If the archive is damaged, Metax stops;
4. the managed NCBI dump in `~/.metax/taxdmp` (section 2.1).

A database built by `metax index` (or given its dump with `metax attach-taxonomy`, §4.6) records the scientific name, rank and lineage of every reference taxid in the dump it was built with (`metax_taxonomy.tsv.gz`, §4.1). `metax profile` compares the dump it uses with this record; a dump other than the database's own is compared before mapping starts. When the dump merges database taxids into others, or gives them another name, rank or lineage, as a later NCBI dump can, Metax warns with the number of such taxa and profiles with the names and lineages of that dump. When database taxids are missing from the dump, Metax stops and lists up to ten of them: such a dump is usually not the taxonomy the database was built for, for example an NCBI dump with a GTDB database. Pass `--skip-missing-taxa` to continue anyway: the references whose taxonomy the dump cannot resolve are then ignored, as if they were not in the database, so reads that map only to them stay unassigned and reads of these organisms may be assigned to relatives. A database built without a record is not compared with the dump; a reference taxid missing from the dump then stops the run after mapping, unless `--skip-missing-taxa` is given.

### 2.1 Managed taxonomy

When neither `--dmp-dir` nor a dump stored with the database is available, Metax uses `~/.metax/taxdmp`. On the first run, or when that copy is missing or invalid, Metax downloads the official NCBI `taxdump.tar.gz`, validates it, and keeps `nodes.dmp`, `names.dmp`, `merged.dmp`, and `delnodes.dmp`. Concurrent Metax processes share an installation lock, so a partially written dump is never used.

### 2.2 Explicit taxonomy directory

```bash
mkdir -p taxonomy
curl -L -o taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz -C taxonomy
metax profile --dmp-dir taxonomy ...
```

An explicit directory, like a dump stored with a database, is never downloaded, repaired, or modified by Metax. Use it for offline installations, pinned taxonomy snapshots, or alternative taxonomies (e.g. GTDB or ICTV) converted to the same dump-file format. For a database that stores its own dump, `--dmp-dir` replaces that dump and is checked against the database's record before mapping (§2). The directory must contain `nodes.dmp`, `names.dmp`, `merged.dmp`, and `delnodes.dmp`; each may also be gzipped (`nodes.dmp.gz` and so on).

## 3. Reference Databases

### 3.1 Pre-built databases

Three databases are available, each as one archive:

| Database | Archive | Size | Index cache |
| --- | --- | --- | --- |
| RefSeq | [metax2.1_refseq_202607_k16w37.tar.gz](https://s3.bi.denbi.de/swift/v1/cami3/metax/metax2.1_refseq_202607_k16w37.tar.gz) | 40 GiB | about 42 GiB |
| GTDB marker genes | [metax2.1_gtdb_r232_reps_markers_k16w37.tar.gz](https://s3.bi.denbi.de/swift/v1/cami3/metax/metax2.1_gtdb_r232_reps_markers_k16w37.tar.gz) | 12 GiB | about 15 GiB |
| GTDB fractional | [metax2.1_gtdb_r232_reps_f010_k16w37.tar.gz](https://s3.bi.denbi.de/swift/v1/cami3/metax/metax2.1_gtdb_r232_reps_f010_k16w37.tar.gz) | 32 GiB | about 34 GiB |

- **RefSeq** is based on the July 2026 NCBI RefSeq snapshot and covers bacteria, archaea, viruses, fungi, protozoa and the human reference genome (GRCh38), with up to three genomes per species for bacteria and archaea. Unnamed "sp." isolates are included only where GTDB (release 232) places them in a species of their own. In total, it contains 29,112 genomes. Taxa are named by the NCBI taxonomy.
- **GTDB marker genes** holds the marker genes that GTDB builds its trees from, of the 199,923 species representatives of GTDB release 232: the bac120 genes of 189,801 bacteria and the ar53 genes of 10,122 archaea, 20,078,456 genes of 22.7 Gbp in total, one reference per gene (§4.3). Read counts are the reads on the marker genes, and depth is that of the whole genome. Species and the levels above them are reported (§7.8), but no genome rows (§7.9). Every genome declares its whole-genome size, so `--estimate-unclassified` is available (§7.12).
- **GTDB fractional** is a fractional index (§4.2) of the same 199,923 representative genomes: genomes of at least 500 kb are sampled to a tenth in segments of 50 kb, and smaller genomes are indexed whole, 69.1 of their 654.7 Gbp in total. Read counts are full-index equivalents. Reads not placed at a species are not reported above species, the identity rules of §7.10 apply, and per-read classification is not written.

Both GTDB databases name taxa by the GTDB taxonomy of release 232, converted to NCBI-style dump files with TaxonKit, so their taxids are not NCBI taxids. Each species is represented by one genome, its GTDB representative.

To use a database, download its archive, unpack it, and pass the unpacked folder to `-d`:

```bash
curl -LO https://s3.bi.denbi.de/swift/v1/cami3/metax/metax2.1_gtdb_r232_reps_markers_k16w37.tar.gz
tar -xzf metax2.1_gtdb_r232_reps_markers_k16w37.tar.gz
(cd metax2.1_gtdb_r232_reps_markers_k16w37 && md5sum -c md5sums.txt)
metax profile -d metax2.1_gtdb_r232_reps_markers_k16w37 \
  -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
  -o sample_pe -t 32
```

Each archive unpacks into one folder of the same name that holds `metax_db.mxi`, `metax_db.tsv.gz`, `metax_taxonomy.tsv.gz`, the taxonomy dump `taxdmp.tar.gz` that the database was built with, and `md5sums.txt` with the checksums of these four files, which `md5sum -c` checks. Unpacking needs about as much space again as the archive, and the archive is not needed afterwards. Profiling reads the dump in place unless `--dmp-dir` is given, and checks the dump it uses against the database (§2). The first run writes the index cache into the unpacked folder, which must therefore be writable; to keep the folder read-only, place the cache elsewhere with the mapper option `--index-cache PATH`, given after `--` (§5.4).

Databases built by earlier Metax 2 releases are read as they are, and the log warns that their index format is older. Build such a database again with `metax index` for the current format, which keeps the same seeds in identical sequence of different genomes and so tells close genomes and strains apart more reliably. Databases of Metax 0.9 and earlier cannot be read; build them again with `metax index`.

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
- The label must lead to a taxon that the dump ranks as species, genus, family, order, class, phylum, kingdom or superkingdom, at the label or among its ancestors. A taxon without one, such as a taxon ranked `domain` directly below the root, is accepted by `metax index` but stops `metax profile` after mapping, unless `--skip-missing-taxa` skips its references.
- The sequences of one genome do not need to be adjacent in the FASTA, but all of them must carry the same taxonomy and, when given, the same `genome_size`.
- Malformed headers, unsupported prefixes, duplicate contigs within a genome, and inconsistent taxonomies or genome sizes are rejected with an explanatory error.
- If a provided `genome_size` is larger than the summed contig lengths, indexing warns that only part of the genome is represented. If it is smaller, indexing fails.

Metax keeps the most specific supported `(rank, taxid)` of every genome as its reporting label. A genome labelled `t:1234;s:9606` is reported at species 9606 and, with `--strain`, can get a genome row with strain taxid 1234 (§7.9). A genome labelled `g:561` is reported at taxon 561, with the rank that the taxonomy dump gives it (§7.8).

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
- `metax_reference.fasta.gz`: only for fractional indexes, the exact sampled FASTA used to build the index (not required for profiling);
- `taxdmp.tar.gz`: the taxonomy dump the database was built with, from `--dmp-dir` or, without it, from the managed NCBI dump (§2.1), as one archive. `metax profile` uses it when `--dmp-dir` is omitted;
- `metax_taxonomy.tsv.gz`: a hash of the scientific name, rank and lineage of every reference taxid in that dump, against which `metax profile` checks the dump it uses (§2).

**`metax_db.tsv.gz`** opens with lines starting with `#`: `#metax_db_version`, the layout of the table (below); `#mapping_index`, the index file; `#fractional`, `true` or `false`; `#reference_type`, `markers` or `mixed`, in layouts 3 and 4 only; and `#columns`, the column names. Then it has one row per indexed sequence, which is a contig, a marker gene or, in a fractional index, a sampled segment:

| Column | Meaning |
| --- | --- |
| `reference_name` | name of the sequence; the mapping file names references by it (§6.4) |
| `reference_length` | its length in bases |
| `full_genome_size` | size of its genome; for a marker genome, the summed length of its genes |
| `effective_indexed_size` | the indexed part of the genome, smaller than `full_genome_size` only in a fractional index |
| `genome` | genome identifier |
| `contig` | sequence identifier; in a fractional index, the contig and the coordinates of the segment |
| `taxonomy` | the taxonomy labels of the FASTA header (§3.2) |
| `reference_type` | `genome` or `markers`; layouts 3 and 4 only |
| `declared_genome_size` | whole-genome size that the FASTA header of a marker genome declares, `NA` without one; layout 4 only |

Layout 2 is written for databases of genomes only, so that earlier Metax 2 releases can read them; layout 3 when some genome is of type `markers`, and layout 4 when some marker genome declares its whole-genome size. The database fingerprint in the output files (§6) is a hash of these rows as read, so recompressing the file does not change it.

**`metax_taxonomy.tsv.gz`** opens with `#metax_taxonomy_version`; then, for a dump stored with the database, `#metax_version` and `#created`, the Metax version and UTC time that recorded it, and one `#dump_file` line per dump file with its size in bytes and SHA-256 checksum; then `#columns`. It has one row per taxid that the references name, their strain and reporting taxids, with `lineage_hash`, a 64-bit hash of the taxid, rank and name of every taxon in its lineage. It holds no names: these come from the dump, and the record only tells whether a dump gives these taxa the same names, ranks and lineages (§2).

**`metax_db.mxi`** is binary. Its header holds the index format, 9 for indexes built by this release, the seed settings (§4.4) and the numbers of references, bases and seeds; `metax info -i` prints them (§4.5).

**`metax_db.mxi.metamcache`**, the index cache, is written by `metax profile` beside the index on first use. It records the size and modification time of the index it was made from, and it is made again when these or its format no longer match.

The strain taxid and the reporting taxid of every reference must be in the taxonomy dump. `metax index` checks them before it builds the index and stops with the missing taxids otherwise, so index with the dump the FASTA headers were written for, for example a GTDB dump for GTDB taxids:

```bash
metax index gtdb_genomes.fa -o db_gtdb --dmp-dir gtdb_taxdump -t 20
```

A built database directory is complete: to publish it, share the directory as it is, without the index cache `metax_db.mxi.metamcache` that profiling writes beside the index on first use.

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
- Larger genomes are sampled reproducibly using genome-specific seeds, with segments distributed across contigs approximately uniformly. A genome whose contigs cannot supply the target from segments of `--segment-length` bases takes its longest contigs instead, from their start.
- The full genome size and the effective indexed size are stored separately.

During profiling with a fractional index, coverage and depth use the effective indexed sequence, and assigned-read support is scaled per genome by `full_size / effective_size` before taxonomic aggregation. Fully indexed genomes keep a factor of 1. Per-read classification output is unavailable for fractional databases.

### 4.3 Marker-gene and mixed databases

A genome can be represented by its assembled sequences (reference type `genome`, the default) or by a set of marker genes with one FASTA record per gene (reference type `markers`), such as the single-copy genes that GTDB uses for its trees. Set the type of every genome with `--reference-type` (default `genome`), or of single genomes with a two-column table (`--reference-type-table`):

```bash
metax index marker_genes.fa -o db_markers --reference-type markers -t 20
metax index genomes_and_markers.fa -o db_mixed --reference-type-table types.tsv --dmp-dir taxonomy -t 20
```

`types.tsv` has one `genome_id<TAB>genome|markers` row per genome whose type differs from `--reference-type`; lines starting with `#` are ignored.

Marker genomes:

- are indexed whole, also with `--fraction`;
- have the summed marker length as genome size: read counts (profile column 4) are the reads mapped to the marker genes, not full-genome equivalents, and coverage refers to the marker genes. Depth (profile column 5) is that of the whole genome: a fragment on a marker gene adds the length of its read over the placements that the genome's genes offer (§7.11). Marker genes are largely single-copy, so relative abundances of marker and genome taxa are comparable even though their read counts are not;
- keep a `genome_size` from the FASTA header that is larger than their summed marker length as the size of the whole genome. It is used only for the unclassified estimate (§7.12) and never scales read counts; `metax info` reports how many marker genomes have one;
- have as chunk breadth the share of their genes that reads hit (§7.4);
- take part in the levels above species, where a row whose fragments fit a marker genome best is judged by the genes they hit (§7.8).

A read counts on a marker gene only when enough of it lies inside the gene; a long read also counts on a gene when it aligns to more than half of the gene (§7.4).

In a mixed database, reads from the part of a marker-only genome that is not indexed can map to the genomes of related species of type `genome`. `metax index` lists up to ten genera that hold species of both types; give the species of one genus the same type where possible.

### 4.4 Index options

Seed scheme:

- `-k, --seed-length`: length of each syncmer component. Default `16`.
- `-w, --seed-window`: seed density/window parameter; larger values produce fewer seeds and a smaller index. Default `37`.
- `--max-occ`: discard seed keys occurring at more than N reference positions. Default `500` for an index of marker genes only (`--reference-type markers` without `--reference-type-table` or `--fraction`), where conserved genes share many seeds, and `100` otherwise.
- `--seed-mode`: `randstrobe` (default) or `minimizer`.
- `--hash-bits`: high hash bits stored as the lookup key (24–64). Default `32`.
- `--strobe-start-dist-min` / `--strobe-start-dist-max`: start-to-start distance range between the two strobes. Defaults `30` and `90`.
- `--first-strobe-only`: compatibility mode that matches only the first strobe; more sensitive but less specific than the default full-pair verification.

Resources:

- `-t, --threads`: worker threads for sampling and parallel seed generation. Default `1`.
- `--seed-chunk-mib`: approximate reference bases held per parallel seeding batch. Default `1024`. A single record larger than this is kept intact.
- `--temp-dir`: existing directory for large temporary seed partitions and sort stores. Fast local scratch is recommended. Without it, temporaries are created inside the hidden staging area beside `--output`. They are removed after indexing.

Progress:

- `-q, --quiet` suppresses progress and summary messages; `--progress` forces live progress output; `--no-progress` disables it; `--progress-interval N` sets the reference records between count-based updates. Default `1000`.

### 4.5 Inspect a database

```bash
metax info -d db
metax info -i db/metax_db.mxi
```

`-d, --db` takes the database directory or its metadata file, and `-i, --index` a mapping index. `info -d` summarizes the database metadata, including its reference type, number of marker genomes, whether it stores its taxonomy dump and how many taxa its taxonomy record holds, and, when the `.mxi` file is present, the mapping-index header. `info -i` inspects an index directly. Neither command loads the full index into memory.

### 4.6 Add the taxonomy dump to an existing database

A database built without a taxonomy dump can be given one without rebuilding the index:

```bash
metax attach-taxonomy -d db --dmp-dir taxonomy
```

Pass the dump the database was built with. Metax checks every reference taxid against it, as `metax index` does, and then adds `taxdmp.tar.gz` and `metax_taxonomy.tsv.gz` to the database (§4.1). Nothing is added when a taxid is missing, and files already in the database are never replaced.

## 5. Run Profiling

### 5.1 Read types

```bash
# Illumina paired-end reads: two read files mean paired-end reads
metax profile -d db -i sample_R1.fastq.gz -I sample_R2.fastq.gz -o sample_pe -t 32

# Illumina single-end reads
metax profile -d db -i sample.fastq.gz -o sample_se -t 32

# Nanopore reads
metax profile -d db -i ont_reads.fastq.gz -s Nanopore -o sample_ont -t 32

# PacBio reads
metax profile -d db -i pacbio_reads.fastq.gz -s PacBio -o sample_pacbio -t 32
```

A *fragment* is a read, or a read pair counted once.

### 5.2 Profile modes and read-level gates

`--mode` selects an evidence preset: `default` (balanced), `recall` (more permissive), or `precision` (stricter). A mapping can support a species when it passes three read-level gates:

| Reads | Gate | Recall | Default | Precision |
| --- | --- | --- | --- | --- |
| Illumina | identity | 0.96 | 0.97 | 0.98 |
| Illumina | mapped bases **and** mapped fraction | 50; 0.60 | 50; 0.70 | 50; 0.80 |
| Nanopore, PacBio | identity | 0.90 | 0.92 | 0.94 |
| Nanopore, PacBio | mapped bases **or** mapped fraction | 300; 0.20 | 500; 0.20 | 500; 0.30 |

Set any of them with `--min-identity`, `--mapped-len` and `--min-mapped-fraction`:

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample -t 32 --min-identity 0.95 --mapped-len 100 --min-mapped-fraction 0.65
```

A pair is judged by its identity over the bases of both mates, and each mate must pass the length and fraction gates. The gates decide whether a mapping can support a species. They do not discard mapper output: a read below them can still count above species (§7.8).

**Identity cutoff from the reads' quality.** Metax samples the first 10,000 reads of each input file (at most 250 million bases) and takes their mean base accuracy from the Phred scores, averaged over both files for pairs. Long reads are sampled in every run, because the accuracy decides whether they count as noisy; short reads are sampled only without `--min-identity`. Without `--min-identity` the accuracy sets the identity cutoff. The room left for divergence from the reference is `0.03` (default), `0.04` (recall) or `0.02` (precision):

- Illumina: the cutoff is `accuracy − room`, never below `0.95`.
- Nanopore and PacBio with an accuracy below `0.95` (*noisy long reads*): the cutoff is `accuracy − room`, never below `0.84`. Also with `--min-identity`, the mapper then looks up a seed by its first strobe alone (§5.4), and reads below the species gates can count for a kept species (§7.8).
- Nanopore and PacBio with an accuracy of at least `0.95`: the preset cutoff.
- FASTA input has no qualities and keeps the preset cutoff.

The accuracy is stored in a saved mapping: `--reuse-mapping` in the same mode applies the same cutoff, and in another mode applies that mode's room to the stored accuracy.

**Qualities that overstate accuracy.** For an automatic long-read cutoff, Metax takes each read's best alignment among those passing the length or fraction gate, and the upper quartile of their identities. When this quartile lies below the sampled accuracy, it replaces the accuracy for the identity cutoff; whether the reads count as noisy still follows the sampled accuracy. The cutoff is only ever lowered.

### 5.3 Profile options

Options before `--` configure the profiler. Rows without a default have none or are resolved as described.

| Option | Default | Meaning |
| --- | --- | --- |
| `-d, --db` | | Database directory, or its `metax_db.tsv.gz`. |
| `--dmp-dir` | | Taxonomy dump directory (§2). |
| `--skip-missing-taxa` | off | Skip references whose taxa the dump lacks (§2). |
| `-i, --reads` | | Single-end or long-read file, or mate 1. |
| `-I, --reads2` | | Mate 2; gives paired-end mode. |
| `-o, --output` | | Output prefix. |
| `-s, --sequencer` | `Illumina` | `Illumina`, `Nanopore` or `PacBio`. |
| `--mode` | `default` | `default`, `recall` or `precision` (§5.2). |
| `-S, --strain` | off | Report genomes within each kept species (§7.9). |
| `-l, --low-biomass` | off | Shorthand for `--min-chunk-breadth 0` (§7.6). |
| `-t, --threads` | `20` | Worker threads for mapping and profiling. |
| `--min-identity` | §5.2 | Identity gate for species evidence. |
| `-m, --mapped-len` | §5.2 | Mapped-bases gate for species evidence. |
| `-f, --min-mapped-fraction` | §5.2 | Mapped-fraction gate for species evidence. |
| `--strict-coverage` | off | Report only species kept by the species tests (§7.8). |
| `--estimate-unclassified` | off | Estimate the share of the input reads that the profile does not explain, written on its first line (§6.1, §7.12). Not with `--host`. |
| `--identity-margin` | §7.10 | Identity margin on fractional databases; without a value `0.025`. |
| `--no-identity-margin` | off | Turn the identity margin off. |
| `--genus-identity` | `0.8` | Identity a read without a species-level match needs to count above species (§7.8). |
| `--low-map-rate-threshold` | `0.3` | Mapping rate of the low-map rule (§7.6). |
| `-b, --min-breadth` | | Drop a taxon below this breadth (§7.11). |
| `-c, --min-chunk-breadth` | automatic | Chunk-breadth cutoff of the amount test (§7.6). |
| `--chunk-breadth-scale` | `1.0` | Multiply the automatic cutoff (§7.6). |
| `--by-mapped` | off | Read the automatic cutoff at the mapped reads (§7.6). |
| `-n, --min-reads` | | Drop a taxon with fewer reads (§7.11). |
| `--min-oebr` | §7.5 | Lower OEBR bound, for species and genomes. |
| `--min-coebr` | §7.5 | Lower COEBR bound, for species and genomes. |
| `--min-marker-genes` | `3` | Marker genes a marker genome must hit (§7.5). |
| `--relative-identity-gap` | §7.10 | Identity gap to a deeper species of the genus. |
| `--no-relative-identity-gap` | off | Turn the relative identity gap off. |
| `--host` | | Host taxid (§5.6). |
| `--pathogen-table` | | Pathogen table for `--host` (§5.6). |
| `-r, --resume` | off | Reuse `<prefix>.mapping.tsv.gz` when present (§5.5). |
| `--reuse-mapping` | | Profile a saved mapping (§5.5). |
| `--save-mapping` | off | Write `<prefix>.mapping.tsv.gz`. |
| `-k, --save-decisions` | off | Write `<prefix>.decisions.tsv` (§6.2). |
| `--save-classification` | off | Write `<prefix>.classify.txt.gz` (§6.3). |
| `--stream-mappings` | off | Pass mappings to the profiler without a temporary file; can raise peak memory. Not with `--save-mapping` or `--temp-dir`. |
| `--temp-dir` | automatic | Directory of the temporary mapping file. Without it, Metax prefers suitable local temporary storage and falls back to network temporary storage or the output directory. Not with `--save-mapping`. |
| `--max-em-iterations` | `200` | Limit on the iterations of each EM. |
| `--lenient-mapping-input` | off | Skip malformed rows of a saved mapping with a warning instead of stopping. |
| `-v, --verbose` | off | Show resolved settings and progress. |
| `--very-verbose` | off | Also show the command and the complete mapper settings. |
| `--progress` | off | Force live progress displays. |
| `--no-progress` | off | Disable live progress and periodic mapper messages. |
| `--progress-interval` | `100000` | Records between progress refreshes. |
| `--mapper-help` | | Show every mapper option and exit (§5.4). |

### 5.4 Mapper options

Options before `--` configure the profiler; options after `--` configure the mapper, i.e. which candidate mappings reach the profiler. Show the complete reference with:

```bash
metax profile --mapper-help
```

Resolved mapper defaults:

| Input | top-n | score ratio | min seeds | mapper identity | mapper fraction |
| --- | ---: | ---: | ---: | ---: | ---: |
| Illumina single-end | 5 | 0.90 (0.80 in recall mode; 0.50 on a database with marker genes) | 2 | 0.80 | 0.50 |
| Illumina paired-end | 50 | 0.50 | 1 | 0.80 | 0.50 |
| Nanopore / PacBio | 10 | 0.70 | 3 | 0.80 | 0.30, or 300 mapped bases |

Near-identical marker genes of sister species share most seeds of a short read, and seed sampling can give the gene a read came from fewer seed hits than one of them, so on a database with marker genes an unpaired short read keeps candidates down to half of the best seed score, as a pair does, and alignment identity decides. Illumina reads and noisy long reads (§5.2) find a seed by its first strobe as well; accurate long reads require both strobes of the seed. A short gene holds few seeds, so on a database with marker genes long reads also find a seed by its first strobe, and on a marker gene a chain needs at most two seeds instead of `min seeds`. `--allow-single-strobe-matches` and `--match-both-strobes` override the seed lookup.

The paired-end candidate set is deliberately broad so that the profiler can see competing species and keep conserved regions ambiguous; the modes differ in the profiler's gates and tests (§5.2), not in candidate generation. Example overrides:

```bash
metax profile ... -- --top-n 20 --score-ratio 0.80
```

With `--strain`, `top-n` is at least 50, also when given. Faster paired-end mapping: `--top-n 10` retains at most ten placements per mate at the default score ratio, apart from mates recovered by mate rescue. It is faster, and can lose recall above species level in communities the database covers poorly.

```bash
metax profile ... -- --top-n 10
```

Index cache: the first run against a database writes a memory-mapped cache of the index beside it (`metax_db.mxi.metamcache`), which later runs open directly. `--index-cache PATH` places the cache at `PATH` instead, for example when the database directory is shared or read-only.

Note that `--min-identity` appears on both sides of `--` with different meanings: before it, the profiler's species/strain threshold; after it, the mapper's emission cutoff.

Paired-end emission: `--paired-emission concordant-only` (default) emits same-reference pairs and, for a fragment without any pair, its mates whose mate may lie beyond the end of their reference (§7.1); `concordant-first` emits same-reference pairs when present and otherwise discordant rows; `all` emits every mate mapping. The latter two are mainly useful together with `--save-mapping` for diagnostics. A pair is written once, also when mate rescue finds again a placement that a mate already has.

Mate rescue: a mate that seeding did not place is searched for in the insert window of its partner's placements. `--mate-rescue-taxa N` (default `8`) extends this search to fragments that already have a proper pair: for up to `N` other taxa per mate placed on that mate only, the missing mate is searched for and kept when that pair matches at least as many bases as the best seeded pair (§7.1). `--mate-rescue-taxa 0` restricts mate rescue to fragments without any proper pair.

### 5.5 Saved and reused data

Normal profiling writes the mappings to a temporary compressed file, releases the mapping index, profiles the records, and removes the file.

**Keep the mappings** with `--save-mapping`, which writes `<prefix>.mapping.tsv.gz` (§6.4):

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample -t 20 --save-mapping
```

**Profile saved mappings again** with other profiler settings, without the reads. The same database must be given:

```bash
metax profile --reuse-mapping sample.mapping.tsv.gz -d db -o sample_again -t 20
```

The header of the mapping file supplies the platform, whether the reads are paired, the mapper mode and the read-level gates, so `--reads`, `--reads2`, `--sequencer`, `--resume`, `--save-mapping`, `--stream-mappings`, `--temp-dir` and mapper options after `--` cannot be given. Profiler options override the recorded gates. The recorded gates are used only when the run has the mode of the saved run and, without `--strain`, when the mapping was saved without `--strain`; otherwise the length and fraction gates are the presets of the run, and the identity cutoff is computed from the recorded accuracy with the room of the run's mode, or is the preset when no accuracy is recorded. A mapping saved without `--strain` keeps fewer genomes per read: profiling it with `--strain` gives a warning and can miss genome rows.

**Resume** with `--resume`: Metax profiles `<prefix>.mapping.tsv.gz` if it exists, and otherwise warns and maps the given reads.

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample -t 20 --resume
```

### 5.6 Host filter

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample -t 20 --host 9606 --pathogen-table pathogen_hosts.tsv
```

`--host` and `--pathogen-table` are given together. The table is tab-separated, has a header line, and has four columns; lists are separated by semicolons, and rows with fewer than four columns are skipped:

```text
pathogen_taxid	host_taxids	host_names	diseases
```

Each profile row is judged by its own taxid:

- A taxon the table lists for the given host, or for a host below it in the taxonomy, is kept and annotated with the table's values; an empty or `unknown` host name or disease field is written as `NA`.
- A taxon the table lists for other hosts only is dropped. Its reads go to no other taxon.
- A taxon the table does not list, or lists with an empty or `unknown` `host_taxids` field, is kept and annotated `unknown`.
- A strain row whose taxid is unlisted in this sense follows its species.

The filter runs before abundance is normalised. The profile and decisions files take the prefix `<prefix>.pathogen`, and the profile gains three columns (§6.1). `--host` turns off the low-map rule (§7.6).

### 5.7 Removed options

These options stop the run with a message that says what to use instead.

| Option | Use instead |
| --- | --- |
| `-p, --is-paired` | Nothing: two read files (`-i R1 -I R2`) mean paired reads, and a saved mapping records it. |
| `--exact-ties`, `--pe-exact-ties` | Nothing: the EM shares reads that align equally well to several species. |
| `-g, --genus-fallback` | Nothing: reads no kept species takes are placed above species (§7.8); `--genus-identity` sets their identity minimum and `--strict-coverage` turns this off. |
| `--strict` | Nothing: the species tests use no p-value (§7.5). |
| `--batch-size` | The mapper option `--batch-size` after `--`. |
| `--save-raw-profile`, `--keep-raw` | `--save-decisions`. |
| `--reuse-sam` | `--save-mapping`, then `--reuse-mapping`. |

## 6. Output Files

For an output prefix `sample`, Metax writes:

- `sample.profile.txt`: the taxonomic profile;
- `sample.log`: the run log;
- `sample.decisions.tsv`: how every judged taxon was decided, with `--save-decisions` (`-k`);
- `sample.classify.txt.gz`: per-read assignments, with `--save-classification`;
- `sample.mapping.tsv.gz`: the coordinate mappings, with `--save-mapping`.

With `--host`, the profile and decisions files use the prefix `sample.pathogen`.

### 6.1 Profile

The profile is tab-separated. It opens with two lines starting with `#`: the run's metadata and the column names.

The first line starts with `#metax_profile` and holds these fields:

| Field | Meaning |
| --- | --- |
| `version` | Metax version |
| `database_fingerprint` | content fingerprint of the database, as in the mapping file |
| `sequencer` | `Illumina`, `Nanopore` or `PacBio` |
| `paired` | `true` or `false` |
| `host` | the `--host` taxid, only with `--host` |
| `input_reads` | the input reads, a read pair counting once, only with `--estimate-unclassified` |
| `estimated_unclassified_percent` | the estimated share of the input reads that the profile does not explain, in percent, only with `--estimate-unclassified`; `NA` when it cannot be estimated (§7.12) |

Columns, with `NA` where a value does not apply. Numbers have five decimals, and p-values below `1e-4` are in scientific notation:

| Column | Name | Meaning |
| ---: | --- | --- |
| 1 | `name` | taxon name |
| 2 | `taxid` | taxid; a genome row without a strain taxid carries its species' taxid |
| 3 | `rank` | rank of the row |
| 4 | `reads` | fragments, as full-index equivalents (on a fractional database scaled up by the indexed share of the row's genome; for a marker genome, the fragments on its marker genes) |
| 5 | `depth` | depth from abundance estimation, each fragment counted once (§7.11) |
| 6 | `abundance` | relative abundance in percent, over the rows of the profile |
| 7 | `breadth` | observed breadth on the row's genome |
| 8 | `expected_breadth` | breadth expected from the row's fragments |
| 9 | `breadth_pvalue` | p-value of the breadth against its expectation |
| 10 | `fixed_chunk_breadth` | share of the genome within 1,000 bases of a fragment's first read; for a marker genome, the share of genes hit |
| 11 | `flex_chunk_breadth` | breadth in chunks scaled to the genome; for a marker genome, the share of genes hit |
| 12 | `expected_flex_chunk_breadth` | expected value of column 11 |
| 13 | `chunk_pvalue` | p-value of column 11 against column 12 |
| 14 | `evidence_depth` | depth of the coverage evidence on the row's genome, before abundance estimation |
| 15 | `best_fit_genome` | accession of the genome that columns 7 to 14 describe |
| 16 | `mean_identity` | mean identity of the taxon's fragments to its best-fit genome (§7.9 says which fragments count) |
| 17 | `estimated_ani` | estimated average nucleotide identity of the sample's strain to that genome, at most 1: the mean identity relative to the sample's baseline identity |

Every profile has these 17 columns, whatever the reads, the database and the options; only `--host` adds columns. The options and the database change the rows instead: `--strain` adds genome rows of rank `strain` (§7.9), `--estimate-unclassified` adds two fields to the first line, and on a marker-gene or fractional database columns 4, 10 and 11 count as stated above.

What columns 7 to 17 describe depends on the row:

| Row | Columns 7 to 15 | Columns 16 and 17 |
| --- | --- | --- |
| species, or a unit labelled above species (§7.8) | its best-fit genome | its identity to that genome; `NA` when none of its fragments has it as best match (§7.9) |
| genome row (rank `strain`) | the genome itself | its identity to the genome |
| row above species, built from the fragments that no kept unit takes (§7.8), or a species reported at its genus (§7.10) | the genome its evidence is measured on | `NA` |

With `--host`, three columns follow:

| Column | Name | Meaning |
| ---: | --- | --- |
| 18 | `host_names` | host names from the pathogen table, `unknown` for an unlisted taxon |
| 19 | `host_taxids` | host taxids from the pathogen table, `unknown` for an unlisted taxon |
| 20 | `diseases` | diseases from the pathogen table, `unknown` for an unlisted taxon |

Rows are sorted by the depth of their taxon. A species is followed by its genome rows. Abundance is normalised once over all rows, not per rank, and a row above species does not include the species listed separately.

### 6.2 Decisions table

`<prefix>.decisions.tsv` (`--save-decisions`, `-k`) has one row per taxon judged at each rank, kept or not. It opens with a `#metax_decisions` line with the fields of §6.1 from `version` to `host`, then the column names.

Row kinds:

- `species`: one row per measured genome of every candidate species. A unit labelled above species has these rows at its rank.
- `strain` (with `--strain`): one row per genome the genome level considered.
- `genus`: one row per kind of fragment (without a species-level match, or of species not kept) whose fragments unique to the genus meet the amount of their evidence model (at least 200 fragments on a genome reference, at least 10 genes hit on a marker genome, §7.8), and one for species reported at the genus when no such row holds the evidence of the genus row.
- `family`, `order`, `class`, `phylum`, `domain`: one row per taxon whose fragments meet that amount.

Columns, with `NA` where a value does not apply:

| Column | Meaning |
| --- | --- |
| `rank`, `taxid`, `name` | the taxon; a strain row as in the profile |
| `parent_taxid` | the species of a strain, the genus of a species, the family of a genus |
| `genome` | the genome the row's evidence is measured on |
| `decided` | `true` on the row the taxon's decision rests on: for a species, the genome it is reported on when kept, otherwise its genome with most reads; for strains, each genome with a row in the profile; for a genus, the row whose evidence the profile shows, otherwise the one with most unique reads; `true` on rows above genus |
| `reads` | observed fragments behind the evidence, not scaled. A strain row counts all fragments resolved to the genome, and a row for species reported at the genus their summed EM reads |
| `depth`, `breadth`, `expected_breadth` | coverage on `genome` |
| `oebr` | observed over expected breadth |
| `fixed_chunk_breadth`, `flex_chunk_breadth`, `expected_flex_chunk_breadth` | chunk coverage on `genome` |
| `coebr` | observed over expected flexible chunk breadth |
| `chunk_cutoff` | chunk-breadth cutoff of the tests |
| `breadth_pvalue`, `chunk_pvalue` | p-values as in the profile |
| `candidate_reads` | fragments whose read group in the first EM lists the species; for a genus, fragments of the kind whose tie set lists it; above genus, fragments reaching the taxon |
| `won_reads` | fragments resolved to the species; for a genus, those unique to it |
| `second_look_reads` | fragments the species would take on the second look (§7.7), whether or not it was then kept |
| `own_reads` | fragments it won that list no kept species, counted for species that reach this step of the own-evidence look (§7.7), `0` for the others |
| `em_share` | a strain's share of its species' reads in the genome EM |
| `strain_identity_reads` | a strain's resolved fragments at strain-level identity |
| `window_spread` | share of the genome's windows that the fragments of a genus or higher row reach; on a marker genome the windows cut its genes laid end to end, and the share is shown but not tested |
| `genome_tests` | `passed`, each failed test with its value and bound, or why the row was not tested; a genus row names the kind of fragment first |
| `decision` | `kept` when the taxon has rows of its own in the profile, otherwise `removed` |
| `reason` | the outcome in words; a genus or higher taxon gives its reads, their mean identity, genes hit on a marker genome, chunk ratio or window spread against the bounds, the spill bound of kept units beneath it, and where the reads of a failed taxon went |
| `reported_as` | rank of the profile rows holding the taxon's reads, or `none` |
| `final_reads`, `final_abundance` | columns 4 and 6 of the taxon's profile row; for a species, with its genome rows; `NA` for a genus whose reads are in the row of a unit labelled at that genus |

Every decisions table has these 31 columns, whatever the reads, the database and the options. `--strain` adds rows of kind `strain`, and with `--host` the table is named `<prefix>.pathogen.decisions.tsv`.

### 6.3 Per-read classification

`<prefix>.classify.txt.gz` is gzip-compressed and tab-separated. It opens with a `#metax_classify` line with the fields of §6.1 from `version` to `host`, then the column names. It has one row per fragment and five columns: `read_id`; `taxid` and `rank` of its largest share; `taxids`, the taxids of its shares; and `shares`, in the same order, largest first, with five decimals. Several taxids or shares are separated by `;`, and a fragment with a single share has share `1`.

The shares are the ones the profile's rows count the fragment with. Shares on taxa without a profile row are dropped and the others renormalised, and a fragment without any share is not written. For a species represented by marker genes, the file lists only the reads placed on those genes, and the log says so. The file is not written for a fractional database.

### 6.4 Mapping file

`<prefix>.mapping.tsv.gz` is tab-separated. It opens with three header lines, for the run, the database and the column names, then has one row per retained mapping and a final line of sample counts:

```text
#metax_mapping  version=5  sequencer=Illumina  paired=true  mapper_mode=default  strain=false  profiler_identity=0.97  profiler_mapped_len=50  profiler_fraction=0.7  sampled_mean_base_accuracy=NA
#metax_mapping_database  fractional=false  sequences=71442  genomes=33533  fingerprint=0123456789abcdef0123456789abcdef
#read_id  reference_name  start  end  identity  read_length  mapped_bases  mate_reference  mate  pair_id  candidate_boundary
S0R0  GCF_000001405.40|NC_000001.11|s:9606|3298430636  1872848  1872998  0.99333  150  150  =  1  1  0
#metax_mapping_sample  total_reads=33237808  mapped_reads=23098244  emitted_records=85758065
```

Mapping columns, eleven for every kind of read: `read_id`; `reference_name`, as in `metax_db.tsv.gz`; `start` and `end`; `identity`; `read_length`; `mapped_bases`; `mate_reference`, where the mate maps (`=` on the same reference, `>` possibly beyond the end of this reference, `*` unmapped, `NA` for a single read, or the mate's reference); `mate`, the mate number (`0`, `1` or `2`); `pair_id`, a pair id local to the read; and `candidate_boundary`, a mask (`0` to `3`) that shows, for single-end and long reads, whether further endpoints may have been omitted by the candidate cap, `0` for pairs. The profiler does not read the mask. Reference, genome, taxonomy and size tables are not copied into the file; they come from the database, whose fingerprint is checked on reuse. The first line records the run: `version`, the layout of the file; `sequencer` and `paired`; `mapper_mode`, the `--mode` the reads were mapped in (§5.2); `strain`, whether the run kept the extra mappings of `--strain` (§7.9); `profiler_identity`, `profiler_mapped_len` and `profiler_fraction`, the read-level gates of the run (§5.2); and `sampled_mean_base_accuracy`, the mean base accuracy measured from the reads' quality scores, `NA` when it was not measured. `--reuse-mapping` profiles with these settings unless the command line sets others; it warns when `--mode` differs, since the saved mappings cannot be made again, and when `--strain` is given for mappings made without it. The second line describes the database: whether it is fractional, its numbers of sequences and genomes, and its fingerprint, which must equal that of the database given with `--reuse-mapping`. The last line counts the input: `total_reads` and `mapped_reads` count every read as sequenced, both mates of a pair separately, unlike the profile, where a read pair counts once; `emitted_records` is the number of mapping rows.

A file without the run and database lines is rejected; a file without the line of column names, as written by earlier releases, is read as well.

## 7. Method

### 7.1 Mapping and candidate generation

The built-in mapper indexes references with syncmer-based randstrobe seeds and retains up to `top-n` candidate mappings per read or mate whose seed score is at least `score-ratio` of the best. For long reads, when no mapping reaches a high identity, a rescue pass also evaluates candidates below the score ratio (`--no-long-read-rescue` turns it off). Candidates are refined to reference intervals with a per-base identity: an ungapped comparison when read and interval have the same length of at most 512 bases, otherwise a gapped identity computed by exact edit distance. Mapper-level acceptance requires identity of at least `0.80` and the platform's mapped-fraction gate (§5.4).

A long read on a marker gene is also accepted when it aligns to more than half of the gene. A long read can run over several genes of a genome, so when any of its candidates lies on a marker gene, candidates on different parts of the read (read intervals overlapping by less than half of the shorter one, directly or through other candidates) do not compete: the score ratio and `top-n` apply within each part. A candidate on a genome reference joins the parts it overlaps, so a read placed along a genome stays one part. Both strobes of a seed must lie inside a gene, so the seeds of a long read stop short of the gene ends; an alignment on a marker gene whose identity lies close to the cutoff is therefore extended into the gene with a banded x-drop alignment for as long as the read keeps matching, each flank being kept only when its own identity reaches the cutoff, and its identity, mapped bases and reference interval are those of the extended alignment.

Only the candidates with the highest seed scores are evaluated (`top-n` × 4, at least 16 for short reads). Before applying `top-n`, the mapper reserves the strongest accepted mapping for each represented profile endpoint (normally species) and each most-specific endpoint (strain when known), then fills the remaining slots with the next best mappings. Among the evaluated candidates, repeated placements from one taxon therefore do not hide another eligible taxon while slots remain. For single-end and long reads, a candidate-boundary mask in the mapping file records whether endpoints may have been omitted.

Paired-end mapping is pair-aware: mates are scored together, preferring pairs within the insert-size window, and by default only pairs on one reference are emitted, whatever their insert size or orientation. A mate without a seeded placement, or a fragment without any proper pair, triggers mate rescue: the read is searched for in the insert window around its partner's placements. A placement found there is evaluated and accepted like a seeded one, so its identity and mapped bases describe the aligned interval. For a short-read mate without any placement of its own, a longer window in which it has no seed hit at all is scanned without gaps for its best position within the sample's species identity cutoff (§5.2), and that position is evaluated like a seeded one. When no placement is found, the read is aligned as a whole at every position of a window of at most 256 bases.

A mate that remains without a mate on its reference after rescue is kept as a lone mate when the insert window on the side where its mate must lie (downstream of a forward read, upstream of a reverse one) runs past the end of the reference: its mate may come from sequence that the reference does not contain. This applies to marker genes and to every reference of a fractional index, whose ends cut the sequence a fragment came from; in a mixed database each reference follows its own type. On genome references of a full index no mate is kept this way, since the end of an assembled sequence mostly does not cut a fragment. Lone mates are emitted only for fragments without any pair and are marked `>` in mapping files. Other lone mates are dropped.

A proper pair on one taxon does not end the search. Seeding can place a second taxon on one mate only, for example a close relative for which the other mate has too few seed hits, and concordant-only emission would then drop that taxon although it may explain the fragment equally well or better. For each such taxon (profile endpoint, up to `--mate-rescue-taxa` per mate) the missing mate is searched for beside the taxon's best placement. The recovered pair is kept only if it matches at least as many bases as the best pair that seeding placed on one reference, so rescue can add an equally good or better explanation of the fragment but never a worse one. Kept mates are appended after the `top-n` cut and therefore never take a `top-n` slot from a seeded placement.

### 7.2 Profiling workflow

1. **Candidates.** Each fragment gets its candidate units, with their errors (§7.3).
2. **First EM.** An EM over all candidates shares the fragments. A fragment *resolves* to a unit when its posterior there is at least `0.9`, and then adds to the coverage of each genome of that unit it aligns to.
3. **Species tests.** A unit is kept when one of its genomes passes the tests of §7.5.
4. **Second judgment.** Removed units are judged again and some kept units are checked on their own reads (§7.7).
5. **Even coverage.** A removed species that covers its genome evenly but too thinly for the amount test can be kept (§7.7).
6. **Second EM.** It shares the fragments among the kept units.
7. **Rows above species.** The fragments that no kept unit takes are placed at a genus or a higher rank (§7.8).
8. **Genomes.** The fragments of each kept species are shared among its genomes. The genome with the largest share is its best-fit genome, and with `--strain` a genome can get a row of its own (§7.9).
9. **Identity rules.** On a fractional database, species whose reads are far less identical than expected are left out or reported at their genus (§7.10).
10. **Filters and abundance.** `--host`, `--min-reads` and `--min-breadth` drop rows, and abundance is normalised over the rows left (§7.11).
11. **Outputs.** The profile and the optional tables are written (§6).

### 7.3 Fragments, candidates and the EM

A *unit* is a species, or one taxon above species holding the genomes whose references are labelled only at that taxon (§7.8). A fragment's candidates are the units whose genomes it aligns to within the read-level gates (§5.2).

The errors of a candidate are the estimated errors of its best placement, `(1 − identity) × mapped bases` over both mates of a pair, plus `0.5` for each base that this placement aligns less than the fragment's longest alignment. Candidates whose errors exceed the best by less than the tie margin count as tied and carry the best's errors. The tie margin is `0.5` errors for short reads and `0.5 %` of the best candidate's aligned bases for long reads.

The EM shares a fragment among its candidates in proportion to `share(u) × exp(−λ × errors(u))`, with `λ = ln((1 − e) / (e / 3))` for the sample's error rate `e`, measured on the three quarters of the fragments that match their best candidate best. A candidate whose likelihood is below `1e-6` of the fragment's best is dropped. Fragments with the same candidates and the same errors above their best, rounded to `0.5`, form one read group. The EM stops when the largest change of a single unit's weight is at most `1e-7` of the total weight, or at the limit of `--max-em-iterations`.

### 7.4 Coverage measures

Each genome is measured on its own; coordinates from different genomes are never combined. For a genome of size `G` with `n` fragments, with `span_i` the bases of the first read and `read_bases_i` the bases both mates cover, counted once where they overlap:

```text
observed_breadth   = merged bases of all reads, both mates of a pair / G
depth              = sum(span_i) / G
mean_span_fraction = clamp(sum(read_bases_i) / G / n, 0, 1)
expected_breadth   = 1 − (1 − mean_span_fraction)^n
OEBR               = observed_breadth / expected_breadth
```

**OEBR** is the observed over expected breadth: reads piled up in a few places give a value well below 1.

The *fixed chunk breadth* is the share of the genome within 1,000 bases of a fragment's read, the first read of a pair. The *flexible chunk breadth* is the share within `√G` bases of that read's midpoint, clamped to its contig:

```text
expected_flex_chunk = 1 − (1 − 1 / ⌊√G / 2⌋)^n
COEBR               = flex_chunk_breadth / expected_flex_chunk
```

**COEBR** is the observed over expected flexible chunk breadth: it tests whether the fragments spread over the genome as a whole.

For a fractional index, `G` is the indexed size of the genome.

The p-values of the profile test the breadth and the flexible chunk breadth against their expectation. A ratio within `[0.9, 1.2]` gives `1`; otherwise the deviation is tested one-sided with an exact occupancy model or its binomial or normal approximation, multiplied by the number of genomes tested and capped at 1. The species tests do not use them; only the even-coverage rule does (§7.7).

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

where `k_b` is the number of the `s_i` possible read starts on gene `i` whose read covers base `b`, so bases near gene ends are expected to be covered less often. The p-value of the genes hit is `1.0` when COEBR lies within `[0.9, 1.2]`, as for the genome p-values (§7.4); outside it, it is the exact probability of hitting at most `genes_hit` genes (multiplied by the number of tests). The observed breadth of a marker genome merges the first read of each fragment only. The breadth p-value is computed as for genomes. A marker species must hit at least `--min-marker-genes` genes (default `3`).

### 7.5 Species tests

A unit is kept when one of its genomes passes all of:

| Test | Rule | Default |
| --- | --- | --- |
| Amount | fixed chunk breadth above the cutoff | automatic (§7.6) |
| OEBR | between the lower bound and `1.50` | lower bound `0.60` |
| COEBR | between the lower bound and `1.50` | lower bound `0.75` for short reads, `0.70` for long reads |
| Won share | at least this share of the fragments whose read group in the first EM lists the unit resolve to it | `0.30` |
| Marker genes | a marker genome hits at least this many genes | `3` |

`--min-oebr`, `--min-coebr` and `--min-marker-genes` set the bounds. A marker genome none of whose genes can hold a read of the sample's length is measured as one sequence, without the marker-gene test. The OEBR bound allows for amplified libraries, such as those made with PCR-based kits, which cover a genome unevenly base by base while their chunks stay evenly covered.

A kept species is reported on its best-fit genome (§7.9). Its coverage columns describe that genome.

### 7.6 Chunk-breadth cutoff and the low-map rule

Without `--min-chunk-breadth`, the cutoff of the amount test is estimated from the read count. With `m` the basis in millions of reads or fragments:

| Reads | Cutoff |
| --- | --- |
| Illumina | `0.66·m/(m+1) + 0.24·m/(m+100)` |
| Nanopore, PacBio | `0.85·m/(m+0.78)` |

Both are capped at `0.95`. The basis is all input fragments.

**Low-map rule.** A run is a low-map run when its genome-equivalent mapping rate is below `--low-map-rate-threshold` (default `0.30`). Each mapped fragment counts as the sampling scale of the genome of its first mapping record that is not on a marker genome (`1` for a genome indexed whole, `full size / indexed size` for a sampled one), fragments placed only on marker genes are left out, and the sum over the input fragments is capped at 1. On a low-map run the basis is the mapped fragments, and on a database that is not fractional it is 1.5 times their number when no `--min-chunk-breadth` is given. `--host` turns this rule off.

- `--by-mapped` uses the mapped fragments as the basis on every run.
- `--chunk-breadth-scale F` multiplies the automatic cutoff, capped at `1`. It cannot be combined with `--min-chunk-breadth` or `--low-biomass`.
- In `recall` mode, above one million fragments in the basis, the automatic cutoff is also multiplied by `0.70 + 0.30/m` for Illumina (at least `0.70`) and by `0.80` for long reads.
- `-c, --min-chunk-breadth` sets the cutoff directly and is not scaled.
- `-l, --low-biomass` is the shorthand for `--min-chunk-breadth 0`, which removes the amount test's cutoff. It cannot be combined with `--min-chunk-breadth` or `--chunk-breadth-scale`.

### 7.7 Second judgment of the species

**Second look.** Sister species that tie on most fragments can split them so that neither resolves them, and both are removed. Removed units are judged again, largest first-EM weight first, and named units before a unit labelled above species of their taxon. A removed unit takes the read groups on which no other removed unit is more likely and that it holds with a posterior of at least `0.9` against the kept units. It is kept if it then passes the tests, and competes with the units judged after it.

**Own-evidence look.** A removed unit that covers a genome well enough but wins too small a share is examined next. If the fragments it shares with a kept species alone span one of its genomes above the cutoff, it stays removed. Otherwise it is kept when the fragments it won with no kept species among their candidates pass the coverage tests on one of its genomes alone. It is reported at its genus when the fragments where it ties with another species for the best match alone span one of its genomes.

**Prophage check.** A virus integrated in a host genome collects the reads of that region. A kept virus is removed when its own fragments, those that list no kept bacterium or archaeon and no taxon labelled above species in their lineages, do not pass the coverage tests over the rest of its genome. Its reads go to the kept units in the second EM.

**Named units first.** A kept unit labelled above species is judged the same way beside the kept, more specific units of its lineage. A species missing from the database whose genomes are labelled only at the genus covers the rest of its genome and stays; an unnamed copy of a named species is removed, and its reads go to that species.

**Even coverage.** A removed species is kept at species when no taxon of its genus is kept and its fragments on its strongest genome meet all of: at least `200` fragments, OEBR and COEBR within `[0.9, 1.1]`, reads in at least `90 %` of the genome's windows, both p-values at least `1e-5`, and a fixed chunk breadth below the cutoff. The genome is cut into `min(100, n/5)` equal windows for its `n` fragments; with fewer than 4 windows the rule does not apply. This keeps large genomes at low depth. Species represented only by marker genes do not take part. The rule is off with `--strict-coverage` and on a fractional database.

### 7.8 Levels above species

The fragments that no kept unit takes, placed on genome references or on marker genomes, go to the lowest rank their matches and evidence support. `--strict-coverage` turns this off, and these fragments are then dropped.

Two kinds of fragments come here:

- **Without a species-level match.** They fail the species gates (§5.2) but pass the floor gates: identity of at least `--genus-identity` (default `0.80`), a mapped fraction of at least `0.60`, and at least 50 mapped bases on short reads. On long reads they need that fraction or the species length gate.
- **Of species not kept.** They pass the species gates, but none of their candidates is kept.

`--genus-identity` is a per-read minimum for fragments without a species-level match; fragments of species not kept have already passed the species gates. A genus row has its own, stricter tests on the mean of its fragments.

A fragment's *tie set* holds the taxa whose placements lie within the tie margin of its best one.

A candidate row is measured on its *evidence genome*, the genome holding most of its fragments' best placements, with that genome's model (§7.4), and its tests depend on the model:

| Test | Evidence on a genome reference | Evidence on a marker genome |
| --- | --- | --- |
| Amount | at least `200` fragments | at least `10` genes hit (a fragment counts on every gene it covers) |
| Genus row: mean identity | at least the species identity cutoff minus `0.09` | at least the species identity cutoff minus `0.09` |
| Genus row: COEBR | at least `0.40` | at least `0.40` over the genes |
| Row above genus: spread | reads in at least `50 %` of the genome's windows | COEBR of at least `0.40` over the genes |

1. **Genus rows.** A fragment whose tie set lies in one genus is unique to it. The two kinds of fragments unique to a genus are judged separately, never added together: the genus gets a row when its unique fragments without a species-level match, or its unique fragments of species not kept, pass the amount, mean identity and COEBR tests on their own, each on its own evidence genome. Once the genus has a row, its unique fragments of both kinds join it, also those of a kind that failed or was too few to judge.
2. **Ties between genera.** A fragment tied across several genera, at least one with a genus row, is shared by an EM over the genera with a row among all its candidates.
3. **Lowest common ancestor rows.** The other fragments go to the lowest common ancestor of their tie set, or the first taxon above it at family, order, class, phylum or domain rank. A taxon gets a row when the fragments that reach it pass the amount and spread tests. Otherwise they move up one of these ranks, and they stay unexplained above domain.

**Kept species on marker genes.** The marker genes of a kept species miss most of its genome and clip reads at their ends, and those reads land on relatives. A taxon holding kept species represented by marker genes therefore gets a row only when its fragments number at least `0.15` times the fragments of those species; otherwise its fragments move up like those of any taxon without a row. For a genus this counts its unique fragments of both kinds together, unlike the other tests. When the row's evidence genome is a genome reference, each fragment of those species counts as `placements(G) / Σ placements(L_i)` whole-genome fragments for the declared size `G` of the marker genome the species is reported on and its genes of lengths `L_i` (§7.12), or as one fragment when that genome declares no size. A kept unit labelled at a genus or higher taxon bounds the row of that taxon too.

A row's depth is that of its evidence genome; for a genus row, that genome comes from the unique fragments of the kind its row shows, the passing kind with most of them. On a genome reference it is the first-read bases of all its fragments over the indexed length; on a marker genome each fragment adds its read length on every gene it covers over the placements of the genes, as for a species (§7.11). Genomes tied exactly on that number share the depth.

**Noisy long reads.** A fragment without a species-level match whose best match, or one tied with it, is a kept species counts for that species and does not go above species. `--strict-coverage` turns this off.

**Limitation.** The reads of a genus missing from the database land on its closest genera. When a sister genus takes them as unique fragments and passes the tests, the missing genus is reported as that sister genus. This holds for rows on marker genomes as for rows on genome references.

**References labelled above species.** A genome whose most specific label is a genus, family or higher taxon belongs to a unit of that taxon, judged like a species and reported at that rank and taxid. It has no identity rules of §7.10. A row above species of the same taxon is merged with it. When the genomes of a species are only partly labelled at species, and the sample's strain is closest to a genome labelled only at the genus, its reads can be reported at the genus.

### 7.9 Genome level (`--strain`)

The fragments the second EM gives to a species are shared among its genomes by the same EM, with a tie margin of `0.1 %` of the longest placement for long reads. The genome with the largest share is the species' best-fit genome. If that genome was not measured, the species is reported on its best-fit genome of the species level: columns 7 to 15 then describe that genome and columns 16 and 17 the genome EM's. The mean identity is taken over the fragments whose best species match, within the tie margin, is this species, at any identity the mapper reports, and that align to the genome. Relative to the sample's baseline (the median over the species, weighted by fragments), it gives `estimated_ani`.

With `--strain`, the mapper keeps at least 50 mappings per read or mate, so that the genomes of a species can be compared, and a genome gets a row of rank `strain` when:

1. the fragments that the genome EM gives to it (posterior of at least `0.9`) and that match it at strain-level identity, the sample's identity to its own genomes less `0.005`, pass the coverage tests: chunk breadth above the cutoff, OEBR from `0.70` to `1.50`, and COEBR from `0.75` (short reads) or `0.70` (long reads) to `1.50`; and
2. in a second genome EM over the genomes passing step 1, its estimated average nucleotide identity, its fragments' mean identity relative to the median over strain-level fragments, is at least `0.999`. A genome below it is dropped and this EM runs again over the rest.

`--min-oebr` and `--min-coebr` set the bounds of this level too. A marker genome never gets a genome row.

A genome row is named by its strain taxid. Without one, it is named `<species name> <accession>` and carries the species' taxid.

The rows of a species are exclusive: a genome row counts each of its fragments by the species' share of it times the genome EM's posterior, and the species row keeps the rest. **When genome rows exist and less than one read is left, the species has no row of its own**, and its reads are all in its genome rows.

### 7.10 Identity rules on fractional databases

On a fractional database, the reads of a genome's unindexed part land on the indexed segments of its relatives, and such a relative can pass every coverage test. Two rules compare the mean identity of the kept species' fragments to their best-fit genome. Units labelled above species are not judged.

**Identity margin** (`--identity-margin F`). The baseline is the median of the species' mean identities, weighted by their fragments. A species more than `F` below the baseline is left out when a species of its genus within the margin is kept, and otherwise reported at its genus; a species without a genus stays. The rule needs at least 3 species with a measured identity. It applies only when some genome is indexed in part, and has no effect on a full database. The default is `0.025` for Illumina reads, paired or single-end; other inputs use it only when it is given. `--no-identity-margin` turns it off.

**Relative identity gap** (`--relative-identity-gap F`). A species is left out when a species of the same genus with more depth has a mean identity at least `F` higher. This also removes a species that is present but much less abundant than a close relative. It is off by default, except `0.01` for single-end Illumina reads on a fractional database; `--no-relative-identity-gap` turns it off. It can be given for any database.

On a fractional database, fragments not placed at a species are not reported above species, and no species is kept for even coverage: the reads of a genome missing from the database land on the indexed segments of many relatives and cannot be counted as full-index equivalents.

### 7.11 Abundance and report filters

Each fragment adds depth to its candidates: the mapped bases of its read (mate 1 of a pair) over the genome size, which for a fractional index is the effective indexed size. For a taxon kept by the coverage evidence this is the size of its best-fit genome (§7.9), whichever of its genomes the fragment placed best on, so that its abundance comes from the genome its presence is reported on. A fragment on a marker gene adds the length of its read, also when the read is clipped at the end of the gene, over `Σ placements(L_i)` for the genome's genes of lengths `L_i` (§7.4); for a long read these are the placements of its own length rather than of the median, since long reads differ widely in length, and a fragment adds this for every gene of the genome it covers. A genome of `G` bases offers about `G` placements, so this is the depth of the whole genome, comparable with the depth of a genome taxon. When no gene of a marker genome can hold a read of the sample's length, its fragments add mapped bases over the summed gene length.

Abundance is the depth of each row, normalised over the rows left after the filters. Read counts are full-index equivalents.

`--min-reads N` drops a taxon whose rows hold fewer than `N` reads, as full-index equivalents, a species counted with its genome rows. `--min-breadth F` drops a taxon whose breadth is below `F`. The reads of a dropped taxon go to no other taxon.

### 7.12 Unclassified estimate

With `--estimate-unclassified`, Metax estimates how much of the input the profile does not explain. Every row explains the fragments in its column 4, which for a fractional database are already full-index equivalents. A species represented only by marker genomes explains more than its fragments on the marker genes: a genome of declared size `G` offers `placements(G)` fragment placements against `Σ placements(L_i)` for its genes of lengths `L_i`, with the placements per sequence of §7.4. The placements of the genes count a fragment on every gene it covers, so the ratio is multiplied by the gene hits per fragment of the fragments resolved to the genome (above 1 where the two mates of a pair lie on two genes). Its fragments are multiplied by that ratio, averaged over its marker genomes. Long reads differ in length and can cover several genes, so for them the ratio is measured on the reads resolved to the species' marker genomes. Only genes that need all `a` bases from a read (`⌊L/2⌋ + 1 ≥ a`, §7.4) count: a read's hits on them are found reliably, while hits on shorter genes, which need only more than half of the gene, are often missed. Each hit of a read of length `R` on such a gene stands for `placements(G) / Σ placements(L_i)` whole-genome reads at that length, summed over these genes, and the ratio is their sum over the number of reads. A species that also has genome references counts its fragments as they are. Then

```text
unclassified share = 1 − min(1, explained fragments / input fragments)
```

The estimate is taken over the rows of the profile, after every filter; a species counts once with its genome rows. Rows above species (§7.8) explain their fragments, those on a marker evidence genome times that genome's whole-genome ratio, computed as for a species on it; other fragments of removed species, fragments of dropped rows, and fragments that fail the read-level gates count as unexplained. For short reads the estimate is unavailable when any marker genome of a reported marker species has no declared size or no gene that can hold a read of the sample's length; for long reads, when one of its genomes with resolved reads has no declared size or no read is resolved to its marker genomes. It is also unavailable when the evidence genome of a row above species is a marker genome without a whole-genome ratio: it declares no size, no gene can hold a read of the sample's length, or, for long reads, no read placed there hits a gene that needs all `a` bases. The log names the taxon, and `estimated_unclassified_percent` is `NA`. The share is also written to the log. The option cannot be combined with `--host`.

### 7.13 Fixed values

| Value | Setting |
| --- | --- |
| posterior at which a fragment resolves to a unit | `0.9` |
| lowest won share of a kept unit | `0.3` |
| upper bound of OEBR and COEBR | `1.5` |
| lower OEBR bound of a species | `0.6` |
| lower OEBR bound of a genome row | `0.7` |
| lower COEBR bound, short reads | `0.75` |
| lower COEBR bound, long reads | `0.7` |
| tie margin, short reads, in errors | `0.5` |
| tie margin, long reads, share of aligned bases | `0.005` |
| fixed chunk radius in bases | `1000` |
| fewest fragments of an even-coverage species | `200` |
| lowest window spread of an even-coverage species | `0.9` |
| fewest fragments of a row above species on a genome reference | `200` |
| fewest genes hit of a row above species on a marker genome | `10` |
| distance of a genus row's mean identity below the species cutoff | `0.09` |
| lowest COEBR of a genus row | `0.4` |
| lowest window spread of a row above genus on a genome reference | `0.5` |
| lowest COEBR of a row above genus on a marker genome | `0.4` |
| share of the fragments of kept species on marker genes that a row above them needs | `0.15` |
| floor mapped bases of a short read | `50` |
| floor mapped fraction | `0.6` |
| strain-level identity margin | `0.005` |
| lowest estimated identity of a genome row | `0.999` |
| species needed for the identity margin | `3` |
| default identity margin | `0.025` |
| relative identity gap, single-end reads on a fractional database | `0.01` |
| noisy long-read accuracy limit | `0.95` |
| lowest automatic long-read identity | `0.84` |
| lowest automatic short-read identity | `0.95` |

## 8. Common Recipes

### 8.1 High-sensitivity Illumina run

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample_recall -t 32 --mode recall
```

### 8.2 Stricter reporting

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample_strict -t 32 --min-reads 10 --min-oebr 0.85 --min-coebr 0.85
```

### 8.3 Fixed chunk-breadth cutoff and minimum reads

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample_filtered -t 32 --min-reads 10 --min-chunk-breadth 0.20
```

### 8.4 Nanopore run with set thresholds

```bash
metax profile -d db -i ont.fastq.gz -s Nanopore -o ont_custom -t 32 --min-identity 0.92 --mapped-len 200 --min-mapped-fraction 0.55
```

### 8.5 Map once, profile several ways

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample -t 32 --save-mapping
metax profile --reuse-mapping sample.mapping.tsv.gz -d db -o sample_precision --mode precision
metax profile --reuse-mapping sample.mapping.tsv.gz -d db -o sample_min10 --min-reads 10
```

### 8.6 Sample with few reads

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample_lowbiomass -t 32 --low-biomass
```

### 8.7 Marker-gene database

```bash
metax index marker_genes.fa -o db_markers --reference-type markers -t 20
metax profile -d db_markers -i R1.fastq.gz -I R2.fastq.gz -o sample_markers -t 32
```

### 8.8 Genome rows and the decisions table

```bash
metax profile -d db -i R1.fastq.gz -I R2.fastq.gz -o sample_strain -t 32 --strain --save-decisions
```

### 8.9 Fractional database with the relative identity gap

Single-end Illumina reads on a fractional database use the gap by default (§7.10). For paired-end reads, add it for samples in which close relatives are not expected together:

```bash
metax profile -d db_frac_10 -i R1.fastq.gz -I R2.fastq.gz -o sample_frac -t 32 --relative-identity-gap 0.01
```

## 9. Limitations

- The built-in mapper is a coordinate mapper, not a full base-level aligner, and normal runs do not produce SAM. Metax profiles only its own coordinate mappings.
- Per-read classification is not written for fractional databases.
- A reference that consists mostly of foreign DNA can produce a false row above species (§7.8).
- A genus missing from the database can be reported as a close sister genus (§7.8).
- The identity margin also reports at genus, or leaves out, a species present only as a strain that has diverged clearly from every reference of its species (§7.10).
- The unclassified estimate counts reads failing the read-level gates, and reads of removed species that no row above species takes, as unexplained, so it errs high. For marker databases it depends on the genome sizes declared in the reference FASTA (§7.12).
- Windows is supported only through WSL.

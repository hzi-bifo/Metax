#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${TMPDIR:-/tmp}/metax_smoke.$$"

cleanup() {
  rm -rf "$WORK"
}
trap cleanup EXIT

run_metax() {
  if [[ -n "${METAX_BIN:-}" ]]; then
    "$METAX_BIN" "$@"
  else
    cargo run --quiet --manifest-path "$ROOT/Cargo.toml" -- "$@"
  fi
}

make_taxdump() {
  local dir="$1"
  mkdir -p "$dir"
  cat >"$dir/nodes.dmp" <<'EOF'
1	|	1	|	no rank	|
2	|	1	|	superkingdom	|
10	|	2	|	genus	|
20	|	2	|	genus	|
100	|	10	|	species	|
101	|	100	|	no rank	|
200	|	20	|	species	|
300	|	2	|	genus	|
EOF
  cat >"$dir/names.dmp" <<'EOF'
1	|	root	|		|	scientific name	|
2	|	Bacteria	|		|	scientific name	|
10	|	GenusA	|		|	scientific name	|
20	|	GenusB	|		|	scientific name	|
100	|	SpeciesA	|		|	scientific name	|
101	|	StrainA	|		|	scientific name	|
200	|	SpeciesB	|		|	scientific name	|
300	|	GenusOnly	|		|	scientific name	|
EOF
  : >"$dir/merged.dmp"
  : >"$dir/delnodes.dmp"
}

dna_repeat() {
  local n="$1"
  local motif="ACGTGCAATGCCGTAAGTCCGATCGTACGATGCTAGCTAGGCTA"
  local out=""
  while ((${#out} < n)); do
    out="${out}${motif}"
  done
  printf '%s' "${out:0:n}"
}

fastq() {
  local id="$1"
  local seq="$2"
  printf '@%s\n%s\n+\n' "$id" "$seq"
  printf '%*s\n' "${#seq}" '' | tr ' ' 'I'
}

assert_exists() {
  [[ -e "$1" ]] || {
    echo "missing expected file: $1" >&2
    exit 1
  }
}

assert_not_exists() {
  [[ ! -e "$1" ]] || {
    echo "unexpected file exists: $1" >&2
    exit 1
  }
}

PROFILE_COLUMNS=$'name\ttaxid\trank\treads\tdepth\tabundance\tbreadth\texpected_breadth\tbreadth_pvalue\tfixed_chunk_breadth\tflex_chunk_breadth\texpected_flex_chunk_breadth\tchunk_pvalue\tevidence_depth\tbest_fit_genome\tmean_identity\testimated_ani'
HOST_COLUMNS=$'host_names\thost_taxids\tdiseases'

# A profile opens with its metadata line and its column names (by default the
# 17 columns of a profile without --host); every row below them has one field
# per column.
assert_profile_format() {
  local profile="$1"
  local columns="${2:-$PROFILE_COLUMNS}"
  awk -F '\t' -v columns="$columns" '
    NR == 1 { ok = $1 == "#metax_profile" && $2 ~ /^version=/ && /\tsequencer=/; next }
    NR == 2 { ok = ok && $0 == "#" columns; width = split(columns, names, "\t"); next }
    /^#/ || NF != width { ok = 0 }
    END { exit !(NR >= 2 && ok) }
  ' "$profile" || {
    echo "malformed profile $profile" >&2
    cat "$profile" >&2
    exit 1
  }
}

# A well-formed profile whose only row is the given taxon.
assert_single_row_profile() {
  local profile="$1" name="$2" taxid="$3" rank="$4"
  assert_profile_format "$profile"
  awk -F '\t' -v name="$name" -v taxid="$taxid" -v rank="$rank" '
    /^#/ { next }
    { rows++ }
    $1 == name && $2 == taxid && $3 == rank { matched++ }
    END { exit !(rows == 1 && matched == 1) }
  ' "$profile" || {
    echo "unexpected profile contents in $profile" >&2
    cat "$profile" >&2
    exit 1
  }
}

assert_single_species_profile() {
  assert_single_row_profile "$1" SpeciesA 100 species
}

assert_single_genus_profile() {
  assert_single_row_profile "$1" GenusOnly 300 genus
}

# The decisions table opens with its metadata line and column names, and every
# row of the profile has a decided row of the same rank and taxid that was
# kept.
assert_decisions_cover_profile() {
  local decisions="$1" profile="$2"
  awk -F '\t' '
    FILENAME == ARGV[1] && FNR == 1 { ok = $1 == "#metax_decisions" && $2 ~ /^version=/; next }
    FILENAME == ARGV[1] && FNR == 2 {
      for (i = 1; i <= NF; i++) column[$i] = i
      ok = ok && column["#rank"] == 1 && column["taxid"] == 2 && ("decided" in column) && ("decision" in column)
      next
    }
    FILENAME == ARGV[1] {
      if ($(column["decided"]) == "true" && $(column["decision"]) == "kept") kept[$1 SUBSEP $2] = 1
      next
    }
    /^#/ { next }
    { rows++; if (!(($3 SUBSEP $2) in kept)) ok = 0 }
    END { exit !(ok && rows > 0) }
  ' "$decisions" "$profile" || {
    echo "the decisions table $decisions does not cover the rows of $profile" >&2
    cat "$decisions" "$profile" >&2
    exit 1
  }
}

# A run with the given profile options fails, and its error output contains
# the given text.
assert_profile_error() {
  local name="$1" expected="$2"
  shift 2
  if run_metax profile --db "$WORK/db_one" --dmp-dir "$WORK/tax" -i "$WORK/se.fq" \
    -o "$WORK/$name" -t 2 "$@" 2>"$WORK/$name.err"; then
    echo "profile with $* unexpectedly succeeded" >&2
    exit 1
  fi
  grep -q -F -- "$expected" "$WORK/$name.err" || {
    echo "unexpected error for $*:" >&2
    cat "$WORK/$name.err" >&2
    exit 1
  }
  assert_not_exists "$WORK/$name.profile.txt"
}

mkdir -p "$WORK"
mkdir -p "$WORK/stage"
make_taxdump "$WORK/tax"

REF_SEQ="$(dna_repeat 1200)"
READ_120="${REF_SEQ:100:120}"
READ_900="${REF_SEQ:200:900}"

cat >"$WORK/one.fa" <<EOF
>G1|contig1|t:101;s:100|1200
$REF_SEQ
EOF

cat >"$WORK/multicontig.fa" <<EOF
>G2|contigA|s:200
${REF_SEQ:0:600}
>G2|contigB|s:200
${REF_SEQ:600:600}
EOF

cat >"$WORK/genus.fa" <<EOF
>G3|contig1|g:300
$REF_SEQ
EOF

cat >"$WORK/bad.fa" <<EOF
>bad|missing_taxonomy
$REF_SEQ
EOF

fastq r1 "$READ_120" >"$WORK/se.fq"
fastq r1/1 "$READ_120" >"$WORK/r1.fq"
fastq r1/2 "$READ_120" >"$WORK/r2.fq"
fastq long1 "$READ_900" >"$WORK/long.fq"
{
  printf '@long_adaptive\n%s\n+\n' "$READ_900"
  printf '%*s\n' "${#READ_900}" '' | tr ' ' '+'
} >"$WORK/long_adaptive.fq"

COMMON_INDEX_ARGS=(
  -k 11
  -w 11
  --max-occ 100
  --strobe-start-dist-min 20
  --strobe-start-dist-max 60
  --temp-dir "$WORK/stage"
  --dmp-dir "$WORK/tax"
  -t 2
)

run_metax index "$WORK/one.fa" -o "$WORK/db_one" "${COMMON_INDEX_ARGS[@]}"
assert_exists "$WORK/db_one/metax_db.tsv.gz"
assert_exists "$WORK/db_one/metax_db.mxi"
assert_not_exists "$WORK/db_one/metax_db.json"
assert_not_exists "$WORK/db_one/metax_db.tsv"
assert_not_exists "$WORK/db_one/metax_db.mmi"
assert_not_exists "$WORK/db_one/metax_reference.fasta"
assert_not_exists "$WORK/db_one/metax_reference.fasta.gz"
assert_exists "$WORK/db_one/taxdmp.tar.gz"
assert_not_exists "$WORK/db_one/taxdmp"
assert_exists "$WORK/db_one/metax_taxonomy.tsv.gz"

run_metax index "$WORK/multicontig.fa" -o "$WORK/db_multi" "${COMMON_INDEX_ARGS[@]}"
run_metax index "$WORK/genus.fa" -o "$WORK/db_genus" "${COMMON_INDEX_ARGS[@]}"
run_metax index "$WORK/one.fa" -o "$WORK/db_frac" -f 0.5 -s 7 -l 200 -m 50 "${COMMON_INDEX_ARGS[@]}"
assert_exists "$WORK/db_frac/metax_db.tsv.gz"
assert_exists "$WORK/db_frac/metax_db.mxi"
assert_exists "$WORK/db_frac/metax_reference.fasta.gz"
assert_not_exists "$WORK/db_frac/metax_reference.fasta"

run_metax info --db "$WORK/db_one" > "$WORK/info.tsv"
grep -q '^info_type	metax_database$' "$WORK/info.tsv"
grep -q '^metadata_reference_count	1$' "$WORK/info.tsv"
grep -q '^index_present	true$' "$WORK/info.tsv"
grep -q '^taxonomy_dump	taxdmp.tar.gz$' "$WORK/info.tsv"
grep -q '^taxonomy_record_taxa	2$' "$WORK/info.tsv"
grep -q '^taxonomy_record_metax_version	' "$WORK/info.tsv"
[[ "$(grep -c '^taxonomy_dump_file	' "$WORK/info.tsv")" == 4 ]]

if run_metax index "$WORK/bad.fa" -o "$WORK/db_bad" "${COMMON_INDEX_ARGS[@]}" 2>"$WORK/bad.err"; then
  echo "malformed header unexpectedly succeeded" >&2
  exit 1
fi
assert_not_exists "$WORK/db_bad"

run_metax profile \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/se.fq" \
  -o "$WORK/se" \
  -t 2 \
  --temp-dir "$WORK/stage" \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_exists "$WORK/se.profile.txt"
assert_single_species_profile "$WORK/se.profile.txt"
assert_not_exists "$WORK/se.classify.txt.gz"
assert_not_exists "$WORK/se.mapping.tsv"
assert_not_exists "$WORK/se.mapping.tsv.gz"
assert_not_exists "$WORK/se.mapping.profile.json"
if find "$WORK/stage" -mindepth 1 -print -quit | grep -q .; then
  echo "temporary staged mapping output was not cleaned up" >&2
  exit 1
fi

# Without --dmp-dir, the dump stored with the database is used.
run_metax profile \
  --db "$WORK/db_one" \
  -i "$WORK/se.fq" \
  -o "$WORK/se_stored_taxonomy" \
  -t 2 \
  --temp-dir "$WORK/stage" \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_single_species_profile "$WORK/se_stored_taxonomy.profile.txt"

# The stored archive unpacked by hand into taxdmp/ is used as well.
cp -r "$WORK/db_one" "$WORK/db_packed"
mkdir "$WORK/db_packed/taxdmp"
tar -xzf "$WORK/db_packed/taxdmp.tar.gz" -C "$WORK/db_packed/taxdmp"
run_metax profile \
  --db "$WORK/db_packed" \
  -i "$WORK/se.fq" \
  -o "$WORK/se_packed_taxonomy" \
  -t 2 \
  --temp-dir "$WORK/stage" \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_single_species_profile "$WORK/se_packed_taxonomy.profile.txt"

# A read-only database works when the index cache is placed elsewhere.
mkdir "$WORK/db_readonly" "$WORK/cache"
cp "$WORK/db_one/metax_db.mxi" "$WORK/db_one/metax_db.tsv.gz" "$WORK/db_one/metax_taxonomy.tsv.gz" \
  "$WORK/db_one/taxdmp.tar.gz" "$WORK/db_readonly/"
chmod a-w "$WORK/db_readonly"
run_metax profile --db "$WORK/db_readonly" -i "$WORK/se.fq" -o "$WORK/se_readonly" -t 2 \
  --temp-dir "$WORK/stage" --min-identity 0.90 --mapped-len 50 --min-mapped-fraction 0.50 \
  -- --index-cache "$WORK/cache/db_readonly.metamcache"
chmod u+w "$WORK/db_readonly"
assert_single_species_profile "$WORK/se_readonly.profile.txt"
assert_exists "$WORK/cache/db_readonly.metamcache"
[[ "$(ls -A "$WORK/db_readonly" | wc -l)" == 4 ]]

# A database built without a dump can be given one.
mkdir "$WORK/db_attach"
cp "$WORK/db_one/metax_db.mxi" "$WORK/db_one/metax_db.tsv.gz" "$WORK/db_attach/"
run_metax attach-taxonomy -d "$WORK/db_attach" --dmp-dir "$WORK/tax"
assert_exists "$WORK/db_attach/taxdmp.tar.gz"
assert_exists "$WORK/db_attach/metax_taxonomy.tsv.gz"
if run_metax attach-taxonomy -d "$WORK/db_attach" --dmp-dir "$WORK/tax" 2>"$WORK/attach.err"; then
  echo "a second attach-taxonomy unexpectedly succeeded" >&2
  exit 1
fi
run_metax profile --db "$WORK/db_attach" -i "$WORK/se.fq" -o "$WORK/se_attach" -t 2 \
  --temp-dir "$WORK/stage" --min-identity 0.90 --mapped-len 50 --min-mapped-fraction 0.50
assert_single_species_profile "$WORK/se_attach.profile.txt"

# A dump that names a database taxon differently gives a warning; one that
# lacks database taxa stops the run before mapping.
cp -r "$WORK/tax" "$WORK/tax_renamed"
sed -i 's/SpeciesA/SpeciesRenamed/' "$WORK/tax_renamed/names.dmp"
run_metax profile --db "$WORK/db_one" --dmp-dir "$WORK/tax_renamed" \
  -i "$WORK/se.fq" -o "$WORK/se_renamed" -t 2 --temp-dir "$WORK/stage" \
  --min-identity 0.90 --mapped-len 50 --min-mapped-fraction 0.50
grep -q "SpeciesRenamed" "$WORK/se_renamed.profile.txt"
grep -q "differs from the one the database was built with" "$WORK/se_renamed.log"
cp -r "$WORK/tax" "$WORK/tax_missing"
sed -i '/^10[01]\t/d' "$WORK/tax_missing/nodes.dmp" "$WORK/tax_missing/names.dmp"
if run_metax profile --db "$WORK/db_one" --dmp-dir "$WORK/tax_missing" \
  -i "$WORK/se.fq" -o "$WORK/se_missing" -t 2 --temp-dir "$WORK/stage" \
  2>"$WORK/missing.err"; then
  echo "profiling with a dump that lacks database taxa unexpectedly succeeded" >&2
  exit 1
fi
grep -q "2 are not in it" "$WORK/missing.err"
grep -q -- "--skip-missing-taxa" "$WORK/missing.err"
assert_not_exists "$WORK/se_missing.profile.txt"
run_metax profile --db "$WORK/db_one" --dmp-dir "$WORK/tax_missing" --skip-missing-taxa \
  -i "$WORK/se.fq" -o "$WORK/se_missing" -t 2 --temp-dir "$WORK/stage" \
  --min-identity 0.90 --mapped-len 50 --min-mapped-fraction 0.50
grep -q "Ignored 1 reference sequence(s) of 1 genome(s)" "$WORK/se_missing.log"

run_metax profile \
  --db "$WORK/db_genus" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/se.fq" \
  -o "$WORK/genus" \
  -t 2 \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_exists "$WORK/genus.profile.txt"
assert_single_genus_profile "$WORK/genus.profile.txt"

run_metax profile \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/se.fq" \
  -o "$WORK/se_stream" \
  -t 2 \
  --stream-mappings \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_exists "$WORK/se_stream.profile.txt"
assert_single_species_profile "$WORK/se_stream.profile.txt"
assert_not_exists "$WORK/se_stream.classify.txt.gz"
assert_not_exists "$WORK/se_stream.mapping.tsv.gz"
cmp "$WORK/se.profile.txt" "$WORK/se_stream.profile.txt"

run_metax profile \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/se.fq" \
  -o "$WORK/keep" \
  -t 2 \
  --save-mapping \
  --save-classification \
  --save-decisions \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_exists "$WORK/keep.mapping.tsv.gz"
assert_exists "$WORK/keep.classify.txt.gz"
[[ -n "$(gzip -dc "$WORK/keep.classify.txt.gz")" ]]
assert_single_species_profile "$WORK/keep.profile.txt"
assert_decisions_cover_profile "$WORK/keep.decisions.tsv" "$WORK/keep.profile.txt"
assert_not_exists "$WORK/keep.mapping.profile.json"
if ! gzip -cd "$WORK/keep.mapping.tsv.gz" | awk -F '\t' '
  $1 == "#metax_mapping" && $2 == "version=5" && $3 == "sequencer=Illumina" && $4 == "paired=false" && $5 == "mapper_mode=default" && $6 == "strain=false" && $7 == "profiler_identity=0.9" && $8 == "profiler_mapped_len=50" && $9 == "profiler_fraction=0.5" && $10 == "sampled_mean_base_accuracy=NA" { file_meta = 1 }
  $1 == "#metax_mapping_database" && $2 == "fractional=false" && $3 == "sequences=1" && $4 == "genomes=1" { db_meta = 1 }
  $1 == "#metax_mapping_reference" { exit 1 }
  $1 == "#metax_mapping_sample" && $2 == "total_reads=1" && $3 == "mapped_reads=1" && $4 ~ /^emitted_records=[1-9][0-9]*$/ { sample_meta = 1 }
  !/^#/ {
    mapping = 1
    if (NF != 11 || $9 != 0 || $10 != 0 || $11 != 0) exit 1
    if ($5 ~ /0$/ || $5 ~ /\.$/) exit 1
  }
  END { if (!(file_meta && db_meta && sample_meta && mapping)) exit 1 }
'; then
  echo "saved mapping output is not a clean coordinate mapping file" >&2
  gzip -cd "$WORK/keep.mapping.tsv.gz" >&2
  exit 1
fi
cmp "$WORK/se.profile.txt" "$WORK/keep.profile.txt"

run_metax profile \
  --reuse-mapping "$WORK/keep.mapping.tsv.gz" \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -o "$WORK/reused" \
  -t 2
assert_exists "$WORK/reused.profile.txt"
assert_single_species_profile "$WORK/reused.profile.txt"
assert_not_exists "$WORK/reused.classify.txt.gz"
assert_not_exists "$WORK/reused.decisions.tsv"
cmp "$WORK/keep.profile.txt" "$WORK/reused.profile.txt"

if find "$WORK" -maxdepth 1 -name '.metax-mapping-*.partial.tsv.gz' -print -quit | grep -q .; then
  echo "partial saved mapping output was not cleaned up" >&2
  exit 1
fi

run_metax profile \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/r1.fq" \
  -I "$WORK/r2.fq" \
  -o "$WORK/pe" \
  -t 2 \
  --save-mapping \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_exists "$WORK/pe.profile.txt"
assert_single_species_profile "$WORK/pe.profile.txt"
assert_not_exists "$WORK/pe.classify.txt.gz"
assert_exists "$WORK/pe.mapping.tsv.gz"
if ! gzip -cd "$WORK/pe.mapping.tsv.gz" | awk -F '\t' '
  $1 == "#metax_mapping" && $2 == "version=5" && $4 == "paired=true" { file_meta = 1 }
  !/^#/ {
    if (NF != 11 || $8 != "=" || ($9 != 1 && $9 != 2) || $10 !~ /^[1-9][0-9]*$/ || $11 != 0) exit 1
    key = $1 SUBSEP $10
    if ($9 == 1) mate1[key]++
    if ($9 == 2) mate2[key]++
    pair[key] = 1
  }
  END {
    if (!file_meta) exit 1
    for (key in pair) {
      if (mate1[key] != 1 || mate2[key] != 1) exit 1
    }
  }
'; then
  echo "paired mapping output does not contain explicit one-to-one pair IDs" >&2
  gzip -cd "$WORK/pe.mapping.tsv.gz" >&2
  exit 1
fi

run_metax profile \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/r1.fq" \
  -I "$WORK/r2.fq" \
  -o "$WORK/pe_strain" \
  -t 2 \
  --strain \
  --save-decisions \
  --very-verbose \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
# The genome level reports the read pair on the species' genome, a strain row
# named after its strain taxid, and the species keeps no reads of its own.
assert_exists "$WORK/pe_strain.profile.txt"
assert_single_row_profile "$WORK/pe_strain.profile.txt" StrainA 101 strain
assert_exists "$WORK/pe_strain.decisions.tsv"
assert_decisions_cover_profile "$WORK/pe_strain.decisions.tsv" "$WORK/pe_strain.profile.txt"
awk -F '\t' '$1 == "species" && $2 == "100" && $6 == "true" { found = 1 } END { exit !found }' \
  "$WORK/pe_strain.decisions.tsv"
assert_not_exists "$WORK/pe_strain.rprofile.txt"
assert_not_exists "$WORK/pe_strain.classify.txt.gz"

# Options that no longer exist stop the run and say what to use instead.
assert_profile_error removed_engine "the level-wise engine always runs" --level-engine
assert_profile_error renamed_decisions "Use --save-decisions instead." --save-raw-profile

# --host filters the profile by the pathogen table, which it needs.
assert_profile_error host_alone "--host needs --pathogen-table" --host 9606
printf 'pathogen_taxid\thost_taxids\thost_names\tdiseases\n100\t9606\tHomo sapiens\tdisease A\n200\t9913\tBos taurus\tdisease B\n' \
  >"$WORK/pathogens.tsv"
run_metax profile \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/se.fq" \
  -o "$WORK/host" \
  -t 2 \
  --host 9606 \
  --pathogen-table "$WORK/pathogens.tsv" \
  --save-decisions \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
# A pathogen of the host is kept and annotated; the files take the prefix
# <prefix>.pathogen and the profile has 20 columns.
assert_not_exists "$WORK/host.profile.txt"
assert_profile_format "$WORK/host.pathogen.profile.txt" "$PROFILE_COLUMNS"$'\t'"$HOST_COLUMNS"
head -n 1 "$WORK/host.pathogen.profile.txt" | grep -q $'\thost=9606$'
awk -F '\t' '
  /^#/ { next }
  { rows++ }
  $1 == "SpeciesA" && $2 == "100" && $18 == "Homo sapiens" && $19 == "9606" && $20 == "disease A" { matched++ }
  END { exit !(rows == 1 && matched == 1) }
' "$WORK/host.pathogen.profile.txt" || {
  echo "unexpected host profile contents" >&2
  cat "$WORK/host.pathogen.profile.txt" >&2
  exit 1
}
assert_decisions_cover_profile "$WORK/host.pathogen.decisions.tsv" "$WORK/host.pathogen.profile.txt"
# A species the table lists for another host only is dropped.
run_metax profile --db "$WORK/db_multi" --dmp-dir "$WORK/tax" -i "$WORK/se.fq" -o "$WORK/multi" -t 2 \
  --min-identity 0.90 --mapped-len 50 --min-mapped-fraction 0.50
assert_single_row_profile "$WORK/multi.profile.txt" SpeciesB 200 species
run_metax profile --db "$WORK/db_multi" --dmp-dir "$WORK/tax" -i "$WORK/se.fq" -o "$WORK/multi_host" -t 2 \
  --host 9606 --pathogen-table "$WORK/pathogens.tsv" \
  --min-identity 0.90 --mapped-len 50 --min-mapped-fraction 0.50
assert_profile_format "$WORK/multi_host.pathogen.profile.txt" "$PROFILE_COLUMNS"$'\t'"$HOST_COLUMNS"
[[ "$(grep -vc '^#' "$WORK/multi_host.pathogen.profile.txt")" == 0 ]]

run_metax profile \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/long.fq" \
  --sequencer Nanopore \
  -o "$WORK/long" \
  -t 2 \
  --min-identity 0.90 \
  --mapped-len 100 \
  --min-mapped-fraction 0.30
assert_exists "$WORK/long.profile.txt"
assert_single_species_profile "$WORK/long.profile.txt"
assert_not_exists "$WORK/long.classify.txt.gz"

run_metax profile \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/long_adaptive.fq" \
  --sequencer Nanopore \
  -o "$WORK/long_saved" \
  -t 2 \
  --save-mapping
assert_exists "$WORK/long_saved.profile.txt"
assert_exists "$WORK/long_saved.mapping.tsv.gz"
assert_single_species_profile "$WORK/long_saved.profile.txt"
if ! gzip -cd "$WORK/long_saved.mapping.tsv.gz" | awk -F '\t' '
  $1 == "#metax_mapping" {
    for (i = 2; i <= NF; i++) {
      split($i, field, "=")
      metadata[field[1]] = field[2]
    }
  }
  END {
    if (metadata["version"] != "5") exit 1
    if (metadata["sequencer"] != "Nanopore") exit 1
    if (metadata["paired"] != "false") exit 1
    if (metadata["mapper_mode"] != "default") exit 1
    if (metadata["strain"] != "false") exit 1
    if (metadata["profiler_identity"] < 0.8699 || metadata["profiler_identity"] > 0.8701) exit 1
    if (metadata["profiler_mapped_len"] != 500) exit 1
    if (metadata["profiler_fraction"] != 0.2) exit 1
    if (metadata["sampled_mean_base_accuracy"] < 0.8999 || metadata["sampled_mean_base_accuracy"] > 0.9001) exit 1
  }
'; then
  echo "saved long-read mapping output is missing adaptive profiler metadata" >&2
  gzip -cd "$WORK/long_saved.mapping.tsv.gz" >&2
  exit 1
fi

run_metax profile \
  --reuse-mapping "$WORK/long_saved.mapping.tsv.gz" \
  --db "$WORK/db_one" \
  --dmp-dir "$WORK/tax" \
  -o "$WORK/long_reused" \
  -t 2
assert_exists "$WORK/long_reused.profile.txt"
assert_not_exists "$WORK/long_reused.classify.txt.gz"
cmp "$WORK/long_saved.profile.txt" "$WORK/long_reused.profile.txt"

run_metax profile \
  --db "$WORK/db_genus" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/long.fq" \
  --sequencer Nanopore \
  -o "$WORK/long_genus" \
  -t 2 \
  --min-identity 0.90 \
  --mapped-len 100 \
  --min-mapped-fraction 0.30
assert_exists "$WORK/long_genus.profile.txt"
assert_single_genus_profile "$WORK/long_genus.profile.txt"
assert_not_exists "$WORK/long_genus.classify.txt.gz"

run_metax profile \
  --db "$WORK/db_frac" \
  --dmp-dir "$WORK/tax" \
  -i "$WORK/se.fq" \
  -o "$WORK/frac" \
  -t 2 \
  --save-mapping \
  --save-classification \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_exists "$WORK/frac.profile.txt"
assert_exists "$WORK/frac.mapping.tsv.gz"
assert_not_exists "$WORK/frac.classify.txt.gz"
if ! gzip -cd "$WORK/frac.mapping.tsv.gz" | awk -F '\t' '
  $1 == "#metax_mapping_database" && $2 == "fractional=true" && $3 ~ /^sequences=[1-9][0-9]*$/ && $4 ~ /^genomes=[1-9][0-9]*$/ { db_meta = 1 }
  !/^#/ {
    if (NF != 11 || $9 != 0 || $10 != 0 || $11 != 0) exit 1
    if (split($2, fields, "|") != 6) exit 1
    if (fields[5] !~ /^[1-9][0-9]*$/ || fields[6] !~ /^[1-9][0-9]*$/) exit 1
    mapping = 1
  }
  $1 == "#metax_mapping_reference" { exit 1 }
  END { if (!(db_meta && mapping)) exit 1 }
'; then
  echo "fractional mapping output is missing required scaling metadata" >&2
  gzip -cd "$WORK/frac.mapping.tsv.gz" >&2
  exit 1
fi

run_metax profile \
  --reuse-mapping "$WORK/frac.mapping.tsv.gz" \
  --db "$WORK/db_frac" \
  --dmp-dir "$WORK/tax" \
  -o "$WORK/frac_reused" \
  -t 2 \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_exists "$WORK/frac_reused.profile.txt"
assert_not_exists "$WORK/frac_reused.classify.txt.gz"
cmp "$WORK/frac.profile.txt" "$WORK/frac_reused.profile.txt"

echo "Metax smoke tests passed: $WORK"

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

assert_single_species_profile() {
  local profile="$1"
  awk -F '\t' '
    $1 == "SpeciesA" && $2 == "100" && $3 == "species" { matched++ }
    END { exit !(NR == 1 && matched == 1) }
  ' "$profile" || {
    echo "unexpected profile contents in $profile" >&2
    cat "$profile" >&2
    exit 1
  }
}

assert_single_genus_profile() {
  local profile="$1"
  awk -F '\t' '
    $1 == "GenusOnly" && $2 == "300" && $3 == "genus" { matched++ }
    END { exit !(NR == 1 && matched == 1) }
  ' "$profile" || {
    echo "unexpected genus profile contents in $profile" >&2
    cat "$profile" >&2
    exit 1
  }
}

assert_single_strain_raw_profile() {
  local profile="$1"
  awk -F '\t' '
    $1 == "StrainA" && $2 == "101" && $3 == "no rank" { matched++ }
    END { exit !(NR == 1 && matched == 1) }
  ' "$profile" || {
    echo "unexpected raw strain profile contents in $profile" >&2
    cat "$profile" >&2
    exit 1
  }
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
assert_not_exists "$WORK/se.classify.txt"
assert_not_exists "$WORK/se.mapping.tsv"
assert_not_exists "$WORK/se.mapping.tsv.gz"
assert_not_exists "$WORK/se.mapping.profile.json"
if find "$WORK/stage" -mindepth 1 -print -quit | grep -q .; then
  echo "temporary staged mapping output was not cleaned up" >&2
  exit 1
fi

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
assert_not_exists "$WORK/se_stream.classify.txt"
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
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_exists "$WORK/keep.mapping.tsv.gz"
assert_exists "$WORK/keep.classify.txt"
[[ -s "$WORK/keep.classify.txt" ]]
assert_single_species_profile "$WORK/keep.profile.txt"
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
assert_not_exists "$WORK/reused.classify.txt"
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
assert_not_exists "$WORK/pe.classify.txt"
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
  --save-raw-profile \
  --very-verbose \
  --min-identity 0.90 \
  --mapped-len 50 \
  --min-mapped-fraction 0.50
assert_exists "$WORK/pe_strain.profile.txt"
[[ ! -s "$WORK/pe_strain.profile.txt" ]]
assert_exists "$WORK/pe_strain.rprofile.txt"
assert_single_strain_raw_profile "$WORK/pe_strain.rprofile.txt"
assert_not_exists "$WORK/pe_strain.classify.txt"

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
assert_not_exists "$WORK/long.classify.txt"

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
assert_not_exists "$WORK/long_reused.classify.txt"
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
assert_not_exists "$WORK/long_genus.classify.txt"

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
assert_not_exists "$WORK/frac.classify.txt"
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
assert_not_exists "$WORK/frac_reused.classify.txt"
cmp "$WORK/frac.profile.txt" "$WORK/frac_reused.profile.txt"

echo "Metax smoke tests passed: $WORK"

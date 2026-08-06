#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_project_config

usage() {
  cat <<'EOF'
Usage:
  PROJECT_ROOT=/path/to/nanopore_16s_mixed \
    bash scripts/09_explore_single_sample.sh SAMPLE_ID [MIN_LEN] [MAX_LEN] [MIN_Q]

Example:
  PROJECT_ROOT=/absolute/path/to/nanopore_16s_mixed \
    bash scripts/09_explore_single_sample.sh 2PJGZM_2_sample_2 1000 2000 10

The script creates a timestamped directory under results/exploratory_single_sample/.
It never modifies the raw FASTQ or primary custom/broad Emu outputs.
EOF
}

[[ "${1:-}" != "-h" && "${1:-}" != "--help" ]] || {
  usage
  exit 0
}
[[ "$#" -ge 1 && "$#" -le 4 ]] || {
  usage >&2
  exit 2
}

sample_id="$1"
explore_min_length="${2:-1000}"
explore_max_length="${3:-2000}"
explore_min_quality="${4:-10}"

[[ "${sample_id}" =~ ^[A-Za-z0-9._-]+$ ]] ||
  die "Unsafe sample ID: ${sample_id}"
[[ "${explore_min_length}" =~ ^[0-9]+$ ]] ||
  die "MIN_LEN must be a positive integer"
[[ "${explore_max_length}" =~ ^[0-9]+$ ]] ||
  die "MAX_LEN must be a positive integer"
[[ "${explore_min_quality}" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
  die "MIN_Q must be numeric"
(( explore_min_length > 0 )) || die "MIN_LEN must be greater than zero"
(( explore_max_length > explore_min_length )) ||
  die "MAX_LEN must be greater than MIN_LEN"

manifest="${PROJECT_ROOT}/metadata/fastq_manifest.tsv"
custom_db="${PROJECT_ROOT}/databases/custom_emu/consortium5_custom_emu"
broad_db="${BROAD_EMU_DB:-${PROJECT_ROOT}/databases/broad_emu/emu_broad_db}"

[[ -f "${manifest}" ]] || die "Missing FASTQ manifest: ${manifest}"
[[ -d "${custom_db}" ]] || die "Missing custom Emu database: ${custom_db}"
[[ -f "${custom_db}/species_taxid.fasta" ]] ||
  die "Invalid custom Emu database: ${custom_db}"
[[ -d "${broad_db}" ]] || die "Missing broad Emu database: ${broad_db}"
[[ -f "${broad_db}/species_taxid.fasta" && -f "${broad_db}/taxonomy.tsv" ]] ||
  die "Invalid broad Emu database: ${broad_db}"

have_cmd emu || die "emu is not available on PATH"
have_cmd seqkit || die "seqkit is not available on PATH"
have_cmd gzip || die "gzip is not available on PATH"

mapfile -t manifest_rows < <(
  awk -F'\t' -v sample="${sample_id}" 'NR > 1 && $1 == sample { print }' \
    "${manifest}"
)
[[ "${#manifest_rows[@]}" -eq 1 ]] ||
  die "Expected exactly one manifest row for ${sample_id}; found ${#manifest_rows[@]}"

IFS=$'\t' read -r \
  manifest_sample original_filename raw_fastq compression \
  file_size_bytes raw_md5 raw_read_count raw_mean_length \
  raw_median_length raw_min_length raw_max_length \
  <<< "${manifest_rows[0]}"

[[ "${manifest_sample}" == "${sample_id}" ]] ||
  die "Manifest parsing failed for ${sample_id}"
[[ -f "${raw_fastq}" ]] || die "Raw FASTQ is missing: ${raw_fastq}"

run_stamp="$(timestamp)"
run_label="len${explore_min_length}_${explore_max_length}_q${explore_min_quality}"
outdir="${PROJECT_ROOT}/results/exploratory_single_sample/${sample_id}/${run_stamp}_${run_label}"
filtered_fastq="${outdir}/${sample_id}.${run_label}.fastq.gz"
custom_out="${outdir}/custom_emu"
broad_out="${outdir}/broad_emu"
log_file="${PROJECT_ROOT}/logs/${run_stamp}_09_explore_${sample_id}.log"

[[ ! -e "${outdir}" ]] || die "Output directory already exists: ${outdir}"
mkdir -p "${outdir}" "${custom_out}" "${broad_out}" "${PROJECT_ROOT}/logs"
exec > >(tee -a "${log_file}") 2>&1

printf 'Exploratory single-sample analysis\n'
printf 'Started: %s\n' "$(date -Is)"
printf 'Project root: %s\n' "${PROJECT_ROOT}"
printf 'Sample ID: %s\n' "${sample_id}"
printf 'Raw FASTQ: %s\n' "${raw_fastq}"
printf 'Raw FASTQ MD5 from manifest: %s\n' "${raw_md5}"
printf 'Primary window: %s-%s bp; Q%s when a quality tool was available\n' \
  "${MIN_LENGTH}" "${MAX_LENGTH}" "${MIN_QUALITY}"
printf 'Exploratory window: %s-%s bp; Q%s when a quality tool is available\n' \
  "${explore_min_length}" "${explore_max_length}" "${explore_min_quality}"
printf 'Output directory: %s\n\n' "${outdir}"

emu --version > "${outdir}/emu_version.txt" 2>&1 || true
emu abundance --help > "${outdir}/emu_abundance_help.txt" 2>&1 || true
seqkit version > "${outdir}/seqkit_version.txt" 2>&1 || true

seqkit stats -a "${raw_fastq}" |
  tee "${outdir}/raw_seqkit_stats.txt"

filter_tool="seqkit_length_only"
if have_cmd chopper; then
  chopper --help > "${outdir}/chopper_help.txt" 2>&1 || true
  filter_tool="chopper"
elif have_cmd NanoFilt; then
  NanoFilt --help > "${outdir}/NanoFilt_help.txt" 2>&1 || true
  filter_tool="NanoFilt"
fi

stream_raw_fastq() {
  if [[ "${compression}" == "gzip" || "${raw_fastq}" == *.gz ]]; then
    gzip -dc "${raw_fastq}"
  else
    command cat "${raw_fastq}"
  fi
}

if [[ "${filter_tool}" == "chopper" ]]; then
  stream_raw_fastq |
    chopper \
      -l "${explore_min_length}" \
      --maxlength "${explore_max_length}" \
      -q "${explore_min_quality}" |
    gzip -c > "${filtered_fastq}"
  filtering_note="${explore_min_length}-${explore_max_length} bp; mean Q >= ${explore_min_quality} with Chopper"
elif [[ "${filter_tool}" == "NanoFilt" ]]; then
  stream_raw_fastq |
    NanoFilt \
      -l "${explore_min_length}" \
      --maxlength "${explore_max_length}" \
      -q "${explore_min_quality}" |
    gzip -c > "${filtered_fastq}"
  filtering_note="${explore_min_length}-${explore_max_length} bp; mean Q >= ${explore_min_quality} with NanoFilt"
else
  warn "Chopper and NanoFilt are unavailable; applying length filtering only."
  stream_raw_fastq |
    seqkit seq \
      -m "${explore_min_length}" \
      -M "${explore_max_length}" |
    gzip -c > "${filtered_fastq}"
  filtering_note="${explore_min_length}-${explore_max_length} bp; no FASTQ quality threshold applied"
fi

gzip -t "${filtered_fastq}"
seqkit stats -a "${filtered_fastq}" |
  tee "${outdir}/exploratory_filtered_seqkit_stats.txt"

retained_reads="$(
  seqkit stats -T "${filtered_fastq}" |
    awk 'NR == 2 { print $4 }'
)"
[[ "${retained_reads:-0}" -gt 0 ]] ||
  die "No reads passed the exploratory filter"

retained_pct="$(
  awk \
    -v retained="${retained_reads}" \
    -v raw="${raw_read_count}" \
    'BEGIN { if (raw > 0) printf "%.3f", 100 * retained / raw; else print "NA" }'
)"

{
  printf 'sample_id\traw_reads\tretained_reads\tretained_pct\tfilter_tool'
  printf '\tminimum_length\tmaximum_length\tminimum_quality\tfiltering_note\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${sample_id}" \
    "${raw_read_count}" \
    "${retained_reads}" \
    "${retained_pct}" \
    "${filter_tool}" \
    "${explore_min_length}" \
    "${explore_max_length}" \
    "${explore_min_quality}" \
    "${filtering_note}"
} > "${outdir}/exploratory_filter_summary.tsv"

help_text="$(< "${outdir}/emu_abundance_help.txt")"
has_output_dir=0
grep -q -- '--output-dir' <<< "${help_text}" && has_output_dir=1

run_emu() {
  local method="$1"
  local database="$2"
  local destination="$3"
  local basename="${sample_id}_${method}_${run_label}"
  local -a args=(
    emu abundance
    "${filtered_fastq}"
    --db "${database}"
    --threads "${THREADS}"
  )

  grep -q -- '--type' <<< "${help_text}" &&
    args+=(--type map-ont)
  grep -q -- '--min-align-len' <<< "${help_text}" &&
    args+=(--min-align-len "${explore_min_length}")
  grep -q -- '--max-align-len' <<< "${help_text}" &&
    args+=(--max-align-len "${explore_max_length}")
  if [[ "${method}" == "broad" ]] &&
    grep -q -- '--min-pid' <<< "${help_text}"; then
    args+=(--min-pid "${EMU_MIN_PID}")
  fi
  grep -q -- '--keep-read-assignments' <<< "${help_text}" &&
    args+=(--keep-read-assignments)
  grep -q -- '--keep-counts' <<< "${help_text}" &&
    args+=(--keep-counts)
  grep -q -- '--output-unclassified' <<< "${help_text}" &&
    args+=(--output-unclassified)
  grep -q -- '--output-basename' <<< "${help_text}" &&
    args+=(--output-basename "${basename}")

  if [[ "${has_output_dir}" -eq 1 ]]; then
    args+=(--output-dir "${destination}")
  fi

  printf '%q ' "${args[@]}" > "${destination}/emu_command.txt"
  printf '\n' >> "${destination}/emu_command.txt"

  if [[ "${has_output_dir}" -eq 1 ]]; then
    "${args[@]}" \
      > "${destination}/emu.stdout.log" \
      2> "${destination}/emu.stderr.log"
  else
    (
      cd "${destination}"
      "${args[@]}"
    ) > "${destination}/emu.stdout.log" 2> "${destination}/emu.stderr.log"
  fi
}

printf '\nRunning custom Emu...\n'
run_emu "custom" "${custom_db}" "${custom_out}"

printf 'Running broad Emu...\n'
run_emu "broad" "${broad_db}" "${broad_out}"

mapfile -t custom_tables < <(
  find "${custom_out}" -type f -name '*rel-abundance*.tsv' | sort
)
mapfile -t broad_tables < <(
  find "${broad_out}" -type f -name '*rel-abundance*.tsv' | sort
)
[[ "${#custom_tables[@]}" -eq 1 ]] ||
  die "Expected one custom Emu abundance table; found ${#custom_tables[@]}"
[[ "${#broad_tables[@]}" -eq 1 ]] ||
  die "Expected one broad Emu abundance table; found ${#broad_tables[@]}"

summarize_emu_table() {
  local method="$1"
  local table="$2"
  awk \
    -F'\t' \
    -v OFS='\t' \
    -v method="${method}" \
    -v total="${retained_reads}" '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          key = tolower($i)
          if (key == "tax_id" || key == "tax id") tax_col = i
          if (key == "species") species_col = i
          if (key == "genus") genus_col = i
          if (key ~ /estimated.*count/) count_col = i
          if (key == "abundance") abundance_col = i
        }
        if (!tax_col || !count_col) {
          print "Unable to identify Emu output columns in " FILENAME > "/dev/stderr"
          exit 1
        }
        next
      }
      {
        tax_id = $tax_col
        species = species_col ? $species_col : ""
        genus = genus_col ? $genus_col : ""
        count = $count_col + 0
        abundance = abundance_col ? $abundance_col + 0 : 0
        printf "%s\t%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
          method,
          tax_id,
          species,
          genus,
          count,
          total > 0 ? 100 * count / total : 0,
          100 * abundance
      }
    ' "${table}"
}

{
  printf 'method\ttax_id\tspecies\tgenus\testimated_counts'
  printf '\tpercentage_all_exploratory_reads\temu_abundance_pct\n'
  summarize_emu_table "custom" "${custom_tables[0]}"
  summarize_emu_table "broad" "${broad_tables[0]}"
} > "${outdir}/exploratory_emu_summary.tsv"

comparison_file="${outdir}/primary_vs_exploratory_emu.tsv"
{
  printf 'analysis_window\tmethod\ttax_id\tspecies\tgenus'
  printf '\testimated_counts\tpercentage_all_filtered_reads\n'

  primary_custom="${PROJECT_ROOT}/results/combined/custom_emu_expected_species_abundance_renormalized.tsv"
  if [[ -f "${primary_custom}" ]]; then
    awk \
      -F'\t' \
      -v OFS='\t' \
      -v sample="${sample_id}" '
        NR == 1 {
          for (i = 1; i <= NF; i++) {
            if ($i == "sample_id") sample_col = i
            if ($i == "species") species_col = i
            if ($i == "estimated_counts") count_col = i
            if ($i == "emu_abundance_pct_all_usable_reads") pct_col = i
          }
          next
        }
        $sample_col == sample {
          print "primary_1200_1800", "custom", "", $species_col, "",
            $count_col, $pct_col
        }
      ' "${primary_custom}"
  fi

  primary_broad="${PROJECT_ROOT}/results/combined/broad_emu_species_abundance.tsv"
  if [[ -f "${primary_broad}" ]]; then
    awk \
      -F'\t' \
      -v OFS='\t' \
      -v sample="${sample_id}" '
        NR == 1 {
          for (i = 1; i <= NF; i++) {
            if ($i == "sample_id") sample_col = i
            if ($i == "tax_id") tax_col = i
            if ($i == "species") species_col = i
            if ($i == "genus") genus_col = i
            if ($i == "estimated_counts") count_col = i
            if ($i == "abundance_pct_all_usable_reads") pct_col = i
          }
          next
        }
        $sample_col == sample {
          print "primary_1200_1800", "broad", $tax_col, $species_col,
            $genus_col, $count_col, $pct_col
        }
      ' "${primary_broad}"
  fi

  awk \
    -F'\t' \
    -v OFS='\t' \
    -v window="exploratory_${run_label}" '
      NR > 1 {
        print window, $1, $2, $3, $4, $5, $6
      }
    ' "${outdir}/exploratory_emu_summary.tsv"
} > "${comparison_file}"

{
  printf 'sample_id\t%s\n' "${sample_id}"
  printf 'run_started\t%s\n' "${run_stamp}"
  printf 'raw_fastq\t%s\n' "${raw_fastq}"
  printf 'raw_fastq_md5\t%s\n' "$(md5_file "${raw_fastq}")"
  printf 'filtered_fastq\t%s\n' "${filtered_fastq}"
  printf 'filtered_fastq_md5\t%s\n' "$(md5_file "${filtered_fastq}")"
  printf 'custom_database\t%s\n' "$(abs_path "${custom_db}")"
  printf 'broad_database\t%s\n' "$(abs_path "${broad_db}")"
  printf 'filter_tool\t%s\n' "${filter_tool}"
  printf 'filtering_note\t%s\n' "${filtering_note}"
  printf 'raw_reads\t%s\n' "${raw_read_count}"
  printf 'retained_reads\t%s\n' "${retained_reads}"
  printf 'retained_pct\t%s\n' "${retained_pct}"
  printf 'emu_version\t%s\n' "$(emu --version 2>&1 || true)"
  printf 'seqkit_version\t%s\n' "$(seqkit version 2>&1 || true)"
} > "${outdir}/run_manifest.tsv"

printf '\nExploratory Emu summary\n'
if have_cmd column; then
  column -t -s $'\t' "${outdir}/exploratory_emu_summary.tsv"
else
  command cat "${outdir}/exploratory_emu_summary.tsv"
fi

printf '\nPrimary versus exploratory comparison\n'
if have_cmd column; then
  column -t -s $'\t' "${comparison_file}"
else
  command cat "${comparison_file}"
fi

printf '\nCompleted: %s\n' "$(date -Is)"
printf 'Results: %s\n' "${outdir}"
printf 'Log: %s\n' "${log_file}"
printf '\nInterpret this run as exploratory only. Partial reads admitted by the wider\n'
printf 'window may be less taxonomically specific than the primary 1200-1800 bp run.\n'

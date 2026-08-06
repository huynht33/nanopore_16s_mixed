#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_project_config
make_log "04_qc_and_filter"

manifest="${PROJECT_ROOT}/metadata/fastq_manifest.tsv"
exclusions="${PROJECT_ROOT}/metadata/sample_exclusions.tsv"
[[ -f "${manifest}" ]] || die "Missing ${manifest}; run 03_inventory_fastq.sh first"
have_cmd seqkit || die "seqkit is required for length filtering"
seqkit version > "${PROJECT_ROOT}/logs/seqkit_version.txt" 2>&1 || true

mkdir -p "${PROJECT_ROOT}/fastq/filtered" "${PROJECT_ROOT}/qc/raw" "${PROJECT_ROOT}/qc/filtered"
qc_table="${PROJECT_ROOT}/qc/filtered/filter_summary.tsv"
printf 'sample_id\traw_reads\treads_after_length_filter\treads_after_quality_filter\tpercent_retained\tmean_retained_read_length\tmedian_retained_read_length\tfiltering_note\tanalysis_status\n' > "${qc_table}"

exclusion_reason() {
  local sample="$1"
  [[ -f "${exclusions}" ]] || return 0
  awk -F'\t' -v sample="${sample}" 'NR > 1 && $1 == sample { print $2; exit }' "${exclusions}"
}

filter_tool="seqkit_length_only"
if command -v chopper >/dev/null 2>&1; then
  chopper --help > "${PROJECT_ROOT}/logs/chopper_help.txt" 2>&1 || true
  filter_tool="chopper"
elif command -v NanoFilt >/dev/null 2>&1; then
  NanoFilt --help > "${PROJECT_ROOT}/logs/NanoFilt_help.txt" 2>&1 || true
  filter_tool="NanoFilt"
else
  warn "No Chopper or NanoFilt found; applying length filtering only with seqkit."
fi

tail -n +2 "${manifest}" | while IFS=$'\t' read -r sample original abs_path compression bytes md5 raw_reads mean_read median_read min_read max_read; do
  reason="$(exclusion_reason "${sample}")"
  if [[ -n "${reason}" ]]; then
    warn "Excluding ${sample}: ${reason}"
    printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\t%s\texcluded\n' \
      "${sample}" "${raw_reads}" "${reason}" >> "${qc_table}"
    continue
  fi

  sample_dir_raw="${PROJECT_ROOT}/qc/raw/${sample}"
  sample_dir_filt="${PROJECT_ROOT}/qc/filtered/${sample}"
  mkdir -p "${sample_dir_raw}" "${sample_dir_filt}"
  seqkit stats -a "${abs_path}" > "${sample_dir_raw}/seqkit_stats.txt"
  if command -v NanoPlot >/dev/null 2>&1; then
    NanoPlot --fastq "${abs_path}" -o "${sample_dir_raw}/NanoPlot" --threads "${THREADS}" > "${sample_dir_raw}/NanoPlot.stdout.log" 2> "${sample_dir_raw}/NanoPlot.stderr.log" || warn "NanoPlot failed for ${sample}; see logs"
  fi

  length_only_stats="$(
    seqkit seq -m "${MIN_LENGTH}" -M "${MAX_LENGTH}" "${abs_path}" |
      seqkit stats -T - 2>/dev/null |
      awk 'NR==2 {print $4}'
  )"
  length_retained="${length_only_stats:-0}"

  filtered="${PROJECT_ROOT}/fastq/filtered/${sample}.full_length_qc.fastq.gz"
  note=""
  if [[ "${filter_tool}" == "chopper" ]]; then
    if [[ "${compression}" == "gzip" ]]; then
      gzip -dc "${abs_path}" | chopper -l "${MIN_LENGTH}" --maxlength "${MAX_LENGTH}" -q "${MIN_QUALITY}" | gzip -c > "${filtered}"
    else
      chopper -l "${MIN_LENGTH}" --maxlength "${MAX_LENGTH}" -q "${MIN_QUALITY}" < "${abs_path}" | gzip -c > "${filtered}"
    fi
    note="length ${MIN_LENGTH}-${MAX_LENGTH}; mean quality >= Q${MIN_QUALITY} with chopper"
  elif [[ "${filter_tool}" == "NanoFilt" ]]; then
    if [[ "${compression}" == "gzip" ]]; then
      gzip -dc "${abs_path}" | NanoFilt -l "${MIN_LENGTH}" --maxlength "${MAX_LENGTH}" -q "${MIN_QUALITY}" | gzip -c > "${filtered}"
    else
      NanoFilt -l "${MIN_LENGTH}" --maxlength "${MAX_LENGTH}" -q "${MIN_QUALITY}" "${abs_path}" | gzip -c > "${filtered}"
    fi
    note="length ${MIN_LENGTH}-${MAX_LENGTH}; mean quality >= Q${MIN_QUALITY} with NanoFilt"
  else
    seqkit seq -m "${MIN_LENGTH}" -M "${MAX_LENGTH}" "${abs_path}" | gzip -c > "${filtered}"
    note="length ${MIN_LENGTH}-${MAX_LENGTH}; no FASTQ quality threshold applied because Chopper/NanoFilt unavailable"
  fi

  seqkit stats -a "${filtered}" > "${sample_dir_filt}/seqkit_stats.txt"
  filt_stats="$(seqkit stats -T -a "${filtered}" | awk 'NR==2 {print $4"\t"$7"\t"$10}')"
  retained="$(printf '%s' "${filt_stats}" | cut -f1)"
  mean_retained="$(printf '%s' "${filt_stats}" | cut -f2)"
  median_retained="$(printf '%s' "${filt_stats}" | cut -f3)"
  percent="$(awk -v r="${retained}" -v raw="${raw_reads}" 'BEGIN { if (raw > 0) printf "%.3f", 100*r/raw; else print "NA" }')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tincluded\n' "${sample}" "${raw_reads}" "${length_retained}" "${retained}" "${percent}" "${mean_retained}" "${median_retained}" "${note}" >> "${qc_table}"
done

printf 'Wrote %s\n' "${qc_table}"

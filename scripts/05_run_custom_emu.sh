#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_project_config
make_log "05_run_custom_emu"

db="${PROJECT_ROOT}/databases/custom_emu/consortium5_custom_emu"
exclusions="${PROJECT_ROOT}/metadata/sample_exclusions.tsv"
[[ -f "${PROJECT_ROOT}/${VALIDATED_CUSTOM_DB_FLAG}" ]] || die "Custom Emu database validation flag missing; run 02_build_custom_emu_database.sh"
[[ -d "${db}" ]] || die "Missing custom Emu database directory: ${db}"
have_cmd emu || die "emu not found on PATH"

emu abundance --help > "${PROJECT_ROOT}/logs/emu_abundance_help.txt" 2>&1 || true
help_text="$(cat "${PROJECT_ROOT}/logs/emu_abundance_help.txt")"
common_args=(--db "${db}" --threads "${THREADS}")
grep -q -- '--type' <<< "${help_text}" && common_args+=(--type map-ont)
grep -q -- '--min-align-len' <<< "${help_text}" && common_args+=(--min-align-len "${MIN_LENGTH}")
grep -q -- '--max-align-len' <<< "${help_text}" && common_args+=(--max-align-len "${MAX_LENGTH}")
grep -q -- '--output-dir' <<< "${help_text}" && has_output_dir=1 || has_output_dir=0
grep -q -- '--keep-read-assignments' <<< "${help_text}" && common_args+=(--keep-read-assignments)
grep -q -- '--keep-counts' <<< "${help_text}" && common_args+=(--keep-counts)
grep -q -- '--output-unclassified' <<< "${help_text}" && common_args+=(--output-unclassified)

mkdir -p "${PROJECT_ROOT}/results/custom_emu" "${PROJECT_ROOT}/results/combined"
if find "${PROJECT_ROOT}/results/custom_emu" -type f -name '*rel-abundance*.tsv' -print -quit | grep -q .; then
  die "Existing custom Emu abundance outputs detected. This script will not overwrite them; run 05_summarize_custom_emu.py to regenerate summaries."
fi
status_table="${PROJECT_ROOT}/results/combined/custom_emu_run_status.tsv"
printf 'sample_id\tusable_reads\tstatus\toutput_directory\n' > "${status_table}"

mapfile -t fastqs < <(find "${PROJECT_ROOT}/fastq/filtered" -type f -name '*.fastq.gz' | sort)
[[ "${#fastqs[@]}" -gt 0 ]] || die "No filtered FASTQ files found; run 04_qc_and_filter.sh first"

for fq in "${fastqs[@]}"; do
  sample="$(basename "${fq}" .full_length_qc.fastq.gz)"
  outdir="${PROJECT_ROOT}/results/custom_emu/${sample}"
  mkdir -p "${outdir}"
  if [[ -f "${exclusions}" ]] && awk -F'\t' -v sample="${sample}" 'NR > 1 && $1 == sample { found=1 } END { exit !found }' "${exclusions}"; then
    warn "Skipping excluded sample ${sample}"
    printf '%s\tNA\texcluded_by_metadata\t%s\n' "${sample}" "${outdir}" >> "${status_table}"
    continue
  fi
  emu --version > "${outdir}/emu_version.txt" 2>&1 || true
  usable_reads="$(seqkit stats -T "${fq}" | awk 'NR==2 {print $4}')"
  usable_reads="${usable_reads:-0}"
  if [[ "${usable_reads}" == "0" ]]; then
    warn "Skipping ${sample}: no reads passed QC"
    printf '%s\t0\tskipped_no_reads\t%s\n' "${sample}" "${outdir}" >> "${status_table}"
    continue
  fi
  cmd=(emu abundance "${fq}" "${common_args[@]}")
  if [[ "${has_output_dir}" == "1" ]]; then
    cmd+=(--output-dir "${outdir}")
    run_dir="${PROJECT_ROOT}"
  else
    run_dir="${outdir}"
  fi
  if grep -q -- '--output-basename' <<< "${help_text}"; then
    cmd+=(--output-basename "${sample}")
  fi
  printf '%q ' "${cmd[@]}" > "${outdir}/emu_command.txt"; printf '\n' >> "${outdir}/emu_command.txt"
  (cd "${run_dir}" && "${cmd[@]}") > "${outdir}/emu.stdout.log" 2> "${outdir}/emu.stderr.log"
  printf '%s\t%s\tcompleted\t%s\n' "${sample}" "${usable_reads}" "${outdir}" >> "${status_table}"
done

python_cmd="python3"
command -v python3 >/dev/null 2>&1 || python_cmd="python"
"${python_cmd}" "${SCRIPT_DIR}/05_summarize_custom_emu.py" --project-root "${PROJECT_ROOT}"

printf 'Custom Emu results written under results/custom_emu and results/combined\n'

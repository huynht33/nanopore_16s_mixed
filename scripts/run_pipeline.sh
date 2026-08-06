#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/run_pipeline.sh --project-root /absolute/path/to/nanopore_16s_mixed [--run-id NAME]

This controlled runner requires reviewed metadata and validated references/database.
It never deletes previous outputs.
EOF
}

PROJECT_ROOT=""
RUN_ID=""
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --run-id) RUN_ID="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage; exit 2 ;;
  esac
done
[[ -n "${PROJECT_ROOT}" ]] || { usage >&2; exit 2; }
export PROJECT_ROOT

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_project_config

RUN_ID="${RUN_ID:-$(timestamp)}"
run_dir="${PROJECT_ROOT}/results/runs/${RUN_ID}"
if [[ -e "${run_dir}" ]]; then
  die "Run directory already exists: ${run_dir}. Choose a new --run-id; previous outputs are never replaced automatically."
fi
mkdir -p "${run_dir}" "${PROJECT_ROOT}/logs"
master_log="${PROJECT_ROOT}/logs/${RUN_ID}_run_pipeline.master.log"
exec >> "${master_log}" 2>&1

printf 'Master log: %s\n' "${master_log}"
printf 'Run directory: %s\n' "${run_dir}"

[[ -f "${PROJECT_ROOT}/${VALIDATED_METADATA_FLAG}" ]] || die "Metadata review flag missing. After editing sample metadata and expected ratios, run: touch ${PROJECT_ROOT}/${VALIDATED_METADATA_FLAG}"
[[ -f "${PROJECT_ROOT}/${VALIDATED_REFERENCES_FLAG}" ]] || die "Reference validation flag missing; run scripts/01_prepare_references.sh"
[[ -f "${PROJECT_ROOT}/${VALIDATED_CUSTOM_DB_FLAG}" ]] || die "Custom database validation flag missing; run scripts/02_build_custom_emu_database.sh"

"${PROJECT_ROOT}/scripts/03_inventory_fastq.sh"
manifest="${PROJECT_ROOT}/metadata/fastq_manifest.tsv"
[[ -s "${manifest}" ]] || die "FASTQ manifest missing or empty"
duplicates="$(tail -n +2 "${manifest}" | cut -f1 | sort | uniq -d)"
[[ -z "${duplicates}" ]] || die "Duplicate sample IDs remain in ${manifest}: ${duplicates}"

"${PROJECT_ROOT}/scripts/04_qc_and_filter.sh"
"${PROJECT_ROOT}/scripts/05_run_custom_emu.sh"
"${PROJECT_ROOT}/scripts/06_run_broad_emu.sh"
"${PROJECT_ROOT}/scripts/07_run_minimap2_crosscheck.sh"
Rscript "${PROJECT_ROOT}/scripts/08_analyze_results.R" "${PROJECT_ROOT}"

printf 'Pipeline complete for run %s\n' "${RUN_ID}"

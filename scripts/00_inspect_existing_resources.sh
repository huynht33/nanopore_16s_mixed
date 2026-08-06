#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_project_config
make_log "00_inspect_existing_resources"

inventory="${PROJECT_ROOT}/reports/existing_resources_inventory.txt"
mkdir -p "${PROJECT_ROOT}/reports" "${PROJECT_ROOT}/databases/broad_emu" "${PROJECT_ROOT}/config"

{
  printf 'Nanopore 16S existing resources inventory\n'
  printf 'Created: %s\n' "$(date -Is)"
  printf 'New project: %s\n' "${PROJECT_ROOT}"
  printf 'Old project inspected read-only: %s\n\n' "${OLD_PROJECT_ROOT}"

  printf '## Host\n'
  hostname || true
  uname -a || true
  printf '\n## Shell\n%s\n\n' "${SHELL:-unknown}"

  printf '## Conda\n'
  if command -v conda >/dev/null 2>&1; then
    conda info || true
    printf '\n### Conda environments\n'
    conda env list || true
  else
    printf 'conda not found on PATH\n'
  fi

  printf '\n## Modules\n'
  if command -v module >/dev/null 2>&1; then
    module avail 2>&1 | sed -n '1,200p' || true
  else
    printf 'module command not available in non-interactive shell\n'
  fi

  printf '\n## Tool discovery\n'
  for tool in emu minimap2 seqkit samtools R Rscript NanoPlot chopper NanoFilt python python3; do
    printf '\n### %s\n' "${tool}"
    if command -v "${tool}" >/dev/null 2>&1; then
      command -v "${tool}"
      case "${tool}" in
        emu)
          emu --help 2>&1 | sed -n '1,80p' || true
          emu abundance --help 2>&1 | sed -n '1,120p' || true
          emu build-database --help 2>&1 | sed -n '1,120p' || true
          ;;
        minimap2) minimap2 --help 2>&1 | sed -n '1,120p' || true ;;
        seqkit) seqkit version 2>&1 || seqkit --help 2>&1 | sed -n '1,60p' || true ;;
        samtools) samtools --version 2>&1 | sed -n '1,40p' || true ;;
        R|Rscript) "${tool}" --version 2>&1 | sed -n '1,20p' || true ;;
        NanoPlot|chopper|NanoFilt) "${tool}" --help 2>&1 | sed -n '1,80p' || true ;;
        python|python3) "${tool}" --version 2>&1 || true ;;
      esac
    else
      printf 'not found\n'
    fi
  done

  printf '\n## Old project listing\n'
  if [[ -d "${OLD_PROJECT_ROOT}" ]]; then
    find "${OLD_PROJECT_ROOT}" -maxdepth 4 -print | sed -n '1,1000p'
  else
    printf 'Old project directory not found: %s\n' "${OLD_PROJECT_ROOT}"
  fi

  printf '\n## Candidate Emu databases\n'
  if [[ -d "${OLD_PROJECT_ROOT}" ]]; then
    find "${OLD_PROJECT_ROOT}" -type f \( -name 'species_taxid.fasta' -o -name 'taxonomy.tsv' \) -print | sort
  fi

  printf '\n## Candidate NCBI taxonomy files\n'
  if [[ -d "${OLD_PROJECT_ROOT}" ]]; then
    find "${OLD_PROJECT_ROOT}" -type f \( -name 'names.dmp' -o -name 'nodes.dmp' \) -print | sort
  fi

  printf '\n## Candidate FASTQ files in old project\n'
  if [[ -d "${OLD_PROJECT_ROOT}" ]]; then
    find "${OLD_PROJECT_ROOT}" -type f \( -iname '*.fastq' -o -iname '*.fastq.gz' -o -iname '*.fq' -o -iname '*.fq.gz' \) -print | sort
  fi
} > "${inventory}"

printf 'Wrote %s\n' "${inventory}"

mapfile -t db_dirs < <(
  if [[ -d "${OLD_PROJECT_ROOT}" ]]; then
    find "${OLD_PROJECT_ROOT}" -type f -name 'species_taxid.fasta' -print |
      while read -r f; do
        d="$(dirname "${f}")"
        [[ -f "${d}/taxonomy.tsv" ]] && printf '%s\n' "${d}"
      done | sort -u
  fi
)

if [[ "${#db_dirs[@]}" -ge 1 ]]; then
  chosen="${db_dirs[0]}"
  link="${PROJECT_ROOT}/databases/broad_emu/emu_broad_db"
  if [[ ! -e "${link}" ]]; then
    ln -s "${chosen}" "${link}"
    printf 'Created broad database symlink: %s -> %s\n' "${link}" "${chosen}"
  else
    printf 'Broad database link/path already exists, leaving unchanged: %s\n' "${link}"
  fi
  {
    printf 'BROAD_EMU_DB="%s"\n' "${link}"
    tax_dir="$(find "${OLD_PROJECT_ROOT}" -type f -name names.dmp -print -quit 2>/dev/null || true)"
    if [[ -n "${tax_dir}" ]]; then
      tax_dir="$(dirname "${tax_dir}")"
      if [[ -f "${tax_dir}/nodes.dmp" ]]; then
        printf 'NCBI_TAXONOMY_DIR="%s"\n' "${tax_dir}"
      fi
    fi
  } > "${PROJECT_ROOT}/config/local_resources.sh"
  printf 'Wrote %s\n' "${PROJECT_ROOT}/config/local_resources.sh"
else
  warn "No valid broad Emu database with species_taxid.fasta and taxonomy.tsv was found under ${OLD_PROJECT_ROOT}"
fi

if command -v conda >/dev/null 2>&1 && [[ -n "${CONDA_PREFIX:-}" ]]; then
  conda env export > "${PROJECT_ROOT}/environments/emu16s_environment.yml" || true
  printf 'Exported active conda environment to environments/emu16s_environment.yml\n'
fi

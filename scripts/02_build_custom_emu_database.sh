#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_project_config
make_log "02_build_custom_emu_database"

combined="${PROJECT_ROOT}/references/combined/consortium5_reference.fasta"
seq2tax="${PROJECT_ROOT}/references/combined/seq2tax.map.tsv"
tax_verify="${PROJECT_ROOT}/references/combined/taxonomy_verification.tsv"
db_parent="${PROJECT_ROOT}/databases/custom_emu"
db_name="consortium5_custom_emu"
db_dir="${db_parent}/${db_name}"
manifest="${db_parent}/custom_emu_database_manifest.tsv"

database_has_expected_files() {
  local candidate="$1"
  [[ -f "${candidate}/species_taxid.fasta" ]] || return 1
  if [[ -f "${candidate}/taxonomy.tsv" ]]; then
    return 0
  fi
  [[ -f "${candidate}/names_df.tsv" && -f "${candidate}/nodes_df.tsv" ]] || return 1
  [[ -f "${candidate}/unique_taxids.tsv" || -f "${candidate}/unqiue_taxids.tsv" ]]
}

[[ -f "${combined}" ]] || die "Missing ${combined}; run 01_prepare_references.sh"
[[ -f "${seq2tax}" ]] || die "Missing ${seq2tax}; run 01_prepare_references.sh with taxonomy enabled"
[[ -f "${tax_verify}" ]] || die "Missing ${tax_verify}"
grep -q $'\tOK$' "${tax_verify}" || die "Taxonomy verification has no OK records"
[[ "$(wc -l < "${seq2tax}" | tr -d ' ')" == "5" ]] || die "seq2tax must contain exactly five headerless mappings"
awk -F'\t' 'NF != 2 || $1 == "sequence_id" || $2 !~ /^[0-9]+$/ { exit 1 }' "${seq2tax}" ||
  die "seq2tax must be a headerless two-column sequence-ID/tax-ID file"
have_cmd emu || die "emu not found on PATH; activate the Emu environment first"

emu --help > "${PROJECT_ROOT}/logs/emu_help.txt" 2>&1 || true
emu abundance --help > "${PROJECT_ROOT}/logs/emu_abundance_help.txt" 2>&1 || true
emu build-database --help > "${PROJECT_ROOT}/logs/emu_build_database_help.txt" 2>&1 || true
emu_version="$(emu --version 2>&1 || emu --help 2>&1 | sed -n '1p')"
printf 'Emu version/help line: %s\n' "${emu_version}"

input_md5="$(md5_file "${combined}")"
map_md5="$(md5_file "${seq2tax}")"
if [[ -f "${PROJECT_ROOT}/${VALIDATED_CUSTOM_DB_FLAG}" && -f "${manifest}" ]]; then
  if grep -q "${input_md5}" "${manifest}" && grep -q "${map_md5}" "${manifest}" && database_has_expected_files "${db_dir}"; then
    printf 'Existing validated custom Emu database matches current reference and seq2tax checksums; not rebuilding.\n'
    exit 0
  fi
fi

if [[ -e "${db_dir}" ]]; then
  die "Custom database path exists but is not a validated checksum match: ${db_dir}. Inspect it manually; this script will not overwrite it."
fi

mkdir -p "${db_parent}"

cmd=(emu build-database "${db_dir}" --sequences "${combined}" --seq2tax "${seq2tax}")
if [[ -n "${NCBI_TAXONOMY_DIR}" && -f "${NCBI_TAXONOMY_DIR}/names.dmp" && -f "${NCBI_TAXONOMY_DIR}/nodes.dmp" ]]; then
  cmd+=(--ncbi-taxonomy "${NCBI_TAXONOMY_DIR}")
elif [[ -n "${TAXONOMY_LIST}" && -f "${TAXONOMY_LIST}" ]]; then
  cmd+=(--taxonomy-list "${TAXONOMY_LIST}")
elif [[ -n "${BROAD_EMU_DB}" && -f "${BROAD_EMU_DB}/taxonomy.tsv" ]]; then
  cmd+=(--taxonomy-list "${BROAD_EMU_DB}/taxonomy.tsv")
elif [[ -f "${PROJECT_ROOT}/databases/broad_emu/emu_broad_db/taxonomy.tsv" ]]; then
  cmd+=(--taxonomy-list "${PROJECT_ROOT}/databases/broad_emu/emu_broad_db/taxonomy.tsv")
else
  die "No NCBI taxonomy directory or taxonomy.tsv available for emu build-database"
fi

printf '%q ' "${cmd[@]}" > "${db_parent}/custom_emu_build_command.txt"
printf '\n' >> "${db_parent}/custom_emu_build_command.txt"
printf 'Running: '
printf '%q ' "${cmd[@]}"
printf '\n'
"${cmd[@]}"

database_has_expected_files "${db_dir}" ||
  die "Built database does not contain species_taxid.fasta plus a supported Emu taxonomy metadata layout"

{
  printf 'item\tpath\tmd5\tbytes\n'
  for f in "${combined}" "${seq2tax}"; do
    printf '%s\t%s\t%s\t%s\n' "$(basename "${f}")" "${f}" "$(md5_file "${f}")" "$(wc -c < "${f}" | tr -d ' ')"
  done
  find "${db_dir}" -maxdepth 1 -type f | sort | while read -r f; do
    printf '%s\t%s\t%s\t%s\n' "$(basename "${f}")" "${f}" "$(md5_file "${f}")" "$(wc -c < "${f}" | tr -d ' ')"
  done
} > "${manifest}"

touch "${PROJECT_ROOT}/${VALIDATED_CUSTOM_DB_FLAG}"
printf 'Validated custom Emu database: %s\n' "${db_dir}"

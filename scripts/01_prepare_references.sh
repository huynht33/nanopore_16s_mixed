#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_project_config
make_log "01_prepare_references"

skip_taxonomy=0
if [[ "${1:-}" == "--skip-taxonomy" ]]; then
  skip_taxonomy=1
fi

orig_dir="${PROJECT_ROOT}/references/original"
clean_dir="${PROJECT_ROOT}/references/clean"
combined_dir="${PROJECT_ROOT}/references/combined"
mkdir -p "${orig_dir}" "${clean_dir}" "${combined_dir}"

summary="${combined_dir}/reference_summary.tsv"
combined="${combined_dir}/consortium5_reference.fasta"
: > "${combined}"
printf 'reference_id\tspecies\tsource_filename\tsequence_length\tgc_percentage\tn_bases\tmd5\n' > "${summary}"

find_source() {
  local prefix="$1" clean_file="$2"
  local src=""
  while IFS= read -r candidate; do
    if [[ -n "${src}" ]]; then
      die "Multiple original FASTA files match prefix '${prefix}'. Resolve ambiguity in ${orig_dir}."
    fi
    src="${candidate}"
  done < <(find "${orig_dir}" -maxdepth 1 -type f \( -name "${prefix}*.fasta" -o -name "${prefix}*.fa" -o -name "${prefix}*.fna" \) | sort)
  if [[ -z "${src}" && -f "${clean_dir}/${clean_file}" ]]; then
    src="${clean_dir}/${clean_file}"
    warn "No original FASTA found for prefix ${prefix}; validating existing clean FASTA ${src}"
  fi
  [[ -n "${src}" ]] || die "Missing reference FASTA beginning with '${prefix}' in ${orig_dir}"
  printf '%s\n' "${src}"
}

clean_one_reference() {
  local source="$1" output="$2" record_id="$3"
  awk -v id="${record_id}" '
    BEGIN { printed=0; seq="" }
    /^>/ {
      if (printed == 1) {
        printf "ERROR: source contains more than one FASTA record\n" > "/dev/stderr"
        exit 2
      }
      printed=1
      next
    }
    {
      gsub(/[[:space:]]/, "", $0)
      seq = seq toupper($0)
    }
    END {
      if (printed == 0 || length(seq) == 0) {
        printf "ERROR: no FASTA sequence found\n" > "/dev/stderr"
        exit 3
      }
      printf ">%s\n", id
      for (i = 1; i <= length(seq); i += 80) print substr(seq, i, 80)
    }
  ' "${source}" > "${output}.tmp"
  mv "${output}.tmp" "${output}"
}

summarize_reference() {
  local file="$1" species="$2" source="$3"
  awk -v species="${species}" -v source="${source}" -v md5="$(md5_file "${file}")" '
    BEGIN { id=""; seq="" }
    /^>/ { id=substr($0,2); next }
    { gsub(/[[:space:]]/, "", $0); seq=seq toupper($0) }
    END {
      len=length(seq); gc=0; n=0
      if (id == "") { printf "Missing FASTA ID in %s\n", FILENAME > "/dev/stderr"; exit 2 }
      if (seq !~ /^[ACGTRYSWKMBDHVN]+$/) { printf "Invalid nucleotide characters in %s\n", FILENAME > "/dev/stderr"; exit 3 }
      for (i=1; i<=len; i++) {
        base=substr(seq,i,1)
        if (base == "G" || base == "C") gc++
        if (base == "N") n++
      }
      if (len < 1200 || len > 1800) {
        printf "WARNING: %s length %d is outside 1200-1800 bp and requires inspection\n", id, len > "/dev/stderr"
      }
      if (n / len > 0.01) {
        printf "WARNING: %s has %.3f fraction N bases and requires inspection\n", id, n / len > "/dev/stderr"
      }
      printf "%s\t%s\t%s\t%d\t%.3f\t%d\t%s\n", id, species, source, len, 100*gc/len, n, md5
    }
  ' "${file}"
}

for spec in "${REFERENCE_SPECS[@]}"; do
  IFS='|' read -r clean_name species prefix clean_file record_id <<< "${spec}"
  source="$(find_source "${prefix}" "${clean_file}")"
  out="${clean_dir}/${clean_file}"
  clean_one_reference "${source}" "${out}" "${record_id}"
  summarize_reference "${out}" "${species}" "$(basename "${source}")" >> "${summary}"
  cat "${out}" >> "${combined}"
done

printf 'Wrote %s\n' "${summary}"
printf 'Wrote %s\n' "${combined}"

if [[ "${skip_taxonomy}" == "1" ]]; then
  warn "Skipping taxonomy resolution because --skip-taxonomy was requested"
  touch "${combined_dir}/.references_sequence_validated_only"
  exit 0
fi

seq2tax="${combined_dir}/seq2tax.map.tsv"
tax_verify="${combined_dir}/taxonomy_verification.tsv"
: > "${seq2tax}"
printf 'species\tproposed_tax_id\tmatch_source\tstatus\n' > "${tax_verify}"

resolve_from_names_dmp() {
  local species="$1" dir="$2"
  awk -v sp="${species}" '
    {
      split($0, a, "\t\\|\t")
      taxid=a[1]; name=a[2]; class=a[4]
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", taxid)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", name)
      if (name == sp && class ~ /scientific name/) print taxid
    }
  ' "${dir}/names.dmp" | sort -u
}

resolve_from_taxonomy_tsv() {
  local species="$1" file="$2"
  awk -F'\t' -v sp="${species}" '
    NR == 1 { next }
    {
      for (i=1; i<=NF; i++) {
        if ($i == sp) {
          print $1
          break
        }
      }
    }
  ' "${file}" | sort -u
}

taxonomy_source=""
if [[ -n "${NCBI_TAXONOMY_DIR}" && -f "${NCBI_TAXONOMY_DIR}/names.dmp" && -f "${NCBI_TAXONOMY_DIR}/nodes.dmp" ]]; then
  taxonomy_source="ncbi:${NCBI_TAXONOMY_DIR}"
elif [[ -n "${TAXONOMY_LIST}" && -f "${TAXONOMY_LIST}" ]]; then
  taxonomy_source="taxonomy-list:${TAXONOMY_LIST}"
elif [[ -n "${BROAD_EMU_DB}" && -f "${BROAD_EMU_DB}/taxonomy.tsv" ]]; then
  taxonomy_source="taxonomy-list:${BROAD_EMU_DB}/taxonomy.tsv"
elif [[ -f "${PROJECT_ROOT}/databases/broad_emu/emu_broad_db/taxonomy.tsv" ]]; then
  taxonomy_source="taxonomy-list:${PROJECT_ROOT}/databases/broad_emu/emu_broad_db/taxonomy.tsv"
else
  die "No taxonomy source found. Run 00_inspect_existing_resources.sh or set NCBI_TAXONOMY_DIR/TAXONOMY_LIST/BROAD_EMU_DB."
fi

printf 'Using taxonomy source: %s\n' "${taxonomy_source}"
for spec in "${REFERENCE_SPECS[@]}"; do
  IFS='|' read -r clean_name species prefix clean_file record_id <<< "${spec}"
  if [[ "${taxonomy_source}" == ncbi:* ]]; then
    mapfile -t taxids < <(resolve_from_names_dmp "${species}" "${taxonomy_source#ncbi:}")
    source_label="names.dmp scientific name"
  else
    mapfile -t taxids < <(resolve_from_taxonomy_tsv "${species}" "${taxonomy_source#taxonomy-list:}")
    source_label="taxonomy.tsv exact field"
  fi
  if [[ "${#taxids[@]}" -ne 1 ]]; then
    printf '%s\t%s\t%s\t%s\n' "${species}" "$(IFS=,; echo "${taxids[*]:-none}")" "${source_label}" "ERROR_exact_match_unavailable_or_ambiguous" >> "${tax_verify}"
    die "Could not resolve an unambiguous exact taxid for '${species}'. Review ${tax_verify}."
  fi
  printf '%s\t%s\n' "${record_id}" "${taxids[0]}" >> "${seq2tax}"
  printf '%s\t%s\t%s\tOK\n' "${species}" "${taxids[0]}" "${source_label}" >> "${tax_verify}"
done

cat "${tax_verify}"
[[ "$(wc -l < "${seq2tax}" | tr -d ' ')" == "5" ]] || die "seq2tax must contain exactly five headerless mappings"
[[ "$(cut -f1 "${seq2tax}" | sort -u | wc -l | tr -d ' ')" == "5" ]] || die "seq2tax sequence IDs are not unique"
[[ "$(cut -f2 "${seq2tax}" | sort -u | wc -l | tr -d ' ')" == "5" ]] || die "seq2tax tax IDs are not unique"
touch "${PROJECT_ROOT}/${VALIDATED_REFERENCES_FLAG}"
printf 'Wrote %s and %s\n' "${seq2tax}" "${tax_verify}"

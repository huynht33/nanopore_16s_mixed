#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_project_config
make_log "03_inventory_fastq"

raw_dir="${PROJECT_ROOT}/fastq/raw"
manifest="${PROJECT_ROOT}/metadata/fastq_manifest.tsv"
list_file="$(fastq_list_file)"
mkdir -p "${PROJECT_ROOT}/metadata"

mapfile -t fastqs < <(find "${raw_dir}" -type f \( -iname '*.fastq' -o -iname '*.fastq.gz' -o -iname '*.fq' -o -iname '*.fq.gz' \) | sort)
printf '%s\n' "${fastqs[@]:-}" | sed '/^$/d' > "${list_file}"

count="${#fastqs[@]}"
if [[ "${count}" -eq 0 ]]; then
  die "No FASTQ files found under ${raw_dir}"
fi
if [[ "${count}" -ne "${EXPECTED_FASTQ_COUNT}" ]]; then
  warn "Expected ${EXPECTED_FASTQ_COUNT} FASTQ files but found ${count}."
  if [[ "${ALLOW_NON_20_FASTQ}" != "1" ]]; then
    die "Review the FASTQ count, then rerun with ALLOW_NON_20_FASTQ=1 if intentional."
  fi
fi

have_cmd seqkit || die "seqkit is required for read-count/read-length inventory"
seqkit version > "${PROJECT_ROOT}/logs/seqkit_version.txt" 2>&1 || true

printf 'sample_id\toriginal_filename\tabsolute_path\tcompression\tfile_size_bytes\tmd5\tread_count\tmean_read_length\tmedian_read_length\tminimum_read_length\tmaximum_read_length\n' > "${manifest}"
tmp_ids="$(mktemp)"

for fq in "${fastqs[@]}"; do
  base="$(basename "${fq}")"
  sample="${base}"
  sample="${sample%.fastq.gz}"; sample="${sample%.fq.gz}"; sample="${sample%.fastq}"; sample="${sample%.fq}"
  printf '%s\n' "${sample}" >> "${tmp_ids}"
  compression="none"
  [[ "${base}" == *.gz ]] && compression="gzip"
  stats="$(seqkit stats -T -a "${fq}" | awk 'NR==2 {print $4"\t"$7"\t"$10"\t"$6"\t"$8}')"
  read_count="$(printf '%s' "${stats}" | cut -f1)"
  mean_len="$(printf '%s' "${stats}" | cut -f2)"
  median_len="$(printf '%s' "${stats}" | cut -f3)"
  min_len="$(printf '%s' "${stats}" | cut -f4)"
  max_len="$(printf '%s' "${stats}" | cut -f5)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${sample}" "${base}" "$(abs_path "${fq}")" "${compression}" "$(wc -c < "${fq}" | tr -d ' ')" "$(md5_file "${fq}")" \
    "${read_count}" "${mean_len}" "${median_len}" "${min_len}" "${max_len}" >> "${manifest}"
done

duplicates="$(sort "${tmp_ids}" | uniq -d)"
rm -f "${tmp_ids}"
if [[ -n "${duplicates}" ]]; then
  printf 'Duplicate sample IDs:\n%s\n' "${duplicates}" >&2
  die "Two or more FASTQ files resolve to the same sample_id. Rename or edit manifest after review."
fi

printf 'Wrote %s\n' "${manifest}"

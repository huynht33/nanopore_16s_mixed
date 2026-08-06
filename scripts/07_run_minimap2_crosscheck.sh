#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_project_config
make_log "07_run_minimap2_crosscheck"

ref="${PROJECT_ROOT}/references/combined/consortium5_reference.fasta"
[[ -f "${ref}" ]] || die "Missing ${ref}; run 01_prepare_references.sh"
have_cmd minimap2 || die "minimap2 not found on PATH"
have_cmd gzip || die "gzip not found on PATH"
minimap2 --help > "${PROJECT_ROOT}/logs/minimap2_help.txt" 2>&1 || true
minimap2 --version > "${PROJECT_ROOT}/logs/minimap2_version.txt" 2>&1

mkdir -p "${PROJECT_ROOT}/results/minimap2" "${PROJECT_ROOT}/results/combined"
fastqs=()
while IFS= read -r fq; do
  fastqs+=("${fq}")
done < <(find "${PROJECT_ROOT}/fastq/filtered" -type f -name '*.fastq.gz' | sort)
[[ "${#fastqs[@]}" -gt 0 ]] || die "No filtered FASTQ files found; run 04_qc_and_filter.sh first"

combined_assignment="${PROJECT_ROOT}/results/combined/minimap2_assignment_summary.tsv"
combined_status="${PROJECT_ROOT}/results/combined/minimap2_read_status_summary.tsv"
combined_confident="${PROJECT_ROOT}/results/combined/minimap2_confident_reference_abundance.tsv"

for output in "${combined_assignment}" "${combined_status}" "${combined_confident}"; do
  [[ ! -e "${output}" ]] || die "Output already exists; refusing to overwrite: ${output}"
done

for fq in "${fastqs[@]}"; do
  sample="$(basename "${fq}" .full_length_qc.fastq.gz)"
  outdir="${PROJECT_ROOT}/results/minimap2/${sample}"
  [[ ! -e "${outdir}" ]] || die "Output directory already exists; refusing to overwrite: ${outdir}"
  mkdir -p "${outdir}"
  paf="${outdir}/${sample}.consortium5.map_ont.paf"
  read_lengths="${outdir}/${sample}.input_read_lengths.tsv"
  assignments="${outdir}/${sample}.per_read_assignments.tsv"

  gzip -cd "${fq}" |
    awk '
      BEGIN { OFS="\t"; print "read_id", "read_length" }
      NR % 4 == 1 {
        header=substr($0, 2)
        split(header, fields, /[[:space:]]+/)
        read_id=fields[1]
        next
      }
      NR % 4 == 2 {
        if (read_id == "") {
          print "FASTQ sequence encountered without a read header" > "/dev/stderr"
          exit 1
        }
        print read_id, length($0)
        read_id=""
      }
    ' > "${read_lengths}"

  input_reads="$(awk 'END { print NR - 1 }' "${read_lengths}")"
  [[ "${input_reads}" -gt 0 ]] || die "No reads found in ${fq}"
  duplicate_read_id="$(awk -F'\t' 'NR > 1 { n[$1]++ } n[$1] == 2 { print $1; exit }' "${read_lengths}")"
  [[ -z "${duplicate_read_id}" ]] ||
    die "Duplicate read ID in ${fq}: ${duplicate_read_id}"

  cmd=(
    minimap2
    -t "${THREADS}"
    -x map-ont
    -c
    --secondary=yes
    -N 10
    "${ref}"
    "${fq}"
  )
  printf '%q ' "${cmd[@]}" > "${outdir}/minimap2_command.txt"; printf '\n' >> "${outdir}/minimap2_command.txt"
  "${cmd[@]}" > "${paf}" 2> "${outdir}/minimap2.stderr.log"

  awk \
    -v mincov="${MINIMAP_MIN_QUERY_COVERAGE}" \
    -v minid="${MINIMAP_MIN_IDENTITY}" \
    -v minmapq="${MINIMAP_MIN_MAPQ}" \
    -v mindelta="${MINIMAP_MIN_AS_DELTA}" \
    -v minaln="${MINIMAP_MIN_ALIGNMENT_LENGTH}" \
    -v maxaln="${MINIMAP_MAX_ALIGNMENT_LENGTH}" '
    function tagval(prefix,   i) {
      for (i=13; i<=NF; i++) if ($i ~ "^" prefix) { split($i,a,":"); return a[3] }
      return ""
    }
    BEGIN { FS=OFS="\t" }
    FNR == NR {
      if (FNR == 1) next
      read_order[++read_count]=$1
      read_length[$1]=$2
      next
    }
    {
      q=$1
      qlen=$2
      query_span=$4-$3
      ref=$6
      alignment_length=$11
      mapq=$12
      coverage=(qlen>0 ? query_span/qlen : 0)
      identity=(alignment_length>0 ? $10/alignment_length : 0)
      alignment_score=tagval("AS:i:")
      if (alignment_score == "") alignment_score=$10

      key=q SUBSEP ref
      if (!(key in reference_score) ||
          alignment_score+0 > reference_score[key]+0) {
        reference_score[key]=alignment_score
        reference_alignment_length[key]=alignment_length
        reference_coverage[key]=coverage
        reference_identity[key]=identity
        reference_mapq[key]=mapq
      }
    }
    END {
      for (key in reference_score) {
        split(key, pair, SUBSEP)
        q=pair[1]
        score=reference_score[key]+0
        if (!(q in best_key) || score > best_score[q]+0) {
          if (q in best_key) {
            second_key[q]=best_key[q]
            second_score[q]=best_score[q]
          }
          best_key[q]=key
          best_score[q]=score
        } else if (!(q in second_key) || score > second_score[q]+0) {
          second_key[q]=key
          second_score[q]=score
        }
      }

      print "read_id\tread_length\tbest_reference\tbest_alignment_length\tquery_coverage\tsequence_identity\tmapping_quality\tsecond_best_reference\tbest_minus_second_alignment_score\tassignment_status"

      for (i=1; i<=read_count; i++) {
        q=read_order[i]
        if (!(q in best_key)) {
          print q, read_length[q], "NA", "NA", "NA", "NA", "NA", "NA", "NA", "unassigned"
          continue
        }

        best=best_key[q]
        split(best, best_pair, SUBSEP)
        best_ref=best_pair[2]
        second_ref="NA"
        delta="NA"
        if (q in second_key) {
          split(second_key[q], second_pair, SUBSEP)
          second_ref=second_pair[2]
          delta=best_score[q]-second_score[q]
        }

        status="confident"
        if (reference_alignment_length[best]+0 < minaln ||
            reference_alignment_length[best]+0 > maxaln ||
            reference_coverage[best]+0 < mincov ||
            reference_identity[best]+0 < minid) {
          status="unassigned"
        } else if (second_ref != "NA" && delta+0 < mindelta) {
          status="ambiguous"
        } else if (reference_mapq[best]+0 < minmapq) {
          status="unassigned"
        }

        printf "%s\t%s\t%s\t%s\t%.5f\t%.5f\t%s\t%s\t%s\t%s\n",
          q,
          read_length[q],
          best_ref,
          reference_alignment_length[best],
          reference_coverage[best],
          reference_identity[best],
          reference_mapq[best],
          second_ref,
          delta,
          status
      }
    }
  ' "${read_lengths}" "${paf}" > "${assignments}"

  assigned_rows="$(awk 'END { print NR - 1 }' "${assignments}")"
  [[ "${assigned_rows}" -eq "${input_reads}" ]] ||
    die "Per-read output count mismatch for ${sample}: ${assigned_rows} versus ${input_reads}"

  awk -F'\t' '
    NR > 1 && $10 != "confident" && $10 != "ambiguous" && $10 != "unassigned" {
      print $10
      exit
    }
  ' "${assignments}" |
    while read -r invalid_status; do
      [[ -z "${invalid_status}" ]] ||
        die "Unexpected assignment status for ${sample}: ${invalid_status}"
    done

  printf '\n%s\n' "${sample}"
  awk -F'\t' '
    NR > 1 { status[$10]++; total++ }
    END {
      printf "  total=%d confident=%d ambiguous=%d unassigned=%d\n",
        total, status["confident"], status["ambiguous"], status["unassigned"]
    }
  ' "${assignments}"
done

{
  printf 'sample_id\treference\tassignment_status\treads\tpercentage_all_input_reads\n'
  find "${PROJECT_ROOT}/results/minimap2" -name '*.per_read_assignments.tsv' | sort | while read -r tbl; do
    sample="$(basename "$(dirname "${tbl}")")"
    awk -F'\t' -v sample="${sample}" '
      BEGIN { OFS="\t" }
      NR > 1 { key=$3 OFS $10; n[key]++; total++ }
      END {
        for (key in n) {
          print sample, key, n[key], 100*n[key]/total
        }
      }
    ' "${tbl}"
  done
} > "${combined_assignment}"

{
  printf 'sample_id\ttotal_reads\tconfident_reads\tambiguous_reads\tunassigned_reads\tconfident_pct\tambiguous_pct\tunassigned_pct\n'
  find "${PROJECT_ROOT}/results/minimap2" -name '*.per_read_assignments.tsv' | sort | while read -r tbl; do
    sample="$(basename "$(dirname "${tbl}")")"
    awk -F'\t' -v sample="${sample}" '
      BEGIN { OFS="\t" }
      NR > 1 { status[$10]++; total++ }
      END {
        print sample,
          total,
          status["confident"]+0,
          status["ambiguous"]+0,
          status["unassigned"]+0,
          100*(status["confident"]+0)/total,
          100*(status["ambiguous"]+0)/total,
          100*(status["unassigned"]+0)/total
      }
    ' "${tbl}"
  done
} > "${combined_status}"

{
  printf 'sample_id\treference\tconfident_reads\tpercentage_all_input_reads\tpercentage_among_confident_reads\n'
  find "${PROJECT_ROOT}/results/minimap2" -name '*.per_read_assignments.tsv' | sort | while read -r tbl; do
    sample="$(basename "$(dirname "${tbl}")")"
    awk -F'\t' -v sample="${sample}" '
      BEGIN { OFS="\t" }
      NR > 1 {
        total++
        if ($10 == "confident") {
          reference[$3]++
          confident++
        }
      }
      END {
        for (ref in reference) {
          print sample,
            ref,
            reference[ref],
            100*reference[ref]/total,
            100*reference[ref]/confident
        }
      }
    ' "${tbl}"
  done
} > "${combined_confident}"

printf '\nMinimap2 cross-check outputs written under results/minimap2\n'
printf 'Combined assignment summary: %s\n' "${combined_assignment}"
printf 'Combined read-status summary: %s\n' "${combined_status}"
printf 'Combined confident-reference summary: %s\n' "${combined_confident}"

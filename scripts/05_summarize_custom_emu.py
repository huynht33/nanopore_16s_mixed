#!/usr/bin/env python3
import argparse
import csv
import glob
import os


SPECIES = [
    ("Bacteroides thetaiotaomicron", "Bacteroides_thetaiotaomicron"),
    ("Cutibacterium acnes", "Cutibacterium_acnes"),
    ("Fusobacterium nucleatum", "Fusobacterium_nucleatum"),
    ("Prevotella melaninogenica", "Prevotella_melaninogenica"),
    ("Veillonella parvula", "Veillonella_parvula"),
]


def number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def formatted(value):
    return f"{value:.6f}"


def write_table(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize existing five-species custom Emu outputs without rerunning Emu."
    )
    parser.add_argument("--project-root", required=True)
    args = parser.parse_args()
    root = os.path.abspath(args.project_root)

    qc_path = os.path.join(root, "qc/filtered/filter_summary.tsv")
    with open(qc_path, newline="") as handle:
        qc_rows = list(csv.DictReader(handle, delimiter="\t"))
    usable_by_sample = {
        row["sample_id"]: int(number(row["reads_after_quality_filter"]))
        for row in qc_rows
        if row.get("analysis_status", "included") == "included"
    }

    rows_all = []
    rows_target = []
    accounting = []
    wide_rows = []

    for sample, usable_reads in usable_by_sample.items():
        sample_dir = os.path.join(root, "results/custom_emu", sample)
        files = sorted(
            glob.glob(os.path.join(sample_dir, "**", "*rel-abundance*.tsv"), recursive=True)
        )
        preferred = [
            path for path in files if "threshold" not in os.path.basename(path).lower()
        ]
        path = preferred[0] if preferred else (files[0] if files else "")
        if not path:
            for display, _ in SPECIES:
                rows_all.append([sample, display, "0.000000", "0.000000", ""])
                rows_target.append(
                    [sample, display, "0.000000", "0.000000", "0.000000", ""]
                )
            accounting.append(
                [
                    sample,
                    usable_reads,
                    "0.000000",
                    "0.000000",
                    "0.000000",
                    formatted(usable_reads),
                    formatted(usable_reads),
                    formatted(100.0 if usable_reads else 0.0),
                    "no_emu_output",
                ]
            )
            wide_rows.append([sample] + ["0.000000"] * len(SPECIES))
            continue

        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            data = list(reader)
            fields = reader.fieldnames or []

        count_col = next((field for field in fields if "count" in field.lower()), None)
        if count_col is None:
            raise RuntimeError(
                f"No estimated-count column found in {path}; Emu must be run with --keep-counts."
            )

        classified_count = sum(number(row.get(count_col)) for row in data)
        values = {}
        for display, key in SPECIES:
            match = next(
                (
                    row
                    for row in data
                    if display in "\t".join(str(value) for value in row.values())
                    or key in "\t".join(str(value) for value in row.values())
                ),
                None,
            )
            estimated_count = number(match.get(count_col)) if match else 0.0
            all_usable_pct = (
                100.0 * estimated_count / usable_reads if usable_reads else 0.0
            )
            values[display] = (estimated_count, all_usable_pct)

        expected_count = sum(value[0] for value in values.values())
        wide_row = [sample]
        for display, (estimated_count, all_usable_pct) in values.items():
            renormalized_pct = (
                100.0 * estimated_count / expected_count if expected_count else 0.0
            )
            rows_all.append(
                [
                    sample,
                    display,
                    formatted(all_usable_pct),
                    formatted(estimated_count),
                    path,
                ]
            )
            rows_target.append(
                [
                    sample,
                    display,
                    formatted(renormalized_pct),
                    formatted(all_usable_pct),
                    formatted(estimated_count),
                    path,
                ]
            )
            wide_row.append(formatted(renormalized_pct))
        wide_rows.append(wide_row)

        emu_other = max(0.0, classified_count - expected_count)
        unmapped_or_filtered = max(0.0, usable_reads - classified_count)
        total_unaccounted = max(0.0, usable_reads - expected_count)
        total_unaccounted_pct = (
            100.0 * total_unaccounted / usable_reads if usable_reads else 0.0
        )
        accounting.append(
            [
                sample,
                usable_reads,
                formatted(classified_count),
                formatted(expected_count),
                formatted(emu_other),
                formatted(unmapped_or_filtered),
                formatted(total_unaccounted),
                formatted(total_unaccounted_pct),
                "completed",
            ]
        )

    out = os.path.join(root, "results/combined")
    os.makedirs(out, exist_ok=True)
    write_table(
        os.path.join(out, "custom_emu_expected_species_abundance_all_usable_reads.tsv"),
        [
            "sample_id",
            "species",
            "emu_abundance_pct_all_usable_reads",
            "estimated_counts",
            "source_file",
        ],
        rows_all,
    )
    write_table(
        os.path.join(out, "custom_emu_expected_species_abundance_renormalized.tsv"),
        [
            "sample_id",
            "species",
            "emu_abundance_pct_renormalized_among_expected",
            "emu_abundance_pct_all_usable_reads",
            "estimated_counts",
            "source_file",
        ],
        rows_target,
    )
    write_table(
        os.path.join(out, "custom_emu_expected_species_abundance_renormalized_wide.tsv"),
        ["sample_id"] + [display for display, _ in SPECIES],
        wide_rows,
    )
    write_table(
        os.path.join(out, "custom_emu_read_accounting.tsv"),
        [
            "sample_id",
            "usable_reads",
            "emu_classified_estimated_reads",
            "expected_species_estimated_reads",
            "emu_other_or_unclassified_estimated_reads",
            "unmapped_or_filtered_reads",
            "total_unassigned_or_unclassified_reads",
            "total_unassigned_or_unclassified_pct",
            "status",
        ],
        accounting,
    )


if __name__ == "__main__":
    main()

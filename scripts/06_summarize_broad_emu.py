#!/usr/bin/env python3
import argparse
import csv
import glob
import os


EXPECTED_SPECIES = {
    "Bacteroides thetaiotaomicron",
    "Cutibacterium acnes",
    "Fusobacterium nucleatum",
    "Prevotella melaninogenica",
    "Veillonella parvula",
}
EXPECTED_GENERA = {
    "Bacteroides",
    "Cutibacterium",
    "Fusobacterium",
    "Prevotella",
    "Veillonella",
}
READ_STATUS_IDS = {"mapped_filtered", "mapped_unclassified", "unmapped"}


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


def classify_taxon(species, genus):
    if species in EXPECTED_SPECIES:
        return "expected_species"
    text = f"{species} {genus}".lower()
    if not species or species.lower() in {"na", "nan", "unclassified"} or "unclassified" in text:
        return "unclassified"
    if genus in EXPECTED_GENERA:
        return "related_expected_genus"
    return "unexpected_taxon"


def main():
    parser = argparse.ArgumentParser(
        description="Summarize existing broad Emu outputs at species and genus levels."
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

    species_rows = []
    unexpected_rows = []
    genus_rows = []
    accounting_rows = []
    read_status_rows = []

    for sample, usable_reads in usable_by_sample.items():
        sample_dir = os.path.join(root, "results/broad_emu", sample)
        files = sorted(
            glob.glob(os.path.join(sample_dir, "**", "*rel-abundance*.tsv"), recursive=True)
        )
        preferred = [
            path for path in files if "threshold" not in os.path.basename(path).lower()
        ]
        path = preferred[0] if preferred else (files[0] if files else "")
        if not path:
            accounting_rows.append(
                [
                    sample,
                    usable_reads,
                    "0.000000",
                    "0.000000",
                    "0.000000",
                    "0.000000",
                    "0.000000",
                    formatted(usable_reads),
                    "0.000000",
                    formatted(100.0 if usable_reads else 0.0),
                    "no_emu_output",
                ]
            )
            continue

        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            data = list(reader)
            fields = reader.fieldnames or []

        count_col = next((field for field in fields if "count" in field.lower()), None)
        abundance_col = next(
            (
                field
                for field in fields
                if "abundance" in field.lower() or "percent" in field.lower()
            ),
            None,
        )
        species_col = next((field for field in fields if field.lower() == "species"), None)
        genus_col = next((field for field in fields if field.lower() == "genus"), None)
        taxid_col = next(
            (
                field
                for field in fields
                if field.lower().replace("_", " ") in {"tax id", "taxid"}
            ),
            None,
        )
        if count_col is None:
            raise RuntimeError(
                f"No estimated-count column found in {path}; Emu must be run with --keep-counts."
            )

        abundance_total = (
            sum(number(row.get(abundance_col)) for row in data) if abundance_col else 0.0
        )
        abundance_scale = 100.0 if abundance_total <= 1.01 else 1.0
        genus_counts = {}
        classified_total = 0.0
        expected_total = 0.0
        unexpected_total = 0.0
        status_counts = {
            "mapped_filtered": 0.0,
            "mapped_unclassified": 0.0,
            "unmapped": 0.0,
            "taxonomy_unclassified": 0.0,
        }

        for row in data:
            species = (row.get(species_col) or "").strip() if species_col else ""
            genus = (row.get(genus_col) or "").strip() if genus_col else ""
            tax_id = (row.get(taxid_col) or "").strip() if taxid_col else ""
            count = number(row.get(count_col))
            classified_pct = (
                number(row.get(abundance_col)) * abundance_scale if abundance_col else 0.0
            )
            all_usable_pct = 100.0 * count / usable_reads if usable_reads else 0.0
            category = classify_taxon(species, genus)
            if tax_id in READ_STATUS_IDS:
                status_counts[tax_id] += count
                read_status_rows.append(
                    [
                        sample,
                        tax_id,
                        formatted(count),
                        formatted(all_usable_pct),
                    ]
                )
                continue
            if category == "unclassified":
                status_counts["taxonomy_unclassified"] += count
                read_status_rows.append(
                    [
                        sample,
                        "taxonomy_unclassified",
                        formatted(count),
                        formatted(all_usable_pct),
                    ]
                )
                continue

            classified_total += count
            if category == "expected_species":
                expected_total += count
            else:
                unexpected_total += count

            output_row = [
                sample,
                tax_id,
                species,
                genus,
                category,
                formatted(count),
                formatted(classified_pct),
                formatted(all_usable_pct),
                path,
            ]
            species_rows.append(output_row)
            if category in {"related_expected_genus", "unexpected_taxon"}:
                unexpected_rows.append(output_row)
            genus_name = genus or "unclassified"
            genus_counts[genus_name] = genus_counts.get(genus_name, 0.0) + count

        for genus, count in sorted(genus_counts.items()):
            genus_rows.append(
                [
                    sample,
                    genus,
                    "expected_genus" if genus in EXPECTED_GENERA else "unexpected_genus",
                    formatted(count),
                    formatted(100.0 * count / usable_reads if usable_reads else 0.0),
                ]
            )

        technical_unassigned = sum(status_counts.values())
        total_non_expected = unexpected_total + technical_unassigned
        accounting_rows.append(
            [
                sample,
                usable_reads,
                formatted(classified_total),
                formatted(expected_total),
                formatted(unexpected_total),
                formatted(status_counts["mapped_unclassified"]),
                formatted(status_counts["mapped_filtered"]),
                formatted(status_counts["unmapped"]),
                formatted(status_counts["taxonomy_unclassified"]),
                formatted(100.0 * total_non_expected / usable_reads if usable_reads else 0.0),
                "completed",
            ]
        )

    unexpected_rows.sort(key=lambda row: (row[0], -number(row[7]), row[2]))
    out = os.path.join(root, "results/combined")
    os.makedirs(out, exist_ok=True)
    species_header = [
        "sample_id",
        "tax_id",
        "species",
        "genus",
        "taxon_category",
        "estimated_counts",
        "emu_abundance_pct_among_classified",
        "abundance_pct_all_usable_reads",
        "source_file",
    ]
    write_table(
        os.path.join(out, "broad_emu_species_abundance.tsv"),
        species_header,
        species_rows,
    )
    write_table(
        os.path.join(out, "broad_emu_unexpected_species.tsv"),
        species_header,
        unexpected_rows,
    )
    write_table(
        os.path.join(out, "broad_emu_genus_abundance.tsv"),
        [
            "sample_id",
            "genus",
            "genus_category",
            "estimated_counts",
            "abundance_pct_all_usable_reads",
        ],
        genus_rows,
    )
    write_table(
        os.path.join(out, "broad_emu_read_status.tsv"),
        [
            "sample_id",
            "read_status",
            "estimated_counts",
            "abundance_pct_all_usable_reads",
        ],
        read_status_rows,
    )
    write_table(
        os.path.join(out, "broad_emu_read_accounting.tsv"),
        [
            "sample_id",
            "usable_reads",
            "emu_classified_estimated_reads",
            "expected_species_estimated_reads",
            "unexpected_species_estimated_reads",
            "mapped_unclassified_reads",
            "mapped_filtered_reads",
            "unmapped_reads",
            "taxonomy_unclassified_reads",
            "total_non_expected_pct_all_usable_reads",
            "status",
        ],
        accounting_rows,
    )


if __name__ == "__main__":
    main()

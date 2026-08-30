#!/usr/bin/env python3
"""Create a self-contained, shareable CSV from the normalized gene-set catalog."""

import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--catalog-dir",
        type=Path,
        default=Path("resources/gene-sets/processed/2026-08-29"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("resources/gene-sets/share/2026-08-29/all_gene_sets.csv"),
    )
    parser.add_argument(
        "--wide-output",
        type=Path,
        default=Path("resources/gene-sets/share/2026-08-29/gene_sets_wide.csv"),
    )
    args = parser.parse_args()

    memberships = read_tsv(args.catalog_dir / "gene_set_membership.tsv")
    sources = {row["source_id"]: row for row in read_tsv(args.catalog_dir / "gene_set_sources.tsv")}
    summaries = {row["gene_set_id"]: row for row in read_tsv(args.catalog_dir / "gene_set_summary.tsv")}

    fields = [
        "catalog_release",
        "gene_set_id",
        "gene_set_size",
        "phenotype",
        "definition",
        "gene_symbol",
        "ensembl_gene_id",
        "evidence",
        "effect_direction",
        "effect_beta",
        "p_value",
        "source_id",
        "source_title",
        "source_release",
        "source_url",
        "source_retrieved_date",
        "source_usage_note",
        "membership_notes",
    ]
    output_rows = []
    for membership in memberships:
        source = sources[membership["source_id"]]
        summary = summaries[membership["gene_set_id"]]
        output_rows.append(
            {
                "catalog_release": "2026-08-29",
                "gene_set_id": membership["gene_set_id"],
                "gene_set_size": int(summary["n_genes"]),
                "phenotype": membership["phenotype"],
                "definition": membership["definition"],
                "gene_symbol": membership["gene_symbol"],
                "ensembl_gene_id": membership["ensembl_gene_id"],
                "evidence": membership["evidence"],
                "effect_direction": membership["effect_direction"],
                "effect_beta": membership["effect_beta"],
                "p_value": membership["p_value"],
                "source_id": membership["source_id"],
                "source_title": source["title"],
                "source_release": source["release"],
                "source_url": source["url"],
                "source_retrieved_date": source["retrieved_date"],
                "source_usage_note": source["usage_note"],
                "membership_notes": membership["notes"],
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    gene_sets = {}
    for membership in memberships:
        gene_sets.setdefault(membership["gene_set_id"], set()).add(membership["gene_symbol"])
    set_ids = sorted(gene_sets)
    sorted_genes = {set_id: sorted(gene_sets[set_id]) for set_id in set_ids}
    max_size = max(map(len, sorted_genes.values()))
    args.wide_output.parent.mkdir(parents=True, exist_ok=True)
    with args.wide_output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(set_ids)
        for index in range(max_size):
            writer.writerow([
                sorted_genes[set_id][index] if index < len(sorted_genes[set_id]) else ""
                for set_id in set_ids
            ])

    print(f"Wrote {len(output_rows):,} memberships to {args.output}")
    print(f"Wrote {len(set_ids)} gene-set columns to {args.wide_output}")


if __name__ == "__main__":
    main()

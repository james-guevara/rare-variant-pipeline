#!/usr/bin/env python3
"""Select synonymous alleles in GeneBayes LoF-tier genes."""

import argparse
import csv
from pathlib import Path

from join_genebayes_lof_tiers import GENEBAYES_FIELDS, load_genebayes, stable_id


def gene_tier(post_mean):
    if post_mean == "":
        return ""
    value = float(post_mean)
    if value >= 0.18:
        return "lof_t1"
    if value >= 0.03:
        return "lof_t2"
    return ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--picked", required=True, type=Path)
    parser.add_argument("--genebayes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    constraints = load_genebayes(args.genebayes)
    written = 0
    with args.picked.open(newline="") as source, args.output.open("w", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or ()) + [
            "genebayes_" + field for field in GENEBAYES_FIELDS
        ] + ["lof_tier"]
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            if "synonymous_variant" not in set(row.get("Consequence", "").split("&")):
                continue
            constraint = constraints.get(stable_id(row.get("Gene", "")))
            for field in GENEBAYES_FIELDS:
                row["genebayes_" + field] = constraint.get(field, "") if constraint else ""
            row["lof_tier"] = gene_tier(row["genebayes_post_mean"])
            if row["lof_tier"]:
                writer.writerow(row)
                written += 1
    print(f"synonymous_tiered={written:,}")


if __name__ == "__main__":
    main()

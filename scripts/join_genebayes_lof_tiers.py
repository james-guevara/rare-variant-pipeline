#!/usr/bin/env python3
"""Join GeneBayes constraint and assign canonical/legacy HC pLoF tiers."""

import argparse
import csv
from collections import Counter
from pathlib import Path


GENEBAYES_FIELDS = (
    "obs_lof", "exp_lof", "prior_mean", "post_mean",
    "post_lower_95", "post_upper_95",
)


def stable_id(value):
    return value.split(".", 1)[0]


def load_genebayes(path):
    with path.open(newline="") as handle:
        return {
            stable_id(row["ensg"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def assign_tier(lof, post_mean, t1_min=0.18, t2_min=0.03):
    if lof != "HC" or post_mean == "":
        return ""
    value = float(post_mean)
    if value >= t1_min:
        return "lof_t1"
    if value >= t2_min:
        return "lof_t2"
    return ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--genebayes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path,
                        help="All input pLoFs with GeneBayes annotations")
    parser.add_argument("--qualifying-output", required=True, type=Path,
                        help="HC pLoFs in either GeneBayes tier")
    parser.add_argument("--hc-output", required=True, type=Path,
                        help="all HC pLoFs, including genes outside GeneBayes tiers")
    parser.add_argument("--lof-t1-min", default=0.18, type=float)
    parser.add_argument("--lof-t2-min", default=0.03, type=float)
    args = parser.parse_args()
    if not 0 <= args.lof_t2_min < args.lof_t1_min <= 1:
        parser.error("thresholds must satisfy 0 <= lof-t2-min < lof-t1-min <= 1")

    genebayes = load_genebayes(args.genebayes)
    counts = Counter()
    rows = 0
    matched = 0
    qualifying = 0

    with (
        args.input.open(newline="") as source,
        args.output.open("w", newline="") as all_target,
        args.qualifying_output.open("w", newline="") as qualifying_target,
        args.hc_output.open("w", newline="") as hc_target,
    ):
        reader = csv.DictReader(source, delimiter="\t")
        input_fields = list(reader.fieldnames or ())
        fields = input_fields + [
            "genebayes_" + field for field in GENEBAYES_FIELDS
        ] + ["lof_tier"]
        all_writer = csv.DictWriter(
            all_target, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        qualifying_writer = csv.DictWriter(
            qualifying_target, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        all_writer.writeheader()
        qualifying_writer.writeheader()
        hc_writer = csv.DictWriter(
            hc_target, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        hc_writer.writeheader()

        for row in reader:
            if row.get("LoF", "") not in ("HC", "LC"):
                continue
            rows += 1
            gene = stable_id(row.get("Gene", ""))
            constraint = genebayes.get(gene)
            if constraint is not None:
                matched += 1
            for field in GENEBAYES_FIELDS:
                row["genebayes_" + field] = (
                    constraint.get(field, "") if constraint is not None else ""
                )
            tier = assign_tier(
                row.get("LoF", ""), row["genebayes_post_mean"],
                args.lof_t1_min, args.lof_t2_min,
            )
            row["lof_tier"] = tier
            all_writer.writerow(row)
            if row.get("LoF") == "HC":
                hc_writer.writerow(row)
            counts[tier or "untiered"] += 1
            if tier:
                qualifying_writer.writerow(row)
                qualifying += 1

    print("plof_rows={:,}".format(rows))
    print("genebayes_matched={:,}/{:,}".format(matched, rows))
    for label in ("lof_t1", "lof_t2", "untiered"):
        print("{}={:,}".format(label, counts[label]))
    print("qualifying={:,}".format(qualifying))


if __name__ == "__main__":
    main()

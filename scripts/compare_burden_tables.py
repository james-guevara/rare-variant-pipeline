#!/usr/bin/env python3
"""Compare participant-level burden counts from two pipeline runs."""

import argparse
import csv
from pathlib import Path


DEFAULT_TIERS = "lof_t1,lof_t2,miss_t1,miss_t2,miss_t3,miss_t4"


def read_by_id(path: Path, id_column: str):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or id_column not in rows[0]:
        raise ValueError(f"{path} does not contain {id_column}")
    indexed = {row[id_column]: row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"{path} contains duplicate {id_column} values")
    return indexed


def write_tsv(path: Path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", required=True, type=Path)
    parser.add_argument("--new", required=True, type=Path)
    parser.add_argument("--old-id", default="SAMPLE")
    parser.add_argument("--new-id", default="IID")
    parser.add_argument("--tiers", default=DEFAULT_TIERS)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    old = read_by_id(args.old, args.old_id)
    new = read_by_id(args.new, args.new_id)
    tiers = [value for value in args.tiers.split(",") if value]
    shared = sorted(set(old) & set(new))
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary = []
    discrepancies = []
    for tier in tiers:
        old_total = new_total = exact = 0
        for sample in shared:
            old_value = int(old[sample][tier])
            new_value = int(new[sample][tier])
            old_total += old_value
            new_total += new_value
            exact += old_value == new_value
            if old_value != new_value:
                discrepancies.append({
                    "IID": sample,
                    "tier": tier,
                    "old_count": old_value,
                    "new_count": new_value,
                    "difference": new_value - old_value,
                })
        summary.append({
            "tier": tier,
            "old_total": old_total,
            "new_total": new_total,
            "total_difference": new_total - old_total,
            "samples_exact": exact,
            "samples_different": len(shared) - exact,
        })

    write_tsv(
        args.output_dir / "tier-summary.tsv",
        ["tier", "old_total", "new_total", "total_difference", "samples_exact", "samples_different"],
        summary,
    )
    write_tsv(
        args.output_dir / "sample-tier-discrepancies.tsv",
        ["IID", "tier", "old_count", "new_count", "difference"],
        discrepancies,
    )
    write_tsv(
        args.output_dir / "sample-set-summary.tsv",
        ["set", "count"],
        [
            {"set": "old", "count": len(old)},
            {"set": "new", "count": len(new)},
            {"set": "shared", "count": len(shared)},
            {"set": "old_only", "count": len(set(old) - set(new))},
            {"set": "new_only", "count": len(set(new) - set(old))},
        ],
    )
    print(f"shared_samples={len(shared)} discrepancies={len(discrepancies)}")


if __name__ == "__main__":
    main()

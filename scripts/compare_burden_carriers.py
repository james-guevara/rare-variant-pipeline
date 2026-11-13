#!/usr/bin/env python3
"""Compare carrier-level burden rows from legacy and targeted pipeline outputs."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow.parquet as pq


TIERS = ("lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4")


def chrom(value):
    value = str(value)
    return value[3:] if value.startswith("chr") else value


def key(row, *, sample, chromosome, tier):
    return (
        str(row[sample]), chrom(row[chromosome]), str(row["POS"]),
        str(row["REF"]), str(row["ALT"]), str(row[tier]),
    )


def selected(row, names):
    return {name: row.get(name, "") for name in names}


def load_legacy(directory):
    rows = []
    columns = ["SAMPLE", "#CHROM", "POS", "REF", "ALT", "tier", "Gene", "SYMBOL", "Feature", "Consequence", "LoF", "LoF_filter"]
    for path in sorted(directory.glob("chr*.filtered_annotated.parquet")):
        for row in pq.read_table(path, columns=columns).to_pylist():
            if row["tier"] in TIERS:
                rows.append((key(row, sample="SAMPLE", chromosome="#CHROM", tier="tier"), selected(row, columns[6:])))
    return rows


def load_targeted(directory):
    rows = []
    columns = ["sample_id", "CHROM", "POS", "REF", "ALT", "Gene", "SYMBOL", "Consequence", "LoF", "lof_tier", "miss_tier"]
    for path in sorted(directory.glob("chr*/11.*-burden-eligible.parquet")):
        tier_column = "lof_tier" if ".plof-" in path.name else "miss_tier"
        for row in pq.read_table(path, columns=columns).to_pylist():
            if row[tier_column] in TIERS:
                rows.append((key(row, sample="sample_id", chromosome="CHROM", tier=tier_column), selected(row, columns[5:])))
    return rows


def write_difference(path, differences, metadata):
    fields = ["IID", "CHROM", "POS", "REF", "ALT", "tier", "multiplicity", "Gene", "SYMBOL", "Feature", "Consequence", "LoF", "LoF_filter"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for item in sorted(differences):
            sample, chromosome, pos, ref, alt, tier = item
            row = {"IID": sample, "CHROM": chromosome, "POS": pos, "REF": ref, "ALT": alt, "tier": tier, "multiplicity": differences[item]}
            row.update(metadata[item][0])
            writer.writerow({field: row.get(field, "") for field in fields})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-dir", required=True, type=Path)
    parser.add_argument("--targeted-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    old_rows = load_legacy(args.legacy_dir)
    new_rows = load_targeted(args.targeted_dir)
    old, new = Counter(k for k, _ in old_rows), Counter(k for k, _ in new_rows)
    old_meta, new_meta = defaultdict(list), defaultdict(list)
    for item, meta in old_rows:
        old_meta[item].append(meta)
    for item, meta in new_rows:
        new_meta[item].append(meta)

    old_only, new_only = old - new, new - old
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_difference(args.output_dir / "carrier-rows-old-only.tsv", old_only, old_meta)
    write_difference(args.output_dir / "carrier-rows-new-only.tsv", new_only, new_meta)
    with (args.output_dir / "carrier-tier-summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tier", "old_rows", "new_rows", "exact_shared", "old_only", "new_only"])
        for tier in TIERS:
            old_n = sum(n for k, n in old.items() if k[-1] == tier)
            new_n = sum(n for k, n in new.items() if k[-1] == tier)
            old_only_n = sum(n for k, n in old_only.items() if k[-1] == tier)
            new_only_n = sum(n for k, n in new_only.items() if k[-1] == tier)
            writer.writerow([tier, old_n, new_n, old_n - old_only_n, old_only_n, new_only_n])
    print(f"old_rows={sum(old.values())} new_rows={sum(new.values())} old_only={sum(old_only.values())} new_only={sum(new_only.values())}")


if __name__ == "__main__":
    main()

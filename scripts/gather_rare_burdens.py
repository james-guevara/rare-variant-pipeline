#!/usr/bin/env python3
"""Gather checksummed chromosome packages into participant-level rare burdens."""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

TIERS = ("lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4")


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_counts(target, path, tiers, suffix=""):
    for row in read_tsv(path):
        sample = row.get("SAMPLE") or row.get("sample_id")
        if not sample:
            raise ValueError(f"missing sample identifier in {path}")
        for tier in tiers:
            column = tier + suffix
            if row.get(column, "") != "":
                target[sample][tier] += int(row[column])


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-manifest", required=True, type=Path)
    parser.add_argument("--package", action="append", required=True, type=Path)
    parser.add_argument("--expected-chromosomes", default=",".join([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--strata-output", required=True, type=Path)
    args = parser.parse_args()

    samples = read_tsv(args.sample_manifest)
    if not samples or "IID" not in samples[0]:
        raise ValueError("sample manifest must contain IID")
    ids = [row["IID"] for row in samples]
    if len(ids) != len(set(ids)):
        raise ValueError("sample manifest contains duplicate IID values")

    expected = set(filter(None, args.expected_chromosomes.split(",")))
    packages = {}
    for root in args.package:
        document = json.loads((root / "targeted-output-manifest.json").read_text())
        chrom = document["chromosome"]
        if document["status"] != "SUCCEEDED" or chrom in packages:
            raise ValueError(f"invalid or duplicate chromosome package: {chrom}")
        for entry in document["files"].values():
            path = root / entry["file"]
            if not path.is_file() or digest(path) != entry["sha256"]:
                raise ValueError(f"package checksum failure: {path}")
        packages[chrom] = (root, document)
    if set(packages) != expected:
        raise ValueError(f"chromosome packages differ from expected: missing={sorted(expected-set(packages))} extra={sorted(set(packages)-expected)}")

    partitions = {
        "autosomal": sorted(expected - {"chrX", "chrY"}),
        "sex_chromosome_primary": sorted(expected & {"chrX", "chrY"}),
    }
    partition_counts = {name: defaultdict(lambda: defaultdict(int)) for name in partitions}
    partition_complete = {
        name: {tier: bool(chroms) for tier in TIERS}
        for name, chroms in partitions.items()
    }
    sensitivity = defaultdict(lambda: defaultdict(int))
    sensitivity_complete = bool(expected & {"chrX", "chrY"})

    for partition, chroms in partitions.items():
        for chrom in chroms:
            root, document = packages[chrom]
            files = document["files"]
            add_counts(partition_counts[partition], root / files["plof_counts"]["file"], TIERS[:2])
            if "missense_counts" in files:
                add_counts(partition_counts[partition], root / files["missense_counts"]["file"], TIERS[2:])
            else:
                for tier in TIERS[2:]: partition_complete[partition][tier] = False
            if chrom in {"chrX", "chrY"}:
                if "plof_sensitivity_burden" in files:
                    add_counts(sensitivity, root / files["plof_sensitivity_burden"]["file"], TIERS[:2], "_variants")
                else:
                    sensitivity_complete = False

    primary_rows, strata_rows = [], []
    for sample in samples:
        iid = sample["IID"]
        base = {"FID": sample.get("FID", ""), "IID": iid, "SEX": sample.get("SEX", "")}
        primary = dict(base)
        for tier in TIERS:
            active = [p for p, chroms in partitions.items() if chroms]
            complete = bool(active) and all(partition_complete[p][tier] for p in active)
            primary[tier] = sum(partition_counts[p][iid][tier] for p in active) if complete else ""
        primary_rows.append(primary)
        for partition in partitions:
            row = dict(base, burden_partition=partition, included_in_primary_total="true", chromosomes_complete="true")
            for tier in TIERS:
                row[tier] = partition_counts[partition][iid][tier] if partition_complete[partition][tier] else ""
            strata_rows.append(row)
        row = dict(base, burden_partition="sex_chromosome_sensitivity", included_in_primary_total="false", chromosomes_complete=str(sensitivity_complete).lower())
        for tier in TIERS:
            row[tier] = sensitivity[iid][tier] if sensitivity_complete and tier.startswith("lof_") else ""
        strata_rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["FID", "IID", "SEX", *TIERS]
    write_tsv(args.output, fields, primary_rows)
    write_tsv(args.strata_output, ["FID", "IID", "SEX", "burden_partition", "included_in_primary_total", "chromosomes_complete", *TIERS], strata_rows)


if __name__ == "__main__": main()

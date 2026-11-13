#!/usr/bin/env python3
"""Validate gathered cohort burdens against a pinned release contract."""

import argparse
import csv
import hashlib
import json
from pathlib import Path


TIERS = ("lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4")
PRIMARY_FIELDS = ("FID", "IID", "SEX", *TIERS)
STRATA_FIELDS = (
    "FID", "IID", "SEX", "burden_partition", "included_in_primary_total",
    "chromosomes_complete", *TIERS,
)
STRATA = {
    "autosomal", "sex_chromosome_primary", "sex_chromosome_sensitivity"
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_rows(path: Path, expected_fields) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(expected_fields):
            raise ValueError(
                f"unexpected columns in {path}: {reader.fieldnames}"
            )
        return list(reader)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--burdens", required=True, type=Path)
    parser.add_argument("--strata", required=True, type=Path)
    args = parser.parse_args()

    contract = json.loads(args.contract.read_text())
    expected = contract["expected_outputs"]
    observed_hashes = {
        "rare_burdens_sha256": digest(args.burdens),
        "chromosome_strata_sha256": digest(args.strata),
    }
    for name, observed in observed_hashes.items():
        if observed != expected[name]:
            raise ValueError(
                f"{name} differs: expected={expected[name]} observed={observed}"
            )

    burdens = read_rows(args.burdens, PRIMARY_FIELDS)
    strata = read_rows(args.strata, STRATA_FIELDS)
    sample_ids = [row["IID"] for row in burdens]
    if len(burdens) != expected["sample_count"]:
        raise ValueError(
            f"sample count differs: expected={expected['sample_count']} "
            f"observed={len(burdens)}"
        )
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("primary burden table contains duplicate IID values")
    if len(strata) != expected["strata_rows"]:
        raise ValueError(
            f"strata rows differ: expected={expected['strata_rows']} "
            f"observed={len(strata)}"
        )
    strata_by_sample = {}
    for row in strata:
        strata_by_sample.setdefault(row["IID"], set()).add(row["burden_partition"])
    if set(strata_by_sample) != set(sample_ids):
        raise ValueError("strata and primary sample sets differ")
    invalid = {
        sample: partitions
        for sample, partitions in strata_by_sample.items()
        if partitions != STRATA
    }
    if invalid:
        raise ValueError(f"samples have incomplete burden strata: {invalid}")

    observed_sums = {}
    for tier in TIERS:
        if any(row[tier] == "" for row in burdens):
            raise ValueError(f"primary burden table has missing {tier} values")
        observed_sums[tier] = sum(int(row[tier]) for row in burdens)
    if observed_sums != expected["primary_tier_sums"]:
        raise ValueError(
            f"tier sums differ: expected={expected['primary_tier_sums']} "
            f"observed={observed_sums}"
        )

    print(
        f"release=PASS id={contract['release_id']} samples={len(burdens)} "
        f"strata={len(strata)} hashes=2 tiers={len(TIERS)}"
    )


if __name__ == "__main__":
    main()

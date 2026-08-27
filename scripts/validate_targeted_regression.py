#!/usr/bin/env python3
"""Validate a completed targeted run against pinned expected counts."""

import argparse
import hashlib
import json
from pathlib import Path

import duckdb


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regression", required=True, type=Path)
    parser.add_argument("--bindings", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()

    expected_document = json.loads(args.regression.read_text())
    expected = expected_document["expected_counts"]
    bindings = json.loads(args.bindings.read_text())["resources"]
    target_bed = Path(bindings["target_bed"])
    candidates = Path(bindings["missense_candidates"])
    run = args.run_root
    if not (run / "_SUCCESS").exists():
        raise ValueError(f"run is not complete: {run}")
    for name, path in (("target_bed", target_bed), ("missense_candidates", candidates)):
        observed_digest = digest(path)
        expected_digest = expected_document[f"{name}_sha256"]
        if observed_digest != expected_digest:
            raise ValueError(
                f"{name} checksum differs: expected={expected_digest} observed={observed_digest}"
            )

    con = duckdb.connect()
    parquet = lambda name: str(run / name)
    tsv = lambda name: str(run / name)
    scalar = lambda query, values=(): con.execute(query, values).fetchone()[0]
    observed = {
        "target_intervals": sum(
            1 for line in target_bed.read_text().splitlines()
            if line and not line.startswith("#")
        ),
        "targeted_alleles": scalar(
            "SELECT count(*) FROM read_parquet(?)", [parquet("01.target-alleles.parquet")]
        ),
        "plof_rows": scalar(
            "SELECT count(*) FROM read_csv(?, delim='\\t', header=true)",
            [tsv("06.plof-genebayes.tsv")],
        ),
        "qualifying_lof": scalar(
            "SELECT count(*) FROM read_csv(?, delim='\\t', header=true)",
            [tsv("06.plof-tiered.tsv")],
        ),
        "lof_carrier_rows": scalar(
            "SELECT count(*) FROM read_parquet(?)",
            [parquet("07.plof-tiered.carriers.parquet")],
        ),
        "lof_carrier_samples": scalar(
            "SELECT count(DISTINCT sample_id) FROM read_parquet(?)",
            [parquet("07.plof-tiered.carriers.parquet")],
        ),
        "observed_missense_candidates": scalar(
            "SELECT count(*) FROM read_parquet(?)", [str(candidates)]
        ),
        "selected_missense": scalar(
            "SELECT count(*) FROM read_parquet(?)",
            [parquet("06.missense-tiered.parquet")],
        ),
        "missense_carrier_rows": scalar(
            "SELECT count(*) FROM read_parquet(?)",
            [parquet("07.missense-tiered.carriers.parquet")],
        ),
        "region_filtered_rows": scalar(
            "SELECT count(*) FROM read_parquet(?)",
            [parquet("08.missense-region-filtered.parquet")],
        ),
        "genotype_qc_rows": scalar(
            "SELECT count(*) FROM read_parquet(?)",
            [parquet("09.missense-genotype-qc.parquet")],
        ),
        "population_af_eligible_rows": scalar(
            "SELECT count(*) FROM read_parquet(?)",
            [parquet("10.missense-population-af-eligible.parquet")],
        ),
        "burden_eligible_rows": scalar(
            "SELECT count(*) FROM read_parquet(?)",
            [parquet("11.missense-burden-eligible.parquet")],
        ),
        "burden_eligible_alleles": scalar(
            "SELECT count(DISTINCT record_id) FROM read_parquet(?)",
            [parquet("11.missense-burden-eligible.parquet")],
        ),
        "burden_eligible_samples": scalar(
            "SELECT count(DISTINCT sample_id) FROM read_parquet(?)",
            [parquet("11.missense-burden-eligible.parquet")],
        ),
    }
    for tier in ("lof_t1", "lof_t2"):
        observed[tier] = scalar(
            "SELECT count(*) FROM read_csv(?, delim='\\t', header=true) WHERE lof_tier = ?",
            [tsv("06.plof-tiered.tsv"), tier],
        )
    for tier in ("miss_t1", "miss_t2", "miss_t3", "miss_t4"):
        observed[f"final_{tier}_rows"] = scalar(
            "SELECT count(*) FROM read_parquet(?) WHERE miss_tier = ?",
            [parquet("11.missense-burden-eligible.parquet"), tier],
        )

    differences = {
        key: {"expected": value, "observed": observed.get(key)}
        for key, value in expected.items()
        if observed.get(key) != value
    }
    if differences:
        raise ValueError(f"regression counts differ: {json.dumps(differences, sort_keys=True)}")
    print(f"regression=PASS chromosome={args.regression.stem} counts={len(expected)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Validate completed sharded Zarr stores and advance a cohort preparation plan."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import zarr


REQUIRED_ARRAYS = {
    "sample_id", "variant_position", "variant_allele", "call_genotype",
    "call_GQ", "call_DP",
}


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"missing TSV header: {path}")
        return reader.fieldnames, list(reader)


def text_values(array):
    return [value.decode() if isinstance(value, bytes) else str(value) for value in array[:]]


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def positions_sorted(array, block_rows=1_000_000):
    previous = None
    for start in range(0, array.shape[0], block_rows):
        values = array[start:min(start + block_rows, array.shape[0])]
        if values.size and previous is not None and values[0] < previous:
            return False
        if values.size > 1 and np.any(np.diff(values) < 0):
            return False
        if values.size:
            previous = values[-1]
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-plan", required=True, type=Path)
    parser.add_argument("--sample-manifest", required=True, type=Path)
    parser.add_argument("--output-plan", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    fields, plan = read_tsv(args.preparation_plan)
    _, sample_rows = read_tsv(args.sample_manifest)
    expected_samples = [row["IID"] for row in sample_rows]
    reports = []
    for row in plan:
        if row["preparation_state"] != "READY_FOR_ZARR":
            raise ValueError(
                f"{row['chromosome']} is {row['preparation_state']}, expected READY_FOR_ZARR"
            )
        store_path = Path(row["zarr_store"])
        checkpoint = Path(f"{store_path}.complete")
        if not checkpoint.is_file():
            raise ValueError(f"missing conversion checkpoint: {checkpoint}")
        if not (store_path / "zarr.json").is_file():
            raise ValueError(f"missing Zarr v3 metadata: {store_path / 'zarr.json'}")
        store = zarr.open_group(str(store_path), mode="r", zarr_format=3)
        names = {name for name, _ in store.arrays()}
        missing = sorted(REQUIRED_ARRAYS - names)
        if missing:
            raise ValueError(f"{row['chromosome']} Zarr is missing arrays: {missing}")
        if "call_AD" in names:
            allele_depth_arrays = ["call_AD"]
            allele_depth_encoding = "AD"
        elif {"call_LAD", "call_LAA"}.issubset(names):
            allele_depth_arrays = ["call_LAD", "call_LAA"]
            allele_depth_encoding = "LAD+LAA"
        else:
            raise ValueError(
                f"{row['chromosome']} Zarr requires call_AD or call_LAD+call_LAA"
            )
        observed_samples = text_values(store["sample_id"])
        if observed_samples != expected_samples:
            raise ValueError(
                f"{row['chromosome']} Zarr sample order differs from sample manifest"
            )
        positions = store["variant_position"]
        n_variants = positions.shape[0]
        if not positions_sorted(positions):
            raise ValueError(f"{row['chromosome']} variant positions are not sorted")
        for name in ("call_genotype", "call_GQ", "call_DP", *allele_depth_arrays):
            array = store[name]
            if array.shape[0] != n_variants or array.shape[1] != len(expected_samples):
                raise ValueError(
                    f"{row['chromosome']} {name} dimensions disagree with variants/samples"
                )
        row["preparation_state"] = "READY_FOR_DERIVED_RESOURCES"
        reports.append({
            "chromosome": row["chromosome"],
            "zarr_store": str(store_path),
            "variants": n_variants,
            "samples": len(observed_samples),
            "arrays": len(names),
            "allele_depth_encoding": allele_depth_encoding,
            "status": "PASS",
        })

    args.output_plan.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_plan, fields, plan)
    args.report.write_text(json.dumps({
        "schema_version": 1,
        "status": "PASS",
        "chromosomes": reports,
        "next_stage": "derive_cohort_resources",
    }, indent=2, sort_keys=True) + "\n")
    print("vcz_validation=PASS")
    print(f"chromosomes={len(reports)}")


if __name__ == "__main__":
    main()

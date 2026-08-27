#!/usr/bin/env python3
"""Fail unless a targeted workflow run exactly matches its validated reference."""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import duckdb


TSV_OUTPUT = "06.plof-tiered.tsv"
PARQUET_OUTPUTS = (
    "07.plof-tiered.carriers.parquet",
    "07.plof-tiered.genotype-summary.parquet",
)


def require_file(root: Path, name: str) -> Path:
    path = root / name
    if not path.is_file():
        raise AssertionError(f"missing output: {path}")
    return path


def tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AssertionError(f"missing TSV header: {path}")
        return reader.fieldnames, list(reader)


def canonical_rows(columns: list[str], rows: list[dict[str, str]]) -> Counter:
    return Counter(tuple(row[column] for column in columns) for row in rows)


def parquet_columns(connection: duckdb.DuckDBPyConnection, path: Path) -> list[str]:
    return [
        row[0]
        for row in connection.execute(
            "DESCRIBE SELECT * FROM read_parquet(?)", [str(path)]
        ).fetchall()
    ]


def parquet_difference(
    connection: duckdb.DuckDBPyConnection, left: Path, right: Path
) -> int:
    return connection.execute(
        "SELECT COUNT(*) FROM ("
        "SELECT * FROM read_parquet(?) EXCEPT ALL SELECT * FROM read_parquet(?)"
        ")",
        [str(left), str(right)],
    ).fetchone()[0]


def check_equal(new_root: Path, reference_root: Path) -> dict[str, int]:
    new_tsv = require_file(new_root, TSV_OUTPUT)
    reference_tsv = require_file(reference_root, TSV_OUTPUT)
    new_columns, new_rows = tsv_rows(new_tsv)
    reference_columns, reference_rows = tsv_rows(reference_tsv)
    if new_columns != reference_columns:
        raise AssertionError(
            f"{TSV_OUTPUT}: schema differs: new={new_columns}, reference={reference_columns}"
        )
    new_counter = canonical_rows(new_columns, new_rows)
    reference_counter = canonical_rows(reference_columns, reference_rows)
    if new_counter != reference_counter:
        raise AssertionError(
            f"{TSV_OUTPUT}: rows differ: "
            f"new_only={sum((new_counter-reference_counter).values())}, "
            f"reference_only={sum((reference_counter-new_counter).values())}"
        )

    connection = duckdb.connect()
    counts = {"qualifying_variants": len(new_rows)}
    for name in PARQUET_OUTPUTS:
        new_path = require_file(new_root, name)
        reference_path = require_file(reference_root, name)
        new_schema = parquet_columns(connection, new_path)
        reference_schema = parquet_columns(connection, reference_path)
        if new_schema != reference_schema:
            raise AssertionError(
                f"{name}: schema differs: new={new_schema}, reference={reference_schema}"
            )
        new_only = parquet_difference(connection, new_path, reference_path)
        reference_only = parquet_difference(connection, reference_path, new_path)
        if new_only or reference_only:
            raise AssertionError(
                f"{name}: rows differ: new_only={new_only}, reference_only={reference_only}"
            )

    counts["targeted_alleles"] = connection.execute(
        "SELECT COUNT(*) FROM read_parquet(?)",
        [str(require_file(new_root, "01.target-alleles.parquet"))],
    ).fetchone()[0]
    carriers = require_file(new_root, PARQUET_OUTPUTS[0])
    counts["carrier_rows"] = connection.execute(
        "SELECT COUNT(*) FROM read_parquet(?)", [str(carriers)]
    ).fetchone()[0]
    counts["carrier_samples"] = connection.execute(
        "SELECT COUNT(DISTINCT sample_id) FROM read_parquet(?)", [str(carriers)]
    ).fetchone()[0]
    counts["lof_t1"] = sum(row.get("lof_tier") == "lof_t1" for row in new_rows)
    counts["lof_t2"] = sum(row.get("lof_tier") == "lof_t2" for row in new_rows)
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--new", required=True, type=Path, help="new workflow run root")
    parser.add_argument(
        "--reference", required=True, type=Path, help="validated workflow run root"
    )
    parser.add_argument("--expectations", type=Path)
    args = parser.parse_args()

    counts = check_equal(args.new, args.reference)
    if args.expectations:
        expected = json.loads(args.expectations.read_text())["expected_counts"]
        mismatches = {
            key: {"expected": value, "observed": counts.get(key)}
            for key, value in expected.items()
            if counts.get(key) != value
        }
        if mismatches:
            raise AssertionError(f"expected counts differ: {json.dumps(mismatches)}")

    for key, value in counts.items():
        print(f"{key}={value}")
    print("regression=PASS")


if __name__ == "__main__":
    main()

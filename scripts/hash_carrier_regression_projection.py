#!/usr/bin/env python3
"""Hash a canonical carrier projection for cross-system regression checks."""

import argparse
import hashlib
import json
from pathlib import Path

import duckdb


def quote(value: str) -> str:
    return '"{}"'.format(value.replace('"', '""'))


def normalized(value):
    if value is None:
        return None
    if isinstance(value, float):
        return format(value, ".17g")
    return str(value)


def digest(rows) -> str:
    payload = json.dumps(
        sorted(tuple(normalized(value) for value in row) for row in rows),
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--sample-column", required=True)
    parser.add_argument("--tier-column", required=True)
    parser.add_argument("--lof-only", action="store_true")
    args = parser.parse_args()

    fixed = [
        "#CHROM", "POS", "REF", "ALT", args.sample_column, "GT", "GQ", "DP",
        "AD", "FILTER", "Gene", "SYMBOL", "Consequence", "LoF",
        args.tier_column, "genebayes_post_mean",
    ]
    connection = duckdb.connect()
    where = (
        " WHERE {} IN ('lof_t1', 'lof_t2')".format(quote(args.tier_column))
        if args.lof_only else ""
    )
    rows = connection.execute(
        "SELECT {} FROM read_parquet(?){}".format(
            ",".join(map(quote, fixed)), where
        ),
        [str(args.input)],
    ).fetchall()
    variants = {tuple(row[index] for index in (0, 1, 2, 3, 10, 11, 12, 13, 14, 15)) for row in rows}
    coordinates = {tuple(row[index] for index in (1, 2, 3)) for row in rows}
    variant_core = {tuple(row[index] for index in (1, 2, 3, 10, 11, 12, 13, 14)) for row in rows}
    carrier_core = {tuple(row[index] for index in (1, 2, 3, 4, 5, 14)) for row in rows}
    print(f"rows={len(rows)}")
    print(f"variants={len(variants)}")
    print(f"carrier_sha256={digest(rows)}")
    print(f"variant_sha256={digest(variants)}")
    print(f"coordinate_sha256={digest(coordinates)}")
    print(f"variant_core_sha256={digest(variant_core)}")
    print(f"carrier_core_sha256={digest(carrier_core)}")
    for index, label in (
        (10, "Gene"), (11, "SYMBOL"), (12, "Consequence"),
        (13, "LoF"), (14, "tier"), (15, "genebayes_post_mean"),
    ):
        values = {tuple(row[item] for item in (1, 2, 3, index)) for row in rows}
        print(f"variant_field_sha256.{label}={digest(values)}")


if __name__ == "__main__":
    main()

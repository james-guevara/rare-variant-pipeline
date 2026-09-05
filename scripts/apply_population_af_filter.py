#!/usr/bin/env python3
"""Emit population-AF-eligible rows without altering the annotated input."""

import argparse
from pathlib import Path

import duckdb


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--column", default="gnomAD4.1_joint_AF")
    parser.add_argument("--max-af", default=0.001, type=float)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect()
    columns = {
        row[0]
        for row in connection.execute(
            "DESCRIBE SELECT * FROM read_parquet(?)", [str(args.input)]
        ).fetchall()
    }
    if args.column not in columns:
        raise ValueError(f"population AF column is absent: {args.column}")
    quoted = '"{}"'.format(args.column.replace('"', '""'))
    connection.read_parquet(str(args.input)).create_view("annotated_input")
    connection.execute(
        f"CREATE TEMP TABLE eligible AS SELECT * FROM annotated_input "
        f"WHERE TRY_CAST({quoted} AS DOUBLE) < ? "
        f"OR TRY_CAST({quoted} AS DOUBLE) IS NULL",
        [args.max_af],
    )
    connection.execute(
        "COPY eligible TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
        [str(args.output)],
    )
    rows = connection.execute(
        "SELECT COUNT(*) FROM read_parquet(?)", [str(args.input)]
    ).fetchone()[0]
    kept = connection.execute(
        "SELECT COUNT(*) FROM read_parquet(?)", [str(args.output)]
    ).fetchone()[0]
    print(f"rows={rows:,}->{kept:,}")
    print(f"population_af_column={args.column}")
    print(f"population_af_max={args.max_af}")


if __name__ == "__main__":
    main()

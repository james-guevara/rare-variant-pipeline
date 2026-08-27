#!/usr/bin/env python3
"""Compare new and reference Parquet outputs over their shared schema."""

import argparse
from pathlib import Path

import duckdb


def quoted(name: str) -> str:
    return '"{}"'.format(name.replace('"', '""'))


def columns(connection: duckdb.DuckDBPyConnection, path: Path) -> list[str]:
    return [
        row[0]
        for row in connection.execute(
            "DESCRIBE SELECT * FROM read_parquet(?)", [str(path)]
        ).fetchall()
    ]


def difference_count(connection, left, right, selected):
    projection = ", ".join(map(quoted, selected))
    query = (
        "SELECT COUNT(*) FROM ("
        f"SELECT {projection} FROM read_parquet(?) EXCEPT ALL "
        f"SELECT {projection} FROM read_parquet(?)"
        ")"
    )
    return connection.execute(query, [str(left), str(right)]).fetchone()[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--new", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    args = parser.parse_args()

    connection = duckdb.connect()
    new_columns = columns(connection, args.new)
    reference_columns = columns(connection, args.reference)
    shared = [name for name in reference_columns if name in new_columns]
    missing = [name for name in reference_columns if name not in new_columns]
    extra = [name for name in new_columns if name not in reference_columns]
    print(f"shared_columns={len(shared)}")
    print("missing_columns={}".format(",".join(missing) or "."))
    print("extra_columns={}".format(",".join(extra) or "."))
    print(
        "new_minus_reference={}".format(
            difference_count(connection, args.new, args.reference, shared)
        )
    )
    print(
        "reference_minus_new={}".format(
            difference_count(connection, args.reference, args.new, shared)
        )
    )
    for name in extra:
        column = quoted(name)
        nonempty = connection.execute(
            f"SELECT COUNT(*) FROM read_parquet(?) WHERE {column} IS NOT NULL "
            f"AND CAST({column} AS VARCHAR) NOT IN ('', '.')",
            [str(args.new)],
        ).fetchone()[0]
        print(f"extra_nonempty.{name}={nonempty}")


if __name__ == "__main__":
    main()

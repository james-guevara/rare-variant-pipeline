#!/usr/bin/env python3
"""Consolidate chromosome-partitioned carrier Parquets into one genome-wide file."""

import argparse
import hashlib
import json
from pathlib import Path

import duckdb


def quote(path: Path) -> str:
    return "'{}'".format(str(path).replace("'", "''"))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--validation-output", required=True, type=Path)
    args = parser.parse_args()

    if len({path.resolve() for path in args.input}) != len(args.input):
        parser.error("duplicate input paths")
    missing = [str(path) for path in args.input if not path.is_file()]
    if missing:
        parser.error("missing inputs: " + ", ".join(missing))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.validation_output.parent.mkdir(parents=True, exist_ok=True)
    paths = "[{}]".format(", ".join(quote(path) for path in args.input))
    source = f"read_parquet({paths}, union_by_name=true)"
    con = duckdb.connect()
    input_rows = con.execute(f"SELECT count(*) FROM {source}").fetchone()[0]

    output = quote(args.output)
    con.execute(f"""
        COPY (
          SELECT * FROM {source}
          ORDER BY
            CASE
              WHEN regexp_extract("#CHROM", '[0-9]+') <> ''
                THEN CAST(regexp_extract("#CHROM", '[0-9]+') AS INTEGER)
              WHEN "#CHROM" IN ('chrX', 'X') THEN 23
              WHEN "#CHROM" IN ('chrY', 'Y') THEN 24
              ELSE 99
            END,
            POS, REF, ALT, sample_id, record_id
        ) TO {output} (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    output_rows = con.execute(
        "SELECT count(*) FROM read_parquet(?)", [str(args.output)]
    ).fetchone()[0]
    if output_rows != input_rows:
        raise RuntimeError(f"row-count mismatch: inputs={input_rows} output={output_rows}")

    # Order-independent fingerprints prove that every logical row was retained.
    columns = [column[0] for column in con.execute(f"SELECT * FROM {source} LIMIT 0").description]
    packed = ", ".join(f'"{column}"' for column in columns)
    input_fingerprint = con.execute(
        f"SELECT count(*), sum(hash({packed})), bit_xor(hash({packed})) FROM {source}"
    ).fetchone()
    output_fingerprint = con.execute(
        f"SELECT count(*), sum(hash({packed})), bit_xor(hash({packed})) "
        "FROM read_parquet(?)", [str(args.output)]
    ).fetchone()
    if input_fingerprint != output_fingerprint:
        raise RuntimeError("row-content fingerprint mismatch")

    report = {
        "status": "PASS",
        "input_files": len(args.input),
        "rows": output_rows,
        "columns": len(columns),
        "output": str(args.output),
        "output_size_bytes": args.output.stat().st_size,
        "output_sha256": digest(args.output),
        "row_fingerprint": [str(value) for value in output_fingerprint],
    }
    args.validation_output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()

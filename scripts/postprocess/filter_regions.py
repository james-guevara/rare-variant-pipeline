#!/usr/bin/env python3
"""Filter out variants overlapping problematic region BEDs.

Anti-join the input parquet against:
  - genomicSuperDups.bed
  - simpleRepeat.bed
  - rmsk.bed restricted to repClass IN ('Simple_repeat', 'Low_complexity')

A variant is dropped if its POS (1-based) falls inside ANY interval (BED
start < POS <= end). All carrier rows for a dropped variant are removed.

Usage:
    python filter_regions.py --cohort ssc --chrom chr22 --resources config/resources.yaml
"""
import argparse
import sys
import json
import duckdb
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--resources", required=True)
    ap.add_argument("--input", help="override input parquet path "
                    "(default: derived from --resources)")
    ap.add_argument("--output", help="override output parquet path "
                    "(default: derived from --resources)")
    args = ap.parse_args()

    cfg = json.loads(Path(args.resources).read_text())
    regions_dir = Path(cfg["regions_dir"].format(chrom=args.chrom))
    if not regions_dir.is_absolute():
        regions_dir = Path(args.resources).resolve().parent / regions_dir
    if args.input:
        in_parquet = Path(args.input)
    else:
        in_dir = Path(cfg["cohorts"][args.cohort]["input_dir"])
        in_parquet = in_dir / f"{args.chrom}.merged.parquet"
    if args.output:
        out_parquet = Path(args.output)
        out_parquet.parent.mkdir(parents=True, exist_ok=True)
    else:
        out_dir = Path(cfg["output_base"]) / args.cohort / "region_filtered"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_parquet = out_dir / f"{args.chrom}.parquet"

    rmsk_bed = regions_dir / cfg["rmsk_bed"].format(chrom=args.chrom)
    rmsk_classes = "(" + ", ".join(repr(c) for c in cfg["rmsk_drop_classes"]) + ")"

    print(f"[{args.cohort} {args.chrom}] filter_regions: {in_parquet} -> {out_parquet}", file=sys.stderr)

    con = duckdb.connect()

    # Build a single combined drop-interval table for this chrom.
    drop_table_parts = []
    for bed_name in cfg["filter_beds"]:
        bed_name = bed_name.format(chrom=args.chrom)
        bed_path = regions_dir / bed_name
        drop_table_parts.append(f"""
            SELECT CAST(column1 AS BIGINT) AS start_, CAST(column2 AS BIGINT) AS end_
            FROM read_csv('{bed_path}', delim='\t', header=false, all_varchar=true)
            WHERE column0 = '{args.chrom}'
        """)
    drop_table_parts.append(f"""
        SELECT CAST(column1 AS BIGINT) AS start_, CAST(column2 AS BIGINT) AS end_
        FROM read_csv('{rmsk_bed}', delim='\t', header=false, all_varchar=true)
        WHERE column0 = '{args.chrom}' AND column6 IN {rmsk_classes}
    """)
    union = " UNION ALL ".join(drop_table_parts)
    con.execute(f"CREATE TABLE _drop AS {union}")
    n_intervals = con.execute("SELECT COUNT(*) FROM _drop").fetchone()[0]
    print(f"  drop intervals: {n_intervals:,}", file=sys.stderr)

    rows_in = con.execute(f"SELECT COUNT(*) FROM read_parquet('{in_parquet}')").fetchone()[0]

    con.execute(f"""
        COPY (
            SELECT c.*
            FROM read_parquet('{in_parquet}') c
            WHERE NOT EXISTS (
                SELECT 1 FROM _drop d
                WHERE CAST(c.POS AS BIGINT) > d.start_
                  AND CAST(c.POS AS BIGINT) <= d.end_
            )
        ) TO '{out_parquet}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    rows_out = con.execute(f"SELECT COUNT(*) FROM read_parquet('{out_parquet}')").fetchone()[0]
    dropped = rows_in - rows_out
    pct = 100 * dropped / max(rows_in, 1)
    print(f"  rows: {rows_in:,} -> {rows_out:,} (dropped {dropped:,} = {pct:.1f}%)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Left-join population AF columns from dbNSFP onto the qc_filtered parquets.

Adds 14 population AF columns (no scores — those are in a later step):
  gnomAD4.1_joint_AF, gnomAD4.1_joint_POPMAX_AF, gnomAD4.1_joint_nhomalt,
  gnomAD4.1_joint_flag, AllofUs_ALL_AF, AllofUs_POPMAX_AF, 1000Gp3_AF,
  ALFA_Total_AF, RegeneronME_ALL_AF, TOPMed_frz8_AC, TOPMed_frz8_AF,
  TOPMed_frz8_AN, dbNSFP_POPMAX_AC, dbNSFP_POPMAX_AF

This stage is annotation-only: it never removes input rows. Population AF fields
are retained for analysis and later, explicitly named eligibility steps.

Join keys: dbNSFP uses '#chr' (no 'chr' prefix), 'pos(1-based)', 'ref', 'alt'.

Usage:
    python join_pop_af.py --cohort ssc --chrom chr22 --resources config/resources.json
"""
import argparse
import json
import sys
import duckdb
from pathlib import Path

AF_COLS = [
    "gnomAD4.1_joint_AF",
    "gnomAD4.1_joint_POPMAX_AF",
    "gnomAD4.1_joint_nhomalt",
    "gnomAD4.1_joint_flag",
    "AllofUs_ALL_AF",
    "AllofUs_POPMAX_AF",
    "1000Gp3_AF",
    "ALFA_Total_AF",
    "RegeneronME_ALL_AF",
    "TOPMed_frz8_AC",
    "TOPMed_frz8_AF",
    "TOPMed_frz8_AN",
    "dbNSFP_POPMAX_AC",
    "dbNSFP_POPMAX_AF",
]


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
    dbnsfp_dir = Path(cfg["dbnsfp_af_dir"])
    if not dbnsfp_dir.is_absolute():
        dbnsfp_dir = Path(args.resources).resolve().parent / dbnsfp_dir
    base = Path(cfg.get("output_base", ".")) / args.cohort
    in_parquet = Path(args.input) if args.input else base / "qc_filtered" / f"{args.chrom}.parquet"
    if args.output:
        out_parquet = Path(args.output)
        out_parquet.parent.mkdir(parents=True, exist_ok=True)
    else:
        out_dir = base / "with_pop_af"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_parquet = out_dir / f"{args.chrom}.parquet"
    dbnsfp_parquet = Path(dbnsfp_dir) / f"{args.chrom}.parquet"

    print(f"[{args.cohort} {args.chrom}] join_pop_af: {in_parquet} -> {out_parquet}", file=sys.stderr)
    print(f"  dbNSFP: {dbnsfp_parquet}", file=sys.stderr)

    con = duckdb.connect()
    rows_in = con.execute(f"SELECT COUNT(*) FROM read_parquet('{in_parquet}')").fetchone()[0]

    # Build dbNSFP table: dedup on (chrom, pos, ref, alt), pre-cast AF cols.
    af_cols_quoted = ", ".join(f'"{c}"' for c in AF_COLS)
    con.execute(f"""
        CREATE TABLE _dbnsfp AS
        SELECT DISTINCT
               '{args.chrom}' AS "#CHROM",
               CAST("pos(1-based)" AS BIGINT) AS POS,
               ref AS REF,
               alt AS ALT,
               {af_cols_quoted}
        FROM read_parquet('{dbnsfp_parquet}')
    """)
    n_dbnsfp = con.execute("SELECT COUNT(*) FROM _dbnsfp").fetchone()[0]
    print(f"  dbNSFP rows for {args.chrom}: {n_dbnsfp:,}", file=sys.stderr)

    con.execute(f"""
        COPY (
            SELECT c.*, {af_cols_quoted}
            FROM read_parquet('{in_parquet}') c
            LEFT JOIN _dbnsfp d
              ON c."#CHROM" = d."#CHROM"
             AND CAST(c.POS AS BIGINT) = d.POS
             AND c.REF = d.REF
             AND c.ALT = d.ALT
        ) TO '{out_parquet}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    rows_out = con.execute(f"SELECT COUNT(*) FROM read_parquet('{out_parquet}')").fetchone()[0]
    if rows_out != rows_in:
        raise RuntimeError(
            f"population-AF annotation changed row count: {rows_in:,} -> {rows_out:,}"
        )
    n_hit = con.execute(f"""
        SELECT COUNT(*) FROM read_parquet('{out_parquet}')
        WHERE "gnomAD4.1_joint_AF" IS NOT NULL
    """).fetchone()[0]
    pct_hit = 100 * n_hit / max(rows_out, 1)
    print(f"  rows: {rows_in:,} -> {rows_out:,} (annotation only; dropped 0)", file=sys.stderr)
    print(f"  gnomAD AF non-null: {n_hit:,} = {pct_hit:.1f}%", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

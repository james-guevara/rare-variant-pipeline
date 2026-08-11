#!/usr/bin/env python3
"""Final filter step: keep FILTER==PASS rows, plus confirmed-DNM rescue.

The upstream stages do not pre-filter FILTER!=PASS, so every per-stage parquet
carries non-PASS variants through. This step applies the FILTER cut and the DNM
rescue in one pass. It requires a DNM callset for the cohort (`dnm_paths` in the
resources JSON), so it is skipped for cohorts that have none.

Logic:
  keep rows where:
      (FILTER passes: NULL/'PASS'/'.'/'')
    OR
      (FILTER fails AND (SAMPLE, chrom, pos, ref, alt) is in cohort DNM list
       with synthdnm_class=1 AND synthdnm_prob >= min_prob)

The DNM list is already pre-filtered to class=1 + prob>=0.5 in the per-cohort
parquet referenced by `dnm_paths` in resources.json, so we just need to
intersect on the variant key + SAMPLE.

Usage:
    python filter_pass_and_rescue.py --cohort ssc --chrom chr22 \\
        --resources config/resources.json
"""
import argparse
import json
import sys
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
    base = Path(cfg["output_base"]) / args.cohort
    in_parquet = base / "with_gene_constraint" / f"{args.chrom}.parquet"
    out_dir = base / "filtered_with_rescue"
    if not args.output:
        out_dir.mkdir(parents=True, exist_ok=True)
    out_parquet = out_dir / f"{args.chrom}.parquet"

    # Explicit-path overrides (used by the Nextflow POSTPROCESS subworkflow).
    if args.input:
        in_parquet = Path(args.input)
    if args.output:
        out_parquet = Path(args.output)
        out_parquet.parent.mkdir(parents=True, exist_ok=True)

    dnm_path = cfg["dnm_paths"][args.cohort]

    print(f"[{args.cohort} {args.chrom}] filter_pass_and_rescue: {in_parquet} -> {out_parquet}", file=sys.stderr)
    print(f"  DNM list: {dnm_path}", file=sys.stderr)

    con = duckdb.connect()
    con.execute("SET threads=2; SET memory_limit='8GB'")

    rows_in = con.execute(f"SELECT COUNT(*) FROM read_parquet('{in_parquet}')").fetchone()[0]

    # FILTER distribution before cut
    filt_dist = con.execute(f"""
        SELECT FILTER, COUNT(*) AS n
        FROM read_parquet('{in_parquet}')
        GROUP BY FILTER ORDER BY n DESC
    """).fetchall()
    print("  FILTER distribution (input):", file=sys.stderr)
    for f, n in filt_dist[:6]:
        print(f"    {str(f):<40} {n:>12,}", file=sys.stderr)

    # Restrict DNM list to this chrom; class+prob already pre-filtered.
    con.execute(f"""
        CREATE TABLE _dnm AS
        SELECT DISTINCT "#CHROM" AS dnm_chrom,
                        CAST(POS AS BIGINT) AS dnm_pos,
                        REF AS dnm_ref,
                        ALT_specific AS dnm_alt,
                        SAMPLE AS dnm_sample
        FROM read_parquet('{dnm_path}')
        WHERE "#CHROM" = '{args.chrom}'
    """)
    n_dnm = con.execute("SELECT COUNT(*) FROM _dnm").fetchone()[0]
    print(f"  DNM list rows for {args.chrom}: {n_dnm:,}", file=sys.stderr)

    # Diagnostic: how many non-PASS rows match a DNM entry?
    n_nonpass = con.execute(f"""
        SELECT COUNT(*) FROM read_parquet('{in_parquet}')
        WHERE FILTER IS NOT NULL
          AND FILTER NOT IN ('PASS', '.', '')
    """).fetchone()[0]
    n_rescued = con.execute(f"""
        SELECT COUNT(*)
        FROM read_parquet('{in_parquet}') c
        JOIN _dnm d
          ON c."#CHROM" = d.dnm_chrom
         AND CAST(c.POS AS BIGINT) = d.dnm_pos
         AND c.REF = d.dnm_ref
         AND c.ALT = d.dnm_alt
         AND c.SAMPLE = d.dnm_sample
        WHERE c.FILTER IS NOT NULL
          AND c.FILTER NOT IN ('PASS', '.', '')
    """).fetchone()[0]
    print(f"  non-PASS rows in input: {n_nonpass:,}", file=sys.stderr)
    print(f"  non-PASS rows matching DNM list (rescued): {n_rescued:,}", file=sys.stderr)

    # Final write: keep PASS OR (non-PASS AND in DNM list).
    con.execute(f"""
        COPY (
            SELECT c.*
            FROM read_parquet('{in_parquet}') c
            LEFT JOIN _dnm d
              ON c."#CHROM" = d.dnm_chrom
             AND CAST(c.POS AS BIGINT) = d.dnm_pos
             AND c.REF = d.dnm_ref
             AND c.ALT = d.dnm_alt
             AND c.SAMPLE = d.dnm_sample
            WHERE
                (c.FILTER IS NULL OR c.FILTER IN ('PASS', '.', ''))
              OR
                (d.dnm_sample IS NOT NULL)
        ) TO '{out_parquet}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    rows_out = con.execute(f"SELECT COUNT(*) FROM read_parquet('{out_parquet}')").fetchone()[0]
    dropped = rows_in - rows_out
    pct = 100 * dropped / max(rows_in, 1)
    print(f"  rows: {rows_in:,} -> {rows_out:,} (dropped {dropped:,} = {pct:.1f}%)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

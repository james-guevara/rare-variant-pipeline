#!/usr/bin/env python3
"""Drop carrier rows failing genotype QC. Strictly per-row.

A row is KEPT only if all of:
  GQ >= min_gq
  DP >= min_dp
  diploid alt dosage 1 AND AB within [het_ab_min, het_ab_max]
  diploid alt dosage 2 OR haploid alt dosage 1 AND AB >= hom_ab_min

The decision uses mask-derived ``called_ploidy`` and ``alt_dosage`` rather than
the formatting of GT. Thus haploid calls encoded as ``1``, ``1/.``, or ``./1``
are equivalent. Any other state is dropped, as is any row where GQ, DP or AB
is missing -- the keep condition is wrapped in COALESCE(..., FALSE) so absent
values fail CLOSED rather than slipping through on SQL three-valued logic.

AB is computed from AD = "ad_ref,ad_alt".

NOTE there is deliberately no family-level propagation. An earlier version dropped
every row of a (variant x family) group when any member's row failed. That made
effective stringency depend on how many family members happened to carry the same
variant, which varies by cohort and by family, and its purpose was inheritance
inference that this pipeline no longer does. Nothing here reads a pedigree.

Usage:
    python qc_genotype.py --cohort ssc --chrom chr22 --resources config/resources.json
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
    ap.add_argument(
        "--passthrough-chroms",
        default="",
        help="Comma-separated chroms to copy through unchanged (no QC applied). "
             "Useful for chrX/chrY where ploidy assumptions break the AB filter.",
    )
    args = ap.parse_args()

    cfg = json.loads(Path(args.resources).read_text())
    base = Path(cfg.get("output_base", ".")) / args.cohort
    in_parquet = Path(args.input) if args.input else base / "region_filtered" / f"{args.chrom}.parquet"
    if args.output:
        out_parquet = Path(args.output)
        out_parquet.parent.mkdir(parents=True, exist_ok=True)
    else:
        out_dir = base / "qc_filtered"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_parquet = out_dir / f"{args.chrom}.parquet"

    passthrough = {c.strip() for c in args.passthrough_chroms.split(",") if c.strip()}
    if args.chrom in passthrough:
        print(f"[{args.cohort} {args.chrom}] qc_genotype: PASSTHROUGH (no QC) -> {out_parquet}", file=sys.stderr)
        con = duckdb.connect()
        rows_in = con.execute(f"SELECT COUNT(*) FROM read_parquet('{in_parquet}')").fetchone()[0]
        con.execute(f"""
            COPY (SELECT * FROM read_parquet('{in_parquet}'))
            TO '{out_parquet}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """)
        print(f"  rows: {rows_in:,} -> {rows_in:,} (passthrough)", file=sys.stderr)
        return 0

    q = cfg["qc"]
    print(f"[{args.cohort} {args.chrom}] qc_genotype: {in_parquet} -> {out_parquet}", file=sys.stderr)

    con = duckdb.connect()
    columns = {
        row[0]
        for row in con.execute(
            f"DESCRIBE SELECT * FROM read_parquet('{in_parquet}')"
        ).fetchall()
    }
    has_mask_derived_state = {"called_ploidy", "alt_dosage"}.issubset(columns)
    if has_mask_derived_state:
        genotype_state = """
                    (called_ploidy = 2 AND alt_dosage = 1
                       AND AB IS NOT NULL
                       AND AB BETWEEN {het_ab_min} AND {het_ab_max})
                    OR
                    ((called_ploidy = 2 AND alt_dosage = 2)
                       OR (called_ploidy = 1 AND alt_dosage = 1))
                       AND AB IS NOT NULL
                       AND AB >= {hom_ab_min}
        """.format(**q)
    else:
        # Compatibility for carrier Parquets produced before called_ploidy was
        # added. New production output always takes the branch above.
        genotype_state = """
                    (GT IN ('0/1','1/0','0|1','1|0')
                       AND AB IS NOT NULL
                       AND AB BETWEEN {het_ab_min} AND {het_ab_max})
                    OR
                    (GT IN ('1/1','1|1','1','1/.','./1','1|.','.|1')
                       AND AB IS NOT NULL
                       AND AB >= {hom_ab_min})
        """.format(**q)
    rows_in = con.execute(f"SELECT COUNT(*) FROM read_parquet('{in_parquet}')").fetchone()[0]

    # Per-row keep flag. COALESCE(..., FALSE) makes missing GQ/DP/AB fail closed.
    con.execute(f"""
        COPY (
            WITH parsed AS (
                SELECT *,
                       TRY_CAST(GQ AS DOUBLE) AS GQ_n,
                       TRY_CAST(DP AS DOUBLE) AS DP_n,
                       TRY_CAST(SPLIT_PART(AD, ',', 1) AS DOUBLE) AS AD_ref,
                       TRY_CAST(SPLIT_PART(AD, ',', 2) AS DOUBLE) AS AD_alt
                FROM read_parquet('{in_parquet}')
            ),
            with_ab AS (
                SELECT *,
                       CASE WHEN (AD_ref + AD_alt) > 0
                            THEN AD_alt / (AD_ref + AD_alt)
                            ELSE NULL END AS AB
                FROM parsed
            )
            SELECT * EXCLUDE (GQ_n, DP_n, AD_ref, AD_alt, AB)
            FROM with_ab
            WHERE COALESCE(
                GQ_n >= {q['min_gq']}
                AND DP_n >= {q['min_dp']}
                AND ({genotype_state}), FALSE)
        ) TO '{out_parquet}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)

    rows_out = con.execute(f"SELECT COUNT(*) FROM read_parquet('{out_parquet}')").fetchone()[0]
    dropped = rows_in - rows_out
    pct = 100 * dropped / max(rows_in, 1)
    print(f"  rows: {rows_in:,} -> {rows_out:,} (dropped {dropped:,} = {pct:.1f}%)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

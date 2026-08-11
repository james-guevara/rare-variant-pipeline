#!/usr/bin/env python3
"""Annotate carrier rows with LoF and missense tiers.

Single `tier` column (string), naming convention: t1 = highest tier across both.

  lof_t1   LoF = 'HC' AND s_het >= 0.18         (highest LoF tier, most constrained)
  lof_t2   LoF = 'HC' AND s_het in [0.03, 0.18)
  miss_t1  missense AND n_flag = 4              (highest missense tier, all 4 metrics)
  miss_t2  missense AND n_flag = 3
  miss_t3  missense AND n_flag = 2
  miss_t4  missense AND n_flag = 1              (lowest missense tier)

n_flag = sum of 4 indicators at v2 forward-selection T*'s:
    ClinPred_rankscore         >= 0.4298
    AlphaMissense_rankscore    >= 0.9603
    popEVE_converted_rankscore >= 0.9209
    MPC_rankscore              >= 0.8947

CONVENTION: t1 is always the MOST severe tier, uniformly for LoF and missense.
Anything consuming these labels inherits that ordering; do not renumber downstream.
Note that older per-sample counters outside this pipeline keyed their columns on
n_flag directly (so their t1 was the LEAST severe) and relabelled at report time —
if joining against those outputs, check which convention they used first.

ANNOTATES, does not filter. Every input row is emitted; `tier` is NULL for rows
that match no tier. That keeps this the pipeline's analysis-ready deliverable —
tiering is one stratification among several (gene sets, consequence class, s_het
bins), so dropping untiered rows here would prevent all the others. Consumers that
want only tiered rows should add `WHERE tier IS NOT NULL`.

`miss_n_flag` kept for inspection.

Usage:
    python tier_variants.py --cohort <cohort> --chrom chr22 --resources resources.json
"""
import argparse
import json
import sys
import duckdb
from pathlib import Path


T_STARS = {
    "ClinPred_rankscore":         0.4298,
    "AlphaMissense_rankscore":    0.9603,
    "popEVE_converted_rankscore": 0.9209,
    "MPC_rankscore":              0.8947,
}

LOF_T1 = 0.18
LOF_T2 = 0.03


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
    out_dir = base / "tiered"
    if not args.output:
        out_dir.mkdir(parents=True, exist_ok=True)
    out_parquet = out_dir / f"{args.chrom}.parquet"

    # Explicit-path overrides (used by the Nextflow POSTPROCESS subworkflow).
    if args.input:
        in_parquet = Path(args.input)
    if args.output:
        out_parquet = Path(args.output)
        out_parquet.parent.mkdir(parents=True, exist_ok=True)

    print(f"[{args.cohort} {args.chrom}] tier_variants: {in_parquet} -> {out_parquet}", file=sys.stderr)

    con = duckdb.connect()
    rows_in = con.execute(f"SELECT COUNT(*) FROM read_parquet('{in_parquet}')").fetchone()[0]

    flag_exprs = [
        f'CAST(COALESCE(TRY_CAST("{col}" AS DOUBLE) >= {t}, FALSE) AS INTEGER)'
        for col, t in T_STARS.items()
    ]
    n_flag_sql = " + ".join(flag_exprs)

    con.execute(f"""
        COPY (
            WITH scored AS (
                SELECT *,
                       CASE
                           WHEN Consequence LIKE '%missense_variant%' THEN ({n_flag_sql})
                           ELSE NULL
                       END AS miss_n_flag
                FROM read_parquet('{in_parquet}')
            ),
            tiered AS (
                SELECT *,
                       CASE
                           -- LoF tiers (t1 = highest, s_het >= 0.18)
                           WHEN LoF = 'HC' AND TRY_CAST(genebayes_post_mean AS DOUBLE) >= {LOF_T1} THEN 'lof_t1'
                           WHEN LoF = 'HC'
                                AND TRY_CAST(genebayes_post_mean AS DOUBLE) >= {LOF_T2}
                                AND TRY_CAST(genebayes_post_mean AS DOUBLE) <  {LOF_T1} THEN 'lof_t2'
                           -- Missense tiers (t1 = highest, n_flag = 4)
                           WHEN miss_n_flag = 4 THEN 'miss_t1'
                           WHEN miss_n_flag = 3 THEN 'miss_t2'
                           WHEN miss_n_flag = 2 THEN 'miss_t3'
                           WHEN miss_n_flag = 1 THEN 'miss_t4'
                           ELSE NULL
                       END AS tier
                FROM scored
            )
            SELECT *
            FROM tiered
        ) TO '{out_parquet}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    rows_out = con.execute(f"SELECT COUNT(*) FROM read_parquet('{out_parquet}')").fetchone()[0]

    breakdown = con.execute(f"""
        SELECT COALESCE(tier, '(untiered)') AS tier, COUNT(*) AS n
        FROM read_parquet('{out_parquet}')
        GROUP BY 1
        ORDER BY 1
    """).fetchall()

    print(f"  rows: {rows_in:,} -> {rows_out:,}", file=sys.stderr)
    for tier, n in breakdown:
        print(f"  {tier}: {n:,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

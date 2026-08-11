#!/usr/bin/env python3
"""Per-sample carrier counts from the filtered, annotated variant parquets.

Deliberately minimal, and deliberately different from
rv_postprocessing_v3/scripts/count_{lof_shet,miss_tier_split,syn_shet}.py:

  * ONE cohort. No UNION across cohorts, so a sample present in more than one
    callset cannot be double-counted.
  * NO covariates. No PCs, ancestry, case status or FID, and therefore no
    dependency on the analysis manifest. Attach covariates downstream.
  * NO ancestry filter. Every sample in the input is counted.
  * The grouping COLUMN is a parameter (--group-col, default `tier`) and its
    VALUES are read from the data, not hardcoded. Tiering is only one way to
    stratify: --group-col Consequence, or a gene-set membership column, work the
    same way with no code change. Because tier values come straight from
    tier_variants.py, the t1 = most-severe convention is inherited by
    construction, with no relabelling step that could invert severity.

Emits one row per sample with a column per group value, so a sample carrying
nothing in a group gets 0 rather than being absent. Rows whose group value is NULL
(e.g. untiered variants) contribute to no group column but the sample still appears.

Usage:
    count_carriers.py --input chr1.parquet chr2.parquet ... \
        --out-counts per_sample_counts.tsv --out-totals group_totals.tsv
    count_carriers.py --input ... --group-col Consequence --out-counts ... --out-totals ...
"""
import argparse
import sys

import duckdb


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", nargs="+", required=True,
                    help="tiered parquet file(s), one or more chromosomes")
    ap.add_argument("--out-counts", required=True, help="per-sample counts TSV")
    ap.add_argument("--out-totals", required=True, help="variant-level group totals TSV")
    ap.add_argument("--sample-col", default="SAMPLE")
    ap.add_argument("--group-col", default="tier",
                    help="column to stratify counts by (default: tier)")
    args = ap.parse_args()

    files = ", ".join(f"'{p}'" for p in args.input)
    src = f"read_parquet([{files}])"
    con = duckdb.connect()

    cols = {d[0] for d in con.execute(f"SELECT * FROM {src} LIMIT 0").description}
    for needed in (args.sample_col, args.group_col):
        if needed not in cols:
            print(f"ERROR: column {needed!r} not in input; found {sorted(cols)[:20]}...",
                  file=sys.stderr)
            return 1

    n_rows = con.execute(f"SELECT COUNT(*) FROM {src}").fetchone()[0]
    print(f"input files: {len(args.input)}  rows: {n_rows:,}", file=sys.stderr)

    # Tier vocabulary straight from the data, so a new tier in tier_variants.py
    # shows up here without a code change.
    tiers = [r[0] for r in con.execute(
        f"SELECT DISTINCT {args.group_col} FROM {src} "
        f"WHERE {args.group_col} IS NOT NULL ORDER BY 1").fetchall()]
    if not tiers:
        print(f"ERROR: no non-null {args.group_col} values found", file=sys.stderr)
        return 1
    print(f"group values found: {', '.join(tiers)}", file=sys.stderr)

    # Variant-level totals (sanity check / provenance)
    con.execute(f"""
        COPY (
            SELECT {args.group_col} AS "{args.group_col}",
                   COUNT(*)                                  AS n_carrier_rows,
                   COUNT(DISTINCT {args.sample_col})          AS n_samples_with_any
            FROM {src}
            WHERE {args.group_col} IS NOT NULL
            GROUP BY 1 ORDER BY 1
        ) TO '{args.out_totals}' (FORMAT CSV, DELIMITER '\t', HEADER)
    """)

    # Per-sample counts, one column per tier. Every sample present in the input
    # appears, including those with zero in a given tier.
    per_tier = ",\n               ".join(
        f"COUNT(*) FILTER (WHERE {args.group_col} = '{t}') AS {t}" for t in tiers)
    con.execute(f"""
        COPY (
            SELECT {args.sample_col} AS SAMPLE,
               {per_tier},
               COUNT(*) FILTER (WHERE {args.group_col} IS NOT NULL) AS any_group
            FROM {src}
            GROUP BY 1
            ORDER BY 1
        ) TO '{args.out_counts}' (FORMAT CSV, DELIMITER '\t', HEADER)
    """)

    n_samp = con.execute(f"SELECT COUNT(DISTINCT {args.sample_col}) FROM {src}").fetchone()[0]
    print(f"wrote {args.out_counts} ({n_samp:,} samples) and {args.out_totals}", file=sys.stderr)
    for t, n, s in con.execute(
            f"SELECT {args.group_col}, COUNT(*), COUNT(DISTINCT {args.sample_col}) FROM {src} "
            f"WHERE {args.group_col} IS NOT NULL GROUP BY 1 ORDER BY 1").fetchall():
        print(f"  {t:10s} rows={n:>10,}  samples={s:>6,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

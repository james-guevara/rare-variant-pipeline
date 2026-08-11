#!/usr/bin/env python3
"""Join per-variant pathogenicity scores from dbNSFP onto the with_pop_af parquets.

Source: dbNSFP 5.3.1a `parquet_expanded_mane_select` — one row per variant.

This parquet has mixed column encodings:
  - Scalar (one value per variant): CADD_phred, CADD_raw, ClinPred_score,
    and ALL `*_rankscore` columns (rankscores are a global percentile rank
    computed once per variant — not per-transcript).
  - Per-transcript semicolon-delimited lists: REVEL_score, MPC_score,
    popEVE_score, AlphaMissense_score (these models score each protein
    context separately).

So we use two extraction patterns:
  - SCALAR: TRY_CAST as-is.
  - LIST:   list_max/list_min over string_split + TRY_CAST per element.

After per-row extraction, aggregate to one row per (chrom,pos,ref,alt) with
MAX/MIN (handles the 0.03% of dup keys; values are identical in dups).

Output columns:
  Rankscores (used by v2 tier_split — all higher = worse, MAX agg):
    ClinPred_rankscore, AlphaMissense_rankscore, popEVE_converted_rankscore,
    MPC_rankscore, REVEL_rankscore
  Raw scores:
    ClinPred_score (scalar, MAX), AlphaMissense_score (list, MAX),
    REVEL_score (list, MAX), MPC_score (list, MAX),
    popEVE_score (list, MIN — lower = worse),
    CADD_phred (scalar, MAX)

Join keys: dbNSFP `#chr` (no 'chr' prefix) / `pos(1-based)` / `ref` / `alt`
matched against `#CHROM` / POS / REF / ALT.

Usage:
    python join_scores.py --cohort ssc --chrom chr22 --resources config/resources.json
"""
import argparse
import json
import sys
import duckdb
from pathlib import Path

# (dbnsfp_col, out_col, encoding, agg)
#   encoding: "scalar" or "list"
#   agg:      "MAX" or "MIN"
SCORES = [
    # rankscores — all scalar, MAX (higher = worse)
    ("ClinPred_rankscore",         "ClinPred_rankscore",         "scalar", "MAX"),
    ("AlphaMissense_rankscore",    "AlphaMissense_rankscore",    "scalar", "MAX"),
    ("popEVE_converted_rankscore", "popEVE_converted_rankscore", "scalar", "MAX"),
    ("MPC_rankscore",              "MPC_rankscore",              "scalar", "MAX"),
    ("REVEL_rankscore",            "REVEL_rankscore",            "scalar", "MAX"),
    # raw scores
    ("ClinPred_score",       "ClinPred_score",       "scalar", "MAX"),
    ("AlphaMissense_score",  "AlphaMissense_score",  "list",   "MAX"),
    ("REVEL_score",          "REVEL_score",          "list",   "MAX"),
    ("MPC_score",            "MPC_score",            "list",   "MAX"),
    ("popEVE_score",         "popEVE_score",         "list",   "MIN"),  # lower = worse
    ("CADD_phred",           "CADD_phred",           "scalar", "MAX"),
]


def extract_expr(src: str, encoding: str, agg: str) -> str:
    if encoding == "scalar":
        return f'TRY_CAST("{src}" AS DOUBLE)'
    # list encoding: split, cast each element, then list_max/list_min
    func = "list_max" if agg == "MAX" else "list_min"
    return (
        f'{func}(list_transform(string_split("{src}", \';\'), '
        f'x -> TRY_CAST(x AS DOUBLE)))'
    )


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
    dbnsfp_dir = cfg["dbnsfp_scores_dir"]
    base = Path(cfg["output_base"]) / args.cohort
    in_parquet = base / "with_pop_af" / f"{args.chrom}.parquet"
    out_dir = base / "with_scores"
    if not args.output:
        out_dir.mkdir(parents=True, exist_ok=True)
    out_parquet = out_dir / f"{args.chrom}.parquet"

    # Explicit-path overrides (used by the Nextflow POSTPROCESS subworkflow).
    if args.input:
        in_parquet = Path(args.input)
    if args.output:
        out_parquet = Path(args.output)
        out_parquet.parent.mkdir(parents=True, exist_ok=True)
    dbnsfp_parquet = Path(dbnsfp_dir) / f"{args.chrom}.parquet"

    print(f"[{args.cohort} {args.chrom}] join_scores: {in_parquet} -> {out_parquet}", file=sys.stderr)
    print(f"  dbNSFP: {dbnsfp_parquet}", file=sys.stderr)

    con = duckdb.connect()
    rows_in = con.execute(f"SELECT COUNT(*) FROM read_parquet('{in_parquet}')").fetchone()[0]

    # Per-row extraction (one value per row, per output column).
    per_row_select = []
    for src, out, encoding, agg in SCORES:
        per_row_select.append(f'{extract_expr(src, encoding, agg)} AS "{out}"')
    per_row_sql = ",\n               ".join(per_row_select)

    con.execute(f"""
        CREATE TABLE _per_row AS
        SELECT '{args.chrom}' AS "#CHROM",
               CAST("pos(1-based)" AS BIGINT) AS POS,
               ref AS REF,
               alt AS ALT,
               {per_row_sql}
        FROM read_parquet('{dbnsfp_parquet}')
    """)

    # Aggregate to one row per (chrom,pos,ref,alt) using MAX/MIN per column.
    agg_select = []
    for src, out, encoding, agg in SCORES:
        agg_select.append(f'{agg}("{out}") AS "{out}"')
    agg_sql = ",\n               ".join(agg_select)

    con.execute(f"""
        CREATE TABLE _scores AS
        SELECT "#CHROM", POS, REF, ALT,
               {agg_sql}
        FROM _per_row
        GROUP BY "#CHROM", POS, REF, ALT
    """)
    n_scores = con.execute("SELECT COUNT(*) FROM _scores").fetchone()[0]
    print(f"  dbNSFP score rows for {args.chrom}: {n_scores:,}", file=sys.stderr)

    out_cols = [t[1] for t in SCORES]
    out_cols_quoted = ", ".join(f's."{c}"' for c in out_cols)

    con.execute(f"""
        COPY (
            SELECT c.*, {out_cols_quoted}
            FROM read_parquet('{in_parquet}') c
            LEFT JOIN _scores s
              ON c."#CHROM" = s."#CHROM"
             AND CAST(c.POS AS BIGINT) = s.POS
             AND c.REF = s.REF
             AND c.ALT = s.ALT
        ) TO '{out_parquet}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    rows_out = con.execute(f"SELECT COUNT(*) FROM read_parquet('{out_parquet}')").fetchone()[0]
    print(f"  rows: {rows_in:,} -> {rows_out:,}", file=sys.stderr)

    for col in out_cols:
        n_hit = con.execute(
            f'SELECT COUNT(*) FROM read_parquet(\'{out_parquet}\') WHERE "{col}" IS NOT NULL'
        ).fetchone()[0]
        pct = 100 * n_hit / max(rows_out, 1)
        print(f"  {col:<32}  non-null: {n_hit:>10,}  ({pct:>5.1f}%)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

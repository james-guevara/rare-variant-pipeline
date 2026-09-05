#!/usr/bin/env python3
"""Left-join per-gene constraint metrics onto the with_scores parquets.

Joins (all on SYMBOL == constraint.gene):

  gnomAD v4.1 constraint metrics (per-gene aggregates from per-transcript table):
    gnomad_lof_pLI_worst        = MAX(lof.pLI) across all transcripts
    gnomad_lof_LOEUF_worst      = MIN(lof.oe_ci.upper) across all transcripts
    gnomad_lof_pLI_mane         = lof.pLI where mane_select is TRUE
    gnomad_lof_LOEUF_mane       = lof.oe_ci.upper where mane_select is TRUE
    gnomad_lof_pLI_canonical    = lof.pLI where canonical is TRUE
    gnomad_lof_LOEUF_canonical  = lof.oe_ci.upper where canonical is TRUE

  GeneBayes (one row per gene already; joins by version-normalized Ensembl gene ID):
    genebayes_obs_lof, genebayes_exp_lof, genebayes_prior_mean,
    genebayes_post_mean, genebayes_post_lower_95, genebayes_post_upper_95

Join keys: SYMBOL for gnomAD constraint; version-normalized Gene/ensg for
GeneBayes. Variants without a matching key get NULL for the corresponding
columns.

Usage:
    python join_gene_constraint.py --cohort ssc --chrom chr22 --resources config/resources.json
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
    in_parquet = base / "with_scores" / f"{args.chrom}.parquet"
    out_dir = base / "with_gene_constraint"
    if not args.output:
        out_dir.mkdir(parents=True, exist_ok=True)
    out_parquet = out_dir / f"{args.chrom}.parquet"

    # Explicit-path overrides (used by the Nextflow POSTPROCESS subworkflow).
    if args.input:
        in_parquet = Path(args.input)
    if args.output:
        out_parquet = Path(args.output)
        out_parquet.parent.mkdir(parents=True, exist_ok=True)

    constraint_tsv = cfg["gnomad_constraint_tsv"]
    genebayes_tsv = cfg["genebayes_tsv"]
    print(f"[{args.cohort} {args.chrom}] join_gene_constraint: {in_parquet} -> {out_parquet}", file=sys.stderr)
    print(f"  gnomAD constraint: {constraint_tsv}", file=sys.stderr)
    print(f"  GeneBayes:         {genebayes_tsv}", file=sys.stderr)

    con = duckdb.connect()
    rows_in = con.execute(f"SELECT COUNT(*) FROM read_parquet('{in_parquet}')").fetchone()[0]

    # gnomAD constraint: aggregate per gene
    con.execute(f"""
        CREATE TABLE _gnomad AS
        SELECT gene,
               MAX(TRY_CAST("lof.pLI" AS DOUBLE)) AS gnomad_lof_pLI_worst,
               MIN(TRY_CAST("lof.oe_ci.upper" AS DOUBLE)) AS gnomad_lof_LOEUF_worst,
               MAX(CASE WHEN mane_select = 'true' THEN TRY_CAST("lof.pLI" AS DOUBLE) END) AS gnomad_lof_pLI_mane,
               MIN(CASE WHEN mane_select = 'true' THEN TRY_CAST("lof.oe_ci.upper" AS DOUBLE) END) AS gnomad_lof_LOEUF_mane,
               MAX(CASE WHEN canonical = 'true'    THEN TRY_CAST("lof.pLI" AS DOUBLE) END) AS gnomad_lof_pLI_canonical,
               MIN(CASE WHEN canonical = 'true'    THEN TRY_CAST("lof.oe_ci.upper" AS DOUBLE) END) AS gnomad_lof_LOEUF_canonical
        FROM read_csv('{constraint_tsv}', delim='\t', header=true, null_padding=true, all_varchar=true)
        GROUP BY gene
    """)
    n_g = con.execute("SELECT COUNT(*) FROM _gnomad").fetchone()[0]
    print(f"  gnomAD genes: {n_g:,}", file=sys.stderr)

    # GeneBayes: already per-gene; join on ensg (Ensembl gene ID)
    con.execute(f"""
        CREATE TABLE _gb AS
        SELECT regexp_replace(ensg, '\\.[0-9]+$', '') AS stable_ensg,
               TRY_CAST(obs_lof AS DOUBLE)        AS genebayes_obs_lof,
               TRY_CAST(exp_lof AS DOUBLE)        AS genebayes_exp_lof,
               TRY_CAST(prior_mean AS DOUBLE)     AS genebayes_prior_mean,
               TRY_CAST(post_mean AS DOUBLE)      AS genebayes_post_mean,
               TRY_CAST(post_lower_95 AS DOUBLE)  AS genebayes_post_lower_95,
               TRY_CAST(post_upper_95 AS DOUBLE)  AS genebayes_post_upper_95
        FROM read_csv('{genebayes_tsv}', delim='\t', header=true, all_varchar=true)
    """)
    n_gb = con.execute("SELECT COUNT(*) FROM _gb").fetchone()[0]
    print(f"  GeneBayes genes: {n_gb:,}", file=sys.stderr)

    con.execute(f"""
        COPY (
            SELECT c.*,
                   g.gnomad_lof_pLI_worst,
                   g.gnomad_lof_LOEUF_worst,
                   g.gnomad_lof_pLI_mane,
                   g.gnomad_lof_LOEUF_mane,
                   g.gnomad_lof_pLI_canonical,
                   g.gnomad_lof_LOEUF_canonical,
                   b.genebayes_obs_lof,
                   b.genebayes_exp_lof,
                   b.genebayes_prior_mean,
                   b.genebayes_post_mean,
                   b.genebayes_post_lower_95,
                   b.genebayes_post_upper_95
            FROM read_parquet('{in_parquet}') c
            LEFT JOIN _gnomad g ON c.SYMBOL = g.gene
            LEFT JOIN _gb b
              ON regexp_replace(c.Gene, '\\.[0-9]+$', '') = b.stable_ensg
        ) TO '{out_parquet}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    rows_out = con.execute(f"SELECT COUNT(*) FROM read_parquet('{out_parquet}')").fetchone()[0]
    print(f"  rows: {rows_in:,} -> {rows_out:,}", file=sys.stderr)

    # Hit rates
    for col in ["gnomad_lof_pLI_worst", "gnomad_lof_LOEUF_worst",
                "gnomad_lof_pLI_mane", "gnomad_lof_LOEUF_mane",
                "genebayes_post_mean"]:
        n_hit = con.execute(
            f'SELECT COUNT(*) FROM read_parquet(\'{out_parquet}\') WHERE "{col}" IS NOT NULL'
        ).fetchone()[0]
        pct = 100 * n_hit / max(rows_out, 1)
        print(f"  {col:<35} non-null: {n_hit:>10,} ({pct:>5.1f}%)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

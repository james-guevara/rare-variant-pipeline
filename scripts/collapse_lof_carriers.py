#!/usr/bin/env python3
"""Create explicit variant- and gene-level LoF carrier burden tables."""

import argparse
from pathlib import Path

import duckdb


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--sample-gene-output", required=True, type=Path)
    parser.add_argument("--sample-burden-output", required=True, type=Path)
    args = parser.parse_args()

    con = duckdb.connect()
    source = str(args.input).replace("'", "''")
    sample_gene = str(args.sample_gene_output).replace("'", "''")
    sample_burden = str(args.sample_burden_output).replace("'", "''")
    con.execute(f"""
        COPY (
            SELECT sample_id, Gene, SYMBOL, lof_tier,
                   COUNT(DISTINCT record_id) AS n_variants,
                   string_agg(DISTINCT record_id, ',' ORDER BY record_id) AS record_ids
            FROM read_parquet('{source}')
            GROUP BY sample_id, Gene, SYMBOL, lof_tier
            ORDER BY sample_id, lof_tier, Gene
        ) TO '{sample_gene}' (FORMAT CSV, DELIMITER '\t', HEADER)
    """)
    con.execute(f"""
        COPY (
            SELECT sample_id AS SAMPLE,
                   COUNT(*) FILTER (WHERE lof_tier = 'lof_t1') AS lof_t1_genes,
                   COALESCE(SUM(n_variants) FILTER (WHERE lof_tier = 'lof_t1'), 0) AS lof_t1_variants,
                   COUNT(*) FILTER (WHERE lof_tier = 'lof_t2') AS lof_t2_genes,
                   COALESCE(SUM(n_variants) FILTER (WHERE lof_tier = 'lof_t2'), 0) AS lof_t2_variants,
                   COUNT(*) AS any_tier_genes,
                   SUM(n_variants) AS any_tier_variants
            FROM (
                SELECT sample_id, Gene, lof_tier,
                       COUNT(DISTINCT record_id) AS n_variants
                FROM read_parquet('{source}')
                GROUP BY sample_id, Gene, lof_tier
            )
            GROUP BY sample_id
            ORDER BY sample_id
        ) TO '{sample_burden}' (FORMAT CSV, DELIMITER '\t', HEADER)
    """)
    rows = con.execute(
        f"SELECT COUNT(*) FROM read_csv_auto('{sample_gene}', delim='\\t', header=true)"
    ).fetchone()[0]
    samples = con.execute(
        f"SELECT COUNT(*) FROM read_csv_auto('{sample_burden}', delim='\\t', header=true)"
    ).fetchone()[0]
    print("sample_gene_rows={:,}".format(rows))
    print("samples={:,}".format(samples))


if __name__ == "__main__":
    main()

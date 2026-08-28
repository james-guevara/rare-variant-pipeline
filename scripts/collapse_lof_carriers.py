#!/usr/bin/env python3
"""Create explicit variant- and gene-level LoF carrier burden tables."""

import argparse
import re
from pathlib import Path

import duckdb


def identifier(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", value):
        raise ValueError(f"unsafe column identifier: {value}")
    return '"{}"'.format(value)


def write_burdens(con, source, sample_gene, sample_burden, where_clause):
    con.execute(f"""
        COPY (
            SELECT sample_id, Gene, SYMBOL, lof_tier,
                   COUNT(DISTINCT record_id) AS n_variants,
                   string_agg(DISTINCT record_id, ',' ORDER BY record_id) AS record_ids
            FROM read_parquet('{source}')
            {where_clause}
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
                {where_clause}
                GROUP BY sample_id, Gene, lof_tier
            )
            GROUP BY sample_id
            ORDER BY sample_id
        ) TO '{sample_burden}' (FORMAT CSV, DELIMITER '\t', HEADER)
    """)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--sample-gene-output", required=True, type=Path)
    parser.add_argument("--sample-burden-output", required=True, type=Path)
    parser.add_argument("--eligibility-col")
    parser.add_argument("--burden-available-col")
    parser.add_argument("--sensitivity-sample-gene-output", type=Path)
    parser.add_argument("--sensitivity-sample-burden-output", type=Path)
    args = parser.parse_args()
    if bool(args.sensitivity_sample_gene_output) != bool(args.sensitivity_sample_burden_output):
        parser.error("both sensitivity outputs must be supplied together")
    if args.sensitivity_sample_gene_output and not args.eligibility_col:
        parser.error("sensitivity outputs require --eligibility-col")

    args.sample_gene_output.parent.mkdir(parents=True, exist_ok=True)
    args.sample_burden_output.parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    source = str(args.input).replace("'", "''")
    sample_gene = str(args.sample_gene_output).replace("'", "''")
    sample_burden = str(args.sample_burden_output).replace("'", "''")
    predicates = []
    if args.eligibility_col:
        predicates.append(f"COALESCE({identifier(args.eligibility_col)}, FALSE)")
    if args.burden_available_col:
        predicates.append(f"COALESCE({identifier(args.burden_available_col)}, FALSE)")
    primary_where = "WHERE " + " AND ".join(predicates) if predicates else ""
    write_burdens(con, source, sample_gene, sample_burden, primary_where)

    if args.sensitivity_sample_gene_output:
        sensitivity_gene = str(args.sensitivity_sample_gene_output).replace("'", "''")
        sensitivity_burden = str(args.sensitivity_sample_burden_output).replace("'", "''")
        args.sensitivity_sample_gene_output.parent.mkdir(parents=True, exist_ok=True)
        args.sensitivity_sample_burden_output.parent.mkdir(parents=True, exist_ok=True)
        sensitivity_predicates = [
            f"NOT COALESCE({identifier(args.eligibility_col)}, FALSE)"
        ]
        if args.burden_available_col:
            sensitivity_predicates.append(
                f"COALESCE({identifier(args.burden_available_col)}, FALSE)"
            )
        write_burdens(
            con, source, sensitivity_gene, sensitivity_burden,
            "WHERE " + " AND ".join(sensitivity_predicates),
        )
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

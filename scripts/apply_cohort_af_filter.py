#!/usr/bin/env python3
"""Annotate carrier rows with cohort AF and emit a burden-eligible subset."""

import argparse
import re
from pathlib import Path

import duckdb


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--allele-summary", required=True, type=Path)
    parser.add_argument("--max-af", required=True, type=float)
    parser.add_argument("--annotated-output", required=True, type=Path)
    parser.add_argument("--eligible-output", required=True, type=Path)
    parser.add_argument(
        "--summary-prefix", default="genotype",
        help="allele-summary metric prefix (genotype or primary_genotype)",
    )
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", args.summary_prefix):
        parser.error("--summary-prefix must be a SQL-safe identifier")

    args.annotated_output.parent.mkdir(parents=True, exist_ok=True)
    args.eligible_output.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    prefix = args.summary_prefix
    con.execute(
        f"""
        CREATE TABLE annotated AS
        SELECT c.*,
               s.{prefix}_ac AS cohort_ac,
               s.{prefix}_an AS cohort_an,
               s.{prefix}_af AS cohort_af,
               s.{"primary_carrier_count" if prefix == "primary_genotype" else "carrier_count"} AS cohort_carrier_count,
               s.{"primary_hom_alt_count" if prefix == "primary_genotype" else "hom_alt_count"} AS cohort_hom_alt_count,
               COALESCE(s.{prefix}_af < ?, FALSE) AS cohort_af_eligible
        FROM read_parquet(?) c
        LEFT JOIN read_parquet(?) s USING (record_id)
        """,
        [args.max_af, str(args.input), str(args.allele_summary)],
    )
    missing = con.execute(
        "SELECT COUNT(*) FROM annotated WHERE cohort_af IS NULL"
    ).fetchone()[0]
    if missing:
        raise ValueError("{} carrier rows lack cohort AF".format(missing))

    con.execute(
        "COPY annotated TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
        [str(args.annotated_output)],
    )
    con.execute(
        "COPY (SELECT * FROM annotated WHERE cohort_af_eligible) "
        "TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
        [str(args.eligible_output)],
    )
    rows, alleles = con.execute(
        "SELECT COUNT(*), COUNT(DISTINCT record_id) FROM annotated"
    ).fetchone()
    kept_rows, kept_alleles = con.execute(
        "SELECT COUNT(*), COUNT(DISTINCT record_id) FROM annotated "
        "WHERE cohort_af_eligible"
    ).fetchone()
    print("rows={:,}->{:,}".format(rows, kept_rows))
    print("alleles={:,}->{:,}".format(alleles, kept_alleles))


if __name__ == "__main__":
    main()

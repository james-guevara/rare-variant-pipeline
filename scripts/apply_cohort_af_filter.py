#!/usr/bin/env python3
"""Annotate carrier rows with cohort AF and emit a burden-eligible subset."""

import argparse
from pathlib import Path

import duckdb


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--allele-summary", required=True, type=Path)
    parser.add_argument("--max-af", required=True, type=float)
    parser.add_argument("--annotated-output", required=True, type=Path)
    parser.add_argument("--eligible-output", required=True, type=Path)
    args = parser.parse_args()

    args.annotated_output.parent.mkdir(parents=True, exist_ok=True)
    args.eligible_output.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(
        """
        CREATE TABLE annotated AS
        SELECT c.*,
               s.genotype_ac AS cohort_ac,
               s.genotype_an AS cohort_an,
               s.genotype_af AS cohort_af,
               s.carrier_count AS cohort_carrier_count,
               s.hom_alt_count AS cohort_hom_alt_count,
               COALESCE(s.genotype_af < ?, FALSE) AS cohort_af_eligible
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

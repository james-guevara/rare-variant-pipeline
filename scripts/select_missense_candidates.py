#!/usr/bin/env python3
"""Retain dbNSFP candidates whose VEP-picked consequence is missense."""

import argparse
from pathlib import Path

import duckdb


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--picked", required=True, type=Path)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    connection = duckdb.connect()
    connection.read_csv(
        str(args.picked), delimiter="\t", header=True, all_varchar=True
    ).create_view("picked_input")
    connection.read_parquet(str(args.candidates)).create_view("candidate_input")
    connection.execute(
        """
        COPY (
          WITH picked AS (
            SELECT *,
                   regexp_replace(CHROM, '^chr', '') AS chrom_key,
                   TRY_CAST(POS AS BIGINT) AS pos_key
            FROM picked_input
          ), candidates AS (
            SELECT *,
                   regexp_replace(chrom, '^chr', '') AS chrom_key,
                   CAST(pos AS BIGINT) AS pos_key
            FROM candidate_input
          )
          SELECT
            p.* EXCLUDE (chrom_key, pos_key),
            c.candidate_genes,
            c.miss_n_flag,
            c.miss_tier,
            c.ClinPred_rankscore,
            c.AlphaMissense_rankscore,
            c.popEVE_converted_rankscore,
            c.MPC_rankscore
          FROM picked p
          JOIN candidates c
            ON p.chrom_key = c.chrom_key
           AND p.pos_key = c.pos_key
           AND p.REF = c.ref
           AND p.ALT = c.alt
          WHERE contains(p.Consequence, 'missense_variant')
            AND list_contains(
                  string_split(c.candidate_genes, ','),
                  regexp_replace(p.Gene, '\\.[0-9]+$', '')
                )
        ) TO ? (FORMAT PARQUET, COMPRESSION ZSTD)
        """,
        [str(args.output)],
    )

    total_candidates = connection.execute(
        "SELECT COUNT(*) FROM read_parquet(?)", [str(args.candidates)]
    ).fetchone()[0]
    selected = connection.execute(
        "SELECT COUNT(*) FROM read_parquet(?)", [str(args.output)]
    ).fetchone()[0]
    print(f"observed_candidates={total_candidates:,}")
    print(f"selected_missense={selected:,}")
    print(f"rejected_selected_consequence={total_candidates-selected:,}")
    for tier, count in connection.execute(
        "SELECT miss_tier, COUNT(*) FROM read_parquet(?) GROUP BY 1 ORDER BY 1",
        [str(args.output)],
    ).fetchall():
        print(f"{tier}={count:,}")


if __name__ == "__main__":
    main()

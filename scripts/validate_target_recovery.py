#!/usr/bin/env python3
"""Check whether target BEDs cover known tier-eligible annotated variants."""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    args = parser.parse_args()

    connection = duckdb.connect()
    connection.execute(
        """
        CREATE TEMP TABLE target AS
        SELECT column0 AS chrom, column1::BIGINT AS start0, column2::BIGINT AS end0
        FROM read_csv(?, delim='\t', header=false, columns={
          'column0': 'VARCHAR', 'column1': 'BIGINT', 'column2': 'BIGINT'
        })
        """,
        [str(args.bed)],
    )
    query = """
        WITH scored AS (
          SELECT *,
            CAST(COALESCE(ClinPred_rankscore >= 0.4298, FALSE) AS INTEGER) +
            CAST(COALESCE(AlphaMissense_rankscore >= 0.9603, FALSE) AS INTEGER) +
            CAST(COALESCE(popEVE_converted_rankscore >= 0.9209, FALSE) AS INTEGER) +
            CAST(COALESCE(MPC_rankscore >= 0.8947, FALSE) AS INTEGER) AS n_flag
          FROM read_parquet(?)
        ), eligible AS (
          SELECT DISTINCT "#CHROM" AS chrom, POS::BIGINT AS pos, REF, ALT, Gene,
            CASE
              WHEN LoF = 'HC' AND genebayes_post_mean >= 0.18 THEN 'lof_t1'
              WHEN LoF = 'HC' AND genebayes_post_mean >= 0.03 THEN 'lof_t2'
              WHEN Consequence LIKE '%missense_variant%' AND n_flag = 4 THEN 'miss_t1'
              WHEN Consequence LIKE '%missense_variant%' AND n_flag = 3 THEN 'miss_t2'
              WHEN Consequence LIKE '%missense_variant%' AND n_flag = 2 THEN 'miss_t3'
              WHEN Consequence LIKE '%missense_variant%' AND n_flag = 1 THEN 'miss_t4'
            END AS tier
          FROM scored
          WHERE (LoF = 'HC' AND genebayes_post_mean >= 0.03)
             OR (Consequence LIKE '%missense_variant%' AND n_flag >= 1)
        ), covered AS (
          SELECT e.*,
            EXISTS (
              SELECT 1 FROM target t
              WHERE t.chrom = e.chrom AND e.pos - 1 < t.end0 AND e.pos > t.start0
            ) AS is_covered
          FROM eligible e
        )
        SELECT tier, count(*) AS known_events,
          count(*) FILTER (WHERE is_covered) AS covered_events,
          count(*) FILTER (WHERE NOT is_covered) AS missed_events
        FROM covered
        GROUP BY tier
        ORDER BY tier
    """
    rows = connection.execute(query, [str(args.variants)]).fetchall()
    print("tier\tknown_events\tcovered_events\tmissed_events")
    missed = 0
    for row in rows:
        print(*row, sep="\t")
        missed += row[3]
    return 1 if missed else 0


if __name__ == "__main__":
    raise SystemExit(main())

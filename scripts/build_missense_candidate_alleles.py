#!/usr/bin/env python3
"""Build exact dbNSFP missense-tier candidates and optionally intersect observed alleles."""

import argparse
from pathlib import Path

import duckdb


THRESHOLDS = {
    "ClinPred_rankscore": 0.4298,
    "AlphaMissense_rankscore": 0.9603,
    "popEVE_converted_rankscore": 0.9209,
    "MPC_rankscore": 0.8947,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dbnsfp", required=True, type=Path)
    parser.add_argument("--genebayes", type=Path)
    parser.add_argument("--all-genes", action="store_true")
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--min-post-mean", default=0.03, type=float)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--observed-alleles", type=Path)
    parser.add_argument("--observed-vcf", type=Path)
    parser.add_argument("--observed-output", type=Path)
    for column, default in THRESHOLDS.items():
        parser.add_argument(
            "--" + column.replace("_", "-").lower() + "-min",
            dest=column,
            type=float,
            default=default,
        )
    for tier, default in (("miss_t1", 4), ("miss_t2", 3), ("miss_t3", 2), ("miss_t4", 1)):
        parser.add_argument("--" + tier.replace("_", "-") + "-pass-count", type=int, default=default)
    args = parser.parse_args()
    thresholds = {column: getattr(args, column) for column in THRESHOLDS}
    if not all(0 <= value <= 1 for value in thresholds.values()):
        parser.error("missense score thresholds must be between 0 and 1")
    tier_pass_counts = {tier: getattr(args, tier + "_pass_count") for tier in (
        "miss_t1", "miss_t2", "miss_t3", "miss_t4"
    )}
    if set(tier_pass_counts.values()) != {1, 2, 3, 4}:
        parser.error("missense tier pass counts must use each count from 1 through 4 once")
    observed_inputs = sum(bool(value) for value in (args.observed_alleles, args.observed_vcf))
    if bool(observed_inputs) != bool(args.observed_output) or observed_inputs > 1:
        parser.error("use --observed-output with exactly one of --observed-alleles or --observed-vcf")
    if args.all_genes == bool(args.genebayes):
        parser.error("use exactly one of --all-genes or --genebayes")

    score_values = ", ".join(
        'TRY_CAST("{}" AS DOUBLE) AS "{}"'.format(column, column)
        for column in thresholds
    )
    flags = " + ".join(
        'CAST(COALESCE("{}" >= {}, FALSE) AS INTEGER)'.format(column, threshold)
        for column, threshold in thresholds.items()
    )
    score_max = ", ".join(
        'MAX("{0}") AS "{0}"'.format(column) for column in thresholds
    )
    chrom_without_prefix = args.chrom.removeprefix("chr")
    con = duckdb.connect()
    target_cte = ""
    target_join = ""
    parameters = []
    if args.genebayes:
        target_cte = """
        target_genes AS (
            SELECT regexp_replace(ensg, '\\.[0-9]+$', '') AS Gene
            FROM read_csv(?, delim='\\t', header=true, all_varchar=true)
            WHERE chrom = ? AND TRY_CAST(post_mean AS DOUBLE) >= ?
        ),
        """
        target_join = "JOIN target_genes USING (Gene)"
        parameters.extend([str(args.genebayes), args.chrom, args.min_post_mean])
    parameters.extend([str(args.dbnsfp), chrom_without_prefix, args.chrom])
    con.execute(
        """
        CREATE TABLE candidates AS
        WITH {target_cte} raw AS (
            SELECT CAST("pos(1-based)" AS BIGINT) AS POS, ref AS REF, alt AS ALT,
                   regexp_replace(gene_id, '\\.[0-9]+$', '') AS Gene,
                   {score_values}
            FROM read_parquet(?),
                 UNNEST(string_split(Ensembl_geneid, ';')) AS genes(gene_id)
            WHERE "#chr" = ?
        ), scored AS (
            SELECT *, ({flags}) AS miss_n_flag
            FROM raw
        )
        SELECT ? AS CHROM, POS, REF, ALT,
               string_agg(DISTINCT Gene, ',' ORDER BY Gene) AS candidate_genes,
               MAX(miss_n_flag) AS miss_n_flag,
               CASE MAX(miss_n_flag)
                   WHEN {miss_t1} THEN 'miss_t1'
                   WHEN {miss_t2} THEN 'miss_t2'
                   WHEN {miss_t3} THEN 'miss_t3'
                   WHEN {miss_t4} THEN 'miss_t4'
               END AS miss_tier,
               {score_max}
        FROM scored {target_join}
        WHERE miss_n_flag >= 1
        GROUP BY POS, REF, ALT
        """.format(target_cte=target_cte, target_join=target_join,
                   score_values=score_values, flags=flags, score_max=score_max,
                   **tier_pass_counts),
        parameters,
    )
    if args.output:
        con.execute(
            "COPY candidates TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
            [str(args.output)],
        )
    count = con.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
    print("candidate_alleles={:,}".format(count))
    for tier, n in con.execute(
        "SELECT miss_tier, COUNT(*) FROM candidates GROUP BY 1 ORDER BY 1"
    ).fetchall():
        print("{}={:,}".format(tier, n))

    if args.observed_alleles:
        observed_path = str(args.observed_alleles).replace("'", "''")
        observed_output = str(args.observed_output).replace("'", "''")
        con.execute(
            f"""
            COPY (
                SELECT o.*, c.* EXCLUDE (CHROM, POS, REF, ALT)
                FROM read_parquet('{observed_path}') o
                JOIN candidates c
                  ON o.chrom = c.CHROM AND o.pos = c.POS
                 AND o.ref = c.REF AND o.alt = c.ALT
            ) TO '{observed_output}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
        observed = con.execute(
            "SELECT COUNT(*) FROM read_parquet(?)", [str(args.observed_output)]
        ).fetchone()[0]
        print("observed_candidate_alleles={:,}".format(observed))
    elif args.observed_vcf:
        metadata_lines = 0
        with args.observed_vcf.open() as handle:
            for line in handle:
                if line.startswith("##"):
                    metadata_lines += 1
                elif line.startswith("#CHROM"):
                    break
                else:
                    raise ValueError("VCF header is missing #CHROM")
        vcf_path = str(args.observed_vcf).replace("'", "''")
        con.execute(
            f"""
            CREATE VIEW normalized_observed AS
            SELECT regexp_replace("#CHROM", '^chr', '') AS chrom,
                   CAST(POS AS BIGINT) AS pos, REF AS ref, ALT AS alt
            FROM read_csv('{vcf_path}', delim='\t', header=true,
                          skip={metadata_lines}, all_varchar=true)
            """
        )
        con.execute(
            """
            COPY (
                SELECT o.chrom, o.pos, o.ref, o.alt,
                       c.* EXCLUDE (CHROM, POS, REF, ALT)
                FROM normalized_observed o
                JOIN candidates c
                  ON o.chrom = regexp_replace(c.CHROM, '^chr', '')
                 AND o.pos = c.POS AND o.ref = c.REF AND o.alt = c.ALT
            ) TO ? (FORMAT PARQUET, COMPRESSION ZSTD)
            """,
            [str(args.observed_output)],
        )
        observed = con.execute(
            "SELECT COUNT(*) FROM read_parquet(?)", [str(args.observed_output)]
        ).fetchone()[0]
        print("observed_candidate_alleles={:,}".format(observed))


if __name__ == "__main__":
    main()

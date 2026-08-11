#!/usr/bin/env python3
"""
Filter rare LOF (loss-of-function) variants using DuckDB streaming.

LOF variants include: stop_gained, frameshift, splice_donor, splice_acceptor,
start_lost, stop_lost, transcript_ablation.

Can filter by LOFTEE confidence (HC = high confidence, LC = low confidence).

Usage:
    python filter_lof.py input.tsv.gz -o filtered.tsv --af 0.01
    python filter_lof.py input.tsv.gz -o filtered.tsv --af 0.01 --loftee hc
    python filter_lof.py input.tsv.gz -o filtered.tsv --af 0.01 --loftee hc --canonical
"""

import argparse
from pathlib import Path

import duckdb

GENE_SETS_DIR = Path("/expanse/projects/sebat1/s3/data/sebat/nf_rare_spark_wes/resources/gene_sets")


def load_gene_set(genes_arg: str) -> set[str] | None:
    if genes_arg is None:
        return None

    gene_path = Path(genes_arg)
    if not gene_path.exists():
        gene_path = GENE_SETS_DIR / genes_arg

    if gene_path.exists():
        with open(gene_path) as f:
            genes = {line.strip() for line in f if line.strip()}
        print(f"Loaded {len(genes)} genes from {gene_path}", file=__import__('sys').stderr)
        return genes
    else:
        genes = {g.strip() for g in genes_arg.split(",") if g.strip()}
        print(f"Using {len(genes)} genes from command line", file=__import__('sys').stderr)
        return genes


def build_where_clause(args, gene_set: set[str] | None) -> str:
    conditions = []

    # LOF consequence filter
    lof_consequences = [
        "stop_gained",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "start_lost",
        "stop_lost",
        "transcript_ablation"
    ]
    lof_condition = " OR ".join([f"Consequence LIKE '%{c}%'" for c in lof_consequences])
    conditions.append(f"({lof_condition})")

    # LOFTEE filter
    if args.loftee:
        if args.loftee.lower() == "hc":
            conditions.append("LoF = 'HC'")
        elif args.loftee.lower() == "lc":
            conditions.append("LoF IN ('HC', 'LC')")
        # else "any": no LOFTEE filter

    # Exclude specific LoF_filter values
    if args.exclude_lof_filters:
        filters_to_exclude = [f.strip() for f in args.exclude_lof_filters.split(",")]
        for f in filters_to_exclude:
            conditions.append(f"(LoF_filter IS NULL OR LoF_filter NOT LIKE '%{f}%')")

    # Transcript filters
    if args.canonical:
        conditions.append("CANONICAL = 'YES'")
    if args.mane:
        conditions.append("(MANE_SELECT IS NOT NULL AND MANE_SELECT != '')")

    # AF filters
    if args.af is not None:
        conditions.append(f"(TRY_CAST(AF AS DOUBLE) < {args.af} OR AF IS NULL)")
    if args.gnomad_af is not None:
        conditions.append(
            f'(TRY_CAST("gnomAD4.1_joint_AF" AS DOUBLE) < {args.gnomad_af} OR "gnomAD4.1_joint_AF" IS NULL)'
        )

    # Gene set filter
    if gene_set:
        gene_list = ", ".join(f"'{g}'" for g in gene_set)
        conditions.append(f"SYMBOL IN ({gene_list})")

    # Gene constraint filters
    if args.pli is not None:
        conditions.append(f"TRY_CAST(gnomad_lof_pLI_canonical AS DOUBLE) > {args.pli}")
    if args.loeuf is not None:
        conditions.append(f"TRY_CAST(gnomad_lof_oe_ci_upper_canonical AS DOUBLE) < {args.loeuf}")
    if args.genebayes is not None:
        conditions.append(f"TRY_CAST(genebayes_post_mean AS DOUBLE) > {args.genebayes}")

    # QC filters
    if args.gq is not None:
        conditions.append(f"TRY_CAST(GQ AS INTEGER) >= {args.gq}")
    if args.dp is not None:
        conditions.append(f"TRY_CAST(DP AS INTEGER) >= {args.dp}")
    if args.qual is not None:
        conditions.append(f"TRY_CAST(QUAL AS DOUBLE) >= {args.qual}")

    # Region filters
    if args.exclude_segdups:
        conditions.append("(segDups IS NULL OR segDups = '')")
    if args.exclude_simple_repeats:
        conditions.append("(simpleRepeats IS NULL OR simpleRepeats = '')")

    # Multiallelic filter
    if args.exclude_multiallelic:
        conditions.append("ALT NOT LIKE '%,%'")

    return " AND ".join(conditions)


def main():
    parser = argparse.ArgumentParser(description="Filter rare LOF variants")
    parser.add_argument("input", help="Input TSV file (use - for stdin)")
    parser.add_argument("--output", "-o", required=True, help="Output TSV file (use - for stdout)")
    
    # LOF-specific
    parser.add_argument("--loftee", choices=["hc", "lc", "any"], 
                        help="LOFTEE filter: hc (high confidence only), lc (hc + low confidence), any (no filter)")
    parser.add_argument("--exclude-lof-filters", 
                        help="Comma-separated LoF_filter values to exclude (e.g., END_TRUNC,SINGLE_EXON)")

    # Transcript filters
    parser.add_argument("--canonical", action="store_true",
                        help="Only include variants on canonical transcripts (CANONICAL='YES')")
    parser.add_argument("--mane", action="store_true",
                        help="Only include variants on MANE Select transcripts")
    
    # AF filters
    parser.add_argument("--af", type=float, help="Cohort AF threshold")
    parser.add_argument("--gnomad-af", type=float, help="gnomAD AF threshold")

    # Gene filters
    parser.add_argument("--genes", "-g", help="Gene set file or comma-separated list")

    # Gene constraint filters
    parser.add_argument("--pli", type=float, help="pLI threshold")
    parser.add_argument("--loeuf", type=float, help="LOEUF threshold")
    parser.add_argument("--genebayes", type=float, help="GeneBayes threshold")

    # QC filters
    parser.add_argument("--gq", type=int, help="Min GQ")
    parser.add_argument("--dp", type=int, help="Min DP")
    parser.add_argument("--qual", type=float, help="Min QUAL")

    # Region filters
    parser.add_argument("--exclude-segdups", action="store_true",
                        help="Exclude variants in segmental duplications")
    parser.add_argument("--exclude-simple-repeats", action="store_true",
                        help="Exclude variants in simple repeats")

    # Multiallelic filter
    parser.add_argument("--exclude-multiallelic", action="store_true",
                        help="Exclude multiallelic sites (ALT contains comma)")

    args = parser.parse_args()

    gene_set = load_gene_set(args.genes)
    where_clause = build_where_clause(args, gene_set)

    input_path = "/dev/stdin" if args.input == "-" else args.input
    output_path = "/dev/stdout" if args.output == "-" else args.output

    import sys
    if args.output != "-":
        print(f"Input: {args.input}", file=sys.stderr)
        print(f"Output: {args.output}", file=sys.stderr)

    query = f"""
        COPY (
            SELECT * 
            FROM read_csv_auto('{input_path}', delim='\t', header=true, null_padding=true, all_varchar=true)
            WHERE {where_clause}
        ) TO '{output_path}' (HEADER, DELIMITER '\t', QUOTE '')
    """

    duckdb.sql(query)
    
    if args.output != "-":
        print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()

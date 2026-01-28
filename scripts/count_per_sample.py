#!/usr/bin/env python3
"""
Count variants per sample from a filtered TSV file.

Usage:
    python count_per_sample.py filtered.tsv -o counts.tsv
    python count_per_sample.py - -o - < filtered.tsv
"""

import argparse
import sys
import duckdb


def main():
    parser = argparse.ArgumentParser(description="Count variants per sample")
    parser.add_argument("input", help="Input TSV file (use - for stdin)")
    parser.add_argument("--output", "-o", required=True, help="Output TSV file (use - for stdout)")

    args = parser.parse_args()

    input_path = "/dev/stdin" if args.input == "-" else args.input
    output_path = "/dev/stdout" if args.output == "-" else args.output

    if args.output != "-":
        print(f"Input: {args.input}", file=sys.stderr)
        print(f"Output: {args.output}", file=sys.stderr)

    query = f"""
        COPY (
            SELECT SAMPLE, COUNT(*) as n_variants
            FROM read_csv_auto('{input_path}', delim='\t', header=true, all_varchar=true)
            WHERE carrier = '1'
            GROUP BY SAMPLE
            ORDER BY n_variants DESC
        ) TO '{output_path}' (HEADER, DELIMITER '\t', QUOTE '')
    """

    duckdb.sql(query)

    if args.output != "-":
        print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()

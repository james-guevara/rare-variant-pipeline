#!/usr/bin/env python3
"""
Convert pipeline output TSV(.gz) files to Parquet using DuckDB.

Uses all_varchar=true to avoid type inference issues with genomic data
(mixed types like GT=0/1, AD=10,25, missing=.).

Usage:
    # Convert a single file
    python tsv_to_parquet.py input.tsv.gz -o output.parquet

    # Convert all .tsv.gz files in a directory (output to same dir)
    python tsv_to_parquet.py /path/to/indexed/

    # Convert to a different output directory
    python tsv_to_parquet.py /path/to/indexed/ -o /path/to/parquet/

    # Convert specific chromosomes
    python tsv_to_parquet.py /path/to/indexed/ --chroms chr1 chr2 chr22
"""
import duckdb
import argparse
import sys
import time
from pathlib import Path


def convert_file(input_path: Path, output_path: Path):
    """Convert a single TSV(.gz) to Parquet via DuckDB."""
    t0 = time.time()
    print(f"  {input_path.name} -> {output_path.name}", end="", flush=True, file=sys.stderr)

    duckdb.sql(f"""
        COPY (
            SELECT * FROM read_csv_auto(
                '{input_path}',
                delim='\t',
                header=true,
                null_padding=true,
                all_varchar=true
            )
        ) TO '{output_path}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)

    size_in = input_path.stat().st_size / (1024 * 1024)
    size_out = output_path.stat().st_size / (1024 * 1024)
    elapsed = time.time() - t0
    print(f"  ({size_in:.0f}MB -> {size_out:.0f}MB, {elapsed:.1f}s)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Convert pipeline output TSV(.gz) files to Parquet."
    )
    parser.add_argument("input", help="Input TSV(.gz) file or directory containing them")
    parser.add_argument("-o", "--output",
                        help="Output parquet file or directory (default: same as input)")
    parser.add_argument("--chroms", nargs="+",
                        help="Only convert these chromosomes (e.g., chr1 chr22)")
    parser.add_argument("--glob", default="*.tsv.gz",
                        help="Glob pattern for finding files in input directory (default: *.tsv.gz)")
    args = parser.parse_args()

    input_path = Path(args.input)

    # Single file mode
    if input_path.is_file():
        if args.output:
            output_path = Path(args.output)
            if output_path.suffix != ".parquet":
                output_path = output_path / input_path.stem.split(".")[0] + ".merged.parquet"
        else:
            # Same directory, swap extension
            stem = input_path.name
            for suffix in [".tsv.gz", ".tsv", ".gz"]:
                if stem.endswith(suffix):
                    stem = stem[: -len(suffix)]
                    break
            output_path = input_path.parent / f"{stem}.parquet"

        output_path.parent.mkdir(parents=True, exist_ok=True)

        if output_path.exists():
            print(f"  SKIP {input_path.name} (parquet already exists)", file=sys.stderr)
        else:
            convert_file(input_path, output_path)
        return

    # Directory mode
    if not input_path.is_dir():
        sys.exit(f"ERROR: Not a file or directory: {input_path}")

    output_dir = Path(args.output) if args.output else input_path

    files = sorted(input_path.glob(args.glob))
    if not files:
        sys.exit(f"ERROR: No files matching '{args.glob}' in {input_path}")

    # Filter by chromosome if specified
    if args.chroms:
        chroms = set(args.chroms)
        files = [f for f in files if any(c in f.stem for c in chroms)]
        if not files:
            sys.exit(f"ERROR: No files matching chromosomes {args.chroms}")

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Converting {len(files)} files: {input_path} -> {output_dir}", file=sys.stderr)
    total_start = time.time()
    converted = 0

    for f in files:
        stem = f.name
        for suffix in [".tsv.gz", ".tsv", ".gz"]:
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        out_file = output_dir / f"{stem}.parquet"

        if out_file.exists():
            print(f"  SKIP {f.name} (parquet already exists)", file=sys.stderr)
            continue

        convert_file(f, out_file)
        converted += 1

    elapsed = time.time() - total_start
    print(f"\nDone. Converted {converted}, skipped {len(files) - converted} in {elapsed:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()

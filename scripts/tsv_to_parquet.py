#!/usr/bin/env python3
"""
Convert pipeline output TSV(.gz) files to Parquet.

Usage:
    # Convert a single file
    python tsv_to_parquet.py input.tsv.gz -o output.parquet

    # Convert all .tsv.gz files in a directory
    python tsv_to_parquet.py /path/to/indexed/ -o /path/to/parquet/

    # Convert specific chromosomes
    python tsv_to_parquet.py /path/to/indexed/ -o /path/to/parquet/ \
        --chroms chr1 chr2 chr22
"""
import polars as pl
import argparse
import sys
import time
from pathlib import Path


def convert_file(input_path: Path, output_path: Path):
    """Convert a single TSV(.gz) to Parquet."""
    t0 = time.time()
    print(f"  {input_path.name} -> {output_path.name}", end="", flush=True, file=sys.stderr)

    lf = pl.scan_csv(str(input_path), separator="\t", infer_schema_length=10000)
    lf.sink_parquet(str(output_path))

    size_in = input_path.stat().st_size / (1024 * 1024)
    size_out = output_path.stat().st_size / (1024 * 1024)
    elapsed = time.time() - t0
    print(f"  ({size_in:.0f}MB -> {size_out:.0f}MB, {elapsed:.1f}s)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Convert pipeline output TSV(.gz) files to Parquet."
    )
    parser.add_argument("input", help="Input TSV(.gz) file or directory containing them")
    parser.add_argument("-o", "--output", required=True,
                        help="Output parquet file or directory")
    parser.add_argument("--chroms", nargs="+",
                        help="Only convert these chromosomes (e.g., chr1 chr22)")
    parser.add_argument("--glob", default="*.tsv.gz",
                        help="Glob pattern for finding files in input directory (default: *.tsv.gz)")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    # Single file mode
    if input_path.is_file():
        if output_path.suffix != ".parquet":
            output_path = output_path / input_path.with_suffix(".parquet").name
        output_path.parent.mkdir(parents=True, exist_ok=True)
        convert_file(input_path, output_path)
        return

    # Directory mode
    if not input_path.is_dir():
        sys.exit(f"ERROR: Not a file or directory: {input_path}")

    files = sorted(input_path.glob(args.glob))
    if not files:
        sys.exit(f"ERROR: No files matching '{args.glob}' in {input_path}")

    # Filter by chromosome if specified
    if args.chroms:
        chroms = set(args.chroms)
        files = [f for f in files if any(c in f.stem for c in chroms)]
        if not files:
            sys.exit(f"ERROR: No files matching chromosomes {args.chroms}")

    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Converting {len(files)} files: {input_path} -> {output_path}", file=sys.stderr)
    total_start = time.time()

    for f in files:
        # .tsv.gz -> .parquet, .merged.tsv.gz -> .merged.parquet
        stem = f.name
        for suffix in [".tsv.gz", ".tsv", ".gz"]:
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        out_file = output_path / f"{stem}.parquet"
        convert_file(f, out_file)

    elapsed = time.time() - total_start
    print(f"\nDone. Converted {len(files)} files in {elapsed:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()

"""
Concatenate per-chromosome files into genome-wide tables.

Combines multiple chromosome-specific TSV files (e.g., chr1, chr2, ..., chrX, chrY)
into a single genome-wide file, sorted by chromosome order.
"""
import polars as pl
import argparse
import sys
import re
from pathlib import Path

def chrom_key(p: Path):
    """Sort key for chromosome files (1-22, X=23, Y=24)."""
    m = re.search(r"chr([0-9]+|X|Y)\b", p.name)
    if not m:
        return (999, p.name)
    c = m.group(1)
    if c.isdigit():
        return (int(c), p.name)
    return (23 if c == "X" else 24, p.name)


def main(pattern: str, out_path: str):
    files = sorted(Path(".").glob(pattern), key=chrom_key)
    if not files:
        sys.exit(f"ERROR: no files matched pattern: {pattern}")

    lfs = [pl.scan_csv(str(f), separator="\t", infer_schema_length=0) for f in files]
    lf = pl.concat(lfs, how="vertical_relaxed")

    lf.collect(engine="streaming").write_csv(out_path, separator="\t")
    print(f"Wrote {out_path} from {len(files)} files", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="out/chr*.family_merged.tsv.gz")
    ap.add_argument("--out", default="/dev/stdout")
    args = ap.parse_args()
    main(args.pattern, args.out)


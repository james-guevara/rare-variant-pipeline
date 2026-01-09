"""
Merge family genotype data with variant annotations (INNER join).
Only keeps variants that have matching annotations.
Uses sink_csv for streaming to disk.
"""
import polars as pl
import argparse
import sys


def main(fam_path: str, var_path: str, out_path: str):
    print(f"Scanning {fam_path}...", file=sys.stderr)
    fam = pl.scan_csv(fam_path, separator="\t", infer_schema_length=10000)

    print(f"Scanning {var_path}...", file=sys.stderr)
    var = pl.scan_csv(var_path, separator="\t", infer_schema_length=10000)

    join_cols = ["#CHROM", "POS0", "END", "REF", "ALT_specific"]

    # Verify columns exist
    fam_cols = set(fam.collect_schema().names())
    var_cols = set(var.collect_schema().names())
    
    missing_fam = set(join_cols) - fam_cols
    missing_var = set(join_cols) - var_cols
    
    if missing_fam:
        sys.exit(f"Family missing columns: {sorted(missing_fam)}")
    if missing_var:
        sys.exit(f"Variants missing columns: {sorted(missing_var)}")

    # Cast join columns to consistent types
    fam = fam.with_columns([
        pl.col("POS0").cast(pl.Int64, strict=False),
        pl.col("END").cast(pl.Int64, strict=False),
    ])

    var = var.with_columns([
        pl.col("POS0").cast(pl.Int64, strict=False),
        pl.col("END").cast(pl.Int64, strict=False),
    ])

    print("Joining (inner) and writing...", file=sys.stderr)
    merged = fam.join(var, on=join_cols, how="inner", suffix="_var")
    
    merged.sink_csv(out_path, separator="\t")
    
    print(f"Wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True, help="Resolved family genotypes TSV")
    ap.add_argument("--variants", required=True, help="Variant annotations TSV")
    ap.add_argument("--out", required=True, help="Output TSV")
    args = ap.parse_args()
    main(args.family, args.variants, args.out)

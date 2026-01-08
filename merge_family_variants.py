"""
Merge family genotype data with variant annotations.

Joins resolved family genotypes with annotated variant data.
Supports gzipped input/output.
"""
import polars as pl
import argparse
import sys
import gzip
import io

def main(fam_path, var_path, out_path):
    """
    Join family and variant tables on ALT_specific.
    Supports gzipped input/output (detected by .gz extension).
    """
    fam = pl.scan_csv(fam_path, separator="\t", infer_schema_length=0)
    var = pl.scan_csv(var_path, separator="\t", infer_schema_length=0)

    needed_fam = {"#CHROM", "POS0", "END", "REF", "ALT_specific"}
    needed_var = {"#CHROM", "POS0", "END", "REF", "ALT_specific"}

    if not needed_fam.issubset(set(fam.columns)):
        sys.exit(f"Family missing columns: {sorted(needed_fam - set(fam.columns))}")
    if not needed_var.issubset(set(var.columns)):
        sys.exit(f"Variants missing columns: {sorted(needed_var - set(var.columns))}")

    fam = fam.with_columns([
        pl.col("POS0").cast(pl.Int64, strict=False),
        pl.col("END").cast(pl.Int64, strict=False),
    ])

    var = var.with_columns([
        pl.col("POS0").cast(pl.Int64, strict=False),
        pl.col("END").cast(pl.Int64, strict=False),
    ])

    join_cols = ["#CHROM", "POS0", "END", "REF", "ALT_specific"]

    merged = fam.join(var, on=join_cols, how="left", suffix="_var")
    result = merged.collect()
    
    # Write output (gzipped if .gz extension)
    if out_path.endswith(".gz"):
        # Write to buffer first, then compress (polars write_csv to gzip handle corrupts data)
        buffer = io.BytesIO()
        result.write_csv(buffer, separator="\t")
        with gzip.open(out_path, "wb") as f:
            f.write(buffer.getvalue())
    else:
        result.write_csv(out_path, separator="\t")
    
    print(f"Wrote {result.height} rows -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    main(args.family, args.variants, args.out)

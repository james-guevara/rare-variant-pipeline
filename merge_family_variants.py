"""
Merge family genotype data with variant annotations.

Joins resolved family genotypes (from resolve_family_genotypes.py) with
annotated variant data (from reformat_variants.py) to create a unified
table with both genotype and annotation information.
"""
import polars as pl
import argparse
import sys

def main(fam_path, var_path, out_path):
    """
    Join family and variant tables.

    Args:
        fam_path: Resolved family genotypes TSV
        var_path: Annotated variants TSV
        out_path: Output merged TSV

    Join keys: #CHROM, POS0, END, REF, ALT (where ALT = ALT_resolved from family)
    """
    fam = pl.scan_csv(fam_path, separator="\t", infer_schema_length=0)
    var = pl.scan_csv(var_path, separator="\t", infer_schema_length=0)

    needed_fam = {"#CHROM","POS0","END","REF","ALT_resolved"}
    needed_var = {"#CHROM","POS0","END","REF","ALT"}

    if not needed_fam.issubset(set(fam.columns)):
        sys.exit(f"Family missing columns: {sorted(needed_fam - set(fam.columns))}")
    if not needed_var.issubset(set(var.columns)):
        sys.exit(f"Variants missing columns: {sorted(needed_var - set(var.columns))}")

    fam = fam.rename({"ALT": "ALT_raw"}).with_columns([
        pl.col("POS0").cast(pl.Int64, strict=False),
        pl.col("END").cast(pl.Int64, strict=False),
        pl.col("ALT_resolved").alias("ALT"),   # join key
    ])

    var = var.with_columns([
        pl.col("POS0").cast(pl.Int64, strict=False),
        pl.col("END").cast(pl.Int64, strict=False),
    ])

    join_cols = ["#CHROM","POS0","END","REF","ALT"]

    merged = fam.join(var, on=["#CHROM","POS0","END","REF","ALT"], how="left", suffix="_var")
    merged.collect(engine="streaming").write_csv(out_path, separator="\t")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    main(args.family, args.variants, args.out)


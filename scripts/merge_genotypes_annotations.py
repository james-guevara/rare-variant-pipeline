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
    fam = pl.scan_csv(fam_path, separator="\t", infer_schema_length=10000, null_values=["."])

    print(f"Scanning {var_path}...", file=sys.stderr)
    var = pl.scan_csv(var_path, separator="\t", infer_schema_length=10000, null_values=["."])

    join_cols = ["#CHROM", "POS0", "END", "REF", "ALT"]

    # Verify columns exist
    fam_cols = set(fam.collect_schema().names())
    var_cols = set(var.collect_schema().names())

    # Fall back to ALT_specific if ALT_specific exists in family data
    if "ALT_specific" in fam_cols and "ALT" not in fam_cols:
        fam = fam.rename({"ALT_specific": "ALT"})
    elif "ALT_specific" in fam_cols:
        # Both exist — join on ALT_specific from family side, ALT from variant side
        join_cols_fam = ["#CHROM", "POS0", "END", "REF", "ALT_specific"]
        join_cols_var = ["#CHROM", "POS0", "END", "REF", "ALT"]

        fam = fam.with_columns([
            pl.col("POS0").cast(pl.Int64, strict=False),
            pl.col("END").cast(pl.Int64, strict=False),
        ])
        var = var.with_columns([
            pl.col("POS0").cast(pl.Int64, strict=False),
            pl.col("END").cast(pl.Int64, strict=False),
        ])

        print("Joining (inner) on ALT_specific <-> ALT...", file=sys.stderr)
        merged = fam.join(
            var, left_on=join_cols_fam, right_on=join_cols_var,
            how="inner", suffix="_var"
        )
        merged.sink_csv(out_path, separator="\t")
        print(f"Wrote {out_path}", file=sys.stderr)
        return

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

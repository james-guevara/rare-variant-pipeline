"""
Resolve family genotypes for biallelic sites.

For each variant site in a family:
1. Determines carrier status (does GT contain '1'?)
2. Extracts allele-specific depths (AD_ref, AD_alt)
3. Flags missing genotypes

Input sites are already biallelic (split upstream by bcftools norm -m-).
"""
import polars as pl
import argparse
import sys


def process_families(df: pl.DataFrame, key_cols: list[str]) -> pl.DataFrame:
    """Process all families: add carrier flags and allele-specific depths."""

    # Carrier: does GT contain a '1'?
    df = df.with_columns(
        pl.col("GT").str.contains(r"\.").alias("is_missing_gt"),
    )

    df = df.with_columns(
        pl.when(pl.col("is_missing_gt"))
          .then(None)
          .when(pl.col("GT").str.contains(r"[1-9]"))
          .then(1)
          .otherwise(0)
          .alias("carrier"),
    )

    # ALT_specific = ALT (already biallelic)
    if "ALT" in df.columns:
        df = df.with_columns(
            pl.col("ALT").alias("ALT_specific"),
        )
    else:
        df = df.with_columns(pl.lit(None).alias("ALT_specific"))

    # Allele-specific depths: AD is "ref,alt" for biallelic
    if "AD" in df.columns:
        df = df.with_columns(
            pl.when(pl.col("AD").is_null() | (pl.col("AD") == "."))
              .then(None)
              .otherwise(pl.col("AD").str.split(","))
              .alias("_ad_list"),
        ).with_columns(
            pl.col("_ad_list").list.get(0).alias("AD_ref"),
            pl.col("_ad_list").list.get(1).alias("AD_alt"),
        ).drop("_ad_list")
    else:
        df = df.with_columns(
            pl.lit(None).alias("AD_ref"),
            pl.lit(None).alias("AD_alt"),
        )

    # Family-level carrier QC flag: does any carrier in the family exist?
    # (Used downstream for filtering families with no carriers)
    group_cols = key_cols + ["FAMILY"]
    existing_group_cols = [c for c in group_cols if c in df.columns]

    fam_has_carrier = (
        df.group_by(existing_group_cols)
          .agg(pl.col("carrier").max().alias("_fam_has_carrier"))
    )
    df = df.join(fam_has_carrier, on=existing_group_cols, how="left")

    # Only keep variant-family groups where at least one member carries ALT
    df = df.filter(pl.col("_fam_has_carrier") > 0).drop("_fam_has_carrier")

    return df


def main(tsv_path: str, output_path: str):
    print(f"Reading {tsv_path}...", file=sys.stderr)
    df = pl.read_csv(tsv_path, separator="\t", infer_schema_length=10000)
    print(f"  Loaded {df.height:,} rows", file=sys.stderr)

    required = {"#CHROM", "FAMILY", "SAMPLE", "GT"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    key_cols = [c for c in ["#CHROM", "POS0", "END", "POS", "REF", "ALT"] if c in df.columns]

    result = process_families(df, key_cols)

    print(f"  {df.height:,} -> {result.height:,} rows (kept families with carriers)", file=sys.stderr)
    result.write_csv(output_path, separator="\t")
    print(f"Wrote {output_path}", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="Input family genotype TSV")
    ap.add_argument("output", help="Output TSV")
    args = ap.parse_args()
    main(args.input, args.output)

"""
Resolve family genotypes by exploding on active ALT alleles.

For each variant site in a family:
1. Identifies all active ALT alleles across family members
2. Creates one row per family member per active allele
3. Adds carrier flags and compound het detection
4. Resolves allele-specific depths

Reads file once, processes per-family via group_by iteration.
"""
import polars as pl
import argparse
import sys


def process_family(df: pl.DataFrame, key_cols: list[str]) -> pl.DataFrame:
    """Process a single family's data."""
    group_cols = key_cols + ["FAMILY"]

    df = df.with_columns(
        pl.col("GT").str.extract_all(r"([1-9][0-9]*)").alias("alleles_in_gt"),
        pl.col("GT").str.contains(r"\.").alias("is_missing_gt"),
    )

    df = df.with_columns(
        (
            (pl.col("alleles_in_gt").list.len() > 1)
            & (pl.col("alleles_in_gt").list.n_unique() > 1)
        ).alias("sample_is_multiallelic_het")
    )

    # Family-level flags
    fam_flag = (
        df.group_by(group_cols)
          .agg(pl.col("sample_is_multiallelic_het").any().alias("family_has_multiallelic_het"))
    )

    # Active alleles per variant in this family
    fam_active = (
        df.group_by(group_cols)
          .agg(pl.col("alleles_in_gt").explode().unique().sort().alias("active_alleles"))
          .explode("active_alleles")
          .rename({"active_alleles": "allele_num_str"})
          .filter(pl.col("allele_num_str").is_not_null())
          .with_columns(pl.col("allele_num_str").cast(pl.Int32).alias("allele_num"))
          .drop("allele_num_str")
    )

    # If no active alleles (all ref or missing), return empty
    if fam_active.height == 0:
        return None

    # Join
    out = df.join(fam_active, on=group_cols, how="inner").join(fam_flag, on=group_cols, how="left")

    # Carrier flag
    out = out.with_columns(
        pl.when(pl.col("is_missing_gt"))
          .then(None)
          .when(pl.col("alleles_in_gt").list.contains(pl.col("allele_num").cast(pl.Utf8)))
          .then(1)
          .otherwise(0)
          .alias("carrier")
    )

    # ALT_specific
    if "ALT" in out.columns:
        out = out.with_columns(
            pl.col("ALT")
              .str.split(",")
              .list.get(pl.col("allele_num") - 1)
              .alias("ALT_specific")
        )
    else:
        out = out.with_columns(pl.lit(None).alias("ALT_specific"))

    # Allele-specific depths
    if "AD" in out.columns:
        out = out.with_columns(
            pl.when(pl.col("AD").is_null() | (pl.col("AD") == "."))
              .then(None)
              .otherwise(pl.col("AD").str.split(","))
              .alias("AD_list")
        ).with_columns(
            pl.col("AD_list").list.get(0).alias("AD_ref"),
            pl.col("AD_list").list.get(pl.col("allele_num")).alias("AD_alt"),
        ).drop("AD_list")
    else:
        out = out.with_columns(pl.lit(None).alias("AD_ref"), pl.lit(None).alias("AD_alt"))

    # Convert list to string
    out = out.with_columns(
        pl.col("alleles_in_gt").list.join(",").alias("alleles_in_gt")
    )

    return out


def resolve_family(tsv_path: str, output_path: str, progress_interval: int = 100):
    """
    Process family genotype TSV per-family.
    Reads file once, iterates over family groups.
    """
    print(f"Reading {tsv_path}...", file=sys.stderr)
    df = pl.read_csv(tsv_path, separator="\t", infer_schema_length=10000)
    print(f"  Loaded {df.height} rows", file=sys.stderr)

    columns = df.columns

    required = {"#CHROM", "FAMILY", "SAMPLE", "GT"}
    missing = sorted(required - set(columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Have: {columns}")

    key_cols = [c for c in ["#CHROM", "POS0", "END", "POS", "REF", "ALT"] if c in columns]
    if not key_cols:
        raise ValueError("Could not find any variant key columns.")

    # Get unique families
    families = df.get_column("FAMILY").unique().to_list()
    print(f"Processing {len(families)} families...", file=sys.stderr)

    # Output column order
    new_cols = [
        "allele_num", "ALT_specific", "carrier", "AD_ref", "AD_alt",
        "sample_is_multiallelic_het", "family_has_multiallelic_het",
        "alleles_in_gt", "is_missing_gt",
    ]
    select_cols = columns + [c for c in new_cols if c not in columns]

    total_rows = 0
    first_write = True

    # Iterate over family groups
    for i, (family_tuple, fam_df) in enumerate(df.group_by("FAMILY")):
        family = family_tuple[0] if isinstance(family_tuple, tuple) else family_tuple

        if fam_df.height == 0:
            continue

        # Process
        result = process_family(fam_df, key_cols)

        if result is None or result.height == 0:
            continue

        # Select columns in order
        result = result.select(select_cols)

        # Write (with header only for first chunk)
        if first_write:
            result.write_csv(output_path, separator="\t")
            first_write = False
        else:
            # Append without header
            with open(output_path, "a") as f:
                result.write_csv(f, separator="\t", include_header=False)

        total_rows += result.height

        # Progress update
        if (i + 1) % progress_interval == 0:
            print(f"  Processed {i + 1}/{len(families)} families, {total_rows} rows so far...", file=sys.stderr)

    print(f"Wrote {total_rows} rows -> {output_path}", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="Input family genotype TSV")
    ap.add_argument("output", help="Output TSV")
    ap.add_argument("--progress-interval", type=int, default=100, help="Print progress every N families (default: 100)")
    args = ap.parse_args()
    resolve_family(args.input, args.output, args.progress_interval)

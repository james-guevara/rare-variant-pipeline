"""
Resolve family genotypes by exploding on active ALT alleles.

For each variant site in a family:
1. Identifies all active ALT alleles across family members
2. Creates one row per family member per active allele
3. Adds carrier flags and compound het detection
4. Resolves allele-specific depths
"""
import polars as pl
import argparse
import sys


def resolve_family(tsv_path: str, output_path: str):
    """
    Process family genotype TSV to explode by active alleles.

    Args:
        tsv_path: Input family genotype TSV (from family_query.py)
        output_path: Output TSV path

    Output includes:
        - allele_num: The specific ALT allele number
        - ALT_resolved: Specific ALT sequence
        - carrier: 1 if this sample carries this allele, 0 otherwise
        - sample_is_multiallelic_het: This sample is 1/2, 2/3, etc.
        - family_has_multiallelic_het: Any family member is multiallelic het
    """
    try:
        df = pl.read_csv(tsv_path, separator="\t", infer_schema_length=0)

        required = {"#CHROM", "FAMILY", "SAMPLE", "GT"}
        missing = sorted(required - set(df.columns))
        if missing:
            raise ValueError(f"Missing required columns: {missing}. Have: {df.columns}")

        # choose a robust variant key (use what exists)
        key_cols = [c for c in ["#CHROM", "POS0", "END", "POS", "REF", "ALT"] if c in df.columns]
        if not key_cols:
            raise ValueError("Could not find any variant key columns among #CHROM/POS0/END/POS/REF/ALT.")

        group_cols = key_cols + ["FAMILY"]

        # parse non-ref allele numbers from GT (strings like "1", "2", ...)
        df = df.with_columns(
            pl.col("GT").str.extract_all(r"([1-9][0-9]*)").alias("alleles_in_gt"),
            pl.col("GT").str.contains(r"\.").alias("is_missing_gt"),
        )

        # sample-level multiallelic-het-at-site flag (e.g. GT=1/2)
        df = df.with_columns(
            (
                (pl.col("alleles_in_gt").list.len() > 1)
                & (pl.col("alleles_in_gt").list.n_unique() > 1)
            ).alias("sample_is_multiallelic_het")
        )

        # family-level flag: any sample at this site is multiallelic het
        fam_flag = (
            df.group_by(group_cols)
              .agg(pl.col("sample_is_multiallelic_het").any().alias("family_has_multiallelic_het"))
        )

        fam_active = (
            df.group_by(group_cols)
              .agg(pl.col("alleles_in_gt").explode().unique().sort().alias("active_alleles"))
              .explode("active_alleles")                      # <-- key fix
              .rename({"active_alleles": "allele_num_str"})
              .filter(pl.col("allele_num_str").is_not_null())
              .with_columns(pl.col("allele_num_str").cast(pl.Int32).alias("allele_num"))
              .drop("allele_num_str")
        )


        # join family allele blocks back to all family members (replicate rows)
        out = df.join(fam_active, on=group_cols, how="inner").join(fam_flag, on=group_cols, how="left")

        # carrier per allele block (null if missing GT, else 1/0)
        out = out.with_columns(
            pl.when(pl.col("is_missing_gt"))
              .then(None)
              .when(pl.col("alleles_in_gt").list.contains(pl.col("allele_num").cast(pl.Utf8)))
              .then(1)
              .otherwise(0)
              .alias("carrier")
        )

        # ALT_resolved from ALT (if present)
        if "ALT" in out.columns:
            out = out.with_columns(
                pl.col("ALT")
                  .str.split(",")
                  .list.get(pl.col("allele_num") - 1)
                  .alias("ALT_resolved")
            )
        else:
            out = out.with_columns(pl.lit(None).alias("ALT_resolved"))

        # AD_ref / AD_alt from AD (if present). AD is "ref,alt1,alt2,..."
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

        # keep original columns first, append new columns after
        new_cols = [
            "allele_num",
            "ALT_resolved",
            "carrier",
            "AD_ref",
            "AD_alt",
            "sample_is_multiallelic_het",
            "family_has_multiallelic_het",
            "alleles_in_gt",
            "is_missing_gt",
        ]
        base_cols = df.columns  # original order
        select_cols = base_cols + [c for c in new_cols if c not in base_cols]
        out = out.select(select_cols)

        out = out.with_columns(
            pl.col("alleles_in_gt").list.join(",").alias("alleles_in_gt")
        )

        out.write_csv(output_path, separator="\t")
        print(f"Wrote {out.height} rows -> {output_path}", file=sys.stderr)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Explode family genotypes by active ALT alleles; keep originals; add carrier/AD/compound-het flags.")
    ap.add_argument("input", help="Input family genotype TSV (can be .gz)")
    ap.add_argument("output", help="Output TSV (write to /dev/stdout if you want to bgzip)")
    args = ap.parse_args()
    resolve_family(args.input, args.output)


#!/usr/bin/env python3
"""
Reshape filtered variants to one row per child × variant.
Each row has child_GT, mother_GT, father_GT together.
Adds mendelian_error flag.

Usage:
    python reshape_trios.py filtered.tsv -o trios.tsv
    python reshape_trios.py filtered.tsv -o trios.tsv --ped custom.ped
"""

import argparse
import sys
import polars as pl

DEFAULT_PED = "/expanse/projects/sebat1/s3/data/sebat/nf_rare_spark_wes/resources/SPARK_iWES_v3.ped"


def load_pedigree(ped_path: str) -> pl.DataFrame:
    """Load PED file and identify roles (proband, father, mother)."""
    ped = pl.read_csv(
        ped_path,
        separator="\t",
        has_header=False,
        new_columns=["FID", "IID", "PID", "MID", "SEX", "PHENO"],
        schema={"FID": pl.Utf8, "IID": pl.Utf8, "PID": pl.Utf8, "MID": pl.Utf8, "SEX": pl.Int64, "PHENO": pl.Int64}
    )

    # Find who is a father/mother by looking at PID/MID references
    fathers = ped.filter(pl.col("PID") != "0").select("PID").unique().rename({"PID": "IID"})
    mothers = ped.filter(pl.col("MID") != "0").select("MID").unique().rename({"MID": "IID"})

    ped = ped.with_columns([
        pl.col("IID").is_in(fathers["IID"]).alias("is_father"),
        pl.col("IID").is_in(mothers["IID"]).alias("is_mother"),
    ])

    ped = ped.with_columns(
        pl.when((pl.col("PID") == "0") & (pl.col("MID") == "0") & pl.col("is_father"))
            .then(pl.lit("father"))
        .when((pl.col("PID") == "0") & (pl.col("MID") == "0") & pl.col("is_mother"))
            .then(pl.lit("mother"))
        .when(pl.col("PID") != "0")  # has parents = child
            .then(pl.lit("child"))
        .otherwise(pl.lit("other"))
        .alias("role")
    )

    return ped.select(["FID", "IID", "PID", "MID", "role"])


def reshape_trios(variants: pl.DataFrame, ped: pl.DataFrame) -> pl.DataFrame:
    """Reshape so each row is one child × variant with parent genotypes."""

    # Join variants with pedigree to get roles
    variants = variants.with_columns(pl.col("FAMILY").cast(pl.Utf8))
    variants_with_roles = variants.join(
        ped.select(["IID", "PID", "MID", "role"]),
        left_on="SAMPLE",
        right_on="IID",
        how="left"
    )

    # Variant key for joining
    variant_key = ["#CHROM", "POS", "REF", "ALT", "FAMILY"]

    # Get children rows
    children = variants_with_roles.filter(pl.col("role") == "child")

    # Get parent data
    parents = variants_with_roles.filter(pl.col("role").is_in(["father", "mother"]))

    mother_data = parents.filter(pl.col("role") == "mother").select(
        variant_key + ["SAMPLE", "GT", "GQ", "DP", "carrier"]
    ).rename({
        "SAMPLE": "mother_ID",
        "GT": "mother_GT",
        "GQ": "mother_GQ",
        "DP": "mother_DP",
        "carrier": "mother_carrier"
    })

    father_data = parents.filter(pl.col("role") == "father").select(
        variant_key + ["SAMPLE", "GT", "GQ", "DP", "carrier"]
    ).rename({
        "SAMPLE": "father_ID",
        "GT": "father_GT",
        "GQ": "father_GQ",
        "DP": "father_DP",
        "carrier": "father_carrier"
    })

    # Rename child columns
    result = children.rename({
        "SAMPLE": "child_ID",
        "GT": "child_GT",
        "GQ": "child_GQ",
        "DP": "child_DP",
        "carrier": "child_carrier"
    })

    # Join with parent data
    result = result.join(mother_data, on=variant_key, how="left")
    result = result.join(father_data, on=variant_key, how="left")

    # Add mendelian_error flag
    result = result.with_columns(
        pl.when(pl.col("mother_GT").is_null() | pl.col("father_GT").is_null())
            .then(pl.lit(None))  # incomplete trio
        # Child 1/1, parent 0/0 = error
        .when((pl.col("child_GT") == "1/1") & ((pl.col("mother_GT") == "0/0") | (pl.col("father_GT") == "0/0")))
            .then(pl.lit(1))
        # Child 0/1, both parents 0/0 = error (putative DNM or error)
        .when((pl.col("child_GT") == "0/1") & (pl.col("mother_GT") == "0/0") & (pl.col("father_GT") == "0/0"))
            .then(pl.lit(1))
        # Child 0/0, both parents have alt = error
        .when((pl.col("child_GT") == "0/0") & 
              pl.col("mother_GT").is_in(["0/1", "1/1"]) & 
              pl.col("father_GT").is_in(["0/1", "1/1"]))
            .then(pl.lit(1))
        .otherwise(pl.lit(0))
        .alias("mendelian_error")
    )

    # Add incomplete_trio flag
    result = result.with_columns(
        pl.when(pl.col("mother_GT").is_null() | pl.col("father_GT").is_null())
            .then(pl.lit(1))
            .otherwise(pl.lit(0))
            .alias("incomplete_trio")
    )

    # Select and order columns
    key_cols = ["#CHROM", "POS", "REF", "ALT", "SYMBOL", "Consequence", "FAMILY"]
    child_cols = ["child_ID", "child_GT", "child_GQ", "child_DP", "child_carrier"]
    mother_cols = ["mother_ID", "mother_GT", "mother_GQ", "mother_DP", "mother_carrier"]
    father_cols = ["father_ID", "father_GT", "father_GQ", "father_DP", "father_carrier"]
    flag_cols = ["incomplete_trio", "mendelian_error"]
    
    # Keep other annotation columns
    other_cols = [c for c in result.columns 
                  if c not in key_cols + child_cols + mother_cols + father_cols + flag_cols
                  and c not in ["IID", "PID", "MID", "role", "is_father", "is_mother"]]

    output_cols = key_cols + flag_cols + child_cols + mother_cols + father_cols + other_cols
    output_cols = [c for c in output_cols if c in result.columns]

    return result.select(output_cols).sort(["#CHROM", "POS", "FAMILY"])


def main():
    parser = argparse.ArgumentParser(description="Reshape variants to trio format")
    parser.add_argument("input", help="Input TSV file")
    parser.add_argument("--output", "-o", required=True, help="Output TSV file")
    parser.add_argument("--ped", default=DEFAULT_PED, help="PED file")

    args = parser.parse_args()

    print(f"Loading PED: {args.ped}", file=sys.stderr)
    ped = load_pedigree(args.ped)
    roles = ped["role"].value_counts()
    print(f"  Roles: {dict(roles.iter_rows())}", file=sys.stderr)

    print(f"Loading variants: {args.input}", file=sys.stderr)
    variants = pl.read_csv(args.input, separator="\t", infer_schema_length=10000)
    print(f"  {len(variants):,} rows", file=sys.stderr)

    print("Reshaping to trio format...", file=sys.stderr)
    result = reshape_trios(variants, ped)
    
    complete = result.filter(pl.col("incomplete_trio") == 0).height
    errors = result.filter(pl.col("mendelian_error") == 1).height
    print(f"  {len(result):,} child × variant rows", file=sys.stderr)
    print(f"  Complete trios: {complete:,}", file=sys.stderr)
    print(f"  Mendelian errors: {errors:,}", file=sys.stderr)

    result.write_csv(args.output, separator="\t")
    print(f"Saved to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

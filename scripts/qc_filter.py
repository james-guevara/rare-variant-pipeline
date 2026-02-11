"""
QC filter for merged rare variant output.

Adds a `qc_filter` column with semicolon-separated flag strings for variants
failing quality checks. Variants are NOT dropped — downstream analysis filters
with `WHERE qc_filter IS NULL`.

Two categories:
  - Variant-level: same for all rows at a site (segdup, repeat, AF, etc.)
  - Genotype-level: per-sample, but if ANY family member fails, ALL rows for
    that variant x family get flagged.

Also adds an `AB` (allele balance) column: AD_alt / (AD_ref + AD_alt).
"""
import polars as pl
import argparse
import sys


def compute_ab(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Add AB = AD_alt / (AD_ref + AD_alt)."""
    ad_ref = pl.col("AD_ref").cast(pl.Float64, strict=False)
    ad_alt = pl.col("AD_alt").cast(pl.Float64, strict=False)
    total = ad_ref + ad_alt
    return lf.with_columns(
        pl.when(total > 0)
        .then(ad_alt / total)
        .otherwise(None)
        .alias("AB")
    )


def build_variant_flags(
    lf: pl.LazyFrame,
    max_af: float,
    no_mane_filter: bool,
) -> pl.LazyFrame:
    """Add variant-level flag columns (same value for all rows at a site)."""
    cols = set(lf.collect_schema().names())
    flags = []

    # segDups > 0
    if "segDups" in cols:
        flags.append(
            pl.when(pl.col("segDups").cast(pl.Float64, strict=False) > 0)
            .then(pl.lit("segdup"))
            .otherwise(None)
            .alias("_vf_segdup")
        )

    # simpleRepeats > 0
    if "simpleRepeats" in cols:
        flags.append(
            pl.when(pl.col("simpleRepeats").cast(pl.Float64, strict=False) > 0)
            .then(pl.lit("simpleRepeat"))
            .otherwise(None)
            .alias("_vf_simpleRepeat")
        )

    # rmsk_Simple_repeat > 0
    if "rmsk_Simple_repeat" in cols:
        flags.append(
            pl.when(pl.col("rmsk_Simple_repeat").cast(pl.Float64, strict=False) > 0)
            .then(pl.lit("rmsk_Simple_repeat"))
            .otherwise(None)
            .alias("_vf_rmsk_Simple_repeat")
        )
    else:
        print("WARNING: rmsk_Simple_repeat column not found, skipping filter", file=sys.stderr)

    # rmsk_Low_complexity > 0
    if "rmsk_Low_complexity" in cols:
        flags.append(
            pl.when(pl.col("rmsk_Low_complexity").cast(pl.Float64, strict=False) > 0)
            .then(pl.lit("rmsk_Low_complexity"))
            .otherwise(None)
            .alias("_vf_rmsk_Low_complexity")
        )
    else:
        print("WARNING: rmsk_Low_complexity column not found, skipping filter", file=sys.stderr)

    # gnomAD AF > threshold (missing = pass)
    af_col = "gnomAD4.1_joint_AF"
    if af_col in cols:
        af_label = f"AF>{max_af}"
        flags.append(
            pl.when(pl.col(af_col).cast(pl.Float64, strict=False) > max_af)
            .then(pl.lit(af_label))
            .otherwise(None)
            .alias("_vf_af")
        )

    # non_MANE: MANE_SELECT == "."
    if not no_mane_filter and "MANE_SELECT" in cols:
        flags.append(
            pl.when(pl.col("MANE_SELECT") == ".")
            .then(pl.lit("non_MANE"))
            .otherwise(None)
            .alias("_vf_non_mane")
        )

    # LOFTEE_LC: LoF == "LC" (only when LoF is non-empty)
    if "LoF" in cols:
        flags.append(
            pl.when(pl.col("LoF") == "LC")
            .then(pl.lit("LOFTEE_LC"))
            .otherwise(None)
            .alias("_vf_loftee_lc")
        )

    return lf.with_columns(flags) if flags else lf


def build_genotype_flags(
    lf: pl.LazyFrame,
    min_gq: int,
    min_dp: int,
    het_ab_min: float,
    het_ab_max: float,
    hom_ab_min: float,
) -> pl.LazyFrame:
    """Add per-sample genotype flag columns."""
    cols = set(lf.collect_schema().names())
    flags = []

    # GQ < threshold (missing = fail)
    if "GQ" in cols:
        gq_label = f"GQ<{min_gq}"
        gq_val = pl.col("GQ").cast(pl.Float64, strict=False)
        flags.append(
            pl.when(gq_val.is_null() | (gq_val < min_gq))
            .then(pl.lit(gq_label))
            .otherwise(None)
            .alias("_gf_gq")
        )

    # DP < threshold (missing = fail)
    if "DP" in cols:
        dp_label = f"DP<{min_dp}"
        dp_val = pl.col("DP").cast(pl.Float64, strict=False)
        flags.append(
            pl.when(dp_val.is_null() | (dp_val < min_dp))
            .then(pl.lit(dp_label))
            .otherwise(None)
            .alias("_gf_dp")
        )

    # het AB outside [het_ab_min, het_ab_max]
    # het = GT matches 0/1 or 0|1 pattern (one allele is 0, other is non-0)
    if "GT" in cols and "AB" in cols:
        is_het = pl.col("GT").str.contains(r"^0[/|][1-9]$|^[1-9][/|]0$")
        ab = pl.col("AB")
        flags.append(
            pl.when(is_het & ((ab < het_ab_min) | (ab > het_ab_max)))
            .then(pl.lit("het_AB"))
            .otherwise(None)
            .alias("_gf_het_ab")
        )

        # hom-alt AB < hom_ab_min
        # hom-alt = GT matches N/N or N|N where N >= 1 and both alleles same
        is_hom_alt = pl.col("GT").str.contains(r"^([1-9])[/|]\1$")
        flags.append(
            pl.when(is_hom_alt & (ab < hom_ab_min))
            .then(pl.lit("hom_AB"))
            .otherwise(None)
            .alias("_gf_hom_ab")
        )

    return lf.with_columns(flags) if flags else lf


def propagate_genotype_flags(
    df: pl.DataFrame,
    gf_cols: list[str],
) -> pl.DataFrame:
    """Propagate genotype flags: if ANY family member fails, ALL rows for that
    variant x family get the flag."""
    if not gf_cols:
        return df

    group_cols = ["#CHROM", "POS0", "END", "ALT_specific", "FAMILY"]
    # Check which group cols actually exist
    existing_group_cols = [c for c in group_cols if c in df.columns]

    # For each genotype flag, compute family-level any()
    agg_exprs = [
        pl.col(c).is_not_null().any().alias(f"{c}_any")
        for c in gf_cols
    ]

    family_flags = df.group_by(existing_group_cols).agg(agg_exprs)

    # Join back
    df = df.join(family_flags, on=existing_group_cols, how="left")

    # Replace individual flags with family-propagated ones
    for c in gf_cols:
        any_col = f"{c}_any"
        # Original flag value (non-null means it has the label string)
        original_label = df[c].drop_nulls()[0] if df[c].drop_nulls().len() > 0 else c
        df = df.with_columns(
            pl.when(pl.col(any_col))
            .then(pl.lit(original_label))
            .otherwise(None)
            .alias(c)
        ).drop(any_col)

    return df


def combine_flags(df: pl.DataFrame, flag_cols: list[str]) -> pl.DataFrame:
    """Combine all flag columns into a single semicolon-separated qc_filter column."""
    if not flag_cols:
        return df.with_columns(pl.lit(None).cast(pl.Utf8).alias("qc_filter"))

    combined = pl.concat_list([pl.col(c) for c in flag_cols]).list.drop_nulls()
    return (
        df.with_columns(
            pl.when(combined.list.len() > 0)
            .then(combined.list.join(";"))
            .otherwise(None)
            .alias("qc_filter")
        )
        .drop(flag_cols)
    )


def main():
    parser = argparse.ArgumentParser(description="Add QC filter column to merged variant TSV")
    parser.add_argument("input", help="Input merged TSV")
    parser.add_argument("output", help="Output QC-annotated TSV")
    parser.add_argument("--min-gq", type=int, default=20, help="Min GQ (default: 20)")
    parser.add_argument("--min-dp", type=int, default=10, help="Min DP (default: 10)")
    parser.add_argument("--max-af", type=float, default=0.001, help="Max gnomAD AF (default: 0.001)")
    parser.add_argument("--het-ab-min", type=float, default=0.25, help="Min AB for hets (default: 0.25)")
    parser.add_argument("--het-ab-max", type=float, default=0.75, help="Max AB for hets (default: 0.75)")
    parser.add_argument("--hom-ab-min", type=float, default=0.9, help="Min AB for hom-alt (default: 0.9)")
    parser.add_argument("--no-mane-filter", action="store_true", help="Disable non-MANE filter")
    args = parser.parse_args()

    print(f"Reading {args.input}...", file=sys.stderr)
    lf = pl.scan_csv(args.input, separator="\t", infer_schema_length=10000)

    # Compute AB column
    lf = compute_ab(lf)

    # Add variant-level flags
    lf = build_variant_flags(lf, args.max_af, args.no_mane_filter)

    # Add genotype-level flags
    lf = build_genotype_flags(
        lf, args.min_gq, args.min_dp,
        args.het_ab_min, args.het_ab_max, args.hom_ab_min,
    )

    # Collect for group_by + join (per-chrom fits in memory)
    print("Collecting and propagating genotype flags...", file=sys.stderr)
    df = lf.collect()

    # Identify flag columns
    all_cols = df.columns
    vf_cols = [c for c in all_cols if c.startswith("_vf_")]
    gf_cols = [c for c in all_cols if c.startswith("_gf_")]

    # Propagate genotype flags across families
    df = propagate_genotype_flags(df, gf_cols)

    # Combine all flags into qc_filter
    flag_cols = vf_cols + gf_cols
    df = combine_flags(df, flag_cols)

    print(f"Writing {args.output}...", file=sys.stderr)
    df.write_csv(args.output, separator="\t")

    # Summary
    total = len(df)
    passing = df.filter(pl.col("qc_filter").is_null()).height
    print(f"QC filter: {passing}/{total} rows pass ({100*passing/total:.1f}%)", file=sys.stderr)

    # Per-flag counts
    for flag in ["segdup", "simpleRepeat", "rmsk_Simple_repeat", "rmsk_Low_complexity",
                 f"AF>{args.max_af}", "non_MANE", "LOFTEE_LC",
                 f"GQ<{args.min_gq}", f"DP<{args.min_dp}", "het_AB", "hom_AB"]:
        n = df.filter(
            pl.col("qc_filter").str.split(";").list.contains(flag)
        ).height
        if n > 0:
            print(f"  {flag}: {n} rows", file=sys.stderr)


if __name__ == "__main__":
    main()

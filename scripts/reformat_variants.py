"""
Reformat multi-allelic VCF annotations and add genomic constraint/annotation layers.
Fully lazy/streaming with polars where possible.

Supports two modes:
  --mode coding     (default) Filter to HIGH/MODERATE IMPACT variants
  --mode regulatory  Intersect with regulatory BED tracks (ABC, PsychENCODE, etc.)
"""
import polars as pl
import polars_bio as pb
import argparse
import sys
from pathlib import Path

CONSEQUENTIAL_IMPACTS = {"HIGH", "MODERATE"}

# Expected regulatory BED filenames in --regulatory-beds directory
REGULATORY_TRACKS = {
    "abc": "abc_enhancers.bed",           # ABC enhancer-gene predictions (4+ cols: chrom, start, end, target_gene, score)
    "psychencode": "psychencode.bed",     # PsychENCODE brain enhancers
    "chromhmm": "chromhmm_fetal_brain.bed",  # Roadmap ChromHMM fetal brain active states
    "phastcons": "phastConsElements.bed", # phastCons conserved elements
    "ccre": "encodeCcreCombined.bed",     # ENCODE cCREs
}

# Problematic region BED filenames (relative to resources_dir/repeats/)
PROBLEMATIC_TRACKS = [
    "genomicSuperDups.bed",
    "simpleRepeat.bed",
    "encodeBlacklist.bed",
]


def add_gnomad_constraint_by_gene(
    lf: pl.LazyFrame,
    constraint_tsv: str,
    gene_col: str = "SYMBOL",
) -> pl.LazyFrame:
    cons_gene = (
        pl.scan_csv(constraint_tsv, separator="\t", null_values=["NA"],
            schema_overrides={"lof.pLI": pl.Float64, "lof.oe_ci.upper": pl.Float64})
        .group_by("gene")
        .agg([
            pl.col("lof.pLI").max().alias("gnomad_lof_pLI_max"),
            pl.col("lof.oe_ci.upper").min().alias("gnomad_lof_oe_ci_upper_min"),
            pl.when(pl.col("canonical") == True).then(pl.col("lof.pLI")).otherwise(None)
              .max().alias("gnomad_lof_pLI_canonical"),
            pl.when(pl.col("canonical") == True).then(pl.col("lof.oe_ci.upper")).otherwise(None)
              .min().alias("gnomad_lof_oe_ci_upper_canonical"),
        ])
        .rename({"gene": gene_col})
    )
    return lf.join(cons_gene, on=gene_col, how="left")


def add_gnomad_constraint_by_transcript(
    lf: pl.LazyFrame,
    constraint_tsv: str,
    feature_col: str = "Feature",
) -> pl.LazyFrame:
    cons_tx = (
        pl.scan_csv(constraint_tsv, separator="\t", null_values=["NA"],
            schema_overrides={"lof.pLI": pl.Float64, "lof.oe_ci.upper": pl.Float64})
        .select([
            pl.col("transcript").alias(feature_col),
            pl.col("lof.pLI").alias("gnomad_lof_pLI_tx"),
            pl.col("lof.oe_ci.upper").alias("gnomad_lof_oe_ci_upper_tx"),
        ])
    )
    return lf.join(cons_tx, on=feature_col, how="left")


def add_genebayes(
    lf: pl.LazyFrame,
    genebayes_tsv: str,
    gene_col: str = "Gene",
) -> pl.LazyFrame:
    gb = pl.scan_csv(genebayes_tsv, separator="\t",
        schema_overrides={"obs_lof": pl.Int64, "exp_lof": pl.Float64, "prior_mean": pl.Float64,
            "post_mean": pl.Float64, "post_lower_95": pl.Float64, "post_upper_95": pl.Float64})
    return lf.join(
        gb.select([
            "ensg",
            pl.col("obs_lof").alias("genebayes_obs_lof"),
            pl.col("exp_lof").alias("genebayes_exp_lof"),
            pl.col("prior_mean").alias("genebayes_prior_mean"),
            pl.col("post_mean").alias("genebayes_post_mean"),
            pl.col("post_lower_95").alias("genebayes_post_lower_95"),
            pl.col("post_upper_95").alias("genebayes_post_upper_95"),
        ]),
        left_on=gene_col,
        right_on="ensg",
        how="left",
    )


def add_1kb_constraint(
    lf: pl.LazyFrame,
    cons_gz: str,
    chrom_col: str = "#CHROM",
    start_col: str = "POS0",
    end_col: str = "END",
) -> pl.LazyFrame:
    v_iv = lf.select([
        pl.col(chrom_col).alias("chrom"),
        pl.col(start_col).cast(pl.Int64).alias("start"),
        pl.col(end_col).cast(pl.Int64).alias("end"),
    ])
    cons = pl.scan_csv(cons_gz, separator="\t",
        schema_overrides={"start": pl.Int64, "end": pl.Int64, "z": pl.Float64, "oe": pl.Float64}
    ).select(["chrom", "start", "end", "z", "oe"])

    pairs = pb.overlap(v_iv, cons, use_zero_based=True, output_type="polars.LazyFrame")
    summ = (
        pairs.group_by(["chrom_1", "start_1", "end_1"])
        .agg([
            pl.col("z_2").mean().alias("constraint_z_1kb_mean"),
            pl.col("oe_2").mean().alias("constraint_oe_1kb_mean"),
        ])
        .rename({"chrom_1": chrom_col, "start_1": start_col, "end_1": end_col})
        .with_columns([pl.col(start_col).cast(pl.Utf8), pl.col(end_col).cast(pl.Utf8)])
    )
    return lf.join(summ, on=[chrom_col, start_col, end_col], how="left")


def add_bed_overlap_flag(
    lf: pl.LazyFrame,
    bed_path: str,
    prefix: str,
    chrom_col: str = "#CHROM",
    start_col: str = "POS0",
    end_col: str = "END",
) -> pl.LazyFrame:
    """Add overlap count with BED regions. Only collects unique coordinates (small)."""
    # Only collect unique coordinates - much smaller than full data
    coords = lf.select([chrom_col, start_col, end_col]).unique().collect()

    bed = pl.scan_csv(bed_path, separator="\t", has_header=False).select([
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").cast(pl.Int64).alias("start"),
        pl.col("column_3").cast(pl.Int64).alias("end"),
    ])

    v_iv = coords.select([
        pl.col(chrom_col).alias("chrom"),
        pl.col(start_col).cast(pl.Int64).alias("start"),
        pl.col(end_col).cast(pl.Int64).alias("end"),
    ])

    counts = pb.count_overlaps(v_iv, bed, use_zero_based=True, output_type="polars.DataFrame")

    # Add overlap count to coordinates
    coord_counts = coords.with_columns(counts.select(pl.col("count").alias(prefix)))

    # Join back to lazy frame
    return lf.join(coord_counts.lazy(), on=[chrom_col, start_col, end_col], how="left")


def filter_regulatory(
    lf: pl.LazyFrame,
    beds_dir: str,
    chrom_col: str = "#CHROM",
    start_col: str = "POS0",
    end_col: str = "END",
) -> pl.LazyFrame:
    """Keep only variants overlapping regulatory BED regions.
    Also adds columns identifying which tracks were hit and ABC target genes."""
    beds_path = Path(beds_dir)

    # Collect unique coordinates once for all intersections
    coords = lf.select([chrom_col, start_col, end_col]).unique().collect()
    v_iv = coords.select([
        pl.col(chrom_col).alias("chrom"),
        pl.col(start_col).cast(pl.Int64).alias("start"),
        pl.col(end_col).cast(pl.Int64).alias("end"),
    ])

    # Start with all-zero flag columns
    flag_df = coords.clone()

    for track_name, bed_filename in REGULATORY_TRACKS.items():
        bed_file = beds_path / bed_filename
        col_name = f"in_{track_name}"

        if not bed_file.exists():
            print(f"WARNING: Regulatory BED not found: {bed_file}, skipping {track_name}", file=sys.stderr)
            flag_df = flag_df.with_columns(pl.lit(0).alias(col_name))
            continue

        if track_name == "abc":
            # ABC has extra columns: chrom, start, end, target_gene, score
            bed = pl.scan_csv(str(bed_file), separator="\t", has_header=False)
            bed_cols = bed.collect_schema().names()

            bed_intervals = bed.select([
                pl.col("column_1").alias("chrom"),
                pl.col("column_2").cast(pl.Int64).alias("start"),
                pl.col("column_3").cast(pl.Int64).alias("end"),
            ])

            # Overlap to get flags
            pairs = pb.overlap(
                v_iv, bed_intervals,
                use_zero_based=True, output_type="polars.DataFrame",
            )

            if len(pairs) > 0:
                # Also extract target gene and score from ABC BED
                bed_with_meta = bed.select([
                    pl.col("column_1").alias("chrom"),
                    pl.col("column_2").cast(pl.Int64).alias("start"),
                    pl.col("column_3").cast(pl.Int64).alias("end"),
                    pl.col("column_4").alias("abc_target_gene") if len(bed_cols) > 3 else pl.lit(None).alias("abc_target_gene"),
                    pl.col("column_5").cast(pl.Float64, strict=False).alias("abc_score") if len(bed_cols) > 4 else pl.lit(None).cast(pl.Float64).alias("abc_score"),
                ]).collect()

                # Get overlap pairs with ABC metadata
                abc_pairs = pb.overlap(
                    v_iv, bed_with_meta.lazy(),
                    use_zero_based=True, output_type="polars.DataFrame",
                )

                # Aggregate: pick best ABC score and its target gene per variant
                abc_summ = (
                    abc_pairs.group_by(["chrom_1", "start_1", "end_1"])
                    .agg([
                        pl.col("abc_score_2").max().alias("abc_score"),
                        # target gene from the row with max ABC score
                        pl.col("abc_target_gene_2").sort_by("abc_score_2", descending=True).first().alias("abc_target_gene"),
                    ])
                    .rename({"chrom_1": chrom_col, "start_1": start_col, "end_1": end_col})
                    .with_columns([
                        pl.col(start_col).cast(pl.Utf8),
                        pl.col(end_col).cast(pl.Utf8),
                    ])
                )
                flag_df = flag_df.join(abc_summ, on=[chrom_col, start_col, end_col], how="left")
            else:
                flag_df = flag_df.with_columns([
                    pl.lit(None).cast(pl.Float64).alias("abc_score"),
                    pl.lit(None).cast(pl.Utf8).alias("abc_target_gene"),
                ])

            # Binary flag: did this variant overlap any ABC region?
            counts = pb.count_overlaps(v_iv, bed_intervals, use_zero_based=True, output_type="polars.DataFrame")
            flag_df = flag_df.with_columns(
                (counts.select(pl.col("count")) > 0).cast(pl.Int32).to_series().alias(col_name)
            )
        else:
            # Standard BED: just chrom, start, end
            bed = pl.scan_csv(str(bed_file), separator="\t", has_header=False).select([
                pl.col("column_1").alias("chrom"),
                pl.col("column_2").cast(pl.Int64).alias("start"),
                pl.col("column_3").cast(pl.Int64).alias("end"),
            ])

            counts = pb.count_overlaps(v_iv, bed, use_zero_based=True, output_type="polars.DataFrame")
            flag_df = flag_df.with_columns(
                (counts.select(pl.col("count")) > 0).cast(pl.Int32).to_series().alias(col_name)
            )

    # Build "in any regulatory track" filter
    track_flag_cols = [f"in_{t}" for t in REGULATORY_TRACKS]
    existing_flags = [c for c in track_flag_cols if c in flag_df.columns]
    any_hit = pl.lit(False)
    for col in existing_flags:
        any_hit = any_hit | (pl.col(col) > 0)
    flag_df = flag_df.with_columns(any_hit.alias("in_any_regulatory"))

    # Join flags back to main LazyFrame, then filter
    lf_with_flags = lf.join(flag_df.lazy(), on=[chrom_col, start_col, end_col], how="left")
    return lf_with_flags.filter(pl.col("in_any_regulatory"))


def filter_problematic_regions(
    lf: pl.LazyFrame,
    resources_dir: str,
    chrom_col: str = "#CHROM",
    start_col: str = "POS0",
    end_col: str = "END",
) -> pl.LazyFrame:
    """Exclude variants in segdups, simple repeats, and ENCODE blacklist regions."""
    res = Path(resources_dir)

    for bed_filename in PROBLEMATIC_TRACKS:
        bed_file = res / "repeats" / bed_filename
        flag_name = f"_exclude_{Path(bed_filename).stem}"

        if not bed_file.exists():
            print(f"WARNING: Problematic region BED not found: {bed_file}, skipping", file=sys.stderr)
            continue

        lf = add_bed_overlap_flag(lf, str(bed_file), flag_name, chrom_col, start_col, end_col)
        lf = lf.filter(pl.col(flag_name) == 0).drop(flag_name)

    return lf


def clean_column_names(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Clean up column names from bcftools +split-vep output.

    - Rename '(null)' to 'INFO'
    - Strip 'CSQ' prefix from VEP columns (added by -p CSQ flag)
    """
    cols = lf.collect_schema().names()
    rename_map = {}
    for col in cols:
        if col == "(null)":
            rename_map[col] = "INFO"
        elif col.startswith("CSQ"):
            rename_map[col] = col[3:]  # Strip 'CSQ' prefix
    return lf.rename(rename_map) if rename_map else lf


def strip_csq_from_info(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Remove the CSQ=... blob from INFO column since VEP fields are now separate columns."""
    return lf.with_columns(
        pl.col("INFO").str.replace(r";CSQ.*", "").alias("INFO")
    )


def reformat_variants_lazy(tsv_path: str) -> pl.LazyFrame:
    lf = pl.scan_csv(tsv_path, separator="\t", infer_schema_length=0)

    # Clean column names first (handle (null) -> INFO, strip CSQ prefix)
    lf = clean_column_names(lf)

    # Strip CSQ blob from INFO column
    lf = strip_csq_from_info(lf)

    idx = pl.col("ALLELE_NUM").cast(pl.Int32, strict=False).fill_null(1).sub(1)
    per_allele_fields = ["AC", "AF", "AQ", "MLEAC", "MLEAF"]

    lf = lf.with_columns(
        pl.col("ALT").str.split(",").list.get(idx, null_on_oob=True).alias("ALT_specific"),
        *[pl.col("INFO").str.extract(f"{field}=([^;]+)", 1).str.split(",").list.get(idx, null_on_oob=True).fill_null("NA").alias(field)
          for field in per_allele_fields],
    )

    # Reorder columns
    cols = lf.collect_schema().names()
    if "ALT" in cols and "ALT_specific" in cols:
        cols.remove("ALT_specific")
        cols.insert(cols.index("ALT") + 1, "ALT_specific")
    if "INFO" in cols:
        insert_idx = cols.index("INFO") + 1
        for field in per_allele_fields:
            if field in cols:
                cols.remove(field)
                cols.insert(insert_idx, field)
                insert_idx += 1
    return lf.select(cols)


def create_merged_bed(lf: pl.LazyFrame, output_path: str) -> None:
    merged = pb.merge(
        lf.select([
            pl.col("#CHROM").alias("chrom"),
            pl.col("POS0").cast(pl.Int64).alias("start"),
            pl.col("END").cast(pl.Int64).alias("end"),
        ]),
        use_zero_based=True, min_dist=1, output_type="polars.LazyFrame"
    ).rename({"chrom": "#CHROM", "start": "POS0", "end": "END"})
    merged.sink_csv(output_path, separator="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--bed", required=True)
    parser.add_argument("--consequential", required=True)
    parser.add_argument("--consequential-bed", required=True)
    parser.add_argument("--resources-dir", required=True)
    parser.add_argument("--mode", choices=["coding", "regulatory"], default="coding")
    parser.add_argument("--regulatory-beds", help="Directory containing regulatory BED files")
    parser.add_argument("--filter-repeats", action="store_true",
                        help="Exclude variants in problematic regions (segdups, repeats, blacklist)")
    args = parser.parse_args()

    res = Path(args.resources_dir)
    constraint = res / "gnomAD/constraint/gnomad.v4.1.constraint_metrics.tsv"
    segdups = res / "repeats/genomicSuperDups.bed"
    simple_repeats = res / "repeats/simpleRepeat.bed"
    constraint_1kb = res / "gnomAD/Genomic_constraint/constraint_z_genome_1kb.qc.download.txt.gz"
    genebayes = res / "GeneBayes/output/Supplementary_Table_1.tsv"

    print("Processing variants (lazy)...", file=sys.stderr)
    lf = (
        reformat_variants_lazy(args.input)
        .pipe(add_bed_overlap_flag, str(segdups), "segDups")
        .pipe(add_bed_overlap_flag, str(simple_repeats), "simpleRepeats")
        .pipe(add_1kb_constraint, str(constraint_1kb))
        .pipe(add_genebayes, str(genebayes))
        .pipe(add_gnomad_constraint_by_transcript, str(constraint), feature_col="Feature")
        .pipe(add_gnomad_constraint_by_gene, str(constraint), gene_col="SYMBOL")
    )

    print(f"Writing {args.output}...", file=sys.stderr)
    lf.sink_csv(args.output, separator="\t")

    print(f"Writing {args.bed}...", file=sys.stderr)
    create_merged_bed(lf, args.bed)

    # Branch: coding vs regulatory filtering
    if args.mode == "coding":
        if "IMPACT" not in lf.collect_schema().names():
            print("ERROR: IMPACT column not found", file=sys.stderr)
            sys.exit(1)
        lf_filtered = lf.filter(pl.col("IMPACT").is_in(list(CONSEQUENTIAL_IMPACTS)))
    elif args.mode == "regulatory":
        if not args.regulatory_beds:
            print("ERROR: --regulatory-beds is required when --mode=regulatory", file=sys.stderr)
            sys.exit(1)
        lf_filtered = filter_regulatory(lf, args.regulatory_beds)

    if args.filter_repeats:
        lf_filtered = filter_problematic_regions(lf_filtered, args.resources_dir)

    print(f"Writing {args.consequential} ({args.mode} mode)...", file=sys.stderr)
    lf_filtered.sink_csv(args.consequential, separator="\t")
    print(f"Writing {args.consequential_bed}...", file=sys.stderr)
    create_merged_bed(lf_filtered, args.consequential_bed)

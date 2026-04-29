"""
Prepare variants from split-VEP TSV output.

1. Clean column names (strip CSQ prefix, rename (null) -> INFO)
2. Extract per-allele INFO fields (AC, AF, AQ, MLEAC, MLEAF)
   (Sites are already biallelic — multi-allelics split upstream by bcftools norm -m-)
3. Filter to rare variants (gnomAD4.1_joint_AF < threshold, default 0.01)
4. Output: all rare variants TSV + BED, consequential TSV + BED

Supports three modes for consequential filtering:
  --mode coding     (default) HIGH/MODERATE IMPACT variants
  --mode regulatory  Intersect with regulatory BED tracks
  --mode splicing    Deep intronic variants with SpliceAI delta score >= threshold

Region overlaps, constraint scores, and other annotations are handled downstream
in the post-processing pipeline.
"""
import polars as pl
import polars_bio as pb
import argparse
import sys
from pathlib import Path

CONSEQUENTIAL_IMPACTS = {"HIGH", "MODERATE"}

# SpliceAI delta score columns from VEP plugin
SPLICEAI_COLS = [
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
]

# Expected regulatory BED filenames
REGULATORY_TRACKS = {
    "abc": "abc_enhancers.bed",
    "psychencode": "psychencode.bed",
    "chromhmm": "chromhmm_fetal_brain.bed",
    "phastcons": "phastConsElements.bed",
    "ccre": "encodeCcreCombined.bed",
    "cpg": "cpgIslandExt.bed",
}



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
            rename_map[col] = col[3:]
    return lf.rename(rename_map) if rename_map else lf


def strip_csq_from_info(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Remove the CSQ=... blob from INFO column since VEP fields are now separate columns."""
    return lf.with_columns(
        pl.col("INFO").str.replace(r";CSQ.*", "").alias("INFO")
    )


def extract_info_fields(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Extract per-allele INFO fields. Sites are already biallelic (split upstream)."""
    info_fields = ["AC", "AF", "AQ", "MLEAC", "MLEAF"]

    lf = lf.with_columns(
        *[pl.col("INFO").str.extract(f"{field}=([^;]+)", 1)
          .fill_null("NA").alias(field)
          for field in info_fields],
    )

    # Reorder: info fields after INFO column
    cols = lf.collect_schema().names()
    if "INFO" in cols:
        insert_idx = cols.index("INFO") + 1
        for field in info_fields:
            if field in cols:
                cols.remove(field)
                cols.insert(insert_idx, field)
                insert_idx += 1
    return lf.select(cols)



def filter_regulatory(
    lf: pl.LazyFrame,
    beds_dir: str,
    chrom_col: str = "#CHROM",
    start_col: str = "POS0",
    end_col: str = "END",
) -> pl.LazyFrame:
    """Keep only variants overlapping regulatory BED regions."""
    beds_path = Path(beds_dir)

    coords = lf.select([chrom_col, start_col, end_col]).unique().collect()
    v_iv = coords.select([
        pl.col(chrom_col).alias("chrom"),
        pl.col(start_col).cast(pl.Int64).alias("start"),
        pl.col(end_col).cast(pl.Int64).alias("end"),
    ])

    flag_df = coords.clone()
    for track_name, bed_filename in REGULATORY_TRACKS.items():
        bed_file = beds_path / bed_filename
        col_name = f"in_{track_name}"

        if not bed_file.exists():
            print(f"WARNING: Regulatory BED not found: {bed_file}, skipping {track_name}", file=sys.stderr)
            flag_df = flag_df.with_columns(pl.lit(0).alias(col_name))
            continue

        bed = pl.scan_csv(str(bed_file), separator="\t", has_header=False).select([
            pl.col("column_1").alias("chrom"),
            pl.col("column_2").cast(pl.Int64).alias("start"),
            pl.col("column_3").cast(pl.Int64).alias("end"),
        ])
        counts = pb.count_overlaps(v_iv, bed, use_zero_based=True, output_type="polars.DataFrame")
        flag_df = flag_df.with_columns(
            (counts.select(pl.col("count")) > 0).cast(pl.Int32).to_series().alias(col_name)
        )

    track_flag_cols = [f"in_{t}" for t in REGULATORY_TRACKS]
    existing_flags = [c for c in track_flag_cols if c in flag_df.columns]
    any_hit = pl.lit(False)
    for col in existing_flags:
        any_hit = any_hit | (pl.col(col) > 0)
    flag_df = flag_df.with_columns(any_hit.alias("in_any_regulatory"))

    lf_with_flags = lf.join(flag_df.lazy(), on=[chrom_col, start_col, end_col], how="left")
    return lf_with_flags.filter(pl.col("in_any_regulatory"))


def filter_splicing(lf: pl.LazyFrame, threshold: float = 0.2) -> pl.LazyFrame:
    """Keep variants with max SpliceAI delta score >= threshold,
    excluding those already captured by coding mode (HIGH/MODERATE IMPACT)."""
    spliceai_max = pl.max_horizontal(
        *[pl.col(c).cast(pl.Float64, strict=False) for c in SPLICEAI_COLS]
    )
    return lf.filter(
        (spliceai_max >= threshold)
        & ~pl.col("IMPACT").is_in(list(CONSEQUENTIAL_IMPACTS))
    ).with_columns(spliceai_max.alias("SpliceAI_max_delta"))


def create_merged_bed(lf: pl.LazyFrame, output_path: str) -> None:
    """Create a merged BED file from variant coordinates."""
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
    parser = argparse.ArgumentParser(
        description="Prepare variants: clean, resolve multi-allelics, filter rare, split by impact."
    )
    parser.add_argument("input", help="Input TSV from bcftools +split-vep")
    parser.add_argument("output", help="Output TSV (all rare variants)")
    parser.add_argument("--bed", required=True, help="Output BED (all rare variants)")
    parser.add_argument("--consequential", required=True, help="Output TSV (HIGH/MODERATE impact only)")
    parser.add_argument("--consequential-bed", required=True, help="Output BED (consequential only)")
    parser.add_argument("--mode", choices=["coding", "regulatory", "splicing"], default="coding",
                        help="Consequential filtering mode (default: coding)")
    parser.add_argument("--regulatory-beds", help="Directory containing regulatory BED files (for regulatory mode)")
    parser.add_argument("--spliceai-threshold", type=float, default=0.2,
                        help="Min SpliceAI max delta score for splicing mode (default: 0.2)")
    parser.add_argument("--max-cohort-af", type=float, default=1.0,
                        help="Drop variants with cohort AF >= this value (default: 1.0 = no filter). "
                             "Set to e.g. 0.01 to restrict to rare variants before family processing.")
    args = parser.parse_args()

    print("Preparing variants...", file=sys.stderr)

    # 1. Read and clean
    lf = pl.scan_csv(args.input, separator="\t", infer_schema_length=0)
    lf = clean_column_names(lf)
    lf = strip_csq_from_info(lf)

    # 2. Extract per-allele INFO fields (sites already biallelic from bcftools norm -m-)
    lf = extract_info_fields(lf)

    # 2a. Optional cohort-AF filter (rows with missing/non-numeric AF are dropped)
    if args.max_cohort_af < 1.0:
        lf = lf.filter(
            pl.col("AF").cast(pl.Float64, strict=False) < args.max_cohort_af
        )

    # 3. Write all variants (cohort-AF-filtered if --max-cohort-af set)
    n_total = lf.select(pl.len()).collect().item()
    print(f"  Total variants: {n_total:,}  (max-cohort-af={args.max_cohort_af})", file=sys.stderr)
    print(f"Writing {args.output}...", file=sys.stderr)
    lf.sink_csv(args.output, separator="\t")

    print(f"Writing {args.bed}...", file=sys.stderr)
    create_merged_bed(lf, args.bed)

    # 5. Filter to consequential variants based on mode
    if args.mode == "coding":
        if "IMPACT" not in lf.collect_schema().names():
            print("ERROR: IMPACT column not found", file=sys.stderr)
            sys.exit(1)
        lf_conseq = lf.filter(pl.col("IMPACT").is_in(list(CONSEQUENTIAL_IMPACTS)))
    elif args.mode == "regulatory":
        if not args.regulatory_beds:
            print("ERROR: --regulatory-beds is required for regulatory mode", file=sys.stderr)
            sys.exit(1)
        lf_conseq = filter_regulatory(lf, args.regulatory_beds)
    elif args.mode == "splicing":
        schema_names = lf.collect_schema().names()
        missing = [c for c in SPLICEAI_COLS if c not in schema_names]
        if missing:
            print(f"ERROR: SpliceAI columns not found: {missing}", file=sys.stderr)
            sys.exit(1)
        lf_conseq = filter_splicing(lf, args.spliceai_threshold)

    n_conseq = lf_conseq.select(pl.len()).collect().item()
    print(f"  Consequential ({args.mode} mode): {n_conseq:,} variants", file=sys.stderr)

    print(f"Writing {args.consequential}...", file=sys.stderr)
    lf_conseq.sink_csv(args.consequential, separator="\t")

    print(f"Writing {args.consequential_bed}...", file=sys.stderr)
    create_merged_bed(lf_conseq, args.consequential_bed)

    print("Done.", file=sys.stderr)

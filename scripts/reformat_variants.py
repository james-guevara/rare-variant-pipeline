"""
Reformat multi-allelic VCF annotations and add genomic constraint/annotation layers.
Fully lazy/streaming with polars where possible.
"""
import polars as pl
import polars_bio as pb
import argparse
import sys
from pathlib import Path

CONSEQUENTIAL_IMPACTS = {"HIGH", "MODERATE"}


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
    
    if "IMPACT" in lf.collect_schema().names():
        lf_conseq = lf.filter(pl.col("IMPACT").is_in(list(CONSEQUENTIAL_IMPACTS)))
        print(f"Writing {args.consequential}...", file=sys.stderr)
        lf_conseq.sink_csv(args.consequential, separator="\t")
        print(f"Writing {args.consequential_bed}...", file=sys.stderr)
        create_merged_bed(lf_conseq, args.consequential_bed)
    else:
        print("ERROR: IMPACT column not found", file=sys.stderr)
        sys.exit(1)

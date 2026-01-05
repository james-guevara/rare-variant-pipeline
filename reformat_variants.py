"""
Reformat multi-allelic VCF annotations and add genomic constraint/annotation layers.
"""
import polars as pl
import polars_bio as pb
import argparse
import sys
import gzip
from pathlib import Path

CONSEQUENTIAL_IMPACTS = {"HIGH", "MODERATE"}

def write_output(df, path, separator="\t"):
    """Write dataframe to TSV, gzipped if path ends with .gz"""
    if path.endswith(".gz"):
        with gzip.open(path, "wt") as f:
            df.write_csv(f, separator=separator)
    else:
        df.write_csv(path, separator=separator)

def add_gnomad_constraint_by_gene(df: pl.DataFrame, constraint_tsv: str, gene_col: str = "SYMBOL") -> pl.DataFrame:
    cons_gene = (
        pl.read_csv(constraint_tsv, separator="\t", null_values=["NA"],
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
    return df.join(cons_gene, on=gene_col, how="left")

def add_gnomad_constraint_by_transcript(df: pl.DataFrame, constraint_tsv: str, feature_col: str = "Feature") -> pl.DataFrame:
    cons_tx = (
        pl.read_csv(constraint_tsv, separator="\t", null_values=["NA"],
            schema_overrides={"lof.pLI": pl.Float64, "lof.oe_ci.upper": pl.Float64})
        .select([
            pl.col("transcript").alias(feature_col),
            pl.col("lof.pLI").alias("gnomad_lof_pLI_tx"),
            pl.col("lof.oe_ci.upper").alias("gnomad_lof_oe_ci_upper_tx"),
        ])
    )
    return df.join(cons_tx, on=feature_col, how="left")

def add_genebayes(df: pl.DataFrame, genebayes_tsv: str, gene_col: str = "Gene") -> pl.DataFrame:
    gb = pl.read_csv(genebayes_tsv, separator="\t",
        schema_overrides={"obs_lof": pl.Int64, "exp_lof": pl.Float64, "prior_mean": pl.Float64,
            "post_mean": pl.Float64, "post_lower_95": pl.Float64, "post_upper_95": pl.Float64})
    return df.join(
        gb.select(["ensg",
            pl.col("obs_lof").alias("genebayes_obs_lof"),
            pl.col("exp_lof").alias("genebayes_exp_lof"),
            pl.col("prior_mean").alias("genebayes_prior_mean"),
            pl.col("post_mean").alias("genebayes_post_mean"),
            pl.col("post_lower_95").alias("genebayes_post_lower_95"),
            pl.col("post_upper_95").alias("genebayes_post_upper_95"),
        ]),
        left_on=gene_col, right_on="ensg", how="left")

def add_1kb_constraint(df: pl.DataFrame, cons_gz: str) -> pl.DataFrame:
    """Add 1kb genomic constraint scores via overlap."""
    v_iv = df.select([
        pl.col("#CHROM").alias("chrom"),
        pl.col("POS0").cast(pl.Int64).alias("start"),
        pl.col("END").cast(pl.Int64).alias("end"),
    ])
    cons = pl.read_csv(cons_gz, separator="\t",
        schema_overrides={"start": pl.Int64, "end": pl.Int64, "z": pl.Float64, "oe": pl.Float64}
    ).select(["chrom", "start", "end", "z", "oe"])
    
    pairs = pb.overlap(v_iv.lazy(), cons.lazy(), use_zero_based=True, output_type="polars.LazyFrame").collect()
    summ = (
        pairs.group_by(["chrom_1", "start_1", "end_1"])
        .agg([
            pl.col("z_2").max().alias("constraint_z_max"),
            pl.col("z_2").mean().alias("constraint_z_mean"),
            pl.col("oe_2").min().alias("constraint_oe_min"),
        ])
        .rename({"chrom_1": "chrom", "start_1": "start", "end_1": "end"})
    )
    keyed = v_iv.join(summ, on=["chrom", "start", "end"], how="left")
    return pl.concat([df, keyed.select(["constraint_z_max", "constraint_z_mean", "constraint_oe_min"])], how="horizontal")

def add_bed_overlap_flag(df: pl.DataFrame, bed_path: str, prefix: str) -> pl.DataFrame:
    """Add overlap count with BED regions."""
    bed = pl.read_csv(bed_path, separator="\t", has_header=False).select([
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").cast(pl.Int64).alias("start"),
        pl.col("column_3").cast(pl.Int64).alias("end"),
    ])
    v_iv = df.select([
        pl.col("#CHROM").alias("chrom"),
        pl.col("POS0").cast(pl.Int64).alias("start"),
        pl.col("END").cast(pl.Int64).alias("end"),
    ])
    counts = pb.count_overlaps(v_iv.lazy(), bed.lazy(), use_zero_based=True, output_type="polars.LazyFrame").collect()
    return pl.concat([df, counts.select(pl.col("count").alias(prefix))], how="horizontal")

def create_merged_bed(df: pl.DataFrame, output_path: str) -> None:
    merged = pb.merge(
        df.lazy().select([
            pl.col("#CHROM").alias("chrom"),
            pl.col("POS0").cast(pl.Int64).alias("start"),
            pl.col("END").cast(pl.Int64).alias("end"),
        ]),
        use_zero_based=True, min_dist=1, output_type="polars.LazyFrame"
    ).collect().rename({"chrom": "#CHROM", "start": "POS0", "end": "END"})
    merged.write_csv(output_path, separator="\t")

def reformat_variants(tsv_path: str) -> pl.DataFrame:
    df = pl.read_csv(tsv_path, separator="\t", infer_schema_length=0)
    idx = pl.col("ALLELE_NUM").cast(pl.Int32, strict=False).fill_null(1).sub(1)
    per_allele_fields = ["AC", "AF", "AQ", "MLEAC", "MLEAF"]
    
    df = df.with_columns(
        pl.col("ALT").str.split(",").list.get(idx, null_on_oob=True).alias("ALT_specific"),
        *[pl.col("INFO").str.extract(f"{field}=([^;]+)", 1).str.split(",").list.get(idx, null_on_oob=True).fill_null("NA").alias(field)
          for field in per_allele_fields],
    )
    cols = df.columns
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
    return df.select(cols)

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

    # Process step by step to avoid polars_bio DataFusion issues
    print("Reading variants...", file=sys.stderr)
    df = reformat_variants(args.input)
    
    print("Adding segDups overlap...", file=sys.stderr)
    df = add_bed_overlap_flag(df, str(segdups), "segDups")
    
    print("Adding simpleRepeats overlap...", file=sys.stderr)
    df = add_bed_overlap_flag(df, str(simple_repeats), "simpleRepeats")
    
    print("Adding 1kb constraint...", file=sys.stderr)
    df = add_1kb_constraint(df, str(constraint_1kb))
    
    print("Adding GeneBayes...", file=sys.stderr)
    df = add_genebayes(df, str(genebayes))
    
    print("Adding gnomAD transcript constraint...", file=sys.stderr)
    df = add_gnomad_constraint_by_transcript(df, str(constraint))
    
    print("Adding gnomAD gene constraint...", file=sys.stderr)
    df = add_gnomad_constraint_by_gene(df, str(constraint))
    
    # Write outputs
    write_output(df, args.output)
    print(f"Wrote {len(df)} variants to {args.output}", file=sys.stderr)
    
    create_merged_bed(df, args.bed)
    print(f"Wrote BED to {args.bed}", file=sys.stderr)
    
    if "IMPACT" in df.columns:
        df_conseq = df.filter(pl.col("IMPACT").is_in(list(CONSEQUENTIAL_IMPACTS)))
        write_output(df_conseq, args.consequential)
        print(f"Wrote {len(df_conseq)} consequential variants to {args.consequential}", file=sys.stderr)
        create_merged_bed(df_conseq, args.consequential_bed)
        print(f"Wrote consequential BED to {args.consequential_bed}", file=sys.stderr)
    else:
        print("ERROR: IMPACT column not found", file=sys.stderr)
        sys.exit(1)

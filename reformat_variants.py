"""
Reformat multi-allelic VCF annotations and add genomic constraint/annotation layers.

This script processes VEP-annotated variant tables by:
1. Resolving multi-allelic variants (splitting ALT alleles)
2. Extracting per-allele INFO fields
3. Adding gene and transcript constraint metrics (gnomAD, GeneBayes)
4. Adding genomic constraint scores (1kb windows)
5. Flagging overlaps with repetitive regions
"""
import polars as pl
import polars_bio as pb
import argparse
import sys

def add_gnomad_constraint_by_gene(
    df: pl.LazyFrame,
    constraint_tsv: str,
    *,
    gene_col: str = "SYMBOL",
) -> pl.LazyFrame:
    """
    Add gnomAD gene-level constraint metrics.

    Args:
        df: Variant LazyFrame
        constraint_tsv: Path to gnomAD constraint file
        gene_col: Column name containing gene symbols

    Returns:
        LazyFrame with added constraint columns (pLI, LOEUF)
    """
    cons_gene = (
        pl.scan_csv(
            constraint_tsv,
            separator="\t",
            has_header=True,
            null_values=["NA"],
            schema_overrides={"lof.pLI": pl.Float64, "lof.oe_ci.upper": pl.Float64},
        )
        .group_by("gene")
        .agg([
            # "worst"/most constrained-ish summaries (your choices)
            pl.col("lof.pLI").max().alias("gnomad_lof_pLI_max"),
            pl.col("lof.oe_ci.upper").min().alias("gnomad_lof_oe_ci_upper_min"),
            # canonical transcript values (one per gene)
            pl.when(pl.col("canonical") == True).then(pl.col("lof.pLI")).otherwise(None)
              .max().alias("gnomad_lof_pLI_canonical"),
            pl.when(pl.col("canonical") == True).then(pl.col("lof.oe_ci.upper")).otherwise(None)
              .min().alias("gnomad_lof_oe_ci_upper_canonical"),
        ])
        .rename({"gene": gene_col})
    )
    return df.join(cons_gene, on=gene_col, how="left")


def add_gnomad_constraint_by_transcript(
    df: pl.LazyFrame,
    constraint_tsv: str,
    *,
    feature_col: str = "Feature",
) -> pl.LazyFrame:
    cons_tx = (
        pl.scan_csv(
            constraint_tsv,
            separator="\t",
            has_header=True,
            null_values=["NA"],
            schema_overrides={"lof.pLI": pl.Float64, "lof.oe_ci.upper": pl.Float64},
        )
        .select([
            pl.col("transcript").alias(feature_col),
            pl.col("lof.pLI").alias("gnomad_lof_pLI_tx"),
            pl.col("lof.oe_ci.upper").alias("gnomad_lof_oe_ci_upper_tx"),
        ])
    )
    return df.join(cons_tx, on=feature_col, how="left")


def add_genebayes(
    df: pl.LazyFrame,
    genebayes_tsv: str,
    *,
    gene_col: str = "Gene",
) -> pl.LazyFrame:
    gb = pl.scan_csv(
        genebayes_tsv,
        separator="\t",
        has_header=True,
        schema_overrides={
            "obs_lof": pl.Int64,
            "exp_lof": pl.Float64,
            "prior_mean": pl.Float64,
            "post_mean": pl.Float64,
            "post_lower_95": pl.Float64,
            "post_upper_95": pl.Float64,
        },
    )

    return (
        df.join(
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
    )



def add_1kb_constraint(df: pl.LazyFrame, cons_gz: str,
    *,
    chrom_col: str = "#CHROM",
    start_col: str = "POS0",
    end_col: str = "END",
    use_zero_based: bool = True,
) -> pl.LazyFrame:
    # interval view of variants
    v_iv = df.select([
        pl.col(chrom_col).alias("chrom"),
        pl.col(start_col).cast(pl.Int64).alias("start"),
        pl.col(end_col).cast(pl.Int64).alias("end"),
    ])

    # constraint windows
    cons = (
        pl.scan_csv(cons_gz, separator="\t", has_header=True,
            schema_overrides={
                "start": pl.Int64,
                "end": pl.Int64,
                "z": pl.Float64,
                "oe": pl.Float64,
            },
        )
        .select(["chrom", "start", "end", "z", "oe"])
    )

    # overlaps; polars-bio suffixes cols with _1/_2
    pairs = pb.overlap(v_iv, cons, use_zero_based=use_zero_based, output_type="polars.LazyFrame")
    # summarize constraint per variant interval
    summ = (
        pairs.group_by(["chrom_1", "start_1", "end_1"])
             .agg([
                 pl.col("z_2").max().alias("constraint_z_max"),
                 pl.col("z_2").mean().alias("constraint_z_mean"),
                 pl.col("oe_2").min().alias("constraint_oe_min"),
             ])
             .rename({"chrom_1": "chrom", "start_1": "start", "end_1": "end"})
    )
    # join back using the interval view keys
    keyed = v_iv.join(summ, on=["chrom", "start", "end"], how="left")
    # attach to full df (row-aligned, like your bed flag code)
    add_cols = keyed.select(["constraint_z_max", "constraint_z_mean", "constraint_oe_min"])
    return pl.concat([df, add_cols], how="horizontal")




def add_bed_overlap_flag(df: pl.LazyFrame, bed_path: str, prefix: str,
    *,
    chrom_col: str = "#CHROM",
    start_col: str = "POS0",
    end_col: str = "END",
    use_zero_based: bool = True,
) -> pl.LazyFrame:

    bed = (
        pl.scan_csv(bed_path, separator="\t", has_header=False)
          .select([
              pl.col("column_1").alias("chrom"),
              pl.col("column_2").cast(pl.Int64).alias("start"),
              pl.col("column_3").cast(pl.Int64).alias("end"),
          ])
    )

    v_iv = df.select([
        pl.col(chrom_col).alias("chrom"),
        pl.col(start_col).cast(pl.Int64).alias("start"),
        pl.col(end_col).cast(pl.Int64).alias("end"),
    ])

    counts = pb.count_overlaps(
        v_iv,
        bed,
        use_zero_based=use_zero_based,
        output_type="polars.LazyFrame",
    ).select(pl.col("count").alias(f"{prefix}"))

    return pl.concat([df, counts], how="horizontal")


def reformat_variants_lazy(tsv_path: str) -> pl.LazyFrame:
    try:
        df = pl.scan_csv(tsv_path, separator="\t", has_header=True, infer_schema_length=0)

        cols = df.collect_schema().names()

        idx = pl.col("ALLELE_NUM").cast(pl.Int32, strict=False).sub(1)
        per_allele_fields = ["AC","AF","AQ","MLEAC","MLEAF"]
        
        df = df.with_columns(
            pl.col("ALT").str.split(",").list.get(idx).alias("ALT_specific"),
            *[
                pl.col("INFO")
                  .str.extract(f"{field}=([^;]+)", 1)
                  .str.split(",")
                  .list.get(idx)
                  .fill_null("NA")
                  .alias(field)
                for field in per_allele_fields 
            ],
        )

        # reorder columns (must use schema again after adding cols)
        cols = df.collect_schema().names()
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

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)


def reformat_variants(tsv_path: str):
    try:
        # 1. Load the TSV
        # infer_schema_length=0 forces all columns to be read as Strings initially.
        df = pl.read_csv(tsv_path, separator='\t', infer_schema_length=0)
        
        # 2. Identify the Index Column (ALLELE_NUM)
        if "ALLELE_NUM" not in df.columns:
            print(f"Error: Column 'ALLELE_NUM' not found. Available columns: {df.columns}", file=sys.stderr)
            sys.exit(1)

        # Cast ALLELE_NUM to int and adjust to 0-based index
        df = df.with_columns(
            pl.col("ALLELE_NUM")
            .cast(pl.Int32, strict=False)
            .sub(1)       
            .alias("allele_index")
        )

        # 3. Create Specific ALT Column (Do not drop original yet)
        df = df.with_columns(
            pl.col("ALT")
            .str.split(",")
            .list.get(pl.col("allele_index"))
            .alias("ALT_specific")
        )

        # 4. Fix Per-Allele INFO Fields
        # Add any field here. If it doesn't exist in the VCF, it will be filled with "NA".
        per_allele_fields = ["AC", "AF", "AQ", "MLEAC", "MLEAF"] 

        info_expressions = []
        
        for field in per_allele_fields:
            # Regex to find "AC=..." inside the INFO string
            regex_pattern = f"{field}=([^;]+)"
            
            expr = (
                pl.col("INFO")
                .str.extract(regex_pattern, 1)    # Returns null if tag not found
                .str.split(",")                   # Returns null if input is null
                .list.get(pl.col("allele_index")) # Returns null if list is null
                .fill_null("NA")                  # Replaces null with "NA"
                .alias(field)                     # Create new column "AC", "AF" etc.
            )
            info_expressions.append(expr)

        # Apply the INFO extractions
        df = df.with_columns(info_expressions)
        
        # Drop the helper index column
        df = df.drop("allele_index")

        # 6. Reorder Columns for readability
        # We want: ... REF, ALT, ALT_raw, ... INFO, AC, AF, AQ ...
        cols = df.columns

        # Put ALT_specific right after ALT
        if "ALT" in cols and "ALT_specific" in cols:
            cols.remove("ALT_specific")
            cols.insert(cols.index("ALT") + 1, "ALT_specific")

        # Move extracted INFO fields (AC, AF, AQ) to be immediately after INFO
        if "INFO" in cols:
            insert_idx = cols.index("INFO") + 1
            for field in per_allele_fields:
                if field in cols:
                    cols.remove(field)
                    cols.insert(insert_idx, field)
                    insert_idx += 1
        
        df = df.select(cols)

        return df

        # 7. Save output
        #df.write_csv(output_path, separator='\t')
        #print(f"Successfully processed {len(df)} rows. Saved to {output_path}", file=sys.stderr)

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tidy up multi-allelic VCF fields split by VEP.")
    parser.add_argument("input", help="Path to input TSV file")
    parser.add_argument("output", help="Path to output TSV file")

    args = parser.parse_args()
    
    constraint = "resources/gnomAD/constraint/gnomad.v4.1.constraint_metrics.tsv"

    lf = (
        reformat_variants_lazy(args.input)
        .pipe(add_bed_overlap_flag, "resources/repeats/genomicSuperDups.bed", "segDups")
        .pipe(add_bed_overlap_flag, "resources/repeats/simpleRepeat.bed", "simpleRepeats")
        .pipe(add_1kb_constraint, "resources/gnomAD/Genomic_constraint/constraint_z_genome_1kb.qc.download.txt.gz")
        .pipe(add_genebayes, "resources/GeneBayes/output/Supplementary_Table_1.tsv")
        .pipe(add_gnomad_constraint_by_transcript, constraint, feature_col="Feature")
        .pipe(add_gnomad_constraint_by_gene, constraint, gene_col="SYMBOL")
    )

    lf.sink_csv(args.output, separator="\t")

    merged = (
        pb.merge(
            lf.select([
                pl.col("#CHROM").alias("chrom"),
                pl.col("POS0").cast(pl.Int64).alias("start"),
                pl.col("END").cast(pl.Int64).alias("end"),
            ]),
            use_zero_based=True,
            min_dist=1,
            # cols=... appears ignored/buggy here, so don't rely on it
            output_type="polars.LazyFrame",
        )
        .rename({"chrom": "#CHROM", "start": "POS0", "end": "END"})
    ).sink_csv("tests/chr22.variants.merged.bed", separator='\t')




import polars as pl

variants_tsv = "tests/chr22.variants.segdup.tsv"          # your annotated variants
genebayes_tsv = "resources/GeneBayes/output/Supplementary_Table_1.tsv"

# variants (keep wide)
variants = pl.scan_csv(
    variants_tsv,
    separator="\t",
    has_header=True,
)

# GeneBayes gene-level table
genebayes = pl.scan_csv(
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

# If your variants "Gene" sometimes has version suffix like ENSG... .1, strip it.
variants_keyed = variants.with_columns(
    pl.col("Gene").cast(pl.Utf8).str.split(".").list.first().alias("ensg_key")
)

genebayes_keyed = genebayes.with_columns(
    pl.col("ensg").cast(pl.Utf8).alias("ensg_key")
)

# Left join: add GeneBayes columns to each variant
out = (
    variants_keyed
    .join(
        genebayes_keyed.select([
            "ensg_key",
            pl.col("obs_lof").alias("genebayes_obs_lof"),
            pl.col("exp_lof").alias("genebayes_exp_lof"),
            pl.col("prior_mean").alias("genebayes_prior_mean"),
            pl.col("post_mean").alias("genebayes_post_mean"),
            pl.col("post_lower_95").alias("genebayes_post_lower_95"),
            pl.col("post_upper_95").alias("genebayes_post_upper_95"),
        ]),
        on="ensg_key",
        how="left",
    )
    .drop("ensg_key")
)

# quick check
print(out.select(["Gene", "genebayes_post_mean"]).limit(10).collect())

# write
out.sink_csv("tests/chr22.variants.segdup.genebayes.tsv", separator="\t")
# or parquet:
# out.sink_parquet("tests/chr22.variants.segdup.genebayes.parquet")


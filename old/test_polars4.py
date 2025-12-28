import polars as pl

variants_tsv = "tests/chr22.variants.segdup.tsv"  # or your latest
constraint_tsv = "resources/gnomAD/constraint/gnomad.v4.1.constraint_metrics.tsv"

variants = pl.scan_csv(variants_tsv, separator="\t", has_header=True)

cons_tx = (
    pl.scan_csv(
        constraint_tsv,
        separator="\t",
        has_header=True,
        null_values=["NA"],
        schema_overrides={"lof.pLI": pl.Float64, "lof.oe_ci.upper": pl.Float64},
    )
    .select([
        pl.col("transcript").alias("Feature"),
        pl.col("lof.pLI").alias("gnomad_lof_pLI_tx"),
        pl.col("lof.oe_ci.upper").alias("gnomad_lof_oe_ci_upper_tx"),
    ])
)
print(cons_tx.limit(10).collect())
out = variants.join(cons_tx, on="Feature", how="left")
with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=100):
    print(out.select(["#CHROM","POS","SYMBOL","Feature","gnomad_lof_pLI_tx","gnomad_lof_oe_ci_upper_tx"]).limit(1000).collect())


#print(out.limit(100).collect())


cons_worst = (
    pl.scan_csv(
        constraint_tsv,
        separator="\t",
        has_header=True,
        schema_overrides={"lof.pLI": pl.Float64},
        null_values = ["NA"]
    )
    .group_by("gene")
    .agg(
        pl.col("lof.pLI").max().alias("gnomad_lof_pLI_min"),
        pl.col("lof.oe_ci.upper").min().alias("lof.oe_ci.upper_min"),
        pl.when(pl.col("canonical") == True).then(pl.col("lof.pLI")).otherwise(None)
          .max().alias("lof_pLI_canonical"),
        pl.when(pl.col("canonical") == True).then(pl.col("lof.oe_ci.upper")).otherwise(None)
          .min().alias("lof_oe_ci_upper_canonical"),
    )
)

#print(cons_worst.limit(10).collect())


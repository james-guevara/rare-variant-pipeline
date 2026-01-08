import polars as pl
import polars_bio as pb

variants_path = "out/chr22.variants.tsv.gz"
segdup_bed    = "resources/repeats/genomicSuperDups.bed"

variants = pl.scan_csv(
    variants_path,
    separator="\t",
    has_header=True,
    schema_overrides={"POS0": pl.Int64, "END": pl.Int64},
)

segdup6 = pl.scan_csv(
    segdup_bed,
    separator="\t",
    has_header=False,
    new_columns=["chrom","start","end","name","score","strand"],
    schema_overrides={"start": pl.Int64, "end": pl.Int64},
)

# 3-col views for polars-bio range ops
v_iv = variants.select([
    pl.col("#CHROM").alias("chrom"),
    pl.col("POS0").alias("start"),
    pl.col("END").alias("end"),
])

s_iv = segdup6.select(["chrom", "start", "end"])

counts = pb.count_overlaps(
    v_iv,
    s_iv,
    use_zero_based=True,
    output_type="polars.LazyFrame",
)

# Attach overlap counts back to the full variants table (row-aligned)
out = pl.concat(
    [variants, counts.select(pl.col("count").alias("segdup_n"))],
    how="horizontal",
).with_columns(
    (pl.col("segdup_n") > 0).alias("in_segdup")
)

print(out.select(["#CHROM", "POS", "REF", "ALT", "segdup_n", "in_segdup"]).limit(10).collect())

# Optional writes (pick one)
# out.sink_parquet("tests/chr22.variants.segdup.parquet")
out.sink_csv("tests/chr22.variants.segdup.tsv", separator="\t")


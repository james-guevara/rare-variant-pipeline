import polars as pl
import polars_bio as pb

variants_path = "out/chr22.variants.tsv.gz"
segdup_bed    = "resources/repeats/genomicSuperDups.bed"

variants_iv = (
    pl.scan_csv(variants_path, separator="\t", has_header=True)
      .select([
          pl.col("#CHROM").alias("chrom"),
          pl.col("POS0").cast(pl.Int64).alias("start"),
          pl.col("END").cast(pl.Int64).alias("end"),
      ])
      .drop_nulls()
)

#bed = pb.read_bed(segdup_bed)

segdup = (
    pl.scan_csv(segdup_bed, separator="\t", has_header=False, new_columns=["chrom","start","end","name","score","strand"],
    )
    .with_columns([
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    ])
    .select(["chrom","start","end"])
)


counts = pb.count_overlaps(
    variants_iv, segdup,
    cols1=["chrom","start","end"],
    cols2=["chrom","start","end"],
    use_zero_based=True,
    output_type="polars.LazyFrame",
)

print(counts.limit(5).collect())





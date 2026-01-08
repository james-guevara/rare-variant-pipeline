import polars as pl
import polars_bio as pb

variants_tsv = "tests/chr22.variants.segdup.tsv"
cons_gz = "resources/gnomAD/Genomic_constraint/constraint_z_genome_1kb.qc.download.txt.gz"

variants = pl.scan_csv(variants_tsv, separator="\t", has_header=True).rename({"#CHROM":"chrom"})
v_iv = variants.select([
    pl.col("chrom"),
    pl.col("POS0").cast(pl.Int64).alias("start"),
    pl.col("END").cast(pl.Int64).alias("end"),
])


cons = pl.scan_csv(
    cons_gz,
    separator="\t",
    has_header=True,
    schema_overrides={"start": pl.Int64, "end": pl.Int64, "z": pl.Float64, "oe": pl.Float64},
).select(["chrom","start","end","z","oe"])

pairs = pb.overlap(v_iv, cons, use_zero_based=True, output_type="polars.LazyFrame")

#print(pairs.limit(10).collect())

# summarize constraint per variant interval (1kb windows -> typically 1 overlap, but aggregate anyway)
summ = (
    pairs.group_by(["chrom_1","start_1","end_1"])
         .agg([
             pl.col("z_2").max().alias("constraint_z_max"),
             pl.col("z_2").mean().alias("constraint_z_mean"),
             pl.col("oe_2").min().alias("constraint_oe_min"),
         ])
         .rename({"chrom_1":"chrom", "start_1":"POS0", "end_1":"END"})
)

out = variants.join(summ, on=["chrom","POS0","END"], how="left")

with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=100):
    print(out.select(["chrom","POS","constraint_z_max","constraint_oe_min"]).limit(1000).collect())

out.sink_csv("tests/chr22.variants.segdup.constraint.tsv", separator="\t")

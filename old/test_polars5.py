#import polars_bio as pb
#import polars as pl
#
#vcf = pb.scan_vcf("out_vep/chr22.vep.vcf.gz", info_fields=["CSQ"])
#
#out = (
#    vcf.select(["chrom","start","end","ref","alt","CSQ"])   # drop AC/AF/etc
#       .explode("CSQ")
#       .with_columns(pl.col("CSQ").str.split("|").alias("csq_fields"))
#       .select(["chrom","start","ref","alt","CSQ","csq_fields"])
#)
#
#print(out.limit(5).collect())

import polars_bio as pb
import polars as pl

out = (
    pb.scan_vcf("out_vep/chr22.vep.vcf.gz", info_fields=["CSQ"])
      .select(["chrom","start","ref","alt","CSQ"])
      .explode("CSQ")
      .with_columns(pl.col("CSQ").str.split("|").alias("csq_fields"))
      .select([
          "chrom","start","ref","alt","CSQ",
          pl.col("csq_fields").list.get(0).alias("csq_allele"),
          pl.col("csq_fields").list.get(1).alias("csq_consequence"),
          pl.col("csq_fields").list.get(2).alias("csq_impact"),
          # add more indices as needed
      ])
      .limit(1000)
      .collect()
)

out.write_csv("tests/csq_debug.tsv", separator="\t")
print("wrote tests/csq_debug.tsv")


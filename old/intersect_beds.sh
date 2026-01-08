#GAP="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/beds/merged/gap.sorted.merged.bed"
#SEGDUP="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/beds/merged/genomicSuperDups.sorted.merged.bed"
#RMSK="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/beds/merged/rmsk.sorted.merged.bed"
#REPEATS="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/beds/merged/simpleRepeat.sorted.merged.bed"
#CONSTRAINT="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/gnomAD/Genomic_constraint/constraint_z_genome_1kb.qc.download.txt.gz"
#
#
#bedtools intersect -header -a test.bed -b "${CONSTRAINT}" -loj > test.plus_constraint.tsv


## key = first 3 cols (#CHROM, start, end)
#awk 'BEGIN{FS=OFS="\t"}
#     NR==1{next}
#     {key=$1 OFS $2 OFS $3}
#     key!=prev && NR>2 {if(n>1) print prev, n; n=0}
#     {n++; prev=key}
#     END{if(n>1) print prev, n}' test.plus_constraint.tsv \
#| head


awk 'BEGIN{FS=OFS="\t"}
     NR==1{print; next}
     {key=$1 OFS $2 OFS $3}
     key==prev {
        if(n==1) print prevline   # print first line of the group
        print                     # print current line
        n++
     }
     key!=prev { prev=key; prevline=$0; n=1 }' test.plus_constraint.tsv \
> duplicated_rows.tsv


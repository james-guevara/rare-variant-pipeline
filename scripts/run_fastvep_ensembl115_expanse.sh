#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=g2mh_fastvep115
#SBATCH --output=fastvep115_%j.out
#SBATCH --error=fastvep115_%j.err

set -euo pipefail

RUN=/expanse/projects/sebat1/s3/data/sebat/fastvep/runs/g2mh_tiered_chr22
RES=/expanse/projects/sebat1/s3/data/sebat/fastvep/resources/ensembl-115
FASTVEP=/expanse/projects/sebat1/s3/data/sebat/fastvep/fastVEP/target/release/fastvep

cd "${RUN}"
/usr/bin/time -v "${FASTVEP}" annotate \
    --input tiered-gene-alleles.normalized.vcf \
    --gff3 "${RES}/Homo_sapiens.GRCh38.115.chr22.gff3" \
    --fasta "${RES}/Homo_sapiens.GRCh38.dna.chromosome.22.fa" \
    --hgvs \
    --output-format vcf \
    --output tiered-gene-alleles.fastvep115-all-transcripts.vcf

grep -vc '^#' tiered-gene-alleles.fastvep115-all-transcripts.vcf
ls -lh tiered-gene-alleles.fastvep115-all-transcripts.vcf

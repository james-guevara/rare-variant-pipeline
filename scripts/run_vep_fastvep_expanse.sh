#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=g2mh_vep_fastvep
#SBATCH --output=annotation_%j.out
#SBATCH --error=annotation_%j.err

set -euo pipefail

RUN=/expanse/projects/sebat1/s3/data/sebat/fastvep/runs/g2mh_tiered_chr22
RES=/expanse/projects/sebat1/s3/data/sebat/fastvep/resources/ensembl-115
PIPE=/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline
INPUT_GZ=${RUN}/tiered-gene-alleles.normalized.vcf.gz
INPUT=${RUN}/tiered-gene-alleles.normalized.vcf
GFF3=${RES}/Homo_sapiens.GRCh38.115.gff3.gz
FASTA=${RES}/Homo_sapiens.GRCh38.dna.chromosome.22.fa
FASTVEP=/expanse/projects/sebat1/s3/data/sebat/fastvep/fastVEP/target/release/fastvep

module load cpu/0.17.3b singularitypro
cd "${RUN}"

if [[ ! -s "${INPUT}" ]]; then
    gzip -dc "${INPUT_GZ}" > "${INPUT}"
fi

echo "=== standard VEP 115: one picked consequence per allele ==="
/usr/bin/time -v singularity exec -B /expanse "${PIPE}/ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools" vep \
    --input_file "${INPUT}" \
    --format vcf \
    --output_file tiered-gene-alleles.vep115.vcf.gz \
    --vcf --compress_output bgzip \
    --canonical --mane --symbol --protein --hgvs --hgvsg --domains --biotype \
    --numbers --uploaded_allele --pick_allele --transcript_version \
    --assembly GRCh38 --cache --dir_cache "${PIPE}/VEP_CACHE" --offline \
    --fasta "${FASTA}" \
    --force_overwrite --fork "${SLURM_CPUS_PER_TASK}"

echo "=== fastVEP: all transcripts; do not use its non-equivalent --pick ==="
/usr/bin/time -v "${FASTVEP}" annotate \
    --input "${INPUT}" \
    --gff3 "${GFF3}" \
    --fasta "${FASTA}" \
    --hgvs \
    --output-format vcf \
    --output tiered-gene-alleles.fastvep-all-transcripts.vcf

echo "=== outputs ==="
bcftools index -t -f tiered-gene-alleles.vep115.vcf.gz
bcftools index -n tiered-gene-alleles.vep115.vcf.gz
grep -vc '^#' tiered-gene-alleles.fastvep-all-transcripts.vcf
ls -lh tiered-gene-alleles.vep115.vcf.gz* tiered-gene-alleles.fastvep-all-transcripts.vcf

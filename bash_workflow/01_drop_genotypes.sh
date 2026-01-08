#!/bin/bash
#SBATCH --job-name=drop_gt
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/01_drop_genotypes_%j.out
#SBATCH --error=logs/01_drop_genotypes_%j.err

set -euo pipefail

SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/bash_workflow"
source "${SCRIPT_DIR}/config.sh"

echo "=== Step 1: Drop Genotypes ==="
echo "Input: ${INPUT_VCF}"
echo "Output: ${SITES_VCF}"

mkdir -p "$(dirname "${SITES_VCF}")"
mkdir -p "${SCRIPT_DIR}/logs"

module load singularitypro

singularity exec --bind /expanse/projects/sebat1/ "${BCFTOOLS_CONTAINER}" \
    bcftools view -G --threads 4 -O z -o "${SITES_VCF}" "${INPUT_VCF}"

singularity exec --bind /expanse/projects/sebat1/ "${BCFTOOLS_CONTAINER}" \
    tabix -p vcf "${SITES_VCF}"

echo "Done: $(date)"

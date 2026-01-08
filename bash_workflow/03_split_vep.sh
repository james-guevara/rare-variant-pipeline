#!/bin/bash
#SBATCH --job-name=split_vep
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/03_split_vep_%j.out
#SBATCH --error=logs/03_split_vep_%j.err

set -euo pipefail

SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/bash_workflow"
source "${SCRIPT_DIR}/config.sh"

echo "=== Step 3: Split VEP ==="
echo "Input: ${VEP_VCF}"
echo "Output: ${VARIANTS_TSV}"

mkdir -p "$(dirname "${VARIANTS_TSV}")"

module load singularitypro

singularity exec --bind /expanse/projects/sebat1/ "${BCFTOOLS_CONTAINER}" \
    bcftools +split-vep -p CSQ -HH -A $'\t' -d -s :missense+ \
    -f '%CHROM\t%POS0\t%END\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t%INFO\t%CSQ\n' \
    "${VEP_VCF}" | \
    sed -E '2,$s/;?CSQ[^=]*=[^;\t]*//g' | \
    sed -E '1s/CSQ//g; 1s/\[[0-9]+\]//g; 1s/\(null\)/INFO/g' | \
    bgzip > "${VARIANTS_TSV}"

echo "Done: $(date)"

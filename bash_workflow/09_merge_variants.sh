#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/09_merge_variants_%j.out
#SBATCH --error=logs/09_merge_variants_%j.err

set -euo pipefail

SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/bash_workflow"
source "${SCRIPT_DIR}/config.sh"

echo "=== Step 9: Merge Variants ==="
echo "Input family: ${RESOLVED_GENOTYPES}"
echo "Input variants: ${CONSEQ_TSV}"
echo "Output: ${MERGED_TSV}"

mkdir -p "$(dirname "${MERGED_TSV}")"

micromamba run -n python3.12_env_default python "${MERGE_SCRIPT}" \
    --family "${RESOLVED_GENOTYPES}" \
    --variants "${CONSEQ_TSV}" \
    --out "${MERGED_TSV}"

echo "Done: $(date)"
echo "Final output: ${MERGED_TSV}"

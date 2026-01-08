#!/bin/bash
#SBATCH --job-name=gather
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/07_gather_genotypes_%j.out
#SBATCH --error=logs/07_gather_genotypes_%j.err

set -euo pipefail
set +o pipefail  # Ignore SIGPIPE from head/tail

SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/bash_workflow"
source "${SCRIPT_DIR}/config.sh"

echo "=== Step 7: Gather Genotypes ==="
echo "Input: ${SCATTER_DIR}/chunk_*.genotypes.tsv.gz"
echo "Output: ${FAMILY_GENOTYPES}"

mkdir -p "$(dirname "${FAMILY_GENOTYPES}")"

# Get header from first file
zcat "${SCATTER_DIR}/chunk_1.genotypes.tsv.gz" | head -1 | gzip > "${FAMILY_GENOTYPES}"

# Concatenate all files (skip headers)
for i in $(seq 1 ${N_CHUNKS}); do
    zcat "${SCATTER_DIR}/chunk_${i}.genotypes.tsv.gz" | tail -n +2
done | gzip >> "${FAMILY_GENOTYPES}"

echo "Total lines: $(zcat "${FAMILY_GENOTYPES}" | wc -l)"
echo "Done: $(date)"

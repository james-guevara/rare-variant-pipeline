#!/bin/bash
#SBATCH --job-name=scatter
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=logs/05_scatter_bed_%j.out
#SBATCH --error=logs/05_scatter_bed_%j.err

set -euo pipefail

SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/bash_workflow"
source "${SCRIPT_DIR}/config.sh"

echo "=== Step 5: Scatter BED ==="
echo "Input: ${CONSEQ_BED}"
echo "Output: ${SCATTER_DIR}/chunk_*.bed"

mkdir -p "${SCATTER_DIR}"

# Count lines
total=$(wc -l < "${CONSEQ_BED}")
chunk_size=$(( (total + N_CHUNKS - 1) / N_CHUNKS ))

echo "Total lines: ${total}, Chunk size: ${chunk_size}, N_CHUNKS: ${N_CHUNKS}"

# Split into chunks
for i in $(seq 1 ${N_CHUNKS}); do
    start=$(( (i - 1) * chunk_size + 1 ))
    end=$(( i * chunk_size ))
    out_file="${SCATTER_DIR}/chunk_${i}.bed"
    
    sed -n "${start},${end}p" "${CONSEQ_BED}" > "${out_file}"
    lines=$(wc -l < "${out_file}")
    echo "  chunk_${i}.bed: ${lines} lines"
done

echo "Done: $(date)"

#!/bin/bash
#SBATCH --job-name=reformat
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/04_reformat_variants_%j.out
#SBATCH --error=logs/04_reformat_variants_%j.err

set -euo pipefail

SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/bash_workflow"
source "${SCRIPT_DIR}/config.sh"

echo "=== Step 4: Reformat Variants ==="
echo "Input: ${VARIANTS_TSV}"
echo "Output: ${REFORMATTED_TSV}, ${CONSEQ_TSV}, ${CONSEQ_BED}"

mkdir -p "$(dirname "${REFORMATTED_TSV}")"

micromamba run -n python3.12_env_default python "${REFORMAT_SCRIPT}" \
    "${VARIANTS_TSV}" "${REFORMATTED_TSV}" \
    --bed "${REFORMATTED_BED}" \
    --consequential "${CONSEQ_TSV}" \
    --consequential-bed "${CONSEQ_BED}" \
    --resources-dir "${PYTHON_RESOURCES}"

echo "Done: $(date)"

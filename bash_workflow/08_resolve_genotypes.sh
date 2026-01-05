#!/bin/bash
#SBATCH --job-name=resolve
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/08_resolve_genotypes_%j.out
#SBATCH --error=logs/08_resolve_genotypes_%j.err

set -euo pipefail

SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/bash_workflow"
source "${SCRIPT_DIR}/config.sh"

echo "=== Step 8: Resolve Genotypes ==="
echo "Input: ${FAMILY_GENOTYPES}"
echo "Output: ${RESOLVED_GENOTYPES}"

micromamba run -n python3.12_env_default python "${RESOLVE_SCRIPT}" \
    "${FAMILY_GENOTYPES}" "${RESOLVED_GENOTYPES}"

echo "Done: $(date)"

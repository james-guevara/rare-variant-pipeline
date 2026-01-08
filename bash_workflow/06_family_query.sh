#!/bin/bash
#SBATCH --job-name=fam_query
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-10
#SBATCH --output=logs/06_family_query_%A_%a.out
#SBATCH --error=logs/06_family_query_%A_%a.err

set -euo pipefail

SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/bash_workflow"
source "${SCRIPT_DIR}/config.sh"

CHUNK="${SLURM_ARRAY_TASK_ID}"
BED_FILE="${SCATTER_DIR}/chunk_${CHUNK}.bed"
OUT_FILE="${SCATTER_DIR}/chunk_${CHUNK}.genotypes.tsv.gz"

echo "=== Step 6: Family Query (chunk ${CHUNK}) ==="
echo "Input BED: ${BED_FILE}"
echo "Input VCF: ${INPUT_VCF}"
echo "Output: ${OUT_FILE}"

if [ -s "${BED_FILE}" ]; then
    micromamba run -n python3.12_env_default python "${FAMILY_QUERY_SCRIPT}" \
        --vcf "${INPUT_VCF}" \
        --ped "${PED_FILE}" \
        --region "${BED_FILE}" \
        --out /dev/stdout | gzip > "${OUT_FILE}"
else
    echo "Empty BED file, creating header-only output"
    echo -e "#CHROM\tPOS0\tEND\tREF\tALT\tFID\tIID\tGT\tGQ\tDP\tAD" | gzip > "${OUT_FILE}"
fi

echo "Done: $(date)"

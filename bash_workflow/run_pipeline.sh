#!/bin/bash
# Rare Variant Pipeline - Orchestrator
# Usage: ./run_pipeline.sh [CHROM]
# Example: ./run_pipeline.sh chrY
#          ./run_pipeline.sh chr22

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set chromosome from argument or default
CHROM="${1:-chrY}"
export CHROM

source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "Rare Variant Pipeline"
echo "Chromosome: ${CHROM}"
echo "Output: ${OUTPUT_DIR}"
echo "=============================================="

mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}"

# Submit jobs with dependencies (--export passes CHROM to each job)
echo "Submitting jobs..."

# Step 1: Drop genotypes
JOB1=$(sbatch --parsable --export=CHROM="${CHROM}" "${SCRIPT_DIR}/01_drop_genotypes.sh")
echo "Step 1 (drop_genotypes): Job ${JOB1}"

# Step 2: VEP annotate
JOB2=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB1} "${SCRIPT_DIR}/02_vep_annotate.sh")
echo "Step 2 (vep_annotate): Job ${JOB2}"

# Step 3: Split VEP
JOB3=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB2} "${SCRIPT_DIR}/03_split_vep.sh")
echo "Step 3 (split_vep): Job ${JOB3}"

# Step 4: Reformat variants
JOB4=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB3} "${SCRIPT_DIR}/04_reformat_variants.sh")
echo "Step 4 (reformat_variants): Job ${JOB4}"

# Step 5: Scatter BED
JOB5=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB4} "${SCRIPT_DIR}/05_scatter_bed.sh")
echo "Step 5 (scatter_bed): Job ${JOB5}"

# Step 6: Family query (array job)
JOB6=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB5} "${SCRIPT_DIR}/06_family_query.sh")
echo "Step 6 (family_query): Job ${JOB6} (array 1-10)"

# Step 7: Gather genotypes
JOB7=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB6} "${SCRIPT_DIR}/07_gather_genotypes.sh")
echo "Step 7 (gather_genotypes): Job ${JOB7}"

# Step 8: Resolve genotypes
JOB8=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB7} "${SCRIPT_DIR}/08_resolve_genotypes.sh")
echo "Step 8 (resolve_genotypes): Job ${JOB8}"

# Step 9: Merge variants
JOB9=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB8} "${SCRIPT_DIR}/09_merge_variants.sh")
echo "Step 9 (merge_variants): Job ${JOB9}"

echo ""
echo "=============================================="
echo "Pipeline submitted!"
echo "Monitor with: squeue -u \$USER"
echo "Final output: ${MERGED_TSV}"
echo "=============================================="

# Save job IDs for reference
cat > "${SCRIPT_DIR}/logs/job_ids_${CHROM}.txt" << JOBS
CHROM=${CHROM}
JOB1_drop_genotypes=${JOB1}
JOB2_vep_annotate=${JOB2}
JOB3_split_vep=${JOB3}
JOB4_reformat_variants=${JOB4}
JOB5_scatter_bed=${JOB5}
JOB6_family_query=${JOB6}
JOB7_gather_genotypes=${JOB7}
JOB8_resolve_genotypes=${JOB8}
JOB9_merge_variants=${JOB9}
JOBS

echo "Job IDs saved to: logs/job_ids_${CHROM}.txt"

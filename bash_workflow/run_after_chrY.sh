#!/bin/bash
# Submit chr20, chr21, chr22 after chrY completes successfully
# Usage: ./run_after_chrY.sh <chrY_final_job_id>

CHRY_JOB="${1:-45645589}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Will submit chr20, chr21, chr22 after chrY job ${CHRY_JOB} completes"

for CHROM in chr20 chr21 chr22; do
    echo "Submitting ${CHROM} pipeline with dependency on ${CHRY_JOB}..."
    
    export CHROM
    source "${SCRIPT_DIR}/config.sh"
    
    # Submit step 1 with dependency on chrY final job
    JOB1=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${CHRY_JOB} "${SCRIPT_DIR}/01_drop_genotypes.sh")
    
    # Chain the rest
    JOB2=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB1} "${SCRIPT_DIR}/02_vep_annotate.sh")
    JOB3=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB2} "${SCRIPT_DIR}/03_split_vep.sh")
    JOB4=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB3} "${SCRIPT_DIR}/04_reformat_variants.sh")
    JOB5=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB4} "${SCRIPT_DIR}/05_scatter_bed.sh")
    JOB6=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB5} "${SCRIPT_DIR}/06_family_query.sh")
    JOB7=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB6} "${SCRIPT_DIR}/07_gather_genotypes.sh")
    JOB8=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB7} "${SCRIPT_DIR}/08_resolve_genotypes.sh")
    JOB9=$(sbatch --parsable --export=CHROM="${CHROM}" --dependency=afterok:${JOB8} "${SCRIPT_DIR}/09_merge_variants.sh")
    
    echo "  ${CHROM}: Jobs ${JOB1} -> ${JOB9}"
    
    # Update CHRY_JOB for next chromosome to chain them
    CHRY_JOB="${JOB9}"
done

echo "All chromosomes submitted. Monitor with: squeue -u \$USER"

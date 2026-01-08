#!/bin/bash
#SBATCH --job-name=vep
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/02_vep_annotate_%j.out
#SBATCH --error=logs/02_vep_annotate_%j.err

set -euo pipefail

SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/bash_workflow"
source "${SCRIPT_DIR}/config.sh"

echo "=== Step 2: VEP Annotation ==="
echo "Input: ${SITES_VCF}"
echo "Output: ${VEP_VCF}"

mkdir -p "$(dirname "${VEP_VCF}")"

module load singularitypro

singularity exec --bind /expanse/projects/sebat1/ "${VEP_CONTAINER}" vep \
    --input_file "${SITES_VCF}" \
    --format vcf \
    --output_file "${VEP_VCF}" \
    --vcf \
    --compress_output bgzip \
    --minimal \
    --canonical \
    --mane \
    --symbol \
    --protein \
    --pubmed \
    --nearest symbol \
    --uploaded_allele \
    --numbers \
    --assembly GRCh38 \
    --cache \
    --dir_cache "${VEP_CACHE}" \
    --offline \
    --fasta "${REFERENCE}" \
    --allele_number \
    --pick_allele \
    --regulatory \
    --biotype \
    --domains \
    --force_overwrite \
    --fork 8 \
    --stats_text \
    --dir_plugins "${VEP_PLUGINS}" \
    --plugin dbNSFP,"${RESOURCES_BASE}/dbNSFP/dbNSFP5.3a_grch38.gz",transcript_match=1,MPC_score,MPC_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,ClinPred_score,ClinPred_rankscore,ClinPred_pred,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,1000Gp3_AC,1000Gp3_AF,AllofUs_ALL_AF,AllofUs_POPMAX_AF,AllofUs_POPMAX_POP,RegeneronME_ALL_AF,gnomAD4.1_joint_flag,gnomAD4.1_joint_AF,gnomAD4.1_joint_nhomalt,gnomAD4.1_joint_POPMAX_AF,gnomAD4.1_joint_POPMAX_nhomalt,ALFA_Total_AF,dbNSFP_POPMAX_AF,dbNSFP_POPMAX_AC,dbNSFP_POPMAX_POP,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,GERP_92_mammals,GERP_92_mammals_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP470way_mammalian,phyloP470way_mammalian_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons470way_mammalian,phastCons470way_mammalian_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore \
    --plugin LoF,loftee_path:"${LOFTEE_PATH}",human_ancestor_fa:"${RESOURCES_BASE}/LOFTEE/human_ancestor.fa.gz",gerp_bigwig:"${RESOURCES_BASE}/LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw",conservation_file:"${RESOURCES_BASE}/LOFTEE/loftee.sql" \
    --plugin SpliceAI,snv="${RESOURCES_BASE}/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz",indel="${RESOURCES_BASE}/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz" \
    --plugin MaxEntScan,"${RESOURCES_BASE}/MaxEntScan/fordownload",SWA,NCSS

singularity exec --bind /expanse/projects/sebat1/ "${VEP_CONTAINER}" \
    tabix -p vcf "${VEP_VCF}"

echo "Done: $(date)"

#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=12
#SBATCH --time=16:00:00
#SBATCH -J vep_annotation
#SBATCH -o logs_vep/py_%A_%a.out
#SBATCH -e logs_vep/py_%A_%a.err

SECONDS=0 

module load singularitypro
mkdir -p logs_vep
mkdir -p out_vep

set -euo pipefail

CHR_INDEX=${SLURM_ARRAY_TASK_ID}
# map array index to chromosome label
if   [[ $CHR_INDEX -le 22 ]]; then
    CHR=$CHR_INDEX
elif [[ $CHR_INDEX -eq 23 ]]; then
    CHR="X"
elif [[ $CHR_INDEX -eq 24 ]]; then
    CHR="Y"
fi

VCF="vcfs/chr${CHR}_jointcall_VQSR_combined.vcf.gz"

GRCH38_REF="/expanse/projects/sebat1/j3guevar/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
PLUGIN="VEP_PLUGINS_ALL"
CACHE="VEP_CACHE"
VEP_CONTAINER="ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools"
LOFTEE="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/VEP_PLUGINS/loftee"

# Make VCF with no genotypes (which will hold annotations)
bcftools view --threads $(nproc) -W -I -O z -o out_vep/chr${CHR}.vcf.gz -G $VCF 

# Run VEP
singularity exec --bind /expanse/projects/sebat1/ ${VEP_CONTAINER} \
	vep \
	--input_file out_vep/chr${CHR}.vcf.gz \
	--format vcf \
	--output_file out_vep/chr${CHR}.vep.vcf.gz \
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
	--dir_cache ${CACHE} \
	--offline \
	--fasta ${GRCH38_REF} \
	--allele_number \
	--pick_allele \
    --regulatory \
	--biotype \
	--domains \
	--force_overwrite \
	--fork 4 \
	--stats_text \
	--dir_plugins ${PLUGIN} \
	--plugin dbNSFP,/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/dbNSFP/dbNSFP5.3a_grch38.gz,transcript_match=1,MPC_score,MPC_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,ClinPred_score,ClinPred_rankscore,ClinPred_pred,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,1000Gp3_AC,1000Gp3_AF,AllofUs_ALL_AF,AllofUs_POPMAX_AF,AllofUs_POPMAX_POP,RegeneronME_ALL_AF,gnomAD4.1_joint_flag,gnomAD4.1_joint_AF,gnomAD4.1_joint_nhomalt,gnomAD4.1_joint_POPMAX_AF,gnomAD4.1_joint_POPMAX_nhomalt,ALFA_Total_AF,dbNSFP_POPMAX_AF,dbNSFP_POPMAX_AC,dbNSFP_POPMAX_POP,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,GERP_92_mammals,GERP_92_mammals_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP470way_mammalian,phyloP470way_mammalian_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons470way_mammalian,phastCons470way_mammalian_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore \
    --plugin LoF,loftee_path:${LOFTEE},human_ancestor_fa:/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/LOFTEE/human_ancestor.fa.gz,gerp_bigwig:/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/LOFTEE/loftee.sql \
	--plugin SpliceAI,snv=/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
	--plugin MaxEntScan,/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/MaxEntScan/fordownload,SWA,NCSS

echo "Job finished in $SECONDS seconds."

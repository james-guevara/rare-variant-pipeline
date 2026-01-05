#!/bin/bash
# Rare Variant Pipeline - Shared Configuration

# Chromosome to process (override with: CHROM=chr22 source config.sh)
CHROM="${CHROM:-chrY}"

# Number of chunks for scatter/gather
N_CHUNKS=10

# Base directories
PIPELINE_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline"
SCRIPTS_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline"
RESOURCES_BASE="${SCRIPTS_DIR}/resources"
OUTPUT_DIR="${PIPELINE_DIR}/bash_workflow/output"

# Input VCF
VCF_DIR="/expanse/projects/sebat1/s3/data/sebat/SPARK_iWES_v3/variants/snv/deepvariant/pvcf"
INPUT_VCF="${VCF_DIR}/SPARK.iWES_v3.2024_08.deepvariant.${CHROM}.vcf.gz"

# Containers (sandbox directories)
BCFTOOLS_CONTAINER="${SCRIPTS_DIR}/bcftools:1.22--h3a4d415_1"
VEP_CONTAINER="${SCRIPTS_DIR}/ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools"

# VEP resources
VEP_CACHE="${PIPELINE_DIR}/VEP_CACHE"
VEP_PLUGINS="${SCRIPTS_DIR}/VEP_PLUGINS_ALL"
LOFTEE_PATH="${SCRIPTS_DIR}/VEP_PLUGINS/loftee"
REFERENCE="/expanse/projects/sebat1/j3guevar/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# Python scripts
REFORMAT_SCRIPT="${PIPELINE_DIR}/reformat_variants.py"
FAMILY_QUERY_SCRIPT="${PIPELINE_DIR}/family_query.py"
RESOLVE_SCRIPT="${PIPELINE_DIR}/resolve_family_genotypes.py"
MERGE_SCRIPT="${PIPELINE_DIR}/merge_family_variants.py"

# Other resources
PED_FILE="${PIPELINE_DIR}/resources/SPARK_iWES_v3.ped"
PYTHON_RESOURCES="${PIPELINE_DIR}/resources"

# Output paths (derived)
SITES_VCF="${OUTPUT_DIR}/sites_only/${CHROM}.sites.vcf.gz"
VEP_VCF="${OUTPUT_DIR}/vep/${CHROM}.vep.vcf.gz"
VARIANTS_TSV="${OUTPUT_DIR}/variants/${CHROM}.variants.tsv.gz"
REFORMATTED_TSV="${OUTPUT_DIR}/reformatted/${CHROM}.reformatted.tsv.gz"
REFORMATTED_BED="${OUTPUT_DIR}/reformatted/${CHROM}.bed"
CONSEQ_TSV="${OUTPUT_DIR}/reformatted/${CHROM}.consequential.tsv.gz"
CONSEQ_BED="${OUTPUT_DIR}/reformatted/${CHROM}.consequential.bed"
SCATTER_DIR="${OUTPUT_DIR}/scatter/${CHROM}"
FAMILY_GENOTYPES="${OUTPUT_DIR}/genotypes/${CHROM}.family_genotypes.tsv.gz"
RESOLVED_GENOTYPES="${OUTPUT_DIR}/genotypes/${CHROM}.resolved_genotypes.tsv.gz"
MERGED_TSV="${OUTPUT_DIR}/merged/${CHROM}.merged.tsv.gz"

# SLURM defaults
SLURM_ACCOUNT="ddp195"
SLURM_PARTITION="shared"

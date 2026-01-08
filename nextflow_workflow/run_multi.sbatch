#!/bin/bash
#SBATCH --job-name=nf_multi
#SBATCH --partition=shared
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/nf_multi_%j.out
#SBATCH --error=logs/nf_multi_%j.err

module load singularitypro
eval "$(micromamba shell hook --shell bash)"
micromamba activate nf_env

cd /expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/nextflow_workflow
nextflow run main.nf --chroms chr20,chr21 --outdir output -resume

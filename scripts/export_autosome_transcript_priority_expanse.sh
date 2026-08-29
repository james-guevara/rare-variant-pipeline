#!/usr/bin/env bash
#SBATCH --partition=ind-shared
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=3-20%6
#SBATCH --job-name=rvp_tx_priority
#SBATCH --output=/expanse/projects/sebat1/resources/rare-variant-pipeline/logs/tx-priority-%A_%a.out
set -euo pipefail

chromosome=chr${SLURM_ARRAY_TASK_ID}
cache_root=/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/VEP_CACHE/homo_sapiens/115_GRCh38
vep_image=/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools
registry=/expanse/projects/sebat1/resources/rare-variant-pipeline
repo=/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-controller-41de024
output_dir=$registry/annotation/ensembl-115/$chromosome

module load cpu/0.17.3b singularitypro/4.1.2
mkdir -p "$output_dir" "$registry/logs"
bash "$repo/scripts/export_vep_transcript_priority.sh" \
  "$cache_root/${chromosome#chr}" "$vep_image" \
  "$output_dir/vep115.$chromosome.transcript-priority.tsv" singularity

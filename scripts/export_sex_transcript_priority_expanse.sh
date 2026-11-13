#!/usr/bin/env bash
#SBATCH --job-name=sex-tx-priority
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/expanse/projects/sebat1/j3guevar/sex-resource-prep/logs/%x-%j.out
#SBATCH --error=/expanse/projects/sebat1/j3guevar/sex-resource-prep/logs/%x-%j.err

set -euo pipefail
module load singularitypro/4.1.2

repo=${RVP_REPO:-/expanse/projects/sebat1/j3guevar/sex-resource-prep/repo}
cache_root=/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/VEP_CACHE/homo_sapiens/115_GRCh38
vep_image=/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools
output_root=/expanse/projects/sebat1/j3guevar/sex-resource-prep/output

mkdir -p "$output_root"
for contig in X Y; do
  "$repo/scripts/export_vep_transcript_priority.sh" \
    "$cache_root/$contig" "$vep_image" \
    "$output_root/vep115.chr${contig}.transcript-priority.tsv" singularity
done
sha256sum "$output_root"/*.tsv > "$output_root/SHA256SUMS"

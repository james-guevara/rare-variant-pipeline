#!/usr/bin/env bash
#SBATCH --partition=ind-shared
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=3-20%4
#SBATCH --job-name=rvp_annotation
#SBATCH --output=/expanse/projects/sebat1/resources/rare-variant-pipeline/logs/annotation-%A_%a.out
set -euo pipefail

contig=$SLURM_ARRAY_TASK_ID
chromosome=chr$contig
registry=/expanse/projects/sebat1/resources/rare-variant-pipeline
source_root=$registry/references/ensembl-115-source
annotation_root=$registry/annotation/ensembl-115/$chromosome
image=$registry/containers/rare-variant-targeted-sha-41de024.sif
fasta_gz=$source_root/Homo_sapiens.GRCh38.dna.chromosome.$contig.fa.gz

mkdir -p "$annotation_root" "$registry/logs"
if ! test -s "$fasta_gz"; then
  curl --fail --location --retry 5 --output "$fasta_gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/$(basename "$fasta_gz")"
fi

module load cpu/0.17.3b singularitypro/4.1.2
singularity exec -B /expanse "$image" env \
  ENSEMBL_SOURCE_ROOT="$source_root" \
  ANNOTATION_ROOT="$annotation_root" \
  LOFTEE_ROOT="$registry/loftee" \
  ENSEMBL_GTF="$source_root/Homo_sapiens.GRCh38.115.chr.gtf.gz" \
  FASTVEP=/opt/fastvep/fastvep \
  RVP_PYTHON=/opt/rvp/.venv/bin/python \
  RVP_REPO=/opt/rvp \
  bash /opt/rvp/scripts/prepare_ensembl_chromosome_resources.sh \
    "$chromosome" "$contig"

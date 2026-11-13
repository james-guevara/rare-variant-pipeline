#!/usr/bin/env bash
#SBATCH --partition=ind-shared
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=2-20%4
#SBATCH --job-name=rvp_candidates
#SBATCH --output=/expanse/projects/sebat1/resources/rare-variant-pipeline/logs/candidates-%A_%a.out
set -euo pipefail

contig=$SLURM_ARRAY_TASK_ID
chromosome=chr$contig
registry=/expanse/projects/sebat1/resources/rare-variant-pipeline
bundle=$registry/candidate-bundles/genebayes-dbnsfp-5.3.1a-v1/$chromosome
image=$registry/containers/rare-variant-targeted-sha-41de024.sif
dbnsfp=/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/parquet_expanded_mane_select/$chromosome.parquet
gtf=$registry/references/ensembl-115-source/Homo_sapiens.GRCh38.115.chr.gtf.gz
genebayes=$registry/genebayes/GeneBayes.Supplementary_Table_1.tsv

mkdir -p "$bundle" "$registry/logs"
lof_bed=$bundle/lof-genebayes-ge-0.03.$chromosome.bed
target_bed=$bundle/lof-plus-missense-candidates.$chromosome.bed
candidates=$bundle/$chromosome.all-missense-candidates.parquet

module load cpu/0.17.3b singularitypro/4.1.2
run() {
  singularity exec -B /expanse "$image" /opt/rvp/.venv/bin/python "$@"
}

test -s "$lof_bed" || run /opt/rvp/scripts/build_target_bed.py \
  --genebayes "$genebayes" --gtf "$gtf" --output "$lof_bed" \
  --min-post-mean 0.03 --features exon --padding 8 --chroms "$contig" \
  --add-chr-prefix
test -s "$target_bed" || run /opt/rvp/scripts/build_missense_candidate_bed.py \
  --dbnsfp "$dbnsfp" --union-bed "$lof_bed" --output "$target_bed" \
  --padding 8 --add-chr-prefix
test -s "$candidates" || run /opt/rvp/scripts/build_missense_candidate_alleles.py \
  --dbnsfp "$dbnsfp" --all-genes --chrom "$chromosome" --output "$candidates"

sha256sum "$lof_bed" "$target_bed" "$candidates" > "$bundle/SHA256SUMS"
printf 'candidate_bundle_ready\t%s\tbed=%s\tcandidates=%s\n' \
  "$chromosome" "$(stat -c %s "$target_bed")" "$(stat -c %s "$candidates")"

#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 8 ]]; then
  echo "usage: $0 CHROM CONTIG DBNSFP ZARR GTF GENEBAYES OUTPUT_DIR PYTHON" >&2
  exit 2
fi

chromosome=$1
contig=$2
dbnsfp=$3
zarr_store=$4
gtf=$5
genebayes=$6
output_dir=$7
python=$8
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

mkdir -p "$output_dir"
lof_bed=$output_dir/lof-genebayes-ge-0.03.$chromosome.bed
target_bed=$output_dir/lof-plus-missense-candidates.$chromosome.bed
target_alleles=$output_dir/$chromosome.target-alleles.parquet
all_candidates=$output_dir/$chromosome.all-missense-candidates.parquet
observed_candidates=$output_dir/$chromosome.observed-missense-candidates.parquet

run_if_missing() {
  local output=$1
  shift
  if [[ -s $output ]]; then
    printf 'checkpoint\tSKIP\t%s\n' "$output"
    return
  fi
  "$@"
  [[ -s $output ]]
  printf 'checkpoint\tDONE\t%s\n' "$output"
}

run_if_missing "$lof_bed" "$python" "$repo/scripts/build_target_bed.py" \
  --genebayes "$genebayes" --gtf "$gtf" --output "$lof_bed" \
  --min-post-mean 0.03 --features exon --padding 8 --chroms "$contig" \
  --add-chr-prefix

run_if_missing "$target_bed" "$python" "$repo/scripts/build_missense_candidate_bed.py" \
  --dbnsfp "$dbnsfp" --union-bed "$lof_bed" --output "$target_bed" \
  --padding 8 --add-chr-prefix

run_if_missing "$target_alleles" "$python" "$repo/scripts/extract_zarr_target_alleles.py" \
  --zarr "$zarr_store" --bed "$target_bed" --chrom "$chromosome" \
  --output "$target_alleles"

if [[ ! -s $all_candidates || ! -s $observed_candidates ]]; then
  "$python" "$repo/scripts/build_missense_candidate_alleles.py" \
    --dbnsfp "$dbnsfp" --all-genes --chrom "$chromosome" \
    --output "$all_candidates" --observed-alleles "$target_alleles" \
    --observed-output "$observed_candidates"
fi
[[ -s $all_candidates && -s $observed_candidates ]]
printf 'target_inputs_ready\t%s\n' "$chromosome"

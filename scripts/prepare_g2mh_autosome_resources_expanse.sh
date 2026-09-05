#!/usr/bin/env bash
# Prepare configurable G2MH autosome region tracks and upload a private bundle.
# Submit only to Expanse ind-shared.
set -euo pipefail

upload_manifest=${1:?usage: prepare_g2mh_autosome_resources_expanse.sh PRIVATE_UPLOAD_MANIFEST}
output_root=${RESOURCE_PREP_ROOT:-/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-resources/autosomes}
region_source=${REGION_SOURCE:-/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/resources/regions}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
read -r -a chromosomes <<< "${CHROMOSOMES:-chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20}"

mkdir -p "$output_root/regions"
for chromosome in "${chromosomes[@]}"; do
  for filename in genomicSuperDups.bed simpleRepeat.bed rmsk.bed; do
    output="$output_root/regions/${filename%.bed}.${chromosome}.bed"
    if ! test -s "$output"; then
      awk -v chromosome="$chromosome" '$1 == chromosome' \
        "$region_source/$filename" > "$output"
    fi
    test -s "$output"
  done
done

"$script_dir/stage_resource_bundle_expanse.sh" "$upload_manifest"

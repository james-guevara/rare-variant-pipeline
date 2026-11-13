#!/usr/bin/env bash
# Download and chromosome-subset the three canonical problematic-region tracks.
set -euo pipefail

chromosome=${1:?usage: prepare_chromosome_region_filters.sh CHROMOSOME OUTPUT_DIR}
output_dir=${2:?usage: prepare_chromosome_region_filters.sh CHROMOSOME OUTPUT_DIR}
task_tmp=$(mktemp -d /tmp/rvp-region-filters.XXXXXX)
trap 'rm -rf "$task_tmp"' EXIT
mkdir -p "$output_dir"

curl --fail --location --retry 5 --silent --show-error \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz \
  | gzip -dc | awk -v chrom="$chromosome" 'BEGIN{OFS="\t"} $2==chrom{print $2,$3,$4}' \
  > "$task_tmp/genomicSuperDups.$chromosome.bed"
curl --fail --location --retry 5 --silent --show-error \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz \
  | gzip -dc | awk -v chrom="$chromosome" 'BEGIN{OFS="\t"} $2==chrom{print $2,$3,$4,$5}' \
  > "$task_tmp/simpleRepeat.$chromosome.bed"
curl --fail --location --retry 5 --silent --show-error \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz \
  | gzip -dc | awk -v chrom="$chromosome" \
      'BEGIN{OFS="\t"} $6==chrom{print $6,$7,$8,$2,$10,$11,$12,$13}' \
  > "$task_tmp/rmsk.$chromosome.bed"

for name in genomicSuperDups simpleRepeat rmsk; do
  source="$task_tmp/$name.$chromosome.bed"
  test -s "$source"
  install -m 0644 "$source" "$output_dir/$name.$chromosome.bed"
done
sha256sum "$output_dir"/*.$chromosome.bed

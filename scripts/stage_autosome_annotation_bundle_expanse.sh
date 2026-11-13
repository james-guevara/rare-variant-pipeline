#!/usr/bin/env bash
#SBATCH --partition=ind-shared
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=rvp_stage_annotation
#SBATCH --output=/expanse/projects/sebat1/resources/rare-variant-pipeline/logs/stage-annotation-%j.out
set -euo pipefail

urls=${1:?usage: stage_autosome_annotation_bundle_expanse.sh PRIVATE_URL_FILE}
registry=/expanse/projects/sebat1/resources/rare-variant-pipeline
staging=$registry/staging
archive=$staging/ensembl115-autosomes-chr2-20.tar.zst
checksum=$archive.sha256
mkdir -p "$staging" "$registry/logs"

members=()
for contig in $(seq 2 20); do
  members+=("annotation/ensembl-115/chr$contig")
  members+=("loftee/ensembl-115/chr$contig.transcripts.sqlite")
done
tar -C "$registry" -cf - "${members[@]}" | zstd -T8 -3 -f -o "$archive"
sha256sum "$archive" > "$checksum"

mapfile -t put_urls < "$urls"
test "${#put_urls[@]}" -eq 2
curl --fail --silent --show-error --retry 5 --upload-file "$archive" "${put_urls[0]}"
curl --fail --silent --show-error --retry 5 --upload-file "$checksum" "${put_urls[1]}"
printf 'annotation_bundle_staged\tbytes=%s\tsha256=%s\n' \
  "$(stat -c %s "$archive")" "$(cut -d" " -f1 "$checksum")"

#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 CACHE_CHROM_DIR VEP_IMAGE OUTPUT_TSV CONTAINER_RUNTIME" >&2
  exit 2
fi

cache_dir=$(cd "$1" && pwd)
vep_image=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
output=$3
runtime=$4
vep_perl5lib=${VEP_PERL5LIB:-/usr/local/share/ensembl-vep-115.2-1}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$(dirname "$output")"
output_dir=$(cd "$(dirname "$output")" && pwd)
output_path=$output_dir/$(basename "$output")

mapfile -t chunks < <(
  find "$cache_dir" -maxdepth 1 -type f -name '*.gz' \
    ! -name '*_reg.gz' ! -name 'all_vars.gz' -print | sort -V
)
if [[ ${#chunks[@]} -eq 0 ]]; then
  echo "no VEP transcript cache chunks found in $cache_dir" >&2
  exit 1
fi

"$runtime" exec \
  --bind "$script_dir:/work/scripts:ro" \
  --bind "$cache_dir:$cache_dir:ro" \
  --bind "$output_dir:$output_dir" \
  "$vep_image" \
  env PERL5LIB="$vep_perl5lib${PERL5LIB:+:$PERL5LIB}" \
  perl /work/scripts/export_vep_transcript_priority.pl "${chunks[@]}" \
  > "$output_path"

rows=$(awk 'END { print NR - 1 }' "$output_path")
if [[ $rows -le 0 ]]; then
  echo "empty transcript-priority output: $output_path" >&2
  exit 1
fi
printf 'transcript_priority_ready\t%s\trows=%s\n' "$output_path" "$rows"

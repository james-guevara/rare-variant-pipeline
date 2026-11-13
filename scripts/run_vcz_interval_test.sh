#!/usr/bin/env bash
set -euo pipefail

if (( $# != 4 )); then
    echo "usage: $0 INPUT.vcf.gz REGION OUTPUT_DIR WORKERS" >&2
    exit 2
fi

input_vcf=$1
region=$2
output_dir=$3
workers=$4
venv=${VCZ_VENV:-/home/ubuntu/work/rare-variant-missense-scope/.venv-vcf-zarr}

mkdir -p "$output_dir"

interval_vcf="$output_dir/interval.vcf.gz"
icf="$output_dir/interval.icf"
vcz="$output_dir/interval.vcz"

/usr/bin/time -v -o "$output_dir/extract.time.txt" \
    bcftools view -r "$region" -Oz -o "$interval_vcf" "$input_vcf"
tabix -f -p vcf "$interval_vcf"

/usr/bin/time -v -o "$output_dir/explode.time.txt" \
    "$venv/bin/vcf2zarr" explode -f -Q -p "$workers" "$interval_vcf" "$icf"

"$venv/bin/vcf2zarr" mkschema "$icf" > "$output_dir/schema.json"

/usr/bin/time -v -o "$output_dir/encode.time.txt" \
    "$venv/bin/vcf2zarr" encode -f -Q -p "$workers" -M 24G \
    -l 1000 -w 1000 -s "$output_dir/schema.json" "$icf" "$vcz"

bcftools index -n "$interval_vcf" > "$output_dir/variant-count.txt"
"$venv/bin/vcf2zarr" inspect "$vcz" > "$output_dir/vcz-inspect.txt"

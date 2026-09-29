#!/usr/bin/env bash
set -euo pipefail
unit=$1
chromosome=$2
source_vcf=$3
checker=$4
threads=$5
output="$unit.sites.vcf.gz"

# Read the large genotyped source once. This unfiltered stream contains all
# source records with their eight site columns and existing INFO annotations.
# Subsequent validation/compression reads only this small genotype-free file.
bcftools view --no-version -G -Ov -o input.sites.vcf "$source_vcf"
awk -v prefix=input -v chromosome="$chromosome" -f "$checker" input.sites.vcf | sha256sum > input.sha256
bcftools view --no-version --threads "$threads" -Oz -o "$output" input.sites.vcf
bcftools index --csi --threads "$threads" "$output"
test -s "$output.csi"

bcftools query -l "$output" > output.samples
test ! -s output.samples
bcftools view --no-version -Ov "$output" |
    awk -v prefix=output -v chromosome="$chromosome" -f "$checker" | sha256sum > output.sha256
cmp input.records output.records
cmp input.sha256 output.sha256
cmp input.csq_header output.csq_header
bcftools index --nrecords "$output" > indexed.records
cmp input.records indexed.records

csq_present=false
if test "$(cat input.csq)" = 1; then csq_present=true; fi
cat receipt-prefix.json > "$unit.receipt.json"
printf '"input_records":%s,"output_records":%s,"output_samples":0,"csq_present":%s,"sites_sha256":"%s","status":"PASS"}\n' \
    "$(cat input.records)" "$(cat output.records)" "$csq_present" \
    "$(awk '{print $1}' input.sha256)" >> "$unit.receipt.json"

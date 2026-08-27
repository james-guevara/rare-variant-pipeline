#!/usr/bin/env bash
set -euo pipefail

if (( $# < 1 )); then
    echo "usage: $0 CHROM [CHROM ...]" >&2
    exit 2
fi

source_dir=/expanse/projects/sebat1/g2mh_data/final_download/wgs/jointcall_dragen/genomics_sample03/submission-82512/g2mh_1065_202606/by_chr
destination=s3://sebat-genomics-work/staged/rare-variant-pilot/g2mh
export AWS_SHARED_CREDENTIALS_FILE=${G2MH_AWS_CREDENTIALS_FILE:-$HOME/.aws/g2mh-transfer.credentials}
export AWS_DEFAULT_REGION=us-east-1

for chromosome in "$@"; do
    source_file="$source_dir/g2mh_1065_chr${chromosome}.vcf.gz"
    test -r "$source_file"
    aws s3 cp "$source_file" "$destination/$(basename "$source_file")" \
        --only-show-errors \
        --expected-size "$(stat -c %s "$source_file")"
done

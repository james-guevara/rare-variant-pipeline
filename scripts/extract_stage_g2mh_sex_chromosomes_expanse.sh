#!/usr/bin/env bash
#SBATCH --job-name=g2mh-sex-chroms
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/expanse/projects/sebat1/j3guevar/g2mh_sex_chromosomes/logs/%x-%j.out
#SBATCH --error=/expanse/projects/sebat1/j3guevar/g2mh_sex_chromosomes/logs/%x-%j.err

set -euo pipefail

source_vcf=/expanse/projects/sebat1/g2mh_data/final_download/wgs/jointcall_dragen/genomics_sample03/submission-82512/g2mh_1065_202606/g2mh_1065.vcf.gz
output_dir=/expanse/projects/sebat1/j3guevar/g2mh_sex_chromosomes/by_chr
destination=s3://sebat-genomics-work/staged/rare-variant-pilot/g2mh

export AWS_SHARED_CREDENTIALS_FILE=${G2MH_AWS_CREDENTIALS_FILE:-$HOME/.aws/g2mh-transfer.credentials}
export AWS_DEFAULT_REGION=us-east-1

mkdir -p "$output_dir"

for chromosome in X Y; do
    output_vcf="$output_dir/g2mh_1065_chr${chromosome}.vcf.gz"
    bcftools view --threads 15 --regions "chr${chromosome}" --output-type z \
        --output-file "$output_vcf" "$source_vcf"
    bcftools index --threads 15 --tbi --force "$output_vcf"
    bcftools index --nrecords "$output_vcf"
    bcftools view --header-only "$output_vcf" | grep -q '^#CHROM'

    aws s3 cp "$output_vcf" "$destination/$(basename "$output_vcf")" \
        --only-show-errors --expected-size "$(stat -c %s "$output_vcf")"
    aws s3 cp "$output_vcf.tbi" "$destination/$(basename "$output_vcf.tbi")" \
        --only-show-errors
done

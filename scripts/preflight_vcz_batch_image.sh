#!/usr/bin/env bash
set -euo pipefail

image=${1:?usage: preflight_vcz_batch_image.sh IMAGE [WORK_DIR]}
work_dir=${2:-/tmp/vcz-batch-preflight}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fixture=${VCZ_PREFLIGHT_FIXTURE:-$script_dir/../tests/fixtures/vcz-preflight.vcf}

rm -rf "$work_dir"
mkdir -p "$work_dir"
cp "$fixture" "$work_dir/tiny.vcf"

docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$work_dir:/work" \
    -v "$work_dir:/fsx" \
    --entrypoint /bin/bash \
    "$image" -lc '
        set -euo pipefail
        python -c "import bio2zarr, cyvcf2, numpy, zarr"
        bgzip -f /work/tiny.vcf
        bcftools index -f /work/tiny.vcf.gz
        vcf2zarr explode -f -Q -p 2 /work/tiny.vcf.gz /work/tiny.icf
        vcf2zarr mkschema /work/tiny.icf > /work/schema.json
        vcf2zarr encode -f -Q -p 2 -M 1G \
            -s /work/schema.json /work/tiny.icf /work/tiny.vcz
        python /usr/local/bin/repack_vcz_zarr_v3.py \
            /work/tiny.vcz /work/tiny.sharded-v3.zarr \
            --workers 2 --codec-threads 1 --compression-level 7 \
            --compressed-passthrough \
            --variant-shard-chunks 2 --overwrite --validate
        vcf2zarr convert -f -Q -p 2 --zarr-format=3 \
            /work/tiny.vcf.gz /work/tiny.unsharded-v3.zarr
        python /usr/local/bin/repack_vcz_zarr_v3.py \
            /work/tiny.unsharded-v3.zarr /work/tiny.v3-to-sharded-v3.zarr \
            --workers 2 --codec-threads 1 --compression-level 7 \
            --compressed-passthrough \
            --variant-shard-chunks 2 --overwrite --validate
        mkdir -p /fsx/checkpoint-source
        cp /work/tiny.vcf.gz /fsx/checkpoint-source/g2mh_1065_chr99.vcf.gz
        cp /work/tiny.vcf.gz.csi /fsx/checkpoint-source/g2mh_1065_chr99.vcf.gz.csi
        VCZ_SOURCE_ROOT=/fsx/checkpoint-source \
        VCZ_WORK_ROOT=/fsx/checkpoint-work \
        VCZ_OUTPUT_ROOT=/fsx/checkpoint-output \
        VCZ_WORKERS=2 VCZ_MEMORY=1G VCZ_VALIDATION=sampled \
            /usr/local/bin/run_chromosome.sh 99
        resume_output=$(VCZ_SOURCE_ROOT=/fsx/checkpoint-source \
            VCZ_WORK_ROOT=/fsx/checkpoint-work \
            VCZ_OUTPUT_ROOT=/fsx/checkpoint-output \
            VCZ_WORKERS=2 VCZ_MEMORY=1G VCZ_VALIDATION=sampled \
            /usr/local/bin/run_chromosome.sh 99)
        grep -q "checkpoint_resume.*stage=complete" <<<"$resume_output"
        python - <<"PY"
import zarr
store = zarr.open_group("/work/tiny.vcz", mode="r", zarr_format=2)
assert store["variant_position"][:].tolist() == [10001, 10011]
assert store["sample_id"][:].tolist() == ["SAMPLE_A", "SAMPLE_B"]
assert store["call_DP"][:].tolist() == [[21, 18], [33, 27]]
sharded = zarr.open_group("/work/tiny.sharded-v3.zarr", mode="r", zarr_format=3)
assert sharded["variant_position"][:].tolist() == [10001, 10011]
assert sharded["sample_id"][:].tolist() == ["SAMPLE_A", "SAMPLE_B"]
assert sharded["call_DP"][:].tolist() == [[21, 18], [33, 27]]
direct_v3 = zarr.open_group(
    "/work/tiny.v3-to-sharded-v3.zarr", mode="r", zarr_format=3
)
assert direct_v3["variant_position"][:].tolist() == [10001, 10011]
assert direct_v3["sample_id"][:].tolist() == ["SAMPLE_A", "SAMPLE_B"]
assert direct_v3["call_DP"][:].tolist() == [[21, 18], [33, 27]]
print("VCZ v2/v3 to sharded-v3 image preflight passed")
PY
    '

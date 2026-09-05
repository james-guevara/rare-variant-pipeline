#!/usr/bin/env bash
set -euo pipefail

chromosome=${1:?usage: run_chromosome.sh CHROMOSOME}
workers=${VCZ_WORKERS:-16}
memory=${VCZ_MEMORY:-48G}
codec_threads=${VCZ_CODEC_THREADS:-1}
compression_level=${VCZ_COMPRESSION_LEVEL:-7}
variant_shard_chunks=${VCZ_VARIANT_SHARD_CHUNKS:-25}
validation=${VCZ_VALIDATION:-sampled}
source_root=${VCZ_SOURCE_ROOT:-/fsx/rare-variant-pilot/g2mh}
work_root=${VCZ_WORK_ROOT:-/fsx/vcz-work}
output_root=${VCZ_OUTPUT_ROOT:-/fsx/vcz-output}

source_vcf=${VCZ_SOURCE_VCF:-"$source_root/g2mh_1065_chr${chromosome}.vcf.gz"}
chrom_work="$work_root/chr${chromosome}"
icf="$chrom_work/chr${chromosome}.icf"
unsharded_v3="$chrom_work/chr${chromosome}.unsharded-v3.zarr"
schema="$chrom_work/schema.json"
v3=${VCZ_OUTPUT_STORE:-"$output_root/chr${chromosome}.sharded-v3.zarr"}
metrics=${VCZ_METRICS_FILE:-"$output_root/chr${chromosome}.metrics.tsv"}
index_checkpoint="$chrom_work/index.complete"
explode_checkpoint="$chrom_work/explode.complete"
encode_checkpoint="$chrom_work/encode.complete"
shard_checkpoint="${v3}.complete"

if [[ "$source_vcf $work_root $output_root $v3" == *"/fsx/"* ]]; then
    mountpoint -q /fsx || { echo "ERROR: /fsx is not a mounted filesystem" >&2; exit 1; }
fi
test -r "$source_vcf" || { echo "ERROR: source VCF is not readable: $source_vcf" >&2; exit 1; }
for command in bcftools vcf2zarr python; do
    command -v "$command" >/dev/null || { echo "ERROR: missing command: $command" >&2; exit 1; }
done
case "$validation" in
    full|sampled) ;;
    *) echo "ERROR: VCZ_VALIDATION must be full or sampled" >&2; exit 1 ;;
esac
mkdir -p "$chrom_work" "$output_root" "$(dirname "$v3")" "$(dirname "$metrics")"
test -w "$chrom_work" -a -w "$(dirname "$v3")" || {
    echo "ERROR: work/output paths are not writable" >&2
    exit 1
}

printf 'preflight\tfsx_mounted=true\tsource=%s\tworkers=%s\tmemory=%s\tvalidation=%s\n' \
    "$source_vcf" "$workers" "$memory" "$validation"

if [ ! -f "$metrics" ]; then
    printf 'stage\telapsed_seconds\n' > "$metrics"
fi

run_stage() {
    local stage=$1
    shift
    local started=$SECONDS
    "$@"
    printf '%s\t%s\n' "$stage" "$((SECONDS - started))" >> "$metrics"
}

validate_store() {
python - "$1" "$schema" <<'PY'
import json
import sys
from pathlib import Path

import zarr

store_path, schema_path = sys.argv[1:]
schema = json.loads(Path(schema_path).read_text())
store = zarr.open_group(store_path, mode="r", zarr_format=3)
expected = schema["dimensions"]["variants"]["size"]
if store["variant_position"].shape[0] != expected:
    raise SystemExit("variant count mismatch")
if "sample_id" not in store or len(list(store.arrays())) == 0:
    raise SystemExit("incomplete VCZ store")
PY
}

if [ -f "$shard_checkpoint" ] && validate_store "$v3"; then
    printf 'checkpoint\tshard\n' >> "$metrics"
    echo "checkpoint_resume\tstage=complete"
    exit 0
fi
rm -f "$shard_checkpoint"

if [ -f "$index_checkpoint" ] && bcftools index -n "$source_vcf" >/dev/null; then
    printf 'index\tskipped_checkpoint\n' >> "$metrics"
else
    rm -f "$index_checkpoint"
    run_stage index bcftools index -f --threads "$workers" "$source_vcf"
    bcftools index -n "$source_vcf" >/dev/null
    touch "$index_checkpoint"
fi

if [ -f "$encode_checkpoint" ] && validate_store "$unsharded_v3"; then
    printf 'explode\tskipped_encode_checkpoint\n' >> "$metrics"
    printf 'encode\tskipped_checkpoint\n' >> "$metrics"
else
    rm -f "$encode_checkpoint"
    if [ -f "$explode_checkpoint" ] && [ -d "$icf" ] && [ -s "$schema" ]; then
        vcf2zarr mkschema "$icf" >/dev/null
        printf 'explode\tskipped_checkpoint\n' >> "$metrics"
    else
        rm -f "$explode_checkpoint"
        run_stage explode vcf2zarr explode -f -Q -p "$workers" "$source_vcf" "$icf"
        schema_tmp="${schema}.tmp.$$"
        vcf2zarr mkschema "$icf" > "$schema_tmp"
        mv "$schema_tmp" "$schema"
        touch "$explode_checkpoint"
    fi
    run_stage encode vcf2zarr encode -f -Q -p "$workers" -M "$memory" \
        --zarr-format=3 -s "$schema" "$icf" "$unsharded_v3"
    validate_store "$unsharded_v3"
    touch "$encode_checkpoint"
fi

run_stage shard python /usr/local/bin/repack_vcz_zarr_v3.py \
    "$unsharded_v3" "$v3" --variant-shard-chunks "$variant_shard_chunks" \
    --workers "$workers" --codec-threads "$codec_threads" \
    --compression-level "$compression_level" --compressed-passthrough --overwrite \
    $([ "$validation" = full ] && printf '%s' --validate)

python - "$source_vcf" "$schema" "$v3" >> "$metrics" <<'PY'
import json
import sys
from pathlib import Path

import zarr

source_vcf, schema_path, v3_path = sys.argv[1:]
schema = json.loads(Path(schema_path).read_text())
store = zarr.open_group(v3_path, mode="r", zarr_format=3)
expected = schema["dimensions"]["variants"]["size"]
observed = store["variant_position"].shape[0]
if observed != expected:
    raise SystemExit(f"variant count mismatch: expected {expected}, observed {observed}")
print(f"variant_count\t{observed}")
print(f"sample_count\t{store['sample_id'].shape[0]}")
print(f"array_count\t{len(list(store.arrays()))}")
print(f"source_vcf\t{source_vcf}")
PY

if [ "$validation" = "sampled" ]; then
python - "$unsharded_v3" "$v3" <<'PY'
import sys

import numpy as np
import zarr

source_path, destination_path = sys.argv[1:]
source = zarr.open_group(source_path, mode="r", zarr_format=3)
destination = zarr.open_group(destination_path, mode="r", zarr_format=3)
for name, source_array in source.arrays():
    target_array = destination[name]
    if source_array.shape != target_array.shape or source_array.dtype != target_array.dtype:
        raise SystemExit(f"metadata mismatch for {name}")
    if source_array.ndim and source_array.shape[0]:
        rows = np.unique(np.linspace(0, source_array.shape[0] - 1, min(16, source_array.shape[0]), dtype=int))
        np.testing.assert_array_equal(source_array[rows], target_array[rows], err_msg=name)
print("sampled_validation\tpassed")
PY
fi

touch "$shard_checkpoint"
rm -rf "$icf"

#!/usr/bin/env bash
set -euo pipefail
umask 002
ROOT=/expanse/projects/sebat1/resources/rare-variant-pipeline
WORK=$ROOT/maintenance/repair-20260929
RELEASE=$ROOT/releases/v1
module load cpu/0.17.3b singularitypro/4.1.2
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1
IMAGE=$ROOT/containers/v1/targeted-81e0e00a97f900867808e5f7cf17cdac6f5c5ef7956e2c6641f3188fa335d9ab.sif
singularity exec --cleanenv -B /expanse -B "$RELEASE:$RELEASE:ro" "$IMAGE" python "$WORK/validate.py" "$RELEASE" "$WORK/validation-readonly"
SMOKE_IMAGE=$ROOT/containers/v1/targeted-7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d.sif
SIF_SHA=$(cut -d ' ' -f1 "$SMOKE_IMAGE.sha256")
singularity exec --cleanenv -B /expanse -B "$RELEASE:$RELEASE:ro" "$SMOKE_IMAGE" python "$WORK/benchmark.py" --input "$WORK/validation-readonly/chr22.synthetic.vcf" --resources "$RELEASE/targeted-annotation/ensembl-115" --outdir "$WORK/benchmark7d-readonly" --sif-sha256 "$SIF_SHA"
cmp "$WORK/validation-readonly/chr22.picked.tsv" "$WORK/benchmark7d-readonly/picked.tsv"
echo 'BOTH CONTAINERS: byte-identical chr22 output with pinned caches and read-only release'
python3 "$WORK/post_verify.py"

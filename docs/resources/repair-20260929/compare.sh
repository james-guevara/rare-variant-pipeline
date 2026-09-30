#!/usr/bin/env bash
set -euo pipefail
umask 002
ROOT=/expanse/projects/sebat1/resources/rare-variant-pipeline
WORK=$ROOT/maintenance/repair-20260929
module load cpu/0.17.3b singularitypro/4.1.2
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1
SMOKE_IMAGE=$ROOT/containers/v1/targeted-7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d.sif
SIF_SHA=$(cut -d ' ' -f1 "$SMOKE_IMAGE.sha256")
singularity exec --cleanenv -B /expanse "$SMOKE_IMAGE" python "$WORK/benchmark.py" \
  --input "$WORK/validation/chr22.synthetic.vcf" \
  --resources "$ROOT/releases/v1/targeted-annotation/ensembl-115" \
  --outdir "$WORK/benchmark7d" --sif-sha256 "$SIF_SHA"
cmp "$WORK/validation/chr22.picked.tsv" "$WORK/benchmark7d/picked.tsv"
echo 'BOTH CONTAINERS: byte-identical chr22 synthetic FastVEP/Rust-picker output'

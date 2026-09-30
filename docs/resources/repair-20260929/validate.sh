#!/usr/bin/env bash
set -euo pipefail
umask 002
ROOT=/expanse/projects/sebat1/resources/rare-variant-pipeline
WORK=$ROOT/maintenance/repair-20260929
module load cpu/0.17.3b singularitypro/4.1.2
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1
IMAGE=$ROOT/containers/v1/targeted-81e0e00a97f900867808e5f7cf17cdac6f5c5ef7956e2c6641f3188fa335d9ab.sif
singularity exec --cleanenv -B /expanse "$IMAGE" python "$WORK/validate.py" "$ROOT/releases/v1" "$WORK/validation"

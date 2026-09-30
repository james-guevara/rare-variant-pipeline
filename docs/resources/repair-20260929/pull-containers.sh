#!/usr/bin/env bash
set -euo pipefail
umask 002
ROOT=/expanse/projects/sebat1/resources/rare-variant-pipeline
WORK=$ROOT/maintenance/repair-20260929
module load cpu/0.17.3b singularitypro/4.1.2
export SINGULARITY_DOCKER_USERNAME=AWS
SINGULARITY_DOCKER_PASSWORD=$(cat "$WORK/private-ecr-token")
export SINGULARITY_DOCKER_PASSWORD
export SINGULARITY_CACHEDIR="$WORK/container-cache"
export SINGULARITY_TMPDIR="$WORK/container-tmp"
mkdir -p "$ROOT/containers/v1" "$SINGULARITY_TMPDIR"
trap 'rm -f "$WORK/private-ecr-token"' EXIT
for DIGEST in 7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d 81e0e00a97f900867808e5f7cf17cdac6f5c5ef7956e2c6641f3188fa335d9ab; do
  OUTPUT="$ROOT/containers/v1/targeted-$DIGEST.sif"
  if ! test -s "$OUTPUT"; then
    singularity pull "$OUTPUT.partial" "docker://640838474376.dkr.ecr.us-east-1.amazonaws.com/rare-variant-pipeline-targeted@sha256:$DIGEST"
    mv "$OUTPUT.partial" "$OUTPUT"
  fi
  sha256sum "$OUTPUT" > "$OUTPUT.sha256"
  singularity inspect --json "$OUTPUT" > "$OUTPUT.inspect.json"
  singularity exec "$OUTPUT" fastvep annotate --help >/dev/null
  singularity exec "$OUTPUT" fastvep-picker --help >/dev/null
  echo "CONTAINER VERIFIED $DIGEST"
done

#!/usr/bin/env bash

# This file is sourced by the CodeBuild bootstrap buildspec so IMAGE_URI and
# IMAGE_DIGEST can be returned as exported CodeBuild variables.
set -euo pipefail

: "${AWS_DEFAULT_REGION:?AWS_DEFAULT_REGION is required}"
: "${ECR_REPOSITORY:?ECR_REPOSITORY is required}"
: "${GITHUB_REPOSITORY:?GITHUB_REPOSITORY is required}"
: "${GITHUB_SHA:?GITHUB_SHA is required}"

case "$GITHUB_REPOSITORY" in
  *[!A-Za-z0-9_./-]*|/*|*/|*..*)
    echo "Invalid GITHUB_REPOSITORY: $GITHUB_REPOSITORY" >&2
    return 2
    ;;
esac

if [[ ! "$GITHUB_SHA" =~ ^[0-9a-fA-F]{40}$ ]]; then
  echo "GITHUB_SHA must be a full 40-character Git commit: $GITHUB_SHA" >&2
  return 2
fi

aws sts get-caller-identity >/dev/null
docker info >/dev/null

test "$(git rev-parse HEAD)" = "$GITHUB_SHA"

aws_account_id="$(aws sts get-caller-identity --query Account --output text)"
registry="${aws_account_id}.dkr.ecr.${AWS_DEFAULT_REGION}.amazonaws.com"
image_tag="sha-${GITHUB_SHA:0:12}"
export IMAGE_URI="${registry}/${ECR_REPOSITORY}:${image_tag}"

existing_digest="$(
  aws ecr describe-images \
    --repository-name "$ECR_REPOSITORY" \
    --image-ids "imageTag=$image_tag" \
    --query 'imageDetails[0].imageDigest' \
    --output text 2>/dev/null || true
)"
if [[ -n "$existing_digest" && "$existing_digest" != "None" ]]; then
  export IMAGE_DIGEST="$existing_digest"
  printf 'Reusing %s@%s\n' "${registry}/${ECR_REPOSITORY}" "$IMAGE_DIGEST"
  return 0
fi

aws ecr get-login-password --region "$AWS_DEFAULT_REGION" \
  | docker login --username AWS --password-stdin "$registry"

docker build \
  --file docker/targeted-annotation/Dockerfile \
  --tag "$IMAGE_URI" \
  --label "org.opencontainers.image.revision=$GITHUB_SHA" \
  .

docker run --rm --entrypoint /bin/bash "$IMAGE_URI" -c '
  set -euo pipefail
  fastvep --version
  bcftools --version | head -n1
  python -c "import duckdb, pyarrow, pyBigWig, pysam, zarr"
  python scripts/run_standalone_loftee.py --help >/dev/null
  python scripts/pick_fastvep_consequences.py --help >/dev/null
'

docker push "$IMAGE_URI"
IMAGE_DIGEST="$(
  aws ecr describe-images \
    --repository-name "$ECR_REPOSITORY" \
    --image-ids "imageTag=$image_tag" \
    --query 'imageDetails[0].imageDigest' \
  --output text
)"
export IMAGE_DIGEST
test -n "$IMAGE_DIGEST"
test "$IMAGE_DIGEST" != "None"
printf 'Published %s@%s\n' "${registry}/${ECR_REPOSITORY}" "$IMAGE_DIGEST"

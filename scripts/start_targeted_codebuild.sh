#!/usr/bin/env bash
# Start the targeted container build from the exact checked-out and pushed commit.
set -euo pipefail

region=${AWS_REGION:-us-east-1}
project=${CODEBUILD_PROJECT:-rare-variant-pipeline-targeted-container}
repository=${GITHUB_REPOSITORY:-james-guevara/rare-variant-pipeline}
ecr_repository=${ECR_REPOSITORY:-rare-variant-pipeline-targeted}

git diff --quiet
git diff --cached --quiet
sha=$(git rev-parse HEAD)
upstream=$(git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}')
git merge-base --is-ancestor "$sha" "$upstream" || {
  echo "ERROR: HEAD $sha is not present in $upstream; push it first" >&2
  exit 1
}

aws codebuild start-build \
  --region "$region" \
  --project-name "$project" \
  --environment-variables-override \
    "name=ECR_REPOSITORY,value=$ecr_repository,type=PLAINTEXT" \
    "name=GITHUB_REPOSITORY,value=$repository,type=PLAINTEXT" \
    "name=GITHUB_SHA,value=$sha,type=PLAINTEXT" \
  --query 'build.{id:id,status:buildStatus}' \
  --output json

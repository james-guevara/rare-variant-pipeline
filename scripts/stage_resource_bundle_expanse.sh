#!/usr/bin/env bash
# Upload a resource bundle using a private TSV of local_path<TAB>presigned_put_url.
# Intended for an Expanse ind-shared Slurm job; URLs are never written to stdout.
set -euo pipefail

upload_manifest=${1:?usage: stage_resource_bundle_expanse.sh PRIVATE_UPLOAD_MANIFEST}
test -r "$upload_manifest"

uploaded=0
while IFS=$'\t' read -r source put_url; do
  test -n "$source" || continue
  test -n "$put_url"
  test -r "$source"
  bytes=$(stat -c %s "$source")
  digest=$(sha256sum "$source" | cut -d' ' -f1)
  printf 'upload_start\t%s\t%s\n' "$(basename "$source")" "$bytes"
  curl --fail --silent --show-error \
    --retry 5 --connect-timeout 30 \
    --upload-file "$source" "$put_url"
  printf 'upload_done\t%s\t%s\n' "$(basename "$source")" "$digest"
  uploaded=$((uploaded + 1))
done < "$upload_manifest"

printf 'bundle_succeeded\tfiles=%s\n' "$uploaded"

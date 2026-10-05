# Unfiltered dbNSFP consolidation — 2026-10-05

The unfiltered `dbNSFP/5.3.1a/parquet_expanded/` files for chr1–22/X/Y were
added to the existing canonical Expanse and S3 v1 release roots. These are
24 additional runtime files, not replacements for the historical
`parquet_expanded_mane_select` representation. No pipeline code, score
aggregation, scientific filters, containers, or original resource bytes changed.

## Inventory and verification

- [Current 137-file runtime manifest](../github-runtime-manifest.tsv)
- [Candidate-resource inventory](../candidate-resource-hashes.json)
- [Per-file source, destination, size and verification evidence](verified.json)
- [S3 catalog](../aws-file-catalog.tsv)

Each source file was SHA-256 hashed, copied to a temporary file on Expanse,
independently hashed at the destination, then atomically installed read-only.
S3 uploads supplied the full-object SHA-256 checksum: S3 validated it during
PUT and the returned checksum and content length were checked using HEAD.
S3 objects were created conditionally to avoid overwriting existing objects.
The aggregate added size is 14791424972 bytes.

Original source:
`/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/parquet_expanded/`

Canonical Expanse destination:
`/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1/dbNSFP/5.3.1a/parquet_expanded/`

Canonical S3 destination:
`s3://sebat-genomics-work/resources/rare-variant-pipeline/v1/dbNSFP/5.3.1a/parquet_expanded/`

`RESOURCE_SHA256SUMS` and `DBNSFP_SHA256SUMS` include the additions. Original
manifest snapshots and the dated consolidation receipt preserve the prior
v1 contents; the September repair's 113-file checks are historical evidence,
not a claim of a new full-bundle validation. The new runtime total is 137.
The release version marker remains unchanged because this is an additive
resource installation; pin individual resource hashes for reproducibility.

Expanse maintenance evidence is under
`/expanse/projects/sebat1/resources/rare-variant-pipeline/maintenance/consolidation-20261005/`.
The release also contains `EXPANDED_CONSOLIDATION_20261005.json`.
Resource bindings expose `dbnsfp_expanded` and `DBNSFP_EXPANDED`.

## Operational scripts

The Python scripts here record the operations used for this installation.
They contain fixed deployment paths and are audit records, not generic pipeline
entrypoints. `consolidate_expanse.py` runs on Expanse; `transfer.py` runs locally
with authenticated boto3 and sends short-lived, per-object signed PUT requests
to `upload.py` over SSH. Credentials and signed URLs are not recorded here.
`install_registry.py` installs manifests and updates the Expanse registry.
`publish_manifests.py` archives the original S3 manifests and publishes the
expanded manifests with conditional writes.

## Shared ddp195 deployment

A separate verified copy for ddp195 collaborators is installed at:
`/expanse/lustre/projects/ddp195/j3guevar/rare-variant-pipeline/`.
It includes the full release under `releases/v1/`, both SIFs under
`containers/v1/`, and local resource bindings under `deployments/ddp195-v1/`.
See [the copy receipt](ddp195-copy-receipt.json) and
[the copy verification script](verify_shared.py).

The directories and files are readable by ddp195 members, including `toedwards`.
The account's membership and the full access path were inspected; this is not
an impersonated login test. The copy was verified against all 137 runtime
checksums, auxiliary file bytes, both SIF checksums, and cache timestamp ordering.

Direct access through the canonical Sebat parent remains restricted. A
user-specific traversal ACL was attempted with the owner's authorization,
but the filesystem returned `Operation not supported`. No broader parent
permissions were granted. Direct access there requires administrator assistance.
The ddp195 quota report alone was not treated as definitive: the owner confirmed
that copying should proceed, and the copy/write results establish its outcome.

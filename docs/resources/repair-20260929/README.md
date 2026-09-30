# Expanse repair completed: 2026-09-29

The consolidated v1 resource release is now physically installed and verified on Expanse:

```text
/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1/
```

`/expanse/projects/sebat1/resources/rare-variant-pipeline/current` resolves to it. Use the explicit versioned path for reproducible runs. The release contains regular, read-only files; none depends on a symlink into a pilot directory. Its source is the documented durable bundle at `s3://sebat-genomics-work/resources/rare-variant-pipeline/v1/`.

## What was repaired

- Restored all **113 runtime files**, verified against the existing GitHub manifest by byte size and SHA-256. Each download was checked while streaming and independently reread from disk.
- Restored the eight construction-source archive objects and two transfer-receipt objects, with independent disk hash checks. The S3 checksum indexes were also retained.
- Created **24 real annotation directories**, chr1–22, chrX, and chrY, with **144 canonical file aliases** into the consolidated release. chr22 is no longer a directory symlink into a mixed pilot bundle. Each chromosome exposes its GFF3/cache plus the same genome-wide FASTA/FAI, priority table, and consequence table.
- Repaired the active GeneBayes path and all five LOFTEE support-resource paths. Exposed the consolidated genome-wide transcript SQLite database through the active LOFTEE path.
- Installed the exact two pinned container images as SIFs, with SIF SHA-256 sidecars and OCI-to-SIF identity records.
- Added resource bindings at `deployments/expanse-v1/resources.json` and `resources.env`, and a full README at the Expanse registry root.
- Preserved the former annotation layout, prior LOFTEE-root link, and prior GeneBayes link under `legacy/registry-before-20260929/`. Existing historical filenames remain available for compatibility. No cohort inputs, original pilot files, or scientific run outputs were deleted.

Canonical data and images are read-only. New deployments should use `releases/v1/`, not legacy chromosome-specific FASTA/priority filenames. The alias layout avoids making 24 copies of the genome-wide reference and tables.

## Final verification

Final Slurm job **54530330** completed with exit code **0:0**. It ran both containers with the release bind-mounted read-only, then independently rehashed all 113 runtime files.

| Check | Result |
|---|---|
| SHA-256 after final validation | **113/113 match** the pinned runtime manifest |
| Direct pinned-cache loads | **24/24**; no stale-cache fallback or cache save |
| Synthetic input | **480 sites**, 20 per chromosome |
| FastVEP + Rust picker | **480 picked rows**, all 24 chromosomes passed |
| Standalone LOFTEE | All 24 chromosomes passed; **397 classified rows** across synthetic candidate sites |
| SQLite | Transcript and LOFTEE conservation databases passed integrity checks |
| dbNSFP | All **48** Parquets have readable metadata and nonzero row counts |
| References/tracks | FASTA/FAI, ancestor, GERP, and region BED checks passed |
| Two-container chr22 comparison | Picked TSVs are byte-identical |
| Initial versus final synthetic annotation | Picked TSVs are byte-identical for all 24 chromosomes |
| Active resource aliases | No broken links in annotation, LOFTEE, GeneBayes, containers, regions, or dbNSFP trees |
| Temporary download credentials | Removed |

This validates the restored resources and synthetic annotation path. It is not an ABCD/NBDC execution, a real-cohort regression, or evidence that cohort-specific paths/settings are valid.

## A validation issue was caught and corrected

The first validation run used writable resource mounts. FastVEP's cache freshness check requires `cache_mtime > gff_mtime`. Parallel downloads left the caches for chr3, chr9, and chr22 older than their corresponding GFF3 files. FastVEP rebuilt those three caches during validation.

The organization guard detected the byte differences and stopped its first attempt before changing the annotation tree. The initial interpretation—that old Expanse caches differed from the S3 release—was incorrect: **the original caches matched the S3 manifest; the validation run had modified the restored copies**.

Correction: restore the three original pinned cache files; make all 24 cache mtimes strictly newer than their GFF3 mtimes; remove write permission from the release; rerun all tests with an explicit read-only container bind; independently hash all 113 files afterward. The final hashes match and every FastVEP run loads the pinned cache directly. All 24 synthetic picked outputs remain byte-identical.

The original [organization receipt](organization-receipt.json) is retained as an observation from that attempt. Its `legacy_byte_differences` field records temporarily rebuilt caches, not a genuine discrepancy between the original legacy files and S3. The [cache-pinning receipt](cache-pinning-receipt.json) and [post-validation hashes](post-validation-hashes.json) document the corrected final state. The initial guarded failure, job 54530141, is retained in [Slurm status](slurm-status.tsv); the successful organization retry is 54530320.

When copying resources to NBDC, verify the file hashes, enforce the cache timestamp ordering, and use read-only resource mounts. The [smoke-test handoff](../../operations/abcd-fastvep-smoke/README.md) now documents this explicitly.

## Where everything is

See the [deployed Expanse README](DEPLOYED_README.md) for the complete path map, container filenames, bindings, compatibility paths, and legacy-tree explanation. That same document is installed at:

```text
/expanse/projects/sebat1/resources/rare-variant-pipeline/README.md
```

Scientific-resource subpaths are unchanged from the S3 v1 bundle:

```text
releases/v1/
  targeted-annotation/ensembl-115/     # GFF3/cache, primary FASTA/FAI, picker tables
  targeted-annotation/GeneBayes.Supplementary_Table_1.tsv
  ensembl-115/transcripts.sqlite
  loftee-grch38/                      # ancestor/indexes, GERP, SQL
  dbNSFP/5.3.1a/parquet_scores_af/
  dbNSFP/5.3.1a/parquet_expanded_mane_select/
  problematic-regions/
  sample-qc/
  postprocess/
  source-archive/
  transfer-receipts/
  runtime-manifest.tsv
  RESOURCE_SHA256SUMS
  DBNSFP_SHA256SUMS
  RESOURCE_BUNDLE_VERSION
  RESTORE_RECEIPT.json
```

The earlier [Expanse/S3 audit](../README.md) remains a dated pre-repair record. The [complete S3 file catalog](../aws-file-catalog.tsv) supplies exact object names, sizes, and runtime hashes. The two staged images are recorded in [containers.json](containers.json); their `.sha256` files use portable basenames.

The old full Ensembl VEP cache is not part of the 113-file runtime contract and was not regenerated or restored. This FastVEP/picker/standalone-LOFTEE path uses the consolidated derived resources above and does not rerun Ensembl VEP. Archived broken links may remain in the explicitly marked legacy tree as historical evidence; the active resource trees resolve correctly.

## Evidence and operational scripts

- [Runtime restoration](restore-receipt.json)
- [Source archive restoration](archive-receipt.json)
- [Final read-only resource validation](validation-readonly-receipt.json)
- [Original pinned benchmark image receipt](benchmark-readonly-receipt.json)
- [Final post-execution hashes](post-validation-hashes.json)
- [Final active-path audit](final-audit.json)
- [Cache restoration and timestamp normalization](cache-pinning-receipt.json)
- [Slurm completion and guarded-failure history](slurm-status.tsv)

The Python and shell files in this directory preserve the operational steps. They are one-time repair records with fixed paths, not a general-purpose installer. In particular, the initial `restore.py` was followed by `pin_caches.py` and `validate-readonly.sh`; do not reproduce the initial writable-validation sequence. No signed download URLs or authentication tokens are included.

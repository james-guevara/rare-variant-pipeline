# Expanse rare-variant resource registry

Canonical local release:

```
/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1/
```

`current` points to this release. Use the versioned path for reproducible runs.
The release consists of read-only regular files copied from:

```
s3://sebat-genomics-work/resources/rare-variant-pipeline/v1/
```

All 113 runtime files were checked against the GitHub runtime manifest by size and SHA-256 during download, then independently reread and hashed on Expanse. Source archives and transfer receipts were also restored and verified. `releases/v1/RESTORE_RECEIPT.json` records the runtime results. `maintenance/repair-20260929/` contains the restoration, resource-validation, and organization receipts and the exact operational scripts.

## Resource map

Paths below are relative to `releases/v1/`:

| Resource | Path |
|---|---|
| 24 chromosome GFF3 files and FastVEP caches | `targeted-annotation/ensembl-115/` |
| Indexed GRCh38 primary-assembly FASTA | `targeted-annotation/ensembl-115/Homo_sapiens.GRCh38.dna.primary_assembly.fa` and `.fai` |
| Consolidated VEP115 transcript priority | `targeted-annotation/ensembl-115/vep115.transcript-priority.tsv` |
| VEP115 consequence ranks | `targeted-annotation/ensembl-115/vep115.consequence-ranks.tsv` |
| Consolidated LOFTEE transcript SQLite | `ensembl-115/transcripts.sqlite` |
| LOFTEE ancestor/FAI/GZI, GERP BigWig and SQL | `loftee-grch38/` |
| GeneBayes | `targeted-annotation/GeneBayes.Supplementary_Table_1.tsv` |
| 24 dbNSFP score/AF Parquets | `dbNSFP/5.3.1a/parquet_scores_af/` |
| 24 dbNSFP expanded MANE Parquets | `dbNSFP/5.3.1a/parquet_expanded_mane_select/` |
| Genome-wide problematic-region tracks | `problematic-regions/` |
| GRCh38 sex-region definitions | `sample-qc/grch38-sex-chromosome-regions.json` |
| Generic postprocess rules | `postprocess/config.json` |
| All runtime identities | `runtime-manifest.tsv`, `RESOURCE_SHA256SUMS`, `DBNSFP_SHA256SUMS` |
| Construction sources | `source-archive/` |
| Original S3 transfer evidence | `transfer-receipts/` |

The postprocess JSON uses paths relative to its directory. Consumers must resolve those paths relative to that JSON, not an arbitrary working directory.

`deployments/expanse-v1/resources.json` and `resources.env` give the corresponding absolute Expanse paths. They are resource bindings, not complete cohort/run configurations. They do not select participants, inputs, filters, scheduler settings, or an output directory.

## Containers

`containers/v1/` contains two immutable-digest SIFs, each accompanied by a SIF SHA-256 and image-inspection JSON:

- `targeted-7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d.sif`: documented FastVEP/Rust picker smoke-test image.
- `targeted-81e0e00a97f900867808e5f7cf17cdac6f5c5ef7956e2c6641f3188fa335d9ab.sif`: documented consolidated-LOFTEE resource image.

These filenames identify OCI digests, not SIF hashes. Check each `.sha256` sidecar for the actual SIF hash. Invoke the required executable explicitly; the default container command runs a broader workflow.

## Compatibility paths and history

`annotation/ensembl-115/chr1` through `chr22`, `chrX`, and `chrY` are all real directories with the same six canonical resource-file aliases. chr22 is no longer a directory symlink into a pilot bundle. Individual aliases explicitly point into `releases/v1/`, which avoids duplicating the genome-wide FASTA and tables 24 times.

Existing historical filenames are retained for compatibility. Those old chromosome FASTAs and chromosome-priority tables may point into `legacy/registry-before-20260929/`. They are not the resource entrypoints for new runs. The legacy tree retains the previous annotation layout, previous LOFTEE root link, and the previous broken GeneBayes link as rollback evidence. Broken links in this archived tree do not indicate that the active repaired aliases are broken.

`loftee/GRCh38/` now resolves the five LOFTEE support files and the consolidated `ensembl-115/transcripts.sqlite`. Historical chromosome databases remain accessible. The active GeneBayes path also resolves to the restored release. Original pilot data, cohort inputs, Zarr stores, region tracks, candidate caches, old containers, and old manifests were preserved.

## Verification scope

Final job **54530330** passed: 113/113 runtime SHA-256 values still matched after read-only container validation. All 24 FastVEP runs loaded their pinned caches directly. See `maintenance/repair-20260929/POST_VALIDATION_HASHES.json` and `validation-readonly/VALIDATION_RECEIPT.json`.

Resource validation exercises synthetic variants on all 24 chromosomes through FastVEP, Rust transcript picking, and standalone LOFTEE. It checks indexed reference access, ancestor and GERP resources, integrity of both SQLite databases, all 48 dbNSFP Parquet footers, and region BED coordinates/chromosome coverage. It is not a real ABCD or G2MH cohort analysis and does not establish cohort-level scientific results.

The two staged containers additionally produce byte-identical picked TSVs for the synthetic chr22 fixture. Receipts and Slurm exit status are retained under `maintenance/repair-20260929/`.

For a receiving site, copy the needed files from this release or the exact S3 version and verify the supplied checksums there. FastVEP/Rust picking alone needs only the six chr22 annotation inputs; it does not require cohort genotypes or the downstream LOFTEE/dbNSFP resources.

## Cache timestamp contract and corrected validation attempt

FastVEP considers a cache fresh only when its filesystem mtime is strictly greater than the GFF3 mtime. Parallel downloads initially reversed this order for chr3, chr9 and chr22; the first validation run rebuilt those caches. Their original bytes were restored from the manifest-matching legacy copies, every cache timestamp was set newer than its matching GFF3, and the entire release was made read-only. All tests were rerun with the release also bind-mounted read-only, followed by an independent full 113-file hash verification. Outputs for all 24 synthetic chromosome fixtures were byte-identical before and after restoring the pinned caches.

The initial organization receipt's `legacy_byte_differences` describes the temporarily rebuilt cache copies, not genuine differences between the original Expanse caches and S3. `CACHE_PINNING_RECEIPT.json` and `POST_VALIDATION_HASHES.json` provide the correction and final verification. Preserve or reestablish cache-newer-than-GFF timestamps when copying resources to another site; do not change cache bytes. Keep production resource mounts read-only.

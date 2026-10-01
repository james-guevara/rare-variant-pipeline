# Historical pre-repair resource audit — superseded 2026-09-29

**Historical snapshot only.** The broken links and proposed restoration described below were observed before the repair. They are not the current deployment state. Start with the [current resource guide](README.md).

**Repair report:** [Expanse v1 restoration and final verification](repair-20260929/README.md). The audit below records the **pre-repair** state; consult the repair report and deployed README for current paths.

Audit: 2026-09-29. This document covers reusable rare-variant scientific resources, container identities, source archives, and deployment paths. It does not inventory protected variant records, genotypes, run outputs, or SPARK WES jobs.

## Existing authoritative documentation

The consolidation work was already documented in the integrated repository:

- [Resource locations and cleanup history](https://github.com/james-guevara/integrated_genomics_pipeline/blob/main/docs/operations/rare-variant-resource-locations.md).
- [Exact 113-file runtime manifest: paths, roles, sizes and SHA-256](https://github.com/james-guevara/integrated_genomics_pipeline/blob/main/resources/manifests/rare-resource-files.tsv).

That document records the durable S3 mirror as verified on 2026-09-14. It explicitly classifies Expanse as **“reconcile into one versioned release”**, not as a completed mirror. Earlier claims in this handoff that Expanse was consolidated were incorrect.

The old targeted-branch `docs/portable-targeted-execution.md` describes the August portable registry. It predates the September genome-wide consolidation. Both locations contain useful files, but their layouts are not interchangeable.

## Which root to use

| Environment | Exact root | Observed status |
|---|---|---|
| Durable S3 release | `s3://sebat-genomics-work/resources/rare-variant-pipeline/v1/` | Present; all 113 runtime manifest paths present. Includes additional source archives, checksum files, and transfer receipts. |
| Historical AWS runtime | `/fsx/rare-variant-resources/v1/` | Source of the archived bundle. No FSx filesystems returned in us-east-1 at audit time; do not assume this filesystem is still mounted. |
| Expanse August registry | `/expanse/projects/sebat1/resources/rare-variant-pipeline/` | Present, mixed layout, six broken resource symlinks found. Not a consolidated v1 mirror. |
| Expanse pilot backing files | `/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-resources/` | Still backs annotation, regions, cohort inputs, and some transcript databases. |
| Older S3 staging namespace | `s3://sebat-genomics-work/rare-variant-resources/` | Contains older annotation/dbNSFP archives and workflow/config objects. Not the consolidated v1 root. |
| NBDC | User-selected destination | Not inspected by Codex. Stage from the exact S3 release and validate checksums there. |

The S3 `RESOURCE_BUNDLE_VERSION` identifies `rare-variant-resources-v1`, created 2026-09-04, GRCh38, Ensembl 115, dbNSFP 5.3.1a. The [transfer receipt](aws-fsx-copy-20260913.json) names `/fsx/rare-variant-resources/v1` as source, this exact S3 prefix as destination, and completion on 2026-09-14. Its 115 transferred files include the 113 runtime files plus checksum indexes; later source archives and receipts account for the larger current listing.

Current listing: 125 objects, 37,576,087,859 bytes. The 113 runtime paths are all present. An object listing is presence evidence, not independent content verification. [HEAD comparison](aws-head-verification-20260929.json) records comparisons with stored SHA-256 metadata; download-time SHA-256 remains the receiving site's content verification.

## Consolidated v1 layout

Every path below is relative to the durable S3 release root. The same relative layout is intended for a restored local runtime directory.

| Resource | Relative path | Purpose |
|---|---|---|
| Ensembl 115 chromosome GFF3 | `targeted-annotation/ensembl-115/Homo_sapiens.GRCh38.115.chr{1..22,X,Y}.gff3` | FastVEP transcript models; 24 files |
| FastVEP caches | Same GFF3 paths plus `.fastvep.cache` | Prebuilt annotation caches; 24 files |
| Genome-wide reference | `targeted-annotation/ensembl-115/Homo_sapiens.GRCh38.dna.primary_assembly.fa` | GRCh38 primary assembly, 3,139,742,371 bytes |
| Reference index | Same FASTA path plus `.fai` | FASTA random access |
| Genome-wide transcript priority | `targeted-annotation/ensembl-115/vep115.transcript-priority.tsv` | Rust picker priority table, 22,466,432 bytes |
| Consequence ranks | `targeted-annotation/ensembl-115/vep115.consequence-ranks.tsv` | VEP-compatible consequence ordering |
| LOFTEE transcript database | `ensembl-115/transcripts.sqlite` | Consolidated genome-wide database, 790,970,368 bytes |
| LOFTEE conservation | `loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw` | 12,617,579,354-byte conservation track |
| LOFTEE ancestor reference | `loftee-grch38/human_ancestor.fa.gz`, `.fa.gz.fai`, `.fa.gz.gzi` | Ancestral sequence and indexes |
| LOFTEE SQL | `loftee-grch38/loftee.sql` | LOFTEE support database |
| GeneBayes | `targeted-annotation/GeneBayes.Supplementary_Table_1.tsv` | Gene-level annotations |
| dbNSFP score/AF | `dbNSFP/5.3.1a/parquet_scores_af/chr{1..22,X,Y}.parquet` | 24 chromosome files |
| dbNSFP expanded MANE | `dbNSFP/5.3.1a/parquet_expanded_mane_select/chr{1..22,X,Y}.parquet` | Separate 24-file representation; not interchangeable with score/AF tables |
| Problematic regions | `problematic-regions/{genomicSuperDups,rmsk,simpleRepeat}.bed` | Three consolidated genome-wide tracks |
| Sex-chromosome regions | `sample-qc/grch38-sex-chromosome-regions.json` | GRCh38 region definitions |
| Generic postprocess configuration | `postprocess/config.json` | Shared rules; inspect environment paths before running elsewhere |
| Release identity | `RESOURCE_BUNDLE_VERSION` | Release/build/version labels |
| Runtime checksums | `RESOURCE_SHA256SUMS` | 113 runtime SHA-256 entries |
| dbNSFP checksums | `DBNSFP_SHA256SUMS` | dbNSFP verification index |
| Construction inputs | `source-archive/` | Source GTF/GFF3/selected compressed FASTAs and upstream checksum files; not normal runtime inputs |
| Transfer evidence | `transfer-receipts/` | FSx-to-S3 receipt and prior checksum index |

Braces above describe filename patterns, not literal filenames. The [complete file catalog](aws-file-catalog.tsv) lists every current S3 object, its exact URI, size, stored manifest checksum where available, and modification time. [Raw S3 inventory](aws-inventory-20260929.json) also preserves ETags; ETags must not be interpreted as SHA-256.

## Why chr22 is a symlink on Expanse

The August registry adopted already-existing pilot bundles by linking to them. The old portable documentation explicitly allowed stable canonical paths to be symlinks pending later consolidation. chr22 uses one directory symlink, chr1 uses individual file symlinks, and chr2–20 were subsequently built into regular per-chromosome directories by `prepare_autosome_annotation_expanse.sh` / `prepare_ensembl_chromosome_resources.sh`.

This explains the observed layout; it does **not** establish that consolidation was finished. There is no scientific reason for chr22 alone to require a directory symlink. The filesystem does not identify who chose each linking style; the migration history is inferred from the documented policy, construction scripts, and observed paths.

Let `R=/expanse/projects/sebat1/resources/rare-variant-pipeline` and `P=/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-resources`.

| Chromosomes | Current annotation location | Physical organization |
|---|---|---|
| chr1 | `R/annotation/ensembl-115/chr1/` | Directory of six file symlinks into `P/portable-chr1-v1/annotation/` |
| chr2–20 | `R/annotation/ensembl-115/chrN/` | Regular files: chromosome GFF3, cache, FASTA, FAI, and priority TSV |
| chr22 | `R/annotation/ensembl-115/chr22/` | Directory symlink to `P/portable-chr22-v1/annotation_root/` |
| chr21, chrX, chrY | Files inside that same chr22 backing directory | No corresponding canonical chromosome directories; filenames exist despite absent directory names |

The chr22 backing directory actually contains annotation files for chr1, chr21, chr22, chrX, and chrY. It is a mixed pilot bundle, not a chr22-only resource set.

Expanse filename patterns differ from v1:

```text
Homo_sapiens.GRCh38.115.chr22.gff3
Homo_sapiens.GRCh38.115.chr22.gff3.fastvep.cache
Homo_sapiens.GRCh38.dna.chromosome.22.fa
Homo_sapiens.GRCh38.dna.chromosome.22.fa.fai
vep115.chr22.transcript-priority.tsv
vep115.consequence-ranks.tsv
```

The v1 bundle instead uses the single primary-assembly FASTA and the single genome-wide priority TSV. Do not rename a chromosome FASTA to `primary_assembly.fa` to make filenames pass a check. For the current ABCD benchmark, use the six actual consolidated v1 objects listed in the handoff.

## Remaining Expanse resource families

| Family | Current location relative to R | Backing/status |
|---|---|---|
| LOFTEE chr2–20 transcripts | `loftee/ensembl-115/chrN.transcripts.sqlite` | Regular files |
| LOFTEE chr1/chr22 transcripts | `loftee/GRCh38/ensembl-115/` | `GRCh38` links to `P/portable-chr22-v1/loftee_root`; chr1 file links onward to `P/portable-chr1-v1/loftee/` |
| LOFTEE support tracks | `loftee/GRCh38/loftee-grch38/` | All five resource symlinks are broken; exact targets in inventory |
| GeneBayes | `genebayes/GeneBayes.Supplementary_Table_1.tsv` | Broken symlink |
| dbNSFP score/AF | `dbnsfp/5.3.1a/parquet_scores_af` | Symlink to `/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/parquet_scores_af`; target directory exists |
| Region tracks, chr1 | `regions/GRCh38/chr1/` | File links to `P/chr1-chr21/regions/` |
| Region tracks, chr22 | `regions/GRCh38/chr22/` | Directory link to `P/portable-chr22-v1/regions/` |
| Downloaded Ensembl sources | `references/ensembl-115-source/` | GFF3/GTF, chr2–20 compressed FASTAs, `SHA256SUMS` |
| GRCh38 reference directory | `references/GRCh38/` | Empty in this audit |
| Optional candidate caches | `candidate-bundles/genebayes-dbnsfp-5.3.1a-v1/chr2` through `chr20` | BEDs, candidate Parquets, per-chromosome `SHA256SUMS`; not needed for ABCD FastVEP smoke test |
| G2MH-specific inputs | `cohorts/g2mh/{candidates,qc,zarr}/` | Legacy cohort files/links; do not treat as generic ABCD inputs |
| Scientific/binding manifests | `manifests/` | Historical chr1/chr22 JSON manifests and regressions |
| Containers | `containers/` | Five symlinks to `/expanse/projects/sebat1/j3guevar/containers/` |
| Old transfer archive | `staging/ensembl115-autosomes-chr2-20.tar.zst` | 1,260,102,775 bytes plus checksum sidecar; older subset, not consolidated v1 |
| Logs | `logs/` | Operational logs, not scientific resources; contents excluded from audit |

The broken GeneBayes and five LOFTEE links point under the historical tree `/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/`. The documented VEP cache path `/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/VEP_CACHE/homo_sapiens/115_GRCh38` also does not exist. All six broken resources have corresponding objects in the S3 v1 bundle.

The [Expanse inventory](expanse-inventory-20260929.json) contains 347 metadata entries with exact paths, resolved paths, link targets, existence, sizes, and modes. It follows the relevant resource aliases to expose broken links hidden behind directory symlinks. It intentionally stops at logs, Zarr and dbNSFP partition roots; no protected records were read. Filesystem presence alone is not scientific validation, and byte equality with v1 has not been established for the old Expanse files.

## Container identity

For the FastVEP + Rust picker benchmark, use the image documented by the [Rust picker validation](https://github.com/james-guevara/rare-variant-pipeline/blob/b686ef796ed6706bcf0ab15ff3731021441f1cd1/docs/benchmarks/rust-fastvep-picker-chr22-20260904.md):

```text
640838474376.dkr.ecr.us-east-1.amazonaws.com/rare-variant-pipeline-targeted@sha256:7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d
```

ECR still contains this digest at audit time, tagged `sha-1671c5a76c36`. Documented FastVEP revision: `cb8113d7bab2db42cb06bb2b2a40c57b60ea2561`; Rust picker source: `1671c5a76c369da50c64320d7dc3c719ac2ab95a`.

The Expanse registry contains SIF aliases named `16e505b`, `199d55e`, `41de024`, `89470da`, and `d6bd22e`. These old names do not establish availability of the September Rust-picker image. No matching SIF was identified in that registry. Use a SIF built from the exact digest with its own recorded SHA-256, and invoke the binaries explicitly rather than the image's default full workflow.

## Restoring the versioned release

The existing GitHub cleanup plan calls for a versioned release on Expanse. The audited registry has not undergone that restoration. A suitable destination is `R/releases/v1/`, preserving the S3 relative layout; that destination is a proposal, not an existing verified directory.

Restoration requires copying the runtime files from S3, verifying `RESOURCE_SHA256SUMS`, running the resource preflight against the selected chromosomes, and then updating environment bindings. Preserve legacy roots while existing consumers still reference them. The existing cleanup document also requires checksum equivalence and a chr22 regression before removing old roots. No resources were moved or deleted during this audit.

For NBDC, only the six FastVEP/picker files are needed initially. See the [corrected smoke-test handoff](../operations/abcd-fastvep-smoke/README.md) for exact S3 download and checksum commands. Standalone LOFTEE and downstream resources are separate later stages.

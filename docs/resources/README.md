# Current rare-variant resource guide

**Start here for resource locations.** The consolidated v1 release was installed
and verified on Expanse on **2026-09-29**. Earlier reports of broken active LOFTEE
or GeneBayes links and an unconsolidated chr22 pilot directory describe the
**pre-repair** state. They do not describe this release.

This guide summarizes the completed repair and recorded deployment paths; it does
not claim a new live filesystem audit on every documentation update. For exact
repair evidence, see the [completed repair report](repair-20260929/README.md).

## Canonical locations

| Purpose | Path |
|---|---|
| Expanse scientific release | `/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1/` |
| Expanse convenience alias | `/expanse/projects/sebat1/resources/rare-variant-pipeline/current` → `releases/v1/` |
| Expanse pinned SIF containers | `/expanse/projects/sebat1/resources/rare-variant-pipeline/containers/v1/` |
| Expanse resource bindings | `/expanse/projects/sebat1/resources/rare-variant-pipeline/deployments/expanse-v1/resources.json` and `resources.env` |
| Expanse installed documentation | `/expanse/projects/sebat1/resources/rare-variant-pipeline/README.md` |
| Durable S3 release | `s3://sebat-genomics-work/resources/rare-variant-pipeline/v1/` |

Use the explicit `releases/v1/` path for reproducible runs. The scientific release
contains read-only regular files and does not depend on old pilot-directory links.
S3 and Expanse use the same relative scientific-resource layout. SIF containers
are in the separate Expanse directory above.

## Files within the scientific release

All paths in this table are relative to `releases/v1/` on Expanse or the S3 v1 root.

| Resource | Relative path |
|---|---|
| Ensembl 115 GFF3 and FastVEP cache, chr1–22/X/Y | `targeted-annotation/ensembl-115/Homo_sapiens.GRCh38.115.chrN.gff3` and `.gff3.fastvep.cache` |
| GRCh38 primary-assembly FASTA and index | `targeted-annotation/ensembl-115/Homo_sapiens.GRCh38.dna.primary_assembly.fa` and `.fa.fai` |
| Rust picker transcript priority | `targeted-annotation/ensembl-115/vep115.transcript-priority.tsv` |
| Rust picker consequence ranks | `targeted-annotation/ensembl-115/vep115.consequence-ranks.tsv` |
| Genome-wide LOFTEE transcript database | `ensembl-115/transcripts.sqlite` |
| LOFTEE ancestor and indexes | `loftee-grch38/human_ancestor.fa.gz`, `.fa.gz.fai`, `.fa.gz.gzi` |
| LOFTEE GERP track | `loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw` |
| LOFTEE conservation database | `loftee-grch38/loftee.sql` |
| GeneBayes | `targeted-annotation/GeneBayes.Supplementary_Table_1.tsv` |
| dbNSFP score/AF Parquets | `dbNSFP/5.3.1a/parquet_scores_af/` |
| dbNSFP expanded MANE Parquets | `dbNSFP/5.3.1a/parquet_expanded_mane_select/` |
| Problematic-region tracks | `problematic-regions/{genomicSuperDups,rmsk,simpleRepeat}.bed` |
| Sex-chromosome regions | `sample-qc/grch38-sex-chromosome-regions.json` |
| Postprocess rules | `postprocess/config.json` |
| Runtime checksums | `RESOURCE_SHA256SUMS`, `DBNSFP_SHA256SUMS` |
| Release identity | `RESOURCE_BUNDLE_VERSION` |
| Construction inputs and transfer evidence | `source-archive/`, `transfer-receipts/` |

`chrN` means chr1–22, chrX, or chrY; braces describe filename patterns.
The reference and priority table are genome-wide files shared by all chromosomes.
The two dbNSFP representations are distinct products. Resolve paths in the
postprocess JSON relative to that JSON's directory.

## Containers

Under the separate Expanse `containers/v1/` directory:

| Use | Filename |
|---|---|
| FastVEP + Rust picker | `targeted-7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d.sif` |
| Standalone LOFTEE / carrier Python environment | `targeted-81e0e00a97f900867808e5f7cf17cdac6f5c5ef7956e2c6641f3188fa335d9ab.sif` |

The filenames identify OCI digests. Actual SIF SHA-256 values are recorded in
[containers.json](repair-20260929/containers.json) and adjacent `.sha256` sidecars
on Expanse. Invoke the required executable explicitly; the default image command
runs a broader workflow.

## NBDC deployment used for ABCD

These are the operator-supplied paths used by the validated ABCD workflow:

```text
/home/ood-guevara-james/abcd-fastvep-smoke/
  resources/ensembl-115/     # FastVEP/picker inputs
  resources/loftee/          # LOFTEE transcript and support resources
  containers/               # the two pinned SIFs
```

See the [NBDC annotation runbook](../operations/abcd-fastvep-smoke/NEXTFLOW.md)
for exact filenames, SIF hashes, resource-lock generation, and read-only mounts;
see the [carrier runbook](../operations/abcd-fastvep-smoke/CARRIERS.md) for the
separate indexed genotype stage. NBDC execution and validation were performed by
the operator; Codex does not have access to that filesystem.

The FastVEP cache must remain newer than its corresponding GFF3. Preserve the
validated bytes and timestamp ordering during transfers and use read-only mounts.
The workflow does not repair timestamps automatically.

## Verification and exact inventories

Final Expanse validation job **54530330** completed successfully on 2026-09-29:
113/113 runtime files matched the pinned SHA-256 manifest after read-only validation;
all 24 chromosome caches loaded directly. These are recorded repair results, not
a claim that an uninspected future transfer is already verified.

- [Completed repair and verification report](repair-20260929/README.md)
- [Complete deployed Expanse README](repair-20260929/DEPLOYED_README.md)
- [Post-validation hashes](repair-20260929/post-validation-hashes.json)
- [Exact 113-file scientific manifest](https://github.com/james-guevara/integrated_genomics_pipeline/blob/main/resources/manifests/rare-resource-files.tsv)
- [S3 object catalog with exact paths and sizes](aws-file-catalog.tsv)

## Historical material — not current deployment guidance

The [pre-repair audit](PRE_REPAIR_AUDIT_20260929.md) and its raw inventories are
retained for traceability. The older integrated-repository consolidation report and
August pilot documentation also predate the completed Expanse repair. Statements
there about broken links, an absent consolidated Expanse release, or a proposed
restoration are superseded by the repair report and this guide.

The repaired compatibility tree has 24 real chromosome directories with file aliases
into the release; chr22 is no longer a directory symlink into its pilot bundle.
Historical files and archived links remain under `legacy/registry-before-20260929/`
for compatibility. New runs should use the versioned release above. Historical
`/fsx/rare-variant-resources/v1/` paths do not establish a currently mounted AWS
filesystem; use the durable S3 release for transfers.

## Candidate-scoring resource addition (2026-10-02)

The ABCD candidate stage requires **unfiltered**
`dbNSFP/5.3.1a/parquet_expanded/chr22.parquet`. The historical v1
`parquet_expanded_mane_select` product is not a substitute. The
[candidate inventory](candidate-resource-hashes.json) records its Expanse-verified
identity; this additional file has not been confirmed staged on NBDC or copied
into the consolidated v1 release/S3. See the
[candidate runbook](../operations/abcd-fastvep-smoke/CANDIDATES.md#resource-inspection-and-identity)
for the exact source, destination, hash, and lock migration.

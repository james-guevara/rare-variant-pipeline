# Generic readiness for the remaining autosomes

## Responsibility and scope

The repository provides chromosome-independent stages, input contracts, resource
identities, tests and command templates. The operator and their execution session
supply local manifest paths/rows, locate completed products, stage any missing
resources and launch jobs. This document does not claim NBDC filesystem access,
create an ABCD manifest by guessing block IDs, or launch jobs.

Chr22 gathering has [operator-reported validation](../../validation/abcd-chr22-post-qc-gather-nbdc.md)
at `282a0862cb479deeb284e54b0469f5e1c7e45450`. Leave completed chr22 outputs and
work directories intact. Expansion concerns chr1–21 only. Scientific auditing,
new QC policies and X/Y are deferred.

## Inspection result

No chromosome-specific scientific implementation is needed. All seven existing
entrypoints accept explicit unit IDs and chromosome labels. Multiple units per
chromosome are supported; identity comes from the manifest, not filename parsing.

One actual repository gap was found: the scoring branch's canonical candidate
inventory listed only unfiltered chr22 dbNSFP. The inventory now includes the
already-verified chr1–22/X/Y files from resource-consolidation commit
`924342cc7a48fefac1cea135954750662165e6ae` (PR #14). This imports recorded hashes;
it does not rebuild/revalidate resources or change scoring. Including X/Y in the
inventory does not authorize their execution. The same chr22 and GeneBayes hashes
are retained. No workflow, QC threshold, carrier definition, or MANE policy changes.

The annotation, candidate and pre-carrier lock builders already accept comma-
separated chromosomes. Existing runtime loaders check required selected coverage
and locked file metadata. The existing sites-only, annotation, candidate, filter,
carrier, QC and gather modules require no new chromosome wiring.

## Manifest contracts

All manifests are tab-separated, with explicit unique safe `unit_id` values.
Use the same ID throughout all stages. Chromosome aliases `21` and `chr21` are
accepted; the examples below use `chr21`. Paths are local to the execution host;
relative input paths resolve against the manifest directory.

| Entrypoint / isolated config | Manifest parameter | Required columns |
|---|---|---|
| `sites_catalog.nf` / `sites_catalog.config` | `--sites_manifest` | `unit_id chromosome vcf` (original genotype VCF) |
| `annotation.nf` / `annotation.config` | `--sites_manifest` | `unit_id chromosome vcf` (published sites-only VCF) |
| `candidates.nf` / `candidates.config` | `--candidate_manifest` | `unit_id chromosome picked loftee` |
| `pre_carrier.nf` / `pre_carrier.config` | `--filter_manifest` | `unit_id chromosome missense lof_hc sites` |
| `filtered_carriers.nf` / `filtered_carriers.config` | `--carrier_manifest` | `unit_id chromosome missense lof_hc vcf index` |
| `post_carrier_qc.nf` / `post_carrier_qc.config` | `--qc_manifest` | `unit_id chromosome carriers samples source_receipt` |
| `gather_post_qc.nf` / `gather_post_qc.config` | `--gather_manifest` | `unit_id chromosome carriers samples source_receipt` |

These are separate manifests: the `vcf` column means original genotypes for
sites preparation and carrier extraction, but sites-only input for annotation.
Similarly, the QC input carriers are raw filtered-candidate extraction outputs;
gather input carriers are the post-QC `carriers.qc.tsv.gz` files.

Example shapes (illustrative paths and unit ID, not an asserted ABCD block list):

```text
# Annotation
unit_id\tchromosome\tvcf
unit_A\tchr21\t/PUBLISHED/sites/unit_A.sites.vcf.gz

# Candidate selection
unit_id\tchromosome\tpicked\tloftee
unit_A\tchr21\t/PUBLISHED/fastvep-picker/unit_A/picked.tsv\t/PUBLISHED/loftee/unit_A/loftee.tsv

# Pre-carrier filtering
unit_id\tchromosome\tmissense\tlof_hc\tsites
unit_A\tchr21\t/PUBLISHED/candidates/unit_A/missense.parquet\t/PUBLISHED/candidates/unit_A/lof_hc.parquet\t/PUBLISHED/sites/unit_A.sites.vcf.gz

# Filtered carriers
unit_id\tchromosome\tmissense\tlof_hc\tvcf\tindex
unit_A\tchr21\t/PUBLISHED/pre-carrier/unit_A/missense.filtered.parquet\t/PUBLISHED/pre-carrier/unit_A/lof_hc.filtered.parquet\t/GENOTYPES/unit_A.vcf.gz\t/GENOTYPES/unit_A.vcf.gz.tbi

# Post-extraction QC
unit_id\tchromosome\tcarriers\tsamples\tsource_receipt
unit_A\tchr21\t/PUBLISHED/filtered-carriers/unit_A/carriers.tsv.gz\t/PUBLISHED/filtered-carriers/unit_A/samples.tsv\t/PUBLISHED/filtered-carriers/unit_A/receipt.json

# Gathering
unit_id\tchromosome\tcarriers\tsamples\tsource_receipt
unit_A\tchr21\t/PUBLISHED/post-carrier-qc/unit_A/carriers.qc.tsv.gz\t/PUBLISHED/post-carrier-qc/unit_A/samples.tsv\t/PUBLISHED/post-carrier-qc/unit_A/receipt.json
```

`\t` above denotes a real tab; omit comment lines in actual manifests.

Most entrypoints support `--select_units ID1,ID2` or `all`. **Sites catalog does
not:** supply only its intended units in its manifest. All loaders check file
existence for all supplied manifest rows before selection. Do not supply rows
pointing to future products that do not exist yet. Maintain one authoritative
source inventory, then derive stage/batch manifests from it as products complete.
A batch can contain multiple chromosomes; no separate manual command per
chromosome is required. Gathering groups the supplied units by chromosome.

## Resource coverage: known versus execution-host prerequisites

The existing canonical Expanse/S3 release inventories cover chr1–22/X/Y for:

- Ensembl 115 GFF3 and matching FastVEP caches;
- unfiltered dbNSFP `parquet_expanded` for candidate scores;
- `parquet_scores_af` for the gnomAD POPmax filter.

Shared FASTA/index, picker tables, genome-wide LOFTEE resources, GeneBayes,
problematic-region tracks and validated containers already exist. No additional
reference type is introduced. Canonical roots:

```text
/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1/
s3://sebat-genomics-work/resources/rare-variant-pipeline/v1/
```

The [candidate inventory](../../resources/candidate-resource-hashes.json) now
supports all autosomes. The full consolidated runtime inventory and verification
are in [PR #14](https://github.com/james-guevara/rare-variant-pipeline/pull/14).

The only execution-host prerequisites to establish are whether the next batch's
original VCF/index and existing stage products are available, whether its
chromosome-specific resources are staged, and whether existing locks cover them.
Their absence is **not asserted** here. Check lock coverage and file metadata;
do not repeat the chr22 benchmark or full resource validation. Reuse existing
locks that already cover the batch and unchanged resources.

If coverage is missing, stage only absent chromosome-specific resources and make
new/extended locks without replacing the chr22 locks used by completed runs.
The current lock-builder CLIs build a whole requested lock and rehash shared
resources; they are not incremental extenders. Do not invoke them automatically
for every chromosome. Batch genuinely missing coverage into a single setup
operation, or handle lock extension in the resource-staging session. New file
identities still require verification against the established canonical hashes.
Do not hand-invent hashes or modify FastVEP resource timestamps.

## Rollout and resume

1. Select one unit on the next remaining autosome (chr21 is a reasonable small
   initial batch). Start at its earliest incomplete stage, consuming completed
   durable outputs. Do not infer completion from a directory's existence alone:
   require the corresponding passed receipt and expected products.
2. Take that unit through the existing stages to gather; verify receipt counts,
   sample-universe consistency, distinct/association reconciliation and separate
   star counts. Repeat unchanged commands with `-resume` to check caching.
3. Expand to the remaining units on that chromosome, then a bounded multi-
   chromosome batch, then the remaining chr1–21 units. Reuse the same per-stage
   work directories, launch directories, output roots, code and locked identities.
   The grouping gather task recomputes when its selected block set expands.
4. Keep chr22 out of execution manifests. Preserve its existing gathering result;
   it need not be regenerated. Include only autosomes in the new batches.
5. Run one stage launch at a time under the reported 8-CPU/user QoS. Four 2-CPU
   tasks use that allowance. `queueSize=4` applies per Nextflow invocation and is
   not a global cap across simultaneous launches. Annotation's two processes
   share its invocation's scheduler limit.

Do not substitute partial batch gathers for full-chromosome results. Use distinct
`--gather_id` values for pilots and complete chromosomes and check the expected
unit set. `full_manifest_selected` proves manifest coverage, not completeness of
an independently unknown cohort block inventory.

## First execution command template

The operator fills in local paths and the explicit unit ID. If a sites catalog
already exists for the selected unit, the first compute command is annotation:

```bash
nextflow -C annotation.config run annotation.nf -profile nbdc \
  --sites_manifest "$SITES_MANIFEST" \
  --resource_lock "$ANNOTATION_LOCK" \
  --select_units "$PILOT_UNIT" \
  --annotation_cpus 2 --annotation_memory '8 GB' \
  --annotation_queue_size 4 \
  --outdir "$ANNOTATION_OUT" \
  -work-dir "$ANNOTATION_WORK" \
  -with-trace "$PILOT_TRACE" -resume
```

The existing `nbdc` profile supplies the validated deployment roots/containers;
override those parameters only if the operator's actual deployment differs, and
keep them consistent with the lock. This is a template, not a claim these shell
variables are already defined. For subsequent stages use the existing
[candidate](CANDIDATES.md), [pre-carrier](PRE_CARRIER.md),
[carrier](FILTERED_CARRIERS.md), [QC](POST_CARRIER_QC.md), and
[gather](GATHER_POST_QC.md) runbooks, replacing manifests and selections with the
batch's rows. Never carry block12-specific expected candidate counts into new units.

If sites preparation is missing, use the established local sites config/SIF and
an original-VCF manifest containing only the intended batch:

```bash
nextflow -C "sites_catalog.config,$SITES_SITE_CONFIG" run sites_catalog.nf \
  --sites_manifest "$SOURCE_BATCH_MANIFEST" --sites_publish_mode copy \
  --outdir "$SITES_OUT" -work-dir "$SITES_WORK" \
  -with-trace "$SITES_TRACE" -resume
```

This deliberately does not guess the operator's bcftools SIF or Slurm config.
Original genotypes remain in shared storage; do not make a genome-wide genotype
copy. Retain work and `.nextflow` metadata for resume. Durable publications can be
used downstream independently of work caches; publication alone does not make
an upstream task `CACHED`.

## Storage and compute planning

Recorded canonical resource sizes for **chr1–21 only** (not proof they are absent
on NBDC):

| Resource | Bytes | Decimal GB |
|---|---:|---:|
| Unfiltered expanded dbNSFP | 13,909,057,521 | 13.91 |
| dbNSFP score/AF | 3,302,769,126 | 3.30 |
| GFF3 plus FastVEP caches | 1,493,404,419 | 1.49 |
| Total possible additional chromosome-specific staging | 18,705,231,066 | 18.71 |

Subtract anything already staged. Shared resources are reused and read-only bound,
not copied for every task.

Observed chr22 sites-only output was approximately 287 MB for 23 blocks. A rough
**block-count proxy only** is `287 MB × N_remaining / 23` for additional sites
publications; if the historical 951-autosome-block inventory still applies,
928 remaining blocks gives about 11.6 GB. Block sizes/chromosome density vary, so
this is not a capacity guarantee or a verified current manifest count.

For picked TSVs, candidate/filter audits, carriers, QC and gather tables, use
actual completed chr22 output sizes as the baseline rather than inventing a
compression ratio. Initial planning multiplier is `N_remaining / 23`; improve it
using record counts for completed pilot batches. Budget at least both published
outputs and retained work products, plus sites preparation's temporary
**uncompressed sites VCF**, staging/cache overhead and operational headroom.
The 18.71 GB resource estimate is not the total run-storage estimate.

Compute requests stay at existing defaults: sites 2 CPUs/2 GB, annotation
2 CPUs/8 GB per task, other downstream stages 2 CPUs/4 GB. Four concurrent tasks
therefore reserve at most 8 CPUs and typically 8–32 GB across the active stage.
Sites preparation reads original genotypes once; exact carrier extraction makes
indexed genotype queries. I/O throughput can dominate either stage.

A defensible initial wall-time estimate is
`sum(observed per-task wall seconds for the remaining workload) / 4`, plus queue
and dependency overhead, assuming four runnable tasks. Use existing chr22 traces
scaled by records/blocks as an initial proxy and refine with the first new batch.
Those traces were not supplied here, so no measured hours/days or CPU-time total
is claimed. Configured task time limits are not runtime predictions.

## Repository checks

The inventory import preserves chr22/GeneBayes identities, has unique paths, and
covers all remaining autosomes. No scientific/runtime implementation changed;
the existing 77-test gather/QC result remains recorded rather than rerunning that
unchanged suite. Real remaining-autosome execution belongs to the operator.

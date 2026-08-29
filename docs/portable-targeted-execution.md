# Portable targeted execution

## Authoritative resource locations

Do not interpret an absent chromosome under the curated portable registry as proof
that the underlying annotation resource is unavailable. The resource layers have
different purposes:

- `/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a` is the
  authoritative Expanse dbNSFP source and contains the genome-wide score and AF
  Parquets.
- `/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/VEP_CACHE/homo_sapiens/115_GRCh38`
  is the existing Ensembl VEP 115 cache. It can be used to export the
  VEP-compatible transcript-priority table for chromosomes not yet packaged.
- `/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools`
  is the validated custom VEP 115 container used for that export.
- `/expanse/projects/sebat1/resources/rare-variant-pipeline` is the curated,
  shared portable-resource registry. A chromosome appears here only after its
  derived FASTA/GFF3, FastVEP cache, LOFTEE transcript database, transcript
  priority, targets, checksums, and manifests have been validated.
- `/expanse/projects/sebat1/resources/rare-variant-pipeline/containers/rare-variant-targeted-sha-41de024.sif`
  is the validated targeted science container currently used for portable runs.
- AWS FSx is a temporary execution/staging layer, not the authoritative resource
  registry. Missing FSx files should first be recovered or derived from Expanse,
  then checksum-staged to AWS.

As of 2026-08-28, chr1, chr21, chr22, chrX, and chrY had validated portable
annotation bundles. For chr2–20, the genome-wide source resources already existed
on Expanse; only the chromosome-specific portable derivatives remained to be
built. This distinction avoids unnecessary source downloads and prevents
AWS-only resources from becoming the accidental source of truth.

Threshold-specific dbNSFP/GeneBayes candidate bundles are optional caches, not a
required workflow layer. The default multi-cohort path should broadly extract
observed coding/splice-region sites from Zarr, run FastVEP, join the fast Parquet
annotations, and recover genotypes only after filtering. This keeps changes to
missense thresholds, gene sets, or constraint models from requiring BED rebuilds.

The targeted workflow has three independent layers. Keeping them separate is the
portability contract.

1. **Scientific manifest** (`manifests/*.json`): cohort, chromosome, reference build,
   annotation release, thresholds, logical resource identities, and checksums. It
   contains no scheduler, cloud, module, account, or filesystem-root settings.
2. **Environment binding** (`config/bindings/*.json`): maps logical resource names to
   paths visible in one deployment. A new environment supplies a new binding; it does
   not edit the scientific manifest.
3. **Executor profile** (`targeted.config`): chooses the scheduler and container
   runtime. The chromosome worker itself contains no Batch, Slurm, FSx, or Expanse
   logic.

`scripts/run_targeted_manifest.py` validates both JSON documents, checks resource
sentinels and pinned checksums, performs a run-root write/delete probe, and maps the
resolved contract into the generic `scripts/run_targeted_chromosome.sh` worker. A
completed run root is immutable after `_SUCCESS`. A binding may also name a packaged
`regression` document. When present, the runner checks every pinned count and the
canonical LoF and missense hashes after scientific execution and fails the job on
any drift. Successful validated runs contain both `_SUCCESS` and
`_REGRESSION_SUCCESS`; the latter records how many count and hash assertions passed.

Some helper resources, such as the postprocessing JSON, necessarily contain local
paths and therefore belong to the binding layer. The resolver verifies that its QC
values exactly match the scientific manifest before execution; a site cannot silently
change genotype QC by supplying a different helper config.

## Preflight without scientific processing

Run inside the targeted scientific container so command checks cover the actual
runtime:

```bash
python /opt/rvp/scripts/run_targeted_manifest.py \
  --manifest manifests/g2mh-chr22.json \
  --bindings config/bindings/site.json \
  --run-root /site/path/to/a/new-run \
  --preflight-only
```

The preflight does not annotate variants or read genotype arrays. Use a fresh run root
for execution; `--run-root` lets a submission choose one without rewriting a binding.

## Nextflow controller

The controller is exactly Nextflow 26.04.6. Use uppercase `-C`: lowercase `-c` merges
the legacy repository `nextflow.config` and can leak its Expanse-specific defaults.

```bash
nextflow -C targeted.config run targeted.nf \
  -profile slurm_singularity \
  --manifest manifests/g2mh-chr22.json \
  --bindings config/bindings/site.json \
  --run_root /site/results/g2mh/chr22/run-id \
  --targeted_container docker://registry.example/rvp@sha256:DIGEST \
  --scheduler_queue queue-name \
  --scheduler_account account-name \
  --container_module singularitypro/4.1.2
```

Validated/example adapters are:

- `local_docker`: local executor and Docker;
- `slurm_singularity`: generic Slurm and Singularity/Apptainer-compatible interface;
- `expanse`: Slurm `ind-shared` and `singularitypro/4.1.2`;
- `aws_batch`: AWS Batch executor, with queue and region supplied as parameters.

AWS Batch additionally requires an S3 Nextflow work location and an explicit mount
for any host filesystem used by the bindings. For the existing us-east-1 FSx compute
environment, the controller command is:

```bash
nextflow -C targeted.config run targeted.nf -profile aws_batch \
  -bucket-dir s3://sebat-genomics-work/nextflow/rare-variant-targeted/work \
  --manifest manifests/g2mh-chrY.json \
  --bindings config/bindings/aws-g2mh-chrY.json \
  --run_root /fsx/loftee-parity/workflows/g2mh/chrY-nextflow-run-id \
  --targeted_container \
    640838474376.dkr.ecr.us-east-1.amazonaws.com/rare-variant-pipeline-targeted@sha256:DIGEST \
  --scheduler_queue rare-variant-vcz-fsx \
  --aws_batch_volumes /fsx:/fsx
```

`-bucket-dir` is controller/task scratch and checkpoint state; scientific resources
and outputs remain on the binding's FSx paths. The launch template must already mount
FSx on the Batch host, while `--aws_batch_volumes` exposes that host mount inside the
task container. Other AWS deployments supply their own queue, work bucket, and volume
mapping rather than editing the scientific manifest.

The Expanse profile is an example, not a privileged workflow target. Expanse currently
provides SingularityPRO 4.1.2 and 3.11 modules; it does not expose an Apptainer module.

Expanse's default `java` is Java 8, which is too old for Nextflow 26.04.6. The
validated controller invocation uses the Java supplied by the existing `nf_latest`
Micromamba environment while pinning the Nextflow launcher itself:

```bash
micromamba run -n nf_latest env NXF_VER=26.04.6 "$HOME/nextflow" \
  -C targeted.config run targeted.nf -profile expanse \
  -work-dir "$launch_root/work" \
  --manifest manifests/g2mh-chr22.json \
  --bindings config/bindings/expanse-g2mh-chr22.json \
  --run_root "$run_root" \
  --targeted_container \
    /expanse/projects/sebat1/resources/rare-variant-pipeline/containers/rare-variant-targeted-sha-89470da.sif
```

Run this from a persistent shell (`tmux`, `screen`, or `nohup`) because Nextflow is
the long-lived scheduler controller. The scientific process itself is submitted to
`ind-shared`; it does not execute on the login node.

The GHCR image for commit `199d55e` was pulled as a scheduled `ind-shared` job and
converted to a 189 MB SIF in shared project space. A second `ind-shared` job validated
FastVEP, bcftools, `procps`, DuckDB/Arrow/Zarr/Pysam/PyBigWig, the embedded manifests,
bindings and regressions, and both portable runner CLIs under SingularityPRO 4.1.2.
Set `SINGULARITYENV_OPENBLAS_NUM_THREADS=1` and
`SINGULARITYENV_OMP_NUM_THREADS=1` for scientific runs.

### Shared Expanse resources

Canonical pipeline paths live under:

```text
/expanse/projects/sebat1/resources/rare-variant-pipeline/
```

The hierarchy separates reference builds, Ensembl releases, LOFTEE, GeneBayes,
dbNSFP, problematic regions, cohort Zarr/candidates/QC, containers, and manifests.
Directories are owned by `jsebat-group` with setgid mode `2775`. Stable paths may be
symlinks to validated existing data initially; consumers use only the canonical path,
so backing files can later be consolidated without editing scientific manifests.

G2MH chr22 was validated through this hierarchy on both AWS and Expanse. The current
container produced identical final LoF and missense rows across systems, passing 35
count assertions and 10 canonical hashes. Its combined contract is
`resources/g2mh-chr22-targeted-regression.json`.

A clean Nextflow-controlled run from commit and container `89470da` passed the same
35 count and 10 hash assertions on 2026-08-28. Its Slurm scientific task used 4 CPUs,
peaked at about 250 MB RSS, and ran for 1 minute 25 seconds. The validated result is
under
`/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-runs/g2mh/chr22-nextflow-clean-89470da-v2`;
the corresponding Nextflow trace, report, and timeline are under
`/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-nextflow/g2mh-chr22-clean-89470da-v2`.

G2MH chr1 was subsequently validated with revision `16e505b`, passing 21 count and
10 canonical hash assertions exactly on both AWS FSx and Expanse. The Expanse
scientific task used 4 CPUs and 16 GB, peaked at about 7.8 GB RSS, and completed in
12 minutes 6 seconds. Its main stage timings were 474 seconds for targeted Zarr
extraction, 1 second for normalization, 31 seconds for FastVEP, and 138 seconds for
standalone LOFTEE. The equivalent current-container AWS run took 367, 1, 23, and
32 seconds for those stages. Extra CPUs alone will not materially improve extraction
until that currently single-threaded stage is parallelized.

The validated chr1 result is
`/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-runs/g2mh/chr1-nextflow-clean-d6bd22e`.
Its final controller receipt and reports are under
`/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-nextflow/g2mh-chr1-final-16e505b`.

Commit and container `41de024` subsequently validated the selected-row extraction
optimization from a fresh chr1 run root. It passed the same 21 count and 10 canonical
hash assertions. The scientific task completed in 6 minutes 59 seconds with 5.4 GB
peak RSS; targeted extraction took 158 seconds, compared with 474 seconds in the
pre-optimization Expanse run. The result and Nextflow controller records are under:

```text
/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-runs/g2mh/chr1-nextflow-41de024
/expanse/projects/sebat1/j3guevar/rare-variant-pipeline-targeted-nextflow/g2mh-chr1-41de024
```

The matching Expanse SIF has SHA-256
`c224347845303561fa535c348939573b26c38f8880a3708c8cef63a8e148effd`.
Use `scripts/pull_targeted_container_expanse.sh` to stage future immutable targeted
images on `ind-shared`; it skips an existing SIF and repeats the complete smoke test.

The same commit was published to ECR as digest
`sha256:fc501151b89a7623649a76a886b5bb5548ff3ffd91f764bf182e90687d2400c3`
and validated against chrX through Batch job definition revision 6. Job
`8ac1c425-9cec-4365-ac76-4bdc15b13347` passed all 35 counts and 10 hashes in 72
seconds of container runtime, including the sex-chromosome karyotype/PAR policy. Its
fresh FSx result is `/fsx/loftee-parity/workflows/g2mh/container-41de024-chrX`.

## Adding another environment

1. Make the immutable scientific container accessible to the site's runtime.
2. Stage or mount the resources and create one binding JSON mapping the existing
   logical identities to visible paths.
3. Add a profile containing only executor/runtime settings, or reuse a generic
   profile with command-line parameters.
4. Run preflight in the container.
5. Run a validated small-chromosome manifest and compare its regression outputs.

The production image embeds the versioned manifests, bindings, and regression JSONs
under `/opt/rvp`; only large biological and cohort resources need site-specific
mounts. For a local Docker preflight, bind those resources at the paths in a local
binding and run the same `run_targeted_manifest.py --preflight-only` command. This
tests the same resolver and container runtime used by Batch without requiring AWS.

Do not add site paths to the worker, duplicate the Python environment on the host, or
change scientific thresholds in an executor profile. PBS, LSF, Kubernetes, Google
Batch, another Slurm cluster, or a workstation should differ only at layers 2 and 3.

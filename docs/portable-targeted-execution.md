# Portable targeted execution

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

The Expanse profile is an example, not a privileged workflow target. Expanse currently
provides SingularityPRO 4.1.2 and 3.11 modules; it does not expose an Apptainer module.

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

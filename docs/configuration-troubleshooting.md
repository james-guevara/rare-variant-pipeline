# Configuration and infrastructure troubleshooting

This is a running ledger of operational problems encountered while validating the
targeted AWS workflow. Add an entry when a failure teaches a reusable lesson; do not
record transient scientific-result differences here.

### Resource helper expected `curl` inside science image (2026-08-28)

- **Symptom:** the Ensembl resource helper exited immediately with
  `curl: command not found` in the validated targeted container.
- **Cause:** the helper can download missing sources, but the deliberately minimal
  runtime image does not contain transfer utilities.
- **Impact:** a three-second chr2 preflight failed before producing annotation
  derivatives; no scientific chromosome run was affected.
- **Resolution:** download immutable source archives on the Expanse host, checksum
  them in the shared registry, and perform only deterministic derivation inside the
  container. Do not expand the science image solely to add `curl`.

## Prevention built into the workflow

- Build and runtime dependencies live in one pinned container.
- `uv sync --frozen` prevents unreviewed Python dependency resolution.
- Batch run roots are unique and immutable after `_SUCCESS` is written.
- The runner checks output write/delete access and all mounted inputs before work.
- Stage checkpoints allow a failed run to resume in its own incomplete run root.
- Container builds use an exact Git commit and smoke-test required commands.
- Terraform broad plans are reviewed for drift; unrelated changes are not applied.

## Incident ledger

### CodeBuild could not fetch the requested commit (2026-08-27)

- **Symptom:** build failed in two seconds at `git fetch --depth 1 origin
  $GITHUB_SHA` with exit 128.
- **Cause:** a full SHA was manually expanded from the abbreviated commit and was
  mistyped. The nonexistent object was passed to CodeBuild.
- **Impact:** no container compilation and effectively no compute expense beyond
  builder startup.
- **Avoidable:** yes.
- **Prevention:** obtain the value directly from `git rev-parse HEAD`; never type or
  reconstruct a full SHA. `scripts/start_targeted_codebuild.sh` implements this.

### Persistent EC2 user could not extend a completed Batch run (2026-08-27)

- **Symptom:** `tee` and DuckDB reported permission denied under a prior FSx run root.
- **Cause:** the directory and files were created by the Batch container user. A later
  exploratory command ran as `ubuntu` and attempted to append a new branch.
- **Impact:** failure occurred before annotation or genotype work; no material compute
  was wasted.
- **Avoidable:** yes. The deeper problem was modifying a completed run.
- **Prevention:** every submission uses a fresh `RUN_ROOT`; completed roots are marked
  with `_SUCCESS` and treated as immutable. A write/delete probe now fails immediately
  with an explicit message.

### DuckDB parameters were used in unsupported statement positions (2026-08-27)

- **Symptom:** the new missense-selector unit test failed while creating parameterized
  views; an earlier combined `COPY` statement also bound paths counterintuitively.
- **Cause:** DuckDB does not permit prepared parameters in `CREATE VIEW`, and parameters
  spanning both input functions and a `COPY` target are easy to bind incorrectly.
- **Impact:** local unit-test time only; no AWS job or biological processing ran.
- **Avoidable:** mostly.
- **Prevention:** construct input relations with DuckDB's Python API, give them named
  temporary views, and parameterize only the output `COPY`. The regression test covers
  this path.

### Batch scaled to zero between short submissions (2026-08-27)

- **Symptom:** a 36-second chr22 workflow waited roughly 4.7 minutes before starting.
- **Cause:** the compute environment correctly terminated idle EC2 workers because its
  minimum is zero. Each isolated submission therefore incurred EC2/ECS/FSx cold start.
- **Impact:** elapsed time, not processing time; idle cost is minimized.
- **Avoidable:** it is a policy tradeoff rather than an error.
- **Prevention:** submit dependent validation jobs together or temporarily retain warm
  capacity during an active validation session. Keep scale-to-zero for intermittent
  production use.

### Terraform state did not match manually changed Batch resources (2026-08-27)

- **Symptom:** a broad plan proposed changes to existing VCZ job definitions unrelated
  to the targeted workflow.
- **Cause:** active AWS revisions and a 32-vCPU definition had been created or changed
  outside the current Terraform state.
- **Impact:** none; the broad plan was not applied.
- **Avoidable:** yes, with a single infrastructure change path.
- **Prevention:** review the full plan, apply only an independently reviewed targeted
  resource when necessary, then reconcile/import drift separately. Do not accept
  unrelated replacement merely to deploy a new scientific job definition.

### Stale Terraform state lock (2026-08-27)

- **Symptom:** Terraform reported that the remote state was locked by an old process.
- **Cause:** a lock from an earlier August run remained although its recorded process
  no longer existed.
- **Impact:** diagnostic delay only; no infrastructure changed.
- **Avoidable:** not always, but easy to handle safely.
- **Prevention:** verify the recorded host/process and confirm no active Terraform run
  before force-unlocking. Never force-unlock solely because a lock is inconvenient.

### Polling variable collided with a zsh reserved parameter (2026-08-27)

- **Symptom:** a local status loop stopped with `read-only variable: status`.
- **Cause:** `status` is a reserved parameter in zsh, although the same variable name
  is ordinary in bash.
- **Impact:** the polling command alone stopped; the Batch job continued unaffected.
- **Avoidable:** yes.
- **Prevention:** use task-specific variable names such as `job_state` and specify
  bash explicitly for bash-oriented operational snippets.

### Saved AF configuration disagreed with validated output (2026-08-27)

- **Symptom:** applying the JSON's `0.001` population-AF setting retained 945 rows,
  while the saved prototype contained 1,151.
- **Cause:** the exploratory JSON was stale. The prototype output proves the executed
  rule was `gnomAD4.1_joint_AF < 0.01`, with nonnumeric `"."` values treated as
  missing. A first implementation also recognized SQL null but not `"."`.
- **Impact:** local/FSx comparison time only; no incorrect output was published and no
  Batch job ran.
- **Avoidable:** yes.
- **Prevention:** scientific thresholds are explicit runtime variables and regression
  metadata now pins expected counts. Missingness is tested via `TRY_CAST(...) IS NULL`.
  Prefer validated output evidence over an unexecuted exploratory config file.

### Historical chr22 fixture disagreed with the current container (2026-08-28)

- **Symptom:** the current AWS and Expanse runs both retained 5,291 missense rows
  after genotype QC and 1,111 after cohort eligibility, while the old split fixture
  expected 5,287 and 1,107.
- **Cause:** the split chr22 missense fixture described an earlier downstream
  container. It was not an AWS-versus-Expanse difference.
- **Evidence:** fresh current-container runs matched at every stage and across five
  LoF plus five missense canonical hashes.
- **Prevention:** use the combined `g2mh-chr22-targeted-regression.json`; the older
  split fixture is explicitly marked superseded.

### Historical chr1 fixture used incomplete aliases and stale QC counts (2026-08-28)

- **Symptom:** the complete validator initially reported several legacy chr1 fields
  as null, then reported seven missense QC/final-count differences after alias repair.
- **Cause:** the chr1 fixture predated the complete LoF-plus-missense validator and
  also described an earlier genotype-QC implementation.
- **Evidence:** fresh current-container runs on AWS FSx and Expanse matched exactly:
  8,966 genotype-QC rows and 4,533 final missense carrier rows. Five canonical LoF
  and five canonical missense hashes also matched across systems.
- **Prevention:** the validator maps every legacy chr1 key explicitly, the AWS chr1
  binding now requires regression, and the fixture pins current counts plus ten
  canonical hashes.

### Nextflow was pinned but inherited Java 8 on Expanse (2026-08-28)

- **Symptom:** Nextflow reported that Java 17 or later was required although the
  pinned 26.04.6 framework JAR was already installed.
- **Cause:** `NXF_VER` selects Nextflow, but the login shell still resolved Expanse's
  `/usr/bin/java` (Java 8).
- **Impact:** controller preflight only; no Slurm task was submitted.
- **Avoidable:** yes.
- **Prevention:** invoke the pinned launcher with a supported JVM and verify its
  version before submission:
  `micromamba run -n nf_latest env NXF_VER=26.04.6 "$HOME/nextflow" -version`.

### Binding and SIF came from different revisions (2026-08-28)

- **Symptom:** resource preflight passed, then execution stopped because
  `/opt/rvp/resources/g2mh-chr22-targeted-regression.json` was absent.
- **Cause:** the binding came from commit `89470da`, while the selected SIF came from
  the older `199d55e` commit. Small versioned contracts are embedded in the image.
- **Impact:** two one-second `ind-shared` tasks failed before scientific processing;
  Nextflow's configured retry accounted for the second task.
- **Avoidable:** yes.
- **Prevention:** pin the manifest, binding, and immutable image to the same source
  revision. Do not bypass regression or substitute an older fixture. The clean chr22
  validation used revision `89470da` for all three.

### SSM RunShellScript used `sh`, not Bash (2026-08-28)

- **Symptom:** a chr1 archive command failed immediately with
  `set: Illegal option -o pipefail`.
- **Cause:** AWS Systems Manager's `AWS-RunShellScript` executed the command block
  with Ubuntu's POSIX `sh`; Bash-specific shell options were used without explicitly
  starting Bash.
- **Impact:** controller time only; compression and upload had not started.
- **Avoidable:** yes.
- **Prevention:** either keep SSM command blocks POSIX-compatible (`set -eu`) or
  explicitly invoke a versioned Bash script. Validate the shell before beginning a
  long FSx or S3 operation.

### Portable Nextflow config inherited legacy site defaults (2026-08-27)

- **Symptom:** `nextflow config` for the local profile still showed Expanse Slurm,
  Singularity, and module settings from the repository's legacy configuration.
- **Cause:** lowercase `-c targeted.config` adds a config but still loads the default
  `nextflow.config`.
- **Impact:** configuration inspection only; no task was submitted.
- **Avoidable:** yes.
- **Prevention:** launch the portable controller with uppercase
  `nextflow -C targeted.config`, which replaces the default config set. Documentation
  and validation commands use only `-C`.

### Nextflow manifest channel had the wrong input shape (2026-08-27)

- **Symptom:** Nextflow 26.04.6 preview reported that `TARGETED_CHROMOSOME` expected two
  inputs but received one.
- **Cause:** a channel emitted a two-element tuple while the process declared two
  independent input channels.
- **Impact:** parser/preview time only; no executor or container started.
- **Avoidable:** yes.
- **Prevention:** declare a single `tuple path(...), path(...)` process input and run
  `nextflow ... -preview` as a controller compilation check.

### Generic worker lost its executable mode during rename (2026-08-27)

- **Symptom:** the manifest preflight passed, then `execvpe` returned permission denied
  for the generic worker.
- **Cause:** the file move preserved its content but not the executable Git mode.
- **Impact:** no-op runtime smoke only; `_SUCCESS` prevented scientific processing.
- **Avoidable:** yes.
- **Prevention:** preserve mode `100755`, smoke-test the packaged entry point, and let
  CodeBuild test the command from the produced image.

### Local Docker daemon was unavailable for controller validation (2026-08-27)

- **Symptom:** `docker run nextflow/nextflow:26.04.6` could not connect to the local
  daemon.
- **Cause:** Docker Desktop/daemon was not running.
- **Impact:** none; the host already had exact Nextflow 26.04.6 installed.
- **Avoidable:** not relevant to scientific execution.
- **Prevention:** validate controller syntax with the exact local Nextflow binary; do
  not create another environment or start remote compute solely for a parser check.

### Expanse curl lacked a newer retry option (2026-08-27)

- **Symptom:** the `ind-shared` resource-staging job failed before its first upload
  because `curl` did not recognize `--retry-all-errors`.
- **Cause:** that option was added after the curl release installed on Expanse.
- **Impact:** no object was uploaded and no scientific processing ran; the job failed
  immediately after its region-subset preparation.
- **Avoidable:** yes.
- **Prevention:** staging uses the older portable `--retry 5` interface. Preflight
  future transfer scripts with the destination site's actual command versions.

### FSx import task used an abbreviated operation name (2026-08-27)

- **Symptom:** AWS rejected a metadata import request with `Request failed validation`.
- **Cause:** the FSx API enum is `IMPORT_METADATA_FROM_REPOSITORY`, not the informal
  abbreviation `IMPORT_METADATA`.
- **Impact:** none; validation failed before a task was created.
- **Avoidable:** yes.
- **Prevention:** copy operation enums from `aws fsx create-data-repository-task help`
  or a prior successful task instead of abbreviating them in operational commands.

### Expanse required explicit Slurm task and node counts (2026-08-27)

- **Symptom:** the Expanse bank-limit plugin rejected a small `ind-shared` job first
  for a missing task count and then for a missing node count.
- **Cause:** Expanse requires both `--nodes` and `--ntasks` even when CPU and memory
  requests otherwise imply a single-node, single-task job.
- **Impact:** none; both submissions were rejected before scheduling.
- **Avoidable:** yes.
- **Prevention:** all Expanse submission templates explicitly specify `--nodes=1`
  and `--ntasks=1`, in addition to CPUs and memory.

### Bare Perl bypassed the VEP image's module directory (2026-08-27)

- **Symptom:** the transcript-priority exporter could not load
  `Bio::EnsEMBL::Attribute` inside the VEP 115 image.
- **Cause:** the image exposes `/usr/local/bin/perl`, but helper scripts do not inherit
  `/usr/local/share/ensembl-vep-115.2-1` in `PERL5LIB` automatically.
- **Impact:** a short metadata job failed before producing output; no chromosome
  workflow ran.
- **Avoidable:** yes.
- **Prevention:** the exporter declares the packaged VEP module directory explicitly,
  with `VEP_PERL5LIB` available for images using a different installation path.

### Targeted Terraform refresh did not converge promptly (2026-08-27)

- **Symptom:** a one-resource targeted plan remained in provider refresh for several
  minutes without producing a plan.
- **Cause:** the AWS provider process was active, but the existing pilot state has
  enough drift and remote dependencies that even a target-scoped refresh was slow.
- **Impact:** local waiting only; the plan was read-only and was cancelled before any
  infrastructure mutation.
- **Avoidable:** partly.
- **Prevention:** keep the generic job definition in Terraform, but retain an exact
  checked-in AWS registration document for time-sensitive validation. Reconcile the
  API-created revision into Terraform state separately from scientific execution.

### AWS Batch rejected manual scale-down (2026-08-27)

- **Symptom:** setting an idle managed compute environment's `desiredvCpus` from 16
  back to zero returned `Manually scaling down compute environment is not supported`.
- **Cause:** AWS Batch permits a desired-capacity increase but owns downward scaling.
- **Impact:** none beyond the normal idle scale-down interval; `minvCpus` remained zero.
- **Avoidable:** yes.
- **Prevention:** use a temporary desired-capacity increase only to accelerate cold
  starts, then allow Batch to return to `minvCpus=0` automatically. Do not disconnect
  a healthy queue merely to force immediate scale-down.

## Entry template

For future incidents record: date, symptom, root cause, impact/cost, whether it was
avoidable, the durable prevention, and the commit or test that enforces the fix.

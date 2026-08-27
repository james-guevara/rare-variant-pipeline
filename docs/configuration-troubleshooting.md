# Configuration and infrastructure troubleshooting

This is a running ledger of operational problems encountered while validating the
targeted AWS workflow. Add an entry when a failure teaches a reusable lesson; do not
record transient scientific-result differences here.

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

## Entry template

For future incidents record: date, symptom, root cause, impact/cost, whether it was
avoidable, the durable prevention, and the commit or test that enforces the fix.

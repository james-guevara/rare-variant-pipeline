# Expanse node-local SIF execution

Source `setup-local-container.sh IMAGE SHA256` from a process `beforeScript`, after loading Singularity. It requires Expanse's Slurm job-local `/scratch/$USER/job_$SLURM_JOB_ID` directory.

The helper copies the selected SIF, checks its pinned SHA-256, and prepends a small executable wrapper that substitutes only that exact image argument. It does not change the image identity recorded by Nextflow or the resource lock. Nextflow preserves PATH when launching its otherwise clean container environment. A gzip/zlib import checks startup; `container-localization.txt` records the source, digest, local path, and node in each task work directory. Slurm manages the job-local lifetime.

This mitigates shared-filesystem reads of container internals. It does not establish the underlying cause of transient I/O failures or test scientific results. Any copy, checksum, or import failure stops the task.

On 2026-10-11, SPARK chr7 post-rarity QC failed during `import gzip` with an I/O error reading `libz.so.1`. Resuming its explicit Nextflow session with this localization retained 28 cached tasks and completed the remaining 36; the chromosome inventory passed all 64 units. QC requests remained 1 CPU and 1 GB.

For rebuilt FastVEP images, the annotation resource lock can specify `containers.fastvep.build_manifest`, an absolute path inside the image. The manifest requires schema 1, a 40-hex source commit, and the executable's SHA-256. The adapter verifies that hash and records the build; it no longer attributes a rebuilt binary to the old documented OCI digest. The resource lock still pins the complete SIF.

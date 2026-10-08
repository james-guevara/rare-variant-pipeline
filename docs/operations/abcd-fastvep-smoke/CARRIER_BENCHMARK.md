# Standalone carrier backend benchmark

This experiment compares the unchanged filtered-carrier extractor with a cyvcf2
reader adapter. It does not change production containers, Nextflow entrypoints,
QC, PSAM handling, frequencies, or candidate definitions. Run only block12.
Existing annotation, extraction and QC outputs are not inputs to be overwritten.

## Pinned image and validation

- Code/image revision: `48421f965465da7ca54afb12a299fc2de6dca66d`
- OCI: `ghcr.io/james-guevara/rare-carrier-benchmark@sha256:f6d9253e1d87bcbcc0047f7bef0828cd5d77a18997e77ab87b964bab86c61a28`
- Python 3.12.12; pysam 0.23.3; cyvcf2 0.31.4; NumPy 2.5.2; DuckDB 1.5.5.
- Base image and all Python dependencies are checksum pinned in
  `benchmarks/carrier/Dockerfile` and `requirements.lock`.
- [Linux build, six synthetic tests, and publication](https://github.com/james-guevara/rare-variant-pipeline/actions/runs/37740966154)
  passed. The six tests also passed locally on macOS/Python 3.12.8.
- No real NBDC extraction or performance measurement has been run by Codex.
  No SIF checksum is claimed in advance: record the converted SIF's checksum
  and verify that same file after transfer.

Neither existing targeted SIF contains this benchmark. Use the separate image.

## Stage the container

On a Linux host with Apptainer and access to GHCR (NBDC itself if permitted):

```bash
set -euo pipefail
OCI='ghcr.io/james-guevara/rare-carrier-benchmark@sha256:f6d9253e1d87bcbcc0047f7bef0828cd5d77a18997e77ab87b964bab86c61a28'
mkdir -p carrier-benchmark-staging
cd carrier-benchmark-staging
apptainer pull --disable-cache carrier-benchmark-48421f9.sif "docker://$OCI"
sha256sum carrier-benchmark-48421f9.sif > carrier-benchmark-48421f9.sif.sha256
python3 - "$OCI" <<'PY'
import json, pathlib, sys
p = pathlib.Path('carrier-benchmark-48421f9.sif')
sha = pathlib.Path(str(p) + '.sha256').read_text().split()[0]
pathlib.Path('carrier-benchmark-container.json').write_text(json.dumps({
    'oci': sys.argv[1], 'source_commit': '48421f965465da7ca54afb12a299fc2de6dca66d',
    'sif': p.name, 'sif_sha256': sha, 'sif_bytes': p.stat().st_size
}, indent=2) + '\n')
PY
```

If the package requires authentication, first run
`apptainer registry login --username YOUR_GITHUB_USERNAME docker://ghcr.io`
and enter an authorized read-packages token at its password prompt. Package
visibility/access has not been assumed or changed. Do not paste tokens into logs.

Transfer the SIF, `.sha256`, and container JSON together using your established
NBDC transfer mechanism to:
`/home/ood-guevara-james/abcd-fastvep-smoke/containers/`.
The NBDC SSH/transfer hostname is deliberately not guessed. If pulling on NBDC,
copy those three files from the staging directory to this destination locally.
Then verify the transferred bytes on NBDC:

```bash
cd /home/ood-guevara-james/abcd-fastvep-smoke/containers
sha256sum -c carrier-benchmark-48421f9.sif.sha256
```

## Submit block12 on NBDC

No checkout is needed to execute: benchmark code and the unchanged production
extractor are baked into the pinned image. This job runs eight extractions in
alternating paired order: pysam/cyvcf2, cyvcf2/pysam, repeated. All use exactly the
same inputs and one allocation. It does not submit production Nextflow tasks.

Run the following on NBDC with Apptainer available in the batch environment:

```bash
set -euo pipefail
umask 077
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
RUN_ROOT="$BASE/carrier-benchmark-block12/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$RUN_ROOT"
cp "$BASE/containers/carrier-benchmark-container.json" "$RUN_ROOT/container.json"
printf '%s\n' '{"unit_id":"chr22_block12","chromosome":"chr22"}' > "$RUN_ROOT/metadata.json"
cat > "$RUN_ROOT/run.sh" <<'SH'
#!/bin/bash
#SBATCH --job-name=carrier-bench-block12
#SBATCH --partition=medium
#SBATCH --qos=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail
umask 077
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
SIF="$BASE/containers/carrier-benchmark-48421f9.sif"
SOURCE=/shared/release/abcd/abcd/concatenated/genetics/sequencing/snv_indel/population_vcf
CAND="$BASE/pre-carrier-nextflow/pre-carrier/chr22_block12"
apptainer exec --cleanenv \
  --bind "$BASE/pre-carrier-nextflow:$BASE/pre-carrier-nextflow:ro" \
  --bind "$SOURCE:$SOURCE:ro" \
  --bind "$RUN_ROOT:$RUN_ROOT:rw" \
  "$SIF" python /opt/rvp/benchmarks/carrier/benchmark.py \
  --metadata "$RUN_ROOT/metadata.json" \
  --missense "$CAND/missense.filtered.parquet" \
  --lof-hc "$CAND/lof_hc.filtered.parquet" \
  --vcf "$SOURCE/abcd_cohort_chr22_block12.vcf.gz" \
  --index "$SOURCE/abcd_cohort_chr22_block12.vcf.gz.tbi" \
  --expected-missense 296 --expected-hc 196 \
  --container-receipt "$RUN_ROOT/container.json" \
  --pairs 4 --outdir "$RUN_ROOT/results"
SH
sbatch --export=ALL,RUN_ROOT="$RUN_ROOT" \
  --output="$RUN_ROOT/slurm-%j.log" "$RUN_ROOT/run.sh"
printf 'Benchmark directory: %s\n' "$RUN_ROOT"
```

Use the site's required account/module initialization if your established Slurm
setup requires it. Two CPUs/8 GB/two hours are an initial allocation, not a
measured runtime estimate. The two engines execute sequentially with one reader
thread; numerical thread counts are fixed to one for both.

Expected candidate annotations are 296 missense (sequence) and 196 HC LoF
(181 sequence + 15 spanning deletions). These are not expected carrier counts.
The driver fails if candidate totals differ. Inspect the receipt's category
breakdowns to confirm the allele-class split and 8,877 source samples.

## Agreement and timing checks

Each run retains its carrier table, complete candidate/sample/burden tables,
receipt, and private logs. `results/benchmark.json` contains aggregate agreement,
timings, dependency versions, input/code identities, and supplied container
identity. Candidate Parquets, metadata, and index are SHA-256 checked before and
after; the large source VCF uses resolved path, size, and modification time to
avoid adding a whole-file scan. Keep the source immutable during the benchmark.

All eight tables are compared against the first pysam run, across every field,
with duplicate multiplicities and row order checked. Scientific receipt fields
must also agree. Thus a changed GT, dosage, annotation, sample order, or burden
cannot pass merely by preserving totals. Gzip bytes and timing/path provenance
are not expected to be byte-identical.

After completion, print only the aggregate report:

```bash
python3 - "$RUN_ROOT/results/benchmark.json" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))
assert r['status'] == 'passed', 'Benchmark failed; inspect artifacts locally'
assert len(r['runs']) == 8
assert all(x['passed'] for x in r['comparisons'])
assert r['aggregate_baseline']['samples_in_source'] == 8877
print(json.dumps({k: r[k] for k in (
    'status', 'versions', 'runs', 'timing_summary', 'aggregate_baseline'
)}, indent=2))
PY
```

Wall time covers each child process from launch through completion. CPU time is
child user + system time, not scheduler allocation time. Both include interpreter
startup, candidate loading, validation, extraction, compression, and summary
writing. Scheduler wait, container startup, and comparison are outside the
measurements. No page cache is flushed; no run is claimed to be cold. Inspect
individual timings in order as well as medians/ranges; the reversed order helps
expose cache effects but cannot eliminate shared-filesystem contention.

Only share aggregate agreement and timing. Per-run tables/logs stay local with
restrictive permissions. Use a fresh RUN_ROOT for another benchmark; existing
results are never overwritten or treated as cached extraction runs.

## Compatibility scope

The adapter changes only the reader in an isolated benchmark child. Candidate
selection, exact allele validation, annotations, BGZF writer, and all summaries
reuse production code. The adapter explicitly preserves missing/partial and
haploid calls, dosage, FORMAT sentinels, sample order, FILTER, and star classes.
It also reproduces pysam's missing representation for an out-of-range GT allele
index on a biallelic record; this is baseline parity, not a new genotype policy.
Synthetic cases cover mixed ploidy/phasing, scalar/vector missing FORMAT values,
TBI/CSI indexes, exact-allele rejection, overlapping candidate types, empty
candidates, zero carriers, and failures for malformed source structure.

Successful parity is a prerequisite for interpreting speed. A production backend
change would be a separate decision after the real benchmark; this PR makes none.

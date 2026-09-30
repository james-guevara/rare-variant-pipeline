# Persistent sites catalog → FastVEP/picker → standalone LOFTEE

`annotation.nf` is a separate modular entrypoint. It consumes the published sites catalog;
`sites_catalog.nf`, the legacy workflow, and their configurations are unchanged.
Use **`-C annotation.config`** to exclude the legacy Expanse configuration.

The two processes are `SITES_ANNOTATION:FASTVEP_PICK` and
`SITES_ANNOTATION:STANDALONE_LOFTEE`, tagged with the explicit `unit_id`.
Each selected block runs independently; its LOFTEE task starts as soon as its picker
finishes. There is no whole-catalog VCF inspection or chromosome completion barrier.
The only global input checks are manifest structure, paths, unique unit IDs, selected
IDs, and resource-lock metadata. Both stages durably publish copies.

## Validation status and rollout gate

The user reports manual NBDC validation (Nextflow 26.04.6, Slurm 24.11.4):

| Unit | Input / picked rows | LOFTEE rows |
|---|---:|---:|
| chr22_block0 | 713,983 | 409 |
| chr22_block12 | 48,001 | 332 |
| chr22_block19 | 372,375 | 430 |

The **new Nextflow annotation entrypoint has not yet run on NBDC**. Local tests use
synthetic data; graph tests replace the scientific runner, command tests use native
bcftools with synthetic FastVEP/picker executables and mock the LOFTEE subprocess.
The local suite passed **12 tests** (six annotation tests plus six existing catalog
tests) with Nextflow 26.04.6. These test orchestration, command arguments, old-CSQ removal, contig aliases,
resource checks, failures, selection, cache invalidation, and durable publication;
they do not establish scientific parity. No protected records are committed.

Run block12 first, verify both exact hashes, then run 0/12/19 and compare counts.
Only after those gates pass should all 23 blocks run. Codex cannot execute on NBDC;
the operator runs the commands below and shares receipts/aggregate results only.

## Resource layout and identity

These are the supplied NBDC locations. The LOFTEE subdirectory layout is an explicit
assumption based on the validated container command and portable resource bundle;
the lock command fails if any required file is missing. Resolve that layout before
launch rather than substituting unverified resources.

```text
/home/ood-guevara-james/abcd-fastvep-smoke/
  sites-catalog-smoke/sites/chr22_block{0..22}.sites.vcf.gz
  resources/ensembl-115/
    Homo_sapiens.GRCh38.115.chr22.gff3
    Homo_sapiens.GRCh38.115.chr22.gff3.fastvep.cache
    Homo_sapiens.GRCh38.dna.primary_assembly.fa
    Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
    vep115.transcript-priority.tsv
    vep115.consequence-ranks.tsv
  resources/loftee/
    ensembl-115/transcripts.sqlite
    loftee-grch38/human_ancestor.fa.gz
    loftee-grch38/human_ancestor.fa.gz.fai
    loftee-grch38/human_ancestor.fa.gz.gzi
    loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
    loftee-grch38/loftee.sql
  containers/targeted-7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d.sif
  containers/targeted-81e0e00a97f900867808e5f7cf17cdac6f5c5ef7956e2c6641f3188fa335d9ab.sif
```

The lock generator hashes resources once, verifies the two exact supplied SIF hashes,
and records absolute paths, sizes, nanosecond mtimes, and SHA-256. It requires
`cache_mtime > gff_mtime`; **it never rewrites bytes or timestamps**. Inventory
creation reads large files and warms filesystem caches; this is not a cold-cache
benchmark. Resource hashes are measured identities, not a new validation of their
scientific contents. The SIF pins are checked against the validated images:

| Image | SIF SHA-256 |
|---|---|
| FastVEP/picker | `6eb015d6cb41ae10c64d373b508963496f14a2f57f14d3294cfcd5e06d922d8f` |
| LOFTEE | `a1722375533ad25de58d01635d3cf36667b4ca909149c59404a8deb6491907cd` |

The FastVEP source is `cb8113d7bab2db42cb06bb2b2a40c57b60ea2561`; Rust picker
source is `1671c5a76c369da50c64320d7dc3c719ac2ab95a` (0.1.0). The original wrapper
still invokes the same binaries and flags. Its new optional resource-lock argument
avoids rehashing the shared FASTA in every block; standalone benchmark use retains
its original resource hashing. LOFTEE explicitly invokes
`python /opt/rvp/scripts/run_standalone_loftee.py`, never the image's default command.

## Exact block12 launch

Run from this repository checkout on NBDC. Keep this launch directory, `.nextflow/`,
and `work/` for subsequent `-resume` runs. `python3`, Nextflow, Apptainer, and Slurm
must be available (host Python 3.9 or newer). If NBDC requires runtime modules, load them in your environment
before starting Nextflow; no Expanse module names are assumed.

```bash
set -euo pipefail
BASE=/home/ood-guevara-james/abcd-fastvep-smoke

# All 23 published paths, explicit block IDs; selection happens before VCF inspection.
python3 - "$BASE/sites-catalog-smoke/sites" > abcd-chr22-sites.tsv <<'PY'
import pathlib, sys
root = pathlib.Path(sys.argv[1])
print('unit_id\tchromosome\tvcf')
for block in range(23):
    unit = f'chr22_block{block}'
    print(f'{unit}\tchr22\t{root / (unit + ".sites.vcf.gz")}')
PY

# One-time checksum inventory. Reuse this file while resources remain immutable.
python3 scripts/annotation_resources.py \
  --annotation "$BASE/resources/ensembl-115" \
  --loftee "$BASE/resources/loftee" \
  --chromosomes chr22 \
  --fastvep-sif "$BASE/containers/targeted-7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d.sif" \
  --loftee-sif "$BASE/containers/targeted-81e0e00a97f900867808e5f7cf17cdac6f5c5ef7956e2c6641f3188fa335d9ab.sif" \
  --output abcd-annotation-resources.json

nextflow -C annotation.config run annotation.nf -profile nbdc \
  --sites_manifest "$PWD/abcd-chr22-sites.tsv" \
  --resource_lock "$PWD/abcd-annotation-resources.json" \
  --select_units chr22_block12 \
  --outdir "$BASE/annotation-nextflow" \
  -work-dir "$BASE/annotation-nextflow-work" \
  -with-trace "$BASE/annotation-block12.trace.tsv" \
  -with-report "$BASE/annotation-block12.report.html" \
  -resume

python3 scripts/check_annotation_pilot.py \
  --outdir "$BASE/annotation-nextflow" --units chr22_block12
```

Expected block12 product hashes (TSV file bytes, supplied from the manual pilot):

```text
cbf5c6bb2f469bbb5e6e7bee1f357f325fca6f24b6246f8ea86d262bd4ff6208  fastvep-picker/chr22_block12/picked.tsv
b9fee4e3527371c8ff11b6caaf25de99088f215d20fc063463be737c99a7ca82  loftee/chr22_block12/loftee.tsv
```

If either hash differs, stop and inspect aggregate receipts before expanding.
Do not replace these expected hashes simply to make a check pass.

For the next run, repeat the Nextflow command with
`--select_units chr22_block0,chr22_block12,chr22_block19` and new trace/report paths;
keep the same work/output directories and `-resume`. Block12 should be `CACHED`.
Run the checker with the same three IDs. Blocks0/19 have count baselines, not supplied
byte-level reference hashes. Only then use `--select_units all` for all manifest rows.
No selection is the default: omitting `--select_units` fails instead of launching all blocks.
Subset manifests are also supported; the full 23-row manifest is not required.

## Configuration, outputs, and resume

NBDC settings default to Slurm partition `medium`, QoS `medium`, 2 CPUs, 8 GB, 4h,
and queueSize 4. Both partition and QoS names are assumptions from the brief; if the
site uses a differently named partition, override it. Configurable parameters:
`--annotation_cpus`, `--annotation_memory`, `--annotation_time`,
`--annotation_queue`, `--annotation_qos`, `--annotation_queue_size`.
Override `--annotation_root`, `--loftee_root`, `--fastvep_container`, and
`--loftee_container` when relocating resources, then regenerate the checksum lock.
The profile bind-mounts both shared resource roots **read-only at their absolute
paths**. They are `val` metadata, not staged Nextflow `path` inputs. Only small scripts,
the block VCF, the picked TSV, and receipts are staged. The lock generator resolves symlinks and rejects resources outside these roots; use
canonical root paths in the config if the resource directories themselves are symlinks.

Outputs per block:

```text
<outdir>/fastvep-picker/<unit_id>/picked.tsv
<outdir>/fastvep-picker/<unit_id>/receipt.json
<outdir>/loftee/<unit_id>/loftee.tsv
<outdir>/loftee/<unit_id>/receipt.json
```

Receipts include original source path, block/chromosome, input/output counts and
hashes, bytes, timing, resource identities, and container identity. LOFTEE's receipt
also hashes its actual Python entrypoint and upstream receipt. FastVEP/picker timing
retains the validated wrapper's semantics (startup/resource loading included;
preparation and hashing separate). Linux CPU time is child user+system time;
peak RSS is the largest individual child, **not** summed concurrent pipeline memory.
LOFTEE retains its existing internal LoF selection; its row count is not expected to
equal picked rows. Existing VEP112 CSQ is removed only from temporary FastVEP input.
Empty or multiallelic inputs remain unsupported by the validated benchmark wrapper.

The checksum identities for each unit's chromosome participate in its Nextflow task
key. Adding/removing other selected blocks does not invalidate that block. Small
staged inputs use deep content caching. Launch-time metadata checks cover shared
resources and SIFs; tasks check resource size/mtime before and after execution.
This assumes immutable resources: an external writer that changes bytes while
preserving size and mtime defeats a metadata check. After deliberate resource
replacement, regenerate the lock (rehashing bytes) and rerun. Never mutate resources
in place during a run. The workflow does not alter cache timestamps.

Failure stops the run; Nextflow retains the failed work directory and trace. A
failed adapter writes `receipt.json` with `status: failed` and an error class, removes
partial final TSV, and retains private subprocess logs. **Nextflow does not publish
failed-task outputs**: inspect its reported work directory for `receipt.json`,
`private-failure.log`, and, for FastVEP, `benchmark-output/private-failure.log`.
Keep these with the trace for diagnosis; don't share protected records from logs.
Previously published successful products can remain after a failed rerun, so a
nonzero Nextflow exit must never be interpreted as success from output existence.

The standalone scripts and Nextflow profile do not run genotype extraction,
HIGH/MODERATE prefiltering, Ensembl VEP, or any downstream cohort analysis.

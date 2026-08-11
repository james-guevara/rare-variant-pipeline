# Running the pipeline on G2MH (Expanse)

A runbook. Follow it top to bottom. Every number below was measured on this cohort,
not estimated.

**Cohort:** G2MH WGS, 1,065 samples, DRAGEN 4.4.9 msVCF delivered via NDA, GRCh38.
A single 343 GB whole-genome VCF; per-chromosome extraction happens inside the
pipeline, so you do not pre-split anything.

---

## 0. One-time setup

**SSH with multiplexing** so you authenticate once — see
[sebat-lab-guides/ssh-config.md](https://github.com/james-guevara/sebat-lab-guides/blob/main/ssh-config.md).

**Nextflow environment:**

```bash
micromamba create -n rvp_env nextflow -c bioconda -c conda-forge
micromamba config set use_lockfiles False      # Expanse NFS doesn't handle lock files
```

**Get the Python container.** It is not in the repo and not built on Expanse (the
cluster refuses `singularity build --fakeroot` — no `/etc/subuid` mapping). CI builds
it and publishes to GHCR; you pull it:

```bash
cd <checkout>
sbatch scripts/pull_python_container.sh docker://ghcr.io/james-guevara/rare-variant-pipeline-python:main
```

Run that **as a job**, not on a login node — `singularity pull` shells out to
`mksquashfs`, which dies with `Out of memory (frag_thrd)` under login-node limits.

Then check `params.python_container` in `nextflow.config` points at the `.sif` it wrote.
The bcftools and VEP containers already exist in project space; see the README's
container provenance section for what is and isn't reproducible.

**Clone per cohort.** The convention is one checkout per cohort, so `work/` caches and
launch logs don't collide:

```bash
cd /expanse/projects/sebat1/s3/data/sebat
git clone https://github.com/james-guevara/rare-variant-pipeline rare_variant_pipeline_g2mh
```

---

## 1. Smoke test on chrY first

chrY is the smallest contig, so it exercises the whole chain in minutes. Do not skip
this — it catches container, path and permission problems before you spend hours.

```bash
cd rare_variant_pipeline_g2mh
sbatch -J norm_chrY   launch/run_stage.sh g2mh normalize   chrY
# wait for it, then:
sbatch -J annot_chrY  launch/run_stage.sh g2mh annotate    chrY
sbatch -J carr_chrY   launch/run_stage.sh g2mh carriers    chrY
sbatch -J export_chrY launch/run_stage.sh g2mh export      chrY
sbatch -J post_chrY   launch/run_stage.sh g2mh postprocess chrY
```

Stages are **not** independent — each reads the previous stage's published output, so
wait for one to finish before launching the next. To chain them automatically:

```bash
J=$(sbatch --parsable -J norm_chrY launch/run_stage.sh g2mh normalize chrY)
J=$(sbatch --parsable --dependency=afterok:$J -J annot_chrY launch/run_stage.sh g2mh annotate chrY)
J=$(sbatch --parsable --dependency=afterok:$J -J carr_chrY  launch/run_stage.sh g2mh carriers chrY)
J=$(sbatch --parsable --dependency=afterok:$J -J export_chrY launch/run_stage.sh g2mh export chrY)
     sbatch            --dependency=afterok:$J -J post_chrY  launch/run_stage.sh g2mh postprocess chrY
```

### Expected chrY numbers

| Check | Value |
|---|---|
| normed records | 291,157 → **330,862** (+13.6% from splitting multiallelics) |
| `mismatch_removed` in the norm log | **0** — anything else means a reference mismatch |
| sites-only records | 330,862, **0 samples** |
| final parquet | ~4,100 rows, **59 columns** |

Verify rather than assume:

```bash
B=/expanse/projects/sebat1/s3/data/sebat
bcftools index -n $B/normed_vcfs/g2mh/norm/chrY.norm.vcf.gz     # 330862
grep -o 'total/split/joined/realigned/mismatch_removed[^ ]*' \
     $(ls -t work/*/*/.command.err | head -20) 2>/dev/null | head -2
```

The local-allele fields must be **gone** from the normed VCF, and `AD`/`PL`/`AF` present:

```bash
V=$B/normed_vcfs/g2mh/norm/chrY.norm.vcf.gz
bcftools view -h $V | grep -cE 'ID=(LAD|LPL|LAA|LAF|LF1R2|LF2R1),'   # must be 0
bcftools view -h $V | grep -cE 'ID=(AD|PL|AF),'                      # must be 3
```

If those local tags are still present, the decode did not run and **every multiallelic
site is silently wrong** — stop and check `params.local_alleles` is true for the profile.

---

## 2. Two autosomes next

chrY is haploid and gene-poor, so it does not exercise diploid genotypes, real VEP load,
or the multi-chunk carrier query. chr21 + chr22 do, and are still the smallest autosomes.

```bash
for S in normalize annotate carriers export postprocess; do
  echo "sbatch -J ${S}_test launch/run_stage.sh g2mh $S chr21,chr22"
done
```

### Measured timings (chr21 + chr22 together: 44 min wall)

| Stage | Per chromosome | CPU / memory |
|---|---|---|
| `BCFTOOLS_NORM` | 18m 40s / 19m 04s | ~147% of 4 cpus, RSS < 500 MB (I/O bound) |
| `BCFTOOLS_SITES` | ~3m 45s | ~118%, 24 MB |
| `VEP_ANNOTATE` | 7m 33s / 8m 56s | ~612% of 8 cpus, RSS 1.3–3.3 GB |
| `SPLIT_VEP` | ~12s | 18 MB |
| `PREPARE_VARIANTS` | ~12s | RSS 1.4–1.6 GB |
| `CARRIER_QUERY` | 8s–2m 27s **per chunk** (5 and 9 chunks) | single-threaded, ~15 MB |
| `EXPORT` (3 processes) | ~1s each | ~35 MB |
| `POSTPROCESS` (7 stages) | ~1–4 min total | 4 cpus / 32 GB allocated |

`BCFTOOLS_NORM` dominates and is I/O bound, so more cpus will not help it. chr1 has
~5.65M source records, roughly 5.7× chr22, so it is the long pole for a full run.

### Expected chr22 numbers

Measured before the genotype-QC change (family propagation removed, fail-closed on
missing values), so D2 may differ by a fraction of a percent. D1 counts are exact.


```
source records                      996,414
normed (after splitting)          1,148,541
carriers.tsv rows                   114,226
D1 parquet rows                     113,468   (~0.7% lost to the inner join — expected)
D2 filtered_annotated rows           70,299
```

The 0.7% loss at the merge is normal: `CARRIER_QUERY` uses a region BED, so co-located
non-consequential ALTs leak in and then correctly fail the join against consequential
annotations.

---

## 3. Full run

```bash
CHR=$(printf 'chr%s,' {1..22})chrX,chrY
for S in normalize annotate carriers export postprocess; do
  sbatch -J ${S}_all launch/run_stage.sh g2mh $S "$CHR"
done
```

Chromosomes run in parallel within a stage, so wall-clock is set by the largest
chromosome, not the sum. Driver jobs request 48 h; increase if chr1 needs it.

If a stage fails partway, fix the cause and re-run the same command with `-resume`
appended inside `launch/run_stage.sh` — Nextflow reuses completed tasks. Note that
editing any module's script text invalidates that process's cache for every chromosome.

---

## 4. What you get

| Deliverable | Path |
|---|---|
| **D1** variant × carrier table | `rv_output/g2mh/parquet/<chrom>.merged.parquet` |
| **D2** filtered + fully annotated | `rv_postprocess/g2mh/postprocess/filtered_annotated/<chrom>.filtered_annotated.parquet` |
| **D3** per-sample carrier counts | `rv_postprocess/g2mh/postprocess/counts/per_sample_counts.tsv` |

D2 is the analysis-ready product. It keeps **untiered rows**, so you can stratify by
gene set, consequence class or s_het bin — not only by tier. Add
`WHERE tier IS NOT NULL` if you want tiered rows only.

D3 carries no covariates by design (no PCs, ancestry, case status or family ID). Join
those on downstream. Change what counts are grouped by with `--count_group_col`.

Tier convention: **`t1` is always the most severe**, for both LoF and missense.

---

## 5. Cohort-specific things to know

**Cohort AF thresholds do not transfer from the large cohorts.** At N = 1,065,
`AF < 0.001` means AC ≤ 2 — singletons and doubletons only. The `g2mh` profile uses
`pop_af_max: 0.01` in `scripts/postprocess/resources.json` for this reason, and
`--max_cohort_af 0.01` upstream.

**Site `FILTER` is uninformative here.** DRAGEN's gVCF-genotyper emits `PASS` at every
site (verified: 300,000/300,000 on raw chr22). The real per-sample quality signal is
`FORMAT/FT`, which the pipeline carries into the parquet — about 5% of carrier rows are
flagged `DRAGENHardQUAL` or `LowDepth`.

**Sample IDs are site-prefixed NDA GUIDs** (`886-NDARNY775ZM1`). To join phenotypes,
strip the prefix and match `subjectkey` in
`wgs/jointcall_dragen/genomics_sample03.txt`, which also carries `src_subject_id` (the
original G2MH ID), `sex`, and `library_prep_batch`. That file duplicates its header
row — skip two lines, not one.

**`INFO` counts in the source are partly stale.** `NS`/`GNS`/`GAC`/`GAN` describe the
original 1,066 samples; one was removed afterwards. The pipeline recomputes `AN/AC/AF`
from genotypes, so downstream values are correct — but do not read the raw fields.

**Per-variant missingness** is available as `INFO/F_MISSING`. On chr22 the distribution
is bimodal: median 0, p90 0.199, p99 0.926. Do not reason from the mean (5.9%) — it is
tail-driven. Per-sample missingness is not yet produced by the pipeline.

---

## 6. When something breaks

| Symptom | Cause |
|---|---|
| Task exits 1, empty `.command.out`, `.command.err` says `Command 'ps' required by nextflow` | The container lacks `procps`. Not your code — the script never ran. |
| `tabix: command not found` | A process without a container directive fell through to the host. |
| `Out of memory (frag_thrd)` during `singularity pull` | You ran it on a login node. Use a job with 32 GB. |
| `could not use fakeroot: no mapping entry found in /etc/subuid` | You cannot build images on Expanse. Build via CI, pull the result. |
| `LAD`/`LAA` still in the normed VCF | `params.local_alleles` is not set for the profile. Every multiallelic is wrong. |
| Norm log shows non-zero `mismatch_removed` | `params.reference` disagrees with the VCF's reference. |

Per-stage timings, CPU and memory land in `reports/<prefix>trace.txt` — the right place
to look before resizing anything.

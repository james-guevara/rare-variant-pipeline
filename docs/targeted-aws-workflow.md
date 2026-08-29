# Targeted AWS workflow

Operational failures and their durable fixes are tracked in
[`configuration-troubleshooting.md`](configuration-troubleshooting.md).
The schemas, row grains, identifiers, and recommended PGS integration keys are in
[`output-data-dictionary.md`](output-data-dictionary.md).
The environment-neutral manifest, binding, and executor contract is documented in
[`portable-targeted-execution.md`](portable-targeted-execution.md).

This is the pre-release workflow for querying the lossless G2MH VCZ stores before
annotation. It deliberately keeps scientific annotation separate from cohort-specific
eligibility filtering.

## chr22 validated path

`scripts/run_targeted_chromosome.sh` runs these checkpointed stages. The old
`run_targeted_chr22_aws.sh` name remains only as a compatibility wrapper.

1. Select exact ALT alleles from the sharded Zarr v3 store using target intervals.
2. Emit a sites-only VCF while retaining `variant_index` and `alt_index` pointers.
3. Normalize against the Ensembl 115 chromosome FASTA.
4. Stream FastVEP output directly into the VEP-compatible transcript picker.
5. Run the independently tested standalone LOFTEE implementation.
6. Join GeneBayes and assign `lof_t1` (`post_mean >= 0.18`) or `lof_t2`
   (`0.03 <= post_mean < 0.18`). There is no LoF tier 3.
7. Recover genotypes and carriers from the unchanged Zarr arrays.
8. Apply problematic-region and genotype QC.
9. Annotate population and cohort AF. Filtering, if requested, is a separately named
   final eligibility output.

If `MISSENSE_CANDIDATES` points to the cohort-observed candidate Parquet generated
from dbNSFP, the same run also joins those candidates to the VEP-picked rows, retains
only candidates whose selected transcript is a `missense_variant` in the candidate
gene, assigns `miss_t1` through `miss_t4`, and recovers their carriers. Candidate
scores use the maximum available dbNSFP rank score rather than MANE-only coverage;
MANE remains a transcript-picking priority, not an inclusion filter.

The G2MH chr22 prototype started with 1,443 observed dbNSFP candidates and retained
1,330 (92.2%) after selected-transcript validation: 6 `miss_t1`, 21 `miss_t2`,
203 `miss_t3`, and 1,100 `miss_t4`. Genotype recovery produced 20,425 carrier rows
across all 1,065 samples. The versioned branch reproduced the prototype exactly over
all selected-variant, carrier, and genotype-summary columns. Its pinned counts and
candidate-input checksum are in `resources/g2mh-chr22-missense-regression.json`.

The combined branch was subsequently run from a fresh directory through AWS Batch
job definition revision 2 using container digest
`sha256:d752c78fba93cf104d4ce8a3b49a65702b88c0a7262dbb6f020b5dfd5511a8c0`.
Scientific processing completed in about 31 seconds. The LoF outputs passed the
versioned exact regression gate, and the three missense outputs had zero row or schema
differences from the saved prototype. The run wrote `_SUCCESS` only after both
branches completed.

When `POSTPROCESS_CONFIG` is set, the missense branch continues through problematic
regions, per-genotype QC, population-AF annotation, population eligibility, cohort-AF
annotation, final burden eligibility, and per-sample tier counts. The annotation
tables remain intact. Population eligibility defaults to
`gnomAD4.1_joint_AF < 0.01` (missing passes), while final G2MH cohort eligibility
defaults to `cohort_af < 0.01`; both thresholds are explicit runtime parameters.

AWS Batch revision 3 validated this complete path with container digest
`sha256:01da69fcf6a3aca8518d35ecaeab89a0495f274a1a1495681eadddb710c70637`.
All five downstream Parquet outputs had exact schemas and row multisets relative to
the prototype, and both count TSVs were byte-identical. The final product contains
1,107 carrier rows across 841 alleles and 697 samples. `_SUCCESS` was written after
all upstream and downstream checks completed.

This targeted Batch runner currently invokes no Nextflow runtime. When it is promoted
to chromosome-wide Nextflow orchestration, use the existing environment's exact
controller pin `nextflow/nextflow:26.04.6`; do not rely on the broader repository's
open-ended `>=25.10.0` constraint.

The script accepts environment-variable overrides for every deployment path. Its
defaults describe the validated persistent-EC2/FSx test environment and are not a
portable installation contract.

Each submission must use a fresh `RUN_ROOT`. The runner performs a write/delete
preflight before scientific work, preserves nonempty stage checkpoints for retries,
and writes `_SUCCESS` only after every enabled branch completes. Resubmitting an
already successful run root is an idempotent no-op; completed directories are never
extended with additional branches.

## Container contract

The targeted image contains FastVEP at the pinned `james-guevara/fastVEP` parity
commit `3bdb862b3153a90ef8cca0b07f02f357a15a3eb0`, bcftools, the picker, standalone
LOFTEE, Zarr/Arrow/DuckDB, every pipeline script, and the versioned manifests,
bindings, and regression expectations. Biological
references remain mounted inputs so they can be versioned and cached independently.

The default container paths are:

- `/inputs/chr22.sharded-v3.zarr`
- `/inputs/targets.chr22.bed`
- `/resources/ensembl-115/` (GFF3, transcript cache, FASTA, picker tables)
- `/resources/loftee/` (ancestor FASTA, GERP BigWig, conservation and transcript DBs)
- `/resources/GeneBayes.Supplementary_Table_1.tsv`
- `/output/`

Example Docker invocation:

```bash
docker run --rm \
  -v /path/to/inputs:/inputs:ro \
  -v /path/to/resources:/resources:ro \
  -v /path/to/output:/output \
  ghcr.io/james-guevara/rare-variant-pipeline-targeted:<immutable-tag>
```

Singularity/Apptainer uses the same contract with bind mounts. Remote schedulers need
only an OCI/SIF runtime, sufficient CPU/memory, these mounted resources, and writable
output storage; they do not need Python, Rust, FastVEP, bcftools or `uv` installed on
the host.

## Validation completed

On G2MH chr22, the integrated FastVEP/picker/standalone-LOFTEE branch processed
28,966 targeted alleles. It produced 53 qualifying LoF alleles before region/QC/cohort
eligibility. The qualifying TSV was byte-identical to the earlier AWS canonical result.

Using the historical G2MH final eligibility rule (`cohort AF < 0.01`) solely as a
regression test produced 38 carrier rows across 19 variants, matching the August 2026
Expanse production run exactly for coordinates, samples, genotypes, genes, LOFTEE,
GeneBayes scores and tiers. One consequence label gained the newer
`splice_polypyrimidine_tract_variant` subcategory; it changed no tier or carrier.

### Containerized Batch regression

The 10-gene chr22 pilot (`259` intervals) was rerun through AWS Batch using image
digest `sha256:1cdd3ec46e0e79739f210cc456b733bb366ec90853fa3ab822ac33b6fca491dc`.
It selected `2,551` alleles and produced `13` qualifying LoF variants (`7` `lof_t1`,
`6` `lof_t2`) and `1,262` carrier rows. All annotation fields for those variants and
all carrier/genotype fields were exact matches to the corresponding subset of the
validated full chr22 run.

This regression also established that the forked FastVEP parity commit is required:
the clean pre-patch upstream binary omitted three indel `HGVS_OFFSET` values, lost one
qualifying LoF indel, and returned two fewer carriers in this pilot. The forked source,
not an untracked local binary, is therefore the reproducible annotation dependency.

The full chr22 output can be checked with one command. It requires the new and
validated run roots to be mounted in the same environment:

```bash
python scripts/validate_targeted_workflow_regression.py \
  --new /path/to/new-run \
  --reference /path/to/validated-run \
  --expectations resources/g2mh-chr22-lof-regression.json
```

The check requires exact schemas and row multisets for the tiered-variant TSV,
carrier Parquet, and genotype-summary Parquet. It also enforces the recorded allele,
tier, carrier-row, and distinct-carrier counts. Any difference exits nonzero, making
the command suitable for Batch validation or CI after the mounted reference fixture
has been provisioned.

### Targeted Zarr extraction benchmark

The chr1 G2MH target extraction was profiled on Expanse `ind-shared` against the
sharded Zarr v3 store using 26,160 merged BED intervals. The original implementation
read complete 25,000-record outer shards and issued two scalar `searchsorted` calls
per interval. It took 493 seconds and emitted 123,262 ALT rows from 115,522 records.

Vectorized interval bounds plus selected-row Zarr reads reduced extraction to 153
seconds with one reader (69% faster). Four concurrent readers took 148 seconds, only
five seconds faster, so the portable default remains one worker; `--workers` is
available for storage-specific benchmarking rather than assumed to scale. Both
optimized outputs matched the reference exactly over all nine columns and both
directions of `EXCEPT ALL`. The reproducible Expanse benchmark and validation launchers
are `scripts/benchmark_targeted_extraction_expanse.sh` and
`scripts/validate_extraction_benchmark_expanse.sh`.

## Population and cohort AF contract

- Keep `gnomAD4.1_joint_AF`, `gnomAD4.1_joint_POPMAX_AF`, and the other population
  fields as annotations.
- Do not use MANE Select as an inclusion filter. It is only part of transcript picking.
- Calculate cohort AC/AN/AF and retain the annotation-complete table.
- If an analysis chooses a cohort threshold, emit a second explicitly eligible table.

## chrX and chrY validated policy

chrX and chrY use explicit GRCh38 PAR/non-PAR regions, cohort sample-sex QC, ploidy-
aware genotype interpretation, haploid allele-balance rules, and sex-aware AC/AN/AF
denominators. Calls from samples with ambiguous or discordant sex-chromosome evidence
are retained when burden counting is possible, but are flagged as sensitivity-only;
they never enter primary burden counts silently. This accommodates possible aneuploid
or mosaic samples without discarding their calls or pretending their karyotype is
resolved.

The packaged G2MH regression contracts cover both LoF and missense branches. With
image digest `sha256:edc522a8b6b9dbbdcf25a6fd138a4c5170fca0e47685db95771d9a1a11e908e6`,
AWS Batch revision 3 reproduced every pinned count and all ten canonical hashes:

- chrX: 40 final LoF carrier rows and 847 final missense carrier rows;
- chrY: 1 final LoF carrier row and 14 final missense carrier rows.

After the selected-row Zarr extraction optimization, commit `41de024` was rebuilt as
ECR digest
`sha256:fc501151b89a7623649a76a886b5bb5548ff3ffd91f764bf182e90687d2400c3`.
AWS Batch revision 6 reran chrX from a fresh FSx root and reproduced all 35 counts and
10 hashes exactly. Target extraction selected 39,983 records/44,371 ALT alleles in 25
seconds, and the complete container ran for 72 seconds.

A subsequent full chrX rerun with FastVEP `cb8113d` corrected terminal shifted-
insertion consequence handling. It removed three false pLoF carrier rows at two
alleles and reproduced the earlier full-VEP result exactly: 40 LoF rows (18 `lof_t1`,
22 `lof_t2`). All 847 missense rows and their canonical hashes were unchanged. The
corrected counts and hashes now define the chrX regression contract.

The policy is cohort-neutral; a new cohort supplies its own sample-sex QC binding and
must establish a regression fixture before production use. Unknown or unusual
karyotypes remain reportable sensitivity strata rather than hard-coded exclusions.

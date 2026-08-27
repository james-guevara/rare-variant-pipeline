# Targeted AWS workflow

This is the pre-release workflow for querying the lossless G2MH VCZ stores before
annotation. It deliberately keeps scientific annotation separate from cohort-specific
eligibility filtering.

## chr22 validated path

`scripts/run_targeted_chr22_aws.sh` runs these checkpointed stages:

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

The script accepts environment-variable overrides for every deployment path. Its
defaults describe the validated persistent-EC2/FSx test environment and are not a
portable installation contract.

## Container contract

The targeted image contains FastVEP at the pinned `james-guevara/fastVEP` parity
commit `3bdb862b3153a90ef8cca0b07f02f357a15a3eb0`, bcftools, the picker, standalone
LOFTEE, Zarr/Arrow/DuckDB, and every pipeline script. Biological
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

## Population and cohort AF contract

- Keep `gnomAD4.1_joint_AF`, `gnomAD4.1_joint_POPMAX_AF`, and the other population
  fields as annotations.
- Do not use MANE Select as an inclusion filter. It is only part of transcript picking.
- Calculate cohort AC/AN/AF and retain the annotation-complete table.
- If an analysis chooses a cohort threshold, emit a second explicitly eligible table.

## chrX and chrY are not finalized

The VCF-to-Zarr conversion is lossless for chrX and chrY, but autosomal carrier QC must
not be reused blindly. Before enabling them in production, define and test:

- sample sex/karyotype source and handling of unknown or discordant sex;
- PAR versus non-PAR intervals for the exact GRCh38 reference;
- haploid, diploid and potential aneuploid genotype interpretation;
- sex-aware AC/AN/AF denominators and missingness;
- allele-balance rules for haploid calls;
- chrY inclusion for samples expected to carry Y sequence;
- whether burden models analyze PAR, non-PAR X and Y separately.

Until those decisions are encoded and regression-tested, chrX/chrY stores are valid
query inputs but their variants must not enter the autosomal QC/counting workflow.

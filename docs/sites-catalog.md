# Durable sites-only catalog for preannotated VCF blocks

`sites_catalog.nf` removes genotype/FORMAT columns from population VCF blocks
with existing annotations. It does not normalize, filter, rerun VEP, or launch
any downstream rare-variant stage. Existing cohort entrypoints are unchanged.

## Input and execution

Supply a tab-delimited manifest with explicit, unique, safe `unit_id` values:

```text
unit_id	chromosome	vcf
chr22_block12	chr22	/path/abcd_cohort_chr22_block12.vcf.gz
chr22_block13	chr22	/path/abcd_cohort_chr22_block13.vcf.gz
```

Relative VCF paths resolve beside the manifest. The supplied manifest can
contain one block, several blocks from the same chromosome, or a larger batch.
Unit IDs are never inferred from filenames. Global validation checks only the
manifest and file paths; it never opens VCF contents. Each row launches one
independent `BLOCK_SITES_CATALOG (unit_id)` task, allowing parallel execution
and visible progress. No source index is required.

Use `-C sites_catalog.config` to load the standalone catalog configuration
instead of the full pipeline's Expanse-specific defaults. Native execution
requires Nextflow, bcftools, Bash, awk, sha256sum, base64, and cmp on PATH.
The dedicated config also supports `-profile docker` or `-profile singularity`,
using bcftools 1.22 by default. Override `--bcftools_container` when needed.
For a scheduler, append your own runtime configuration to the `-C` list, e.g.
`-C sites_catalog.config,/path/to/site.config`; this stage assumes no particular
cluster. Task defaults are two CPUs and 2 GB memory per block.

## ABCD chr22_block12 smoke test

Run inside the updated `rare-variant-pipeline` checkout. Replace only the source
path below with the actual path on your execution system:

```bash
abcd_source=/path/abcd_cohort_chr22_block12.vcf.gz
printf 'unit_id\tchromosome\tvcf\nchr22_block12\tchr22\t%s\n' \
  "$abcd_source" > abcd-sites-smoke.tsv

nextflow -C sites_catalog.config run sites_catalog.nf \
  --sites_manifest "$PWD/abcd-sites-smoke.tsv" \
  --outdir "$PWD/abcd-sites-smoke" \
  -work-dir "$PWD/work-abcd-sites" \
  -with-trace "$PWD/abcd-sites-smoke.trace.tsv"

cat abcd-sites-smoke/sites/chr22_block12.receipt.json
bcftools query -l abcd-sites-smoke/sites/chr22_block12.sites.vcf.gz | wc -l
bcftools index --nrecords abcd-sites-smoke/sites/chr22_block12.sites.vcf.gz
```

For the block previously measured by the user, expect **48,001 input/output
records, zero output samples, `csq_present: true`, and `status: PASS`**. Those
real-data observations were provided by the user, not independently validated
by Codex. Inspect the receipt and outputs before choosing a downstream stage.

To add block13, append its explicit row to the manifest and rerun with `-resume`
from the same launch directory and work directory. Use a new trace filename.
Unchanged block12 is reused while block13 runs independently. An expanded or
reduced manifest does not change individual block task keys. Keep each unit ID
tied to its source; changed source content under the same ID replaces that
unit's published product. Use separate output roots for distinct releases.

## Output and validation

```text
<outdir>/sites/
  chr22_block12.sites.vcf.gz
  chr22_block12.sites.vcf.gz.csi
  chr22_block12.receipt.json
```

Publication defaults to `copy`. `--sites_publish_mode link` supports hardlinks
on a shared filesystem; symlink publication is deliberately unsupported. The
published catalog survives task-work cleanup. Retain work and `.nextflow` if
you still want normal Nextflow cache reuse. A subset run leaves other published
units intact.

The task reads the large source once with unfiltered `bcftools view -G -Ov`.
Every input record enters a small, temporary, genotype-free VCF; all subsequent
compression and validation read only that sites representation. No sample
selection, allele trimming, normalization, or filtering options are used.
See the [bcftools view documentation](https://samtools.github.io/bcftools/bcftools.html#view).

Before reporting success, each task checks:

- Readable compressed output with exactly eight VCF columns and zero samples.
- Identical record counts and SHA-256 of canonical CHROM/POS/ID/REF/ALT/QUAL/
  FILTER/INFO records before compression and after reopening the output.
- Unchanged CSQ header definition; full INFO values, including CSQ, are part of
  the record checksum. `csq_present` means the INFO/CSQ header is declared.
- Records agree with the manifest chromosome (optional `chr` prefix and M/MT
  aliases are accepted); the task never silently filters mismatched records.
- A nonempty CSI index whose reported record count matches the output.

Empty valid VCF blocks and inputs without CSQ are supported. An absent CSQ
declaration is reported as false. The per-unit JSON receipt contains `unit_id`,
`chromosome`, `source_vcf`, `input_records`, `output_records`, `output_samples`,
`csq_present`, `sites_sha256`, and `status`. Failed tasks emit no PASS receipt
and do not publish their outputs. Genotypes are not present in the catalog;
preserve source genotyped VCFs for any later carrier analysis.

## Synthetic validation

On 2026-09-29, `pytest -q tests/test_sites_catalog.py` passed all six tests using
Nextflow 26.04.6 and native bcftools 1.17 / htslib 1.24. Tests cover existing
VEP112 CSQ and other INFO fields, QUAL/FILTER/IDs/indels, zero samples, readable
CSI output, empty inputs, absent CSQ, chromosome mismatch, invalid manifests,
an 8,877-sample block, two same-chromosome units, subset expansion/resume, source
paths with spaces/quotes, and durable publication after deleting work.

The configured bcftools 1.22 container and real ABCD data have not been run in
this validation. The next step is the real one/two-block smoke test above.

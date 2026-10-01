# Exact-allele HC-LOFTEE carriers: first block12 gate

This separate `carriers.nf` entrypoint consumes completed LOFTEE TSVs and original
indexed genotype-bearing VCFs. It does not invoke or modify the sites catalog or
FastVEP/picker/LOFTEE annotation workflows. The operator's first real block12 attempt
exposed the ALT=`*` validator bug; the updated implementation has not been rerun on
NBDC by Codex. NBDC execution remains with the operator; start only block12.

## Data flow and definitions

```text
published loftee/<unit_id>/loftee.tsv
  → select LoF == HC, validate unique exact CHROM/POS/REF/ALT keys
original genotype-bearing block VCF + explicit TBI/CSI
  → indexed fetch at candidate positions
  → explicit CHROM/POS/REF/ALT equality
  → called ALT index 1 in GT
  → durable carriers + candidate/matching audit + preliminary HC burden tables
```

The existing `CARRIER_QUERY` is region-based with an optional INFO/AF cap; neither
provides exact allele matching. It and all existing carrier/annotation modules are
unchanged. The new `EXACT_HC_CARRIERS` process uses the `pysam` already present in
the validated LOFTEE SIF. It performs bounded indexed reads and never writes an
intermediate genotype-bearing VCF. The original VCF and index are staged as
symlinks by Nextflow, not copied or rewritten. TSV products are published as copies.

Scientific rules are deliberately explicit:

- Only `LoF == HC`; no additional allele frequency, genotype-quality, site FILTER,
  or FT threshold. This is a raw carrier audit, not a final QC-filtered rare burden.
- One candidate is one exact sequence allele. Duplicate HC keys fail instead of
  double-counting transcripts. Exactly one ALT is required: a sequence allele or
  the literal `*` spanning-deletion allele. REF remains a sequence allele. Other
  symbolic forms, including `<*>`, are not treated as `*` and remain unsupported.
  See the [VCF specification, ALT definition](https://samtools.github.io/hts-specs/VCFv4.5.pdf).
  Carrier, candidate, and unmatched tables include `allele_class` (`sequence` or
  `spanning_deletion`). All HC `*` candidates are retained in the raw audit.
- `22` and `chr22` are allowed naming aliases. If both exist in the source header,
  selection fails rather than silently choosing one. REF/ALT and POS must match
  exactly; no allele normalization or position-only matching is performed here.
- A sample carries a candidate when its GT contains allele index 1. Both `0/1` and
  `1/1` yield one carrier record. `alt_dosage` separately records one versus two ALT
  copies. Phasing, haploid calls, and partial calls (`1/.`) are preserved. Partial
  calls with a called ALT are included and counted separately in the receipt;
  fully missing/reference-only calls do not produce a carrier record.
- GT, GQ, DP, AD, FT, and site FILTER are retained; absent optional fields are `.`.
  A matching record without GT, a multiallelic record at a candidate position, an
  invalid GT allele index, or a duplicate exact source record fails explicitly.
- An HC allele absent from the source is reported as unmatched. A matched allele
  with no observed carriers is a separate category, not an unmatched allele.
- Preliminary sample×gene outputs separate `sequence_hc_variant_count` from
  `spanning_deletion_hc_record_count`; neither is ALT dosage. No combined burden
  column is emitted. A star record is not assumed to be an independent deletion
  event: no upstream-event mapping or deduplication is attempted. Its eventual
  biological burden treatment remains undecided.

## Block12 expected versus observed

Operator-reported annotation baselines:

| Quantity | Expected for real block12 | Available evidence |
|---|---:|---|
| LOFTEE rows | 332 | Not measured on NBDC |
| HC rows / unique candidates | 302 | Operator reports 302 |
| Sequence HC candidates | 203 | Operator reports 203 |
| Spanning-deletion HC candidates (`*`) | 99 | Operator reports 99 |
| LC rows | 30 | Not measured on NBDC |
| Carrier records | Unknown until extraction | Not measured on NBDC |
| Samples / genes with HC carriers | Unknown until extraction | Not measured on NBDC |
| Unmatched exact candidates | Measure and investigate | Not measured on NBDC |

The reported 99 stars comprise 29 `C/*`, 28 `G/*`, 18 `A/*`, 9 `T/*`, and
15 longer-REF `/*` records. This is aggregate operator evidence, not a new Codex
measurement. The old implementation rejected these candidates; the revised
`--expected_hc 302` guard includes both classes (203 + 99), not just sequence alleles.

The operator also reports all 23 chr22 annotation blocks complete (3,502,058 picked
sites; 6,721 LOFTEE rows). That does not establish genotype/carrier counts.

Synthetic fixture observations: six HC candidates, one LC annotation, four exact
matched candidates, two unmatched candidates (wrong REF and absent position), four
carrier records, two carrier samples, one gene with carriers, and two matched
variants without carriers. A different ALT at the same position and an overlapping
deletion do not leak into the carrier table. Five observed ALT copies correspond
to four carrier records because a homozygous ALT counts once for carrier burden.

## Exact NBDC setup and command

From this repository checkout on NBDC, using Python 3 and Nextflow 26.04.6. The
existing annotation resource lock is reused only to identify/check the pinned SIF;
carrier extraction does not mount or consume the large annotation databases.
The SIF checksum verified in that lock must be
`a1722375533ad25de58d01635d3cf36667b4ca909149c59404a8deb6491907cd`.
As for annotation, the lock assumes the SIF remains immutable after verification.

```bash
set -euo pipefail
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
SOURCE=/shared/release/abcd/abcd/concatenated/genetics/sequencing/snv_indel/population_vcf

# Explicit block identity and paths. Resolve only the existing index sidecar.
# CSI takes precedence if both sidecars exist. No VCF-content preflight occurs.
python3 - "$BASE" "$SOURCE" > abcd-block12-carriers.tsv <<'PY'
from pathlib import Path
import sys
base, source = map(Path, sys.argv[1:])
unit = 'chr22_block12'
vcf = source / f'abcd_cohort_{unit}.vcf.gz'
loftee = base / 'annotation-nextflow' / 'loftee' / unit / 'loftee.tsv'
index = next((Path(str(vcf)+ext) for ext in ['.csi','.tbi'] if Path(str(vcf)+ext).is_file()), None)
if not vcf.is_file() or not loftee.is_file() or index is None:
    raise SystemExit('Required source VCF, index, or published LOFTEE TSV is missing')
print('unit_id\tchromosome\tloftee\tvcf\tindex')
print(f'{unit}\tchr22\t{loftee}\t{vcf}\t{index}')
PY

# Use the same resource-lock file used by the validated annotation launch.
# Adjust this argument if the lock lives in a different launch directory.
nextflow -C carriers.config run carriers.nf -profile nbdc \
  --carrier_manifest "$PWD/abcd-block12-carriers.tsv" \
  --resource_lock "$PWD/abcd-annotation-resources.json" \
  --select_units chr22_block12 \
  --expected_hc 302 \
  --outdir "$BASE/carrier-nextflow" \
  -work-dir "$BASE/carrier-nextflow-work" \
  -with-trace "$BASE/carrier-block12.trace.tsv" \
  -resume

# Share aggregate receipt only; do not print the carrier or sample tables to logs.
python3 - "$BASE/carrier-nextflow/carriers/chr22_block12/receipt.json" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))
assert r['status'] == 'passed'
assert r['annotation_rows'] == 332 and r['candidate_hc_variants'] == 302
assert r['non_hc_annotation_rows'] == 30
assert r['by_allele_class']['sequence']['candidate_hc_variants'] == 203
assert r['by_allele_class']['spanning_deletion']['candidate_hc_variants'] == 99
keys = ['candidate_hc_variants','matched_candidate_variants','unmatched_candidate_variants',
        'matched_variants_without_carriers','carrier_records','samples_in_source',
        'samples_with_hc_plof','genes_with_hc_plof','partial_call_carrier_records', 'by_allele_class']
print(json.dumps({k:r[k] for k in keys}, indent=2))
PY
```

Before rerunning, update `feat/exact-allele-carriers` using `git pull --ff-only`.
The implementation change invalidates the carrier task cache; the failed attempt
will execute the corrected script. The sites and annotation workflows are untouched.
Keep `--expected_hc 302`. Use a fresh trace filename for the retry.

Repeat the successful corrected command with another new trace filename and
`-resume` to confirm caching.
Keep the manifest, source paths, launch directory, and work directory stable. The
standard Nextflow path cache uses size/mtime for the large VCF and index; the receipt
also hashes the index, LOFTEE TSV, and all output products. The multi-GB genotype
VCF is deliberately not fully hashed or scanned before every selected task. Source
VCFs and indexes must be immutable. Do not use the annotation workflow's work
folder for this separate launch. **Do not expand to all chr22 blocks yet.**

Configurable resources: `--carrier_cpus` (2), `--carrier_memory` (`4 GB`),
`--carrier_time` (`4h`), `--carrier_queue` (`medium`), `--carrier_qos` (`medium`),
`--carrier_queue_size` (4). Use `-C carriers.config` to isolate legacy Expanse defaults.
The profile assumes the same Apptainer/Slurm environment as annotation. Container
path can be overridden but must match the validated LOFTEE SIF identity in the lock.

## Published files and privacy

Each block publishes under `<outdir>/carriers/<unit_id>/`:

| File | Contents |
|---|---|
| `carriers.tsv.gz` | BGZF TSV; one exact allele × carrier sample, gene/transcript, GT/GQ/DP/AD/FT, ALT dosage, site FILTER |
| `candidates.tsv` | All unique HC allele keys and their picked gene/transcript |
| `unmatched.tsv` | Exact candidate keys not found in the source |
| `samples.tsv` | Full source sample roster, including zero-carrier samples |
| `sample_gene_burden.tsv` | Sample×gene rows with either class nonzero; separate sequence variant and spanning-deletion record counts |
| `sample_burden.tsv` | Separate class counts per sample, including zero counts |
| `gene_burden.tsv` | Carrier records and unique carrier samples separately for each class; includes zero-carrier candidate genes |
| `receipt.json` | Aggregate counts, definitions, paths, identities, hashes, runtime, status |

Receipt schema version 2 includes `by_allele_class.sequence` and
`by_allele_class.spanning_deletion`, each with candidate/matched/unmatched,
carrier-record, unique-sample/gene, partial-call, and ALT-copy counts. Top-level
counts remain raw audit totals across both classes for reconciliation; the legacy
`samples_with_hc_plof` / `genes_with_hc_plof` names describe HC-labelled input records
in that audit, not a decision that each star is an independent biological pLoF.
A sample/gene present in both classes counts only once in the top-level unique total.

**Schema change:** sample and sample×gene tables replace `hc_variant_count` with
`sequence_hc_variant_count` and `spanning_deletion_hc_record_count`. Gene tables use
`sequence_carrier_records`, `sequence_carrier_samples`,
`spanning_deletion_carrier_records`, and `spanning_deletion_carrier_samples`.
Do not add these columns together to infer a deduplicated biological burden.

The TSV products contain protected sample or variant information and belong on
NBDC. Only aggregate receipt fields are printed by the process. Source parser
messages are redirected to `private-input.log` in the task work directory so
Nextflow does not echo record-bearing parser errors into aggregate logs. Do not share
that private log.
Failures remove partial products, retain a failed receipt in the work directory,
and stop the run. As usual, Nextflow does not publish failed-task receipts. A
previously published successful output is not evidence that a later failed run passed.

## Comparison scope

The operator deferred the ABCD-provided matrix comparison; it is not part of this
workflow. No matrix loader, column-encoding assumptions, or agreement/correlation
results are added. The per-sample/per-gene tables remain useful standalone outputs.

Block12-only counts are not whole-chr22 burdens. They omit variants from 22 other
blocks, potentially including other variants in the same genes. Any future
comparison must account for coverage and distinguish Ensembl 115 + FastVEP +
Rust picker + LOFTEE HC from ABCD's VEP112 LoF/HIGH definition. Discrepancies
between these annotation systems are not automatically errors.

## Local validation

`uv run --with pytest --with pysam python -m pytest -q tests/test_exact_carriers.py tests/test_sites_annotation.py tests/test_sites_catalog.py`
passed 30 tests (18 carrier tests plus 12 existing regression tests), using real
pysam indexed VCF reads and real Nextflow execution for the carrier stage. Carrier
tests cover exact REF/ALT matching, same-position ALTs, overlapping records,
HC/LC selection, phased/homozygous/haploid/partial GTs, optional FORMAT values and
headers, TBI and CSI, duplicate/multiallelic failure, empty HC sets, expected-HC
guards, subset/resume caching, and durable outputs. Star-specific tests cover
longer REF/*, phased/homozygous/partial star GTs, different-ALT/different-REF
exclusion, multiallelic `C,*` rejection, `<*>` rejection, and separated burden counts. They use synthetic data;
there has been no carrier run on real ABCD data by Codex. NBDC configuration parsing,
Python compilation, and diff whitespace checks also passed.

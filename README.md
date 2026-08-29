# Rare Variant Pipeline

A modular Nextflow DSL2 pipeline for annotating and analyzing rare genetic variants in family-based sequencing studies (WGS/WES).

## Overview

This pipeline processes VCF files through variant annotation (VEP), filters for consequential variants (HIGH/MODERATE impact), and queries family genotypes to identify rare variants segregating in families.

## Pipeline Structure

Five subworkflows, each with its own entry point so stages can be run and re-run
independently.

```
NORMALIZE            BCFTOOLS_NORM
                     [decode local alleles] -> norm -m- -> fill-tags(AN,AC,AF,F_MISSING)
                          |
                          v  <chrom>.norm.vcf.gz
ANNOTATE             BCFTOOLS_SITES -> VEP_ANNOTATE -> SPLIT_VEP -> PREPARE_VARIANTS
                          |
                          v  <chrom>.consequential.{tsv,bed}
CARRIER_EXTRACTION   SCATTER_BED -> CARRIER_QUERY (per chunk) -> GATHER_CARRIERS
                          |
                          v  <chrom>.carriers.tsv
EXPORT               MERGE_ANNOTATIONS -> SORT_INDEX -> CONVERT_PARQUET
                          |
                          v  D1: <chrom>.merged.parquet
POSTPROCESS          PP_FILTER_REGIONS -> PP_QC_GENOTYPE -> PP_JOIN_POP_AF
                     -> PP_JOIN_SCORES -> PP_JOIN_GENE_CONSTRAINT
                     -> PP_TIER_VARIANTS (annotate)
                          |                        \
                          v                         v
                     D2: filtered_annotated/   D3: PP_COUNT_CARRIERS -> counts/
```

`CARRIER_EXTRACTION` reads no pedigree despite the historical name of its
predecessor — it extracts carriers (`GT="alt"`), nothing more.

## Directory Structure

```
├── main.nf                    # entry points + parameter defaults
├── nextflow.config            # profiles, resources, containers
├── Dockerfile                 # Python-stages image (built by CI)
├── .github/workflows/         # container build + publish to GHCR
├── modules/                   # one process each
│   ├── bcftools_norm.nf       #   decode/split/fill-tags
│   ├── bcftools_sites.nf      #   drop genotypes
│   ├── vep_annotate.nf
│   ├── split_vep.nf
│   ├── prepare_variants.nf
│   ├── scatter_bed.nf
│   ├── carrier_query.nf
│   ├── gather_carriers.nf
│   ├── merge_annotations.nf
│   ├── sort_index.nf
│   └── convert_parquet.nf
├── subworkflows/
│   ├── normalize.nf
│   ├── annotate.nf
│   ├── carrier_extraction.nf
│   ├── export.nf
│   └── postprocess.nf
├── scripts/
│   ├── prepare_variants.py
│   ├── merge_genotypes_annotations.py
│   ├── tsv_to_parquet.py
│   ├── pull_python_container.sh      # one-time container setup (run via sbatch)
│   ├── download_*.sh                 # fetch region/regulatory BED resources
│   └── postprocess/                  # DuckDB stages + resources.json
│       ├── filter_regions.py
│       ├── qc_genotype.py
│       ├── join_pop_af.py
│       ├── join_scores.py
│       ├── join_gene_constraint.py
│       ├── tier_variants.py
│       ├── count_carriers.py
│       └── resources.json
└── archive/                   # superseded code, kept findable (see archive/README.md)
```

## Requirements

- Nextflow >= 25.10.0
- Singularity / SingularityPRO

Nothing else. **Every process runs in a container**, so there is no host Python,
bcftools or tabix dependency. Earlier versions relied on a micromamba env at an
absolute path inside one user's home directory, which meant nobody else could run
the pipeline; that is gone.

### Containers

| Image | Processes |
|---|---|
| `bcftools:1.22` | `BCFTOOLS_NORM`, `BCFTOOLS_SITES`, `CARRIER_QUERY`, `SPLIT_VEP`, `SCATTER_BED`, `GATHER_CARRIERS`, `SORT_INDEX` |
| `ensembl-vep:115.2` (with samtools) | `VEP_ANNOTATE` |
| `rare-variant-pipeline-python` | `PREPARE_VARIANTS`, `MERGE_ANNOTATIONS`, `CONVERT_PARQUET`, all `PP_*` |

The Python image is built from the repo's `Dockerfile` by CI
(`.github/workflows/container.yml`) and published to GHCR. Versions are pinned to
duckdb 1.5.2 / polars 1.39.3 / polars-bio 0.28.0 / pyarrow 22.0.0 — **not** latest,
because `prepare_variants.py` is written against the polars_bio 0.28 API. The build
asserts that API surface exists, so a bad pin fails CI rather than a long run.

#### Container provenance (known gap)

Only one of the three is currently reproducible from this repo.

| Container | Provenance |
|---|---|
| `bcftools:1.22--h3a4d415_1` | Stock biocontainer. Recoverable: `docker://quay.io/biocontainers/bcftools:1.22--h3a4d415_1` |
| `ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools` | **Custom, no build recipe.** Built 2025-11-25 from a local sandbox (`bootstrap: localimage`, `from: vep_sandbox`). Contains VEP 115.2 plus samtools 1.22.1, bgzip, tabix and perl in `/usr/local/bin`. samtools is there because **LOFTEE requires it** — the plugin shells out to samtools for ancestral-allele lookups, so the stock VEP image is not sufficient. If the `.sif` is lost there is currently no documented way to rebuild it. |
| `rare-variant-pipeline-python` | Built from this repo's `Dockerfile` by CI. |

Both older images live in shared project space, so anyone with `sebat1` access can run
the pipeline today. Making the VEP image reproducible — a `Dockerfile` starting from
the `ensembl-vep:115.2--pl5321h2a3209d_1` biocontainer and adding samtools 1.22.1,
built by the same CI workflow — is outstanding work, and may be superseded if the
LOFTEE step is reimplemented (LOFTEE is the only reason samtools is needed). Note a rebuild would be
*equivalent*, not bit-identical, so it should be validated by re-running a
chromosome and comparing VEP output before being swapped in.

Separately, `VEP_CACHE` is ~48 GB of Ensembl data and is not containerized. It is
downloaded from Ensembl, not shipped here.

One-time setup on the cluster:

```bash
sbatch scripts/pull_python_container.sh   # writes rvp_python.sif next to the others
```

Run it **as a job**. `singularity pull` shells out to `mksquashfs`, which dies with
`Out of memory (frag_thrd)` under login-node limits but succeeds with 32 GB in a job.
Building the image on the cluster is not possible: there is no `/etc/subuid` mapping
for user accounts, so `singularity build --fakeroot` is refused — hence CI.

## Usage

### Reusable targeted workflow

`targeted.nf` remains the standalone entrypoint, but is now a thin wrapper around
the named DSL2 workflow in `workflows/targeted_manifest.nf`:

```nextflow
include { TARGETED_MANIFEST_WORKFLOW } from './workflows/targeted_manifest'

workflow {
    manifest_bindings = Channel.of(tuple(file(params.manifest), file(params.bindings)))
    TARGETED_MANIFEST_WORKFLOW(manifest_bindings)
}
```

The reusable workflow takes a channel of `(science_manifest, environment_bindings)`
file tuples. The standalone wrapper creates a one-element channel; a cohort-level
parent can provide one tuple per chromosome and Nextflow will scatter the same
validated process without copying its implementation. Each binding supplies its own
unique run root. `params.targeted_container` remains the single digest-pinned runtime
for every tuple in a launch.

The named workflow uses the same manifest, environment bindings, container,
preflight, and run-root parameters as the standalone command and emits both
`execution_receipt` and a portable `chromosome_outputs` directory. The latter
contains the available validated chromosome-level LoF, missense, primary LoF, and
sensitivity LoF burden TSVs plus a checksummed `targeted-output-manifest.json`.
A future cohort-level parent workflow can therefore invoke
the rare-variant and PGS branches independently and join their participant-level
tables only after both branches complete. Scheduler and container-runtime settings
remain site-specific configuration rather than part of the scientific workflow.

`RARE_BURDEN_GATHER_WORKFLOW` in `workflows/rare_burden_gather.nf` consumes a cohort
sample manifest and the completed chromosome-output packages. It emits
`rare_burdens.tsv` and `rare_burdens_by_chromosome_stratum.tsv`. The thin
`gather.nf` entrypoint allows this aggregation contract to be run and tested on its
own before composition with the PGS workflow.

`cohort.nf` connects those two reusable pieces without adding another scientific
implementation. Its tab-separated run sheet has exactly three columns:

```text
chromosome\tmanifest\tbindings
chr1\t/path/to/g2mh-chr1.json\t/path/to/aws-g2mh-chr1.json
chrX\t/path/to/g2mh-chrX.json\t/path/to/aws-g2mh-chrX.json
```

Each binding must name a unique chromosome run root. The wrapper scatters all rows,
passes the resulting checksummed chromosome packages directly to the gather workflow,
and publishes `rare_burdens.tsv` plus
`rare_burdens_by_chromosome_stratum.tsv` under `${params.outdir}/rare_burdens`.
`--expected_chromosomes` is required and must match the packages exactly, preventing
an incomplete cohort from silently producing zero-filled burdens.

### PGS and rare-burden integration

`INTEGRATED_ANALYSIS_WORKFLOW` in `workflows/integrated_analysis.nf` is the narrow
join boundary between this repository and the reusable `PGS_WORKFLOW` from
`james-guevara/pgs_pipeline`. It consumes the PGS `analysis_dataset` and
`analysis_dictionary` outputs plus this pipeline's gathered `rare_burdens` table.
The standalone `integrate.nf` wrapper accepts the same artifacts as paths.

The join deliberately uses the PGS analysis dataset as its participant universe, in
the original row order. Every PGS IID must be present in the completed wide
rare-burden table; missing participants or blank/non-integer burden values fail the
job instead of being silently converted to zero. Rare-only participants are excluded
but counted in `integrated_analysis_qc.json`. Nonempty `FID` and `SEX` values must
agree if both inputs provide them.

Outputs under `${params.outdir}/integrated_analysis` are:

- `integrated_analysis_dataset.tsv`, containing PGS/PCA variables and `lof_t1`,
  `lof_t2`, and `miss_t1` through `miss_t4`;
- `integrated_analysis_dictionary.tsv`, following the merged variable template;
- `integrated_analysis_qc.json`, recording cohort ID, input/output participant counts,
  excluded rare-only participants, and the analysis-universe rule.

```bash
nextflow -C targeted.config run integrate.nf -profile local_docker \
  --pgs_dataset /path/to/analysis_dataset.tsv \
  --pgs_dictionary /path/to/analysis_dataset_dictionary.tsv \
  --rare_burdens /path/to/rare_burdens.tsv \
  --cohort_id g2mh \
  --targeted_container docker.io/example/targeted@sha256:DIGEST \
  --outdir results/g2mh
```

**New to this pipeline?** [`docs/running-g2mh.md`](docs/running-g2mh.md) is a runbook you
can follow top to bottom — one-time setup, a chrY smoke test with the numbers to check
against, then the full run. The rest of this section is reference material.

### Reference

### Full Pipeline
```bash
nextflow run main.nf -profile <cohort> --chroms <chromosomes> -resume
```

### Individual Subworkflows

Entry points, in dependency order:

```bash
nextflow run main.nf -profile ssc -entry RUN_NORMALIZE              --chroms chr22
nextflow run main.nf -profile ssc -entry RUN_ANNOTATE               --chroms chr22  # sites -> VEP -> split-vep -> prepare
nextflow run main.nf -profile ssc -entry RUN_CARRIER_EXTRACTION     --chroms chr22  # scatter -> carrier query -> gather
nextflow run main.nf -profile ssc -entry RUN_EXPORT                 --chroms chr22  # merge -> sort/index -> parquet
nextflow run main.nf -profile ssc -entry RUN_POSTPROCESS            --chroms chr22  # filter/annotate -> tier -> counts
```

`RUN_ANNOTATE` can also be run one step at a time, which is what the per-step
launch scripts do so each stage gets its own outdir and trace:

```bash
-entry RUN_SITES_ONLY          # normed -> sites-only VCF
-entry RUN_VEP_ONLY            # sites  -> VEP-annotated VCF
-entry RUN_SPLIT_VEP_ONLY      # VEP    -> per-consequence TSV
-entry RUN_PREPARE_VARIANTS_ONLY  # TSV -> reformatted + consequential TSV/BED
```

Renamed 2026-08-10 (old names are gone, not aliased): `RUN_VCF_PROCESSING` →
`RUN_ANNOTATE`, `RUN_FAMILY_PROCESSING` → `RUN_CARRIER_EXTRACTION`,
`RUN_MERGE_INDEX` → `RUN_EXPORT`. The old names were misleading: nothing in the
pipeline reads a pedigree, so "family processing" only ever extracted carriers.

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--chroms` | Comma-separated chromosomes | `chr22` |
| `--outdir` | Output directory | `output` |
| `--tmpdir` | `SORT_INDEX` scratch. Must be on a shared filesystem — **never `/tmp`**, which is node-local on Expanse | `${outdir}/tmp` |
| `--regions_per_chunk` | Regions per scatter chunk | `1000` |
| `--trace_prefix` | Prefix for report files | `""` |
| `--single_vcf` | Path to one whole-genome VCF; per-chrom extraction happens in `BCFTOOLS_NORM`. Leave null for per-chrom inputs | `null` |
| `--target_bed` | Optional BED on the shared site filesystem. `BCFTOOLS_NORM` uses the source VCF index to extract these intervals before decoding, normalization, VEP, and carrier extraction. On AWS, place it on the same FSx mount visible to the controller and Batch workers. | `null` |
| `--local_alleles` | Input stores local-allele FORMAT fields — see below | `false` |
| `--pp_cohort` | Key into `scripts/postprocess/resources.json` `cohorts` | set per profile |
| `--count_group_col` | Column `PP_COUNT_CARRIERS` stratifies by | `tier` |

### Input dialects (DRAGEN msVCF)

### Targeted early extraction

Set `--target_bed` to restrict processing at the first VCF operation. The source
VCF must be bgzip-compressed and indexed. This is materially different from
filtering after VEP: only records overlapping the BED are decoded, normalized,
annotated, queried for carriers, and postprocessed.

The BED is an optimization boundary, not a tier definition. Generate it from a
versioned transcript annotation with enough exon/splice padding to preserve every
variant that could receive the intended consequence, and validate a targeted run
against an existing full run before treating it as equivalent. A BED can represent
GeneBayes-eligible genes, a curated gene set, or their union.

```bash
python scripts/build_target_bed.py \
  --genebayes GeneBayes.Supplementary_Table_1.tsv \
  --gtf Homo_sapiens.GRCh38.115.gtf.gz \
  --output targeted-lof-shet-ge-0.03.bed \
  --min-post-mean 0.03 \
  --padding 8 \
  --add-chr-prefix

nextflow run main.nf -profile <cohort> \
  --chroms chr22 \
  --target_bed /shared/resources/targeted-lof.v1.bed \
  -resume
```

Use `--add-chr-prefix` only when the source VCF contigs are named `chr1`,
`chr2`, and so on while the GTF uses `1`, `2`, and so on. A contig-name
mismatch yields an empty extraction, so the smoke test must assert a nonzero
record count.

Callers differ in what they store per sample. Both flags below default to `false`,
and with both off `BCFTOOLS_NORM` emits the original single-command form, so task
hashes — and therefore `-resume` — are unchanged for existing cohorts.

**`--local_alleles`** — DRAGEN's gVCF-genotyper msVCF (`--generate-msvcf`) stores
per-sample depths and likelihoods as *local-allele* fields `LAD`/`LPL`/`LAA`
instead of `AD`/`PL`. These are `Number=.`, so `bcftools norm -m-` has no rule for
subsetting them and **copies them verbatim onto every split row, silently**. The
`LAA` indices then no longer describe the row they sit on, and a non-carrier
inherits the carrier's depth (`AD=0,32` on a `0/0` row → allele balance 1.0).

So the decode **must** precede the split. This is a property of bcftools, verified
identical on 1.20 and 1.24 — not a version bug to wait out. When the flag is set,
`BCFTOOLS_NORM` runs `+tag2tag -- --LXX-to-XX -r` (only the family form works; the
single-tag `--LAD-to-AD` is rejected in both versions despite its own `--help`),
then strips the remaining allele-indexed fields `LAF`/`LF1R2`/`LF2R1`, and only
then splits.

`BCFTOOLS_NORM` always ends with `+fill-tags -t AN,AC,AF,F_MISSING`. Both parts are
deliberate and unconditional:

- **`AN,AC,AF` are recomputed from `GT`.** Downstream stages filter on `INFO/AF`, so
  if it were the caller's own value the filter would mean something different per
  cohort — callers compute it over different sample sets and with different ploidy
  handling, and it can be stale (e.g. after samples are dropped post-joint-calling).
  Recomputing makes it mean "frequency among the samples in this file", everywhere.
- **`F_MISSING`** is the per-variant missing call rate, which the pipeline otherwise
  has no measure of at all. `prepare_variants.py` extracts it as a column.

Genotypes already stream through this step, so neither costs an extra pass.

### Profiles

- `spark_wgs` — SPARK WGS (DeepVariant, per-chrom)
- `spark_iwes` — SPARK iWES v3 (DeepVariant, per-chrom)
- `ssc` — SSC (GATK, per-chrom)
- `cirm_dragen`, `cirm_gatk` — CIRM
- `g2mh_wo_utrecht` — G2MH, in-house GATK/VQSR joint call (superseded)
- `g2mh` — G2MH finalized DRAGEN msVCF from NDA; sets `single_vcf`,
  `local_alleles`, and the postprocess cohort key

Note `--max_cohort_af` does not transfer across cohort sizes: at N≈1,000 an
AF of 0.001 means AC ≤ 2, i.e. singletons and doubletons only.

## Deliverables

The pipeline produces three things. Naming them explicitly matters because the
middle one is the analysis-ready product and is easy to lose.

**D1 — variant × carrier table.** `RUN_EXPORT` → `<outdir>/parquet/<chrom>.merged.parquet`.
One row per (consequential rare variant × carrier sample) with VEP annotations,
`GT/GQ/DP/AD/PL/FT`, `FILTER`, and recomputed `AC/AN/AF/F_MISSING`.

**D2 — filtered, fully annotated variant × carrier table.**
`RUN_POSTPROCESS` → `<outdir>/postprocess/filtered_annotated/<chrom>.filtered_annotated.parquet`.
D1 after region filtering, genotype QC, population-AF annotation, dbNSFP scores
and gene constraint, plus a `tier` column. Population AF never removes rows in
this stage; cohort-specific eligibility is a separate, explicit final output.

**This is deliberately stratification-agnostic.** `tier_variants` *annotates* and
does not filter, so untiered rows are retained. Tiering is only one way to slice the
data — gene sets, consequence classes and s_het bins are others — and dropping
untiered variants here would foreclose all of them. Consumers wanting only tiered
rows add `WHERE tier IS NOT NULL`.

**D3 — per-sample carrier counts.** `<outdir>/postprocess/counts/`:
`per_sample_counts.tsv` (one row per sample, one column per group value, explicit
zeros) and `group_totals.tsv`. Grouping is set by `--count_group_col`, default
`tier`; any annotated column works.

D3 carries **no covariates** — no PCs, ancestry, case status or family ID — and does
no cross-cohort union. Those belong to analysis, and keeping them out means the
pipeline has no dependency on an analysis manifest. It also means a sample appearing
in two callsets cannot be double-counted here.

### Tier convention

`t1` is always the **most severe** tier, for both LoF and missense:

| Tier | Definition |
|---|---|
| `lof_t1` | LoF = `HC` and s_het ≥ 0.18 |
| `lof_t2` | LoF = `HC` and s_het ∈ [0.03, 0.18) |
| `miss_t1`…`miss_t4` | missense with 4, 3, 2, 1 of {ClinPred, AlphaMissense, popEVE, MPC} rankscores over threshold |

Counts inherit this by construction — `PP_COUNT_CARRIERS` groups on the `tier`
column written by `tier_variants.py`, with no relabelling step that could invert
severity.

### What postprocessing deliberately does NOT do

- **No DNM rescue.** It can only fire where trios exist, so the same variant would be
  filtered differently between cohorts with and without them. De novo status belongs
  in a separate annotation step.
- **No `FILTER` cut.** `FILTER` is carried as a column so the cut can be applied
  downstream where it is visible. (Note some callers emit a constant `PASS` at every
  site and put their real per-sample filtering in `FORMAT/FT`, which is why `FT` is
  carried too.)
- **No family-level propagation in genotype QC.** Dropping a whole (variant × family)
  group when one member fails makes effective stringency depend on how many family
  members happen to carry that variant, which varies by cohort and by family. QC is
  strictly per-row, and fails closed on missing `GQ`/`DP`/`AB`.
- **No pedigree is read anywhere.**

## Output

Final output is a gzipped, tabix-indexed TSV with:
- Variant annotations (VEP consequences, impact, gene, etc.)
- Family genotypes (GT, GQ, DP, AD per sample)
- Family relationships

```
output/
├── sites/          # Filtered VCF sites
├── vep/            # VEP annotated VCFs
├── variants/       # Split VEP TSV
├── reformat/       # Reformatted + consequential variants
├── scatter/        # BED chunks
├── family_query/   # Per-chunk genotypes
├── gather/         # Gathered genotypes per chromosome
├── resolve/        # Resolved family genotypes
├── merged/         # Merged annotations + genotypes
└── indexed/        # Final gzipped + tabix indexed output
```

### Output Columns

Each row in the final TSV represents one **sample x ALT allele x variant** — multi-allelic sites and multi-member families expand into multiple rows. Columns come from three sources, joined on `(#CHROM, POS0, END, REF, ALT_specific)`.

#### Genotype columns (from `family_query.py` + `resolve_family_genotypes.py`)

| Column | Description |
|--------|-------------|
| `#CHROM` | Chromosome |
| `POS0` | 0-based start position |
| `END` | End position (exclusive) |
| `POS` | 1-based position |
| `REF` | Reference allele |
| `ALT` | All ALT alleles (comma-separated) |
| `SAMPLE` | Sample ID |
| `GT` | Genotype (e.g., `0/1`) |
| `GQ` | Genotype quality |
| `DP` | Read depth |
| `AD` | Allelic depths (comma-separated) |
| `FAMILY` | Family ID from PED file |
| `allele_num` | Which ALT allele this row represents (1-based index) |
| `ALT_specific` | The specific ALT allele for this row |
| `carrier` | 1 if sample carries this allele, 0 if not, null if missing GT |
| `AD_ref` | Reference allele depth |
| `AD_alt` | ALT allele depth (allele-specific) |
| `sample_is_multiallelic_het` | True if sample is het for >1 different ALT allele |
| `family_has_multiallelic_het` | True if any family member is multiallelic het |
| `alleles_in_gt` | Comma-separated ALT allele indices in this sample's GT |
| `is_missing_gt` | True if genotype contains `.` |

#### VEP annotation columns (from `split_vep.nf` + `reformat_variants.py`)

VCF-level fields:

| Column | Description |
|--------|-------------|
| `ID` | Variant ID |
| `QUAL` | Variant quality score |
| `INFO` | INFO field (CSQ blob stripped) |

Standard VEP fields:

| Column | Description |
|--------|-------------|
| `Allele` | VEP allele |
| `Consequence` | VEP consequence type(s) |
| `IMPACT` | HIGH, MODERATE, LOW, or MODIFIER |
| `SYMBOL` | Gene symbol |
| `Gene` | Ensembl gene ID |
| `Feature_type` | Type of feature (Transcript, RegulatoryFeature, etc.) |
| `Feature` | Ensembl transcript/feature ID |
| `BIOTYPE` | Transcript biotype |
| `EXON` | Exon number/total |
| `INTRON` | Intron number/total |
| `HGVSc` | HGVS coding notation |
| `HGVSp` | HGVS protein notation |
| `cDNA_position` | Position in cDNA |
| `CDS_position` | Position in CDS |
| `Protein_position` | Position in protein |
| `Amino_acids` | Reference/variant amino acids |
| `Codons` | Reference/variant codons |
| `Existing_variation` | Known variant IDs (e.g., rsIDs) |
| `ALLELE_NUM` | Allele number from VCF |
| `DISTANCE` | Distance to transcript |
| `STRAND` | Transcript strand |
| `FLAGS` | Transcript quality flags |
| `SYMBOL_SOURCE` | Source of gene symbol |
| `HGNC_ID` | HGNC identifier |
| `CANONICAL` | YES if canonical transcript |
| `MANE_SELECT` | MANE Select transcript ID |
| `MANE_PLUS_CLINICAL` | MANE Plus Clinical transcript ID |
| `NEAREST` | Nearest gene symbol |
| `DOMAINS` | Protein domains |
| `PUBMED` | PubMed IDs |
| `Uploaded_allele` | Allele as uploaded |

LOFTEE plugin:

| Column | Description |
|--------|-------------|
| `LoF` | LOFTEE verdict: HC, LC, or empty |
| `LoF_filter` | LOFTEE filter reason |
| `LoF_flags` | LOFTEE warning flags |
| `LoF_info` | LOFTEE detailed info |

SpliceAI plugin:

| Column | Description |
|--------|-------------|
| `SpliceAI_pred` | SpliceAI prediction scores |

MaxEntScan plugin (with SWA, NCSS):

| Column | Description |
|--------|-------------|
| `MaxEntScan_ref` | MaxEntScan reference score |
| `MaxEntScan_alt` | MaxEntScan alternate score |
| `MaxEntScan_diff` | MaxEntScan score difference |

dbNSFP plugin — pathogenicity scores:

| Column | Description |
|--------|-------------|
| `MPC_score` | MPC missense constraint score |
| `MPC_rankscore` | MPC rank score |
| `PrimateAI_score` | PrimateAI pathogenicity score |
| `PrimateAI_rankscore` | PrimateAI rank score |
| `PrimateAI_pred` | PrimateAI prediction |
| `ClinPred_score` | ClinPred pathogenicity score |
| `ClinPred_rankscore` | ClinPred rank score |
| `ClinPred_pred` | ClinPred prediction |
| `AlphaMissense_score` | AlphaMissense pathogenicity score |
| `AlphaMissense_rankscore` | AlphaMissense rank score |
| `AlphaMissense_pred` | AlphaMissense prediction |
| `CADD_raw` | CADD raw score |
| `CADD_raw_rankscore` | CADD raw rank score |
| `CADD_phred` | CADD phred-scaled score |

dbNSFP plugin — population frequencies:

| Column | Description |
|--------|-------------|
| `1000Gp3_AC` | 1000 Genomes allele count |
| `1000Gp3_AF` | 1000 Genomes allele frequency |
| `AllofUs_ALL_AF` | All of Us overall AF |
| `AllofUs_POPMAX_AF` | All of Us population max AF |
| `AllofUs_POPMAX_POP` | All of Us population max population |
| `RegeneronME_ALL_AF` | Regeneron ME overall AF |
| `gnomAD4.1_joint_flag` | gnomAD v4.1 joint flag |
| `gnomAD4.1_joint_AF` | gnomAD v4.1 joint AF |
| `gnomAD4.1_joint_nhomalt` | gnomAD v4.1 joint homozygous alt count |
| `gnomAD4.1_joint_POPMAX_AF` | gnomAD v4.1 population max AF |
| `gnomAD4.1_joint_POPMAX_nhomalt` | gnomAD v4.1 population max homozygous alt count |
| `ALFA_Total_AF` | ALFA total allele frequency |
| `dbNSFP_POPMAX_AF` | dbNSFP cross-database population max AF |
| `dbNSFP_POPMAX_AC` | dbNSFP cross-database population max AC |
| `dbNSFP_POPMAX_POP` | dbNSFP population max population |

dbNSFP plugin — ClinVar:

| Column | Description |
|--------|-------------|
| `clinvar_id` | ClinVar variant ID |
| `clinvar_clnsig` | ClinVar clinical significance |
| `clinvar_trait` | ClinVar trait name |
| `clinvar_review` | ClinVar review status |
| `clinvar_hgvs` | ClinVar HGVS notation |
| `clinvar_var_source` | ClinVar variant source |
| `clinvar_MedGen_id` | ClinVar MedGen ID |
| `clinvar_OMIM_id` | ClinVar OMIM ID |
| `clinvar_Orphanet_id` | ClinVar Orphanet ID |

dbNSFP plugin — conservation scores:

| Column | Description |
|--------|-------------|
| `GERP++_NR` | GERP neutral rate |
| `GERP++_RS` | GERP rejected substitutions |
| `GERP++_RS_rankscore` | GERP RS rank score |
| `GERP_92_mammals` | GERP score across 92 mammals |
| `GERP_92_mammals_rankscore` | GERP 92 mammals rank score |
| `phyloP100way_vertebrate` | phyloP 100-way vertebrate conservation |
| `phyloP100way_vertebrate_rankscore` | phyloP 100-way rank score |
| `phyloP470way_mammalian` | phyloP 470-way mammalian conservation |
| `phyloP470way_mammalian_rankscore` | phyloP 470-way rank score |
| `phyloP17way_primate` | phyloP 17-way primate conservation |
| `phyloP17way_primate_rankscore` | phyloP 17-way primate rank score |
| `phastCons100way_vertebrate` | phastCons 100-way vertebrate |
| `phastCons100way_vertebrate_rankscore` | phastCons 100-way rank score |
| `phastCons470way_mammalian` | phastCons 470-way mammalian |
| `phastCons470way_mammalian_rankscore` | phastCons 470-way rank score |
| `phastCons17way_primate` | phastCons 17-way primate |
| `phastCons17way_primate_rankscore` | phastCons 17-way primate rank score |

Per-allele INFO fields (extracted by `reformat_variants.py`):

| Column | Description |
|--------|-------------|
| `AC` | Allele count (allele-specific) |
| `AF` | Allele frequency (allele-specific) |
| `AQ` | Allele quality (allele-specific) |
| `MLEAC` | Maximum likelihood allele count |
| `MLEAF` | Maximum likelihood allele frequency |

#### Annotation overlay columns (from `reformat_variants.py`)

| Column | Description |
|--------|-------------|
| `segDups` | Overlap count with segmental duplications |
| `simpleRepeats` | Overlap count with simple repeats |
| `constraint_z_1kb_mean` | Mean genomic constraint z-score (1kb windows) |
| `constraint_oe_1kb_mean` | Mean genomic constraint observed/expected (1kb windows) |
| `genebayes_obs_lof` | GeneBayes observed LoF count |
| `genebayes_exp_lof` | GeneBayes expected LoF count |
| `genebayes_prior_mean` | GeneBayes prior mean |
| `genebayes_post_mean` | GeneBayes posterior mean |
| `genebayes_post_lower_95` | GeneBayes posterior 95% lower bound |
| `genebayes_post_upper_95` | GeneBayes posterior 95% upper bound |
| `gnomad_lof_pLI_tx` | gnomAD pLI for this specific transcript |
| `gnomad_lof_oe_ci_upper_tx` | gnomAD LOEUF for this specific transcript |
| `gnomad_lof_pLI_max` | gnomAD pLI max across all transcripts for gene |
| `gnomad_lof_oe_ci_upper_min` | gnomAD LOEUF min across all transcripts for gene |
| `gnomad_lof_pLI_canonical` | gnomAD pLI for canonical transcript |
| `gnomad_lof_oe_ci_upper_canonical` | gnomAD LOEUF for canonical transcript |

## Post-Pipeline Region Annotation

> **Superseded.** This section describes an earlier, partial postprocessing layer.
> The canonical postprocessing now lives outside this repo at
> `/expanse/projects/sebat1/s3/data/sebat/rv_postprocessing_v3/` (DuckDB + SLURM):
> `filter_regions → qc_genotype → join_pop_af → join_scores` (dbNSFP) `→
> join_gene_constraint → filter_pass_and_rescue`, then `tier_variants.py` and the
> `count_*.py` scripts. `annotate_regions.py` and friends were moved to
> `archive/scripts/` — see `archive/README.md`. Kept here for reference only.


After the pipeline produces per-chromosome parquet files, you can add region-based annotations (repeat overlaps, blacklist flags, regulatory elements, etc.) using `annotate_regions.py`.

### Setup

1. Define tracks in the manifest (`resources/region_tracks.tsv`)
2. Download the BED files:
   ```bash
   bash scripts/download_region_tracks.sh
   ```

### Usage

```bash
# Annotate with all downloaded tracks
python archive/scripts/annotate_regions.py output/indexed/chr22.parquet \
    -o output/annotated/chr22.parquet \
    --manifest resources/region_tracks.tsv \
    --regions-dir resources/regions/

# Annotate with specific BEDs only
python archive/scripts/annotate_regions.py output/indexed/chr22.parquet \
    -o output/annotated/chr22.parquet \
    --regions-dir resources/regions/ \
    --beds rmsk.bed encodeBlacklist.bed windowmaskerSdust.bed

# Re-run after adding new tracks (skips already-annotated columns)
python archive/scripts/annotate_regions.py output/annotated/chr22.parquet \
    -o output/annotated/chr22.parquet \
    --manifest resources/region_tracks.tsv \
    --regions-dir resources/regions/ \
    --skip-existing
```

Each track adds an overlap count column (e.g., `rmsk`, `encodeBlacklist`, `windowmaskerSdust`). To add a new annotation, add a row to `resources/region_tracks.tsv`, download it, and re-run with `--skip-existing`.

### Available tracks

See `resources/region_tracks.tsv` for the full list. Key tracks:

| Track | Category | Description |
|-------|----------|-------------|
| `genomicSuperDups` | filter | Segmental duplications |
| `simpleRepeat` | filter | Simple/tandem repeats (TRF) |
| `rmsk` | filter | RepeatMasker (LINE, SINE, LTR, DNA transposons) |
| `windowmaskerSdust` | filter | Low-complexity regions |
| `gap` | filter | Assembly gaps (centromeres, telomeres) |
| `encodeBlacklist` | filter | ENCODE unified blacklist (~910 artifact regions) |
| `encodeCcreCombined` | annotate | ENCODE cis-regulatory elements (promoter, enhancer, CTCF) |
| `cpgIslandExt` | annotate | CpG islands |
| `phastConsElements100way` | annotate | Conserved elements (100 vertebrates) |
| `gwasCatalog` | annotate | NHGRI-EBI GWAS Catalog hits |

## Configuration

Edit `nextflow.config` to customize:
- Container paths
- VEP cache/plugin locations
- Python environment path
- SLURM queue settings
- Resource allocations per process

## License

MIT

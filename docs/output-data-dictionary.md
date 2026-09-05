# Targeted rare-variant pipeline output data dictionary

**Dictionary version:** 1.0  
**Pipeline:** `scripts/run_targeted_chromosome.sh`  
**Genome build:** GRCh38  
**Scope:** LoF and optional missense branches, including autosomes and the validated
chrX/chrY policy.

This document describes the logical schema. Parquet physical integer widths and the
presence of optional columns can vary with the source VCZ and enabled branches. Empty
strings in TSV annotation fields mean unavailable/not applicable; Parquet annotation
fields generally use null for unavailable numeric values.

## Recommended analysis products

Synonymous negative-control outputs use the same tier-gene definitions and QC/AF
eligibility stack as LoF, but are named `*.synonymous-*` and remain outside the
primary rare-burden variables. They are intended for calibration and sensitivity
analyses, not as additional LoF events.

| File pattern | Grain | Natural key | Purpose |
|---|---|---|---|
| `11.plof-burden-eligible.parquet` | qualifying allele × carrier sample | `record_id`, `sample_id` | Final carrier-level LoF analysis table after region, genotype, population-AF, and cohort-AF eligibility |
| `11.missense-burden-eligible.parquet` | qualifying allele × carrier sample | `record_id`, `sample_id` | Equivalent final missense table |
| `12.plof-primary-sample-gene.tsv` | sample × gene × LoF tier | `sample_id`, `Gene`, `lof_tier` | Deduplicated gene-level LoF burdens |
| `12.plof-primary-sample-burden.tsv` | sample | `SAMPLE` | Wide per-sample LoF gene and variant counts; preferred LoF-to-PGS integration table |
| `12.*-per-sample-counts.tsv` | sample | `SAMPLE` | Wide per-sample carrier-row counts by tier |
| `12.*-tier-totals.tsv` | tier | `lof_tier` or `miss_tier` | Run-level QC summary, not a participant table |
| `12.plof-sensitivity-sample-*.tsv` | same as corresponding primary table | same as primary | chrX/chrY burdens that are countable but excluded from the primary analysis |

`_SUCCESS` indicates that all enabled branches and output checks completed. Do not
treat a run directory without `_SUCCESS` as a complete dataset.

## Identifiers and genomic coordinates

| Column | Logical type | Definition |
|---|---|---|
| `sample_id` | string | Sample identifier stored in the input VCZ. Unique only within its cohort/callset unless the data owner guarantees otherwise. |
| `SAMPLE` | string | Renamed `sample_id` in wide count/burden TSVs. |
| `record_id` | string | Stable normalized-allele pointer, encoded as `zarr_<variant_index>_<alt_index>`. It links annotations and recovered genotypes within the exact source VCZ. It is not portable across independently generated Zarr stores. |
| `variant_index` | integer | Zero-based variant-record row in the source Zarr store. |
| `alt_index` | integer | One-based ALT allele index within that source variant record. |
| `sample_index` | integer | Zero-based sample row in the source Zarr store. |
| `#CHROM` | string | GRCh38 chromosome with `chr` prefix. Preferred chromosome field for downstream joins. |
| `CHROM` | string | Annotation chromosome, commonly without `chr`; retained for compatibility. |
| `POS` | integer | GRCh38 VCF position, one-based. |
| `REF` | string | Normalized reference allele. |
| `ALT` | string | Normalized alternate allele represented by this row. |

For cross-pipeline variant joins, use the normalized composite key
`genome_build, #CHROM, POS, REF, ALT`; do not use `record_id`. For participant joins,
construct `cohort_id, sample_id` before combining cohorts.

## Carrier-level columns

These columns originate in `07.*.carriers.parquet` and are retained through the
numbered filtering and annotation stages unless stated otherwise.

| Column | Logical type | Definition |
|---|---|---|
| `Gene` | string | Ensembl gene ID from the VEP-picked transcript; may include a version in early annotation output. |
| `SYMBOL` | string | Gene symbol associated with the picked transcript. |
| `Feature` | string | Selected Ensembl transcript ID (present in allele annotation tables, not guaranteed in carrier tables). |
| `Consequence` | string | Selected VEP/FastVEP consequence term(s), `&`-delimited when multiple apply. |
| `LoF` | string | Standalone LOFTEE classification, normally `HC`, `LC`, or empty/not applicable. |
| `LoF_filter` | string | LOFTEE rejection/filter reason; `.` means none. Present in allele annotation TSVs; not currently copied to carrier Parquet. |
| `LoF_flags` | string | LOFTEE warning flags; `.` means none. Present in allele annotation TSVs; not currently copied to carrier Parquet. |
| `LoF_info` | string | Supporting LOFTEE details; `.` means none. Present in allele annotation TSVs; not currently copied to carrier Parquet. |
| `lof_tier` | categorical string | `lof_t1`: HC pLoF with GeneBayes `post_mean >= 0.18`; `lof_t2`: HC pLoF with `0.03 <= post_mean < 0.18`. There is no current `lof_t3`. |
| `genebayes_post_mean` | float | GeneBayes posterior mean LoF intolerance used for LoF tiering. |
| `miss_tier` | categorical string | `miss_t1` through `miss_t4`, where `miss_t1` satisfies four of four score flags and `miss_t4` satisfies one of four. Empty for LoF rows. |
| `miss_n_flag` | integer | Number of qualifying missense score thresholds, 1–4. |
| `candidate_genes` | string | Comma-delimited Ensembl gene IDs for which the normalized allele passed the dbNSFP candidate rules. |
| `ClinPred_rankscore` | float/null | Maximum available dbNSFP ClinPred percentile rank score. |
| `AlphaMissense_rankscore` | float/null | Maximum available dbNSFP AlphaMissense percentile rank score. |
| `popEVE_converted_rankscore` | float/null | Maximum available dbNSFP popEVE converted percentile rank score. |
| `MPC_rankscore` | float/null | Maximum available dbNSFP MPC percentile rank score; not restricted to MANE transcript coverage. |
| `FILTER` | string | Original source-VCF site FILTER labels joined by `;`; `.` means no label. Whether PASS is required is cohort configuration, not implicit here. |
| `GT` | string | Source genotype rendered with `/` or `|`; fixed-width stores may render haploid calls with a masked slot. Do not infer ploidy from this string. |
| `genotype` | string | Compatibility alias of `GT`. |
| `AD` | string | Reference and selected-ALT depths as `ref,alt`; empty if unavailable. |
| `GQ` | integer/null | Genotype quality from the source callset. |
| `DP` | integer/null | Sample-level read depth from the source callset. |
| `called_ploidy` | integer | Count of unmasked, called genotype slots; authoritative for haploid/diploid interpretation. |
| `alt_dosage` | integer | Count of called alleles equal to this row's `alt_index`. |
| `is_heterozygous` | boolean | `called_ploidy == 2` and `alt_dosage == 1`. |
| `is_homozygous_alt` | boolean | `called_ploidy == 2` and `alt_dosage == 2`. |
| `is_hemizygous_alt` | boolean | `called_ploidy == 1` and `alt_dosage == 1`. |
| `GT`, `genotype` | string | Genotype recoded relative to the carrier row's single `ALT` (`0/1`, `1/1`, and haploid equivalents). |
| `source_GT` | string | Original VCZ genotype encoding before biallelic decomposition; this may contain allele indices such as `0/2`. |
| `AD` | string | Reference and selected-ALT depths for this biallelic carrier row. |

The genotype-QC stage computes allele balance from `AD` but does not retain a separate
`AB` column. Its thresholds are recorded in the run's postprocessing configuration.

## Population- and cohort-frequency columns

`10.*-population-af-annotated.parquet` adds the following without removing rows:

- `gnomAD4.1_joint_AF`, `gnomAD4.1_joint_POPMAX_AF`,
  `gnomAD4.1_joint_nhomalt`, `gnomAD4.1_joint_flag`
- `AllofUs_ALL_AF`, `AllofUs_POPMAX_AF`, `1000Gp3_AF`, `ALFA_Total_AF`
- `RegeneronME_ALL_AF`, `TOPMed_frz8_AC`, `TOPMed_frz8_AF`, `TOPMed_frz8_AN`
- `dbNSFP_POPMAX_AC`, `dbNSFP_POPMAX_AF`

`10.*-population-af-eligible.parquet` is the subset passing the configured population
rule. The current default is `gnomAD4.1_joint_AF < 0.01`, with missing population AF
passing. POPMAX and all other fields are retained as annotations, not substituted for
the joint-AF rule.

`11.*-cohort-af-annotated.parquet` then adds:

| Column | Logical type | Definition |
|---|---|---|
| `cohort_ac` | integer | Alternate allele count in the cohort denominator selected for this chromosome. |
| `cohort_an` | integer | Called allele number in that denominator. |
| `cohort_af` | float | `cohort_ac / cohort_an`. |
| `cohort_carrier_count` | integer | Number of carrier samples in the selected denominator. |
| `cohort_hom_alt_count` | integer | Number of diploid homozygous-alt samples in the selected denominator. |
| `cohort_af_eligible` | boolean | Whether `cohort_af` is below the configured cohort threshold. |

`11.*-burden-eligible.parquet` retains rows with `cohort_af_eligible = true`. The
threshold is cohort-specific runtime provenance and should be stored with integrated
analysis data rather than inferred from the filename.

## Genotype-summary Parquet

`07.*.genotype-summary.parquet` has one row per normalized qualifying allele, including
alleles with no carriers.

| Column | Definition |
|---|---|
| `record_id`, `Gene`, `SYMBOL`, tier | Allele identity and branch annotation (`lof_tier`; missense currently has `miss_tier` only in carrier rows) |
| `genotype_ac`, `genotype_an`, `genotype_af` | AC, AN, and AF across all source samples |
| `carrier_count` | Samples with ALT dosage greater than zero |
| `hom_alt_count` | Samples with ALT dosage two |
| `primary_genotype_ac`, `primary_genotype_an`, `primary_genotype_af` | chrX/chrY primary-policy denominator values; absent on autosomes |
| `primary_carrier_count`, `primary_hom_alt_count` | chrX/chrY primary-policy carrier summaries; absent on autosomes |

## chrX/chrY policy columns

The sex-chromosome carrier tables add these columns. Rows remain reportable even when
they are excluded from primary burden counts.

| Column | Definition |
|---|---|
| `inferred_karyotype` | Cohort QC category such as `XX-like`, `XY-like`, or `ambiguous`; it is evidence-based and not a clinical diagnosis. |
| `copy_number_evidence_pattern` | Compact description of the X/Y evidence underlying that category. |
| `sex_chromosome_region` | `PAR` or `nonPAR`. |
| `par_duplicate_excluded` | True for the noncanonical duplicate representation of a PAR call. |
| `burden_count_available` | Whether a burden can be reported for this row at all. |
| `primary_analysis_eligible` | Whether the row enters the primary burden analysis. |
| `frequency_denominator_eligible` | Whether the sample contributes to the primary AC/AN denominator. |
| `sex_chromosome_policy_reason` | Machine-readable reason for the policy decision. |
| `sensitivity_analysis_group` | `primary`, `ambiguous_karyotype`, or a QC-only exclusion category. |

## Sample-level TSVs

### `12.plof-primary-sample-gene.tsv`

| Column | Definition |
|---|---|
| `sample_id` | Carrier sample ID |
| `Gene`, `SYMBOL` | Gene identifiers |
| `lof_tier` | `lof_t1` or `lof_t2` |
| `n_variants` | Distinct qualifying normalized alleles carried by this sample in this gene/tier |
| `record_ids` | Comma-delimited contributing `record_id` values for auditability |

### `12.plof-primary-sample-burden.tsv`

| Column | Definition |
|---|---|
| `SAMPLE` | Sample ID |
| `lof_t1_genes`, `lof_t2_genes` | Distinct genes with at least one qualifying carried allele in the tier |
| `lof_t1_variants`, `lof_t2_variants` | Distinct qualifying carried alleles in the tier, summed across genes |
| `any_tier_genes`, `any_tier_variants` | Combined LoF tier totals |

Sensitivity LoF files have the same schema but include countable chrX/chrY rows that
are not primary-analysis eligible.

### `12.*-per-sample-counts.tsv`

`SAMPLE` is followed by dynamically generated tier columns (`lof_t1`, `lof_t2`, or
`miss_t1`–`miss_t4`) and `any_group`. Values count carrier rows after all configured
eligibility rules. These are variant-carrier counts, not distinct-gene counts.

### `12.*-tier-totals.tsv`

Contains the tier column plus `n_carrier_rows` and `n_samples_with_any`. This is a
run-level QC summary and should not be joined to a participant-level PGS table.

## Checkpoints and intermediate outputs

| Prefix | Contents |
|---|---|
| `01` | Target alleles selected from the Zarr store |
| `02` | Sites-only VCF carrying Zarr allele pointers |
| `03` | Reference-normalized sites-only VCF |
| `04` | FastVEP rows after VEP-compatible transcript picking |
| `05` | Standalone LOFTEE annotations |
| `06` | GeneBayes/LoF or dbNSFP/missense tier assignment |
| `07` | Recovered carriers and allele-level genotype summaries |
| `08` | Carrier rows after problematic-region filtering |
| `09` | Carrier rows after per-genotype QC |
| `10` | Population-AF annotated and population-eligible tables |
| `11` | Cohort-AF annotated and final burden-eligible tables |
| `12` | Per-sample and run-level burden summaries |

The `*.raw.parquet` sex-chromosome carrier files precede sex-policy annotation and are
checkpoints, not recommended analysis products. `workflow.log` records checkpoint
timings and executed paths but is not a structured scientific output.

## Gene-set burden table

`scripts/build_gene_set_burdens.py` creates one row per participant and gene set:

| Column | Definition |
|---|---|
| `FID`, `IID` | Participant identifiers from the supplied sample file |
| `gene_set_id` | Versioned identifier from `gene_set_membership.tsv` |
| `plof` | Distinct eligible HC-pLoF alleles carried in any member gene |
| `miss_t1`–`miss_t4` | Distinct eligible missense alleles carried in member genes, by classifier tier |

These are separate from the genome-wide GeneBayes burdens. In particular, `plof`
does not require `lof_t1` or `lof_t2`; its input must contain eligible HC-pLoF
carriers from all genes represented by the selected gene sets. The script writes
explicit zeros for callset participants without a qualifying event, deduplicates
transcript/repeated carrier rows, and permits an allele to contribute to multiple
overlapping gene sets. It does not calculate GeneBayes-tier × gene-set intersections.

Example:

```bash
uv run scripts/build_gene_set_burdens.py \
  --samples cohort.psam \
  --gene-sets resources/gene-sets/processed/2026-08-29/gene_set_membership.tsv \
  --plof chr*/11.all-gene-plof-burden-eligible.parquet \
  --missense chr*/11.missense-burden-eligible.parquet \
  --output gene_set_burdens.tsv
```

## PGS integration contract

The proposed concrete modeling variables are listed individually in
[`resources/integrated-analysis-variables.tsv`](../resources/integrated-analysis-variables.tsv).
Its primary rare-variant predictors are exactly `lof_t1`, `lof_t2`, and `miss_t1`
through `miss_t4`. `FID` is included as a nullable identifier until a validated
pedigree/sample mapping supplies it; it must not be inferred from the spelling of
`IID`.

The supplemental chromosome-partition schema is
[`resources/rare-burden-chromosome-strata-variables.tsv`](../resources/rare-burden-chromosome-strata-variables.tsv).
It separates `autosomal`, `sex_chromosome_primary`, and
`sex_chromosome_sensitivity` rows per participant. Primary-eligible sex-chromosome
counts contribute to the six main totals; sensitivity-only counts remain separate.

Create one participant table with the key `(cohort_id, sample_id)` and left-join both
PGS and rare-variant summaries to a cohort sample manifest. A left join to the manifest
is essential: absence from a carrier table can mean zero burden, whereas absence due
to incomplete processing or sample-ID mismatch must not silently become zero.

Recommended provenance columns alongside the joined burdens are:

- `genome_build` (`GRCh38` here), cohort/callset ID, and pipeline commit/image digest;
- LoF and missense tier-definition version;
- population-AF and cohort-AF thresholds;
- reference/dbNSFP/GeneBayes/FastVEP/LOFTEE resource versions;
- chromosome completion status and whether chrX/chrY primary or sensitivity burdens
were used;
- PGS identifier/version, scoring genome build, ancestry model, and standardization
  procedure.

Before integration, standardize `SAMPLE` to `sample_id`, assert uniqueness of
`(cohort_id, sample_id)` in each wide table, and explicitly fill burden columns with
zero only for samples proven present in a successfully completed callset.

### Validated G2MH PGS crosswalk

The validated PGS analysis products are regenerated under `<outdir>/07_analysis/`:

- `analysis_dataset.tsv`: one row per PGS participant;
- `analysis_dataset_dictionary.tsv`: machine-readable column dictionary with fields
  `variable`, `data_type`, `nullable`, `description`, and `source`.

For run `g2mh-integrated-analysis-v1`, `analysis_dataset.tsv` has 1,043 rows, 55
columns, 1,043 unique nonempty `IID` values, and no duplicate IDs. The join is:

| PGS | Rare variant | Integrated canonical name |
|---|---|---|
| `IID` | `SAMPLE` in wide TSVs; `sample_id` in long tables | `sample_id` |
| implicit G2MH run context | explicit integration provenance | `cohort_id = g2mh` |

The PGS columns comprise ancestry assignment/probabilities, ten global PCs,
within-ancestry PCA status and up to ten within-ancestry PCs, and twenty `PGS_*`
scores. The scores are centered PLINK `SCORE1_AVG` values; they are neither
residualized nor standardized, so any analysis-specific transformation must be named
and recorded rather than silently applied during the join.

Observed IDs align exactly when present: all 130 chr1 LoF-burden samples and all 26
chr22 LoF-burden samples were found in the PGS dataset. Carrier-only files are sparse
by design. For example, chr1 missense counts contained 1,048 carrier samples, of whom
1,026 were in the 1,043-person PGS analysis set. Therefore the PGS dataset—not a
carrier file—should define the integrated analysis universe when the scientific goal
is analysis among participants with PGS values. Samples outside that universe must be
reported as excluded, while PGS participants absent from a successfully aggregated
rare-variant carrier table receive zero burden only after callset membership and
chromosome completion are verified.

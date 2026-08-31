# Rare-variant and Zarr integration roadmap

This file records intended extensions that are not required for the current burden-count
workflow. They should remain simple optional branches rather than prerequisites for a cohort
run.

## 1. Preserve and use PSAM pedigree data (initializer complete)

The cohort initializer requires a PSAM, checks `IID` against VCF samples, and retains:

- `FID`: family identifier;
- `IID`: sample identifier;
- `PAT`: paternal sample identifier;
- `MAT`: maternal sample identifier;
- `SEX`: reported sex;

Nonzero `PAT` and `MAT` identifiers must exist in the cohort unless
`--allow-external-parents` is selected. Pedigree-aware outputs should join through this
normalized manifest; annotation and ordinary burden counting do not require complete pedigrees.

## 2. Synonymous-variant output for tiered LoF genes (implemented and chr22-validated)

The negative-control branch selects `synonymous_variant` consequences in genes belonging
to the configured GeneBayes LoF tiers. It supports the analysis tiers currently used for
burden counts (`lof_t1`, GeneBayes posterior mean >= 0.18; and `lof_t2`, > 0 and < 0.18).

The branch reuses all-observed FastVEP output and the existing Zarr genotype extractor.
It does not run LOFTEE on synonymous consequences or mix synonymous counts into LoF burdens.
It emits allele-level annotations, carrier-level genotypes, and per-sample counts under
distinct synonymous output names.

Expected cost: modest once all-observed annotation exists. Annotation is reused; extra work is
a gene/consequence filter, sparse Zarr genotype extraction, and additional Parquet output.

The G2MH chr22 AWS validation completed successfully at commit `175774c`; the full
LoF/missense regression remained unchanged.

## 3. Emit genotypes for relatives of carriers (implemented and chr22-validated)

For every qualifying allele with at least one carrier:

1. use `FID` in the PSAM-derived sample manifest to identify each carrier's family;
2. select all sequenced members of those families;
3. emit `GT`, called ploidy, dosage/carrier status, `GQ`, `DP`, and allele depths for every
   selected family member, including reference and missing calls;
4. retain `FID`, an index-carrier flag, and the carrier sample IDs.

The existing Zarr extractor already reads each selected variant chunk across all samples before
discarding noncarriers. Family expansion therefore requires little additional Zarr I/O or
computation. Output size should grow roughly with the number of sequenced family members per
carrier family, rather than with the entire cohort. It is implemented as an optional,
default-off family-genotype output, not as a replacement for the compact carrier table.

The G2MH chr22 AWS validation completed successfully at commit `3907573`. It added 69 LoF,
1,379 missense, and 5,961 synonymous same-FID rows, while all 32 pinned count checks and ten
canonical hashes passed. The immutable validation image was
`rare-variant-pipeline-targeted@sha256:d2f65272e4cbe2e23313b3744fb7a2082eeeaa1e85d0e3f810c353e2428123ad`.
The FSx run root is
`/fsx/loftee-parity/workflows/g2mh/family-genotypes-chr22-3907573`.

## 4. Reuse cohort Zarr stores for CNV SNP-signal extraction (autosomes complete)

The current CNV preparation path uses `bcftools query` to extract selected SNP-marker rows with
`CHROM`, `POS`, `REF`, `ALT`, `SAMPLE`, `GT`, `AD`, `DP`, and `GQ`. The cohort Zarr stores retain
the corresponding variant, sample, genotype, depth, quality, and allele-depth arrays, so a Zarr
reader can implement:

- interval/marker selection;
- exclusion-region filtering;
- SNP and biallelic-site filtering;
- sample selection;
- long-form or Parquet signal output for downstream LRR/BAF calculation.

This is particularly attractive when the expensive VCF-to-Zarr conversion has already been
performed for the rare-variant workflow. Benchmark a chromosome against the existing bcftools
output and require exact marker/sample/genotype/depth parity before switching the default.

Important boundary: a joint-call Zarr store cannot substitute for sample-gVCF operations that
depend on reference blocks, `INFO/END`, or `FORMAT/MIN_DP` unless those records and fields were
included during conversion. Keep the bcftools/gVCF route for those operations or create a
separate lossless gVCF Zarr representation.

### Proposed pVCF-first integration

The CNV workflow already has a validated grouped-signal input contract and downstream
transpose/assembly/PyPennCNV route. Add one chromosome-parallel producer that converts each
cohort pVCF Zarr store directly into that existing grouped-v2 Parquet contract:

```text
pVCF Zarr -> grouped marker Parquet -> existing transpose/assembly -> PyPennCNV
```

The producer should scan Zarr variant chunks, intersect positions and REF/ALT alleles with the
configured BIM panel, apply the problematic-region mask and `DP`/`GQ` thresholds, resolve
multiallelic `AD` or `LAD+LAA`, and emit aligned `sample_id`, `dp`, and `bac` lists per marker.
It must preserve the current conservative pVCF rule for homozygous-reference calls: retain them
only when the panel ALT has nonzero allele-depth evidence, unless an explicit compatibility
option requests otherwise.

Use the Zarr `sample_id` order to build the numeric sample dictionary once. Do not perform sparse
marker-by-marker reads; scan `variant_position` chunks and read genotype/depth arrays only for
matched rows. Keep the existing grouped-v2 Parquet as the single checkpoint because the CNV
workflow already validates and consumes it.

The autosomal implementation and validation are complete in `gVCFToLRRBAF`. Zarr-derived
grouped rows matched the localized-allele pVCF/bcftools path exactly for the pinned G2MH
comparison, and both PyPennCNV and PyQuantiSNP now consume the assembled signals.

## 5. Add chrX/chrY CNV signals and aneuploidy QC

The validated CNV signal path currently covers chr1--22. PyPennCNV now emits per-chromosome and
genome-wide LRR/BAF summaries, but chrX and chrY cannot appear until sex-chromosome marker
signals are supplied.

Add an optional sex-chromosome signal branch that:

1. extracts chrX and chrY markers from the existing cohort Zarr stores;
2. separates pseudoautosomal (PAR) and non-PAR regions;
3. reports X/Y LRR, BAF, marker count, missingness, and depth relative to the autosomal baseline;
4. preserves reported sex and the independent rare-variant karyotype audit as separate evidence;
5. flags possible X0, XXY, XYY, mosaic, and ambiguous patterns without forcing a diagnosis;
6. calibrates thresholds on SPARK and SSC before enabling production classification.

Estimated effort is 2--4 hours for extraction, summaries, and an initial G2MH validation, plus
approximately one additional day for a calibrated production classifier across SPARK and SSC.
The all-marker BAF mean is affected by pVCF marker ascertainment, so detection should emphasize
relative chromosome LRR, depth, and heterozygous BAF structure rather than BAF mean alone.

## Prioritized remaining work

Keep these optional and avoid adding cohort-specific preparation to the ordinary annotation
path.

1. Run end-to-end production-style validations on SPARK and SSC, including their differing
   callset conventions and PSAM/sample manifests.
2. Finalize chrX/chrY rare-variant counting policy for haploid, diploid, and ambiguous
   karyotypes; retain flagged burdens rather than silently discarding ambiguous samples.
3. Add chrX/chrY CNV marker signals and calibrate aneuploidy QC on SPARK and SSC.
4. Add fine-grained timing around each Zarr genotype extraction, burden aggregation,
   regression validation, and output packaging stage.
5. Profile transcript picking/FastVEP and standalone LOFTEE. Optimize only after profiling;
   in the G2MH chr22 validation they took 124 and 122 seconds respectively out of 381 seconds.
6. Exercise container portability on another Slurm/Singularity or Apptainer environment in
   addition to AWS Batch and Expanse.
7. Validate the participant-manifest-centered analysis assembler on G2MH, then expose it as a
   reusable parent-workflow boundary. PGS, rare variants, and CNVs must remain optional
   components with explicit completeness policies.
8. Define the CNV release step that filters/QCs normalized calls, intersects an approved gene
   annotation, and counts unique affected genes per participant separately for DEL and DUP.
   Do not use raw segment counts as the primary CNV burden.
9. Preserve the optional synonymous and family-genotype branches in the output data dictionary
   and in future combined analysis workflows.

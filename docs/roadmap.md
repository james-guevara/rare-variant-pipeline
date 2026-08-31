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

## 2. Optional synonymous-variant output for tiered LoF genes

Add a negative-control/output branch for `synonymous_variant` consequences in genes belonging
to the configured GeneBayes LoF tiers. Initially support the analysis tiers currently used for
burden counts (`lof_t1`, GeneBayes posterior mean >= 0.18; and `lof_t2`, > 0 and < 0.18).

The branch should reuse all-observed FastVEP output and the existing Zarr genotype extractor.
It should not run LOFTEE on synonymous consequences or mix synonymous counts into LoF burdens.
Emit allele-level annotations, carrier-level genotypes, and optional per-sample counts under
distinct synonymous output names.

Expected cost: modest once all-observed annotation exists. Annotation is reused; extra work is
a gene/consequence filter, sparse Zarr genotype extraction, and additional Parquet output.

## 3. Emit genotypes for relatives of carriers

For every qualifying allele with at least one carrier:

1. use the PSAM-derived pedigree table to identify each carrier's family;
2. select all sequenced members of those families;
3. emit `GT`, called ploidy, dosage/carrier status, `GQ`, `DP`, and allele depths for every
   selected relative, including reference and missing calls;
4. retain `FID`, relationship to the index carrier when derivable, and carrier sample IDs;
5. distinguish relatives absent from the VCF from relatives with a missing genotype.

The existing Zarr extractor already reads each selected variant chunk across all samples before
discarding noncarriers. Family expansion therefore requires little additional Zarr I/O or
computation. Output size should grow roughly with the number of sequenced family members per
carrier family (often approximately 3--5 times the carrier-only row count), rather than with the
entire cohort. Implement this as an optional family-genotype output, not as a replacement for the
compact carrier table.

## 4. Reuse cohort Zarr stores for CNV SNP-signal extraction

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

Initial validation should use chr22 and pinned samples, comparing Zarr-derived grouped rows with
the existing bcftools-derived grouped rows for exact marker, sample, DP, BAC, and derived BAF
parity. Then run the existing signal/call regression gate. Only after parity should the Zarr
producer replace bcftools as the default pVCF extraction engine.

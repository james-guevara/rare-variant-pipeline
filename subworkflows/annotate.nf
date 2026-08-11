// Subworkflow: VCF Processing
// Drop genotypes → VEP → Split VEP → Prepare variants
// Expects normalized (split + left-aligned) input VCFs.

include { BCFTOOLS_SITES } from '../modules/bcftools_sites'
include { VEP_ANNOTATE } from '../modules/vep_annotate'
include { SPLIT_VEP } from '../modules/split_vep'
include { PREPARE_VARIANTS } from '../modules/prepare_variants'

workflow VCF_PROCESSING {
    take:
    input_vcfs       // tuple(chrom, vcf, tbi) — normalized genotyped VCFs

    main:
    BCFTOOLS_SITES(input_vcfs)
    VEP_ANNOTATE(BCFTOOLS_SITES.out.sites)
    SPLIT_VEP(VEP_ANNOTATE.out.vep)
    PREPARE_VARIANTS(SPLIT_VEP.out.variants)

    emit:
    reformatted = PREPARE_VARIANTS.out.reformatted              // tuple(chrom, reformatted.tsv)
    consequential = PREPARE_VARIANTS.out.consequential          // tuple(chrom, consequential.tsv)
    consequential_bed = PREPARE_VARIANTS.out.consequential_bed  // tuple(chrom, consequential.bed)
    input_vcfs = input_vcfs                                     // pass through normed VCFs for family query
}

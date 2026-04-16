// Subworkflow: VCF Processing (Parts 1-4)
// Drop genotypes → VEP → Split VEP → Prepare variants (clean, resolve, filter rare)

include { BCFTOOLS_NORM } from '../modules/bcftools_norm'
include { VEP_ANNOTATE } from '../modules/vep_annotate'
include { SPLIT_VEP } from '../modules/split_vep'
include { PREPARE_VARIANTS } from '../modules/prepare_variants'

workflow VCF_PROCESSING {
    take:
    input_vcfs       // tuple(chrom, vcf, tbi)

    main:
    // Drop genotypes (sites-only VCF for annotation)
    BCFTOOLS_NORM(input_vcfs)

    // VEP annotation
    VEP_ANNOTATE(BCFTOOLS_NORM.out.sites)

    // Split VEP to TSV
    SPLIT_VEP(VEP_ANNOTATE.out.vep)

    // Prepare: clean columns, resolve multi-allelics, filter rare, split by impact
    PREPARE_VARIANTS(SPLIT_VEP.out.variants)

    emit:
    reformatted = PREPARE_VARIANTS.out.reformatted              // tuple(chrom, reformatted.tsv)
    consequential = PREPARE_VARIANTS.out.consequential          // tuple(chrom, consequential.tsv)
    consequential_bed = PREPARE_VARIANTS.out.consequential_bed  // tuple(chrom, consequential.bed)
    input_vcfs = BCFTOOLS_NORM.out.input_vcf                  // tuple(chrom, vcf, tbi) - pass through
}

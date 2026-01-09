// Subworkflow: VCF Processing (Parts 1-4)
// Drop genotypes → VEP → Split VEP → Reformat

include { BCFTOOLS_FILTER } from '../modules/bcftools_filter'
include { VEP_ANNOTATE } from '../modules/vep_annotate'
include { SPLIT_VEP } from '../modules/split_vep'
include { REFORMAT_VARIANTS } from '../modules/reformat_variants'

workflow VCF_PROCESSING {
    take:
    input_vcfs       // tuple(chrom, vcf, tbi)

    main:
    // Drop genotypes (sites-only VCF for annotation)
    BCFTOOLS_FILTER(input_vcfs)

    // VEP annotation
    VEP_ANNOTATE(BCFTOOLS_FILTER.out.sites)

    // Split VEP to TSV
    SPLIT_VEP(VEP_ANNOTATE.out.vep)

    // Reformat and add constraint annotations
    REFORMAT_VARIANTS(SPLIT_VEP.out.variants)

    emit:
    reformatted = REFORMAT_VARIANTS.out.reformatted          // tuple(chrom, reformatted.tsv)
    consequential = REFORMAT_VARIANTS.out.consequential      // tuple(chrom, consequential.tsv)
    consequential_bed = REFORMAT_VARIANTS.out.consequential_bed  // tuple(chrom, consequential.bed)
    input_vcfs = BCFTOOLS_FILTER.out.input_vcf               // tuple(chrom, vcf, tbi) - pass through
}

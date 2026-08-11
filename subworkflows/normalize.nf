// Subworkflow: Normalize input VCFs (split multi-allelics + left-align indels)

include { BCFTOOLS_NORM } from '../modules/bcftools_norm'

workflow NORMALIZE {
    take:
    input_vcfs       // tuple(chrom, vcf, tbi)

    main:
    BCFTOOLS_NORM(input_vcfs)

    emit:
    norm_vcfs = BCFTOOLS_NORM.out.norm   // tuple(chrom, norm.vcf.gz, tbi)
}

// Module: Normalize VCF — split multi-allelics and left-align indels
// Output: with genotypes (for family query) + sites-only (for VEP)

process BCFTOOLS_NORM {
    tag "${chrom}"
    cpus 4
    memory '4 GB'
    time '48h'
    container "${params.bcftools_container}"

    input:
    tuple val(chrom), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("${chrom}.sites.vcf.gz"), path("${chrom}.sites.vcf.gz.csi"), emit: sites
    tuple val(chrom), path("${chrom}.norm.vcf.gz"), path("${chrom}.norm.vcf.gz.tbi"), emit: input_vcf

    script:
    """
    # Split multi-allelics + normalize indels
    bcftools norm -m- -f ${params.reference} -r ${chrom} --threads ${task.cpus} -O z -o ${chrom}.norm.vcf.gz ${vcf}
    bcftools index -t --threads ${task.cpus} ${chrom}.norm.vcf.gz

    # Sites-only for VEP
    bcftools view -G --threads ${task.cpus} -O z -o ${chrom}.sites.vcf.gz ${chrom}.norm.vcf.gz
    bcftools index --threads ${task.cpus} ${chrom}.sites.vcf.gz
    """
}

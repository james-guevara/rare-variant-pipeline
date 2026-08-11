// Module: Drop genotypes → sites-only VCF (for VEP)

process BCFTOOLS_SITES {
    tag "${chrom}"
    cpus 4
    memory '4 GB'
    time '8h'
    container "${params.bcftools_container}"

    input:
    tuple val(chrom), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("${chrom}.sites.vcf.gz"), path("${chrom}.sites.vcf.gz.csi"), emit: sites

    script:
    """
    bcftools view -G -r ${chrom} --threads ${task.cpus} -O z -o ${chrom}.sites.vcf.gz ${vcf}
    bcftools index --threads ${task.cpus} ${chrom}.sites.vcf.gz
    """
}

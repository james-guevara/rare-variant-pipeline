// Module: Drop genotypes from VCF (sites-only)

process BCFTOOLS_FILTER {
    tag "${chrom}"
    cpus 4
    memory '16 GB'
    time '8h'
    container "${params.bcftools_container}"

    input:
    tuple val(chrom), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("${chrom}.sites.vcf.gz"), path("${chrom}.sites.vcf.gz.csi"), emit: sites
    tuple val(chrom), path(vcf), path(tbi), emit: input_vcf  // Pass through for family query

    script:
    """
    bcftools view -G --threads ${task.cpus} -O z -o ${chrom}.sites.vcf.gz ${vcf}
    bcftools index --threads ${task.cpus} ${chrom}.sites.vcf.gz
    """
}

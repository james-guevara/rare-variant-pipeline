// Module: Normalize VCF — split multi-allelics and left-align indels

process BCFTOOLS_NORM {
    tag "${chrom}"
    cpus 4
    memory '4 GB'
    time '48h'
    container "${params.bcftools_container}"

    input:
    tuple val(chrom), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("${chrom}.norm.vcf.gz"), path("${chrom}.norm.vcf.gz.tbi"), emit: norm

    script:
    """
    bcftools norm -m- --force -f ${params.reference} -r ${chrom} --threads ${task.cpus} -O z -o ${chrom}.norm.vcf.gz ${vcf}
    bcftools index -t --threads ${task.cpus} ${chrom}.norm.vcf.gz
    """
}

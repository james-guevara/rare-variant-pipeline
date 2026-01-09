// Module: Split VEP annotations to TSV

process SPLIT_VEP {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '4h'
    container "${params.bcftools_container}"
    containerOptions '--env BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools'

    input:
    tuple val(chrom), path(vep_vcf), path(vep_tbi)

    output:
    tuple val(chrom), path("${chrom}.variants.tsv"), emit: variants

    script:
    """
    bcftools +split-vep -p CSQ -HH -d -s :missense+ \
        -f '%CHROM\t%POS0\t%END\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t%INFO\t%CSQ\n' \
        -A '\t' \
        ${vep_vcf} > ${chrom}.variants.tsv
    """
}

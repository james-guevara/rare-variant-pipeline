// Module: Query family genotypes from VCF for BED regions

process FAMILY_QUERY {
    tag "${chrom}_${chunk_bed.baseName}"
    cpus 1
    memory '16 GB'
    time '4h'
    conda "${params.python_env}"

    input:
    tuple val(chrom), path(chunk_bed), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("${chunk_bed.baseName}.genotypes.tsv"), emit: genotypes

    script:
    """
    python ${params.family_query_script} \\
        --vcf ${vcf} \\
        --ped ${params.ped_file} \\
        --region ${chunk_bed} \\
        --out ${chunk_bed.baseName}.genotypes.tsv
    """
}

// Module: Resolve family genotypes per chunk

process RESOLVE_GENOTYPES {
    tag "${chrom}_${genotypes.baseName}"
    cpus 1
    memory '4 GB'
    time '1h'
    conda "${params.python_env}"

    input:
    tuple val(chrom), path(genotypes)

    output:
    tuple val(chrom), path("${genotypes.baseName}.resolved.tsv"), emit: resolved

    script:
    """
    python ${params.resolve_script} \\
        ${genotypes} \\
        ${genotypes.baseName}.resolved.tsv
    """
}

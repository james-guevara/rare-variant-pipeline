// Module: Resolve family genotypes

process RESOLVE_GENOTYPES {
    tag "${chrom}"
    cpus 1
    memory '16 GB'
    time '2h'
    conda "${params.python_env}"

    input:
    tuple val(chrom), path(genotypes)

    output:
    tuple val(chrom), path("${chrom}.resolved_genotypes.tsv"), emit: resolved

    script:
    """
    python ${params.resolve_script} \\
        ${genotypes} \\
        ${chrom}.resolved_genotypes.tsv
    """
}

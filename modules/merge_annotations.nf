// Module: Merge genotypes with variant annotations

process MERGE_ANNOTATIONS {
    tag "${chrom}"
    cpus 1
    memory '16 GB'
    time '1h'
    conda "${params.python_env}"

    input:
    tuple val(chrom), path(resolved_genotypes), path(consequential_tsv)

    output:
    tuple val(chrom), path("${chrom}.merged.tsv"), emit: merged

    script:
    """
    python ${params.merge_script} \\
        --family ${resolved_genotypes} \\
        --variants ${consequential_tsv} \\
        --out ${chrom}.merged.tsv
    """
}

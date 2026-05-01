// Module: Merge genotypes with variant annotations

process MERGE_ANNOTATIONS {
    tag "${chrom}"
    cpus 1
    memory '16 GB'
    time '1h'
    conda "${params.python_env}"

    input:
    tuple val(chrom), path(carriers_tsv), path(consequential_tsv)

    output:
    tuple val(chrom), path("${chrom}.merged.tsv"), emit: merged

    script:
    """
    python ${params.merge_script} \\
        --family ${carriers_tsv} \\
        --variants ${consequential_tsv} \\
        --out ${chrom}.merged.tsv
    """
}

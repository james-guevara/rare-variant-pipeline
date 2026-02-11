// Module: Apply QC filter to merged annotations

process APPLY_QC_FILTER {
    tag "${chrom}"
    cpus 1
    memory '16 GB'
    time '1h'
    conda "${params.python_env}"

    input:
    tuple val(chrom), path(merged_tsv)

    output:
    tuple val(chrom), path("${chrom}.qc_filtered.tsv"), emit: filtered

    script:
    def mane_arg = params.qc_no_mane_filter ? "--no-mane-filter" : ""
    """
    export POLARS_MAX_THREADS=${task.cpus}
    python ${params.qc_filter_script} \\
        ${merged_tsv} \\
        ${chrom}.qc_filtered.tsv \\
        --min-gq ${params.qc_min_gq} \\
        --min-dp ${params.qc_min_dp} \\
        --max-af ${params.qc_max_af} \\
        --het-ab-min ${params.qc_het_ab_min} \\
        --het-ab-max ${params.qc_het_ab_max} \\
        --hom-ab-min ${params.qc_hom_ab_min} \\
        ${mane_arg}
    """
}

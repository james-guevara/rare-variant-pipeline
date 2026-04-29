// Module: Prepare variants - clean VEP output, extract INFO fields, split by impact

process PREPARE_VARIANTS {
    tag "${chrom}"
    cpus 4
    memory '32 GB'
    time '4h'
    conda "${params.python_env}"

    input:
    tuple val(chrom), path(tsv)

    output:
    tuple val(chrom), path("${chrom}.reformatted.tsv"), emit: reformatted
    tuple val(chrom), path("${chrom}.reformatted.bed"), emit: bed
    tuple val(chrom), path("${chrom}.consequential.tsv"), emit: consequential
    tuple val(chrom), path("${chrom}.consequential.bed"), emit: consequential_bed

    script:
    def mode_args = ""
    if (params.mode == "regulatory" && params.regulatory_beds) {
        mode_args = "--mode regulatory --regulatory-beds ${params.regulatory_beds}"
    } else if (params.mode == "splicing") {
        mode_args = "--mode splicing --spliceai-threshold ${params.spliceai_threshold ?: '0.2'}"
    } else {
        mode_args = "--mode coding"
    }
    def af_arg = "--max-cohort-af ${params.max_cohort_af ?: 1.0}"
    """
    python ${params.prepare_script} \\
        ${tsv} \\
        ${chrom}.reformatted.tsv \\
        --bed ${chrom}.reformatted.bed \\
        --consequential ${chrom}.consequential.tsv \\
        --consequential-bed ${chrom}.consequential.bed \\
        ${mode_args} \\
        ${af_arg}
    """
}

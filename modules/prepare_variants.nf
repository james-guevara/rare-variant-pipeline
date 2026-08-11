// Module: Prepare variants - clean VEP output, extract INFO fields, split by impact

process PREPARE_VARIANTS {
    tag "${chrom}"
    cpus 4
    memory '32 GB'
    time '4h'

    container "${params.python_container}"

    input:
    tuple val(chrom), path(tsv)

    output:
    tuple val(chrom), path("${chrom}.reformatted.tsv"), emit: reformatted
    tuple val(chrom), path("${chrom}.reformatted.bed"), emit: bed
    tuple val(chrom), path("${chrom}.consequential.tsv"), emit: consequential
    tuple val(chrom), path("${chrom}.consequential.bed"), emit: consequential_bed

    script:
    def mode = params.mode ?: "coding"
    def mode_args = "--mode ${mode}"
    if (mode == "regulatory" && params.regulatory_beds) {
        mode_args += " --regulatory-beds ${params.regulatory_beds}"
    } else if (mode == "splicing") {
        mode_args += " --spliceai-threshold ${params.spliceai_threshold ?: '0.2'}"
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

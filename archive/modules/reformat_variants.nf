// Module: Reformat variants and add annotations

process REFORMAT_VARIANTS {
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
    def mode_arg = params.mode ? "--mode ${params.mode}" : "--mode coding"
    def reg_arg = params.mode == "regulatory" && params.regulatory_beds
        ? "--regulatory-beds ${params.regulatory_beds}" : ""
    def filter_arg = params.filter_repeats ? "--filter-repeats" : ""
    """
    python ${params.reformat_script} \\
        ${tsv} \\
        ${chrom}.reformatted.tsv \\
        --bed ${chrom}.reformatted.bed \\
        --consequential ${chrom}.consequential.tsv \\
        --consequential-bed ${chrom}.consequential.bed \\
        --resources-dir ${params.resources_dir} \\
        ${mode_arg} ${reg_arg} ${filter_arg}
    """
}

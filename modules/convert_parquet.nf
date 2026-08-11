// Module: Convert merged TSV.gz to Parquet

process CONVERT_PARQUET {
    tag "${chrom}"
    memory '8 GB'
    time '1h'

    container "${params.python_container}"

    input:
    tuple val(chrom), path(merged_tsv_gz), path(tbi)

    output:
    tuple val(chrom), path("${chrom}.merged.parquet"), emit: parquet

    script:
    """
    python ${params.tsv_to_parquet_script} \
        ${merged_tsv_gz} -o ${chrom}.merged.parquet
    """
}

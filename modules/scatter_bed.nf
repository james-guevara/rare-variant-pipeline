// Module: Scatter BED into chunks for parallel family queries

process SCATTER_BED {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '30m'

    input:
    tuple val(chrom), path(bed)

    output:
    tuple val(chrom), path("${chrom}.chunk_*.bed"), emit: chunks

    script:
    """
    n_regions=\$(wc -l < ${bed})

    if [ \$n_regions -le ${params.regions_per_chunk} ]; then
        cp ${bed} ${chrom}.chunk_00.bed
    else
        split -l ${params.regions_per_chunk} -d --additional-suffix=.bed ${bed} ${chrom}.chunk_
    fi

    echo "Chromosome: ${chrom}, Total regions: \$n_regions, Chunks created: \$(ls ${chrom}.chunk_*.bed | wc -l)"
    """
}

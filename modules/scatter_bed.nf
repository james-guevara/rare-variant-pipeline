// Module: Scatter BED into chunks for parallel family queries

process SCATTER_BED {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '30m'

    container "${params.bcftools_container}"

    input:
    tuple val(chrom), path(bed)

    output:
    tuple val(chrom), path("${chrom}.chunk_*.bed"), env(N_CHUNKS), emit: chunks

    script:
    """
    n_regions=\$(wc -l < ${bed})

    # Chunk with awk rather than `split -d --additional-suffix`. Those are GNU-only
    # flags, and this process runs in a container whose coreutils are BusyBox, which
    # rejects them outright (`split: invalid option -- 'd'`). Chromosomes small enough
    # to take the cp branch below hid this: chrY has 2,823 regions and never called
    # split, so the failure would first appear on a real autosome.
    # Splits on raw lines, exactly as `split -l` did, so chunk contents are unchanged.
    if [ \$n_regions -le ${params.regions_per_chunk} ]; then
        cp ${bed} ${chrom}.chunk_00.bed
    else
        awk -v n=${params.regions_per_chunk} -v pre="${chrom}.chunk_" '
            { if ((NR - 1) % n == 0) { f = sprintf("%s%02d.bed", pre, int((NR - 1) / n)) }
              print > f }
        ' ${bed}
    fi

    N_CHUNKS=\$(ls ${chrom}.chunk_*.bed | wc -l)
    echo "Chromosome: ${chrom}, Total regions: \$n_regions, Chunks created: \$N_CHUNKS"
    """
}

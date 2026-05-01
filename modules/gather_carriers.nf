// Module: Concatenate per-chunk carrier TSVs into a single per-chromosome file.

process GATHER_CARRIERS {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '30m'

    input:
    tuple val(chrom), path(carrier_chunks)

    output:
    tuple val(chrom), path("${chrom}.carriers.tsv"), emit: carriers

    script:
    """
    chunks=\$(ls *.chunk_*.carriers.tsv | sort -V)
    head -1 \$(echo \$chunks | awk '{print \$1}') > ${chrom}.carriers.tsv
    for f in \$chunks; do
        tail -n +2 "\$f"
    done >> ${chrom}.carriers.tsv
    """
}

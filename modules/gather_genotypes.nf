// Module: Gather resolved genotype chunks into single file per chromosome

process GATHER_GENOTYPES {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '30m'

    input:
    tuple val(chrom), path(resolved_chunks)

    output:
    tuple val(chrom), path("${chrom}.resolved_genotypes.tsv"), emit: resolved

    script:
    """
    # Header from first chunk
    head -1 \$(ls *.resolved.tsv | sort -V | head -1) > ${chrom}.resolved_genotypes.tsv

    # Data from all chunks (skip headers)
    for f in \$(ls *.resolved.tsv | sort -V); do
        tail -n +2 "\$f"
    done >> ${chrom}.resolved_genotypes.tsv
    """
}

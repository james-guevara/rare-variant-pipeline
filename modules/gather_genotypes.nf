// Module: Gather genotype chunks into single file per chromosome

process GATHER_GENOTYPES {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '30m'

    input:
    tuple val(chrom), path(genotype_chunks)

    output:
    tuple val(chrom), path("${chrom}.family_genotypes.tsv"), emit: genotypes

    script:
    """
    # Header from first chunk
    head -1 \$(ls *.genotypes.tsv | sort -V | head -1) > ${chrom}.family_genotypes.tsv

    # Data from all chunks (skip headers)
    for f in \$(ls *.genotypes.tsv | sort -V); do
        tail -n +2 "\$f"
    done >> ${chrom}.family_genotypes.tsv
    """
}

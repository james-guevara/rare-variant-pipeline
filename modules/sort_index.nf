// Module: Sort, bgzip, and tabix index

process SORT_INDEX {
    tag "${chrom}"
    cpus 8
    memory '32 GB'
    time '2h'

    input:
    tuple val(chrom), path(merged_tsv)

    output:
    tuple val(chrom), path("${chrom}.merged.tsv.gz"), path("${chrom}.merged.tsv.gz.tbi"), emit: indexed

    script:
    """
    # Sort (header first, then data by position), bgzip
    (head -1 ${merged_tsv}; tail -n +2 ${merged_tsv} | sort --parallel=${task.cpus} -S 4G -k1,1 -k2,2n) | bgzip -@ ${task.cpus} > ${chrom}.merged.tsv.gz

    # Tabix index
    tabix -@ ${task.cpus} -s1 -b2 -e3 ${chrom}.merged.tsv.gz
    """
}

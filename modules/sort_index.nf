// Module: Sort, bgzip, and tabix index

process SORT_INDEX {
    tag "${chrom}"
    cpus 8
    memory '32 GB'
    time '2h'

    container "${params.bcftools_container}"

    input:
    tuple val(chrom), path(merged_tsv)

    output:
    tuple val(chrom), path("${chrom}.merged.tsv.gz"), path("${chrom}.merged.tsv.gz.tbi"), emit: indexed

    script:
    """
    # pipefail matters here: this process runs in a container whose coreutils are
    # BusyBox. BusyBox sort silently ignores GNU's --parallel/-S and emits NOTHING,
    # exiting 0 — so bgzip happily wrote a header-only file and the whole stage
    # "succeeded" while producing an empty parquet downstream. Without pipefail a
    # mid-pipe failure here is invisible.
    set -o pipefail

    # No --parallel / -S: BusyBox sort has neither. Sorting is not the bottleneck in
    # this stage anyway; the merged TSV is far smaller than the VCFs upstream.
    (head -1 ${merged_tsv}; tail -n +2 ${merged_tsv} | sort -k1,1 -k2,2n) | bgzip -@ ${task.cpus} > ${chrom}.merged.tsv.gz

    tabix -@ ${task.cpus} -s1 -b2 -e3 ${chrom}.merged.tsv.gz

    # Guard against a silently-empty result reaching CONVERT_PARQUET.
    n_out=\$(zcat ${chrom}.merged.tsv.gz | wc -l)
    n_in=\$(wc -l < ${merged_tsv})
    if [ "\$n_out" -ne "\$n_in" ]; then
        echo "ERROR: sort/bgzip lost rows: in=\$n_in out=\$n_out" >&2
        exit 1
    fi
    """
}

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir = "out_final"
params.outdir = "out_indexed"
params.chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

process SORT_BGZIP_TABIX {
    tag "${chrom}"
    cpus 8
    memory '32 GB'
    time '2h'

    input:
    val chrom

    output:
    tuple path("*.tsv.gz"), path("*.tsv.gz.tbi"), emit: indexed

    script:
    """
    INPUT="${launchDir}/${params.input_dir}/${chrom}.merged.tsv"
    OUT="${chrom}.merged.tsv.gz"
    
    # Sort (header first, then data by position), bgzip
    (head -1 "\$INPUT"; tail -n +2 "\$INPUT" | sort --parallel=${task.cpus} -S 4G -k1,1 -k2,2n) | bgzip -@ ${task.cpus} > "\$OUT"
    
    # Tabix index
    tabix -@ ${task.cpus} -s1 -b2 -e3 "\$OUT"
    """
}

workflow {
    chroms_ch = Channel.fromList(params.chroms.tokenize(','))
    SORT_BGZIP_TABIX(chroms_ch)
}

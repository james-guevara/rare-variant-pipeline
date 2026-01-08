#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reformat_dir = "out_reformat"
params.scatter_dir = "out_scatter"
params.genotypes_dir = "out_family_query"
params.outdir = "out_gather"
params.trace_file = "trace_part6.txt"
params.chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
params.regions_per_chunk = 1000
params.verify_script = "verify_and_gather.py"

process VERIFY_AND_GATHER {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '30m'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    val chrom
    path reformat_dir
    path scatter_dir
    path genotypes_dir
    path trace_file
    path verify_script

    output:
    path "${chrom}.family_genotypes.tsv"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate python3.12_env_default

    python3 ${verify_script} \\
        --chrom ${chrom} \\
        --reformat-dir ${reformat_dir} \\
        --scatter-dir ${scatter_dir} \\
        --genotypes-dir ${genotypes_dir} \\
        --trace-file ${trace_file} \\
        --output ${chrom}.family_genotypes.tsv \\
        --regions-per-chunk ${params.regions_per_chunk}
    """
}

workflow {
    chrom_list = params.chroms.tokenize(',')
    chroms = Channel.fromList(chrom_list)

    reformat_dir = file(params.reformat_dir)
    scatter_dir = file(params.scatter_dir)
    genotypes_dir = file(params.genotypes_dir)
    trace_file = file(params.trace_file)
    verify_script = file(params.verify_script)

    VERIFY_AND_GATHER(chroms, reformat_dir, scatter_dir, genotypes_dir, trace_file, verify_script)
}

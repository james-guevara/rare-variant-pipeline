#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.genotypes_dir = "out_gather"
params.outdir = "out_resolve"
params.chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
params.resolve_script = "resolve_family_genotypes.py"

process RESOLVE_GENOTYPES {
    tag "${chrom}"
    cpus 1
    memory '16 GB'
    time '2h'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    val chrom
    path genotypes_dir
    path resolve_script

    output:
    path "${chrom}.resolved_genotypes.tsv"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate python3.12_env_default

    python3 ${resolve_script} \\
        ${genotypes_dir}/${chrom}.family_genotypes.tsv \\
        ${chrom}.resolved_genotypes.tsv
    """
}

workflow {
    chrom_list = params.chroms.tokenize(',')
    chroms = Channel.fromList(chrom_list)

    genotypes_dir = file(params.genotypes_dir)
    resolve_script = file(params.resolve_script)

    RESOLVE_GENOTYPES(chroms, genotypes_dir, resolve_script)
}

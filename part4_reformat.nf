#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.indir = "output_variants"
params.outdir = "out_reformat"
params.chroms = "chr21,chr22,chrY"
params.resources = "resources"

process REFORMAT_VARIANTS {
    tag "${chrom}"
    cpus 4
    memory '32 GB'
    time '4h'
    
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(chrom), path(tsv)
    path resources_dir
    path reformat_script
    
    output:
    tuple val(chrom), path("${chrom}.reformatted.tsv"), emit: tsv
    tuple val(chrom), path("${chrom}.reformatted.bed"), emit: bed
    tuple val(chrom), path("${chrom}.consequential.tsv"), emit: conseq_tsv
    tuple val(chrom), path("${chrom}.consequential.bed"), emit: conseq_bed
    
    script:
    """
    micromamba run -n python3.12_env_default python ${reformat_script} \
        ${tsv} \
        ${chrom}.reformatted.tsv \
        --bed ${chrom}.reformatted.bed \
        --consequential ${chrom}.consequential.tsv \
        --consequential-bed ${chrom}.consequential.bed \
        --resources-dir ${resources_dir}
    """
}

workflow {
    // Split chroms into separate channel items
    chrom_list = params.chroms.tokenize(',')
    chroms = Channel.fromList(chrom_list)
    
    // Build input channel from split_vep TSV files
    tsv_files = chroms.map { chrom ->
        def tsv = file("${params.indir}/${chrom}.variants.tsv")
        tuple(chrom, tsv)
    }
    
    resources_dir = file(params.resources)
    reformat_script = file("reformat_variants.py")
    
    REFORMAT_VARIANTS(tsv_files, resources_dir, reformat_script)
}

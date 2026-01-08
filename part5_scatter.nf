#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.indir = "out_reformat"
params.outdir = "out_scatter"
params.chroms = "chr21,chr22,chrY"
params.regions_per_chunk = 1000

process SCATTER_BED {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '30m'
    
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(chrom), path(bed)
    
    output:
    tuple val(chrom), path("${chrom}.chunk_*.bed")
    
    script:
    """
    n_regions=\$(wc -l < ${bed})
    
    if [ \$n_regions -le ${params.regions_per_chunk} ]; then
        cp ${bed} ${chrom}.chunk_00.bed
    else
        split -l ${params.regions_per_chunk} -d --additional-suffix=.bed ${bed} ${chrom}.chunk_
    fi
    
    echo "Chromosome: ${chrom}, Total regions: \$n_regions, Chunks created: \$(ls ${chrom}.chunk_*.bed | wc -l)"
    """
}

workflow {
    chrom_list = params.chroms.tokenize(',')
    chroms = Channel.fromList(chrom_list)
    
    bed_files = chroms.map { chrom ->
        def bed = file("${params.indir}/${chrom}.consequential.bed")
        tuple(chrom, bed)
    }
    
    SCATTER_BED(bed_files)
}

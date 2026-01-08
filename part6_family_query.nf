#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.scatter_dir = "out_scatter"
params.vcf_dir = "vcfs"
params.outdir = "out_family_query"
params.chroms = "chr21,chr22,chrY"
params.ped = "resources/SPARK_iWGS_v1.1.ped"
params.family_query_script = "family_query.py"

process FAMILY_QUERY {
    tag "${chrom}_${chunk_bed.baseName}"
    cpus 1
    memory '16 GB'
    time '4h'
    
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(chrom), path(chunk_bed), path(vcf), path(tbi)
    path ped
    path script
    
    output:
    tuple val(chrom), path("${chunk_bed.baseName}.genotypes.tsv")
    
    script:
    """
    micromamba run -n python3.12_env_default python ${script} \\
        --vcf ${vcf} \\
        --ped ${ped} \\
        --region ${chunk_bed} \\
        --out ${chunk_bed.baseName}.genotypes.tsv
    """
}

workflow {
    chrom_list = params.chroms.tokenize(',')
    
    // Get all chunk BED files using Channel.fromPath
    chunk_beds = Channel.fromPath("${params.scatter_dir}/*.chunk_*.bed")
        .map { bed ->
            def name = bed.name.toString()
            def chrom = name.replaceAll(/\.chunk_.*\.bed$/, '')
            tuple(chrom, bed)
        }
        .filter { chrom, bed -> chrom in chrom_list }
    
    // Create VCF channel
    vcfs = Channel.fromList(chrom_list)
        .map { chrom ->
            def vcf = file("${params.vcf_dir}/wgs_12519_genome.deepvariant.${chrom}.vcf.gz")
            def tbi = file("${params.vcf_dir}/wgs_12519_genome.deepvariant.${chrom}.vcf.gz.tbi")
            tuple(chrom, vcf, tbi)
        }
    
    // Join chunks with VCFs by chromosome
    chunks_with_vcf = chunk_beds.combine(vcfs, by: 0)
    
    ped = file(params.ped)
    script = file(params.family_query_script)
    
    FAMILY_QUERY(chunks_with_vcf, ped, script)
}

#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Part 1: Filter WGS VCFs to exonic regions and drop genotypes

// Parameters
params.chroms = "chr22"  // comma-separated list
params.outdir = "${projectDir}/output_sites"

// Input paths - SPARK iWGS
params.vcf_dir = "${projectDir}/vcfs"
params.vcf_pattern = "wgs_12519_genome.deepvariant.{chrom}.vcf.gz"

// Regions file for filtering (exonic regions)
params.regions_bed = "${projectDir}/resources/GENCODE/gencode.v46.basic.exons.merged.bed"

// Container
params.bcftools_container = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/bcftools:1.22--h3a4d415_1"

process FILTER_AND_DROP_GENOTYPES {
    tag "${chrom}"
    cpus 4
    memory '16 GB'
    time '8h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(chrom), path(vcf), path(tbi)
    path regions_bed

    output:
    tuple val(chrom), path("${chrom}.exonic.sites.vcf.gz"), path("${chrom}.exonic.sites.vcf.gz.csi")

    script:
    """
    singularity exec --bind /expanse/projects/sebat1/ ${params.bcftools_container} \\
        bcftools view -G -R ${regions_bed} --threads ${task.cpus} -O z -o ${chrom}.exonic.sites.vcf.gz ${vcf}

    singularity exec --bind /expanse/projects/sebat1/ ${params.bcftools_container} \\
        bcftools index --threads ${task.cpus} ${chrom}.exonic.sites.vcf.gz
    """
}

workflow {
    chroms_ch = Channel.fromList(params.chroms.tokenize(','))

    input_ch = chroms_ch.map { chrom ->
        def vcf_path = params.vcf_pattern.replace('{chrom}', chrom)
        def vcf_file = file("${params.vcf_dir}/${vcf_path}")
        def tbi_file = file("${params.vcf_dir}/${vcf_path}.tbi")
        return tuple(chrom, vcf_file, tbi_file)
    }

    regions_bed = file(params.regions_bed)

    FILTER_AND_DROP_GENOTYPES(input_ch, regions_bed)
}

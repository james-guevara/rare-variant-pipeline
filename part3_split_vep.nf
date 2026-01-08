#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Part 3: Split VEP annotations to TSV

// Parameters
params.chroms = "chr22"
params.outdir = "${projectDir}/output_variants"

// Input from Part 2
params.vep_dir = "${projectDir}/output_vep"
params.vep_pattern = "{chrom}.vep.vcf.gz"

// Container
params.bcftools_container = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/bcftools:1.22--h3a4d415_1"

process SPLIT_VEP {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '4h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(chrom), path(vep_vcf)

    output:
    tuple val(chrom), path("${chrom}.variants.tsv")

    script:
    """
    singularity exec --bind /expanse/projects/sebat1/ ${params.bcftools_container} \
        bcftools +split-vep -p CSQ -HH -d -s :missense+ \
        -f '%CHROM\t%POS0\t%END\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t%INFO\t%CSQ\n' \
        -A '\t' \
        ${vep_vcf} | \
        sed -E '1s/\\[[0-9]+\\]//g' | \
        sed '1s/CSQ//g' | \
        sed '1s/(null)/INFO/g' \
        > ${chrom}.variants.tsv
    """
}

workflow {
    chroms_ch = Channel.fromList(params.chroms.tokenize(','))

    input_ch = chroms_ch.map { chrom ->
        def vcf_path = params.vep_pattern.replace('{chrom}', chrom)
        def vcf_file = file("${params.vep_dir}/${vcf_path}")
        return tuple(chrom, vcf_file)
    }

    SPLIT_VEP(input_ch)
}

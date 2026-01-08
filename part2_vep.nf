#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Part 2: VEP annotation of filtered sites VCFs

// Parameters
params.chroms = "chr22"
params.outdir = "${projectDir}/output_vep"

// Input from Part 1
params.sites_dir = "${projectDir}/output_sites"
params.sites_pattern = "{chrom}.exonic.sites.vcf.gz"

// Container and resources
params.vep_container = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools"
params.vep_cache = "${projectDir}/VEP_CACHE"
params.vep_plugins = "${projectDir}/VEP_PLUGINS_ALL"
params.loftee_path = "${projectDir}/VEP_PLUGINS/loftee"
params.resources_base = "${projectDir}/resources"
params.reference = "/expanse/projects/sebat1/j3guevar/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

process VEP_ANNOTATE {
    tag "${chrom}"
    cpus 8
    memory '32 GB'
    time '8h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(chrom), path(sites_vcf), path(sites_tbi)

    output:
    tuple val(chrom), path("${chrom}.vep.vcf.gz"), path("${chrom}.vep.vcf.gz.tbi")

    script:
    """
    singularity exec --bind /expanse/projects/sebat1/ ${params.vep_container} vep \\
        --input_file ${sites_vcf} \\
        --format vcf \\
        --output_file ${chrom}.vep.vcf.gz \\
        --vcf \\
        --compress_output bgzip \\
        --minimal --canonical --mane --symbol --protein --pubmed \\
        --nearest symbol --uploaded_allele --numbers \\
        --assembly GRCh38 --cache --dir_cache ${params.vep_cache} --offline \\
        --fasta ${params.reference} \\
        --allele_number --pick_allele --regulatory --biotype --domains \\
        --force_overwrite --fork ${task.cpus} --stats_text \\
        --dir_plugins ${params.vep_plugins} \\
        --plugin dbNSFP,${params.resources_base}/dbNSFP/dbNSFP5.3a_grch38.gz,transcript_match=1,MPC_score,MPC_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,ClinPred_score,ClinPred_rankscore,ClinPred_pred,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,1000Gp3_AC,1000Gp3_AF,AllofUs_ALL_AF,AllofUs_POPMAX_AF,AllofUs_POPMAX_POP,RegeneronME_ALL_AF,gnomAD4.1_joint_flag,gnomAD4.1_joint_AF,gnomAD4.1_joint_nhomalt,gnomAD4.1_joint_POPMAX_AF,gnomAD4.1_joint_POPMAX_nhomalt,ALFA_Total_AF,dbNSFP_POPMAX_AF,dbNSFP_POPMAX_AC,dbNSFP_POPMAX_POP,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,GERP_92_mammals,GERP_92_mammals_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP470way_mammalian,phyloP470way_mammalian_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons470way_mammalian,phastCons470way_mammalian_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore \\
        --plugin LoF,loftee_path:${params.loftee_path},human_ancestor_fa:${params.resources_base}/LOFTEE/human_ancestor.fa.gz,gerp_bigwig:${params.resources_base}/LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:${params.resources_base}/LOFTEE/loftee.sql \\
        --plugin SpliceAI,snv=${params.resources_base}/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=${params.resources_base}/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \\
        --plugin MaxEntScan,${params.resources_base}/MaxEntScan/fordownload,SWA,NCSS

    singularity exec --bind /expanse/projects/sebat1/ ${params.vep_container} \\
        tabix -p vcf ${chrom}.vep.vcf.gz
    """
}

workflow {
    chroms_ch = Channel.fromList(params.chroms.tokenize(','))

    input_ch = chroms_ch.map { chrom ->
        def vcf_path = params.sites_pattern.replace('{chrom}', chrom)
        def vcf_file = file("${params.sites_dir}/${vcf_path}")
        def tbi_file = file("${params.sites_dir}/${vcf_path}.csi")
        return tuple(chrom, vcf_file, tbi_file)
    }

    VEP_ANNOTATE(input_ch)
}

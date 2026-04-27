#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Rare Variant Pipeline - Modular Version
// Run full pipeline or individual subworkflows with -entry

include { NORMALIZE } from './subworkflows/normalize'
include { VCF_PROCESSING } from './subworkflows/vcf_processing'
include { FAMILY_PROCESSING } from './subworkflows/family_processing'
include { MERGE_INDEX } from './subworkflows/merge_index'
include { BCFTOOLS_SITES } from './modules/bcftools_sites'
include { VEP_ANNOTATE } from './modules/vep_annotate'
include { SPLIT_VEP } from './modules/split_vep'

// ============================================================================
// Parameters
// ============================================================================

// Chromosomes to process
params.chroms = "chr22"

// Input VCFs
params.vcf_dir = "/expanse/projects/sebat1/s3/data/sebat/SSC_JG/gatk"
params.vcf_pattern = "{chrom}.masked.vcf.gz"
params.single_vcf = null  // null = per-chrom mode (default); set to VCF path for single-VCF mode
params.normed_vcf_dir = null  // if set, skip NORMALIZE and read per-chrom ${chrom}.norm.vcf.gz from this dir
params.sites_vcf_dir = null   // for RUN_VEP_ONLY: read per-chrom ${chrom}.sites.vcf.gz from this dir
params.vep_vcf_dir = null     // for RUN_SPLIT_VEP_ONLY: read per-chrom ${chrom}.vep.vcf.gz from this dir

// Variant filtering
params.mode = "coding"           // "coding", "regulatory", or "splicing"
params.af_threshold = 0.01       // gnomAD AF threshold for rare filtering (in prepare_variants)
params.cohort_af_threshold = 0.01  // Cohort AF threshold for pre-VEP filtering
params.regulatory_beds = "${projectDir}/resources/regulatory"
params.spliceai_threshold = 0.2

// Family query settings
params.regions_per_chunk = 1000
params.ped_file = "${projectDir}/resources/SPARK_iWGS_v1.1.ped"

// Output
params.outdir = "${projectDir}/output"

// Containers
params.bcftools_container = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/bcftools:1.22--h3a4d415_1"
params.vep_container = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools"

// VEP resources
params.vep_cache = "${projectDir}/VEP_CACHE"
params.vep_plugins = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/VEP_PLUGINS_ALL"
params.loftee_path = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/VEP_PLUGINS/loftee"
params.resources_base = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources"
params.reference = "/expanse/projects/sebat1/j3guevar/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

// Python scripts and resources
params.prepare_script = "${projectDir}/scripts/prepare_variants.py"
params.family_query_script = "${projectDir}/family_query.py"
params.resolve_script = "${projectDir}/resolve_family_genotypes.py"
params.merge_script = "${projectDir}/merge_genotypes_annotations.py"
params.tsv_to_parquet_script = "${projectDir}/scripts/tsv_to_parquet.py"
params.qc_filter_script = "${projectDir}/scripts/qc_filter.py"
params.resources_dir = "${projectDir}/resources"

// QC filter thresholds
params.qc_min_gq = 20
params.qc_min_dp = 10
params.qc_max_af = 0.001
params.qc_het_ab_min = 0.25
params.qc_het_ab_max = 0.75
params.qc_hom_ab_min = 0.9
params.qc_no_mane_filter = false

// ============================================================================
// Helper function to build input channel
// ============================================================================

def findIndex(vcf_path) {
    // Check for .tbi first, then .csi
    def tbi = file("${vcf_path}.tbi")
    if (tbi.exists()) return tbi
    def csi = file("${vcf_path}.csi")
    if (csi.exists()) return csi
    // Return .tbi path as default (will fail at runtime with clear error)
    return tbi
}

def buildNormedChannel(chroms_str, normed_dir) {
    Channel.fromList(chroms_str.tokenize(',')).map { chrom ->
        def vcf_file = file("${normed_dir}/${chrom}.norm.vcf.gz")
        def tbi_file = findIndex("${normed_dir}/${chrom}.norm.vcf.gz")
        return tuple(chrom, vcf_file, tbi_file)
    }
}

def buildSitesChannel(chroms_str, sites_dir) {
    Channel.fromList(chroms_str.tokenize(',')).map { chrom ->
        def vcf_file = file("${sites_dir}/${chrom}.sites.vcf.gz")
        def idx_file = findIndex("${sites_dir}/${chrom}.sites.vcf.gz")
        return tuple(chrom, vcf_file, idx_file)
    }
}

def buildVepChannel(chroms_str, vep_dir) {
    Channel.fromList(chroms_str.tokenize(',')).map { chrom ->
        def vcf_file = file("${vep_dir}/${chrom}.vep.vcf.gz")
        def tbi_file = file("${vep_dir}/${chrom}.vep.vcf.gz.tbi")
        return tuple(chrom, vcf_file, tbi_file)
    }
}

def buildInputChannel(chroms_str, vcf_dir, vcf_pattern) {
    if (params.single_vcf) {
        // Single VCF: all chroms point to the same file, extracted by region downstream
        Channel.fromList(chroms_str.tokenize(',')).map { chrom ->
            def vcf_file = file(params.single_vcf)
            def idx_file = findIndex(params.single_vcf)
            return tuple(chrom, vcf_file, idx_file)
        }
    } else {
        // Per-chrom VCFs (current behavior)
        Channel.fromList(chroms_str.tokenize(',')).map { chrom ->
            def vcf_path = vcf_pattern.replace('{chrom}', chrom)
            def vcf_file = file("${vcf_dir}/${vcf_path}")
            def idx_file = findIndex("${vcf_dir}/${vcf_path}")
            return tuple(chrom, vcf_file, idx_file)
        }
    }
}

// ============================================================================
// Default workflow - Full pipeline
// ============================================================================

workflow {
    // Get normalized VCFs: either run NORMALIZE, or read from normed_vcf_dir if provided
    if (params.normed_vcf_dir) {
        normed = buildNormedChannel(params.chroms, params.normed_vcf_dir)
    } else {
        input_vcfs = buildInputChannel(params.chroms, params.vcf_dir, params.vcf_pattern)
        NORMALIZE(input_vcfs)
        normed = NORMALIZE.out.norm_vcfs
    }

    // VCF Processing: Sites-only → VEP → Split → Prepare
    VCF_PROCESSING(normed)

    // Family Processing: Scatter → Query → Gather → Resolve
    FAMILY_PROCESSING(VCF_PROCESSING.out.consequential_bed, VCF_PROCESSING.out.input_vcfs)

    // Merge and Index: Merge → Sort → Index
    MERGE_INDEX(FAMILY_PROCESSING.out.resolved, VCF_PROCESSING.out.consequential)
}

// ============================================================================
// Entry points for individual subworkflows
// ============================================================================

workflow RUN_NORMALIZE {
    input_vcfs = buildInputChannel(params.chroms, params.vcf_dir, params.vcf_pattern)
    NORMALIZE(input_vcfs)
}

workflow RUN_VCF_PROCESSING {
    // Reads normed VCFs from params.normed_vcf_dir, or from ${params.outdir}/norm by default
    def normed_dir = params.normed_vcf_dir ?: "${params.outdir}/norm"
    normed = buildNormedChannel(params.chroms, normed_dir)
    VCF_PROCESSING(normed)
}

workflow RUN_SITES_ONLY {
    // Drop genotypes from normed VCFs to produce sites-only VCFs.
    // Reads normed VCFs from params.normed_vcf_dir, or from ${params.outdir}/norm by default.
    def normed_dir = params.normed_vcf_dir ?: "${params.outdir}/norm"
    normed = buildNormedChannel(params.chroms, normed_dir)
    BCFTOOLS_SITES(normed)
}

workflow RUN_VEP_ONLY {
    // Run VEP annotation on pre-built sites-only VCFs.
    // Reads sites VCFs from params.sites_vcf_dir, or from ${params.outdir}/sites by default.
    def sites_dir = params.sites_vcf_dir ?: "${params.outdir}/sites"
    sites = buildSitesChannel(params.chroms, sites_dir)
    VEP_ANNOTATE(sites)
}

workflow RUN_SPLIT_VEP_ONLY {
    // Expand the CSQ INFO field of VEP-annotated VCFs into a TSV via bcftools +split-vep.
    // Reads VEP VCFs from params.vep_vcf_dir, or from ${params.outdir}/vep by default.
    def vep_dir = params.vep_vcf_dir ?: "${params.outdir}/vep"
    vep = buildVepChannel(params.chroms, vep_dir)
    SPLIT_VEP(vep)
}

workflow RUN_FAMILY_PROCESSING {
    // Requires outputs from NORMALIZE and VCF_PROCESSING
    chroms = Channel.fromList(params.chroms.tokenize(','))

    consequential_bed = chroms.map { chrom ->
        tuple(chrom, file("${params.outdir}/reformat/${chrom}.consequential.bed"))
    }

    def normed_dir = params.normed_vcf_dir ?: "${params.outdir}/norm"
    input_vcfs = buildNormedChannel(params.chroms, normed_dir)

    FAMILY_PROCESSING(consequential_bed, input_vcfs)
}

workflow RUN_MERGE_INDEX {
    // Requires outputs from both VCF_PROCESSING and FAMILY_PROCESSING
    chroms = Channel.fromList(params.chroms.tokenize(','))

    resolved = chroms.map { chrom ->
        tuple(chrom, file("${params.outdir}/resolve/${chrom}.resolved_genotypes.tsv"))
    }

    consequential = chroms.map { chrom ->
        tuple(chrom, file("${params.outdir}/reformat/${chrom}.consequential.tsv"))
    }

    MERGE_INDEX(resolved, consequential)
}

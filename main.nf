#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Rare Variant Pipeline - Modular Version
// Run full pipeline or individual subworkflows with -entry

include { NORMALIZE } from './subworkflows/normalize'
include { ANNOTATE } from './subworkflows/annotate'
include { CARRIER_EXTRACTION } from './subworkflows/carrier_extraction'
include { EXPORT } from './subworkflows/export'
include { POSTPROCESS } from './subworkflows/postprocess'
include { BCFTOOLS_SITES } from './modules/bcftools_sites'
include { VEP_ANNOTATE } from './modules/vep_annotate'
include { SPLIT_VEP } from './modules/split_vep'
include { PREPARE_VARIANTS } from './modules/prepare_variants'

// ============================================================================
// Parameters
// ============================================================================

// Chromosomes to process
params.chroms = "chr22"

// Input VCFs
params.vcf_dir = "/expanse/projects/sebat1/s3/data/sebat/SSC_JG/gatk"
params.vcf_pattern = "{chrom}.masked.vcf.gz"
params.single_vcf = null  // null = per-chrom mode (default); set to VCF path for single-VCF mode
// Optional BED used for indexed early extraction before normalization. This is
// deliberately a direct shared-filesystem path rather than a process input: very
// large cohort VCFs and the BED live on the shared site filesystem (Expanse or
// FSx) and must not be copied into each task work directory.
params.target_bed = null

// Input dialect — what the source caller actually stores. Consumed by BCFTOOLS_NORM.
params.local_alleles = false  // DRAGEN msVCF: decode LAD/LPL/LAA -> AD/PL BEFORE splitting

// Postprocessing (POSTPROCESS subworkflow). Wraps the DuckDB stage scripts.
params.pp_cohort = null       // key into resources.json "cohorts"; defaults to profile below
params.pp_resources = "${projectDir}/scripts/postprocess/resources.json"
params.count_group_col = "tier"  // column PP_COUNT_CARRIERS stratifies by
params.normed_vcf_dir = null  // if set, skip NORMALIZE and read per-chrom ${chrom}.norm.vcf.gz from this dir
params.sites_vcf_dir = null   // for RUN_VEP_ONLY: read per-chrom ${chrom}.sites.vcf.gz from this dir
params.vep_vcf_dir = null     // for RUN_SPLIT_VEP_ONLY: read per-chrom ${chrom}.vep.vcf.gz from this dir
params.parquet_dir = null    // for RUN_POSTPROCESS: read per-chrom ${chrom}.merged.parquet from this dir
params.split_vep_dir = null   // for RUN_PREPARE_VARIANTS_ONLY: read per-chrom ${chrom}.variants.tsv from this dir

// Variant filtering
params.mode = "coding"           // "coding", "regulatory", or "splicing"
params.af_threshold = 0.01       // gnomAD AF threshold for rare filtering (in prepare_variants)
params.cohort_af_threshold = 0.01  // Cohort AF threshold for pre-VEP filtering
params.max_cohort_af = 1.0       // Cohort AF cap in PREPARE_VARIANTS (1.0 = no filter; set to 0.01 for rare-only)
params.regulatory_beds = "${projectDir}/resources/regulatory"
params.spliceai_threshold = 0.2

// Family query settings
params.regions_per_chunk = 1000
params.ped_file = "${projectDir}/resources/SPARK_iWGS_v1.1.ped"

// Output
params.outdir = "${projectDir}/output"
params.tmpdir = "${params.outdir}/tmp"  // SORT_INDEX scratch; shared fs, never /tmp

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
params.merge_script = "${projectDir}/scripts/merge_genotypes_annotations.py"
params.tsv_to_parquet_script = "${projectDir}/scripts/tsv_to_parquet.py"
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

def buildSplitVepChannel(chroms_str, split_vep_dir) {
    Channel.fromList(chroms_str.tokenize(',')).map { chrom ->
        return tuple(chrom, file("${split_vep_dir}/${chrom}.variants.tsv"))
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
    ANNOTATE(normed)

    // Family Processing: Scatter → Query → Gather → Resolve
    CARRIER_EXTRACTION(ANNOTATE.out.consequential_bed, ANNOTATE.out.input_vcfs)

    // Merge and Index: Merge → Sort → Index
    EXPORT(CARRIER_EXTRACTION.out.carriers, ANNOTATE.out.consequential)
}

// ============================================================================
// Entry points for individual subworkflows
// ============================================================================

workflow RUN_NORMALIZE {
    input_vcfs = buildInputChannel(params.chroms, params.vcf_dir, params.vcf_pattern)
    NORMALIZE(input_vcfs)
}

workflow RUN_ANNOTATE {
    // Reads normed VCFs from params.normed_vcf_dir, or from ${params.outdir}/norm by default
    def normed_dir = params.normed_vcf_dir ?: "${params.outdir}/norm"
    normed = buildNormedChannel(params.chroms, normed_dir)
    ANNOTATE(normed)
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

workflow RUN_PREPARE_VARIANTS_ONLY {
    // Run only the PREPARE_VARIANTS step on existing SPLIT_VEP outputs.
    // Reads per-chrom ${chrom}.variants.tsv from params.split_vep_dir, or
    // from ${params.outdir}/split_vep by default.
    def sv_dir = params.split_vep_dir ?: "${params.outdir}/split_vep"
    sv = buildSplitVepChannel(params.chroms, sv_dir)
    PREPARE_VARIANTS(sv)
}

workflow RUN_CARRIER_EXTRACTION {
    // Requires outputs from NORMALIZE and ANNOTATE
    chroms = Channel.fromList(params.chroms.tokenize(','))

    consequential_bed = chroms.map { chrom ->
        tuple(chrom, file("${params.outdir}/reformat/${chrom}.consequential.bed"))
    }

    def normed_dir = params.normed_vcf_dir ?: "${params.outdir}/norm"
    input_vcfs = buildNormedChannel(params.chroms, normed_dir)

    CARRIER_EXTRACTION(consequential_bed, input_vcfs)
}

workflow RUN_EXPORT {
    // Requires outputs from both ANNOTATE and CARRIER_EXTRACTION
    chroms = Channel.fromList(params.chroms.tokenize(','))

    carriers = chroms.map { chrom ->
        tuple(chrom, file("${params.outdir}/carriers/${chrom}.carriers.tsv"))
    }

    consequential = chroms.map { chrom ->
        tuple(chrom, file("${params.outdir}/reformat/${chrom}.consequential.tsv"))
    }

    EXPORT(carriers, consequential)
}

workflow RUN_POSTPROCESS {
    // Postprocessing on the pipeline's parquet output: region + genotype QC, gnomAD
    // AF, dbNSFP scores, gene constraint, tiering, then per-sample tier counts.
    // Reads per-chrom ${chrom}.merged.parquet from params.parquet_dir, or from
    // ${params.outdir}/parquet by default.
    def pq_dir = params.parquet_dir ?: "${params.outdir}/parquet"
    parquets = Channel.fromList(params.chroms.tokenize(',')).map { chrom ->
        tuple(chrom, file("${pq_dir}/${chrom}.merged.parquet"))
    }
    POSTPROCESS(parquets)
}

#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Rare Variant Pipeline - Modular Version
// Run full pipeline or individual subworkflows with -entry

include { VCF_PROCESSING } from './subworkflows/vcf_processing'
include { FAMILY_PROCESSING } from './subworkflows/family_processing'
include { MERGE_INDEX } from './subworkflows/merge_index'

// ============================================================================
// Parameters
// ============================================================================

// Chromosomes to process
params.chroms = "chr22"

// Input VCFs
params.vcf_dir = "/expanse/projects/sebat1/s3/data/sebat/SSC_JG/gatk"
params.vcf_pattern = "{chrom}.masked.vcf.gz"

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
params.reformat_script = "${projectDir}/reformat_variants.py"
params.family_query_script = "${projectDir}/family_query.py"
params.resolve_script = "${projectDir}/resolve_family_genotypes.py"
params.merge_script = "${projectDir}/merge_genotypes_annotations.py"
params.resources_dir = "${projectDir}/resources"

// ============================================================================
// Helper function to build input channel
// ============================================================================

def buildInputChannel(chroms_str, vcf_dir, vcf_pattern) {
    Channel.fromList(chroms_str.tokenize(',')).map { chrom ->
        def vcf_path = vcf_pattern.replace('{chrom}', chrom)
        def vcf_file = file("${vcf_dir}/${vcf_path}")
        def tbi_file = file("${vcf_dir}/${vcf_path}.tbi")
        return tuple(chrom, vcf_file, tbi_file)
    }
}

// ============================================================================
// Default workflow - Full pipeline
// ============================================================================

workflow {
    // Build input channel
    input_vcfs = buildInputChannel(params.chroms, params.vcf_dir, params.vcf_pattern)

    // VCF Processing: Drop genotypes → VEP → Split → Reformat
    VCF_PROCESSING(input_vcfs)

    // Family Processing: Scatter → Query → Gather → Resolve
    FAMILY_PROCESSING(VCF_PROCESSING.out.consequential_bed, VCF_PROCESSING.out.input_vcfs)

    // Merge and Index: Merge → Sort → Index
    MERGE_INDEX(FAMILY_PROCESSING.out.resolved, VCF_PROCESSING.out.consequential)
}

// ============================================================================
// Entry points for individual subworkflows
// ============================================================================

workflow RUN_VCF_PROCESSING {
    input_vcfs = buildInputChannel(params.chroms, params.vcf_dir, params.vcf_pattern)
    VCF_PROCESSING(input_vcfs)
}

workflow RUN_FAMILY_PROCESSING {
    // Requires outputs from VCF_PROCESSING
    chroms = Channel.fromList(params.chroms.tokenize(','))

    consequential_bed = chroms.map { chrom ->
        tuple(chrom, file("${params.outdir}/reformat/${chrom}.consequential.bed"))
    }

    input_vcfs = chroms.map { chrom ->
        def vcf_path = params.vcf_pattern.replace('{chrom}', chrom)
        tuple(chrom, file("${params.vcf_dir}/${vcf_path}"), file("${params.vcf_dir}/${vcf_path}.tbi"))
    }

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

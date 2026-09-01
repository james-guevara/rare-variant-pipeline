#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.participant_manifest = null
params.pgs_dataset = null
params.pgs_dictionary = null
params.rare_burdens = null
params.cnv_dataset = null
params.cnv_dictionary = null
params.variable_template = "${projectDir}/resources/integrated-analysis-variables.tsv"
params.cohort_id = null
params.targeted_container = null
params.outdir = 'results'
params.missing_pgs_policy = 'allow'
params.missing_rare_policy = 'allow'
params.missing_cnv_policy = 'allow'

include { ANALYSIS_DATASET_WORKFLOW } from './workflows/analysis_dataset'

workflow {
    if (!params.participant_manifest || !params.cohort_id || !params.targeted_container) {
        error 'Required: --participant_manifest, --cohort_id, and --targeted_container'
    }
    if ((params.pgs_dataset ? true : false) != (params.pgs_dictionary ? true : false)) {
        error '--pgs_dataset and --pgs_dictionary must be supplied together'
    }
    if ((params.cnv_dataset ? true : false) != (params.cnv_dictionary ? true : false)) {
        error '--cnv_dataset and --cnv_dictionary must be supplied together'
    }
    if (!params.pgs_dataset && !params.rare_burdens && !params.cnv_dataset) {
        error 'At least one PGS, rare-burden, or CNV dataset is required'
    }

    sentinel = file("${projectDir}/resources/empty-analysis-component.tsv", checkIfExists: true)
    pgs_inputs = params.pgs_dataset \
        ? Channel.value(tuple(true,
            file(params.pgs_dataset, checkIfExists: true),
            file(params.pgs_dictionary, checkIfExists: true))) \
        : Channel.value(tuple(false, sentinel, sentinel))
    rare_input = params.rare_burdens \
        ? Channel.value(tuple(true, file(params.rare_burdens, checkIfExists: true))) \
        : Channel.value(tuple(false, sentinel))
    cnv_inputs = params.cnv_dataset \
        ? Channel.value(tuple(true,
            file(params.cnv_dataset, checkIfExists: true),
            file(params.cnv_dictionary, checkIfExists: true))) \
        : Channel.value(tuple(false, sentinel, sentinel))

    ANALYSIS_DATASET_WORKFLOW(
        Channel.value(file(params.participant_manifest, checkIfExists: true)),
        pgs_inputs,
        rare_input,
        cnv_inputs,
        Channel.value(file(params.variable_template, checkIfExists: true))
    )
}

#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.pgs_dataset = null
params.pgs_dictionary = null
params.rare_burdens = null
params.variable_template = "${projectDir}/resources/integrated-analysis-variables.tsv"
params.cohort_id = null
params.targeted_container = null
params.outdir = 'results'

include { INTEGRATED_ANALYSIS_WORKFLOW } from './workflows/integrated_analysis'

workflow {
    if (!params.pgs_dataset || !params.pgs_dictionary || !params.rare_burdens || !params.cohort_id || !params.targeted_container) {
        error 'Required: --pgs_dataset, --pgs_dictionary, --rare_burdens, --cohort_id, and --targeted_container'
    }
    INTEGRATED_ANALYSIS_WORKFLOW(
        Channel.value(file(params.pgs_dataset, checkIfExists: true)),
        Channel.value(file(params.pgs_dictionary, checkIfExists: true)),
        Channel.value(file(params.rare_burdens, checkIfExists: true)),
        Channel.value(file(params.variable_template, checkIfExists: true))
    )
}

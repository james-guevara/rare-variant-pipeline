#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.run_sheet = null
params.participant_manifest = null
params.expected_chromosomes = null
params.pgs_dataset = null
params.pgs_dictionary = null
params.cnv_dataset = null
params.cnv_dictionary = null
params.variable_template = "${projectDir}/resources/integrated-analysis-variables.tsv"
params.cohort_id = null
params.targeted_container = null
params.lof_only = false
params.run_root = null
params.outdir = 'results'
params.missing_pgs_policy = 'allow'
params.missing_rare_policy = 'error'
params.missing_cnv_policy = 'allow'

include { COHORT_ANALYSIS_WORKFLOW } from './workflows/cohort_analysis'

workflow {
    if (!params.run_sheet || !params.participant_manifest ||
        !params.expected_chromosomes || !params.cohort_id || !params.targeted_container) {
        error 'Required: --run_sheet, --participant_manifest, --expected_chromosomes, --cohort_id, and --targeted_container'
    }
    if (params.run_root) {
        error '--run_root is only valid for a single-chromosome launch; each cohort binding must supply a unique run_root'
    }
    if ((params.pgs_dataset ? true : false) != (params.pgs_dictionary ? true : false)) {
        error '--pgs_dataset and --pgs_dictionary must be supplied together'
    }
    if ((params.cnv_dataset ? true : false) != (params.cnv_dictionary ? true : false)) {
        error '--cnv_dataset and --cnv_dictionary must be supplied together'
    }

    manifest_bindings = Channel.fromPath(params.run_sheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row ->
            if (!row.chromosome || !row.manifest || !row.bindings) {
                error 'Run sheet requires chromosome, manifest, and bindings columns'
            }
            tuple(
                file(row.manifest.toString(), checkIfExists: true),
                file(row.bindings.toString(), checkIfExists: true)
            )
        }

    pgs_inputs = params.pgs_dataset \
        ? Channel.value(tuple(true,
            file(params.pgs_dataset, checkIfExists: true),
            file(params.pgs_dictionary, checkIfExists: true))) \
        : Channel.value(tuple(false,
            file("${projectDir}/resources/empty-analysis-pgs-dataset.tsv", checkIfExists: true),
            file("${projectDir}/resources/empty-analysis-pgs-dictionary.tsv", checkIfExists: true)))

    cnv_inputs = params.cnv_dataset \
        ? Channel.value(tuple(true,
            file(params.cnv_dataset, checkIfExists: true),
            file(params.cnv_dictionary, checkIfExists: true))) \
        : Channel.value(tuple(false,
            file("${projectDir}/resources/empty-analysis-cnv-dataset.tsv", checkIfExists: true),
            file("${projectDir}/resources/empty-analysis-cnv-dictionary.tsv", checkIfExists: true)))

    COHORT_ANALYSIS_WORKFLOW(
        manifest_bindings,
        Channel.value(file(params.participant_manifest, checkIfExists: true)),
        Channel.value(params.expected_chromosomes),
        pgs_inputs,
        cnv_inputs,
        Channel.value(file(params.variable_template, checkIfExists: true))
    )
}

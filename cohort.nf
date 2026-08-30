#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.run_sheet = null
params.sample_manifest = null
params.expected_chromosomes = null
params.targeted_container = null
params.lof_only = false
params.run_root = null
params.outdir = 'results'

include { COHORT_RARE_BURDEN_WORKFLOW } from './workflows/cohort_rare_burden'

workflow {
    if (!params.run_sheet || !params.sample_manifest || !params.expected_chromosomes || !params.targeted_container) {
        error 'Required: --run_sheet, --sample_manifest, --expected_chromosomes, and --targeted_container'
    }
    if (params.run_root) {
        error '--run_root is only valid for a single-chromosome launch; each cohort binding must supply a unique run_root'
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

    COHORT_RARE_BURDEN_WORKFLOW(
        manifest_bindings,
        Channel.value(file(params.sample_manifest, checkIfExists: true)),
        Channel.value(params.expected_chromosomes)
    )
}

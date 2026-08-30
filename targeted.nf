#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.manifest = null
params.bindings = null
params.targeted_container = null
params.preflight_only = false
params.lof_only = false
params.run_root = null

include { TARGETED_MANIFEST_WORKFLOW } from './workflows/targeted_manifest'

workflow {
    if (!params.manifest || !params.bindings || !params.targeted_container) {
        error 'Required: --manifest, --bindings, and --targeted_container'
    }
    manifest_bindings = Channel.of(tuple(file(params.manifest), file(params.bindings)))
    TARGETED_MANIFEST_WORKFLOW(manifest_bindings)
}

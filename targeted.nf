#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.manifest = null
params.bindings = null
params.targeted_container = null
params.preflight_only = false
params.run_root = null

include { TARGETED_MANIFEST_WORKFLOW } from './workflows/targeted_manifest'

workflow {
    TARGETED_MANIFEST_WORKFLOW()
}

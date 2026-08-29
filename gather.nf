#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.sample_manifest = null
params.package_glob = null
params.expected_chromosomes = (1..22).collect { "chr${it}" }.plus(['chrX', 'chrY']).join(',')
params.targeted_container = null
params.outdir = 'results'

include { RARE_BURDEN_GATHER_WORKFLOW } from './workflows/rare_burden_gather'

workflow {
    if (!params.sample_manifest || !params.package_glob || !params.targeted_container) {
        error 'Required: --sample_manifest, --package_glob, and --targeted_container'
    }
    sample_manifest = Channel.value(file(params.sample_manifest))
    packages = Channel.fromPath(params.package_glob, type: 'dir', checkIfExists: true)
    expected = Channel.value(params.expected_chromosomes)
    RARE_BURDEN_GATHER_WORKFLOW(sample_manifest, packages, expected)
}

#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.sample_manifest = null
params.package_glob = null
params.expected_chromosomes = (1..22).collect { "chr${it}" }.plus(['chrX', 'chrY']).join(',')
params.targeted_container = null
params.outdir = 'results'
params.gene_set_membership = "${projectDir}/resources/gene-sets/processed/2026-08-29/gene_set_membership.tsv"

include { RARE_BURDEN_GATHER_WORKFLOW } from './workflows/rare_burden_gather'

workflow {
    if (!params.sample_manifest || !params.package_glob || !params.targeted_container) {
        error 'Required: --sample_manifest, --package_glob, and --targeted_container'
    }
    sample_manifest = Channel.value(file(params.sample_manifest))
    packages = Channel.fromPath(params.package_glob, type: 'dir', checkIfExists: true)
    expected = Channel.value(params.expected_chromosomes)
    gene_sets = Channel.value(file(params.gene_set_membership, checkIfExists: true))
    consolidate = Channel.value(file("${projectDir}/scripts/consolidate_carrier_parquets.py", checkIfExists: true))
    gene_burdens = Channel.value(file("${projectDir}/scripts/build_gene_set_burdens.py", checkIfExists: true))
    RARE_BURDEN_GATHER_WORKFLOW(
        sample_manifest, packages, expected, gene_sets, consolidate, gene_burdens
    )
}

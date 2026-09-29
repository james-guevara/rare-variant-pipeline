#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BLOCK_SITES_CATALOG } from './modules/block_sites_catalog'

params.sites_manifest = null
params.sites_publish_mode = 'copy'

workflow {
    if (!params.sites_manifest) error 'Required: --sites_manifest (unit_id, chromosome, vcf TSV)'
    if (!(params.sites_publish_mode in ['copy', 'link'])) error '--sites_publish_mode must be copy or link'
    def manifest = file(params.sites_manifest, checkIfExists: true)
    // Validate the complete supplied manifest before scheduling any blocks.
    // This reads TSV rows and stats paths, never VCF contents.
    def rows = SitesCatalogManifest.load(manifest)
    BLOCK_SITES_CATALOG(Channel.fromList(rows).map { row ->
        tuple(row.unit_id, row.chromosome, row.vcf, file(row.vcf, checkIfExists: true))
    }, Channel.value(file("${projectDir}/scripts/make_sites_catalog.sh")),
       Channel.value(file("${projectDir}/scripts/check_sites_records.awk")))
}

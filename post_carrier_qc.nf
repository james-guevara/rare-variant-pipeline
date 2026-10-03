nextflow.enable.dsl=2
include { POST_CARRIER_QC } from './modules/post_carrier_qc'

workflow {
    if (!params.qc_manifest) error 'Required: --qc_manifest and --select_units'
    def containerIdentity = null
    if (params.qc_container) {
        if (!params.resource_lock) error 'Required with container: --resource_lock (existing annotation checksum lock)'
        def lock = new groovy.json.JsonSlurper().parseText(file(params.resource_lock, checkIfExists: true).text)
        containerIdentity = lock.containers.loftee
        def sif = file(params.qc_container, checkIfExists: true)
        if (containerIdentity.path != sif.toString() ||
            containerIdentity.sha256 != 'a1722375533ad25de58d01635d3cf36667b4ca909149c59404a8deb6491907cd' ||
            java.nio.file.Files.size(sif) != containerIdentity.bytes ||
            java.nio.file.Files.getLastModifiedTime(sif).toMillis() != containerIdentity.mtime_ms)
            error 'Carrier SIF differs from validated annotation lock; verify container identity'
    }
    def rows = AnnotationResources.select(
        CarrierQcManifest.load(file(params.qc_manifest, checkIfExists: true)), params.select_units)
    def entries = Channel.fromList(rows).map { row ->
        def meta = row + [container: params.qc_container, container_identity: containerIdentity]
        tuple(meta, file(row.carriers), file(row.samples), file(row.source_receipt))
    }
    POST_CARRIER_QC(entries, Channel.value(file("${projectDir}/scripts/qc_filtered_carriers.py")))
}

nextflow.enable.dsl=2
include { FINAL_RARITY } from './modules/final_rarity'

workflow {
    if (!params.rarity_manifest) error 'Required: --rarity_manifest and --select_units'
    def containerIdentity = null
    if (params.rarity_container) {
        if (!params.resource_lock) error 'Required with container: --resource_lock (existing annotation checksum lock)'
        def lock = new groovy.json.JsonSlurper().parseText(file(params.resource_lock, checkIfExists: true).text)
        containerIdentity = lock.containers.loftee
        def sif = file(params.rarity_container, checkIfExists: true)
        if (containerIdentity.path != sif.toString() ||
            containerIdentity.sha256 != 'a1722375533ad25de58d01635d3cf36667b4ca909149c59404a8deb6491907cd' ||
            java.nio.file.Files.size(sif) != containerIdentity.bytes ||
            java.nio.file.Files.getLastModifiedTime(sif).toMillis() != containerIdentity.mtime_ms)
            error 'Carrier SIF differs from validated annotation lock; verify container identity'
    }
    def rows = AnnotationResources.select(
        FinalRarityManifest.load(file(params.rarity_manifest, checkIfExists: true)), params.select_units)
    def entries = Channel.fromList(rows).map { row ->
        def meta = row + [container: params.rarity_container, container_identity: containerIdentity]
        tuple(meta, file(row.carriers), file(row.samples), file(row.source_receipt), file(row.candidates), file(row.frequencies))
    }
    FINAL_RARITY(entries, Channel.value(file("${projectDir}/scripts/filter_final_rarity.py")), Channel.value(file("${projectDir}/scripts/qc_filtered_carriers.py")))
}

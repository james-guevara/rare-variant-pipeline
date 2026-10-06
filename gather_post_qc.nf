nextflow.enable.dsl=2
include { GATHER_POST_QC } from './modules/gather_post_qc'
workflow {
    if (!params.gather_manifest) error 'Required: --gather_manifest --select_units --gather_id'
    if (!(params.gather_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Required: safe --gather_id identifying this selection'
    def containerIdentity = null
    if (params.gather_container) {
        if (!params.resource_lock) error 'Required with container: --resource_lock (existing annotation checksum lock)'
        def lock = new groovy.json.JsonSlurper().parseText(file(params.resource_lock, checkIfExists: true).text)
        containerIdentity = lock.containers.loftee
        def sif = file(params.gather_container, checkIfExists: true)
        if (containerIdentity.path != sif.toString() ||
            containerIdentity.sha256 != 'a1722375533ad25de58d01635d3cf36667b4ca909149c59404a8deb6491907cd' ||
            java.nio.file.Files.size(sif) != containerIdentity.bytes ||
            java.nio.file.Files.getLastModifiedTime(sif).toMillis() != containerIdentity.mtime_ms)
            error 'Carrier SIF differs from validated annotation lock; verify container identity'
    }
    def allRows = CarrierQcManifest.load(file(params.gather_manifest,checkIfExists:true))
    def selected = AnnotationResources.select(allRows,params.select_units)
    def groups = selected.groupBy { it.chromosome }.collect { chrom, units ->
        units = units.sort { it.unit_id }
        def meta = [chromosome:chrom,gather_id:params.gather_id,units:units,
                    manifest_units:allRows.findAll{it.chromosome==chrom}.collect{it.unit_id}.sort(),
                    container:params.gather_container,container_identity:containerIdentity]
        tuple(meta,units.collect{file(it.carriers)},units.collect{file(it.samples)},units.collect{file(it.source_receipt)})
    }
    GATHER_POST_QC(Channel.fromList(groups),Channel.value(file("${projectDir}/scripts/gather_post_qc.py")),Channel.value(file("${projectDir}/scripts/qc_filtered_carriers.py")))
}

nextflow.enable.dsl=2
include { FILTERED_CARRIERS } from './modules/filtered_carriers'

workflow {
    if (!params.carrier_manifest) error 'Required: --carrier_manifest and --select_units'
    if (!(params.compute_frequencies in [true,false,'true','false'])) error 'compute_frequencies must be true or false'
    def frequencyEnabled=params.compute_frequencies.toString()=='true'
    if (frequencyEnabled && !params.psam) error 'Corrected frequencies require --psam'
    def containerIdentity = null
    if (params.carrier_container) {
        if (!params.resource_lock) error 'Required with container: --resource_lock (existing annotation checksum lock)'
        def lock = new groovy.json.JsonSlurper().parseText(file(params.resource_lock, checkIfExists: true).text)
        containerIdentity = lock.containers.loftee
        def sif = file(params.carrier_container, checkIfExists: true)
        if (containerIdentity.path != sif.toString() ||
            containerIdentity.sha256 != 'a1722375533ad25de58d01635d3cf36667b4ca909149c59404a8deb6491907cd' ||
            java.nio.file.Files.size(sif) != containerIdentity.bytes ||
            java.nio.file.Files.getLastModifiedTime(sif).toMillis() != containerIdentity.mtime_ms)
            error 'Carrier SIF differs from validated annotation lock; verify container identity'
    }
    def rows = AnnotationResources.select(
        FilteredCarrierManifest.load(file(params.carrier_manifest, checkIfExists: true)), params.select_units)
    if (params.sex_chromosome_policy && params.sex_chromosome_policy != 'grch38_x_only_par')
        error 'Unsupported sex chromosome frequency policy'
    if (frequencyEnabled && rows.any { it.chromosome in ['chrX','chrY'] } && params.sex_chromosome_policy != 'grch38_x_only_par')
        error 'Corrected X/Y frequencies require --sex_chromosome_policy grch38_x_only_par'
    def entries = Channel.fromList(rows).map { row ->
        def meta = row + [container: params.carrier_container, container_identity: containerIdentity, expected_hc: params.expected_hc, expected_missense: params.expected_missense, psam: params.psam, compute_frequencies: frequencyEnabled, sex_chromosome_policy: params.sex_chromosome_policy]
        tuple(meta, file(row.missense), file(row.lof_hc), file(row.vcf), file(row.index))
    }
    FILTERED_CARRIERS(entries, Channel.value(file("${projectDir}/scripts/extract_filtered_carriers.py")),
        Channel.value(file("${projectDir}/scripts/extract_exact_carriers.py")),
        Channel.value(file("${projectDir}/scripts/carrier_psam.py")),
        Channel.value(file("${projectDir}/scripts/carrier_frequencies.py")),
        Channel.value(params.psam ? file(params.psam,checkIfExists:true) : []))
}

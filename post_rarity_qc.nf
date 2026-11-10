nextflow.enable.dsl=2
include { POST_RARITY_QC } from './modules/post_rarity_qc'

workflow {
    if (!params.post_rarity_qc_manifest) error 'Required: --post_rarity_qc_manifest and --select_units'
    if (!(params.site_filter_policy in ['pass','pass_or_missing'])) error 'Unsupported site_filter_policy'
    if (!params.psam) error 'Required: --psam'
    if (params.sex_chromosome_policy && params.sex_chromosome_policy!='grch38_x_only_par') error 'Unsupported sex chromosome policy'
    def containerIdentity = null
    if (params.post_rarity_qc_container) {
        if (!params.resource_lock) error 'Required with container: --resource_lock (existing annotation checksum lock)'
        def lock = new groovy.json.JsonSlurper().parseText(file(params.resource_lock, checkIfExists: true).text)
        containerIdentity = lock.containers.loftee
        def sif = file(params.post_rarity_qc_container, checkIfExists: true)
        if (containerIdentity.path != sif.toString() ||
            containerIdentity.sha256 != 'a1722375533ad25de58d01635d3cf36667b4ca909149c59404a8deb6491907cd' ||
            java.nio.file.Files.size(sif) != containerIdentity.bytes ||
            java.nio.file.Files.getLastModifiedTime(sif).toMillis() != containerIdentity.mtime_ms)
            error 'Carrier SIF differs from validated annotation lock; verify container identity'
    }
    def rows = AnnotationResources.select(
        CarrierQcManifest.load(file(params.post_rarity_qc_manifest, checkIfExists: true)), params.select_units)
    if (rows.any { it.chromosome in ['chrX','chrY'] } && params.sex_chromosome_policy!='grch38_x_only_par') error 'X/Y QC requires --sex_chromosome_policy grch38_x_only_par'
    def entries = Channel.fromList(rows).map { row ->
        def meta = row + [container: params.post_rarity_qc_container, container_identity: containerIdentity, psam: params.psam, sex_chromosome_policy: params.sex_chromosome_policy, site_filter_policy: params.site_filter_policy]
        tuple(meta, file(row.carriers), file(row.samples), file(row.source_receipt))
    }
    POST_RARITY_QC(entries, Channel.value(file("${projectDir}/scripts/qc_post_rarity.py")),
        Channel.value(['qc_filtered_carriers.py','filter_final_rarity.py','carrier_frequencies.py','carrier_psam.py','extract_exact_carriers.py'].collect { file("${projectDir}/scripts/${it}") }),
        Channel.value(file(params.psam,checkIfExists:true)))
}

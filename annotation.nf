nextflow.enable.dsl = 2
include { SITES_ANNOTATION } from './subworkflows/sites_annotation'

workflow {
    if (!params.sites_manifest || !params.resource_lock)
        error 'Required: --sites_manifest, --resource_lock, --select_units (IDs or all)'
    def rows = AnnotationResources.select(
        SitesCatalogManifest.load(file(params.sites_manifest, checkIfExists: true)), params.select_units)
    def resources = AnnotationResources.load(file(params.resource_lock, checkIfExists: true), rows)
    if (params.annotation_root != resources.annotation_root || params.loftee_root != resources.loftee_root)
        error 'Configured annotation_root/loftee_root must match resource lock (read-only bind roots)'
    if (params.fastvep_container != resources.containers.fastvep.path ||
        params.loftee_container != resources.containers.loftee.path)
        error 'Configured SIF paths must match resource lock'
    def entries = Channel.fromList(rows).map { row ->
        def chrom = AnnotationResources.chromosome(row.chromosome)
        def meta = [unit_id: row.unit_id, chromosome: chrom, source_vcf: row.vcf,
                    annotation_root: resources.annotation_root,
                    fastvep: resources.chromosomes[chrom], loftee: resources.loftee,
                    containers: resources.containers]
        tuple(meta, file(row.vcf, checkIfExists: true))
    }
    SITES_ANNOTATION(entries,
        Channel.value(file("${projectDir}/scripts/run_annotation_stage.py")),
        Channel.value(file("${projectDir}/docs/operations/abcd-fastvep-smoke/benchmark.py")))
}

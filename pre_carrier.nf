nextflow.enable.dsl=2
include { FILTER_PRE_CARRIER } from './modules/filter_pre_carrier'
workflow {
    if(!params.filter_manifest || !params.filter_resource_lock)
        error 'Required: --filter_manifest --filter_resource_lock --select_units'
    def rows=AnnotationResources.select(PreCarrierManifest.load(file(params.filter_manifest,checkIfExists:true)),params.select_units)
    def lock=PreCarrierManifest.resources(file(params.filter_resource_lock,checkIfExists:true),rows)
    if(params.filter_resource_root!=lock.resource_root) error 'Resource root must match lock and read-only bind'
    if(params.filter_container && params.filter_container!=lock.container.path) error 'Container path must match lock'
    def entries=Channel.fromList(rows).map { row ->
        def meta=row+[popmax:lock.popmax[row.chromosome],regions:lock.regions,container:lock.container]
        tuple(meta,file(row.missense),file(row.lof_hc),file(row.sites))
    }
    FILTER_PRE_CARRIER(entries,Channel.value(file("${projectDir}/scripts/filter_pre_carrier.py")))
}

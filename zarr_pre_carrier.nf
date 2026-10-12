nextflow.enable.dsl=2

process ZARR_PRE_CARRIER {
 tag meta.unit_id
 cpus 2
 memory params.preliminary_memory
 time '4h'
 container params.carrier_container
 publishDir { "${params.outdir}/pre-carrier/${meta.unit_id}" }, mode:'copy'
 input:
 tuple val(meta),path(missense,stageAs:'input/missense.parquet'),path(lof_hc,stageAs:'input/lof_hc.parquet'),path(sites,stageAs:'input/sites.vcf.gz'),path(store,stageAs:'input/source.zarr')
 path scripts,stageAs:'scripts/*'
 path psam,stageAs:'input/samples.psam'
 output:
 path '*.parquet'
 path '*.tsv'
 path '*.json'
 script:
 def encoded=groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
 def sex=params.sex_chromosome_policy ? "--sex-chromosome-policy '${params.sex_chromosome_policy}'" : ''
 """
 printf '%s' '${encoded}' | base64 -d > unit.json
 python scripts/zarr_preliminary_frequencies.py --zarr '${store}' --chromosome '${meta.chromosome}' --missense '${missense}' --lof-hc '${lof_hc}' --psam '${psam}' ${sex} --output preliminary-frequencies.tsv --receipt preliminary-frequencies.json
 python scripts/filter_pre_carrier.py --missense '${missense}' --lof-hc '${lof_hc}' --sites '${sites}' --metadata unit.json --threads 2 --memory 3GB --cohort-frequencies preliminary-frequencies.tsv --frequency-receipt preliminary-frequencies.json --outdir .
 """
}


workflow {
    if(!params.filter_manifest || !params.filter_resource_lock)
        error 'Required: --filter_manifest --filter_resource_lock --select_units'
    def rows=AnnotationResources.select(PreCarrierManifest.load(file(params.filter_manifest,checkIfExists:true)),params.select_units)
    def lock=PreCarrierManifest.resources(file(params.filter_resource_lock,checkIfExists:true),rows)
    if(params.filter_resource_root!=lock.resource_root) error 'Resource root must match lock and read-only bind'
    if(params.filter_container && params.filter_container!=lock.container.path) error 'Container path must match lock'
    def entries=Channel.fromList(rows).map { row ->
        def meta=row+[popmax:lock.popmax[row.chromosome],regions:lock.regions,container:lock.container]
        if(!row.zarr) error 'Zarr path required'
        tuple(meta,file(row.missense),file(row.lof_hc),file(row.sites),file(row.zarr,checkIfExists:true))
    }
    ZARR_PRE_CARRIER(entries,Channel.value(['filter_pre_carrier.py','zarr_preliminary_frequencies.py','zarr_carrier_source.py','zarr_frequencies.py','carrier_frequencies.py','carrier_psam.py','extract_exact_carriers.py'].collect {file("${projectDir}/scripts/${it}")}),Channel.value(file(params.psam,checkIfExists:true)))
}

nextflow.enable.dsl=2
process ZARR_SITES {
    tag unit
    cpus 1
    memory '2 GB'
    container params.carrier_container
    publishDir { "${params.outdir}/sites/${unit}" }, mode:'copy'
    input:
    tuple val(unit),val(chromosome),path(store,stageAs:'source.zarr')
    path scripts,stageAs:'scripts/*'
    output:
    tuple val(unit),val(chromosome),path('sites.vcf.gz'),emit:sites
    path 'receipt.json',emit:receipt
    script:
    """
    python scripts/zarr_sites.py --zarr '${store}' --chromosome '${chromosome}' --output sites.vcf.gz --receipt receipt.json
    """
}
workflow {
    if(!params.zarr_manifest) error 'Require --zarr_manifest: unit_id, chromosome, zarr'
    def manifest=file(params.zarr_manifest,checkIfExists:true)
    def seen=[]
    rows=Channel.fromPath(manifest).splitCsv(header:true,sep:'\t').map {r ->
        if(!(r.unit_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || r.unit_id in seen) error 'Unsafe or duplicate unit ID'
        if(!(r.chromosome ==~ /chr([1-9]|1[0-9]|2[0-2]|X|Y)/) || !r.zarr) error 'Missing/invalid chromosome or Zarr path'
        seen.add(r.unit_id)
        tuple(r.unit_id,r.chromosome,file(manifest.parent.resolve(r.zarr).normalize(),checkIfExists:true))
    }
    ZARR_SITES(rows,Channel.value(['zarr_sites.py','zarr_carrier_source.py','extract_exact_carriers.py'].collect {file("${projectDir}/scripts/${it}",checkIfExists:true)}))
}

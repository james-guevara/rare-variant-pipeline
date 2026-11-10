nextflow.enable.dsl=2

process ZARR_FILTERED_CARRIERS {
    tag meta.unit_id
    cpus params.carrier_cpus
    memory params.carrier_memory
    time params.carrier_time
    container params.carrier_container
    publishDir { "${params.outdir}/filtered-carriers/${meta.unit_id}" }, mode:'copy'
    input:
    tuple val(meta), path(missense,stageAs:'input/missense.parquet'), path(lof_hc,stageAs:'input/lof_hc.parquet'), path(store,stageAs:'input/source.zarr')
    path scripts, stageAs:'scripts/*'
    path psam, stageAs:'input/samples.psam'
    output:
    tuple val(meta), path('carriers.tsv.gz'), emit:carriers
    path '*.tsv', emit:summaries
    path 'receipt.json', emit:receipt
    script:
    def encoded=groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    def sex=params.sex_chromosome_policy ? "--sex-chromosome-policy '${params.sex_chromosome_policy}'" : ''
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python scripts/extract_filtered_carriers.py --metadata unit.json --missense '${missense}' --lof-hc '${lof_hc}' --zarr '${store}' --psam '${psam}' --compute-frequencies ${sex} --outdir . 2> private-input.log
    """
}

workflow {
    if (!params.carrier_manifest || !params.psam) error 'Require carrier_manifest and psam; corrected frequencies are mandatory'
    def manifest=file(params.carrier_manifest,checkIfExists:true)
    def rows=manifest.readLines().findAll {it.trim()}
    def header=rows.remove(0).split('\t',-1).toList()
    if (!header.containsAll(['unit_id','chromosome','zarr','missense','lof_hc'])) error 'Require unit_id, chromosome, zarr, missense, lof_hc columns'
    def ids=[]
    def entries=rows.collect {line ->
        def values=line.split('\t',-1).toList()
        if(values.size()!=header.size()) error 'Malformed manifest row'
        def row=[header,values].transpose().collectEntries()
        if (!(row.unit_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || row.unit_id in ids) error 'Unsafe or duplicate unit_id'
        ids.add(row.unit_id)
        if (!(row.chromosome ==~ /chr([1-9]|1[0-9]|2[0-2]|X|Y)/)) error 'Invalid chromosome'
        if(row.chromosome in ['chrX','chrY'] && params.sex_chromosome_policy!='grch38_x_only_par') error 'X/Y require explicit ploidy policy'
        ['zarr','missense','lof_hc'].each {k ->
            if(!row[k]) error "Missing ${k}"
            row[k]=manifest.parent.resolve(row[k]).normalize().toString()
        }
        tuple(row,file(row.missense,checkIfExists:true),file(row.lof_hc,checkIfExists:true),file(row.zarr,checkIfExists:true))
    }
    if(!entries) error 'No units selected'
    ZARR_FILTERED_CARRIERS(Channel.fromList(entries),Channel.value([
        'extract_filtered_carriers.py','extract_exact_carriers.py','carrier_psam.py','carrier_frequencies.py',
        'zarr_carrier_source.py','zarr_frequencies.py'].collect {file("${projectDir}/scripts/${it}",checkIfExists:true)}),
        Channel.value(file(params.psam,checkIfExists:true)))
}

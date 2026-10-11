process FILTERED_CARRIERS {
    tag "${meta.unit_id}"
    cpus { params.carrier_cpus }
    memory { params.carrier_memory }
    time { params.carrier_time }
    container params.carrier_container
    publishDir { "${params.outdir}/filtered-carriers/${meta.unit_id}" }, mode:'copy',overwrite:true
    input:
    tuple val(meta),path(missense,stageAs:'input/missense.filtered.parquet'),path(lof_hc,stageAs:'input/lof_hc.filtered.parquet'),path(vcf,stageAs:'input/source.vcf.gz'),path(index,stageAs:'input/source.index')
    path runner
    path core,stageAs:'extract_exact_carriers.py'
    path psamHelper,stageAs:'carrier_psam.py'
    path psam,stageAs:'input/samples.psam'
    output:
    tuple val(meta),path('carriers.tsv.gz'),emit:carriers
    path '*.tsv',emit:summaries
    path 'receipt.json',emit:receipt
    script:
    def encoded=groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    def psamArg=psam ? "--psam '${psam}'" : ''
    def expectedHc=params.expected_hc==null ? '' : "--expected-hc ${params.expected_hc as Integer}"
    def expectedMiss=params.expected_missense==null ? '' : "--expected-missense ${params.expected_missense as Integer}"
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python '${runner}' --metadata unit.json --missense '${missense}' --lof-hc '${lof_hc}' --vcf '${vcf}' --index '${index}' --outdir . --batch-bp ${params.carrier_batch_bp as Integer} ${expectedHc} ${expectedMiss} ${psamArg} 2> private-input.log
    """
}

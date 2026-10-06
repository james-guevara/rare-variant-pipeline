process POST_CARRIER_QC {
    tag "${meta.unit_id}"
    cpus { params.qc_cpus }
    memory { params.qc_memory }
    time { params.qc_time }
    container params.qc_container
    publishDir { "${params.outdir}/post-carrier-qc/${meta.unit_id}" }, mode:'copy',overwrite:true
    input:
    tuple val(meta),path(carriers,stageAs:'input/carriers.tsv.gz'),path(samples,stageAs:'input/samples.tsv'),path(source_receipt,stageAs:'input/extraction-receipt.json')
    path runner
    output:
    tuple val(meta),path('carriers.qc.tsv.gz'),emit:carriers
    path 'qc_audit.tsv.gz',emit:audit
    path '*.tsv',emit:summaries
    path 'receipt.json',emit:receipt
    script:
    def encoded=groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python '${runner}' --metadata unit.json --carriers '${carriers}' --samples '${samples}' --source-receipt '${source_receipt}' --outdir . 2> private-input.log
    """
}

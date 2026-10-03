process FILTER_PRE_CARRIER {
    tag "${meta.unit_id}"
    cpus { params.filter_cpus }
    memory { params.filter_memory }
    time { params.filter_time }
    container params.filter_container
    publishDir { "${params.outdir}/pre-carrier/${meta.unit_id}" }, mode:'copy',overwrite:true
    input:
    tuple val(meta),path(missense,stageAs:'input/missense.parquet'),path(lof_hc,stageAs:'input/lof_hc.parquet'),path(sites,stageAs:'input/sites.vcf.gz')
    path runner
    output:
    tuple val(meta),path('missense.filtered.parquet'),path('lof_hc.filtered.parquet'),emit:filtered
    path 'filter_audit.parquet',emit:audit
    path 'receipt.json',emit:receipt
    script:
    def encoded=groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    def mem=(task.memory.toBytes()*0.70).toLong()
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python '${runner}' --missense '${missense}' --lof-hc '${lof_hc}' --sites '${sites}' --metadata unit.json --threads ${task.cpus} --memory '${mem}B' --outdir . 2> private-input.log
    """
}

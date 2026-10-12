process POST_RARITY_QC {
    tag "${meta.unit_id}"
    cpus 1
    memory { params.post_rarity_qc_memory }
    time { params.post_rarity_qc_time }
    container params.post_rarity_qc_container
    publishDir { "${params.outdir}/post-rarity-qc/${meta.unit_id}" }, mode:'copy',overwrite:true
    input:
    tuple val(meta),path(carriers,stageAs:'input/carriers.rare.tsv.gz'),path(samples,stageAs:'input/samples.tsv'),path(source_receipt,stageAs:'input/rarity-receipt.json')
    path runner
    path helpers
    path psam,stageAs:'input/samples.psam'
    output:
    tuple val(meta),path('carriers.qc.tsv.gz'),emit:carriers
    path 'qc_audit.tsv.gz',emit:audit
    path '*.tsv',emit:summaries
    path 'receipt.json',emit:receipt
    script:
    def encoded=groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    def policy=meta.sex_chromosome_policy ? "--sex-chromosome-policy ${meta.sex_chromosome_policy}" : ''
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python '${runner}' --metadata unit.json --carriers '${carriers}' --samples '${samples}' --source-receipt '${source_receipt}' --psam '${psam}' --outdir . --site-filter-policy ${meta.site_filter_policy} ${policy} 2> private-input.log
    """
}

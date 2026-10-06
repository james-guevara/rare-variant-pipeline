process GATHER_POST_QC {
    tag "${meta.chromosome}:${meta.gather_id}"
    cpus { params.gather_cpus }
    memory { params.gather_memory }
    time { params.gather_time }
    container params.gather_container
    publishDir { "${params.outdir}/post-qc-gather/${meta.chromosome}/${meta.gather_id}" }, mode:'copy',overwrite:true
    input:
    tuple val(meta),path(carriers,stageAs:'units??/carriers.qc.tsv.gz'),path(samples,stageAs:'units??/samples.tsv'),path(receipts,stageAs:'units??/receipt.json')
    path runner
    path helper
    output:
    path 'carriers.qc.tsv.gz',emit:carriers
    path '*.tsv',emit:summaries
    path 'receipt.json',emit:receipt
    script:
    def encoded=groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    def carrierArgs=(carriers instanceof List ? carriers : [carriers]).collect{"'${it}'"}.join(' ')
    def sampleArgs=(samples instanceof List ? samples : [samples]).collect{"'${it}'"}.join(' ')
    def receiptArgs=(receipts instanceof List ? receipts : [receipts]).collect{"'${it}'"}.join(' ')
    """
    printf '%s' '${encoded}' | base64 -d > selection.json
    python '${runner}' --metadata selection.json --carriers ${carrierArgs} --samples ${sampleArgs} --receipts ${receiptArgs} --outdir . 2> private-input.log
    """
}

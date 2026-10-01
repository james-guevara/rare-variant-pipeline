process EXACT_HC_CARRIERS {
    tag "${meta.unit_id}"
    cpus { params.carrier_cpus }
    memory { params.carrier_memory }
    time { params.carrier_time }
    container params.carrier_container
    publishDir { "${params.outdir}/carriers/${meta.unit_id}" }, mode: 'copy', overwrite: true
    input:
    tuple val(meta), path(loftee, stageAs: 'input/loftee.tsv'), path(vcf, stageAs: 'input/source.vcf.gz'), path(index, stageAs: 'input/source.index')
    path runner
    output:
    tuple val(meta), path('carriers.tsv.gz'), emit: carriers
    path '*.tsv', emit: summaries
    path 'receipt.json', emit: receipt
    script:
    def encoded = groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    def expected = params.expected_hc == null ? '' : "--expected-hc ${params.expected_hc as Integer}"
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python '${runner}' --metadata unit.json --loftee '${loftee}' --vcf '${vcf}' --index '${index}' --outdir . ${expected} 2> private-input.log
    """
}

process FASTVEP_PICK {
    tag "${meta.unit_id}"
    label 'sites_annotation'
    container params.fastvep_container
    publishDir { "${params.outdir}/fastvep-picker/${meta.unit_id}" }, mode: 'copy', overwrite: true
    input:
    tuple val(meta), path(vcf, stageAs: 'input/sites.vcf.gz')
    path runner
    path benchmark
    output:
    tuple val(meta), path('picked.tsv'), path('receipt.json'), emit: picked
    script:
    def encoded = groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python '${runner}' --stage fastvep --metadata unit.json --input '${vcf}' --benchmark '${benchmark}'
    """
}

process STANDALONE_LOFTEE {
    tag "${meta.unit_id}"
    label 'sites_annotation'
    container params.loftee_container
    publishDir { "${params.outdir}/loftee/${meta.unit_id}" }, mode: 'copy', overwrite: true
    input:
    tuple val(meta), path(picked), path(pick_receipt, stageAs: 'upstream-receipt.json')
    path runner
    output:
    tuple val(meta), path('loftee.tsv'), path('receipt.json'), emit: annotated
    script:
    def encoded = groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python '${runner}' --stage loftee --metadata unit.json --input '${picked}' --upstream '${pick_receipt}'
    """
}

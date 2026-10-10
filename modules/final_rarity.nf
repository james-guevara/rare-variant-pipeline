process FINAL_RARITY {
    tag "${meta.unit_id}"
    cpus 1
    memory { params.rarity_memory }
    time { params.rarity_time }
    container params.rarity_container
    publishDir { "${params.outdir}/final-rarity/${meta.unit_id}" }, mode:'copy',overwrite:true
    input:
    tuple val(meta),path(carriers,stageAs:'input/carriers.tsv.gz'),path(samples,stageAs:'input/samples.tsv'),path(source_receipt,stageAs:'input/extraction-receipt.json'),path(candidates,stageAs:'input/candidates.tsv'),path(frequencies,stageAs:'input/variant_frequencies.tsv')
    path runner
    path summaryHelper,stageAs:'qc_filtered_carriers.py'
    output:
    tuple val(meta),path('carriers.rare.tsv.gz'),emit:carriers
    path '*.tsv',emit:tables
    path 'receipt.json',emit:receipt
    script:
    def encoded=groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python '${runner}' --metadata unit.json --carriers '${carriers}' --samples '${samples}' --source-receipt '${source_receipt}' --candidates '${candidates}' --frequencies '${frequencies}' --outdir . 2> private-input.log
    """
}

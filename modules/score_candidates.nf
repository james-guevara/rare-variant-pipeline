process SCORE_CANDIDATES {
    tag "${meta.unit_id}"
    cpus { params.candidate_cpus }
    memory { params.candidate_memory }
    time { params.candidate_time }
    container params.candidate_container
    publishDir { "${params.outdir}/candidates/${meta.unit_id}" }, mode:'copy', overwrite:true
    input:
    tuple val(meta), path(picked,stageAs:'input/picked.tsv'), path(loftee,stageAs:'input/loftee.tsv')
    path runner
    path tiers, stageAs:'postprocess/tier_variants.py'
    path scores, stageAs:'postprocess/join_scores.py'
    output:
    tuple val(meta), path('missense.parquet'), path('lof_hc.parquet'), emit:candidates
    path 'receipt.json', emit:receipt
    script:
    def encoded=groovy.json.JsonOutput.toJson(meta).bytes.encodeBase64().toString()
    def mem=(task.memory.toBytes()*0.70).toLong()
    """
    printf '%s' '${encoded}' | base64 -d > unit.json
    python '${runner}' --picked '${picked}' --loftee '${loftee}' --metadata unit.json --postprocess-dir postprocess --threads ${task.cpus} --memory '${mem}B' --outdir . 2> private-input.log
    """
}

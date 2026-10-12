nextflow.enable.dsl=2
params.chromosome_manifest = null
params.expected_chromosomes = (1..22).join(',')
params.outdir = 'genome_burdens_output'
params.gather_publish_mode = 'copy'
params.gather_container = null

process GATHER_GENOME_BURDENS {
    cpus 1
    memory '4 GB'
    time '2h'
    container params.gather_container
    publishDir params.outdir, mode:params.gather_publish_mode
    input:
    val chromosomes
    path directories, stageAs:'chrom??'
    path runner
    output:
    path 'genome', emit: summaries
    script:
    def dirs = directories instanceof List ? directories : [directories]
    def manifest = 'chromosome\tdirectory\n' + chromosomes.withIndex().collect { c,i -> "${c}\t${dirs[i]}" }.join('\n') + '\n'
    def encoded = manifest.bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded}' | base64 -d > chromosomes.tsv
    python '${runner}' --manifest chromosomes.tsv --chromosomes '${params.expected_chromosomes}' --outdir genome
    """
}
workflow {
    if(!params.chromosome_manifest) error '--chromosome_manifest is required'
    def rows = file(params.chromosome_manifest,checkIfExists:true).readLines().findAll { it.trim() }
    if(rows[0] != 'chromosome\tdirectory') error 'Manifest header must be chromosome<TAB>directory'
    def entries = rows.drop(1).collect { it.split('\t',-1) }
    if(!entries || entries.any { it.size()!=2 }) error 'Invalid chromosome manifest'
    GATHER_GENOME_BURDENS(Channel.value(entries.collect { it[0] }),
        Channel.value(entries.collect { file(it[1],checkIfExists:true) }),
        Channel.value(file("${projectDir}/scripts/gather_genome_burdens.py")))
}

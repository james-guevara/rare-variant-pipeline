process BLOCK_SITES_CATALOG {
    tag unit_id
    cpus 2
    memory '2 GB'
    time '8h'
    container params.bcftools_container
    publishDir "${params.outdir}/sites", mode: params.sites_publish_mode

    input:
    tuple val(unit_id), val(chromosome), val(source_vcf), path(vcf, stageAs: 'source/input.vcf.gz')
    path runner
    path checker

    output:
    tuple val(unit_id), val(chromosome), path("${unit_id}.sites.vcf.gz"),
        path("${unit_id}.sites.vcf.gz.csi"), emit: sites
    path "${unit_id}.receipt.json", emit: receipt

    script:
    def metadata = groovy.json.JsonOutput.toJson([
        unit_id:unit_id, chromosome:chromosome, source_vcf:source_vcf
    ]).replaceFirst(/}$/, ',').bytes.encodeBase64().toString()
    """
    printf '%s' '${metadata}' | base64 -d > receipt-prefix.json
    bash '${runner}' '${unit_id}' '${chromosome}' '${vcf}' '${checker}' ${task.cpus}
    """
}

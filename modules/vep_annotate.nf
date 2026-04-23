// Module: VEP annotation (core + LOFTEE only)

process VEP_ANNOTATE {
    tag "${chrom}"
    cpus 8
    memory '32 GB'
    time '8h'
    container "${params.vep_container}"

    input:
    tuple val(chrom), path(sites_vcf), path(sites_idx)

    output:
    tuple val(chrom), path("${chrom}.vep.vcf.gz"), path("${chrom}.vep.vcf.gz.tbi"), emit: vep

    script:
    """
    vep \
        --input_file ${sites_vcf} \
        --format vcf \
        --output_file ${chrom}.vep.vcf.gz \
        --vcf \
        --compress_output bgzip \
        --canonical --mane --symbol --protein --hgvs --hgvsg --domains --biotype \
        --numbers --uploaded_allele --pick \
        --assembly GRCh38 --cache --dir_cache ${params.vep_cache} --offline \
        --fasta ${params.reference} \
        --force_overwrite --fork ${task.cpus} --stats_text \
        --dir_plugins ${params.vep_plugins} \
        --plugin LoF,loftee_path:${params.loftee_path},human_ancestor_fa:${params.resources_base}/LOFTEE/human_ancestor.fa.gz,gerp_bigwig:${params.resources_base}/LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:${params.resources_base}/LOFTEE/loftee.sql

    tabix -@ ${task.cpus} -p vcf ${chrom}.vep.vcf.gz
    """
}

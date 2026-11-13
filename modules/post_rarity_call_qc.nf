process POST_RARITY_CALL_QC {
    tag "${meta.chromosome}"
    cpus 1
    memory { params.call_qc_memory }
    time { params.call_qc_time }
    container params.call_qc_container
    publishDir { "${params.outdir}/post-rarity-call-qc/${meta.chromosome}" }, mode:params.call_qc_publish_mode,overwrite:true
    input:
    tuple val(meta),path(gather,stageAs:'gather')
    path runner
    path helpers
    path rescue,stageAs:'input/dnm_rescue.tsv.gz'
    output:
    path 'out/carriers.qc.tsv.gz',emit:carriers
    path 'out/*',emit:products
    script:
    def zarr=meta.zarr_dir ? "--zarr-dir '${meta.zarr_dir}'" : ''
    def rescueArg=meta.dnm_rescue ? "--dnm-rescue '${rescue}'" : ''
    """
    python '${runner}' --gather-dir gather --outdir out ${zarr} ${rescueArg} --snv-site-ab-min ${meta.snv_site_ab_min} --indel-site-ab-min ${meta.indel_site_ab_min} 2> private-input.log
    """
}

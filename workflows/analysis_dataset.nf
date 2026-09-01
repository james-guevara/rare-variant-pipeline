process BUILD_ANALYSIS_DATASET {
    tag params.cohort_id
    container params.targeted_container
    publishDir "${params.outdir}/analysis", mode: 'copy'
    cpus 1
    memory '2 GB'
    time '30m'

    input:
    path participant_manifest
    tuple val(pgs_enabled), path(pgs_dataset), path(pgs_dictionary)
    tuple val(rare_enabled), path(rare_burdens)
    tuple val(cnv_enabled), path(cnv_dataset), path(cnv_dictionary)
    path variable_template

    output:
    path 'analysis_dataset.tsv', emit: dataset
    path 'analysis_dictionary.tsv', emit: dictionary
    path 'analysis_qc.json', emit: qc
    path 'analysis_exclusions.tsv', emit: exclusions

    script:
    def pgsArgs = pgs_enabled ? "--pgs-dataset ${pgs_dataset} --pgs-dictionary ${pgs_dictionary}" : ''
    def rareArgs = rare_enabled ? "--rare-burdens ${rare_burdens}" : ''
    def cnvArgs = cnv_enabled ? "--cnv-dataset ${cnv_dataset} --cnv-dictionary ${cnv_dictionary}" : ''
    """
    python /opt/rvp/scripts/build_analysis_dataset.py \
      --participant-manifest ${participant_manifest} \
      ${pgsArgs} \
      ${rareArgs} \
      ${cnvArgs} \
      --variable-template ${variable_template} \
      --cohort-id '${params.cohort_id}' \
      --missing-pgs-policy '${params.missing_pgs_policy}' \
      --missing-rare-policy '${params.missing_rare_policy}' \
      --missing-cnv-policy '${params.missing_cnv_policy}' \
      --output analysis_dataset.tsv \
      --dictionary analysis_dictionary.tsv \
      --qc analysis_qc.json \
      --exclusions analysis_exclusions.tsv
    """
}

workflow ANALYSIS_DATASET_WORKFLOW {
    take:
    participant_manifest
    pgs_inputs
    rare_input
    cnv_inputs
    variable_template

    main:
    BUILD_ANALYSIS_DATASET(
        participant_manifest, pgs_inputs, rare_input, cnv_inputs, variable_template
    )

    emit:
    dataset = BUILD_ANALYSIS_DATASET.out.dataset
    dictionary = BUILD_ANALYSIS_DATASET.out.dictionary
    qc = BUILD_ANALYSIS_DATASET.out.qc
    exclusions = BUILD_ANALYSIS_DATASET.out.exclusions
}

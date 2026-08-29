process BUILD_INTEGRATED_ANALYSIS_DATASET {
    tag params.cohort_id
    container params.targeted_container
    publishDir "${params.outdir}/integrated_analysis", mode: 'copy'
    cpus 1
    memory '2 GB'
    time '30m'

    input:
    path pgs_dataset
    path pgs_dictionary
    path rare_burdens
    path variable_template

    output:
    path 'integrated_analysis_dataset.tsv', emit: dataset
    path 'integrated_analysis_dictionary.tsv', emit: dictionary
    path 'integrated_analysis_qc.json', emit: qc

    script:
    """
    python /opt/rvp/scripts/build_integrated_analysis_dataset.py \
      --pgs-dataset ${pgs_dataset} \
      --pgs-dictionary ${pgs_dictionary} \
      --rare-burdens ${rare_burdens} \
      --variable-template ${variable_template} \
      --cohort-id '${params.cohort_id}' \
      --output integrated_analysis_dataset.tsv \
      --dictionary integrated_analysis_dictionary.tsv \
      --qc integrated_analysis_qc.json
    """
}

workflow INTEGRATED_ANALYSIS_WORKFLOW {
    take:
    pgs_dataset
    pgs_dictionary
    rare_burdens
    variable_template

    main:
    BUILD_INTEGRATED_ANALYSIS_DATASET(
        pgs_dataset, pgs_dictionary, rare_burdens, variable_template
    )

    emit:
    dataset = BUILD_INTEGRATED_ANALYSIS_DATASET.out.dataset
    dictionary = BUILD_INTEGRATED_ANALYSIS_DATASET.out.dictionary
    qc = BUILD_INTEGRATED_ANALYSIS_DATASET.out.qc
}

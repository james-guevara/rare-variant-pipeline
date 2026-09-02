process BUILD_ANALYSIS_DATASET {
    tag params.cohort_id
    container params.targeted_container ?: params.python_container
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
    path analysis_builder

    output:
    path 'analysis_dataset.tsv', emit: dataset
    path 'analysis_dataset_dictionary.tsv', emit: dictionary
    path 'analysis_qc.json', emit: qc
    path 'analysis_exclusions.tsv', emit: exclusions

    script:
    def pgsArgs = pgs_enabled ? "--pgs-dataset ${pgs_dataset} --pgs-dictionary ${pgs_dictionary}" : ''
    def rareArgs = rare_enabled ? "--rare-burdens ${rare_burdens}" : ''
    def cnvArgs = cnv_enabled ? "--cnv-dataset ${cnv_dataset} --cnv-dictionary ${cnv_dictionary}" : ''
    """
    python3 ${analysis_builder} \
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
      --dictionary analysis_dataset_dictionary.tsv \
      --qc analysis_qc.json \
      --exclusions analysis_exclusions.tsv
    """
}

process BUILD_ANALYSIS_REPORTS {
    tag params.cohort_id
    container params.targeted_container ?: params.python_container
    publishDir "${params.outdir}/analysis", mode: 'copy'
    cpus 1
    memory '1 GB'
    time '15m'

    input:
    path dataset
    path dictionary
    path qc
    path report_builder
    val pgs_mode
    val rare_mode
    val cnv_mode

    output:
    path 'missingness_report.tsv', emit: missingness
    path 'run_manifest.json', emit: run_manifest

    script:
    """
    python3 ${report_builder} \
      --dataset ${dataset} \
      --dictionary ${dictionary} \
      --qc ${qc} \
      --cohort-id '${params.cohort_id}' \
      --pgs-mode '${pgs_mode}' \
      --rare-mode '${rare_mode}' \
      --cnv-mode '${cnv_mode}' \
      --missingness missingness_report.tsv \
      --run-manifest run_manifest.json
    """
}

workflow ANALYSIS_DATASET_WORKFLOW {
    take:
    participant_manifest
    pgs_inputs
    rare_input
    cnv_inputs
    variable_template
    analysis_builder
    report_builder
    pgs_mode
    rare_mode
    cnv_mode

    main:
    BUILD_ANALYSIS_DATASET(
        participant_manifest, pgs_inputs, rare_input, cnv_inputs, variable_template,
        analysis_builder
    )
    BUILD_ANALYSIS_REPORTS(
        BUILD_ANALYSIS_DATASET.out.dataset,
        BUILD_ANALYSIS_DATASET.out.dictionary,
        BUILD_ANALYSIS_DATASET.out.qc,
        report_builder, pgs_mode, rare_mode, cnv_mode
    )

    emit:
    dataset = BUILD_ANALYSIS_DATASET.out.dataset
    dictionary = BUILD_ANALYSIS_DATASET.out.dictionary
    qc = BUILD_ANALYSIS_DATASET.out.qc
    exclusions = BUILD_ANALYSIS_DATASET.out.exclusions
    missingness = BUILD_ANALYSIS_REPORTS.out.missingness
    run_manifest = BUILD_ANALYSIS_REPORTS.out.run_manifest
}

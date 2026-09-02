include { COHORT_RARE_BURDEN_WORKFLOW } from './cohort_rare_burden'
include { ANALYSIS_DATASET_WORKFLOW } from './analysis_dataset'

workflow COHORT_ANALYSIS_WORKFLOW {
    take:
    manifest_bindings
    participant_manifest
    expected_chromosomes
    pgs_inputs
    cnv_inputs
    variable_template

    main:
    COHORT_RARE_BURDEN_WORKFLOW(
        manifest_bindings,
        participant_manifest,
        expected_chromosomes
    )

    rare_input = COHORT_RARE_BURDEN_WORKFLOW.out.rare_burdens.map { burdens ->
        tuple(true, burdens)
    }

    ANALYSIS_DATASET_WORKFLOW(
        participant_manifest,
        pgs_inputs,
        rare_input,
        cnv_inputs,
        variable_template,
        Channel.value(file("${projectDir}/scripts/build_analysis_dataset.py", checkIfExists: true)),
        Channel.value(file("${projectDir}/scripts/build_analysis_reports.py", checkIfExists: true)),
        params.pgs_dataset ? 'precomputed' : 'disabled',
        'computed',
        params.cnv_dataset ? 'precomputed' : 'disabled'
    )

    emit:
    analysis_dataset = ANALYSIS_DATASET_WORKFLOW.out.dataset
    analysis_dictionary = ANALYSIS_DATASET_WORKFLOW.out.dictionary
    analysis_qc = ANALYSIS_DATASET_WORKFLOW.out.qc
    analysis_exclusions = ANALYSIS_DATASET_WORKFLOW.out.exclusions
    rare_burdens = COHORT_RARE_BURDEN_WORKFLOW.out.rare_burdens
    chromosome_strata = COHORT_RARE_BURDEN_WORKFLOW.out.chromosome_strata
    chromosome_outputs = COHORT_RARE_BURDEN_WORKFLOW.out.chromosome_outputs
}

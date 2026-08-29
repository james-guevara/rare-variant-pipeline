include { TARGETED_MANIFEST_WORKFLOW } from './targeted_manifest'
include { RARE_BURDEN_GATHER_WORKFLOW } from './rare_burden_gather'

workflow COHORT_RARE_BURDEN_WORKFLOW {
    take:
    manifest_bindings
    sample_manifest
    expected_chromosomes

    main:
    TARGETED_MANIFEST_WORKFLOW(manifest_bindings)
    RARE_BURDEN_GATHER_WORKFLOW(
        sample_manifest,
        TARGETED_MANIFEST_WORKFLOW.out.chromosome_outputs,
        expected_chromosomes
    )

    emit:
    rare_burdens = RARE_BURDEN_GATHER_WORKFLOW.out.rare_burdens
    chromosome_strata = RARE_BURDEN_GATHER_WORKFLOW.out.chromosome_strata
    chromosome_outputs = TARGETED_MANIFEST_WORKFLOW.out.chromosome_outputs
}

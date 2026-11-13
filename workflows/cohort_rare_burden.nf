include { TARGETED_MANIFEST_WORKFLOW } from './targeted_manifest'
include { RARE_BURDEN_GATHER_WORKFLOW } from './rare_burden_gather'

workflow COHORT_RARE_BURDEN_WORKFLOW {
    take:
    manifest_bindings
    sample_manifest
    expected_chromosomes
    gene_set_membership
    consolidate_script
    gene_set_script

    main:
    TARGETED_MANIFEST_WORKFLOW(manifest_bindings)
    RARE_BURDEN_GATHER_WORKFLOW(
        sample_manifest,
        TARGETED_MANIFEST_WORKFLOW.out.chromosome_outputs,
        expected_chromosomes,
        gene_set_membership,
        consolidate_script,
        gene_set_script
    )

    emit:
    rare_burdens = RARE_BURDEN_GATHER_WORKFLOW.out.rare_burdens
    chromosome_strata = RARE_BURDEN_GATHER_WORKFLOW.out.chromosome_strata
    plof_carriers = RARE_BURDEN_GATHER_WORKFLOW.out.plof_carriers
    missense_carriers = RARE_BURDEN_GATHER_WORKFLOW.out.missense_carriers
    gene_set_burdens = RARE_BURDEN_GATHER_WORKFLOW.out.gene_set_burdens
    chromosome_outputs = TARGETED_MANIFEST_WORKFLOW.out.chromosome_outputs
}

process GATHER_RARE_BURDENS {
    tag 'cohort'
    container params.targeted_container
    cpus 1
    memory '4 GB'
    time '1h'

    input:
    path sample_manifest
    path chromosome_packages, stageAs: 'packages/package??'
    val expected_chromosomes

    output:
    path 'rare_burdens.tsv', emit: burdens
    path 'rare_burdens_by_chromosome_stratum.tsv', emit: strata

    script:
    def packageArgs = chromosome_packages.collect { "--package ${it}" }.join(' ')
    """
    python /opt/rvp/scripts/gather_rare_burdens.py \
      --sample-manifest ${sample_manifest} \
      ${packageArgs} \
      --expected-chromosomes '${expected_chromosomes}' \
      --output rare_burdens.tsv \
      --strata-output rare_burdens_by_chromosome_stratum.tsv
    """
}

workflow RARE_BURDEN_GATHER_WORKFLOW {
    take:
    sample_manifest
    chromosome_packages
    expected_chromosomes

    main:
    packages = chromosome_packages.collect()
    GATHER_RARE_BURDENS(sample_manifest, packages, expected_chromosomes)

    emit:
    rare_burdens = GATHER_RARE_BURDENS.out.burdens
    chromosome_strata = GATHER_RARE_BURDENS.out.strata
}

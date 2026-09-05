process GATHER_RARE_BURDENS {
    tag 'cohort'
    container params.rare_gather_container ?: params.targeted_container
    publishDir "${params.outdir}/rare_burdens", mode: 'copy', pattern: 'rare_burdens/*', saveAs: { it.tokenize('/').last() }
    publishDir "${params.outdir}/eligible_carriers", mode: 'copy', pattern: 'eligible_carriers/*', saveAs: { it.tokenize('/').last() }
    publishDir "${params.outdir}/gene_set_burdens", mode: 'copy', pattern: 'gene_set_burdens/*', saveAs: { it.tokenize('/').last() }
    publishDir "${params.outdir}/validation", mode: 'copy', pattern: 'validation/*', saveAs: { it.tokenize('/').last() }
    cpus 1
    memory '4 GB'
    time '1h'

    input:
    path sample_manifest
    path chromosome_packages, stageAs: 'packages/package??'
    val expected_chromosomes
    path gene_set_membership
    path consolidate_script
    path gene_set_script

    output:
    path 'rare_burdens/rare_burdens.tsv', emit: burdens
    path 'rare_burdens/rare_burdens_by_chromosome_stratum.tsv', emit: strata
    path 'eligible_carriers/plof-burden-eligible.parquet', emit: plof_carriers
    path 'eligible_carriers/missense-burden-eligible.parquet', emit: missense_carriers
    path 'gene_set_burdens/gene_set_burdens.tsv', emit: gene_set_burdens
    path 'validation/plof-carrier-consolidation-validation.json', emit: plof_validation
    path 'validation/missense-carrier-consolidation-validation.json', emit: missense_validation

    script:
    """
    mkdir -p rare_burdens eligible_carriers gene_set_burdens validation
    python /opt/rvp/scripts/gather_rare_burdens.py \
      --sample-manifest ${sample_manifest} \
      --package-root packages \
      --expected-chromosomes '${expected_chromosomes}' \
      --output rare_burdens/rare_burdens.tsv \
      --strata-output rare_burdens/rare_burdens_by_chromosome_stratum.tsv
    python ${consolidate_script} \
      --input packages/package*/11.plof-burden-eligible.parquet \
      --output eligible_carriers/plof-burden-eligible.parquet \
      --validation-output validation/plof-carrier-consolidation-validation.json
    python ${consolidate_script} \
      --input packages/package*/11.missense-burden-eligible.parquet \
      --output eligible_carriers/missense-burden-eligible.parquet \
      --validation-output validation/missense-carrier-consolidation-validation.json
    python ${gene_set_script} \
      --samples ${sample_manifest} \
      --gene-sets ${gene_set_membership} \
      --plof eligible_carriers/plof-burden-eligible.parquet \
      --missense eligible_carriers/missense-burden-eligible.parquet \
      --output gene_set_burdens/gene_set_burdens.tsv
    """
}

workflow RARE_BURDEN_GATHER_WORKFLOW {
    take:
    sample_manifest
    chromosome_packages
    expected_chromosomes
    gene_set_membership
    consolidate_script
    gene_set_script

    main:
    packages = chromosome_packages.collect()
    GATHER_RARE_BURDENS(
        sample_manifest, packages, expected_chromosomes, gene_set_membership,
        consolidate_script, gene_set_script
    )

    emit:
    rare_burdens = GATHER_RARE_BURDENS.out.burdens
    chromosome_strata = GATHER_RARE_BURDENS.out.strata
    plof_carriers = GATHER_RARE_BURDENS.out.plof_carriers
    missense_carriers = GATHER_RARE_BURDENS.out.missense_carriers
    gene_set_burdens = GATHER_RARE_BURDENS.out.gene_set_burdens
    plof_validation = GATHER_RARE_BURDENS.out.plof_validation
    missense_validation = GATHER_RARE_BURDENS.out.missense_validation
}

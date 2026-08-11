// Subworkflow: POSTPROCESS
//
// Wraps the DuckDB postprocessing stages (previously run by
// rv_postprocessing_v3/run_chrom.sh as a per-chrom SLURM job) as Nextflow
// processes, so they get the same -resume, tracing and fan-out as the rest of the
// pipeline. The scripts themselves are unchanged apart from optional
// --input/--output overrides; with those omitted they behave exactly as before,
// which is what makes this migration verifiable against the existing outputs.
//
// Deliverables this produces:
//   D2  filtered + fully annotated variant x carrier table (per chrom), the
//       terminus of the chain. Stratification-agnostic: carries every score,
//       constraint metric and a `tier` column, so tiers, gene sets, consequence
//       classes etc. can all be selected downstream.
//   D3  per-sample carrier counts, grouped by params.count_group_col.
//   (D1 is the upstream EXPORT parquet, before any of this filtering.)
//
// Per-chrom chain:
//   filter_regions -> qc_genotype -> join_pop_af -> join_scores
//     -> join_gene_constraint -> tier_variants(annotate) -> count_carriers
//
// Notes on what this chain deliberately does NOT do:
//   * tier_variants reads with_gene_constraint directly, which is the existing
//     behaviour and is unaffected by the removal above.
//   * No DNM rescue and no FILTER cut. Both were deliberately removed: a DNM
//     rescue can only fire where trios exist, so the same variant would be
//     filtered differently between cohorts with and without them, which biases
//     cross-cohort comparison. FILTER is carried as a column instead, so the cut
//     can be applied downstream where it is visible. De novo status belongs in a
//     later annotation step, not a filter here.
//
// All PP_* processes run in params.python_container (they need duckdb).

def pp (name) { "python ${projectDir}/scripts/postprocess/${name}" }

process PP_FILTER_REGIONS {
    tag "${chrom}"
    container "${params.python_container}"
    input:  tuple val(chrom), path(in_parquet)
    output: tuple val(chrom), path("${chrom}.region_filtered.parquet"), emit: out
    script:
    """
    ${pp('filter_regions.py')} \\
        --cohort ${params.pp_cohort} --chrom ${chrom} --resources ${params.pp_resources} \\
        --input ${in_parquet} --output ${chrom}.region_filtered.parquet
    """
}

process PP_QC_GENOTYPE {
    tag "${chrom}"
    container "${params.python_container}"
    input:  tuple val(chrom), path(in_parquet)
    output: tuple val(chrom), path("${chrom}.qc_filtered.parquet"), emit: out
    script:
    """
    ${pp('qc_genotype.py')} \\
        --cohort ${params.pp_cohort} --chrom ${chrom} --resources ${params.pp_resources} \\
        --input ${in_parquet} --output ${chrom}.qc_filtered.parquet
    """
}

process PP_JOIN_POP_AF {
    tag "${chrom}"
    container "${params.python_container}"
    input:  tuple val(chrom), path(in_parquet)
    output: tuple val(chrom), path("${chrom}.with_pop_af.parquet"), emit: out
    script:
    """
    ${pp('join_pop_af.py')} \\
        --cohort ${params.pp_cohort} --chrom ${chrom} --resources ${params.pp_resources} \\
        --input ${in_parquet} --output ${chrom}.with_pop_af.parquet
    """
}

process PP_JOIN_SCORES {
    tag "${chrom}"
    container "${params.python_container}"
    input:  tuple val(chrom), path(in_parquet)
    output: tuple val(chrom), path("${chrom}.with_scores.parquet"), emit: out
    script:
    """
    ${pp('join_scores.py')} \\
        --cohort ${params.pp_cohort} --chrom ${chrom} --resources ${params.pp_resources} \\
        --input ${in_parquet} --output ${chrom}.with_scores.parquet
    """
}

process PP_JOIN_GENE_CONSTRAINT {
    tag "${chrom}"
    container "${params.python_container}"
    input:  tuple val(chrom), path(in_parquet)
    output: tuple val(chrom), path("${chrom}.with_gene_constraint.parquet"), emit: out
    script:
    """
    ${pp('join_gene_constraint.py')} \\
        --cohort ${params.pp_cohort} --chrom ${chrom} --resources ${params.pp_resources} \\
        --input ${in_parquet} --output ${chrom}.with_gene_constraint.parquet
    """
}

process PP_TIER_VARIANTS {
    tag "${chrom}"
    container "${params.python_container}"
    input:  tuple val(chrom), path(in_parquet)
    output: tuple val(chrom), path("${chrom}.filtered_annotated.parquet"), emit: out
    script:
    """
    ${pp('tier_variants.py')} \\
        --cohort ${params.pp_cohort} --chrom ${chrom} --resources ${params.pp_resources} \\
        --input ${in_parquet} --output ${chrom}.filtered_annotated.parquet
    """
}

// Cohort-level: per-sample carrier counts, stratified by params.count_group_col
// (default `tier`; any annotated column works, e.g. Consequence or a gene-set
// membership flag). Deliberately carries NO covariates (no PCs / ancestry / case
// status / FID) and does no cross-cohort union — those belong to analysis.
process PP_COUNT_CARRIERS {
    tag "${params.pp_cohort}:${params.count_group_col}"
    container "${params.python_container}"
    input:  path(annotated_parquets)
    output: path("per_sample_counts.tsv"), emit: counts
            path("group_totals.tsv"),      emit: totals
    script:
    """
    ${pp('count_carriers.py')} \\
        --input ${annotated_parquets.join(' ')} \\
        --group-col ${params.count_group_col} \\
        --out-counts per_sample_counts.tsv \\
        --out-totals group_totals.tsv
    """
}

workflow POSTPROCESS {
    take:
    variant_parquets    // tuple(chrom, parquet) — the pipeline's EXPORT output

    main:
    PP_FILTER_REGIONS(variant_parquets)
    PP_QC_GENOTYPE(PP_FILTER_REGIONS.out.out)
    PP_JOIN_POP_AF(PP_QC_GENOTYPE.out.out)
    PP_JOIN_SCORES(PP_JOIN_POP_AF.out.out)
    PP_JOIN_GENE_CONSTRAINT(PP_JOIN_SCORES.out.out)

    // Tier ANNOTATION (adds a `tier` column, drops nothing) — the terminus of the
    // filter/annotate chain and the pipeline's analysis-ready deliverable.
    PP_TIER_VARIANTS(PP_JOIN_GENE_CONSTRAINT.out.out)

    PP_COUNT_CARRIERS(PP_TIER_VARIANTS.out.out.map { chrom, pq -> pq }.collect())

    emit:
    annotated = PP_TIER_VARIANTS.out.out       // D2: filtered + fully annotated
    counts    = PP_COUNT_CARRIERS.out.counts   // D3: per-sample carrier counts
    totals    = PP_COUNT_CARRIERS.out.totals
}

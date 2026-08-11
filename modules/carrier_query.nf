// Module: Query carrier genotypes from a normed cohort VCF.
//
// Per-sample carrier filter via bcftools' GT="alt" predicate, plus an
// optional INFO/AF cap to drop common ALT records that share a position
// with rare targets (bcftools -R is region-based, so co-located common
// ALTs leak in without an AF check). One row per (variant × carrier sample).

// Resource sizing: bcftools query streams the VCF and is single-threaded;
// observed RSS ~50 MB for 9k-sample SSC chunks, independent of BED size.
process CARRIER_QUERY {
    tag "${chrom}_${chunk_bed.baseName}"
    cpus 1
    memory '2 GB'
    time '4h'
    container "${params.bcftools_container}"

    input:
    tuple val(chrom), path(chunk_bed), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("${chunk_bed.baseName}.carriers.tsv"), emit: carriers

    script:
    def af_cap = (params.max_cohort_af as Double) < 1.0 ? "INFO/AF<${params.max_cohort_af} && " : ""
    def filter_expr = "${af_cap}ALT!=\"*\" && GT=\"alt\""
    // FORMAT/FT is the per-sample filter. Some callers put their genotype-level
    // filtering there and leave the site-level FILTER constant, in which case FT is
    // the only quality filter that varies and dropping it loses that information.
    // Safe where the field is absent: bcftools prints '.'.
    """
    bcftools query \\
        -HH \\
        -R ${chunk_bed} \\
        -i '${filter_expr}' \\
        -f '[%CHROM\\t%POS0\\t%END\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\t%GQ\\t%DP\\t%PL\\t%AD\\t%FT\\n]' \\
        ${vcf} \\
        > ${chunk_bed.baseName}.carriers.tsv
    """
}

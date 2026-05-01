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
    """
    bcftools query \\
        -HH \\
        -R ${chunk_bed} \\
        -i '${af_cap}GT="alt"' \\
        -f '[%CHROM\\t%POS0\\t%END\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\t%GQ\\t%DP\\t%AD{0}\\t%AD{1}\\n]' \\
        ${vcf} \\
        | sed '1s/\\tAD\\tAD\$/\\tAD_ref\\tAD_alt/' \\
        > ${chunk_bed.baseName}.carriers.tsv
    """
}

// Subworkflow: Merge and Index (Parts 9-12)
// Merge annotations → QC Filter → Sort/Bgzip/Tabix → Parquet

include { MERGE_ANNOTATIONS } from '../modules/merge_annotations'
include { APPLY_QC_FILTER } from '../modules/qc_filter'
include { SORT_INDEX } from '../modules/sort_index'
include { CONVERT_PARQUET } from '../modules/convert_parquet'

workflow MERGE_INDEX {
    take:
    resolved        // tuple(chrom, resolved.tsv.gz)
    consequential   // tuple(chrom, consequential.tsv)

    main:
    // Join resolved genotypes with consequential annotations
    merge_input = resolved.join(consequential)

    // Merge (inner join)
    MERGE_ANNOTATIONS(merge_input)

    // Apply QC filter (adds qc_filter + AB columns)
    APPLY_QC_FILTER(MERGE_ANNOTATIONS.out.merged)

    // Sort, bgzip, tabix
    SORT_INDEX(APPLY_QC_FILTER.out.filtered)

    // Convert to Parquet
    CONVERT_PARQUET(SORT_INDEX.out.indexed)

    emit:
    indexed = SORT_INDEX.out.indexed      // tuple(chrom, tsv.gz, tbi)
    parquet = CONVERT_PARQUET.out.parquet  // tuple(chrom, parquet)
}

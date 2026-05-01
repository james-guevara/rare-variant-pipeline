// Subworkflow: Merge and Index
// Merge annotations → Sort/Bgzip/Tabix → Parquet
// QC filtering is intentionally deferred to downstream post-processing.

include { MERGE_ANNOTATIONS } from '../modules/merge_annotations'
include { SORT_INDEX } from '../modules/sort_index'
include { CONVERT_PARQUET } from '../modules/convert_parquet'

workflow MERGE_INDEX {
    take:
    carriers        // tuple(chrom, carriers.tsv)
    consequential   // tuple(chrom, consequential.tsv)

    main:
    merge_input = carriers.join(consequential)

    MERGE_ANNOTATIONS(merge_input)

    SORT_INDEX(MERGE_ANNOTATIONS.out.merged)

    CONVERT_PARQUET(SORT_INDEX.out.indexed)

    emit:
    indexed = SORT_INDEX.out.indexed
    parquet = CONVERT_PARQUET.out.parquet
}

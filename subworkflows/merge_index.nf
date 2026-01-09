// Subworkflow: Merge and Index (Parts 9-11)
// Merge annotations → Sort/Bgzip/Tabix

include { MERGE_ANNOTATIONS } from '../modules/merge_annotations'
include { SORT_INDEX } from '../modules/sort_index'

workflow MERGE_INDEX {
    take:
    resolved        // tuple(chrom, resolved.tsv.gz)
    consequential   // tuple(chrom, consequential.tsv)

    main:
    // Join resolved genotypes with consequential annotations
    merge_input = resolved.join(consequential)

    // Merge (inner join)
    MERGE_ANNOTATIONS(merge_input)

    // Sort, bgzip, tabix
    SORT_INDEX(MERGE_ANNOTATIONS.out.merged)

    emit:
    indexed = SORT_INDEX.out.indexed  // tuple(chrom, tsv.gz, tbi)
}

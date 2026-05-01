// Subworkflow: Family Processing
// Scatter consequential BED into chunks, query carriers from the normed VCF,
// gather per-chrom. RESOLVE_GENOTYPES is no longer needed because CARRIER_QUERY
// now emits carriers-only rows with AD pre-split and FAMILY attached.

include { SCATTER_BED } from '../modules/scatter_bed'
include { CARRIER_QUERY } from '../modules/carrier_query'
include { GATHER_GENOTYPES } from '../modules/gather_genotypes'

workflow FAMILY_PROCESSING {
    take:
    consequential_bed   // tuple(chrom, bed) - regions to query
    input_vcfs          // tuple(chrom, vcf, tbi) - normed cohort VCF with genotypes

    main:
    SCATTER_BED(consequential_bed)

    chunk_counts = SCATTER_BED.out.chunks
        .map { chrom, chunks, n_chunks -> tuple(chrom, n_chunks.toInteger()) }

    chunks_with_vcf = SCATTER_BED.out.chunks
        .join(input_vcfs)
        .flatMap { chrom, chunks, n_chunks, vcf, tbi ->
            def chunkList = chunks instanceof List ? chunks : [chunks]
            chunkList.collect { chunk_file ->
                return tuple(chrom, chunk_file, vcf, tbi)
            }
        }

    CARRIER_QUERY(chunks_with_vcf)

    // Gather chunks back into per-chrom files
    carriers_with_count = CARRIER_QUERY.out.carriers
        .combine(chunk_counts, by: 0)

    gathered = carriers_with_count
        .map { chrom, carriers_file, n_chunks -> tuple( groupKey(chrom, n_chunks), carriers_file ) }
        .groupTuple()
        .map { chrom_key, carrier_files -> tuple(chrom_key.toString(), carrier_files) }

    GATHER_GENOTYPES(gathered)

    emit:
    resolved = GATHER_GENOTYPES.out.resolved
}

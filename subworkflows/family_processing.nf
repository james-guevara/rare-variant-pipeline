// Subworkflow: Family Processing (Parts 5-8)
// Scatter → Family Query → Gather → Resolve

include { SCATTER_BED } from '../modules/scatter_bed'
include { FAMILY_QUERY } from '../modules/family_query'
include { GATHER_GENOTYPES } from '../modules/gather_genotypes'
include { RESOLVE_GENOTYPES } from '../modules/resolve_genotypes'

workflow FAMILY_PROCESSING {
    take:
    consequential_bed   // tuple(chrom, bed) - regions to query
    input_vcfs          // tuple(chrom, vcf, tbi) - original VCFs with genotypes

    main:
    // Scatter BED into chunks - now emits (chrom, chunks, n_chunks)
    SCATTER_BED(consequential_bed)

    // Extract chunk counts per chromosome for groupTuple
    chunk_counts = SCATTER_BED.out.chunks
        .map { chrom, chunks, n_chunks -> tuple(chrom, n_chunks.toInteger()) }

    // Join scatter output with input VCFs and flatten chunks
    chunks_with_vcf = SCATTER_BED.out.chunks
        .join(input_vcfs)
        .flatMap { chrom, chunks, n_chunks, vcf, tbi ->
            def chunkList = chunks instanceof List ? chunks : [chunks]
            chunkList.collect { chunk_file ->
                return tuple(chrom, chunk_file, vcf, tbi)
            }
        }

    // Query family genotypes per chunk
    FAMILY_QUERY(chunks_with_vcf)

    // Add chunk count back to FAMILY_QUERY output for groupTuple sizing
    family_query_with_count = FAMILY_QUERY.out.genotypes
        .combine(chunk_counts, by: 0)  // (chrom, genotypes_file, n_chunks)

    // Gather chunks by chromosome with known size - enables immediate completion
    gathered = family_query_with_count
        .map { chrom, genotypes, n_chunks -> tuple( groupKey(chrom, n_chunks), genotypes ) }
        .groupTuple()
        .map { chrom_key, genotypes -> tuple(chrom_key.toString(), genotypes) }

    GATHER_GENOTYPES(gathered)

    // Resolve genotypes (one row per variant per family)
    RESOLVE_GENOTYPES(GATHER_GENOTYPES.out.genotypes)

    emit:
    resolved = RESOLVE_GENOTYPES.out.resolved
}

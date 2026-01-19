// Subworkflow: Family Processing (Parts 5-8)
// Scatter → Family Query → Resolve (per chunk) → Gather

include { SCATTER_BED } from '../modules/scatter_bed'
include { FAMILY_QUERY } from '../modules/family_query'
include { RESOLVE_GENOTYPES } from '../modules/resolve_genotypes'
include { GATHER_GENOTYPES } from '../modules/gather_genotypes'

workflow FAMILY_PROCESSING {
    take:
    consequential_bed   // tuple(chrom, bed) - regions to query
    input_vcfs          // tuple(chrom, vcf, tbi) - original VCFs with genotypes

    main:
    // Scatter BED into chunks - emits (chrom, chunks, n_chunks)
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

    // Resolve genotypes per chunk (low memory - each chunk is small)
    RESOLVE_GENOTYPES(FAMILY_QUERY.out.genotypes)

    // Add chunk count to resolved output for groupTuple sizing
    resolved_with_count = RESOLVE_GENOTYPES.out.resolved
        .combine(chunk_counts, by: 0)  // (chrom, resolved_file, n_chunks)

    // Gather resolved chunks by chromosome with known size
    gathered = resolved_with_count
        .map { chrom, resolved, n_chunks -> tuple( groupKey(chrom, n_chunks), resolved ) }
        .groupTuple()
        .map { chrom_key, resolved_files -> tuple(chrom_key.toString(), resolved_files) }

    GATHER_GENOTYPES(gathered)

    emit:
    resolved = GATHER_GENOTYPES.out.resolved
}

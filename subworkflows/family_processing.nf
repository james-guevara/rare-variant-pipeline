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
    // Scatter BED into chunks
    SCATTER_BED(consequential_bed)

    // Join scatter output with input VCFs
    chunks_with_vcf = SCATTER_BED.out.chunks
        .join(input_vcfs)
        .flatMap { chrom, chunks, vcf, tbi ->
            def chunkList = chunks instanceof List ? chunks : [chunks]; chunkList.collect { chunk_file ->
                def chunk_id = chunk_file.name.replace('chunk_', '').replace('.bed', '')
                return tuple(chrom, chunk_file, vcf, tbi)
            }
        }

    // Query family genotypes per chunk
    FAMILY_QUERY(chunks_with_vcf)

    // Gather chunks by chromosome (wait for all chunks per chrom)
    gathered = FAMILY_QUERY.out.genotypes.groupTuple()
    GATHER_GENOTYPES(gathered)

    // Resolve genotypes (one row per variant per family)
    RESOLVE_GENOTYPES(GATHER_GENOTYPES.out.genotypes)

    emit:
    resolved = RESOLVE_GENOTYPES.out.resolved  // tuple(chrom, resolved.tsv.gz)
}

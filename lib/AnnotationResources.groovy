import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Paths

/** Only metadata is checked at launch; lock creation hashes immutable resources once. */
class AnnotationResources {
    static Map load(def lock, List rows) {
        def data = new JsonSlurper().parseText(lock.text)
        if (data.schema != 1) throw new IllegalArgumentException('Unsupported annotation resource lock')
        def identities = data.containers.values() + data.loftee.values()
        rows.collect { chromosome(it.chromosome) }.unique().each { chrom ->
            if (!data.chromosomes[chrom]) throw new IllegalArgumentException('No resources locked for '+chrom)
            identities.addAll(data.chromosomes[chrom].values())
            if (data.chromosomes[chrom].cache.mtime_ns <= data.chromosomes[chrom].gff3.mtime_ns)
                throw new IllegalArgumentException('Cache must be newer than GFF3: '+chrom)
        }
        identities.unique { it.path }.each { item ->
            def path = Paths.get(item.path)
            if (!Files.isRegularFile(path) || Files.size(path) != item.bytes ||
                Files.getLastModifiedTime(path).toMillis() != item.mtime_ms)
                throw new IllegalArgumentException('Resource changed/missing; verify and regenerate lock: '+item.path)
            if (!(item.sha256 ==~ /[a-f0-9]{64}/)) throw new IllegalArgumentException('Invalid resource checksum')
        }
        data
    }
    static String chromosome(String value) { value.startsWith('chr') ? value : 'chr'+value }
    static List select(List rows, def selection) {
        if (!selection) throw new IllegalArgumentException('Required: --select_units explicit IDs or all')
        if (selection.toString() == 'all') return rows
        def ids = selection.toString().split(',', -1).toList()
        if (ids.toSet().size() != ids.size() || ids.any { !it || !(it in rows*.unit_id) })
            throw new IllegalArgumentException('Unknown, empty, or duplicate --select_units ID')
        rows.findAll { it.unit_id in ids }
    }
}

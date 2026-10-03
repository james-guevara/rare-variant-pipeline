import java.nio.file.Files
import java.nio.file.Paths

class FilteredCarrierManifest {
    static List load(def manifest) {
        def rows = SitesCatalogManifest.load(manifest)
        def lines = manifest.readLines().findAll { it.trim() }
        def header = lines[0].split('\t', -1).toList()
        if (!header.containsAll(['missense', 'lof_hc', 'index']))
            throw new IllegalArgumentException('Filtered carrier manifest requires unit_id, chromosome, missense, lof_hc, vcf, index')
        rows.eachWithIndex { row, i ->
            def values = lines[i+1].split('\t', -1).toList()
            ['missense', 'lof_hc', 'index'].each { name ->
                def value = values[header.indexOf(name)].trim()
                if (!value) throw new IllegalArgumentException('Missing '+name+' for '+row.unit_id)
                def path = Paths.get(value)
                if (!path.isAbsolute()) path = manifest.parent.resolve(path)
                path = path.toAbsolutePath().normalize()
                if (!Files.isRegularFile(path)) throw new IllegalArgumentException('Missing '+name+' for '+row.unit_id)
                row[name] = path.toString()
            }
        }
        rows
    }
}

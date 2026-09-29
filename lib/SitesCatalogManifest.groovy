import java.nio.file.Files

/** Block identities are explicit. Subset manifests are valid catalog inputs. */
class SitesCatalogManifest {
    static List load(def manifest) {
        def lines = manifest.readLines().findAll { it.trim() }
        if (lines.size() < 2) throw new IllegalArgumentException('Sites manifest has no block rows')
        def header = lines[0].split('\t', -1).toList()
        if (header.toSet().size() != header.size() || !header.containsAll(['unit_id','chromosome','vcf']))
            throw new IllegalArgumentException('Sites manifest requires unique columns: unit_id, chromosome, vcf')
        def seen = new HashSet<String>()
        lines.drop(1).collect { line ->
            def fields = line.split('\t', -1).toList()
            if (fields.size() != header.size()) throw new IllegalArgumentException('Malformed sites manifest row')
            def row = [header, fields].transpose().collectEntries { [(it[0]): it[1].trim()] }
            if (!(row.unit_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
                throw new IllegalArgumentException('Explicit safe unit_id required in every sites manifest row')
            if (!seen.add(row.unit_id)) throw new IllegalArgumentException('Duplicate unit_id: '+row.unit_id)
            if (!(row.chromosome ==~ /(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|XY|M|MT|PAR1|PAR2)/))
                throw new IllegalArgumentException('Invalid chromosome for '+row.unit_id+': '+row.chromosome)
            if (!row.vcf) throw new IllegalArgumentException('Missing VCF path for '+row.unit_id)
            def source = java.nio.file.Paths.get(row.vcf)
            if (!source.isAbsolute()) source = manifest.parent.resolve(source)
            source = source.toAbsolutePath().normalize()
            if (!Files.isRegularFile(source)) throw new IllegalArgumentException('Missing VCF for '+row.unit_id+': '+source)
            [unit_id:row.unit_id, chromosome:row.chromosome, vcf:source.toString()]
        }
    }
}

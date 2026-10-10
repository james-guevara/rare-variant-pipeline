import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Paths

class FinalRarityManifest {
    static List load(def manifest) {
        def lines=manifest.readLines().findAll{it.trim()}
        if(lines.size()<2) throw new IllegalArgumentException('Empty final rarity manifest')
        def header=lines[0].split('\t',-1).toList()
        if(header.toSet().size()!=header.size() || !header.containsAll(['unit_id','chromosome','carriers','samples','source_receipt','candidates','frequencies']))
            throw new IllegalArgumentException('Required: unit_id chromosome carriers samples source_receipt candidates frequencies')
        def seen=[] as Set
        lines.drop(1).collect { line ->
            def values=line.split('\t',-1).toList()
            if(values.size()!=header.size()) throw new IllegalArgumentException('Invalid manifest row')
            def row=[header,values].transpose().collectEntries{[(it[0]):it[1].trim()]}
            if(!(row.unit_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || !seen.add(row.unit_id))
                throw new IllegalArgumentException('Unique safe unit_id required')
            row.chromosome=AnnotationResources.chromosome(row.chromosome)
            if(!(row.chromosome ==~ /chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)/)) throw new IllegalArgumentException('Invalid chromosome')
            ['carriers','samples','source_receipt','candidates','frequencies'].each { key ->
                if(!row[key]) throw new IllegalArgumentException('Missing '+key)
                def path=Paths.get(row[key]); if(!path.isAbsolute()) path=manifest.parent.resolve(path)
                if(!Files.isRegularFile(path)) throw new IllegalArgumentException('Missing '+key+' for '+row.unit_id)
                row[key]=path.toAbsolutePath().normalize().toString()
            }
            row
        }
    }
}

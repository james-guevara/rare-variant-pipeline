import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Paths

class CandidateManifest {
    static List load(def manifest) {
        def lines = manifest.readLines().findAll { it.trim() }
        if (lines.size()<2) throw new IllegalArgumentException('Candidate manifest has no rows')
        def header = lines[0].split('\t', -1).toList()
        if (header.toSet().size()!=header.size() || !header.containsAll(['unit_id','chromosome','picked','loftee']))
            throw new IllegalArgumentException('Required manifest columns: unit_id, chromosome, picked, loftee')
        def seen = [] as Set
        lines.drop(1).collect { line ->
            def values=line.split('\t',-1).toList()
            if (values.size()!=header.size()) throw new IllegalArgumentException('Malformed candidate row')
            def row=[header,values].transpose().collectEntries { [(it[0]):it[1].trim()] }
            if (!(row.unit_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || !seen.add(row.unit_id))
                throw new IllegalArgumentException('Explicit unique safe unit_id required')
            row.chromosome=AnnotationResources.chromosome(row.chromosome)
            if (!(row.chromosome ==~ /chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)/))
                throw new IllegalArgumentException('Unsupported chromosome')
            ['picked','loftee'].each { key ->
                if (!row[key]) throw new IllegalArgumentException('Missing '+key)
                def path=Paths.get(row[key])
                if (!path.isAbsolute()) path=manifest.parent.resolve(path)
                if (!Files.isRegularFile(path)) throw new IllegalArgumentException('Missing '+key+' file for '+row.unit_id)
                row[key]=path.toAbsolutePath().normalize().toString()
            }
            row
        }
    }
    static Map resources(def path, List rows) {
        def data=new JsonSlurper().parseText(path.text)
        if (data.schema!=2 || data.dbnsfp_representation!='parquet_expanded')
            throw new IllegalArgumentException('Rebuild candidate lock with unfiltered parquet_expanded resources; old MANE-filtered locks are unsupported')
        def items=[data.genebayes,data.container]
        rows*.chromosome.unique().each { chrom ->
            if (!data.dbnsfp[chrom]) throw new IllegalArgumentException('No locked dbNSFP for '+chrom)
            items.add(data.dbnsfp[chrom])
        }
        items.each { item ->
            def file=Paths.get(item.path)
            if (!Files.isRegularFile(file) || Files.size(file)!=item.bytes || Files.getLastModifiedTime(file).toMillis()!=item.mtime_ms)
                throw new IllegalArgumentException('Resource changed since lock creation: '+item.path)
        }
        data
    }
}

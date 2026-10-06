nextflow.enable.dsl=2
include { SCORE_CANDIDATES } from './modules/score_candidates'

workflow {
    if (!params.candidate_manifest || !params.candidate_resource_lock)
        error 'Required: --candidate_manifest, --candidate_resource_lock, --select_units'
    def rows=AnnotationResources.select(CandidateManifest.load(file(params.candidate_manifest,checkIfExists:true)),params.select_units)
    def lock=CandidateManifest.resources(file(params.candidate_resource_lock,checkIfExists:true),rows)
    if (params.candidate_container && params.candidate_container != lock.container.path)
        error 'Container path must match candidate resource lock'
    if (params.candidate_resource_root != lock.resource_root)
        error 'Resource root must match lock and read-only bind'
    def entries=Channel.fromList(rows).map { row ->
        def meta=row+[dbnsfp_representation:lock.dbnsfp_representation,dbnsfp:lock.dbnsfp[row.chromosome],genebayes:lock.genebayes,container:lock.container]
        tuple(meta,file(row.picked),file(row.loftee))
    }
    SCORE_CANDIDATES(entries,
        Channel.value(file("${projectDir}/scripts/select_scored_candidates.py")),
        Channel.value(file("${projectDir}/scripts/postprocess/tier_variants.py")),
        Channel.value(file("${projectDir}/scripts/postprocess/join_scores.py")))
}

nextflow.enable.dsl=2
include { POST_RARITY_CALL_QC } from './modules/post_rarity_call_qc'
// Call-level QC of gathered chromosomes (scripts/post_rarity_call_qc.py). Manifest: chromosome<TAB>directory (per-chromosome gather outputs).
workflow {
    if (!params.call_qc_manifest) error 'Required: --call_qc_manifest (chromosome<TAB>directory of GATHER_POST_RARITY_QC outputs)'
    def rows = file(params.call_qc_manifest,checkIfExists:true).readLines().findAll { it.trim() }
    if (rows[0] != 'chromosome\tdirectory') error 'Manifest header must be chromosome<TAB>directory'
    def entries = rows.drop(1).collect { it.split('\t',-1) }
    if (!entries || entries.any { it.size()!=2 } || entries.collect { it[0] }.unique().size()!=entries.size()) error 'Invalid call-QC manifest'
    if (params.zarr_dir) file(params.zarr_dir,checkIfExists:true)
    def ch = Channel.fromList(entries).map { c, d ->
        tuple([chromosome:c, zarr_dir:params.zarr_dir, dnm_rescue:params.dnm_rescue ? true : false,
               snv_site_ab_min:params.snv_site_ab_min, indel_site_ab_min:params.indel_site_ab_min], file(d,checkIfExists:true))
    }
    POST_RARITY_CALL_QC(ch, Channel.value(file("${projectDir}/scripts/post_rarity_call_qc.py")),
        Channel.value(['gather_post_qc.py','qc_filtered_carriers.py','qc_post_rarity.py','carrier_frequencies.py','carrier_psam.py','filter_final_rarity.py','extract_exact_carriers.py'].collect { file("${projectDir}/scripts/${it}") }),
        Channel.value(params.dnm_rescue ? file(params.dnm_rescue,checkIfExists:true) : file("${projectDir}/assets/NO_DNM_RESCUE.tsv.gz",checkIfExists:true)))
}

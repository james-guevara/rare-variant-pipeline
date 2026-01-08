#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.resolved_dir = "out_resolve"
params.annotations_dir = "out_reformat"
params.outdir = "out_merge"
params.chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
params.merge_script = "merge_genotypes_annotations.py"

process MERGE_GENOTYPES {
    tag "${chrom}"
    cpus 1
    memory '16 GB'
    time '1h'


    input:
    val chrom

    output:
    path "*.merged.tsv", emit: merged

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate python3.12_env_default
    python3 ${launchDir}/${params.merge_script} \
        --family ${launchDir}/${params.resolved_dir}/${chrom}.resolved_genotypes.tsv \
        --variants ${launchDir}/${params.annotations_dir}/${chrom}.consequential.tsv \
        --out ${chrom}.merged.tsv
    """
}

workflow {
    chroms_ch = Channel.fromList(params.chroms.tokenize(','))
    MERGE_GENOTYPES(chroms_ch)
}

#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.chroms = "chrY"  // comma-separated list, e.g., "chr20,chr21,chr22"
params.n_chunks = 10
params.outdir = "${projectDir}/output"

// Input paths
params.vcf_dir = "/expanse/projects/sebat1/s3/data/sebat/SPARK_iWES_v3/variants/snv/deepvariant/pvcf"
params.vcf_pattern = "SPARK.iWES_v3.2024_08.deepvariant.{chrom}.vcf.gz"

// Container paths (sandbox directories)
params.scripts_dir = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline"
params.bcftools_container = "${params.scripts_dir}/bcftools:1.22--h3a4d415_1"
params.vep_container = "${params.scripts_dir}/ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools"

// Resource paths
params.pipeline_dir = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline"
params.resources_base = "${params.scripts_dir}/resources"
params.vep_cache = "${params.pipeline_dir}/VEP_CACHE"
params.vep_plugins = "${params.scripts_dir}/VEP_PLUGINS_ALL"
params.loftee_path = "${params.scripts_dir}/VEP_PLUGINS/loftee"
params.reference = "/expanse/projects/sebat1/j3guevar/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

// Python scripts
params.reformat_script = "${params.pipeline_dir}/reformat_variants.py"
params.family_query_script = "${params.pipeline_dir}/family_query.py"
params.resolve_script = "${params.pipeline_dir}/resolve_family_genotypes.py"
params.merge_script = "${params.pipeline_dir}/merge_family_variants.py"
params.ped_file = "${params.pipeline_dir}/resources/SPARK_iWES_v3.ped"
params.python_resources = "${params.pipeline_dir}/resources"

// Derive input VCF path
def getInputVcf(chrom) {
    return file("${params.vcf_dir}/${params.vcf_pattern.replace('{chrom}', chrom)}")
}

process DROP_GENOTYPES {
    tag "${chrom}"
    cpus 4
    memory '16 GB'
    time '4h'
    
    input:
    tuple val(chrom), path(vcf), path(tbi)
    
    output:
    val chrom, emit: chrom
    path "*.sites.vcf.gz", emit: vcf
    path "*.sites.vcf.gz.tbi", emit: tbi
    tuple val(chrom), path(vcf), path(tbi), emit: input_vcf

    
    script:
    """
    singularity exec --bind /expanse/projects/sebat1/ ${params.bcftools_container} \\
        bcftools view -G --threads ${task.cpus} -O z -o ${chrom}.sites.vcf.gz ${vcf}
    
    singularity exec --bind /expanse/projects/sebat1/ ${params.bcftools_container} \\
        tabix -p vcf ${chrom}.sites.vcf.gz
    """
}

process VEP_ANNOTATE {
    tag "${chrom}"
    cpus 8
    memory '32 GB'
    time '8h'
    
    input:
    val chrom
    path sites_vcf
    
    output:
    val chrom, emit: chrom
    path "*.vep.vcf.gz", emit: vcf
    path "*.vep.vcf.gz.tbi", emit: tbi

    
    script:
    """
    singularity exec --bind /expanse/projects/sebat1/ ${params.vep_container} vep \\
        --input_file ${sites_vcf} \\
        --format vcf \\
        --output_file ${chrom}.vep.vcf.gz \\
        --vcf \\
        --compress_output bgzip \\
        --minimal --canonical --mane --symbol --protein --pubmed \\
        --nearest symbol --uploaded_allele --numbers \\
        --assembly GRCh38 --cache --dir_cache ${params.vep_cache} --offline \\
        --fasta ${params.reference} \\
        --allele_number --pick_allele --regulatory --biotype --domains \\
        --force_overwrite --fork ${task.cpus} --stats_text \\
        --dir_plugins ${params.vep_plugins} \\
        --plugin dbNSFP,${params.resources_base}/dbNSFP/dbNSFP5.3a_grch38.gz,transcript_match=1,MPC_score,MPC_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,ClinPred_score,ClinPred_rankscore,ClinPred_pred,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,1000Gp3_AC,1000Gp3_AF,AllofUs_ALL_AF,AllofUs_POPMAX_AF,AllofUs_POPMAX_POP,RegeneronME_ALL_AF,gnomAD4.1_joint_flag,gnomAD4.1_joint_AF,gnomAD4.1_joint_nhomalt,gnomAD4.1_joint_POPMAX_AF,gnomAD4.1_joint_POPMAX_nhomalt,ALFA_Total_AF,dbNSFP_POPMAX_AF,dbNSFP_POPMAX_AC,dbNSFP_POPMAX_POP,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,GERP_92_mammals,GERP_92_mammals_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP470way_mammalian,phyloP470way_mammalian_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons470way_mammalian,phastCons470way_mammalian_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore \\
        --plugin LoF,loftee_path:${params.loftee_path},human_ancestor_fa:${params.resources_base}/LOFTEE/human_ancestor.fa.gz,gerp_bigwig:${params.resources_base}/LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:${params.resources_base}/LOFTEE/loftee.sql \\
        --plugin SpliceAI,snv=${params.resources_base}/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=${params.resources_base}/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \\
        --plugin MaxEntScan,${params.resources_base}/MaxEntScan/fordownload,SWA,NCSS
    
    singularity exec --bind /expanse/projects/sebat1/ ${params.vep_container} \\
        tabix -p vcf ${chrom}.vep.vcf.gz
    """
}

process SPLIT_VEP {
    tag "${chrom}"
    cpus 1
    memory "4 GB"
    time "4h"
    
    input:
    val chrom
    path vep_vcf
    
    output:
    val chrom, emit: chrom
    path "*.variants.tsv.gz", emit: tsv
    
    script:
    def bcftools = params.bcftools_container
    """
    singularity exec --bind /expanse/projects/sebat1/ ${bcftools} \
        bcftools +split-vep -p CSQ -HH -d -s :missense+ \
        -f '%CHROM	%POS0	%END	%POS	%REF	%ALT	%ID	%QUAL	%INFO	%CSQ
' \
        -A '	' \
        ${vep_vcf} | \
        sed -E '1s/[[][0-9]+[]]//g' | \
        sed '1s/CSQ//g' | \
        sed '1s/(null)/INFO/g' | \
        bgzip > ${chrom}.variants.tsv.gz
    """
}

process REFORMAT_VARIANTS {
    tag "${chrom}"
    cpus 4
    memory '16 GB'
    time '4h'
    
    input:
    val chrom
    path variants_tsv
    
    output:
    val chrom, emit: chrom
    path "*.reformatted.tsv.gz", emit: reformatted
    path "*.bed", emit: bed
    tuple val(chrom), path("*.consequential.tsv.gz"), emit: conseq_tsv
    path "*.consequential.bed", emit: conseq_bed
    
    script:
    """
    micromamba run -n python3.12_env_default python ${params.reformat_script} \\
        ${variants_tsv} ${chrom}.reformatted.tsv.gz \\
        --bed ${chrom}.bed \\
        --consequential ${chrom}.consequential.tsv.gz \\
        --consequential-bed ${chrom}.consequential.bed \\
        --resources-dir ${params.python_resources}
    """
}

process SCATTER_BED {
    tag "${chrom}"
    cpus 1
    memory '1 GB'
    time '5m'
    
    input:
    val chrom
    path conseq_bed
    
    output:
    tuple val(chrom), path("chunk_*.bed"), emit: chunks
    
    script:
    """
    total=\$(wc -l < ${conseq_bed})
    chunk_size=\$(( (total + ${params.n_chunks} - 1) / ${params.n_chunks} ))
    
    for i in \$(seq 1 ${params.n_chunks}); do
        start=\$(( (i - 1) * chunk_size + 1 ))
        end=\$(( i * chunk_size ))
        sed -n "\${start},\${end}p" ${conseq_bed} > chunk_\${i}.bed || touch chunk_\${i}.bed
    done
    """
}

process FAMILY_QUERY {
    tag "${chrom}_chunk${chunk_id}"
    cpus 1
    memory '16 GB'
    time '4h'
    
    input:
    tuple val(chrom), val(chunk_id), path(bed_chunk), path(input_vcf), path(input_tbi)
    
    output:
    tuple val(chrom), path("*.genotypes.tsv.gz"), emit: genotypes
    
    script:
    """
    if [ -s ${bed_chunk} ]; then
        micromamba run -n python3.12_env_default python ${params.family_query_script} \\
            --vcf ${input_vcf} \\
            --ped ${params.ped_file} \\
            --region ${bed_chunk} \\
            --out /dev/stdout | gzip > chunk_${chunk_id}.genotypes.tsv.gz
    else
        echo -e "#CHROM\\tPOS0\\tEND\\tREF\\tALT\\tFID\\tIID\\tGT\\tGQ\\tDP\\tAD" | gzip > chunk_${chunk_id}.genotypes.tsv.gz
    fi
    """
}

process GATHER_GENOTYPES {
    tag "${chrom}"
    cpus 1
    memory '4 GB'
    time '10m'
    
    input:
    tuple val(chrom), path(genotype_chunks)
    
    output:
    tuple val(chrom), path("*.family_genotypes.tsv.gz"), emit: genotypes
    
    script:
    """
    set +o pipefail
    
    zcat \$(ls chunk_*.genotypes.tsv.gz | sort -V | head -1) | head -1 | gzip > ${chrom}.family_genotypes.tsv.gz
    
    for f in \$(ls chunk_*.genotypes.tsv.gz | sort -V); do
        zcat "\$f" | tail -n +2
    done | gzip >> ${chrom}.family_genotypes.tsv.gz
    """
}

process RESOLVE_GENOTYPES {
    tag "${chrom}"
    cpus 1
    memory '16 GB'
    time '4h'
    
    input:
    tuple val(chrom), path(family_genotypes)
    
    output:
    tuple val(chrom), path("*.resolved_genotypes.tsv.gz"), emit: resolved
    
    script:
    """
    micromamba run -n python3.12_env_default python ${params.resolve_script} \\
        ${family_genotypes} ${chrom}.resolved_genotypes.tsv.gz
    """
}

process MERGE_VARIANTS {
    tag "${chrom}"
    cpus 1
    memory '16 GB'
    time '4h'
    publishDir "${params.outdir}/merged", mode: 'copy'
    
    input:
    tuple val(chrom), path(resolved_genotypes), path(conseq_tsv)
    
    output:
    path "*.merged.tsv.gz", emit: merged
    
    script:
    """
    micromamba run -n python3.12_env_default python ${params.merge_script} \\
        --family ${resolved_genotypes} \\
        --variants ${conseq_tsv} \\
        --out ${chrom}.merged.tsv.gz
    """
}
workflow {
    // Create channel from comma-separated chromosome list
    chroms_ch = Channel.fromList(params.chroms.tokenize(','))
    
    // Create input channel with chrom, vcf path, and tbi path
    input_ch = chroms_ch.map { chrom -> 
        def vcf_path = params.vcf_pattern.replace('{chrom}', chrom)
        def vcf_file = file("${params.vcf_dir}/${vcf_path}")
        def tbi_file = file("${params.vcf_dir}/${vcf_path}.tbi")
        return tuple(chrom, vcf_file, tbi_file)
    }
    
    DROP_GENOTYPES(input_ch)
    VEP_ANNOTATE(DROP_GENOTYPES.out.chrom, DROP_GENOTYPES.out.vcf)
    SPLIT_VEP(VEP_ANNOTATE.out.chrom, VEP_ANNOTATE.out.vcf)
    REFORMAT_VARIANTS(SPLIT_VEP.out.chrom, SPLIT_VEP.out.tsv)
    SCATTER_BED(REFORMAT_VARIANTS.out.chrom, REFORMAT_VARIANTS.out.conseq_bed)
    
    // Join scatter output with input VCF+TBI (propagated from DROP_GENOTYPES)
    vcf_by_chrom = DROP_GENOTYPES.out.input_vcf
    
    // Create chunks with VCF+TBI for family query
    chunks_with_vcf = SCATTER_BED.out.chunks
        .join(vcf_by_chrom)
        .flatMap { chrom, chunks, vcf, tbi ->
            chunks.collect { chunk_file ->
                def chunk_id = chunk_file.name.replace('chunk_', '').replace('.bed', '')
                return tuple(chrom, chunk_id, chunk_file, vcf, tbi)
            }
        }
    
    FAMILY_QUERY(chunks_with_vcf)
    
    // Group by chrom for gather
    gathered = FAMILY_QUERY.out.genotypes.groupTuple(size: params.n_chunks)
    
    GATHER_GENOTYPES(gathered)
    RESOLVE_GENOTYPES(GATHER_GENOTYPES.out.genotypes)
    
    // Combine resolved with conseq_tsv by chrom
    merge_input = RESOLVE_GENOTYPES.out.resolved.join(REFORMAT_VARIANTS.out.conseq_tsv)
    
    MERGE_VARIANTS(merge_input)
}

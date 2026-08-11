// Module: Normalize VCF — split multi-allelics and left-align indels

process BCFTOOLS_NORM {
    tag "${chrom}"
    cpus 4
    memory '16 GB'
    time '48h'
    container "${params.bcftools_container}"

    input:
    tuple val(chrom), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("${chrom}.norm.vcf.gz"), path("${chrom}.norm.vcf.gz.tbi"), emit: norm

    script:
    // Stages are assembled rather than hardcoded because callers differ in what
    // they store per sample. One optional stage, off by default:
    //
    // params.local_alleles — some callers store per-sample depths and likelihoods
    //   as LOCAL-allele fields (LAD/LPL/LAA) rather than AD/PL. Those are Number=.,
    //   so `norm -m-` has no subsetting rule and copies them VERBATIM onto every
    //   split row, silently. The LAA indices then no longer describe the row they
    //   sit on, and a non-carrier can inherit a carrier's depth. The decode
    //   therefore MUST precede the split. This is bcftools behaviour, not a bug in a
    //   particular version, so it will not go away with an upgrade. `-r` drops the
    //   local tags once decoded so nothing downstream can read a stale index.
    //   Only the family form --LXX-to-XX works; the single-tag --LAD-to-AD is
    //   rejected despite appearing in tag2tag's own --help.
    //
    def stages = []
    stages << "bcftools view -r ${chrom} --threads ${task.cpus} -O u ${vcf}"
    if (params.local_alleles) {
        stages << "bcftools +tag2tag -O u -- --LXX-to-XX -r"
        // tag2tag -r removes LAD/LPL/LAA, but NOT the other local-allele-indexed
        // fields. LF1R2/LF2R1 ("...supporting each local allele") are Number=. and
        // so get copied verbatim onto every split row exactly like LAD did; LAF is
        // header-declared but absent from records. Nothing downstream reads them
        // today, but leaving allele-indexed fields whose indices no longer match
        // their row is how the original bug happened, so strip them here.
        stages << "bcftools annotate -x FORMAT/LAF,FORMAT/LF1R2,FORMAT/LF2R1 -O u"
    }
    stages << "bcftools norm -m- --force -f ${params.reference} --threads ${task.cpus} -O u"
    // AN/AC/AF are always recomputed from GT, and F_MISSING always added. Both are
    // deliberate and unconditional:
    //   * AF: downstream filters on INFO/AF, so if it were the caller's own value the
    //     filter would mean something different per cohort — different callers compute
    //     it over different sample sets and with different ploidy handling, and it can
    //     be stale (e.g. after samples are removed post-joint-calling). Recomputing
    //     makes it mean "frequency among the samples in this file", everywhere.
    //   * F_MISSING: per-variant missing call rate, which the pipeline otherwise has
    //     no measure of at all.
    // Genotypes are already streaming past here, so neither costs an extra pass.
    stages << "bcftools +fill-tags -O z -o ${chrom}.norm.vcf.gz -- -t AN,AC,AF,F_MISSING"
    // pipefail is essential here: Nextflow runs .command.sh under `bash -ue`,
    // which does NOT fail a pipeline when a NON-final stage dies. Without it a
    // crashed tag2tag/norm would leave the last stage exiting 0 on truncated
    // input -> a silently short normed VCF. Continuations are indented to match
    // the surrounding block so Nextflow's stripIndent() removes them uniformly.
    """
    set -o pipefail
    ${stages.join(" \\\n        | ")}
    bcftools index -t --threads ${task.cpus} ${chrom}.norm.vcf.gz
    """
}

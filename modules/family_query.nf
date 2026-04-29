// Module: Query carrier genotypes from a normed cohort VCF.
//
// Replaces the previous cyvcf2-based family-expansion implementation.
// Outputs ONE row per (variant × carrier sample) — non-carrier family members
// are NOT emitted. De novo / pedigree-QC use cases that need non-carrier
// parent rows should run a separate pass on the normed VCF for the relevant
// sites; for rare-variant burden testing, carriers are sufficient.

process FAMILY_QUERY {
    tag "${chrom}_${chunk_bed.baseName}"
    cpus 1
    memory '4 GB'
    time '4h'
    container "${params.bcftools_container}"

    input:
    tuple val(chrom), path(chunk_bed), path(vcf), path(tbi)

    output:
    tuple val(chrom), path("${chunk_bed.baseName}.carriers.tsv"), emit: carriers

    script:
    """
    # Header for downstream consumers
    printf '#CHROM\\tPOS0\\tEND\\tPOS\\tREF\\tALT\\tSAMPLE\\tGT\\tGQ\\tDP\\tAD_ref\\tAD_alt\\tFAMILY\\n' > ${chunk_bed.baseName}.carriers.tsv

    # bcftools query: emit one line per (record × sample), then awk filters to
    # carrier samples, splits AD, joins FAMILY from the PED file.
    bcftools query \\
        -R ${chunk_bed} \\
        -i 'GT[*]~"[1-9]"' \\
        -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\t%AD\\t%GQ\\t%DP\\n]' \\
        ${vcf} | \\
    awk -F'\\t' -v OFS='\\t' -v ped='${params.ped_file}' '
        BEGIN {
            while ((getline line < ped) > 0) {
                if (line == "" || substr(line,1,1) == "#") continue
                n = split(line, f, /[ \\t]+/)
                if (n < 2) continue
                lf = tolower(f[1]); ls = tolower(f[2])
                if ((lf == "fam" || lf == "famid" || lf == "family" || lf == "familyid") && (ls == "iid" || ls == "sample" || ls == "sampleid")) continue
                fam[f[2]] = f[1]
            }
            close(ped)
        }
        {
            # \$1=CHROM \$2=POS \$3=REF \$4=ALT \$5=SAMPLE \$6=GT \$7=AD \$8=GQ \$9=DP
            if (\$6 !~ /[1-9]/) next
            pos0 = \$2 - 1
            end_ = pos0 + length(\$3)
            ad_ref = "."; ad_alt = "."
            if (\$7 != "" && \$7 != ".") {
                m = split(\$7, ad_parts, ",")
                if (m >= 1) ad_ref = ad_parts[1]
                if (m >= 2) ad_alt = ad_parts[2]
            }
            family_id = (\$5 in fam) ? fam[\$5] : "."
            print \$1, pos0, end_, \$2, \$3, \$4, \$5, \$6, \$8, \$9, ad_ref, ad_alt, family_id
        }
    ' >> ${chunk_bed.baseName}.carriers.tsv
    """
}

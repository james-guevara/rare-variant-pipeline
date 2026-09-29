# Validate and hash only variant columns; header provenance lines are excluded.
BEGIN {
    FS = OFS = "\t"
    count = 0; csq = 0; header = 0
    sub(/^chr/, "", chromosome)
    if (chromosome == "M") chromosome = "MT"
    printf "" > (prefix ".csq_header")
}
/^##INFO=<ID=CSQ[,>]/ {
    csq = 1
    print > (prefix ".csq_header")
}
/^#CHROM\t/ {
    if (NF != 8) { print "Sites header must have exactly eight columns" > "/dev/stderr"; exit 1 }
    header = 1
}
/^#/ { next }
{
    if (!header || NF != 8) { print "Malformed sites record or FORMAT/sample columns present" > "/dev/stderr"; exit 1 }
    observed = $1; sub(/^chr/, "", observed)
    if (observed == "M") observed = "MT"
    if (observed != chromosome) { print "Record chromosome disagrees with manifest" > "/dev/stderr"; exit 1 }
    count++
    print
}
END {
    if (!header) { print "Missing sites VCF header" > "/dev/stderr"; exit 1 }
    print count > (prefix ".records")
    print csq > (prefix ".csq")
}

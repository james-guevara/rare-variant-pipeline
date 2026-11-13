"""Pure Python port of the pinned LOFTEE decision layer.

Compatibility target: konradjk/loftee commit
a46b502a68c812c8ae0c5a5721c0603fe81cae8d with its default GRCh38 settings.
Resource-dependent calculations belong in the context builder, not here.
"""

from __future__ import annotations

from .model import Classification, TranscriptContext


OTHER_LOF = frozenset({"stop_gained", "frameshift_variant"})
SPLICE_LOF = frozenset({"splice_acceptor_variant", "splice_donor_variant"})


def _fmt(value: float | int) -> str:
    # Perl's default scalar stringification uses approximately 15 significant
    # digits. Match it so LoF_info remains byte-comparable to the plugin output.
    return format(value, ".15g") if isinstance(value, float) else str(value)


def classify(
    context: TranscriptContext,
    *,
    min_intron_size: int = 15,
    gerp_end_trunc_cutoff: float = -58,
    use_gerp_end_trunc: bool = True,
    filter_position: float = 0.05,
    check_complete_cds: bool = False,
    use_conservation: bool = True,
    use_human_ancestor: bool = True,
) -> Classification:
    """Classify one variant-transcript-allele context.

    Filter, flag, and info ordering intentionally follows LoF.pm because those
    strings are part of the chr22 parity contract.
    """

    if context.biotype != "protein_coding":
        return Classification()

    consequences = set(context.consequences)
    other_lof = bool(consequences & OTHER_LOF)
    splice_lof = bool(consequences & SPLICE_LOF)
    if not other_lof and not splice_lof:
        return Classification()

    filters: list[str] = []
    flags: list[str] = []
    info: list[str] = []

    # The pinned configuration has run_splice_predictions=0, so VEP splice
    # donor/acceptor consequences enter as HC before structural filters.
    confidence = "HC"

    if other_lof:
        if context.percentile is not None:
            info.append(f"PERCENTILE:{_fmt(context.percentile)}")
        if use_gerp_end_trunc:
            if context.gerp_weighted_distance is not None:
                info.append(f"GERP_DIST:{_fmt(context.gerp_weighted_distance)}")
            if context.distance_to_stop is not None:
                info.append(f"BP_DIST:{context.distance_to_stop}")
            if context.exon_number is None:
                flags.append("NO_EXON_NUMBER")
            elif (
                context.distance_to_stop is not None
                and context.last_exon_coding_length is not None
                and context.gerp_weighted_distance is not None
            ):
                distance_from_last_exon = (
                    context.distance_to_stop - context.last_exon_coding_length
                )
                info.append(f"DIST_FROM_LAST_EXON:{distance_from_last_exon}")
                info.append(
                    "50_BP_RULE:"
                    + ("FAIL" if distance_from_last_exon <= 50 else "PASS")
                )
                if (
                    distance_from_last_exon <= 50
                    and context.gerp_weighted_distance <= gerp_end_trunc_cutoff
                ):
                    filters.append("END_TRUNC")
        elif context.percentile is not None and context.percentile >= 1 - filter_position:
            filters.append("END_TRUNC")

    if other_lof and context.exon_number is not None:
        if context.exon_intron_undefined:
            filters.append("EXON_INTRON_UNDEF")
        elif context.exon_number[1] == 1:
            flags.append("SINGLE_EXON")
        elif check_complete_cds and context.incomplete_cds:
            filters.append("INCOMPLETE_CDS")

        if use_conservation:
            if context.conservation_too_short:
                info.append("PHYLOCSF_TOO_SHORT")
            elif (
                context.conservation_orf_score is not None
                and context.conservation_max_score is not None
            ):
                info.append(f"ANN_ORF:{_fmt(context.conservation_orf_score)}")
                info.append(f"MAX_ORF:{_fmt(context.conservation_max_score)}")
                if context.conservation_orf_score < 0:
                    flags.append(
                        "PHYLOCSF_UNLIKELY_ORF"
                        if context.conservation_max_score > 0
                        else "PHYLOCSF_WEAK"
                    )

    if context.intron_number is not None:
        if context.exon_intron_undefined:
            filters.append("EXON_INTRON_UNDEF")
        else:
            if context.intron_size is not None:
                info.append(f"INTRON_SIZE:{context.intron_size}")
                if context.intron_size < min_intron_size:
                    filters.append("SMALL_INTRON")
            if splice_lof:
                if "splice_donor_variant" in consequences and context.gc_to_gt_donor:
                    filters.append("GC_TO_GT_DONOR")
                if context.noncanonical_intron:
                    flags.append("NON_CAN_SPLICE")
                # VEP's explicit consequence is authoritative for indels that
                # span a CDS boundary; coordinate-only inference can place the
                # event inside the CDS even though its splice consequence is UTR.
                if context.five_prime_utr or "5_prime_UTR_variant" in consequences:
                    filters.append("5UTR_SPLICE")
                if context.three_prime_utr or "3_prime_UTR_variant" in consequences:
                    filters.append("3UTR_SPLICE")
            if "splice_acceptor_variant" in consequences and context.nagnag_site:
                flags.append("NAGNAG_SITE")

    # Some VEP splice/UTR consequences do not carry an INTRON value. Preserve
    # LOFTEE's UTR filter even when structural intron context is unavailable.
    if splice_lof and "5_prime_UTR_variant" in consequences and "5UTR_SPLICE" not in filters:
        filters.append("5UTR_SPLICE")
    if splice_lof and "3_prime_UTR_variant" in consequences and "3UTR_SPLICE" not in filters:
        filters.append("3UTR_SPLICE")

    if use_human_ancestor and context.ancestral_allele:
        filters.append("ANC_ALLELE")

    if filters:
        confidence = "LC"
    return Classification(confidence, tuple(filters), tuple(flags), tuple(info))

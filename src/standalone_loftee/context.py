"""Derive classifier inputs from base-VEP fields and pinned resources."""

from __future__ import annotations

import re

from .model import TranscriptContext
from .resources import LofteeResources
from .transcripts import Transcript


COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(COMPLEMENT)[::-1]


def parse_number(value: str) -> tuple[int, int] | None:
    if not value or value == ".":
        return None
    match = re.fullmatch(r"(\d+)/(\d+)", value)
    return (int(match[1]), int(match[2])) if match else None


def cds_coordinate(value: str, *, insertion: bool) -> int | None:
    if not value or value == ".":
        return None
    numerator = value.split("/", 1)[0]
    positions = numerator.split("-", 1)
    selected = positions[0] if insertion else positions[-1]
    return int(selected) if selected.isdigit() else None


def feature_interval(position: int, ref: str, alt: str) -> tuple[int, int]:
    """Convert normalized VCF alleles to VEP's affected-feature interval."""
    prefix = 0
    while prefix < len(ref) and prefix < len(alt) and ref[prefix] == alt[prefix]:
        prefix += 1
    ref_remainder = ref[prefix:]
    start = position + prefix
    end = start + max(len(ref_remainder), 1) - 1
    return start, end


def deletion_shift_maps_to_exon(
    transcript: Transcript, start: int, end: int, hgvs_offset: int
) -> bool:
    direction = transcript.strand
    shifted_start = start + direction * hgvs_offset
    shifted_end = end + direction * hgvs_offset
    shifted_start, shifted_end = sorted((shifted_start, shifted_end))
    return any(
        exon.start <= shifted_start and shifted_end <= exon.end
        for exon in transcript.exons
    )


def build_context(
    row: dict[str, str], transcript: Transcript, resources: LofteeResources
) -> TranscriptContext:
    consequences = tuple(row["Consequence"].split("&"))
    position, end = feature_interval(int(row["POS"]), row["REF"], row["ALT"])
    exon_number = parse_number(row.get("EXON", ""))
    intron_number = parse_number(row.get("INTRON", ""))
    insertion = len(row["ALT"]) > len(row["REF"])
    coding_position = cds_coordinate(row.get("CDS_position", ""), insertion=insertion)
    hgvs_offset_value = row.get("HGVS_OFFSET", "")
    hgvs_offset = int(hgvs_offset_value) if hgvs_offset_value not in {"", "."} else 0
    shifted_maps_to_exon = False
    if hgvs_offset and not insertion:
        shifted_maps_to_exon = deletion_shift_maps_to_exon(
            transcript, position, end, hgvs_offset
        )
    if coding_position is not None and (insertion or shifted_maps_to_exon):
        # VEP exposes the 3-prime HGVS repeat shift separately from its printed
        # CDS_position. Shifted deletions that land in a mapper gap fall back to
        # their original coordinates; exon-contained shifts retain the offset.
        coding_position += hgvs_offset

    percentile = coding_position / transcript.cds_length if coding_position and transcript.cds_length else None
    gerp_distance = bp_distance = last_exon_length = None
    if consequences and {"stop_gained", "frameshift_variant"}.intersection(consequences):
        gerp_distance, bp_distance = resources.gerp_weighted_distance(transcript, position)
        last_exon_length = transcript.last_exon_coding_length()

    intron_size = None
    intron_sequence = ""
    if intron_number is not None and 1 <= intron_number[0] <= len(transcript.introns):
        intron = transcript.introns[intron_number[0] - 1]
        intron_size = intron.length
        intron_sequence = resources.sequence(transcript.seqname, intron.start, intron.end)
        if transcript.strand == -1:
            intron_sequence = reverse_complement(intron_sequence)

    ref, alt = row["REF"], row["ALT"]
    if transcript.strand == -1:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    gc_to_gt = (
        "splice_donor_variant" in consequences
        and intron_sequence.startswith("GC")
        and ref == "C"
        and alt == "T"
    )
    noncanonical = bool(intron_sequence) and not (
        intron_sequence.startswith("GT") and intron_sequence.endswith("AG")
    )

    nagnag = False
    if "splice_acceptor_variant" in consequences and len(row["REF"]) == len(row["ALT"]) == 1:
        splice_context = resources.sequence(transcript.seqname, position - 4, end + 4)
        if transcript.strand == -1:
            splice_context = reverse_complement(splice_context)
        nagnag = re.search(r"AG.AG", splice_context) is not None

    five_utr = (
        end < transcript.coding_start if transcript.strand == 1 else position > transcript.coding_end
    ) if transcript.coding_start is not None and transcript.coding_end is not None else False
    three_utr = (
        position > transcript.coding_end if transcript.strand == 1 else end < transcript.coding_start
    ) if transcript.coding_start is not None and transcript.coding_end is not None else False

    ancestral = False
    if len(row["REF"]) == len(row["ALT"]) == 1:
        ancestral = resources.ancestral_allele(transcript.seqname, position, end) == row["ALT"].upper()

    conservation = resources.phylocsf(transcript.transcript_id, exon_number[0]) if exon_number else None
    return TranscriptContext(
        consequences=consequences,
        biotype=row.get("BIOTYPE", transcript.biotype),
        exon_number=exon_number,
        intron_number=intron_number,
        five_prime_utr=five_utr,
        three_prime_utr=three_utr,
        intron_size=intron_size,
        noncanonical_intron=noncanonical,
        gc_to_gt_donor=gc_to_gt,
        nagnag_site=nagnag,
        ancestral_allele=ancestral,
        gerp_weighted_distance=gerp_distance,
        distance_to_stop=bp_distance,
        last_exon_coding_length=last_exon_length,
        percentile=percentile,
        conservation_orf_score=(conservation.corresponding_orf_score if conservation else None),
        conservation_max_score=(conservation.max_score if conservation else None),
        conservation_too_short=exon_number is not None and conservation is None,
    )

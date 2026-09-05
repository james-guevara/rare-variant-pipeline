"""Typed boundary between transcript/resource lookups and LOFTEE rules."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class TranscriptContext:
    """All information consumed by the pure classification layer.

    A separate context builder is responsible for deriving these values from
    base-VEP output, the pinned Ensembl transcript model, FASTA, GERP, the
    ancestral sequence, and the LOFTEE conservation database.
    """

    consequences: tuple[str, ...]
    biotype: str = "protein_coding"
    exon_number: tuple[int, int] | None = None
    intron_number: tuple[int, int] | None = None
    five_prime_utr: bool = False
    three_prime_utr: bool = False
    exon_intron_undefined: bool = False
    incomplete_cds: bool = False
    intron_size: int | None = None
    noncanonical_intron: bool = False
    gc_to_gt_donor: bool = False
    nagnag_site: bool = False
    ancestral_allele: bool = False
    gerp_weighted_distance: float | None = None
    distance_to_stop: int | None = None
    last_exon_coding_length: int | None = None
    percentile: float | None = None
    conservation_orf_score: float | None = None
    conservation_max_score: float | None = None
    conservation_too_short: bool = False


@dataclass(frozen=True)
class Classification:
    lof: str = ""
    filters: tuple[str, ...] = field(default_factory=tuple)
    flags: tuple[str, ...] = field(default_factory=tuple)
    info: tuple[str, ...] = field(default_factory=tuple)

    @property
    def lof_filter(self) -> str:
        return ",".join(self.filters)

    @property
    def lof_flags(self) -> str:
        return ",".join(self.flags)

    @property
    def lof_info(self) -> str:
        return ",".join(self.info)

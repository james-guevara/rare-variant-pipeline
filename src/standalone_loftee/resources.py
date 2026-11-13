"""Version-locked sequence, GERP, and PhyloCSF resource access."""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .transcripts import Transcript


@dataclass(frozen=True)
class Conservation:
    corresponding_orf_score: float
    max_score: float


class LofteeResources:
    def __init__(
        self,
        reference_fasta: Path,
        ancestor_fasta: Path,
        gerp_bigwig: Path,
        conservation_db: Path,
    ):
        # Imported lazily so the pure classifier and mocked resource tests do
        # not require platform-specific native extensions at import time.
        import pyBigWig
        import pysam

        self.reference = pysam.FastaFile(str(reference_fasta))
        self.ancestor = pysam.FastaFile(str(ancestor_fasta))
        self.gerp = pyBigWig.open(str(gerp_bigwig))
        self.conservation = sqlite3.connect(conservation_db)

    def close(self) -> None:
        self.reference.close()
        self.ancestor.close()
        self.gerp.close()
        self.conservation.close()

    @lru_cache(maxsize=65_536)
    def ancestral_allele(self, chrom: str, start: int, end: int) -> str:
        """Return uppercase ancestral bases for a 1-based inclusive interval."""
        contig = chrom.removeprefix("chr")
        return self.ancestor.fetch(contig, start - 1, end).upper()

    @lru_cache(maxsize=65_536)
    def sequence(self, chrom: str, start: int, end: int) -> str:
        contig = chrom.removeprefix("chr")
        return self.reference.fetch(contig, start - 1, end).upper()

    @lru_cache(maxsize=262_144)
    def interval_gerp(self, chrom: str, start: int, end: int) -> float:
        """Sum BigWig scores over a 1-based inclusive interval.

        Missing BigWig positions contribute zero, matching the Perl summary
        feature loop, which only accumulates covered features.
        """
        contig = chrom.removeprefix("chr")
        left, right = sorted((start, end))
        # Bio::DB::BigWig's `features(type => 'summary')` uses the BigWig zoom
        # summary rather than exact base values. pyBigWig's default stats call
        # uses the same summary path; mean * covered interval length reproduces
        # LOFTEE's length * binMean(score).
        mean = self.gerp.stats(contig, left - 1, right, type="mean", exact=False)[0]
        return 0.0 if mean is None else mean * (right - left + 1)

    @lru_cache(maxsize=65_536)
    def gerp_weighted_distance(
        self, transcript: Transcript, variant_position: int
    ) -> tuple[float, int]:
        """Port gerp_dist.pl::get_gerp_weighted_dist."""
        stop = transcript.coding_end if transcript.strand == 1 else transcript.coding_start
        if stop is None:
            raise ValueError(f"transcript {transcript.transcript_id} has no CDS")
        weighted_distance = 0.0
        distance = 0
        for exon in transcript.exons:
            if transcript.strand == -1 and variant_position < exon.start:
                continue
            if transcript.strand == 1 and variant_position > exon.end:
                continue
            last_exon = exon.start < stop <= exon.end if transcript.strand == 1 else exon.start <= stop < exon.end
            in_affected_exon = exon.start <= variant_position <= exon.end
            if last_exon:
                start = variant_position if in_affected_exon else (exon.start if transcript.strand == 1 else exon.end)
                end = stop
            elif in_affected_exon:
                start = variant_position
                end = exon.end if transcript.strand == 1 else exon.start
            else:
                start, end = exon.start, exon.end
            weighted_distance += self.interval_gerp(transcript.seqname, start, end)
            distance += abs(end - start)
        return weighted_distance, distance

    @lru_cache(maxsize=65_536)
    def phylocsf(self, transcript_id: str, exon_number: int) -> Conservation | None:
        stable_id = transcript_id.split(".", 1)[0]
        row = self.conservation.execute(
            "SELECT corresponding_orf_score, max_score FROM phylocsf_summary "
            "WHERE transcript = ? AND exon = ?",
            (stable_id, exon_number),
        ).fetchone()
        return Conservation(*map(float, row)) if row is not None else None

    def __enter__(self) -> "LofteeResources":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

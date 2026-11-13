"""Build and query a compact transcript model from a version-locked GTF."""

from __future__ import annotations

import gzip
import re
import sqlite3
from functools import lru_cache
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


ATTRIBUTE_RE = re.compile(r'([^ ]+) "([^"]*)";')


@dataclass(frozen=True)
class Interval:
    start: int
    end: int
    rank: int | None = None
    phase: int | None = None

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    gene_id: str
    seqname: str
    strand: int
    biotype: str
    exons: tuple[Interval, ...]
    cds: tuple[Interval, ...]
    start_codon: tuple[Interval, ...]
    stop_codon: tuple[Interval, ...]

    @property
    def coding_start(self) -> int | None:
        return min((part.start for part in self.cds + self.stop_codon), default=None)

    @property
    def coding_end(self) -> int | None:
        return max((part.end for part in self.cds + self.stop_codon), default=None)

    @property
    def cds_length(self) -> int:
        # Ensembl GTF CDS features exclude the stop codon, while the VEP API's
        # cdna_coding_start/end span used by LOFTEE includes it.
        return sum(part.length for part in self.cds + self.stop_codon)

    @property
    def introns(self) -> tuple[Interval, ...]:
        genomic = sorted(self.exons, key=lambda exon: exon.start)
        introns = tuple(
            Interval(left.end + 1, right.start - 1, rank=index + 1)
            for index, (left, right) in enumerate(zip(genomic, genomic[1:]))
        )
        return introns if self.strand == 1 else tuple(reversed(introns))

    def last_exon_coding_length(self) -> int | None:
        """Match LoF.pm get_last_exon_coding_length (exclusive length)."""
        stop = self.coding_end if self.strand == 1 else self.coding_start
        if stop is None:
            return None
        for exon in reversed(self.exons):
            if exon.start <= stop <= exon.end:
                return stop - exon.start if self.strand == 1 else exon.end - stop
        return None


def parse_attributes(value: str) -> dict[str, str]:
    return dict(ATTRIBUTE_RE.findall(value))


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def build_transcript_db(gtf: Path, output: Path, *, seqname: str | None = None) -> None:
    """Create a deterministic SQLite transcript database from an Ensembl GTF."""
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        output.unlink()
    connection = sqlite3.connect(output)
    try:
        connection.executescript(
            """
            PRAGMA journal_mode = OFF;
            PRAGMA synchronous = OFF;
            CREATE TABLE transcript (
                transcript_id TEXT PRIMARY KEY,
                gene_id TEXT NOT NULL,
                seqname TEXT NOT NULL,
                strand INTEGER NOT NULL,
                biotype TEXT NOT NULL
            );
            CREATE TABLE feature (
                transcript_id TEXT NOT NULL,
                feature TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL,
                rank INTEGER,
                phase INTEGER,
                PRIMARY KEY (transcript_id, feature, start, end)
            );
            CREATE INDEX feature_transcript_idx ON feature(transcript_id, feature);
            """
        )
        transcripts: dict[str, tuple[str, str, int, str]] = {}
        features: list[tuple[str, str, int, int, int | None, int | None]] = []

        def flush_features() -> None:
            if not features:
                return
            connection.executemany(
                "INSERT OR IGNORE INTO feature VALUES (?, ?, ?, ?, ?, ?)",
                features,
            )
            features.clear()
        with _open_text(gtf) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"invalid GTF row with {len(fields)} fields")
                row_seqname, _, feature, start, end, _, strand, phase, raw_attrs = fields
                if seqname is not None and row_seqname != seqname:
                    continue
                attrs = parse_attributes(raw_attrs)
                transcript_id = attrs.get("transcript_id")
                if transcript_id is None:
                    continue
                gene_id = attrs["gene_id"]
                biotype = attrs.get("transcript_biotype", attrs.get("gene_biotype", ""))
                strand_value = 1 if strand == "+" else -1
                transcripts[transcript_id] = (gene_id, row_seqname, strand_value, biotype)
                if feature in {"exon", "CDS", "start_codon", "stop_codon"}:
                    rank = int(attrs["exon_number"]) if attrs.get("exon_number", "").isdigit() else None
                    phase_value = int(phase) if phase in {"0", "1", "2"} else None
                    features.append(
                        (transcript_id, feature, int(start), int(end), rank, phase_value)
                    )
                    if len(features) >= 50_000:
                        flush_features()
        flush_features()
        connection.executemany(
            "INSERT INTO transcript VALUES (?, ?, ?, ?, ?)",
            ((tid, *values) for tid, values in sorted(transcripts.items())),
        )
        connection.commit()
    finally:
        connection.close()


class TranscriptStore:
    def __init__(self, path: Path):
        self.connection = sqlite3.connect(path)

    def close(self) -> None:
        self.connection.close()

    @lru_cache(maxsize=16_384)
    def get(self, transcript_id: str) -> Transcript | None:
        stable_id = transcript_id.split(".", 1)[0]
        row = self.connection.execute(
            "SELECT transcript_id, gene_id, seqname, strand, biotype "
            "FROM transcript WHERE transcript_id = ?",
            (stable_id,),
        ).fetchone()
        if row is None:
            return None
        feature_rows = self.connection.execute(
            "SELECT feature, start, end, rank, phase FROM feature "
            "WHERE transcript_id = ? ORDER BY start, end",
            (stable_id,),
        ).fetchall()
        grouped: dict[str, list[Interval]] = {
            "exon": [], "CDS": [], "start_codon": [], "stop_codon": []
        }
        for feature, start, end, rank, phase in feature_rows:
            grouped[feature].append(Interval(start, end, rank, phase))
        order = lambda item: (item.rank if item.rank is not None else item.start)
        for values in grouped.values():
            values.sort(key=order)
        return Transcript(*row, *(tuple(grouped[name]) for name in grouped))

    def __enter__(self) -> "TranscriptStore":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

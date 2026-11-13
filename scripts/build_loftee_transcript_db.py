#!/usr/bin/env python3
"""Build a compact SQLite transcript model for standalone LOFTEE."""

from __future__ import annotations

import argparse
from pathlib import Path

from standalone_loftee.transcripts import build_transcript_db


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seqname", help="restrict to one GTF sequence name, e.g. 22")
    args = parser.parse_args()
    build_transcript_db(args.gtf, args.output, seqname=args.seqname)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

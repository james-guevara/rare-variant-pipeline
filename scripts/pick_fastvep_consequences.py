#!/usr/bin/env python3
"""Apply Ensembl VEP's default --pick_allele hierarchy to fastVEP CSQ rows."""

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path


CSQ_FORMAT = re.compile(r'Format: ([^"]+)')


def open_text(path):
    if str(path) == "-":
        return sys.stdin
    return gzip.open(str(path), "rt") if path.suffix == ".gz" else path.open()


def bare_transcript(value):
    return value.split(".", 1)[0]


def iter_csq_records(path):
    """Yield one VCF record and its decoded CSQ rows without retaining the VCF."""
    fields = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = CSQ_FORMAT.search(line)
                if not match:
                    raise ValueError("Could not parse CSQ schema in {}".format(path))
                fields = match.group(1).split("|")
                continue
            if line.startswith("#"):
                continue
            if fields is None:
                raise ValueError("No CSQ schema found before records in {}".format(path))

            columns = line.rstrip("\n").split("\t")
            csq_value = next(
                (value[4:] for value in columns[7].split(";")
                 if value.startswith("CSQ=")),
                None,
            )
            if csq_value is None:
                continue

            fixed = {
                "record_id": columns[2],
                "CHROM": columns[0],
                "POS": columns[1],
                "REF": columns[3],
                "ALT": columns[4],
            }
            rows = []
            for value in csq_value.split(","):
                parts = value.split("|")
                parts.extend([""] * (len(fields) - len(parts)))
                row = dict(zip(fields, parts))
                row.update(fixed)
                rows.append(row)
            yield fixed["record_id"], rows


def load_table(path, key):
    with path.open(newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def integer(value, default):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def picker_key(row, priority, ranks):
    transcript = bare_transcript(row.get("Feature", ""))
    meta = priority.get(transcript, {})
    consequence_rank = min(
        (ranks.get(term, 1000) for term in row.get("Consequence", "").split("&")),
        default=1000,
    )
    return (
        0 if integer(meta.get("mane_select"), 0) else 1,
        0 if integer(meta.get("mane_plus_clinical"), 0) else 1,
        0 if integer(meta.get("canonical"), 0) else 1,
        integer(meta.get("appris"), 100),
        integer(meta.get("tsl"), 100),
        0 if integer(meta.get("protein_coding"), 0) else 1,
        0 if integer(meta.get("ccds"), 0) else 1,
        consequence_rank,
        integer(meta.get("length"), 0),
        integer(meta.get("ensembl"), 1),
        integer(meta.get("refseq"), 1),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastvep", required=True, type=Path)
    parser.add_argument("--transcript-priority", required=True, type=Path)
    parser.add_argument("--consequence-ranks", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--vep-oracle", type=Path)
    args = parser.parse_args()

    priority = load_table(args.transcript_priority, "transcript_id")
    rank_rows = load_table(args.consequence_ranks, "consequence")
    ranks = {key: int(value["rank"]) for key, value in rank_rows.items()}

    fields = [
        "record_id", "CHROM", "POS", "REF", "ALT", "Feature", "Gene",
        "SYMBOL", "BIOTYPE", "Consequence", "EXON", "INTRON",
        "CDS_position", "Protein_position", "HGVS_OFFSET", "MANE_SELECT",
        "MANE_PLUS_CLINICAL", "CANONICAL", "TSL", "APPRIS", "CCDS",
    ]
    oracle = None
    if args.vep_oracle:
        oracle = {
            record_id: bare_transcript(rows[0].get("Feature", ""))
            for record_id, rows in iter_csq_records(args.vep_oracle)
        }

    record_count = 0
    differences = []
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record_id, rows in iter_csq_records(args.fastvep):
            row = min(rows, key=lambda candidate: picker_key(candidate, priority, ranks))
            writer.writerow({field: row.get(field, "") for field in fields})
            record_count += 1
            if oracle is not None:
                expected_transcript = oracle.get(record_id, "")
                actual_transcript = bare_transcript(row.get("Feature", ""))
                if actual_transcript != expected_transcript:
                    differences.append((record_id, expected_transcript, actual_transcript))

    print("records={:,}".format(record_count))
    if args.vep_oracle:
        print(
            "transcript_match={:,}/{:,} ({:.1%})".format(
                record_count - len(differences), record_count,
                (record_count - len(differences)) / record_count,
            )
        )
        counts = Counter((expected, actual) for _, expected, actual in differences)
        for (expected, actual), count in counts.most_common(20):
            print("difference\t{}\t{}\t{}".format(count, expected, actual))
        mismatch_path = args.output.with_suffix(".transcript-mismatches.tsv")
        with mismatch_path.open("w") as handle:
            handle.write("record_id\tvep_transcript\tfastvep_picked_transcript\n")
            for row in differences:
                handle.write("\t".join(row) + "\n")
        print("transcript_mismatches={}".format(mismatch_path))


if __name__ == "__main__":
    main()

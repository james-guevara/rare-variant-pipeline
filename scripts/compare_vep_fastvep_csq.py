#!/usr/bin/env python3
"""Compare VEP-picked CSQ rows with the same transcript emitted by fastVEP."""

import argparse
import gzip
import re
from collections import Counter
from pathlib import Path


CSQ_FORMAT = re.compile(r"Format: ([^\"]+)")


def open_text(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_csq(path):
    fields = []
    records = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = CSQ_FORMAT.search(line)
                if not match:
                    raise ValueError(f"Could not parse CSQ schema in {path}")
                fields = match.group(1).split("|")
            elif line.startswith("#"):
                continue
            else:
                columns = line.rstrip().split("\t")
                record_id = columns[2]
                csq_value = next(
                    value[4:]
                    for value in columns[7].split(";")
                    if value.startswith("CSQ=")
                )
                rows = []
                for value in csq_value.split(","):
                    parts = value.split("|")
                    parts.extend([""] * (len(fields) - len(parts)))
                    row = dict(zip(fields, parts))
                    row.update({
                        "record_id": record_id,
                        "CHROM": columns[0],
                        "POS": columns[1],
                        "REF": columns[3],
                        "ALT": columns[4],
                    })
                    rows.append(row)
                records[record_id] = rows
    return fields, records


def bare_transcript(value):
    return value.split(".", 1)[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep", required=True, type=Path)
    parser.add_argument("--fastvep", required=True, type=Path)
    args = parser.parse_args()

    _, vep = read_csq(args.vep)
    _, fast = read_csq(args.fastvep)
    total = len(vep)
    matched = 0
    unmatched = Counter()
    comparison_fields = (
        "Gene", "BIOTYPE", "STRAND", "Consequence", "EXON", "INTRON",
        "CDS_position", "Protein_position", "HGVS_OFFSET",
    )
    agreements = Counter()
    denominators = Counter()
    consequence_differences = Counter()
    hgvs_offset_differences = Counter()
    mismatch_rows = []

    for record_id, vep_rows in vep.items():
        if len(vep_rows) != 1:
            unmatched[f"vep_rows={len(vep_rows)}"] += 1
            continue
        selected = vep_rows[0]
        transcript = bare_transcript(selected.get("Feature", ""))
        candidates = [
            row for row in fast.get(record_id, [])
            if bare_transcript(row.get("Feature", "")) == transcript
        ]
        if not candidates:
            unmatched["selected_transcript_absent"] += 1
            continue
        matched += 1
        candidate = candidates[0]
        for field in comparison_fields:
            denominators[field] += 1
            if selected.get(field, "") == candidate.get(field, ""):
                agreements[field] += 1
        if selected.get("Consequence", "") != candidate.get("Consequence", ""):
            consequence_differences[(
                selected.get("Consequence", ""),
                candidate.get("Consequence", ""),
            )] += 1
            mismatch_rows.append((
                record_id,
                selected.get("Feature", ""),
                selected.get("Gene", ""),
                selected.get("Consequence", ""),
                candidate.get("Consequence", ""),
                selected.get("HGVS_OFFSET", ""),
                candidate.get("HGVS_OFFSET", ""),
            ))
        if selected.get("HGVS_OFFSET", "") != candidate.get("HGVS_OFFSET", ""):
            hgvs_offset_differences[(
                selected.get("HGVS_OFFSET", "") or "<blank>",
                candidate.get("HGVS_OFFSET", "") or "<blank>",
            )] += 1

    print(f"records={total:,}")
    print(f"selected_transcript_matched={matched:,}/{total:,} ({matched/total:.1%})")
    for reason, count in unmatched.most_common():
        print(f"unmatched_{reason}={count:,}")
    for field in comparison_fields:
        numerator = agreements[field]
        denominator = denominators[field]
        print(f"{field}={numerator:,}/{denominator:,} ({numerator/denominator:.1%})")
    print("consequence_difference_pairs:")
    for pair, count in consequence_differences.most_common(20):
        print(f"  {count:,}\t{pair[0]}\t{pair[1]}")
    print("hgvs_offset_difference_pairs:")
    for pair, count in hgvs_offset_differences.most_common(20):
        print(f"  {count:,}\t{pair[0]}\t{pair[1]}")
    mismatch_path = args.fastvep.with_suffix(".consequence-mismatches.tsv")
    with mismatch_path.open("w") as handle:
        handle.write(
            "record_id\tvep_feature\tgene\tvep_consequence\t"
            "fastvep_consequence\tvep_hgvs_offset\tfastvep_hgvs_offset\n"
        )
        for row in mismatch_rows:
            handle.write("\t".join(row) + "\n")
    print(f"consequence_mismatches={mismatch_path}")


if __name__ == "__main__":
    main()

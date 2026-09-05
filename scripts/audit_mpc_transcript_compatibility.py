#!/usr/bin/env python3
"""Audit whether dbNSFP MPC transcript contexts match selected missense contexts."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import duckdb


def values(value):
    return str(value or "").split(";")


def numeric(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", required=True, type=Path)
    parser.add_argument("--dbnsfp", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    con = duckdb.connect()
    rows = con.execute(
        """
        SELECT s.record_id, s.Feature AS selected_transcript,
               d.Ensembl_transcriptid, d.aaref, d.aaalt,
               d.MPC_score, d.MPC_rankscore
        FROM read_parquet(?) s
        JOIN read_parquet(?) d
          ON CAST(s.POS AS BIGINT) = CAST(d."pos(1-based)" AS BIGINT)
         AND s.REF = d.ref AND s.ALT = d.alt
        """,
        [str(args.selected), str(args.dbnsfp)],
    ).fetchall()

    audit = defaultdict(lambda: {
        "selected_transcript": "", "transcripts": set(), "aa_pairs": set(),
        "selected_aa_pairs": set(), "mpc_aa_pairs": set(),
        "mpc_values": [], "selected_mpc_values": [], "rankscore_values": [],
    })
    for record_id, selected, transcript_text, ref_text, alt_text, mpc_text, rank in rows:
        item = audit[record_id]
        item["selected_transcript"] = str(selected).split(".", 1)[0]
        transcripts = values(transcript_text)
        refs = values(ref_text)
        alts = values(alt_text)
        mpcs = values(mpc_text)
        width = max(len(transcripts), len(refs), len(alts), len(mpcs))
        for series in (transcripts, refs, alts, mpcs):
            series.extend([""] * (width - len(series)))
        for transcript, ref, alt, mpc in zip(
            transcripts, refs, alts, mpcs
        ):
            transcript = transcript.split(".", 1)[0]
            pair = (ref, alt) if ref not in ("", ".") and alt not in ("", ".") else None
            if transcript not in ("", "."):
                item["transcripts"].add(transcript)
            if pair:
                item["aa_pairs"].add(pair)
                if transcript == item["selected_transcript"]:
                    item["selected_aa_pairs"].add(pair)
            score = numeric(mpc)
            if score is not None:
                item["mpc_values"].append(score)
                if transcript == item["selected_transcript"]:
                    item["selected_mpc_values"].append(score)
                if pair:
                    item["mpc_aa_pairs"].add(pair)
        rank_value = numeric(rank)
        if rank_value is not None:
            item["rankscore_values"].append(rank_value)

    fields = [
        "record_id", "selected_transcript", "selected_transcript_in_dbnsfp",
        "n_dbnsfp_aa_pairs", "selected_aa_known", "mpc_raw_available",
        "selected_transcript_mpc_available", "mpc_rankscore_available",
        "mpc_context_compatible",
    ]
    counts = defaultdict(int)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record_id in sorted(audit):
            item = audit[record_id]
            selected_in = item["selected_transcript"] in item["transcripts"]
            selected_known = bool(item["selected_aa_pairs"])
            raw_available = bool(item["mpc_values"])
            selected_mpc_available = bool(item["selected_mpc_values"])
            rank_available = bool(item["rankscore_values"])
            compatible = (
                not raw_available
                or len(item["aa_pairs"]) <= 1
                or bool(item["selected_aa_pairs"] & item["mpc_aa_pairs"])
            )
            writer.writerow({
                "record_id": record_id,
                "selected_transcript": item["selected_transcript"],
                "selected_transcript_in_dbnsfp": int(selected_in),
                "n_dbnsfp_aa_pairs": len(item["aa_pairs"]),
                "selected_aa_known": int(selected_known),
                "mpc_raw_available": int(raw_available),
                "selected_transcript_mpc_available": int(selected_mpc_available),
                "mpc_rankscore_available": int(rank_available),
                "mpc_context_compatible": int(compatible),
            })
            counts["records"] += 1
            counts["selected_in"] += selected_in
            counts["one_aa_pair"] += len(item["aa_pairs"]) <= 1
            counts["raw"] += raw_available
            counts["selected_mpc"] += selected_mpc_available
            counts["rank"] += rank_available
            counts["compatible"] += compatible

    for key in (
        "records", "selected_in", "one_aa_pair", "raw", "selected_mpc",
        "rank", "compatible",
    ):
        print("{}={:,}".format(key, counts[key]))


if __name__ == "__main__":
    main()

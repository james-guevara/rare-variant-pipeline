#!/usr/bin/env python3
"""Build the integration sample manifest from participant IDs and a PLINK PED."""
import argparse, csv
from pathlib import Path

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--participants", required=True, type=Path)
    p.add_argument("--ped", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    with a.participants.open(newline="") as h:
        participants = list(csv.DictReader(h, delimiter="\t"))
    pedigree = {}
    with a.ped.open(newline="") as h:
        for row in csv.reader(h, delimiter="\t"):
            if len(row) < 6: raise ValueError("PED rows must have six columns")
            fid, iid, _, _, sex, _ = row[:6]
            if iid in pedigree: raise ValueError(f"duplicate PED IID: {iid}")
            pedigree[iid] = (fid, {"1": "M", "2": "F", "0": ""}.get(sex, ""))
    missing = [r["IID"] for r in participants if r["IID"] not in pedigree]
    if missing: raise ValueError(f"participants absent from PED: {len(missing)}")
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open("w", newline="") as h:
        writer = csv.writer(h, delimiter="\t", lineterminator="\n")
        writer.writerow(["FID", "IID", "SEX"])
        for row in participants:
            iid = row["IID"]; fid, sex = pedigree[iid]
            writer.writerow([fid, iid, sex])

if __name__ == "__main__": main()

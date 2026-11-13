import csv, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/build_cohort_sample_manifest.py"

def test_builds_manifest_in_participant_order(tmp_path: Path):
    participants = tmp_path / "participants.tsv"
    participants.write_text("IID\nS2\nS1\n")
    ped = tmp_path / "cohort.ped"
    ped.write_text("F1\tS1\t0\t0\t1\t-9\nF2\tS2\t0\t0\t2\t-9\n")
    output = tmp_path / "manifest.tsv"
    subprocess.run([sys.executable, str(SCRIPT), "--participants", str(participants), "--ped", str(ped), "--output", str(output)], check=True)
    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert rows == [{"FID": "F2", "IID": "S2", "SEX": "F"}, {"FID": "F1", "IID": "S1", "SEX": "M"}]

import csv
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "build_analysis_reports.py"


def test_reports_capture_missingness_modes_and_hashes(tmp_path):
    dataset = tmp_path / "analysis_dataset.tsv"
    dictionary = tmp_path / "analysis_dataset_dictionary.tsv"
    qc = tmp_path / "analysis_qc.json"
    missingness = tmp_path / "missingness_report.tsv"
    manifest = tmp_path / "run_manifest.json"

    dataset.write_text("FID\tIID\tPGS_trait\tlof_t1\nF1\tS1\t0.5\t1\nF2\tS2\t\t0\n")
    dictionary.write_text(
        "variable\tdata_type\nFID\tstring\nIID\tstring\nPGS_trait\tfloat\nlof_t1\tinteger\n"
    )
    qc.write_text(json.dumps({
        "cohort_id": "pilot",
        "components": {
            "pgs": {"missing_policy": "allow"},
            "rare": {"missing_policy": "error"},
            "cnv": {"missing_policy": "allow"},
        },
    }))

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--dataset", str(dataset),
            "--dictionary", str(dictionary),
            "--qc", str(qc),
            "--cohort-id", "pilot",
            "--pgs-mode", "precomputed",
            "--rare-mode", "computed",
            "--cnv-mode", "disabled",
            "--missingness", str(missingness),
            "--run-manifest", str(manifest),
        ],
        check=True,
    )

    with missingness.open() as handle:
        rows = {row["variable"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["PGS_trait"]["missing_participants"] == "1"
    assert rows["lof_t1"]["missing_participants"] == "0"

    report = json.loads(manifest.read_text())
    assert report["participants"] == 2
    assert report["variables"] == 4
    assert report["component_modes"] == {
        "pgs": "precomputed", "rare": "computed", "cnv": "disabled"
    }
    assert set(report["artifacts"]) == {
        "analysis_dataset.tsv",
        "analysis_dataset_dictionary.tsv",
        "analysis_qc.json",
    }

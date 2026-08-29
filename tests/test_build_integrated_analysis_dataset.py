import csv
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "build_integrated_analysis_dataset.py"
TEMPLATE = REPO / "resources" / "integrated-analysis-variables.tsv"


def rows(path):
    return list(csv.DictReader(path.open(), delimiter="\t"))


def test_integrates_on_pgs_universe_and_reports_rare_only(tmp_path):
    pgs = tmp_path / "pgs.tsv"
    pgs.write_text("IID\tANCESTRY\tPGS_trait\nS1\tEUR\t0.5\nS2\tAFR\t-0.2\n")
    dictionary = tmp_path / "pgs_dictionary.tsv"
    dictionary.write_text(
        "variable\tdata_type\tnullable\tdescription\tsource\n"
        "IID\tstring\tfalse\tID\tPGS\n"
        "ANCESTRY\tcategorical\tfalse\tAncestry\tPGS\n"
        "PGS_trait\tfloat\tfalse\tTrait score\tPGS\n"
    )
    rare = tmp_path / "rare.tsv"
    rare.write_text(
        "FID\tIID\tSEX\tlof_t1\tlof_t2\tmiss_t1\tmiss_t2\tmiss_t3\tmiss_t4\n"
        "F1\tS1\tF\t1\t0\t0\t1\t0\t2\n"
        "F2\tS2\tM\t0\t1\t0\t0\t1\t0\n"
        "F3\tS3\tF\t0\t0\t0\t0\t0\t0\n"
    )
    output, merged_dictionary, qc = (
        tmp_path / "integrated.tsv",
        tmp_path / "integrated_dictionary.tsv",
        tmp_path / "qc.json",
    )
    subprocess.run([
        sys.executable, str(SCRIPT),
        "--pgs-dataset", str(pgs),
        "--pgs-dictionary", str(dictionary),
        "--rare-burdens", str(rare),
        "--variable-template", str(TEMPLATE),
        "--cohort-id", "test",
        "--output", str(output),
        "--dictionary", str(merged_dictionary),
        "--qc", str(qc),
    ], check=True)

    integrated = rows(output)
    assert [row["IID"] for row in integrated] == ["S1", "S2"]
    assert integrated[0]["FID"] == "F1" and integrated[0]["lof_t1"] == "1"
    assert integrated[1]["SEX"] == "M" and integrated[1]["miss_t3"] == "1"
    assert [row["variable"] for row in rows(merged_dictionary)] == list(integrated[0])
    report = json.loads(qc.read_text())
    assert report["integrated_participants"] == 2
    assert report["rare_only_participants_excluded"] == 1


def test_rejects_pgs_participant_without_completed_rare_burdens(tmp_path):
    pgs = tmp_path / "pgs.tsv"; pgs.write_text("IID\nS1\n")
    dictionary = tmp_path / "dictionary.tsv"
    dictionary.write_text("variable\tdata_type\tnullable\tdescription\tsource\nIID\tstring\tfalse\tID\tPGS\n")
    rare = tmp_path / "rare.tsv"
    rare.write_text("FID\tIID\tSEX\tlof_t1\tlof_t2\tmiss_t1\tmiss_t2\tmiss_t3\tmiss_t4\n")
    result = subprocess.run([
        sys.executable, str(SCRIPT), "--pgs-dataset", str(pgs),
        "--pgs-dictionary", str(dictionary), "--rare-burdens", str(rare),
        "--variable-template", str(TEMPLATE), "--cohort-id", "test",
        "--output", str(tmp_path / "out"), "--dictionary", str(tmp_path / "dict"),
        "--qc", str(tmp_path / "qc"),
    ], capture_output=True, text=True)
    assert result.returncode != 0
    assert "PGS participants are absent" in result.stderr

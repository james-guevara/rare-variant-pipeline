import csv
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "build_analysis_dataset.py"
TEMPLATE = REPO / "resources" / "integrated-analysis-variables.tsv"


def read_rows(path):
    return list(csv.DictReader(path.open(), delimiter="\t"))


def run(tmp_path, extra):
    outputs = {
        "output": tmp_path / "analysis.tsv",
        "dictionary": tmp_path / "dictionary.tsv",
        "qc": tmp_path / "qc.json",
        "exclusions": tmp_path / "exclusions.tsv",
    }
    command = [
        sys.executable, str(SCRIPT), "--participant-manifest", str(tmp_path / "samples.tsv"),
        "--variable-template", str(TEMPLATE), "--cohort-id", "test",
    ]
    for name, path in outputs.items():
        command.extend([f"--{name}", str(path)])
    subprocess.run([*command, *extra], check=True)
    return outputs


def test_manifest_universe_combines_optional_components(tmp_path):
    (tmp_path / "samples.tsv").write_text(
        "#FID\tIID\tSEX\nF1\tS1\t2\nF1\tS2\t1\nF2\tS3\t2\n"
    )
    (tmp_path / "pgs.tsv").write_text("IID\tANCESTRY\tPGS_trait\nS1\tEUR\t0.5\nS2\tAFR\t-0.2\n")
    (tmp_path / "pgs.dict.tsv").write_text(
        "variable\tdata_type\tnullable\tdescription\tsource\n"
        "ANCESTRY\tcategorical\ttrue\tAncestry\tPGS\n"
        "PGS_trait\tfloat\ttrue\tTrait score\tPGS\n"
    )
    (tmp_path / "rare.tsv").write_text(
        "FID\tIID\tSEX\tlof_t1\tlof_t2\tmiss_t1\tmiss_t2\tmiss_t3\tmiss_t4\n"
        "F1\tS1\tF\t1\t0\t0\t1\t0\t2\nF1\tS2\tM\t0\t1\t0\t0\t1\t0\n"
        "F2\tS3\tF\t0\t0\t0\t0\t0\t0\n"
    )
    (tmp_path / "cnv.tsv").write_text(
        "IID\tCNV_DEL_GENE_COUNT\tCNV_DUP_GENE_COUNT\nS1\t2\t1\nS3\t0\t3\n"
    )
    (tmp_path / "cnv.dict.tsv").write_text(
        "variable\trole\tdata_type\tnullable\tdefault_analysis_use\tsource\tdescription\n"
        "CNV_DEL_GENE_COUNT\tpredictor\tinteger\ttrue\tCNV burden\tCNV\tDeleted genes\n"
        "CNV_DUP_GENE_COUNT\tpredictor\tinteger\ttrue\tCNV burden\tCNV\tDuplicated genes\n"
    )
    result = run(tmp_path, [
        "--pgs-dataset", str(tmp_path / "pgs.tsv"),
        "--pgs-dictionary", str(tmp_path / "pgs.dict.tsv"),
        "--rare-burdens", str(tmp_path / "rare.tsv"),
        "--cnv-dataset", str(tmp_path / "cnv.tsv"),
        "--cnv-dictionary", str(tmp_path / "cnv.dict.tsv"),
    ])
    rows = read_rows(result["output"])
    assert [row["IID"] for row in rows] == ["S1", "S2", "S3"]
    assert [row["SEX"] for row in rows] == ["F", "M", "F"]
    assert rows[1]["CNV_DEL_GENE_COUNT"] == ""
    assert rows[2]["PGS_trait"] == ""
    assert [row["variable"] for row in read_rows(result["dictionary"])] == list(rows[0])
    qc = json.loads(result["qc"].read_text())
    assert qc["analysis_universe"] == "participant_manifest"
    assert qc["components"]["pgs"]["manifest_participants_missing"] == 1
    assert qc["components"]["cnv"]["manifest_participants_missing"] == 1


def test_exclude_policy_is_component_specific(tmp_path):
    (tmp_path / "samples.tsv").write_text("FID\tIID\tSEX\nF1\tS1\tF\nF2\tS2\tM\n")
    (tmp_path / "pgs.tsv").write_text("IID\tPGS_trait\nS1\t1.0\n")
    (tmp_path / "pgs.dict.tsv").write_text(
        "variable\tdata_type\tnullable\tdescription\tsource\n"
        "PGS_trait\tfloat\ttrue\tTrait\tPGS\n"
    )
    result = run(tmp_path, [
        "--pgs-dataset", str(tmp_path / "pgs.tsv"),
        "--pgs-dictionary", str(tmp_path / "pgs.dict.tsv"),
        "--missing-pgs-policy", "exclude",
    ])
    assert [row["IID"] for row in read_rows(result["output"])] == ["S1"]
    assert read_rows(result["exclusions"]) == [{"IID": "S2", "reason": "missing_from_pgs"}]

import csv
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "initialize_cohort.py"


def test_initializes_from_joint_vcf_and_psam_without_claiming_derived_data(tmp_path):
    psam = tmp_path / "cohort.psam"
    psam.write_text("#FID IID SEX\nF1 S1 2\n0 S2 1\n")
    output = tmp_path / "initialized"
    subprocess.run([
        sys.executable, str(SCRIPT), "--cohort", "demo",
        "--joint-vcf", "s3://bucket/cohort.vcf.gz", "--psam", str(psam),
        "--chromosomes", "chr1,X,Y", "--shared-resources-root", "/shared/rvp",
        "--cohort-root", "/cohorts/demo", "--output-dir", str(output),
    ], check=True)

    samples = list(csv.DictReader((output / "sample_manifest.tsv").open(), delimiter="\t"))
    assert samples == [
        {"FID": "F1", "IID": "S1", "SEX": "F"},
        {"FID": "", "IID": "S2", "SEX": "M"},
    ]
    plan = list(csv.DictReader((output / "chromosome_preparation.tsv").open(), delimiter="\t"))
    assert [row["chromosome"] for row in plan] == ["chr1", "chrX", "chrY"]
    assert all(row["preparation_state"] == "PENDING_DERIVED_RESOURCES" for row in plan)
    assert plan[0]["input_vcf"] == "s3://bucket/cohort.vcf.gz"
    assert plan[1]["zarr_store"] == "/cohorts/demo/zarr/chrX.sharded-v3.zarr"
    assert plan[1]["vcz_work_root"] == "/cohorts/demo/vcz-work/chrX"
    qc = json.loads((output / "initialization_qc.json").read_text())
    assert qc["derived_resources_ready"] is False
    assert qc["sex_counts"] == {"F": 1, "M": 1}


def test_per_chromosome_template_must_have_placeholder(tmp_path):
    psam = tmp_path / "cohort.psam"; psam.write_text("#IID\nS1\n")
    result = subprocess.run([
        sys.executable, str(SCRIPT), "--cohort", "demo", "--vcf-template", "input.vcf.gz",
        "--psam", str(psam), "--shared-resources-root", "/shared",
        "--cohort-root", "/cohort", "--output-dir", str(tmp_path / "out"),
    ], capture_output=True, text=True)
    assert result.returncode != 0
    assert "must contain {chrom}" in result.stderr

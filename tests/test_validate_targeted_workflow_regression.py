import csv
import json
import subprocess
import sys
from pathlib import Path

import duckdb


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/validate_targeted_workflow_regression.py"


def _write_run(root: Path, *, symbol: str = "GENE1") -> None:
    root.mkdir()
    with (root / "06.plof-tiered.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["CHROM", "POS", "REF", "ALT", "SYMBOL", "lof_tier"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {"CHROM": "22", "POS": 100, "REF": "A", "ALT": "G", "SYMBOL": symbol, "lof_tier": "lof_t1"}
        )
    con = duckdb.connect()
    con.execute(
        "COPY (SELECT 1::BIGINT AS variant_index) TO ? (FORMAT PARQUET)",
        [str(root / "01.target-alleles.parquet")],
    )
    con.execute(
        "COPY (SELECT 'S1' AS sample_id, 100::BIGINT AS POS) TO ? (FORMAT PARQUET)",
        [str(root / "07.plof.carriers.parquet")],
    )
    con.execute(
        "COPY (SELECT 'v1' AS record_id, 1::BIGINT AS carrier_count) TO ? (FORMAT PARQUET)",
        [str(root / "07.plof.genotype-summary.parquet")],
    )


def test_exact_workflow_regression_passes(tmp_path: Path):
    reference, new = tmp_path / "reference", tmp_path / "new"
    _write_run(reference)
    _write_run(new)
    expected = tmp_path / "expected.json"
    expected.write_text(json.dumps({"expected_counts": {
        "targeted_alleles": 1, "qualifying_variants": 1, "lof_t1": 1,
        "lof_t2": 0, "carrier_rows": 1, "carrier_samples": 1,
    }}))
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--new", str(new), "--reference", str(reference),
         "--expectations", str(expected)],
        check=True, capture_output=True, text=True,
    )
    assert "regression=PASS" in result.stdout


def test_annotation_drift_fails(tmp_path: Path):
    reference, new = tmp_path / "reference", tmp_path / "new"
    _write_run(reference)
    _write_run(new, symbol="CHANGED")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--new", str(new), "--reference", str(reference)],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "rows differ" in result.stderr

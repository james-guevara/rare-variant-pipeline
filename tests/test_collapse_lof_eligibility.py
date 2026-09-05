import csv
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


REPO = Path(__file__).resolve().parents[1]


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_primary_and_sensitivity_lof_burdens_are_separate(tmp_path: Path):
    source = tmp_path / "carriers.parquet"
    pq.write_table(pa.Table.from_pylist([
        {"sample_id": "primary", "Gene": "G1", "SYMBOL": "ONE", "lof_tier": "lof_t1", "record_id": "v1", "primary_analysis_eligible": True, "burden_count_available": True},
        {"sample_id": "primary", "Gene": "G1", "SYMBOL": "ONE", "lof_tier": "lof_t1", "record_id": "v2", "primary_analysis_eligible": True, "burden_count_available": True},
        {"sample_id": "ambiguous", "Gene": "G2", "SYMBOL": "TWO", "lof_tier": "lof_t2", "record_id": "v3", "primary_analysis_eligible": False, "burden_count_available": True},
        {"sample_id": "qc_only", "Gene": "G3", "SYMBOL": "THREE", "lof_tier": "lof_t2", "record_id": "v4", "primary_analysis_eligible": False, "burden_count_available": False},
        {"sample_id": "primary", "Gene": "G4", "SYMBOL": "FOUR", "lof_tier": "", "record_id": "v5", "primary_analysis_eligible": True, "burden_count_available": True},
    ]), source)
    outputs = {name: tmp_path / f"{name}.tsv" for name in (
        "primary_gene", "primary_burden", "sensitivity_gene", "sensitivity_burden"
    )}
    subprocess.run([
        sys.executable, str(REPO / "scripts/collapse_lof_carriers.py"),
        "--input", str(source),
        "--sample-gene-output", str(outputs["primary_gene"]),
        "--sample-burden-output", str(outputs["primary_burden"]),
        "--eligibility-col", "primary_analysis_eligible",
        "--burden-available-col", "burden_count_available",
        "--sensitivity-sample-gene-output", str(outputs["sensitivity_gene"]),
        "--sensitivity-sample-burden-output", str(outputs["sensitivity_burden"]),
    ], check=True)
    primary = read_tsv(outputs["primary_burden"])
    sensitivity = read_tsv(outputs["sensitivity_burden"])
    assert primary == [{
        "SAMPLE": "primary", "lof_t1_genes": "1", "lof_t1_variants": "2",
        "lof_t2_genes": "0", "lof_t2_variants": "0", "any_tier_genes": "1",
        "any_tier_variants": "2",
    }]
    assert [row["SAMPLE"] for row in sensitivity] == ["ambiguous"]

import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


REPO = Path(__file__).resolve().parents[1]


def test_primary_counts_exclude_but_input_retains_sensitivity_rows(tmp_path: Path):
    source = tmp_path / "carriers.parquet"
    pq.write_table(pa.Table.from_pylist([
        {"sample_id": "primary", "miss_tier": "miss_t1", "primary_analysis_eligible": True},
        {"sample_id": "ambiguous", "miss_tier": "miss_t1", "primary_analysis_eligible": False},
    ]), source)
    counts = tmp_path / "counts.tsv"
    totals = tmp_path / "totals.tsv"
    subprocess.run([
        sys.executable, str(REPO / "scripts/postprocess/count_carriers.py"),
        "--input", str(source), "--sample-col", "sample_id",
        "--group-col", "miss_tier", "--eligibility-col", "primary_analysis_eligible",
        "--out-counts", str(counts), "--out-totals", str(totals),
    ], check=True)
    assert "primary" in counts.read_text()
    assert "ambiguous" not in counts.read_text()
    assert pq.read_table(source).num_rows == 2

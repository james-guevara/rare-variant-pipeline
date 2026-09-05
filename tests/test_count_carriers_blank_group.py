import csv
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).parents[1]


def test_blank_tier_is_not_a_count_group(tmp_path):
    source = tmp_path / "carriers.parquet"
    pq.write_table(pa.Table.from_pylist([
        {"sample_id": "S1", "lof_tier": "lof_t1", "record_id": "v1"},
        {"sample_id": "S1", "lof_tier": "", "record_id": "v2"},
    ]), source)
    counts = tmp_path / "counts.tsv"
    totals = tmp_path / "totals.tsv"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/postprocess/count_carriers.py"),
        "--input", str(source), "--sample-col", "sample_id", "--group-col", "lof_tier",
        "--out-counts", str(counts), "--out-totals", str(totals),
    ], check=True)
    with counts.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row == {"SAMPLE": "S1", "lof_t1": "1", "any_group": "1"}
    with totals.open() as handle:
        total_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert total_rows == [{"lof_tier": "lof_t1", "n_carrier_rows": "1", "n_samples_with_any": "1"}]

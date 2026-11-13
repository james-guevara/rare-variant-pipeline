import json
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).parents[1]


def test_consolidates_and_validates_chromosome_partitions(tmp_path):
    chr2 = tmp_path / "chr2.parquet"
    chrx = tmp_path / "chrX.parquet"
    pq.write_table(pa.Table.from_pylist([
        {"#CHROM": "chr2", "POS": 20, "REF": "A", "ALT": "G", "sample_id": "S2", "record_id": "v2"},
        {"#CHROM": "chr2", "POS": 10, "REF": "C", "ALT": "T", "sample_id": "S1", "record_id": "v1"},
    ]), chr2)
    pq.write_table(pa.Table.from_pylist([
        {"#CHROM": "chrX", "POS": 1, "REF": "G", "ALT": "A", "sample_id": "S1", "record_id": "vx", "primary_analysis_eligible": True},
    ]), chrx)
    output = tmp_path / "genome.parquet"
    validation = tmp_path / "validation.json"
    subprocess.run([
        sys.executable, str(ROOT / "scripts" / "consolidate_carrier_parquets.py"),
        "--input", str(chrx), str(chr2), "--output", str(output),
        "--validation-output", str(validation),
    ], check=True)
    table = pq.read_table(output).to_pylist()
    assert [(row["#CHROM"], row["POS"]) for row in table] == [
        ("chr2", 10), ("chr2", 20), ("chrX", 1)
    ]
    report = json.loads(validation.read_text())
    assert report["status"] == "PASS"
    assert report["input_files"] == 2
    assert report["rows"] == 3

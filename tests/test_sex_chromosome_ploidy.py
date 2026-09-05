import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


REPO = Path(__file__).resolve().parents[1]


def load_audit_module():
    path = REPO / "scripts/audit_sex_chromosome_ploidy.py"
    spec = importlib.util.spec_from_file_location("audit_sex_chromosome_ploidy", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_inferred_karyotype_thresholds():
    classify = load_audit_module().classify
    assert classify(0.95, 0.01, 0.99) == "XY-like"
    assert classify(0.01, 0.95, 0.02) == "XX-like"
    assert classify(0.75, 0.20, 0.60) == "ambiguous"


def test_qc_uses_mask_derived_state_not_gt_spelling(tmp_path: Path):
    rows = []
    for sample_id, gt, ploidy, dosage, ad in (
        ("s1", "1", 1, 1, "0,20"),
        ("s2", "1/.", 1, 1, "0,20"),
        ("s3", "./1", 1, 1, "0,20"),
        ("s4", "0/1", 2, 1, "10,10"),
        ("s5", "1/.", 1, 1, "10,10"),
        # Deliberately misleading text: numeric mask-derived state wins.
        ("s6", "0/1", 1, 1, "0,20"),
    ):
        rows.append({
            "sample_id": sample_id,
            "GT": gt,
            "called_ploidy": ploidy,
            "alt_dosage": dosage,
            "AD": ad,
            "GQ": 40,
            "DP": 20,
        })
    source = tmp_path / "source.parquet"
    output = tmp_path / "output.parquet"
    pq.write_table(pa.Table.from_pylist(rows), source)
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "cohorts": {"g2mh": {"input_dir": "."}},
        "output_base": ".",
        "qc": {
            "min_gq": 20,
            "min_dp": 10,
            "het_ab_min": 0.25,
            "het_ab_max": 0.75,
            "hom_ab_min": 0.9,
        },
    }))
    subprocess.run([
        sys.executable,
        str(REPO / "scripts/postprocess/qc_genotype.py"),
        "--cohort", "g2mh",
        "--chrom", "chrX",
        "--resources", str(config),
        "--input", str(source),
        "--output", str(output),
    ], check=True)
    kept = set(pq.read_table(output, columns=["sample_id"])["sample_id"].to_pylist())
    assert kept == {"s1", "s2", "s3", "s4", "s6"}

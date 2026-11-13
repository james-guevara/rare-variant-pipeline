import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import zarr


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "validate_vcz_plan.py"


def test_validates_store_and_advances_plan(tmp_path):
    store_path = tmp_path / "chr1.sharded-v3.zarr"
    store = zarr.open_group(store_path, mode="w", zarr_format=3)
    store.create_array("sample_id", data=np.asarray(["S1", "S2"], dtype="U2"))
    store.create_array("variant_position", data=np.asarray([10, 20], dtype="i4"))
    store.create_array("variant_allele", data=np.asarray([["A", "C"], ["G", "T"]], dtype="U1"))
    store.create_array("call_genotype", data=np.zeros((2, 2, 2), dtype="i1"))
    store.create_array("call_GQ", data=np.ones((2, 2), dtype="i2"))
    store.create_array("call_DP", data=np.ones((2, 2), dtype="i2"))
    store.create_array("call_AD", data=np.ones((2, 2, 2), dtype="i2"))
    Path(f"{store_path}.complete").touch()
    plan = tmp_path / "plan.tsv"
    plan.write_text(
        "chromosome\tzarr_store\tpreparation_state\n"
        f"chr1\t{store_path}\tREADY_FOR_ZARR\n"
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text("FID\tIID\tSEX\n\tS1\t\n\tS2\t\n")
    output, report = tmp_path / "validated.tsv", tmp_path / "report.json"
    subprocess.run([
        sys.executable, str(SCRIPT), "--preparation-plan", str(plan),
        "--sample-manifest", str(samples), "--output-plan", str(output),
        "--report", str(report),
    ], check=True)
    row = next(csv.DictReader(output.open(), delimiter="\t"))
    assert row["preparation_state"] == "READY_FOR_DERIVED_RESOURCES"
    assert json.loads(report.read_text())["chromosomes"][0]["variants"] == 2


def test_accepts_localized_allele_depth_arrays(tmp_path):
    store_path = tmp_path / "chrY.sharded-v3.zarr"
    store = zarr.open_group(store_path, mode="w", zarr_format=3)
    store.create_array("sample_id", data=np.asarray(["S1"], dtype="U2"))
    store.create_array("variant_position", data=np.asarray([10], dtype="i4"))
    store.create_array("variant_allele", data=np.asarray([["A", "C"]], dtype="U1"))
    store.create_array("call_genotype", data=np.zeros((1, 1, 2), dtype="i1"))
    store.create_array("call_GQ", data=np.ones((1, 1), dtype="i2"))
    store.create_array("call_DP", data=np.ones((1, 1), dtype="i2"))
    store.create_array("call_LAD", data=np.ones((1, 1, 2), dtype="i2"))
    store.create_array("call_LAA", data=np.ones((1, 1, 1), dtype="i1"))
    Path(f"{store_path}.complete").touch()
    plan = tmp_path / "plan.tsv"
    plan.write_text(
        "chromosome\tzarr_store\tpreparation_state\n"
        f"chrY\t{store_path}\tREADY_FOR_ZARR\n"
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text("FID\tIID\tSEX\n\tS1\t\n")
    output, report = tmp_path / "validated.tsv", tmp_path / "report.json"
    subprocess.run([
        sys.executable, str(SCRIPT), "--preparation-plan", str(plan),
        "--sample-manifest", str(samples), "--output-plan", str(output),
        "--report", str(report),
    ], check=True)
    chromosome = json.loads(report.read_text())["chromosomes"][0]
    assert chromosome["allele_depth_encoding"] == "LAD+LAA"

import hashlib
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "package_targeted_outputs.py"


def science_manifest(tmp_path: Path) -> Path:
    path = tmp_path / "science.json"
    path.write_text(json.dumps({"cohort": "test", "chromosome": "chr22"}))
    return path


def test_packages_available_counts_and_records_hashes(tmp_path: Path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    (run_root / "_SUCCESS").touch()
    lof = run_root / "12.plof-per-sample-counts.tsv"
    lof.write_text("SAMPLE\tlof_t1\tlof_t2\nS1\t1\t0\n")
    miss = run_root / "12.missense-per-sample-counts.tsv"
    miss.write_text("SAMPLE\tmiss_t1\tmiss_t4\nS1\t0\t2\n")
    bindings = tmp_path / "bindings.json"
    bindings.write_text(json.dumps({"run_root": str(run_root)}))
    output = tmp_path / "package"

    subprocess.run(
        [sys.executable, str(SCRIPT), "--science-manifest", str(science_manifest(tmp_path)),
         "--bindings", str(bindings),
         "--output-dir", str(output)],
        check=True,
    )

    manifest = json.loads((output / "targeted-output-manifest.json").read_text())
    assert manifest["status"] == "SUCCEEDED"
    assert (manifest["cohort"], manifest["chromosome"]) == ("test", "chr22")
    assert set(manifest["files"]) == {"plof_counts", "missense_counts"}
    assert manifest["files"]["plof_counts"]["sha256"] == hashlib.sha256(
        lof.read_bytes()
    ).hexdigest()
    assert (output / lof.name).read_bytes() == lof.read_bytes()


def test_rejects_completed_run_without_required_lof_counts(tmp_path: Path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    (run_root / "_SUCCESS").touch()
    bindings = tmp_path / "bindings.json"
    bindings.write_text(json.dumps({"run_root": str(run_root)}))

    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--science-manifest", str(science_manifest(tmp_path)),
         "--bindings", str(bindings),
         "--output-dir", str(tmp_path / "package")],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "missing required LoF counts" in result.stderr


def test_preflight_package_contains_no_scientific_outputs(tmp_path: Path):
    bindings = tmp_path / "bindings.json"
    bindings.write_text(json.dumps({"run_root": str(tmp_path / "absent")}))
    output = tmp_path / "package"

    subprocess.run(
        [sys.executable, str(SCRIPT), "--science-manifest", str(science_manifest(tmp_path)),
         "--bindings", str(bindings),
         "--output-dir", str(output), "--preflight-only"],
        check=True,
    )
    manifest = json.loads((output / "targeted-output-manifest.json").read_text())
    assert manifest == {
        "schema_version": 1, "cohort": "test", "chromosome": "chr22",
        "status": "PREFLIGHT_ONLY", "files": {},
    }

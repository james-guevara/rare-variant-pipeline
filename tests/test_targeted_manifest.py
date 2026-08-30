import hashlib
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/run_targeted_manifest.py"


def _fixture(tmp_path: Path) -> tuple[Path, Path]:
    resources = {}
    specs = {}
    for name in ("target_bed", "genebayes"):
        path = tmp_path / name
        path.write_text(name)
        resources[name] = str(path)
        specs[name] = {"kind": "file"}
    for name, sentinel in (
        ("zarr_store", "zarr.json"),
        ("annotation_root", "cache"),
        ("loftee_root", "transcripts.sqlite"),
    ):
        path = tmp_path / name
        path.mkdir()
        (path / sentinel).write_text(name)
        resources[name] = str(path)
        specs[name] = {"kind": "directory", "sentinel": sentinel}
    digest = hashlib.sha256(Path(resources["target_bed"]).read_bytes()).hexdigest()
    specs["target_bed"]["sha256"] = digest
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({
        "schema_version": 1, "cohort": "test", "chromosome": "chr22",
        "contig": "22", "contig_length": 50818468,
        "thresholds": {"population_af_max": 0.01, "cohort_af_max": 0.02},
        "resources": specs,
    }))
    bindings = tmp_path / "bindings.json"
    bindings.write_text(json.dumps({
        "schema_version": 1, "run_root": str(tmp_path / "run"), "resources": resources,
    }))
    return manifest, bindings


def test_manifest_preflight_validates_bindings_and_checksums(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--bindings", str(bindings), "--preflight-only", "--skip-runtime-checks"],
        check=True, capture_output=True, text=True,
    )
    assert "preflight=PASS" in result.stdout
    assert "checksum.target_bed=" in result.stdout


def test_manifest_preflight_fails_on_checksum_drift(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    document = json.loads(bindings.read_text())
    Path(document["resources"]["target_bed"]).write_text("changed")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--bindings", str(bindings), "--preflight-only", "--skip-runtime-checks"],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "resource checksum differs" in result.stderr


def test_manifest_preflight_rejects_site_qc_drift(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    postprocess = tmp_path / "postprocess.json"
    postprocess.write_text(json.dumps({"qc": {"min_gq": 10}}))
    science = json.loads(manifest.read_text())
    science["qc"] = {"min_gq": 20}
    science["resources"]["postprocess_config"] = {"kind": "file"}
    manifest.write_text(json.dumps(science))
    deployment = json.loads(bindings.read_text())
    deployment["resources"]["postprocess_config"] = str(postprocess)
    bindings.write_text(json.dumps(deployment))
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--bindings", str(bindings), "--preflight-only", "--skip-runtime-checks"],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "postprocess QC differs" in result.stderr


def test_sex_chromosome_manifest_requires_qc_and_par_resources(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    science = json.loads(manifest.read_text())
    science["chromosome"] = "chrX"
    science["contig"] = "X"
    science["contig_length"] = 156040895
    manifest.write_text(json.dumps(science))
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--bindings", str(bindings), "--preflight-only", "--skip-runtime-checks"],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "sex-chromosome manifest requires resources" in result.stderr


def test_lof_only_does_not_require_declared_missense_binding(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    science = json.loads(manifest.read_text())
    science["resources"]["missense_candidates"] = {"kind": "file"}
    manifest.write_text(json.dumps(science))
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--bindings", str(bindings), "--lof-only", "--preflight-only",
         "--skip-runtime-checks"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "preflight=PASS" in result.stdout

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import duckdb

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/run_targeted_manifest.py"
SPEC = importlib.util.spec_from_file_location("run_targeted_manifest", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
resolve = MODULE.resolve


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


def test_manifest_preflight_accepts_all_observed_mode_without_target_bed(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    science = json.loads(manifest.read_text())
    science["resources"].pop("target_bed")
    manifest.write_text(json.dumps(science))
    deployment = json.loads(bindings.read_text())
    deployment["resources"].pop("target_bed")
    bindings.write_text(json.dumps(deployment))
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--bindings", str(bindings), "--preflight-only", "--skip-runtime-checks"],
        check=True, capture_output=True, text=True,
    )
    assert "preflight=PASS" in result.stdout


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
    science["resources"]["missense_dbnsfp"] = {"kind": "file"}
    manifest.write_text(json.dumps(science))
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--bindings", str(bindings), "--lof-only", "--preflight-only",
         "--skip-runtime-checks"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "preflight=PASS" in result.stdout


def test_synonymous_controls_default_on_and_can_be_disabled(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    science = json.loads(manifest.read_text())
    deployment = json.loads(bindings.read_text())
    environment, _ = resolve(science, deployment)
    assert environment["SYNONYMOUS_TIERED_CONTROLS"] == "1"
    science["optional_outputs"] = {"synonymous_tiered_controls": False}
    environment, _ = resolve(science, deployment)
    assert environment["SYNONYMOUS_TIERED_CONTROLS"] == "0"


def test_tier_configuration_is_validated_and_exported(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    science = json.loads(manifest.read_text())
    deployment = json.loads(bindings.read_text())
    science["tiers"] = {
        "lof": {
            "lof_t1_min_genebayes_post_mean": 0.2,
            "lof_t2_min_genebayes_post_mean": 0.04,
        },
        "missense": {
            "score_thresholds": {
                "ClinPred_rankscore": 0.5,
                "AlphaMissense_rankscore": 0.9,
                "popEVE_converted_rankscore": 0.8,
                "MPC_rankscore": 0.7,
            },
            "tier_pass_counts": {
                "miss_t1": 4, "miss_t2": 3, "miss_t3": 2, "miss_t4": 1,
            },
        },
    }
    environment, _ = resolve(science, deployment)
    assert environment["LOF_T1_MIN_GENEBAYES_POST_MEAN"] == "0.2"
    assert environment["LOF_T2_MIN_GENEBAYES_POST_MEAN"] == "0.04"
    assert environment["MISSENSE_MPC_RANKSCORE_MIN"] == "0.7"
    assert environment["MISS_T1_PASS_COUNT"] == "4"


def test_preflight_rejects_missense_file_without_score_columns(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    missense = tmp_path / "missense.parquet"
    duckdb.connect().execute("COPY (SELECT 1 AS position) TO ? (FORMAT PARQUET)", [str(missense)])
    science = json.loads(manifest.read_text())
    deployment = json.loads(bindings.read_text())
    science["resources"]["missense_dbnsfp"] = {"kind": "file"}
    deployment["resources"]["missense_dbnsfp"] = str(missense)
    try:
        resolve(science, deployment)
        assert False, "missing score columns should fail"
    except ValueError as error:
        assert "lacks required score columns" in str(error)


def test_family_genotypes_are_optional_and_require_sample_manifest(tmp_path: Path):
    manifest, bindings = _fixture(tmp_path)
    science = json.loads(manifest.read_text())
    deployment = json.loads(bindings.read_text())
    science["optional_outputs"] = {"family_genotypes": True}
    try:
        resolve(science, deployment)
        assert False, "missing sample manifest should fail"
    except ValueError as error:
        assert "sample_manifest" in str(error)

    sample_manifest = tmp_path / "samples.tsv"
    sample_manifest.write_text("FID\tIID\tSEX\nF1\tS1\tF\n")
    science["resources"]["sample_manifest"] = {"kind": "file"}
    deployment["resources"]["sample_manifest"] = str(sample_manifest)
    environment, _ = resolve(science, deployment)
    assert environment["FAMILY_GENOTYPES"] == "1"
    assert environment["SAMPLE_MANIFEST"] == str(sample_manifest)

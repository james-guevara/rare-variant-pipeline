import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def test_targeted_wrapper_compiles_and_exposes_receipt(tmp_path: Path):
    manifest = tmp_path / "manifest.json"
    bindings = tmp_path / "bindings.json"
    manifest.write_text("{}")
    bindings.write_text("{}")

    result = subprocess.run(
        [
            "nextflow", "run", "targeted.nf", "-preview",
            "--manifest", str(manifest),
            "--bindings", str(bindings),
            "--targeted_container", "example.invalid/targeted:test",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    workflow = (REPO / "workflows" / "targeted_manifest.nf").read_text()
    assert "workflow TARGETED_MANIFEST_WORKFLOW" in workflow
    assert "execution_receipt = TARGETED_CHROMOSOME.out[0]" in workflow
    assert "chromosome_outputs = TARGETED_CHROMOSOME.out[1]" in workflow

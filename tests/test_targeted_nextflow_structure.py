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
    assert "take:\n    manifest_bindings" in workflow
    assert "TARGETED_CHROMOSOME(manifest_bindings)" in workflow
    assert "execution_receipt = TARGETED_CHROMOSOME.out[0]" in workflow
    assert "chromosome_outputs = TARGETED_CHROMOSOME.out[1]" in workflow


def test_cohort_wrapper_compiles_with_scatter_and_gather(tmp_path: Path):
    manifest = tmp_path / "manifest.json"
    bindings = tmp_path / "bindings.json"
    samples = tmp_path / "samples.tsv"
    sheet = tmp_path / "runs.tsv"
    manifest.write_text("{}")
    bindings.write_text("{}")
    samples.write_text("FID\tIID\tSEX\nF1\tS1\tXX\n")
    sheet.write_text(
        "chromosome\tmanifest\tbindings\n"
        f"chrY\t{manifest}\t{bindings}\n"
    )

    result = subprocess.run(
        [
            "nextflow", "run", "cohort.nf", "-preview",
            "--run_sheet", str(sheet),
            "--sample_manifest", str(samples),
            "--expected_chromosomes", "chrY",
            "--targeted_container", "example.invalid/targeted:test",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    source = (REPO / "workflows" / "cohort_rare_burden.nf").read_text()
    assert "workflow COHORT_RARE_BURDEN_WORKFLOW" in source
    assert "TARGETED_MANIFEST_WORKFLOW(manifest_bindings)" in source
    assert "RARE_BURDEN_GATHER_WORKFLOW(" in source
    assert "rare_burdens = RARE_BURDEN_GATHER_WORKFLOW.out.rare_burdens" in source

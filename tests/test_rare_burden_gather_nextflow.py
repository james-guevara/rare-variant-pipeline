import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

def test_gather_wrapper_compiles(tmp_path: Path):
    samples = tmp_path / "samples.tsv"
    samples.write_text("IID\nS1\n")
    package = tmp_path / "chr1"
    package.mkdir()
    result = subprocess.run([
        "nextflow", "run", "gather.nf", "-preview",
        "--sample_manifest", str(samples), "--package_glob", str(package),
        "--expected_chromosomes", "chr1",
        "--targeted_container", "example.invalid/targeted:test",
    ], cwd=REPO, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    workflow = (REPO / "workflows/rare_burden_gather.nf").read_text()
    assert "workflow RARE_BURDEN_GATHER_WORKFLOW" in workflow
    assert "rare_burdens = GATHER_RARE_BURDENS.out.burdens" in workflow
    assert "gene_set_burdens = GATHER_RARE_BURDENS.out.gene_set_burdens" in workflow
    assert "plof_carriers = GATHER_RARE_BURDENS.out.plof_carriers" in workflow
    assert "--package-root packages" in workflow

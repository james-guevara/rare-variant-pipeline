import json
import importlib.util
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "submit_vcz_plan.py"


def load_module():
    spec = importlib.util.spec_from_file_location("submit_vcz_plan", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_dry_run_builds_one_job_per_ready_plan_row(tmp_path):
    plan = tmp_path / "plan.tsv"
    plan.write_text(
        "chromosome\tcontig\tinput_vcf\tzarr_store\tvcz_work_root\tpreparation_state\n"
        "chr1\t1\t/fsx/input/cohort.vcf.gz\t/fsx/cohort/zarr/chr1.sharded-v3.zarr\t/fsx/cohort/vcz-work/chr1\tREADY_FOR_ZARR\n"
        "chrX\tX\t/fsx/input/cohort.vcf.gz\t/fsx/cohort/zarr/chrX.sharded-v3.zarr\t/fsx/cohort/vcz-work/chrX\tREADY_FOR_ZARR\n"
    )
    output = tmp_path / "submission.json"
    subprocess.run([
        sys.executable, str(SCRIPT), "--preparation-plan", str(plan),
        "--cohort", "demo", "--queue", "queue", "--job-definition", "definition:1",
        "--output", str(output),
    ], check=True)
    document = json.loads(output.read_text())
    assert document["mode"] == "DRY_RUN"
    assert [job["chromosome"] for job in document["jobs"]] == ["chr1", "chrX"]
    assert all(job["status"] == "PLANNED" for job in document["jobs"])


def test_rejects_row_not_ready_for_zarr(tmp_path):
    plan = tmp_path / "plan.tsv"
    plan.write_text(
        "chromosome\tcontig\tinput_vcf\tzarr_store\tvcz_work_root\tpreparation_state\n"
        "chr1\tauto\tinput.vcf.gz\t/out/chr1.sharded-v3.zarr\t/work/chr1\tPENDING_DERIVED_RESOURCES\n"
    )
    result = subprocess.run([
        sys.executable, str(SCRIPT), "--preparation-plan", str(plan),
        "--cohort", "demo", "--queue", "queue", "--job-definition", "definition:1",
        "--output", str(tmp_path / "out.json"),
    ], capture_output=True, text=True)
    assert result.returncode != 0
    assert "expected READY_FOR_ZARR" in result.stderr


def test_adapter_resolves_prefixed_contigs_and_rejects_empty_subsets():
    adapter = load_module().ADAPTER
    assert 'for candidate in "$CONTIG" "$CHROMOSOME" "$token"' in adapter
    assert 'test "$record_count" -gt 0' in adapter


def test_adapter_invalidates_conversion_when_staged_source_changes():
    adapter = load_module().ADAPTER
    assert 'source_fingerprint=$(sha256sum "$adapted"' in adapter
    assert 'previous_fingerprint' in adapter
    assert 'rm -rf "$conversion_work/chr${token}" "$ZARR_STORE"' in adapter

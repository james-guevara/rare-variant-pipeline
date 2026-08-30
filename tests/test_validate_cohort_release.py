import hashlib
import json
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "validate_cohort_release.py"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_validates_pinned_gathered_outputs(tmp_path):
    burdens = tmp_path / "rare.tsv"
    strata = tmp_path / "strata.tsv"
    burdens.write_text(
        "FID\tIID\tSEX\tlof_t1\tlof_t2\tmiss_t1\tmiss_t2\tmiss_t3\tmiss_t4\n"
        "F1\tS1\tF\t1\t2\t3\t4\t5\t6\n"
    )
    strata.write_text(
        "FID\tIID\tSEX\tburden_partition\tincluded_in_primary_total\tchromosomes_complete\tlof_t1\tlof_t2\tmiss_t1\tmiss_t2\tmiss_t3\tmiss_t4\n"
        "F1\tS1\tF\tautosomal\ttrue\ttrue\t1\t2\t3\t4\t5\t6\n"
        "F1\tS1\tF\tsex_chromosome_primary\ttrue\ttrue\t0\t0\t0\t0\t0\t0\n"
        "F1\tS1\tF\tsex_chromosome_sensitivity\tfalse\ttrue\t0\t0\t\t\t\t\n"
    )
    contract = tmp_path / "contract.json"
    contract.write_text(json.dumps({
        "release_id": "fixture",
        "expected_outputs": {
            "sample_count": 1,
            "strata_rows": 3,
            "primary_tier_sums": {
                "lof_t1": 1, "lof_t2": 2, "miss_t1": 3,
                "miss_t2": 4, "miss_t3": 5, "miss_t4": 6,
            },
            "rare_burdens_sha256": sha256(burdens),
            "chromosome_strata_sha256": sha256(strata),
        },
    }))
    result = subprocess.run(
        ["python", SCRIPT, "--contract", contract, "--burdens", burdens,
         "--strata", strata], capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "release=PASS id=fixture" in result.stdout


def test_fails_if_gathered_output_changes(tmp_path):
    burdens = tmp_path / "rare.tsv"
    strata = tmp_path / "strata.tsv"
    burdens.write_text("changed\n")
    strata.write_text("changed\n")
    contract = tmp_path / "contract.json"
    contract.write_text(json.dumps({
        "release_id": "fixture",
        "expected_outputs": {
            "rare_burdens_sha256": "0" * 64,
            "chromosome_strata_sha256": "0" * 64,
        },
    }))
    result = subprocess.run(
        ["python", SCRIPT, "--contract", contract, "--burdens", burdens,
         "--strata", strata], capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "rare_burdens_sha256 differs" in result.stderr

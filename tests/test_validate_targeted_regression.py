import importlib.util
from pathlib import Path


path = Path(__file__).resolve().parents[1] / "scripts/validate_targeted_regression.py"
spec = importlib.util.spec_from_file_location("validate_targeted_regression", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
add_legacy_missense_aliases = module.add_legacy_missense_aliases


def test_existing_missense_expectation_names_are_populated():
    observed = {
        "missense_burden_eligible_rows": 10,
        "missense_burden_eligible_alleles": 8,
        "missense_burden_eligible_samples": 7,
        "missense_primary_burden_rows": 9,
        "missense_sensitivity_only_rows": 1,
        "missense_primary_miss_t1_rows": 1,
        "missense_primary_miss_t2_rows": 2,
        "missense_primary_miss_t3_rows": 3,
        "missense_primary_miss_t4_rows": 3,
    }
    add_legacy_missense_aliases(observed)
    assert observed["sensitivity_only_burden_rows"] == 1
    assert observed["primary_miss_t4_rows"] == 3
    assert observed["burden_eligible_rows"] == 10

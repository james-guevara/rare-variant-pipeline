import importlib.util
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def load_module():
    path = REPO / "scripts/annotate_sex_chromosome_carriers.py"
    spec = importlib.util.spec_from_file_location("annotate_sex_carriers", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


REGIONS = {
    "pseudoautosomal_regions": {
        "chrX": [{"start": 10_001, "end": 2_781_479}],
        "chrY": [{"start": 10_001, "end": 2_781_479}],
    },
    "par_canonical_representation": "chrX",
}


def test_ambiguous_nonpar_counts_are_retained_but_not_primary():
    result = load_module().policy("chrX", 50_000_000, "ambiguous", REGIONS)
    assert result["burden_count_available"] is True
    assert result["primary_analysis_eligible"] is False
    assert result["frequency_denominator_eligible"] is False
    assert result["sensitivity_analysis_group"] == "ambiguous_karyotype"


def test_par_is_autosomal_but_chrY_duplicate_is_not_counted():
    policy = load_module().policy
    x = policy("chrX", 20_000, "ambiguous", REGIONS)
    y = policy("chrY", 20_000, "XY-like", REGIONS)
    assert x["burden_count_available"] is True
    assert x["primary_analysis_eligible"] is True
    assert y["par_duplicate_excluded"] is True
    assert y["burden_count_available"] is False
    assert y["primary_analysis_eligible"] is False


def test_y_nonpar_primary_requires_xy_like_pattern():
    policy = load_module().policy
    assert policy("chrY", 5_000_000, "XY-like", REGIONS)["primary_analysis_eligible"]
    assert not policy("chrY", 5_000_000, "XX-like", REGIONS)["primary_analysis_eligible"]
    xx = policy("chrY", 5_000_000, "XX-like", REGIONS)
    assert xx["burden_count_available"] is False
    assert xx["sensitivity_analysis_group"] == "qc_only_y_ineligible_karyotype"

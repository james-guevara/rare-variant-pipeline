import importlib.util
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def load_module():
    path = REPO / "scripts/interpret_sex_chromosome_qc.py"
    spec = importlib.util.spec_from_file_location("interpret_sex_chromosome_qc", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


THRESHOLDS = {
    "one_x_dp_ratio_min": 0.3,
    "one_x_dp_ratio_max": 0.7,
    "two_x_dp_ratio_min": 0.75,
    "two_x_dp_ratio_max": 1.25,
    "y_present_call_rate_min": 0.5,
    "y_absent_call_rate_max": 0.2,
    "excess_y_dp_ratio_min": 0.75,
}


def row(x, y, calls):
    return {
        "x_autosome_dp_ratio": str(x),
        "y_autosome_dp_ratio": str(y),
        "y_call_rate": str(calls),
    }


def test_non_diagnostic_evidence_patterns():
    interpret = load_module().interpret
    assert interpret(row(0.94, 0.46, 0.98), THRESHOLDS) == "two-X-plus-Y-compatible"
    assert interpret(row(0.47, 0.91, 0.98), THRESHOLDS) == "one-X-with-excess-Y-signal"
    assert interpret(row(0.52, 0.36, 0.96), THRESHOLDS) == "one-X-plus-Y-with-GT-ploidy-discordance"
    assert interpret(row(0.95, 0.05, 0.28), THRESHOLDS) == "two-X-with-uncertain-Y-signal"

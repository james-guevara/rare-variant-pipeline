import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_standalone_loftee.py"
SPEC = importlib.util.spec_from_file_location("run_standalone_loftee", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_recognizes_rows_without_transcript_annotations():
    assert MODULE.transcript_is_absent("")
    assert MODULE.transcript_is_absent("-")
    assert MODULE.transcript_is_absent(".")
    assert not MODULE.transcript_is_absent("ENST00000335137.4")

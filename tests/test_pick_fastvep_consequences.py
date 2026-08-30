import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "pick_fastvep_consequences.py"
SPEC = importlib.util.spec_from_file_location("pick_fastvep_consequences", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_picker_key_is_orderable_and_prioritizes_mane():
    rows = [
        {"Feature": "ENST_OTHER", "Consequence": "stop_gained"},
        {"Feature": "ENST_MANE", "Consequence": "missense_variant"},
    ]
    priority = {
        "ENST_OTHER": {},
        "ENST_MANE": {"mane_select": "1"},
    }
    ranks = {"stop_gained": 1, "missense_variant": 10}

    assert min(rows, key=lambda row: MODULE.picker_key(row, priority, ranks)) == rows[1]


def test_alternate_lof_rows_requires_picked_missense_and_other_transcript():
    selected = {
        "Feature": "ENST_SELECTED.1",
        "BIOTYPE": "protein_coding",
        "Consequence": "missense_variant",
    }
    alternate = {
        "Feature": "ENST_ALTERNATE.2",
        "BIOTYPE": "protein_coding",
        "Consequence": "stop_gained",
    }
    same_transcript = {
        "Feature": "ENST_SELECTED.2",
        "BIOTYPE": "protein_coding",
        "Consequence": "splice_donor_variant",
    }
    noncoding = {
        "Feature": "ENST_NONCODING",
        "BIOTYPE": "processed_transcript",
        "Consequence": "frameshift_variant",
    }

    assert MODULE.alternate_lof_rows(
        [selected, alternate, same_transcript, noncoding], selected
    ) == [alternate]


def test_alternate_lof_rows_ignores_nonmissense_picked_row():
    selected = {
        "Feature": "ENST_SELECTED",
        "BIOTYPE": "protein_coding",
        "Consequence": "synonymous_variant",
    }
    alternate = {
        "Feature": "ENST_ALTERNATE",
        "BIOTYPE": "protein_coding",
        "Consequence": "stop_gained",
    }

    assert MODULE.alternate_lof_rows([selected, alternate], selected) == []

import importlib.util
from pathlib import Path

import numpy as np


REPO = Path(__file__).parents[1]
SCRIPT = REPO / "scripts" / "extract_zarr_allele_genotypes.py"
SPEC = importlib.util.spec_from_file_location("extract_zarr_genotypes", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_render_genotype_preserves_source_multiallelic_codes():
    genotype = np.array([0, 2])
    called = np.array([True, True])

    assert MODULE.render_genotype(genotype, called, "/") == "0/2"
    assert MODULE.render_genotype(genotype, called, "/", alt_index=2) == "0/1"


def test_render_genotype_recodes_selected_alt_and_preserves_phase_and_missingness():
    genotype = np.array([2, 2, -1])
    called = np.array([True, True, False])

    assert MODULE.render_genotype(genotype, called, "|", alt_index=2) == "1|1|."

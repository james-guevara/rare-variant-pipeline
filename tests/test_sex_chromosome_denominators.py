import importlib.util
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[1]


def load_module():
    path = REPO / "scripts/extract_zarr_allele_genotypes.py"
    spec = importlib.util.spec_from_file_location("extract_zarr_genotypes", path)
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
KARYOTYPES = np.asarray(["XX-like", "XY-like", "ambiguous"])


def test_nonpar_denominators_exclude_ambiguous_and_y_requires_xy():
    mask = load_module().primary_sample_mask
    assert mask("chrX", 50_000_000, KARYOTYPES, REGIONS).tolist() == [True, True, False]
    assert mask("chrY", 5_000_000, KARYOTYPES, REGIONS).tolist() == [False, True, False]


def test_par_denominator_uses_only_canonical_x_representation():
    mask = load_module().primary_sample_mask
    assert mask("chrX", 20_000, KARYOTYPES, REGIONS).tolist() == [True, True, True]
    assert mask("chrY", 20_000, KARYOTYPES, REGIONS).tolist() == [False, False, False]

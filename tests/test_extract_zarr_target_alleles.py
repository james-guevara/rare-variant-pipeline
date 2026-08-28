import importlib.util
from pathlib import Path

import numpy as np


path = Path(__file__).resolve().parents[1] / "scripts/extract_zarr_target_alleles.py"
spec = importlib.util.spec_from_file_location("extract_zarr_target_alleles", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def brute_force(positions, lengths, intervals):
    selected = []
    for index, (position, length) in enumerate(zip(positions, lengths, strict=True)):
        start = int(position) - 1
        end = start + int(length)
        if any(start < bed_end and end > bed_start for bed_start, bed_end in intervals):
            selected.append(index)
    return np.asarray(selected, dtype=np.int64)


def test_vectorized_overlap_matches_bed_semantics():
    positions = np.asarray([1, 5, 10, 14, 20, 30], dtype=np.int32)
    lengths = np.asarray([1, 6, 1, 8, 1, 2], dtype=np.int16)
    intervals = [(0, 1), (9, 10), (15, 20), (29, 31)]
    expected = brute_force(positions, lengths, intervals)
    observed = module.overlapping_indexes(positions, lengths, intervals)
    np.testing.assert_array_equal(observed, expected)


def test_vectorized_overlap_handles_empty_inputs():
    positions = np.asarray([], dtype=np.int32)
    lengths = np.asarray([], dtype=np.int16)
    np.testing.assert_array_equal(
        module.overlapping_indexes(positions, lengths, [(0, 10)]),
        np.asarray([], dtype=np.int64),
    )
    np.testing.assert_array_equal(
        module.overlapping_indexes(
            np.asarray([1], dtype=np.int32),
            np.asarray([1], dtype=np.int16),
            [],
        ),
        np.asarray([], dtype=np.int64),
    )

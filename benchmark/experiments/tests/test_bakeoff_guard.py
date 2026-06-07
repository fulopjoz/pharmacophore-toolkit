"""Guard test: load_bench refuses a dataset that lacks a valid reference SDF, so a
non-CCR2 dataset can never silently use CCR2 references for the ref-based scorers."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import bakeoff


def test_load_bench_raises_on_missing_ref():
    caller = lambda m: (["c1ccccc1"], ["CCCC"], {})
    with pytest.raises(NotImplementedError):
        bakeoff.load_bench("created/load.py", caller, None)
    with pytest.raises(NotImplementedError):
        bakeoff.load_bench("created/load.py", caller, "/no/such/reference.sdf")


def test_dataset_specs_all_carry_a_ref():
    for spec in bakeoff.DATASET_SPECS:
        assert len(spec) == 4, "DATASET_SPECS must be (name, relpath, caller, ref_sdf)"
        assert spec[3] and os.path.exists(spec[3]), f"{spec[0]} has no valid ref SDF"

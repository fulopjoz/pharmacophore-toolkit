"""Smoke test for the prism scorer: registers + runs end-to-end, aligned length."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from harness import BenchData, discover, REGISTRY  # noqa: E402


def test_prism_registers_and_scores_aligned_length():
    discover()
    assert "prism" in REGISTRY
    act = ["c1ccc(CCN)cc1O", "c1ccc(CCN)cc1", "c1ccc2ccccc2c1CN"]
    dec = ["CCCCCCN", "CCCCCCO", "CCCCCCCC"]
    ref = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                       "tutorials", "data", "CCR2_reference_ligands.sdf"))
    data = BenchData.from_lists(act, dec, ref)
    train_idx = np.array([0, 1, 3, 4])
    test_idx = np.array([2, 5])
    scores = REGISTRY["prism"](data, train_idx, test_idx)
    assert scores.shape == (2,)
    assert np.all(np.isfinite(scores))

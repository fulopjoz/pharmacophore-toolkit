import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from harness import BenchData, discover, REGISTRY


def test_fixed_registers_runs_and_differs_from_learned():
    discover()
    assert "prism_fixed" in REGISTRY
    act = ["c1ccc(CCN)cc1O", "c1ccc(CCN)cc1", "c1ccc2ccccc2c1CN"]
    dec = ["CCCCCCN", "CCCCCCO", "CCCCCCCC"]
    ref = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                       "tutorials", "data", "CCR2_reference_ligands.sdf"))
    data = BenchData.from_lists(act, dec, ref)
    tr, te = np.array([0, 1, 3, 4]), np.array([2, 5])
    s_fixed = REGISTRY["prism_fixed"](data, tr, te)
    assert s_fixed.shape == (2,) and np.all(np.isfinite(s_fixed))
    s_learned = REGISTRY["prism"](data, tr, te)
    assert not np.allclose(s_fixed, s_learned)   # fixed mean != learned logistic

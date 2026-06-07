"""Data-readiness test: the committed CCR5 / CXCR4 created datasets load valid
actives + decoys. (Used to power the bake-off's cross-target Friedman test once
per-target references are wired — CCR2 references are not valid for CCR5/CXCR4.)"""
import os
import sys
import importlib.util

import pytest
from rdkit import Chem

_HERE = os.path.dirname(os.path.abspath(__file__))
_CREATED = os.path.dirname(_HERE)


def _load(target):
    spec = importlib.util.spec_from_file_location("created_load", os.path.join(_CREATED, "load.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.load(target)


@pytest.mark.parametrize("target", ["CCR5", "CXCR4"])
def test_created_target_loads_valid(target):
    if not os.path.isdir(os.path.join(_CREATED, target)):
        pytest.skip(f"{target} not present")
    act, dec, meta = _load(target)
    assert len(act) > 0 and len(dec) > 0
    for s in act[:30] + dec[:30]:
        assert Chem.MolFromSmiles(str(s)) is not None, f"unparseable SMILES: {s}"

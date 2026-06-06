"""Toolkit core: 2D pharmacophore-fingerprint scoring (Pharm2DScorer).

Wraps the toolkit's own `pharmacophore.pharm2d_scoring.Pharm2DScorer` — Gobbi
Pharm2D feature-pair fingerprints, Tanimoto vs the reference ligands, max over
refs. This is a faithful adapter (no reimplementation): construct = fit, the
references ARE the model.

UNSUPERVISED — `train_idx` is ignored (the query is the fixed reference set, not
learned from the train split), so it is leakage-free under the held-out split.
2D only: no conformer generation, fast.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from rdkit import Chem

from harness import register, BenchData

# Make the installed toolkit package importable from this experiments dir.
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
from pharmacophore.pharm2d_scoring import Pharm2DScorer  # noqa: E402


def _load_refs(data: BenchData):
    if "_pharm2d_refs" not in data._cache:
        data._cache["_pharm2d_refs"] = [
            m for m in Chem.SDMolSupplier(data.ref_sdf, removeHs=False) if m is not None
        ]
    return data._cache["_pharm2d_refs"]


@register("pharm2d")
def pharm2d(data: BenchData, train_idx, test_idx) -> np.ndarray:
    """Gobbi Pharm2D fingerprint Tanimoto vs reference ligands (toolkit core)."""
    refs = _load_refs(data)
    if not refs:
        return np.zeros(len(test_idx), dtype=float)
    scorer = Pharm2DScorer(refs)
    mols = data.mols()
    test_mols = [mols[i] for i in test_idx]
    scores = scorer.score_all(test_mols)
    return np.asarray(scores, dtype=float)

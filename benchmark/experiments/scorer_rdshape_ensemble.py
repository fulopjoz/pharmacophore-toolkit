"""Toolkit core: 3D shape+color ensemble scoring (ReferenceEnsembleScorer).

Wraps `pharmacophore.rdshape_optimizer.ReferenceEnsembleScorer` — rdShapeAlign
combo Tanimoto (shape + color, opt_param=0.5) of each query (ETKDGv3 conformers)
against the reference ensemble, aggregated as max over refs. This is the
toolkit's own 3D scoring core; a faithful adapter (constructor = fit, stateless
`score_batch` = score, CV-safe per the code-review map).

UNSUPERVISED — `train_idx` ignored (query = fixed reference set).

COST: generates conformers per query → compute-heavy on large sets. Tunable via
env vars so the front node can smoke-test cheaply and PBS can run at full
fidelity:
    BAKEOFF_NCONF   conformers per query   (default 5; ReferenceEnsembleScorer min is 1)
    BAKEOFF_NJOBS   parallel workers       (default 1)
"""
from __future__ import annotations

import os
import sys

import numpy as np
from rdkit import Chem

from harness import register, BenchData, SEED

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
from pharmacophore.rdshape_optimizer import ReferenceEnsembleScorer  # noqa: E402

_NCONF = int(os.environ.get("BAKEOFF_NCONF", "5"))
_NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))


def _load_refs(data: BenchData):
    if "_rdshape_refs" not in data._cache:
        data._cache["_rdshape_refs"] = [
            m for m in Chem.SDMolSupplier(data.ref_sdf, removeHs=False) if m is not None
        ]
    return data._cache["_rdshape_refs"]


@register("rdshape_ensemble")
def rdshape_ensemble(data: BenchData, train_idx, test_idx) -> np.ndarray:
    """rdShapeAlign combo Tanimoto vs reference ensemble (toolkit 3D core)."""
    refs = _load_refs(data)
    if not refs:
        return np.zeros(len(test_idx), dtype=float)
    scorer = ReferenceEnsembleScorer(
        refs, opt_param=0.5, n_conformers=_NCONF, aggregation="max",
        use_colors=True, random_seed=SEED, n_jobs=_NJOBS,
    )
    mols = data.mols()
    test_mols = [mols[i] for i in test_idx]
    scores = scorer.score_batch(test_mols)
    return np.asarray(scores, dtype=float)

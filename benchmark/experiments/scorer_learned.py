"""Toolkit core: supervised learned shape/color weighting (LearnedScorer).

Wraps `pharmacophore.learned_scoring.LearnedScorer` — logistic regression over
the per-reference [shape, color] feature vector produced by the toolkit's
`UnifiedEvaluator` (3D reference-mode alignment). This is the toolkit's flagship
SUPERVISED method: it learns the shape-vs-color balance from labeled actives /
decoys instead of a fixed opt_param.

Held-out adaptation (the key to a leakage-free comparison): the toolkit's
fit()/predict() both score *every* molecule held by their evaluator, so we
  1. build a TRAIN-only UnifiedEvaluator and fit the logistic on it, then
  2. repoint the fitted scorer at a TEST-only UnifiedEvaluator and predict.
The trained scaler+coefficients (from train) are applied to fresh test features.
Test rows are grouped by their true label only to keep both evaluator lists
non-empty; the per-molecule score never uses the label, so nothing leaks.

COST: 3D reference-mode alignment with conformers, fit + predict → the slowest
method here. PBS only for large datasets. Tunable:
    BAKEOFF_NCONF   conformers per molecule (floored at 5 per EvaluationConfig)
    BAKEOFF_NJOBS   parallel workers
"""
from __future__ import annotations

import os
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from harness import register, BenchData, SEED

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
from pharmacophore.evaluation import UnifiedEvaluator, EvaluationConfig  # noqa: E402
from pharmacophore.learned_scoring import LearnedScorer  # noqa: E402

_NCONF = max(5, int(os.environ.get("BAKEOFF_NCONF", "5")))   # EvaluationConfig requires >= 5
_NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))


def _load_refs(data: BenchData):
    if "_learned_refs" not in data._cache:
        data._cache["_learned_refs"] = [
            m for m in Chem.SDMolSupplier(data.ref_sdf, removeHs=False) if m is not None
        ]
    return data._cache["_learned_refs"]


def _embed(mol, nconf, seed):
    """Embed `nconf` conformers (ETKDGv3, with Hs) or return None on failure.

    Pre-embedding here is essential for robust held-out alignment: UnifiedEvaluator
    ._prepare_molecules SILENTLY DROPS any molecule whose embedding yields 0
    conformers (no identity preserved), which would desync per-molecule scores
    from their test indices. A molecule arriving with >= nconf conformers is
    REUSED by the evaluator (no re-embed, no drop), so passing pre-embedded mols
    guarantees the returned breakdowns are 1:1 with what we pass.
    """
    if mol is None:
        return None
    m = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.numThreads = 1                  # match evaluator; avoid RDKit threadpool issues
    AllChem.EmbedMultipleConfs(m, numConfs=nconf, params=p)
    return m if m.GetNumConformers() > 0 else None


@register("learned_scorer")
def learned_scorer(data: BenchData, train_idx, test_idx) -> np.ndarray:
    """Logistic on per-ref shape/color features; fit on train, score held-out test.

    Un-embeddable test molecules are floored to 0.0 (ranked last) rather than
    dropped, so the returned scores stay aligned with `test_idx`."""
    refs = _load_refs(data)
    if not refs:
        return np.zeros(len(test_idx), dtype=float)

    mols = data.mols()
    y = data.y
    cfg = EvaluationConfig(scoring_mode="reference", n_conformers=_NCONF, opt_param=0.5)

    # --- fit on TRAIN only (embeddable mols, so X rows align with y_true) ---
    tr = np.asarray(train_idx)
    train_act = [em for i in tr if y[i] == 1 for em in (_embed(mols[i], _NCONF, SEED),) if em is not None]
    train_dec = [em for i in tr if y[i] == 0 for em in (_embed(mols[i], _NCONF, SEED),) if em is not None]
    train_eval = UnifiedEvaluator(refs, train_act, train_dec, random_state=SEED, n_jobs=_NJOBS)
    ls = LearnedScorer(train_eval, C=1.0, random_state=SEED)
    ls.fit(cfg)

    # --- score held-out TEST: pre-embed, keep only embeddable (1:1), floor the rest ---
    te = np.asarray(test_idx)
    ok_pos, ok_mols = [], []
    for k in range(len(te)):
        em = _embed(mols[te[k]], _NCONF, SEED)
        if em is not None:
            ok_pos.append(k)
            ok_mols.append(em)

    out = np.zeros(len(te), dtype=float)              # floor: un-embeddable -> 0.0 (ranked last)
    if not ok_mols:
        return out

    # Pass all embeddable test mols as "actives" (decoys=[]); predict() only needs
    # the per-molecule breakdowns, not the labels, so this avoids any empty-class issue.
    test_eval = UnifiedEvaluator(refs, ok_mols, [], random_state=SEED, n_jobs=_NJOBS)
    ls.evaluator = test_eval                          # reuse train-fitted scaler + logistic
    preds = np.asarray(ls.predict(cfg), dtype=float)  # 1:1 with ok_mols (conformers reused)

    n = min(len(preds), len(ok_pos))                  # defensive: never broadcast-mismatch
    out[np.asarray(ok_pos[:n])] = preds[:n]
    return out

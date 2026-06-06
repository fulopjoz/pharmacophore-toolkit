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


@register("learned_scorer")
def learned_scorer(data: BenchData, train_idx, test_idx) -> np.ndarray:
    """Logistic on per-ref shape/color features; fit on train, score held-out test."""
    refs = _load_refs(data)
    if not refs:
        return np.zeros(len(test_idx), dtype=float)

    mols = data.mols()
    y = data.y
    cfg = EvaluationConfig(scoring_mode="reference", n_conformers=_NCONF, opt_param=0.5)

    # --- fit on TRAIN only ---
    tr = np.asarray(train_idx)
    train_act = [mols[i] for i in tr if y[i] == 1]
    train_dec = [mols[i] for i in tr if y[i] == 0]
    train_eval = UnifiedEvaluator(refs, train_act, train_dec, random_state=SEED, n_jobs=_NJOBS)
    ls = LearnedScorer(train_eval, C=1.0, random_state=SEED)
    ls.fit(cfg)

    # --- score held-out TEST: repoint the fitted scorer at a test-only evaluator ---
    te = np.asarray(test_idx)
    labels_te = y[te]
    act_pos = [k for k in range(len(te)) if labels_te[k] == 1]
    dec_pos = [k for k in range(len(te)) if labels_te[k] == 0]
    test_act = [mols[te[k]] for k in act_pos]
    test_dec = [mols[te[k]] for k in dec_pos]
    test_eval = UnifiedEvaluator(refs, test_act, test_dec, random_state=SEED, n_jobs=_NJOBS)

    ls.evaluator = test_eval          # reuse the train-fitted scaler + logistic model
    preds = np.asarray(ls.predict(cfg), dtype=float)   # order: [test_act..., test_dec...]

    out = np.zeros(len(te), dtype=float)
    out[act_pos] = preds[: len(act_pos)]
    out[dec_pos] = preds[len(act_pos):]
    return out

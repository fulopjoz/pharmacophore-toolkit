"""prism_esp — PRISM plus an electrostatic-similarity channel.

Same templates + per-(template×type) color features as prism, plus a per-template Carbó
ESP scalar (RDKit Gasteiger charges on the same color-optimized pose) -> (K*7) -> logistic.
Tests whether electrostatics adds discrimination beyond shape+color. Honest scope: an
ESP-similarity FEATURE on a color-optimized pose, not field-based alignment (ShaEP/EON).
Trust only if it lifts the unbiased (MUBD) set over prism with non-overlapping CIs;
the K extra columns are otherwise overfitting capacity."""
from __future__ import annotations

import os
import sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData, SEED
from prism_core import prism_features


@register("prism_esp")
def prism_esp(data: BenchData, train_idx, test_idx) -> np.ndarray:
    X, tr = prism_features(data, train_idx, with_esp=True)
    if X is None:
        return np.zeros(len(test_idx), dtype=float)
    sc = StandardScaler().fit(X[tr])
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
    lr.fit(sc.transform(X[tr]), data.y[tr])
    return lr.predict_proba(sc.transform(X[np.asarray(test_idx)]))[:, 1]

"""PRISM — Pharmacophore-Resolved Importance-weighted Similarity Model (was `s3_3d`).

Discrimination-weighted 3D color, feature-type × template. Per held-out split:
  1. cluster the TRAIN actives -> K templates (leakage-safe), embed;
  2. for every molecule, align to each template and decompose the color into 6
     per-feature-type overlaps -> (K*6) feature matrix (prism_core.prism_features);
  3. logistic regression (StandardScaler + LogisticRegression) fit on train rows
     learns the per-(type × template) discrimination weights and predicts test.

The logistic IS the discrimination weighting; the 3D-overlap features make it 3D.
The prism metaphor: a prism resolves white light into a spectrum, as PRISM resolves
ROCS "color" into its feature types and reweights them by learned importance.

Lineage (scite, 2026-06-07): weight-matrix learning over pharmacophore fingerprints
(Rehioui 2022, 10.1002/minf.202200210), PLS-DA on pharmacophore FPs (Askjaer 2008,
10.1021/ci700356w), PharmRF ML scoring (Kumar 2022, 10.1002/jcc.26840), multi-template
VS (Madzhidov 2020, 10.3390/molecules25020385)."""
from __future__ import annotations

import os
import sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData, SEED  # noqa: E402
from prism_core import prism_features  # noqa: E402


@register("prism")
def prism(data: BenchData, train_idx, test_idx) -> np.ndarray:
    X, tr = prism_features(data, train_idx, with_esp=False)
    if X is None:
        return np.zeros(len(test_idx), dtype=float)
    sc = StandardScaler().fit(X[tr])
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
    lr.fit(sc.transform(X[tr]), data.y[tr])
    return lr.predict_proba(sc.transform(X[np.asarray(test_idx)]))[:, 1]

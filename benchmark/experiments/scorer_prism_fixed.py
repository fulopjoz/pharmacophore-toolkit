"""prism_fixed — fixed-weight ablation of PRISM (isolates the learned weighting).

IDENTICAL per-(template×type) 3D-color features as prism, aggregated with EQUAL fixed
weights (mean) instead of a learned logistic. The mean is taken AFTER the same
StandardScaler(fit on train) that prism feeds its logistic, and ONLY over columns with
non-zero train variance — so the ONLY difference between prism and prism_fixed is
learned-vs-fixed weights on the *same standardized, informative features*.

Why the two guards matter: (1) a raw mean of unstandardized Tanimoto columns would
conflate weighting with scaling (hydrophobe/rings overlaps sit ~0.5-0.8 while
anion/cation sit ~0, so an unscaled mean implicitly up-weights high-baseline types);
(2) a color type absent from ALL train actives has zero train variance -> StandardScaler
sets scale_=1 and passes its raw [0,1] test value through unscaled, contaminating the
mean; those columns carry no train signal, so they are excluded from the equal-weight
mean. Templates still from TRAIN actives -> leakage-safe. Honest label: a fixed-weight
ablation of OUR features ('SHAFTS-style'), not literal SHAFTS (not a Python lib)."""
from __future__ import annotations

import os
import sys

import numpy as np
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData
from prism_core import prism_features


@register("prism_fixed")
def prism_fixed(data: BenchData, train_idx, test_idx) -> np.ndarray:
    X, tr = prism_features(data, train_idx, with_esp=False)
    if X is None:
        return np.zeros(len(test_idx), dtype=float)
    cols = X[tr].std(axis=0) > 1e-9                      # only columns with train signal
    sc = StandardScaler().fit(X[tr])
    Xz = sc.transform(X[np.asarray(test_idx)])
    return Xz[:, cols].mean(axis=1) if cols.any() else Xz.mean(axis=1)

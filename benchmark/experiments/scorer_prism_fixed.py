"""prism_fixed — fixed-weight ablation of PRISM (isolates the learned weighting).

IDENTICAL per-(template×type) 3D-color features as prism, aggregated with EQUAL fixed
weights (mean) instead of a learned logistic. The mean is taken AFTER the same
StandardScaler(fit on train) that prism feeds its logistic — so the ONLY difference
between prism and prism_fixed is learned-vs-fixed weights on the *same standardized
features*. (A raw mean of the unstandardized Tanimoto columns would conflate weighting
with scaling: hydrophobe/rings overlaps sit ~0.5-0.8 while anion/cation sit ~0, so an
unscaled mean implicitly up-weights the high-baseline types — not 'equal weight'.)
Templates still from TRAIN actives -> leakage-safe. Honest label: a fixed-weight ablation
of OUR features ('SHAFTS-style'), not literal SHAFTS (not a Python lib)."""
from __future__ import annotations
import os, sys
import numpy as np
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData, SEED
from prism_core import make_templates, feature_matrix

_NCONF = int(os.environ.get("BAKEOFF_NCONF", "5"))
_NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))
_MAX_TEMPLATES = int(os.environ.get("BAKEOFF_MAX_TEMPLATES", "8"))
_ALPHA = float(os.environ.get("BAKEOFF_COLOR_ALPHA", "0.5"))


@register("prism_fixed")
def prism_fixed(data: BenchData, train_idx, test_idx) -> np.ndarray:
    tr = np.asarray(train_idx)
    train_act = [data.smiles[i] for i in tr if data.y[i] == 1]
    tmpl_mols = make_templates(train_act, nconf=_NCONF, max_templates=_MAX_TEMPLATES, seed=SEED)
    if not tmpl_mols:
        return np.zeros(len(test_idx), dtype=float)
    X = feature_matrix(data, tmpl_mols, with_esp=False, nconf=_NCONF, njobs=_NJOBS, alpha=_ALPHA)
    sc = StandardScaler().fit(X[tr])                      # same standardization as prism
    return sc.transform(X[np.asarray(test_idx)]).mean(axis=1)   # equal fixed weights, no learning

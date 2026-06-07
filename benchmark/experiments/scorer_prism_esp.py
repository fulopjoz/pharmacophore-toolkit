"""prism_esp — PRISM plus an electrostatic-similarity channel. Same templates + per-
(template×type) color features as prism, plus a per-template Carbo ESP scalar (RDKit
Gasteiger charges on the same color-optimized pose) -> (K*7) -> logistic. Tests whether
electrostatics adds discrimination beyond shape+color. Honest scope: an ESP-similarity
FEATURE on a color-optimized pose, not field-based alignment (ShaEP/EON)."""
from __future__ import annotations
import os, sys
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData, SEED
from prism_core import make_templates, feature_matrix

_NCONF = int(os.environ.get("BAKEOFF_NCONF", "5"))
_NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))
_MAX_TEMPLATES = int(os.environ.get("BAKEOFF_MAX_TEMPLATES", "8"))
_ALPHA = float(os.environ.get("BAKEOFF_COLOR_ALPHA", "0.5"))


@register("prism_esp")
def prism_esp(data: BenchData, train_idx, test_idx) -> np.ndarray:
    tr = np.asarray(train_idx)
    train_act = [data.smiles[i] for i in tr if data.y[i] == 1]
    tmpl_mols = make_templates(train_act, nconf=_NCONF, max_templates=_MAX_TEMPLATES, seed=SEED)
    if not tmpl_mols:
        return np.zeros(len(test_idx), dtype=float)
    X = feature_matrix(data, tmpl_mols, with_esp=True, nconf=_NCONF, njobs=_NJOBS, alpha=_ALPHA)
    sc = StandardScaler().fit(X[tr])
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
    lr.fit(sc.transform(X[tr]), data.y[tr])
    return lr.predict_proba(sc.transform(X[np.asarray(test_idx)]))[:, 1]

"""PRISM — Pharmacophore-Resolved Importance-weighted Similarity Model (was `s3_3d`).

Discrimination-weighted 3D color, feature-type × template. Per held-out split:
  1. cluster the TRAIN actives -> K templates (leakage-safe), embed;
  2. for every molecule, align to each template and decompose the color into 6
     per-feature-type overlaps -> (K*6) feature matrix (prism_core);
  3. logistic regression (StandardScaler + LogisticRegression) fit on train rows
     learns the per-(type × template) discrimination weights and predicts test.

The logistic IS the discrimination weighting; the 3D-overlap features make it 3D.
The prism metaphor: a prism resolves white light into a spectrum, as PRISM resolves
ROCS "color" into its feature types and reweights them by learned importance.
Compare on the benchmark vs rdshape_ensemble (unweighted 3D color) and s3_weighted (2D).

Lineage (scite, 2026-06-07): weight-matrix learning over pharmacophore fingerprints
(Rehioui 2022, 10.1002/minf.202200210), PLS-DA on pharmacophore FPs (Askjaer 2008,
10.1021/ci700356w), PharmRF ML scoring (Kumar 2022, 10.1002/jcc.26840), multi-template
VS (Madzhidov 2020, 10.3390/molecules25020385). Novelty: 3D color decomposed by
feature type, validated on unbiased MUBD decoys (overfitting caveat: Wallach 2011
10.1021/ci100374f -> controlled by the held-out scaffold split + MUBD gate)."""
from __future__ import annotations

import os
import sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData, SEED  # noqa: E402
from prism_core import make_templates, feature_matrix  # noqa: E402

_NCONF = int(os.environ.get("BAKEOFF_NCONF", "5"))
_NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))
_MAX_TEMPLATES = int(os.environ.get("BAKEOFF_MAX_TEMPLATES", "8"))
_ALPHA = float(os.environ.get("BAKEOFF_COLOR_ALPHA", "0.5"))


@register("prism")
def prism(data: BenchData, train_idx, test_idx) -> np.ndarray:
    tr = np.asarray(train_idx)
    train_act = [data.smiles[i] for i in tr if data.y[i] == 1]
    tmpl_mols = make_templates(train_act, nconf=_NCONF, max_templates=_MAX_TEMPLATES, seed=SEED)
    if not tmpl_mols:
        return np.zeros(len(test_idx), dtype=float)
    X = feature_matrix(data, tmpl_mols, with_esp=False, nconf=_NCONF, njobs=_NJOBS, alpha=_ALPHA)
    sc = StandardScaler().fit(X[tr])
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
    lr.fit(sc.transform(X[tr]), data.y[tr])
    return lr.predict_proba(sc.transform(X[np.asarray(test_idx)]))[:, 1]

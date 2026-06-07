"""Discrimination-weighted 3D color (S3-in-3D), feature-type x template.

Per held-out split:
  1. cluster the TRAIN actives -> K templates (leakage-safe), embed + annotate;
  2. for every molecule, align to each template and decompose into 6 per-type
     color overlaps -> (K*6) feature matrix;
  3. logistic regression (StandardScaler + LogisticRegression) fit on train rows
     learns the per-(type x template) discrimination weights and predicts test.

The logistic IS the S3 weighting; the 3D-overlap features make it 3D. Compare on
the benchmark vs rdshape_ensemble (unweighted 3D color) and s3_weighted (2D).

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
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData, SEED, P4  # noqa: E402
from templates import cluster_templates  # noqa: E402
from color3d import embed, align_decompose  # noqa: E402

_NCONF = int(os.environ.get("BAKEOFF_NCONF", "5"))
_NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))
_MAX_TEMPLATES = int(os.environ.get("BAKEOFF_MAX_TEMPLATES", "8"))
_ALPHA = float(os.environ.get("BAKEOFF_COLOR_ALPHA", "0.5"))


def _row(smi, tmpl_mols):
    """(K*6) feature row for one molecule: per-type color overlaps vs each template."""
    q = embed(smi, _NCONF, SEED)
    if q is None:
        return np.zeros(len(tmpl_mols) * len(P4), dtype=float)
    return np.concatenate([align_decompose(q, t, _ALPHA) for t in tmpl_mols])


@register("s3_3d")
def s3_3d(data: BenchData, train_idx, test_idx) -> np.ndarray:
    tr = np.asarray(train_idx)
    train_act = [data.smiles[i] for i in tr if data.y[i] == 1]
    templates = cluster_templates(train_act, sim_cutoff=0.65,
                                  max_templates=_MAX_TEMPLATES, seed=SEED)
    # Templates use a few conformers each (kept small to bound cost).
    tmpl_mols = [m for m in (embed(s, max(1, _NCONF // 2), SEED) for s in templates) if m]
    if not tmpl_mols:
        return np.zeros(len(test_idx), dtype=float)

    rows = Parallel(n_jobs=_NJOBS, backend="loky")(
        delayed(_row)(data.smiles[i], tmpl_mols) for i in range(len(data.smiles)))
    X = np.vstack(rows)

    sc = StandardScaler().fit(X[tr])
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
    lr.fit(sc.transform(X[tr]), data.y[tr])
    return lr.predict_proba(sc.transform(X[np.asarray(test_idx)]))[:, 1]

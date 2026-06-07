"""Shared core for the PRISM family (Pharmacophore-Resolved Importance-weighted
Similarity Model; was the internal codename `s3_3d`).

Provides leakage-safe templates + the per-(template×type) 3D-color feature matrix
(optionally with an electrostatic-similarity column per template). Factored out so
`prism` / `prism_fixed` / `prism_esp` share one tested feature builder (DRY).

The `with_esp=True` path calls `color3d.align_decompose(..., with_esp=True)`, which
is added by the ESP work item; it is dormant (never executed) when `with_esp=False`."""
from __future__ import annotations

import os
import sys

import numpy as np
from joblib import Parallel, delayed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import BenchData, SEED, P4
from templates import cluster_templates
from color3d import embed, align_decompose


def make_templates(train_act_smiles, nconf: int, max_templates: int = 8, seed: int = SEED):
    """Cluster TRAIN actives -> embedded template mols (leakage-safe)."""
    smis = cluster_templates(train_act_smiles, sim_cutoff=0.65,
                             max_templates=max_templates, seed=seed)
    return [m for m in (embed(s, max(1, nconf // 2), seed) for s in smis) if m]


def _row(smi, tmpl_mols, with_esp, nconf, alpha):
    """(K*(6|7)) feature row for one molecule: per-type color (+esp) overlaps per template."""
    q = embed(smi, nconf, SEED)
    width = len(P4) + (1 if with_esp else 0)
    if q is None:
        return np.zeros(len(tmpl_mols) * width, dtype=float)
    if with_esp:
        return np.concatenate([align_decompose(q, t, alpha, with_esp=True) for t in tmpl_mols])
    return np.concatenate([align_decompose(q, t, alpha) for t in tmpl_mols])


def feature_matrix(data: BenchData, tmpl_mols, with_esp: bool, nconf: int,
                   njobs: int, alpha: float = 0.5) -> np.ndarray:
    """(n_mols, K*(6 or 7)) per-(template×type) color (+esp) overlap matrix."""
    rows = Parallel(n_jobs=njobs, backend="loky")(
        delayed(_row)(data.smiles[i], tmpl_mols, with_esp, nconf, alpha)
        for i in range(len(data.smiles)))
    return np.vstack(rows)

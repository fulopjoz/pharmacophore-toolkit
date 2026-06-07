"""Shared core for the PRISM family (Pharmacophore-Resolved Importance-weighted
Similarity Model; was the internal codename `s3_3d`).

Single source of truth for: the env-tunable config, leakage-safe templates, the
per-(template×type) 3D-color feature matrix (optionally + an electrostatic-similarity
column per template), and the shared scorer preamble `prism_features`. The three
scorers (prism / prism_fixed / prism_esp) differ ONLY in aggregation, so everything
up to the feature matrix lives here (DRY — avoids config drift between scorers).

The `with_esp=True` path calls `color3d.align_decompose(..., with_esp=True)`."""
from __future__ import annotations

import os
import sys

import numpy as np
from joblib import Parallel, delayed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import BenchData, SEED, P4
from templates import cluster_templates
from color3d import embed, align_decompose

# Single source of truth for the env-tunable config (was duplicated per scorer).
NCONF = int(os.environ.get("BAKEOFF_NCONF", "5"))
NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))
MAX_TEMPLATES = int(os.environ.get("BAKEOFF_MAX_TEMPLATES", "8"))
ALPHA = float(os.environ.get("BAKEOFF_COLOR_ALPHA", "0.5"))


def make_templates(train_act_smiles, nconf: int, max_templates: int = 8, seed: int = SEED):
    """Cluster TRAIN actives -> embedded template mols (leakage-safe)."""
    smis = cluster_templates(train_act_smiles, sim_cutoff=0.65,
                             max_templates=max_templates, seed=seed)
    return [m for m in (embed(s, max(1, nconf // 2), seed) for s in smis) if m]


def _row(smi, tmpl_mols, with_esp, nconf, alpha, seed):
    """(K*(6|7)) feature row for one molecule: per-type color (+esp) overlaps per template."""
    q = embed(smi, nconf, seed)
    width = len(P4) + (1 if with_esp else 0)
    if q is None:
        return np.zeros(len(tmpl_mols) * width, dtype=float)
    return np.concatenate([align_decompose(q, t, alpha, with_esp=with_esp) for t in tmpl_mols])


def feature_matrix(data: BenchData, tmpl_mols, with_esp: bool, nconf: int,
                   njobs: int, alpha: float = ALPHA, seed: int = SEED) -> np.ndarray:
    """(n_mols, K*(6 or 7)) per-(template×type) color (+esp) overlap matrix."""
    rows = Parallel(n_jobs=njobs, backend="loky")(
        delayed(_row)(data.smiles[i], tmpl_mols, with_esp, nconf, alpha, seed)
        for i in range(len(data.smiles)))
    return np.vstack(rows)


def prism_features(data: BenchData, train_idx, with_esp: bool):
    """Shared scorer preamble: TRAIN-derived templates + the (n, K*(6|7)) feature matrix.

    Returns ``(X, tr)`` where ``tr = np.asarray(train_idx)``, or ``(None, tr)`` if no
    templates could be built (caller then returns a zero score vector).

    The (expensive, conformer-heavy) feature matrix is cached on ``data._cache`` keyed by
    ``with_esp`` + the train split, so the three prism scorers built on the SAME data and
    split reuse it: prism + prism_fixed (both with_esp=False) share ONE build — which both
    avoids ~2x redundant 3D work in the bake-off AND guarantees the ablation compares
    byte-identical features."""
    tr = np.asarray(train_idx)
    ckey = "_prism_X_esp" if with_esp else "_prism_X_color"
    cached = data._cache.get(ckey)
    if cached is not None and np.array_equal(cached[0], tr):
        return cached[1], tr

    train_act = [data.smiles[i] for i in tr if data.y[i] == 1]
    tmpl_mols = make_templates(train_act, nconf=NCONF, max_templates=MAX_TEMPLATES, seed=SEED)
    if not tmpl_mols:
        return None, tr
    X = feature_matrix(data, tmpl_mols, with_esp=with_esp, nconf=NCONF, njobs=NJOBS,
                       alpha=ALPHA, seed=SEED)
    data._cache[ckey] = (tr.copy(), X)
    return X, tr

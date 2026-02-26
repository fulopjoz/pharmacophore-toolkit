"""Optimal transport scoring for pharmacophore model comparison.

Computes the earth mover's distance (Wasserstein distance) between two
pharmacophore models treated as discrete measures over typed 3D points.

The cost between two feature points combines spatial distance and type mismatch:
    c(i,j) = ||pos_i - pos_j||^2 + alpha * type_distance(i,j)

Implementation hierarchy:
1. If `pot` (Python Optimal Transport) is installed: use ot.emd() for exact OT
   or ot.sinkhorn() for fast entropy-regularized OT. Handles unequal sizes natively.
2. Fallback: scipy.optimize.linear_sum_assignment with dummy-padded cost matrix.

Source: Zhang et al. 2022 (WIREs), Gachon et al. 2025 (Cytometry).
"""

import numpy as np
from typing import List, Optional
from scipy.optimize import linear_sum_assignment

from .constants import get_type_distance

try:
    import ot as pot_lib
    HAS_POT = True
except ImportError:
    HAS_POT = False


def _build_cost_matrix(
    features_a: List[List],
    features_b: List[List],
    alpha: float = 0.5,
    max_distance: float = 10.0
) -> np.ndarray:
    """Build pairwise cost matrix between feature sets.

    Args:
        features_a: First pharmacophore model features [type, idx, x, y, z].
        features_b: Second pharmacophore model features.
        alpha: Weight for type mismatch (0=spatial only, 1=type only).
        max_distance: Cap on spatial distance (Angstroms).

    Returns:
        Cost matrix of shape (n_a, n_b).
    """
    n_a = len(features_a)
    n_b = len(features_b)
    cost = np.zeros((n_a, n_b), dtype=np.float64)

    for i, fa in enumerate(features_a):
        pos_a = np.array([fa[2], fa[3], fa[4]])
        for j, fb in enumerate(features_b):
            pos_b = np.array([fb[2], fb[3], fb[4]])
            spatial = min(np.sum((pos_a - pos_b) ** 2), max_distance ** 2)
            type_dist = get_type_distance(fa[0], fb[0])
            cost[i, j] = (1.0 - alpha) * spatial + alpha * type_dist

    return cost


def wasserstein_pharmacophore_distance(
    features_a: List[List],
    features_b: List[List],
    alpha: float = 0.5,
    max_distance: float = 10.0,
    use_sinkhorn: bool = False,
    sinkhorn_reg: float = 0.1
) -> float:
    """Compute earth mover's distance between two pharmacophore models.

    Each model is treated as a uniform discrete measure over typed 3D points.

    Args:
        features_a: First pharmacophore model features.
        features_b: Second pharmacophore model features.
        alpha: Weight for type distance vs spatial distance (0-1).
        max_distance: Spatial distance cap (Angstroms).
        use_sinkhorn: Use entropy-regularized OT (faster, approximate).
            Requires `pot` library.
        sinkhorn_reg: Regularization parameter for Sinkhorn (lower = more exact).

    Returns:
        Distance (lower = more similar). Normalized to [0, 1].
    """
    if not features_a and not features_b:
        return 0.0
    if not features_a or not features_b:
        return 1.0

    n_a = len(features_a)
    n_b = len(features_b)

    cost = _build_cost_matrix(features_a, features_b, alpha, max_distance)

    # Normalize cost to [0, 1] range
    max_cost = (1.0 - alpha) * max_distance ** 2 + alpha * 1.0
    if max_cost > 0:
        cost_norm = cost / max_cost
    else:
        cost_norm = cost

    if HAS_POT:
        # Uniform weights (all features equally important)
        a = np.ones(n_a) / n_a
        b = np.ones(n_b) / n_b

        if use_sinkhorn:
            transport_plan = pot_lib.sinkhorn(a, b, cost_norm, reg=sinkhorn_reg)
        else:
            transport_plan = pot_lib.emd(a, b, cost_norm)

        distance = np.sum(transport_plan * cost_norm)
    else:
        # Scipy fallback: pad to equal size, use Hungarian algorithm
        n_max = max(n_a, n_b)
        padded = np.ones((n_max, n_max), dtype=np.float64)  # dummy cost = 1.0

        padded[:n_a, :n_b] = cost_norm

        row_ind, col_ind = linear_sum_assignment(padded)
        total = padded[row_ind, col_ind].sum()
        distance = total / n_max

    return min(1.0, max(0.0, distance))


def wasserstein_similarity(
    features_a: List[List],
    features_b: List[List],
    alpha: float = 0.5,
    max_distance: float = 10.0
) -> float:
    """Compute Wasserstein similarity (1 - distance) between pharmacophore models.

    Convenience wrapper that returns similarity in [0, 1] (higher = more similar),
    compatible with other scoring functions in the toolkit.

    Args:
        features_a: First pharmacophore model features.
        features_b: Second pharmacophore model features.
        alpha: Weight for type distance vs spatial distance (0-1).
        max_distance: Spatial distance cap (Angstroms).

    Returns:
        Similarity in [0, 1]. 1 = identical, 0 = maximally different.
    """
    return 1.0 - wasserstein_pharmacophore_distance(
        features_a, features_b, alpha=alpha, max_distance=max_distance
    )


def wasserstein_similarity_aligned(
    query_mol,
    ref_mol,
    blend_alpha: float = 0.3,
    ot_alpha: float = 0.5,
    max_distance: float = 10.0,
    opt_param: float = 0.5,
    max_preiters: int = 10,
    max_postiters: int = 30,
) -> float:
    """Score query against reference using aligned OT pharmacophore distance.

    Aligns query onto ref via RDKit shape alignment, then computes Wasserstein
    pharmacophore distance on the aligned features. This ensures both feature
    sets share the same coordinate frame, fixing the random-AUC bug of the
    unaligned ``wasserstein_similarity()``.

    Args:
        query_mol: Query molecule with 3D conformer (a copy is made internally).
        ref_mol: Reference molecule with 3D conformer.
        blend_alpha: Blending weight for shape quality vs OT similarity.
            Final = blend_alpha * shape_quality + (1-blend_alpha) * ot_sim.
        ot_alpha: Weight for type distance vs spatial distance in OT cost (0-1).
        max_distance: Spatial distance cap (Angstroms).
        opt_param: AlignMol balance (0=color, 0.5=balanced, 1=shape).
        max_preiters: AlignMol pre-optimization iterations.
        max_postiters: AlignMol post-optimization iterations.

    Returns:
        Similarity in [0, 1]. Higher = more similar.
    """
    from rdkit import Chem
    from .shape_alignment import align_and_extract_features

    # Work on a copy so we don't mutate the caller's molecule
    probe = Chem.RWMol(query_mol)

    try:
        ref_feats, query_feats, shape_tani, color_tani = align_and_extract_features(
            ref_mol, probe,
            opt_param=opt_param,
            max_preiters=max_preiters,
            max_postiters=max_postiters,
        )
    except (ValueError, RuntimeError):
        return 0.0

    if not ref_feats or not query_feats:
        return 0.0

    # OT distance on aligned features
    dist = wasserstein_pharmacophore_distance(
        ref_feats, query_feats,
        alpha=ot_alpha,
        max_distance=max_distance,
    )
    ot_sim = 1.0 - dist

    # Shape quality from alignment (range 0-1 each)
    shape_quality = (shape_tani + color_tani) / 2.0

    # Blend: mostly OT signal, some shape quality
    final = blend_alpha * shape_quality + (1.0 - blend_alpha) * ot_sim
    return max(0.0, min(1.0, final))

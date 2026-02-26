"""Optimal feature assignment between pharmacophore models via Hungarian algorithm.

Uses the Hungarian algorithm (scipy.optimize.linear_sum_assignment) to find the
minimum-cost bijective assignment between features of two pharmacophore models.

The cost between two features combines spatial distance and type compatibility:
    cost(i,j) = spatial_weight * ||pos_i - pos_j||^2 + type_weight * type_dist(i,j)

Source: G3PS algorithm (Seidel et al. 2021), graph theory.
"""

import numpy as np
from typing import List, Tuple, Optional
from scipy.optimize import linear_sum_assignment

from .constants import get_type_distance


def _extract_coords_and_types(features: List[List]) -> Tuple[np.ndarray, List[str]]:
    """Extract coordinates and types from feature list.

    Args:
        features: List of [type, indices, x, y, z].

    Returns:
        Tuple of (coords array (N,3), type list).
    """
    coords = np.array([[f[2], f[3], f[4]] for f in features])
    types = [f[0] for f in features]
    return coords, types


def build_cost_matrix(
    features_a: List[List],
    features_b: List[List],
    spatial_weight: float = 1.0,
    type_weight: float = 1.0,
    max_distance: float = 10.0
) -> np.ndarray:
    """Build cost matrix for feature assignment.

    Args:
        features_a: First pharmacophore model features.
        features_b: Second pharmacophore model features.
        spatial_weight: Weight for spatial distance component.
        type_weight: Weight for type mismatch component.
        max_distance: Maximum spatial distance before capping (Angstroms).

    Returns:
        Cost matrix of shape (max(n_a, n_b), max(n_a, n_b)).
        Padded with dummy entries (cost=max_distance^2 + type_weight)
        when feature sets differ in size.
    """
    coords_a, types_a = _extract_coords_and_types(features_a)
    coords_b, types_b = _extract_coords_and_types(features_b)

    n_a = len(features_a)
    n_b = len(features_b)
    n_max = max(n_a, n_b)

    # Dummy cost for unmatched features
    dummy_cost = spatial_weight * max_distance ** 2 + type_weight

    cost = np.full((n_max, n_max), dummy_cost, dtype=np.float64)

    for i in range(n_a):
        for j in range(n_b):
            spatial_dist_sq = np.sum((coords_a[i] - coords_b[j]) ** 2)
            spatial_dist_sq = min(spatial_dist_sq, max_distance ** 2)
            type_dist = get_type_distance(types_a[i], types_b[j])
            cost[i, j] = spatial_weight * spatial_dist_sq + type_weight * type_dist

    return cost


def match_features(
    features_a: List[List],
    features_b: List[List],
    spatial_weight: float = 1.0,
    type_weight: float = 1.0,
    max_distance: float = 10.0
) -> Tuple[List[Tuple[int, int]], List[int], List[int], float]:
    """Find optimal assignment between two pharmacophore feature sets.

    Uses the Hungarian algorithm for minimum-cost bipartite matching.
    Handles unequal sizes by padding with high-cost dummy features.

    Args:
        features_a: First pharmacophore model features.
        features_b: Second pharmacophore model features.
        spatial_weight: Weight for spatial distance.
        type_weight: Weight for type mismatch.
        max_distance: Maximum spatial distance cutoff (Angstroms).

    Returns:
        Tuple of:
            - matched_pairs: List of (i, j) index pairs.
            - unmatched_a: Indices in features_a with no match.
            - unmatched_b: Indices in features_b with no match.
            - total_cost: Sum of matched pair costs.
    """
    if not features_a or not features_b:
        return (
            [],
            list(range(len(features_a))),
            list(range(len(features_b))),
            0.0
        )

    cost = build_cost_matrix(
        features_a, features_b,
        spatial_weight=spatial_weight,
        type_weight=type_weight,
        max_distance=max_distance
    )

    row_ind, col_ind = linear_sum_assignment(cost)

    n_a = len(features_a)
    n_b = len(features_b)
    dummy_cost = spatial_weight * max_distance ** 2 + type_weight

    matched_pairs = []
    unmatched_a = set(range(n_a))
    unmatched_b = set(range(n_b))
    total_cost = 0.0

    for r, c in zip(row_ind, col_ind):
        if r < n_a and c < n_b:
            pair_cost = cost[r, c]
            # Only count as matched if cost is below dummy threshold
            if pair_cost < dummy_cost * 0.95:
                matched_pairs.append((r, c))
                unmatched_a.discard(r)
                unmatched_b.discard(c)
                total_cost += pair_cost

    return matched_pairs, sorted(unmatched_a), sorted(unmatched_b), total_cost


def pharmacophore_distance(
    features_a: List[List],
    features_b: List[List],
    spatial_weight: float = 1.0,
    type_weight: float = 1.0,
    max_distance: float = 10.0
) -> float:
    """Compute distance between two pharmacophore models via optimal matching.

    Normalized to [0, 1] range for compatibility with other scoring methods.

    Args:
        features_a: First pharmacophore model.
        features_b: Second pharmacophore model.
        spatial_weight: Weight for spatial distance.
        type_weight: Weight for type mismatch.
        max_distance: Maximum spatial distance (Angstroms).

    Returns:
        Normalized distance in [0, 1]. 0 = identical, 1 = maximally different.
    """
    if not features_a and not features_b:
        return 0.0
    if not features_a or not features_b:
        return 1.0

    matched_pairs, unmatched_a, unmatched_b, total_cost = match_features(
        features_a, features_b,
        spatial_weight=spatial_weight,
        type_weight=type_weight,
        max_distance=max_distance
    )

    n_max = max(len(features_a), len(features_b))
    max_possible_cost = n_max * (spatial_weight * max_distance ** 2 + type_weight)

    if max_possible_cost < 1e-12:
        return 0.0

    # Add penalty for unmatched features
    dummy_cost = spatial_weight * max_distance ** 2 + type_weight
    penalty = (len(unmatched_a) + len(unmatched_b)) * dummy_cost
    normalized = (total_cost + penalty) / max_possible_cost

    return min(1.0, normalized)


def pharmacophore_similarity_aligned(
    query_mol,
    ref_mol,
    alpha: float = 0.3,
    spatial_weight: float = 1.0,
    type_weight: float = 1.0,
    max_distance: float = 10.0,
    opt_param: float = 0.5,
    max_preiters: int = 10,
    max_postiters: int = 30,
) -> float:
    """Score query against reference using aligned pharmacophore features.

    Aligns query onto ref via RDKit shape alignment, then computes Hungarian
    pharmacophore distance on the aligned features. This ensures both feature
    sets share the same coordinate frame, fixing the random-AUC bug of the
    unaligned ``pharmacophore_distance()``.

    Args:
        query_mol: Query molecule with 3D conformer (a copy is made internally).
        ref_mol: Reference molecule with 3D conformer.
        alpha: Blending weight for shape quality vs Hungarian similarity.
            Final = alpha * shape_quality + (1-alpha) * hungarian_sim.
        spatial_weight: Weight for spatial distance in cost matrix.
        type_weight: Weight for type mismatch in cost matrix.
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

    # Hungarian distance on aligned features
    dist = pharmacophore_distance(
        ref_feats, query_feats,
        spatial_weight=spatial_weight,
        type_weight=type_weight,
        max_distance=max_distance,
    )
    hungarian_sim = 1.0 - dist

    # Shape quality from alignment (range 0-1 each)
    shape_quality = (shape_tani + color_tani) / 2.0

    # Blend: mostly Hungarian signal, some shape quality
    final = alpha * shape_quality + (1.0 - alpha) * hungarian_sim
    return max(0.0, min(1.0, final))

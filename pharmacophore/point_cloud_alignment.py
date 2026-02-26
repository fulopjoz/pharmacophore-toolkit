"""Point cloud ICP alignment for pharmacophore model comparison.

Implements a colored Iterative Closest Point (ICP) algorithm that combines
spatial distance with pharmacophore feature-type distance for alignment and
scoring of pharmacophore point clouds.

The combined distance follows Zhou, Griffith & Gaeta (BMC Bioinformatics 2014):
    D(i,j) = (1-lambda) * ||pos_i - pos_j|| + lambda * type_distance(i,j)

Key refinements over the original paper:
    - Centroid normalization before ICP (center at origin)
    - Graceful handling of mismatched pharmacophore sizes
    - Gap penalty for unmatched features (PhAST-inspired)
    - RMSD-based convergence rather than fixed iteration count
    - scipy.linalg.orthogonal_procrustes for Kabsch rotation

Dependencies: scipy only (no new packages).
"""

import logging
import numpy as np
from typing import List, Tuple, Optional

from scipy.linalg import orthogonal_procrustes
from scipy.spatial.distance import cdist

from .constants import get_type_distance

logger = logging.getLogger(__name__)


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


def _build_type_distance_matrix(
    types_a: List[str],
    types_b: List[str],
) -> np.ndarray:
    """Build pairwise type distance matrix between two feature type lists.

    Args:
        types_a: Feature types for point set A.
        types_b: Feature types for point set B.

    Returns:
        Type distance matrix of shape (n_a, n_b), values in [0, 1].
    """
    n_a = len(types_a)
    n_b = len(types_b)
    type_dist = np.zeros((n_a, n_b), dtype=np.float64)
    for i in range(n_a):
        for j in range(n_b):
            type_dist[i, j] = get_type_distance(types_a[i], types_b[j])
    return type_dist


def _build_colored_distance_matrix(
    coords_a: np.ndarray,
    types_a: List[str],
    coords_b: np.ndarray,
    types_b: List[str],
    lambda_color: float = 0.5,
) -> np.ndarray:
    """Build colored distance matrix combining spatial and type distances.

    D(i,j) = (1-lambda) * euclidean(i,j) + lambda * type_distance(i,j)

    Spatial distances are normalized by dividing by the maximum observed
    distance to put them on the same scale as type distances (0-1).

    Args:
        coords_a: Coordinates for set A, shape (n_a, 3).
        types_a: Feature types for set A.
        coords_b: Coordinates for set B, shape (n_b, 3).
        types_b: Feature types for set B.
        lambda_color: Weight for chemical/type distance (0=spatial only, 1=type only).

    Returns:
        Combined distance matrix of shape (n_a, n_b).
    """
    spatial = cdist(coords_a, coords_b, metric='euclidean')

    # Normalize spatial distances to [0, 1] for fair combination with type distances
    max_spatial = spatial.max() if spatial.size > 0 and spatial.max() > 0 else 1.0
    spatial_norm = spatial / max_spatial

    type_dist = _build_type_distance_matrix(types_a, types_b)

    return (1.0 - lambda_color) * spatial_norm + lambda_color * type_dist


def colored_icp_align(
    features_a: List[List],
    features_b: List[List],
    lambda_color: float = 0.5,
    max_iter: int = 50,
    convergence_threshold: float = 1e-4,
) -> Tuple[np.ndarray, np.ndarray, float, List[Tuple[int, int]]]:
    """Align two pharmacophore point clouds using colored ICP.

    Centers both point sets at origin, then iteratively:
    1. Builds colored distance matrix (spatial + type)
    2. Finds nearest-neighbor correspondences
    3. Computes optimal rotation via Kabsch (orthogonal Procrustes)
    4. Applies rotation, checks RMSD convergence

    Handles different-sized point sets by matching the smaller set to
    the larger set (only the smaller set's points get assigned).

    Args:
        features_a: First pharmacophore model [type, idx, x, y, z].
        features_b: Second pharmacophore model [type, idx, x, y, z].
        lambda_color: Weight for type distance in combined metric.
        max_iter: Maximum ICP iterations.
        convergence_threshold: RMSD change threshold for convergence.

    Returns:
        Tuple of:
            - aligned_coords_a: Aligned coordinates of set A (N_a, 3).
            - rotation: Rotation matrix R (3, 3).
            - rmsd: Final RMSD of matched points.
            - correspondences: List of (i_a, i_b) matched pairs.
    """
    if not features_a or not features_b:
        return np.empty((0, 3)), np.eye(3), float('inf'), []

    coords_a, types_a = _extract_coords_and_types(features_a)
    coords_b, types_b = _extract_coords_and_types(features_b)

    # Center both sets at origin (Zhou et al. 2014: normalize centroids to (0,0,0))
    centroid_a = coords_a.mean(axis=0)
    centroid_b = coords_b.mean(axis=0)
    coords_a_c = coords_a - centroid_a
    coords_b_c = coords_b - centroid_b

    # Working copy that gets rotated each iteration
    working_a = coords_a_c.copy()
    R_total = np.eye(3)

    prev_rmsd = float('inf')

    for iteration in range(max_iter):
        # Build colored distance matrix
        D = _build_colored_distance_matrix(
            working_a, types_a, coords_b_c, types_b, lambda_color
        )

        # Nearest-neighbor assignment: for each point in A, find closest in B
        nn_b = D.argmin(axis=1)  # shape (n_a,)

        # Build matched point sets for Kabsch
        matched_a = working_a
        matched_b = coords_b_c[nn_b]

        # Kabsch rotation via orthogonal Procrustes
        # Finds R that minimizes ||matched_a @ R - matched_b||
        R, _ = orthogonal_procrustes(matched_a, matched_b)
        R_total = R_total @ R

        # Apply rotation
        working_a = working_a @ R

        # Compute RMSD of matched points
        diffs = working_a - matched_b
        rmsd = np.sqrt((diffs ** 2).sum(axis=1).mean())

        # Check convergence
        if abs(prev_rmsd - rmsd) < convergence_threshold:
            break
        prev_rmsd = rmsd

    # Final correspondences: re-compute with final alignment
    D_final = _build_colored_distance_matrix(
        working_a, types_a, coords_b_c, types_b, lambda_color
    )
    nn_final = D_final.argmin(axis=1)
    correspondences = [(i, int(nn_final[i])) for i in range(len(features_a))]

    # Translate aligned coords back to B's frame
    aligned_a = working_a + centroid_b

    return aligned_a, R_total, rmsd, correspondences


def point_cloud_similarity(
    features_a: List[List],
    features_b: List[List],
    lambda_color: float = 0.5,
    max_iter: int = 50,
    sigma: float = 2.0,
    gap_penalty: float = 1.0,
) -> float:
    """Compute similarity between two pharmacophore models via ICP alignment.

    Runs colored ICP to align the point clouds, then combines alignment
    quality (RMSD-based) with matched fraction to produce a similarity score.

    Handles mismatched sizes via a matched-fraction penalty following
    Zhou et al.'s observation that size ratios > 2:1 degrade performance.

    Args:
        features_a: First pharmacophore model features.
        features_b: Second pharmacophore model features.
        lambda_color: Weight for type distance in ICP.
        max_iter: Maximum ICP iterations.
        sigma: RMSD decay parameter (Angstroms). Controls how quickly
            similarity drops with increasing RMSD. Default 2.0 A.
        gap_penalty: Penalty factor per unmatched feature (0-1 range).

    Returns:
        Similarity in [0, 1]. Higher = more similar.
    """
    if not features_a and not features_b:
        return 1.0
    if not features_a or not features_b:
        return 0.0

    n_a = len(features_a)
    n_b = len(features_b)

    aligned_a, R, rmsd, correspondences = colored_icp_align(
        features_a, features_b,
        lambda_color=lambda_color,
        max_iter=max_iter,
    )

    if not correspondences or rmsd == float('inf'):
        return 0.0

    # Quality from RMSD: Gaussian decay
    quality = np.exp(-(rmsd ** 2) / (2.0 * sigma ** 2))

    # Matched fraction: penalize size mismatch
    # Only min(n_a, n_b) features can be matched; unmatched get gap penalty
    n_matched = min(n_a, n_b)
    n_total = max(n_a, n_b)
    n_unmatched = n_total - n_matched
    matched_frac = n_matched / n_total

    # Type-aware matching quality: check how well types align
    _, types_a = _extract_coords_and_types(features_a)
    _, types_b = _extract_coords_and_types(features_b)
    type_match_score = 0.0
    for i_a, i_b in correspondences:
        td = get_type_distance(types_a[i_a], types_b[i_b])
        type_match_score += 1.0 - td
    type_match_score /= len(correspondences)

    # Combine: alignment quality * matched fraction * type matching
    # Gap penalty reduces score for unmatched features
    gap_factor = 1.0 - gap_penalty * (n_unmatched / n_total)
    gap_factor = max(0.0, gap_factor)

    similarity = quality * type_match_score * gap_factor
    return max(0.0, min(1.0, similarity))


def point_cloud_similarity_aligned(
    query_mol,
    ref_mol,
    alpha: float = 0.3,
    lambda_color: float = 0.5,
    sigma: float = 2.0,
    opt_param: float = 0.5,
    max_preiters: int = 10,
    max_postiters: int = 30,
) -> float:
    """Score query against reference using shape alignment + ICP.

    First aligns query to ref via RDKit AlignMol (establishes coordinate
    frame), then runs colored ICP on the aligned pharmacophore features.

    Args:
        query_mol: Query molecule with 3D conformer (copy made internally).
        ref_mol: Reference molecule with 3D conformer.
        alpha: Blending weight for shape quality vs ICP similarity.
            Final = alpha * shape_quality + (1-alpha) * icp_sim.
        lambda_color: Weight for type distance in ICP combined metric.
        sigma: RMSD decay parameter for ICP similarity.
        opt_param: AlignMol balance (0=color, 0.5=balanced, 1=shape).
        max_preiters: AlignMol pre-optimization iterations.
        max_postiters: AlignMol post-optimization iterations.

    Returns:
        Similarity in [0, 1]. Higher = more similar.
    """
    from rdkit import Chem
    from .shape_alignment import align_and_extract_features

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

    # ICP similarity on aligned features
    icp_sim = point_cloud_similarity(
        ref_feats, query_feats,
        lambda_color=lambda_color,
        sigma=sigma,
    )

    # Shape quality from alignment
    shape_quality = (shape_tani + color_tani) / 2.0

    # Blend
    final = alpha * shape_quality + (1.0 - alpha) * icp_sim
    return max(0.0, min(1.0, final))

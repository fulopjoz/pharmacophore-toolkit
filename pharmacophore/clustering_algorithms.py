"""Alternative clustering algorithms for consensus pharmacophore generation.

This module provides different clustering approaches to compare against the
default hierarchical clustering:
1. K-means clustering (fast, assumes spherical clusters)
2. Grid-based binning (ultra-fast, grid discretization)
"""

import numpy as np
from typing import List, Tuple, Dict
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import warnings


def cluster_kmeans(
    coords: np.ndarray,
    tolerance: float,
    occurrence_threshold: float,
    n_molecules: int
) -> List[np.ndarray]:
    """Cluster features using K-means algorithm.
    
    K-means assumes spherical clusters and requires pre-specifying number
    of clusters. We estimate n_clusters from spatial extent and tolerance.
    
    Args:
        coords: Array of feature coordinates (N, 3)
        tolerance: Spatial radius for cluster membership (Å)
        occurrence_threshold: Minimum fraction of molecules
        n_molecules: Total number of molecules
    
    Returns:
        List of cluster centroids that meet occurrence threshold
    """
    if len(coords) == 0:
        return []
    
    if len(coords) == 1:
        return [coords[0]]
    
    # Estimate number of clusters from spatial extent
    # Assume average cluster spacing = 2 * tolerance
    extent = np.ptp(coords, axis=0).max()  # Max range
    n_clusters_est = max(1, int(extent / (2 * tolerance)))
    n_clusters = min(n_clusters_est, len(coords))
    
    # K-means clustering (deterministic with random_state)
    kmeans = KMeans(
        n_clusters=n_clusters,
        random_state=42,
        n_init=10,
        max_iter=300
    )
    
    try:
        labels = kmeans.fit_predict(coords)
        centroids = kmeans.cluster_centers_
    except Exception as e:
        warnings.warn(f"K-means failed: {e}. Returning all features.")
        return list(coords)
    
    # Count molecules per cluster
    # Assume each coordinate represents one feature from one molecule
    # Simple approach: each cluster needs enough features
    min_features = int(np.ceil(occurrence_threshold * n_molecules))
    
    consensus_centroids = []
    for cluster_id in range(n_clusters):
        cluster_mask = labels == cluster_id
        cluster_size = np.sum(cluster_mask)
        
        if cluster_size >= min_features:
            consensus_centroids.append(centroids[cluster_id])
    
    return consensus_centroids


def cluster_grid_binning(
    coords: np.ndarray,
    tolerance: float,
    occurrence_threshold: float,
    n_molecules: int
) -> List[np.ndarray]:
    """Cluster features using grid-based binning (ultra-fast).
    
    This method discretizes 3D space into cubic bins of size 2*tolerance.
    Features in the same bin are merged. This is the fastest method but
    assumes uniform spatial distribution.
    
    Advantages:
    - O(N) complexity (vs O(N²) for hierarchical)
    - Completely deterministic
    - No iterative optimization
    
    Disadvantages:
    - Sensitive to grid alignment
    - Can split natural clusters at bin boundaries
    
    Args:
        coords: Array of feature coordinates (N, 3)
        tolerance: Half-width of grid bins (Å)
        occurrence_threshold: Minimum fraction of molecules
        n_molecules: Total number of molecules
    
    Returns:
        List of bin centroids that meet occurrence threshold
    """
    if len(coords) == 0:
        return []
    
    if len(coords) == 1:
        return [coords[0]]
    
    # Create grid bins
    bin_size = 2 * tolerance
    
    # Shift coordinates to ensure positive bins
    coords_shifted = coords - coords.min(axis=0) + bin_size
    
    # Assign to bins
    bin_indices = (coords_shifted / bin_size).astype(int)
    
    # Group by bin (use tuple as dict key)
    bin_dict: Dict[Tuple[int, int, int], List[np.ndarray]] = {}
    
    for i, coord in enumerate(coords):
        bin_key = tuple(bin_indices[i])
        if bin_key not in bin_dict:
            bin_dict[bin_key] = []
        bin_dict[bin_key].append(coord)
    
    # Calculate centroids for bins that meet occurrence threshold
    min_features = int(np.ceil(occurrence_threshold * n_molecules))
    
    consensus_centroids = []
    for bin_coords in bin_dict.values():
        if len(bin_coords) >= min_features:
            centroid = np.mean(bin_coords, axis=0)
            consensus_centroids.append(centroid)
    
    return consensus_centroids


def cluster_features_with_labels(
    coords: np.ndarray,
    tolerance: float,
    occurrence_threshold: float,
    n_molecules: int,
    method: str = 'hierarchical',
    linkage: str = 'complete',
    filter_by_occurrence: bool = True,
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Cluster features and return native labels alongside centroids.

    Unlike ``cluster_features_generic`` (which only returns centroids), this
    function preserves the per-point cluster labels produced by each algorithm.
    This avoids post-hoc ``cdist + argmin`` reassignment that erases the
    differences between clustering methods.

    Args:
        coords: Array of feature coordinates (N, 3).
        tolerance: Clustering distance threshold (Angstroms).
        occurrence_threshold: Minimum fraction of molecules (0.0-1.0).
        n_molecules: Total number of input molecules.
        method: Clustering algorithm ('hierarchical', 'kmeans', 'grid').
        linkage: For hierarchical only ('average', 'complete', 'single', 'ward').
        filter_by_occurrence: If True (default), mark points in small clusters
            as -1 and only return centroids for surviving clusters.  If False,
            return ALL cluster labels (contiguous, >= 0) and all centroids —
            the caller is responsible for occurrence filtering.

    Returns:
        Tuple of (labels, centroids) where:
            - labels: int array of length N, cluster id per point
              (-1 for points in clusters below occurrence threshold
               when ``filter_by_occurrence=True``)
            - centroids: list of centroid arrays for surviving clusters
              (or all clusters when ``filter_by_occurrence=False``)
    """
    if len(coords) == 0:
        return np.array([], dtype=int), []

    if len(coords) == 1:
        return np.array([0]), [coords[0]]

    min_features = int(np.ceil(occurrence_threshold * n_molecules))

    if method == 'kmeans':
        extent = np.ptp(coords, axis=0).max()
        n_clusters_est = max(1, int(extent / (2 * tolerance)))
        n_clusters = min(n_clusters_est, len(coords))

        kmeans = KMeans(
            n_clusters=n_clusters, random_state=42,
            n_init=10, max_iter=300,
        )
        try:
            raw_labels = kmeans.fit_predict(coords)
        except Exception as e:
            warnings.warn(f"K-means failed: {e}. Returning all features.")
            return np.zeros(len(coords), dtype=int), list(coords)

        if not filter_by_occurrence:
            # Return raw contiguous labels and all centroids
            centroids = [kmeans.cluster_centers_[cid] for cid in range(n_clusters)]
            return raw_labels.astype(int), centroids

        # Filter by occurrence and remap labels
        labels = np.full(len(coords), -1, dtype=int)
        centroids = []
        new_id = 0
        for cid in range(n_clusters):
            mask = raw_labels == cid
            if mask.sum() >= min_features:
                labels[mask] = new_id
                centroids.append(kmeans.cluster_centers_[cid])
                new_id += 1
        return labels, centroids

    elif method == 'grid':
        bin_size = 2 * tolerance
        coords_shifted = coords - coords.min(axis=0) + bin_size
        bin_indices = (coords_shifted / bin_size).astype(int)

        # Map each point to its bin key
        bin_keys = [tuple(row) for row in bin_indices]

        # Group points by bin
        from collections import OrderedDict
        bin_groups: Dict[Tuple[int, int, int], List[int]] = OrderedDict()
        for i, key in enumerate(bin_keys):
            if key not in bin_groups:
                bin_groups[key] = []
            bin_groups[key].append(i)

        if not filter_by_occurrence:
            # Return all bins as clusters, contiguous labels
            labels = np.zeros(len(coords), dtype=int)
            centroids = []
            for cid, indices in enumerate(bin_groups.values()):
                for idx in indices:
                    labels[idx] = cid
                centroids.append(np.mean(coords[indices], axis=0))
            return labels, centroids

        labels = np.full(len(coords), -1, dtype=int)
        centroids = []
        new_id = 0
        for indices in bin_groups.values():
            if len(indices) >= min_features:
                for idx in indices:
                    labels[idx] = new_id
                centroids.append(np.mean(coords[indices], axis=0))
                new_id += 1
        return labels, centroids

    elif method == 'hierarchical':
        from sklearn.cluster import AgglomerativeClustering

        clustering = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=tolerance,
            linkage=linkage,
            metric='euclidean',
        )
        raw_labels = clustering.fit_predict(coords)

        if not filter_by_occurrence:
            # Return raw labels and all centroids
            centroids = []
            for cid in sorted(set(raw_labels)):
                centroids.append(np.mean(coords[raw_labels == cid], axis=0))
            return raw_labels.astype(int), centroids

        labels = np.full(len(coords), -1, dtype=int)
        centroids = []
        new_id = 0
        for cid in sorted(set(raw_labels)):
            mask = raw_labels == cid
            if mask.sum() >= min_features:
                labels[mask] = new_id
                centroids.append(np.mean(coords[mask], axis=0))
                new_id += 1
        return labels, centroids

    else:
        raise ValueError(
            f"Unknown clustering method: {method}. "
            f"Choose from: 'hierarchical', 'kmeans', 'grid'"
        )


def cluster_features_generic(
    coords: np.ndarray,
    tolerance: float,
    occurrence_threshold: float,
    n_molecules: int,
    method: str = 'hierarchical',
    linkage: str = 'complete'
) -> List[np.ndarray]:
    """Generic clustering dispatcher for different algorithms.
    
    This function routes to the appropriate clustering algorithm and
    returns consensus centroids.
    
    Args:
        coords: Array of feature coordinates (N, 3)
        tolerance: Clustering distance threshold (Å)
        occurrence_threshold: Minimum fraction of molecules (0.0-1.0)
        n_molecules: Total number of input molecules
        method: Clustering algorithm ('hierarchical', 'kmeans', 'grid')
        linkage: For hierarchical only ('average', 'complete', 'single', 'ward')
    
    Returns:
        List of consensus feature centroids
    
    Raises:
        ValueError: If method is not recognized
    """
    if method == 'kmeans':
        return cluster_kmeans(coords, tolerance, occurrence_threshold, n_molecules)
    
    elif method == 'dbscan':
        warnings.warn(
            "DBSCAN was removed for non-determinism (see CLAUDE.md). "
            "Falling back to 'hierarchical'. Update your code to use "
            "'hierarchical', 'kmeans', or 'grid'.",
            DeprecationWarning, stacklevel=2
        )
        return cluster_features_generic(
            coords, tolerance, occurrence_threshold, n_molecules,
            method='hierarchical', linkage=linkage
        )

    elif method == 'grid':
        return cluster_grid_binning(coords, tolerance, occurrence_threshold, n_molecules)
    
    elif method == 'hierarchical':
        # Use existing hierarchical clustering from consensus.py
        from sklearn.cluster import AgglomerativeClustering
        
        if len(coords) == 0:
            return []
        if len(coords) == 1:
            return [coords[0]]
        
        clustering = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=tolerance,
            linkage=linkage,
            metric='euclidean'
        )
        
        labels = clustering.fit_predict(coords)
        
        # Calculate centroids and filter by occurrence
        min_features = int(np.ceil(occurrence_threshold * n_molecules))
        
        consensus_centroids = []
        for label in set(labels):
            cluster_mask = labels == label
            cluster_size = np.sum(cluster_mask)
            
            if cluster_size >= min_features:
                cluster_coords = coords[cluster_mask]
                centroid = np.mean(cluster_coords, axis=0)
                consensus_centroids.append(centroid)
        
        return consensus_centroids
    
    else:
        raise ValueError(
            f"Unknown clustering method: {method}. "
            f"Choose from: 'hierarchical', 'kmeans', 'grid'"
        )


# Algorithm metadata for comparison
ALGORITHM_INFO = {
    'hierarchical': {
        'name': 'Hierarchical Clustering',
        'complexity': 'O(N² log N)',
        'advantages': ['Deterministic', 'Flexible linkage', 'No parameter tuning'],
        'disadvantages': ['Slow for large N', 'Memory intensive'],
        'best_for': 'High accuracy, small-medium datasets'
    },
    'kmeans': {
        'name': 'K-Means Clustering',
        'complexity': 'O(N × K × iter)',
        'advantages': ['Fast', 'Well-understood', 'Scales well'],
        'disadvantages': ['Requires estimating K', 'Assumes spherical clusters'],
        'best_for': 'Large datasets, uniform cluster sizes'
    },
    'grid': {
        'name': 'Grid-Based Binning',
        'complexity': 'O(N)',
        'advantages': ['Ultra-fast', 'Deterministic', 'Simple'],
        'disadvantages': ['Grid alignment artifacts', 'Fixed bin size'],
        'best_for': 'Very large datasets, speed priority'
    }
}

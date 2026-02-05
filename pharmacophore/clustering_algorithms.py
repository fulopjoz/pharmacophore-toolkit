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

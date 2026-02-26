"""Tests for clustering_algorithms module, including cluster_features_with_labels."""

import numpy as np
import pytest

from pharmacophore.clustering_algorithms import (
    cluster_features_with_labels,
    cluster_features_generic,
    cluster_kmeans,
    cluster_grid_binning,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def sample_coords():
    """Two tight clusters at (0,0,0) and (10,10,10), 5 points each."""
    rng = np.random.RandomState(42)
    cluster_a = rng.randn(5, 3) * 0.3  # Near origin
    cluster_b = rng.randn(5, 3) * 0.3 + 10.0  # Near (10,10,10)
    return np.vstack([cluster_a, cluster_b])


@pytest.fixture
def single_cluster_coords():
    """One tight cluster at (5,5,5)."""
    rng = np.random.RandomState(42)
    return rng.randn(8, 3) * 0.2 + 5.0


# ---------------------------------------------------------------------------
# cluster_features_with_labels: shape and semantics
# ---------------------------------------------------------------------------
class TestClusterFeaturesWithLabels:
    """Tests for the label-returning clustering wrapper."""

    def test_empty_input(self):
        labels, centroids = cluster_features_with_labels(
            np.array([]).reshape(0, 3), tolerance=2.0,
            occurrence_threshold=0.3, n_molecules=5,
        )
        assert len(labels) == 0
        assert centroids == []

    def test_single_point(self):
        coords = np.array([[1.0, 2.0, 3.0]])
        labels, centroids = cluster_features_with_labels(
            coords, tolerance=2.0,
            occurrence_threshold=0.3, n_molecules=5,
        )
        assert len(labels) == 1
        assert labels[0] == 0
        assert len(centroids) == 1
        np.testing.assert_array_almost_equal(centroids[0], [1.0, 2.0, 3.0])

    @pytest.mark.parametrize("method", ["hierarchical", "kmeans", "grid"])
    def test_label_shape_matches_input(self, sample_coords, method):
        labels, centroids = cluster_features_with_labels(
            sample_coords, tolerance=2.0,
            occurrence_threshold=0.1, n_molecules=10,
            method=method,
        )
        assert len(labels) == len(sample_coords)
        assert isinstance(labels, np.ndarray)
        assert labels.dtype == int

    @pytest.mark.parametrize("method", ["hierarchical", "kmeans", "grid"])
    def test_finds_two_clusters(self, sample_coords, method):
        """Two well-separated clusters should produce >= 2 centroids."""
        labels, centroids = cluster_features_with_labels(
            sample_coords, tolerance=2.0,
            occurrence_threshold=0.1, n_molecules=10,
            method=method,
        )
        assert len(centroids) >= 2, (
            f"{method} should find >= 2 clusters, got {len(centroids)}"
        )

    @pytest.mark.parametrize("method", ["hierarchical", "kmeans", "grid"])
    def test_labels_reference_valid_centroids(self, sample_coords, method):
        """Every non-negative label should index a valid centroid."""
        labels, centroids = cluster_features_with_labels(
            sample_coords, tolerance=2.0,
            occurrence_threshold=0.1, n_molecules=10,
            method=method,
        )
        n_centroids = len(centroids)
        for lbl in labels:
            assert lbl == -1 or 0 <= lbl < n_centroids

    @pytest.mark.parametrize("method", ["hierarchical", "kmeans", "grid"])
    def test_occurrence_filtering_marks_negative(self, method):
        """Points in small clusters should get label -1."""
        # 8 points near origin, 1 outlier at (100,100,100)
        coords = np.vstack([
            np.random.RandomState(0).randn(8, 3) * 0.2,
            [[100.0, 100.0, 100.0]],
        ])
        labels, centroids = cluster_features_with_labels(
            coords, tolerance=2.0,
            occurrence_threshold=0.3, n_molecules=10,
            method=method,
        )
        # The outlier (last point) should be filtered out
        assert labels[-1] == -1

    def test_unknown_method_raises(self, sample_coords):
        with pytest.raises(ValueError, match="Unknown clustering method"):
            cluster_features_with_labels(
                sample_coords, tolerance=2.0,
                occurrence_threshold=0.3, n_molecules=5,
                method='magic',
            )

    @pytest.mark.parametrize("method", ["hierarchical", "kmeans", "grid"])
    def test_labels_differ_between_clusters(self, sample_coords, method):
        """First 5 points (cluster A) and last 5 (cluster B) should have
        different labels."""
        labels, centroids = cluster_features_with_labels(
            sample_coords, tolerance=2.0,
            occurrence_threshold=0.1, n_molecules=10,
            method=method,
        )
        labels_a = set(labels[:5]) - {-1}
        labels_b = set(labels[5:]) - {-1}
        # The two groups should not share a label
        assert labels_a.isdisjoint(labels_b), (
            f"{method}: clusters share labels: A={labels_a}, B={labels_b}"
        )


# ---------------------------------------------------------------------------
# Consistency: with_labels centroids match generic centroids
# ---------------------------------------------------------------------------
class TestConsistencyWithGeneric:
    """cluster_features_with_labels centroids should match cluster_features_generic."""

    def test_hierarchical_centroids_match(self, sample_coords):
        labels, centroids_wl = cluster_features_with_labels(
            sample_coords, tolerance=2.0,
            occurrence_threshold=0.1, n_molecules=10,
            method='hierarchical',
        )
        centroids_gen = cluster_features_generic(
            sample_coords, tolerance=2.0,
            occurrence_threshold=0.1, n_molecules=10,
            method='hierarchical',
        )
        # Same number of centroids
        assert len(centroids_wl) == len(centroids_gen)

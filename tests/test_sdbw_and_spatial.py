"""Tests for S_Dbw cluster validation and spatial scaling."""

import pytest
import numpy as np

from pharmacophore.evaluation import compute_sdbw
from pharmacophore.constants import SPATIAL_SCALE, FEATURE_COMPATIBILITY, get_type_distance


class TestSpatialScale:
    """Test SPATIAL_SCALE constants."""

    def test_donor_tighter_than_baseline(self):
        assert SPATIAL_SCALE['Donor'] < 1.0

    def test_acceptor_tighter_than_baseline(self):
        assert SPATIAL_SCALE['Acceptor'] < 1.0

    def test_hydrophobe_looser_than_baseline(self):
        assert SPATIAL_SCALE['Hydrophobe'] > 1.0

    def test_aromatic_is_baseline(self):
        assert SPATIAL_SCALE['Aromatic'] == 1.0

    def test_pos_ionizable_tighter(self):
        assert SPATIAL_SCALE['PosIonizable'] < 1.0


class TestFeatureCompatibility:
    """Test the type compatibility matrix."""

    def test_donor_acceptor_compatible(self):
        assert 0.0 < get_type_distance('Donor', 'Acceptor') < 1.0

    def test_aromatic_hydrophobe_compatible(self):
        assert 0.0 < get_type_distance('Aromatic', 'Hydrophobe') < 1.0

    def test_hydrophobe_lumped_identical(self):
        assert get_type_distance('Hydrophobe', 'LumpedHydrophobe') == 0.0

    def test_donor_hydrophobe_incompatible(self):
        assert get_type_distance('Donor', 'Hydrophobe') == 1.0


class TestComputeSdbw:
    """Test S_Dbw cluster validation index."""

    def _make_metadata(self, coords_list, labels_list):
        """Helper to build metadata dict."""
        meta = {}
        for i, (coords, labels) in enumerate(zip(coords_list, labels_list)):
            meta[f'Type_{i}'] = {
                'coordinates': np.array(coords),
                'labels': np.array(labels),
                'cluster_to_mols': {},
                'centroids': {},
                'n_clusters': len(set(labels)),
                'n_valid': len(set(labels)),
            }
        return meta

    def test_well_separated_clusters(self):
        """Well-separated clusters should have low S_Dbw."""
        coords = [
            [0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.1, 0.1, 0.0],
            [10.0, 10.0, 10.0], [10.1, 10.0, 10.0], [10.0, 10.1, 10.0],
        ]
        labels = [0, 0, 0, 1, 1, 1]
        meta = self._make_metadata([coords], [labels])
        sdbw = compute_sdbw(meta)
        assert sdbw < 1.0, f"Expected low S_Dbw for well-separated clusters, got {sdbw}"

    def test_overlapping_clusters_higher(self):
        """Overlapping clusters should have higher S_Dbw."""
        coords = [
            [0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0],
            [0.5, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0],
        ]
        labels = [0, 0, 0, 1, 1, 1]
        meta = self._make_metadata([coords], [labels])
        sdbw_overlap = compute_sdbw(meta)

        # Well-separated
        coords2 = [
            [0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0],
            [100.0, 0.0, 0.0], [100.1, 0.0, 0.0], [100.2, 0.0, 0.0],
        ]
        labels2 = [0, 0, 0, 1, 1, 1]
        meta2 = self._make_metadata([coords2], [labels2])
        sdbw_sep = compute_sdbw(meta2)

        assert sdbw_overlap > sdbw_sep, (
            f"Overlapping ({sdbw_overlap}) should have higher S_Dbw than "
            f"separated ({sdbw_sep})"
        )

    def test_empty_metadata(self):
        sdbw = compute_sdbw({})
        assert sdbw == float('inf')

    def test_single_cluster(self):
        """Single cluster should return inf (need at least 2)."""
        coords = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]
        labels = [0, 0]
        meta = self._make_metadata([coords], [labels])
        sdbw = compute_sdbw(meta)
        assert sdbw == float('inf')

    def test_multiple_feature_types(self):
        """S_Dbw should work across multiple feature types."""
        coords1 = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [5.0, 5.0, 5.0], [5.1, 5.0, 5.0]]
        labels1 = [0, 0, 1, 1]
        coords2 = [[1.0, 1.0, 1.0], [1.1, 1.0, 1.0], [8.0, 8.0, 8.0], [8.1, 8.0, 8.0]]
        labels2 = [0, 0, 1, 1]
        meta = self._make_metadata([coords1, coords2], [labels1, labels2])
        sdbw = compute_sdbw(meta)
        assert sdbw != float('inf')
        assert sdbw >= 0.0

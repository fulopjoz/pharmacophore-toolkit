"""Tests for pharmacophore Hungarian matching module."""

import pytest
import numpy as np

from pharmacophore.hungarian_matching import (
    match_features,
    build_cost_matrix,
    pharmacophore_distance,
)
from pharmacophore.constants import get_type_distance


class TestGetTypeDistance:
    """Test the type compatibility function."""

    def test_same_type_returns_zero(self):
        assert get_type_distance('Donor', 'Donor') == 0.0
        assert get_type_distance('Acceptor', 'Acceptor') == 0.0

    def test_compatible_types(self):
        d = get_type_distance('Donor', 'Acceptor')
        assert 0.0 < d < 1.0, f"Expected 0 < {d} < 1 for compatible types"

    def test_incompatible_types(self):
        d = get_type_distance('Donor', 'Hydrophobe')
        assert d == 1.0

    def test_symmetric(self):
        assert get_type_distance('Donor', 'Acceptor') == get_type_distance('Acceptor', 'Donor')

    def test_lumped_hydrophobe_matches_hydrophobe(self):
        assert get_type_distance('Hydrophobe', 'LumpedHydrophobe') == 0.0


class TestBuildCostMatrix:
    """Test cost matrix construction."""

    def test_identical_features(self):
        features = [['Donor', (), 1.0, 2.0, 3.0]]
        cost = build_cost_matrix(features, features)
        assert cost[0, 0] == 0.0

    def test_same_type_different_position(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Donor', (), 1.0, 0.0, 0.0]]
        cost = build_cost_matrix(fa, fb, spatial_weight=1.0, type_weight=0.0)
        assert cost[0, 0] == pytest.approx(1.0)  # 1^2 = 1

    def test_different_type_same_position(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Hydrophobe', (), 0.0, 0.0, 0.0]]
        cost = build_cost_matrix(fa, fb, spatial_weight=0.0, type_weight=1.0)
        assert cost[0, 0] == 1.0  # incompatible types

    def test_unequal_sizes_padded(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0], ['Acceptor', (), 1.0, 0.0, 0.0]]
        fb = [['Donor', (), 0.0, 0.0, 0.0]]
        cost = build_cost_matrix(fa, fb)
        assert cost.shape == (2, 2)  # padded to square


class TestMatchFeatures:
    """Test optimal feature matching."""

    def test_identical_models(self):
        features = [
            ['Donor', (), 1.0, 2.0, 3.0],
            ['Acceptor', (), 4.0, 5.0, 6.0],
        ]
        pairs, unmatched_a, unmatched_b, cost = match_features(features, features)
        assert len(pairs) == 2
        assert len(unmatched_a) == 0
        assert len(unmatched_b) == 0
        assert cost == pytest.approx(0.0)

    def test_different_sizes(self):
        fa = [
            ['Donor', (), 0.0, 0.0, 0.0],
            ['Acceptor', (), 1.0, 0.0, 0.0],
            ['Aromatic', (), 2.0, 0.0, 0.0],
        ]
        fb = [
            ['Donor', (), 0.1, 0.0, 0.0],
            ['Acceptor', (), 1.1, 0.0, 0.0],
        ]
        pairs, unmatched_a, unmatched_b, cost = match_features(fa, fb)
        assert len(pairs) == 2  # Two matched
        assert len(unmatched_a) == 1  # One extra in A
        assert len(unmatched_b) == 0

    def test_empty_inputs(self):
        pairs, ua, ub, cost = match_features([], [])
        assert pairs == []
        assert cost == 0.0

    def test_one_empty(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        pairs, ua, ub, cost = match_features(fa, [])
        assert len(pairs) == 0
        assert len(ua) == 1

    def test_prefers_same_type(self):
        """Matching should prefer features of the same type."""
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [
            ['Donor', (), 0.5, 0.0, 0.0],     # Same type, close
            ['Hydrophobe', (), 0.1, 0.0, 0.0],  # Different type, closer
        ]
        pairs, _, _, _ = match_features(fa, fb, spatial_weight=1.0, type_weight=5.0)
        assert (0, 0) in pairs  # Should match Donor-Donor despite being farther


class TestPharmacophoreDistance:
    """Test normalized distance function."""

    def test_identical_returns_zero(self):
        features = [['Donor', (), 1.0, 2.0, 3.0]]
        d = pharmacophore_distance(features, features)
        assert d == pytest.approx(0.0)

    def test_very_different_returns_high(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Hydrophobe', (), 100.0, 100.0, 100.0]]
        d = pharmacophore_distance(fa, fb)
        assert d > 0.5

    def test_bounded_zero_one(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0], ['Acceptor', (), 5.0, 5.0, 5.0]]
        fb = [['Aromatic', (), 3.0, 3.0, 3.0]]
        d = pharmacophore_distance(fa, fb)
        assert 0.0 <= d <= 1.0

    def test_empty_same(self):
        assert pharmacophore_distance([], []) == 0.0

    def test_empty_different(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        assert pharmacophore_distance(fa, []) == 1.0
        assert pharmacophore_distance([], fa) == 1.0

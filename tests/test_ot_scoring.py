"""Tests for pharmacophore optimal transport scoring module."""

import pytest
import numpy as np

from pharmacophore.ot_scoring import (
    wasserstein_pharmacophore_distance,
    wasserstein_similarity,
    _build_cost_matrix,
    HAS_POT,
)


class TestBuildCostMatrix:
    """Test internal cost matrix construction."""

    def test_identical_features_zero_cost(self):
        f = [['Donor', (), 1.0, 2.0, 3.0]]
        cost = _build_cost_matrix(f, f)
        assert cost[0, 0] == pytest.approx(0.0)

    def test_spatial_only(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Donor', (), 1.0, 0.0, 0.0]]
        cost = _build_cost_matrix(fa, fb, alpha=0.0)  # spatial only
        assert cost[0, 0] == pytest.approx(1.0)  # 1^2 = 1

    def test_type_only(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Hydrophobe', (), 0.0, 0.0, 0.0]]
        cost = _build_cost_matrix(fa, fb, alpha=1.0)  # type only
        assert cost[0, 0] == pytest.approx(1.0)  # incompatible

    def test_shape(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0], ['Acceptor', (), 1.0, 0.0, 0.0]]
        fb = [['Donor', (), 0.0, 0.0, 0.0]]
        cost = _build_cost_matrix(fa, fb)
        assert cost.shape == (2, 1)


class TestWassersteinDistance:
    """Test Wasserstein pharmacophore distance."""

    def test_identical_is_zero(self):
        features = [
            ['Donor', (), 1.0, 2.0, 3.0],
            ['Acceptor', (), 4.0, 5.0, 6.0],
        ]
        d = wasserstein_pharmacophore_distance(features, features)
        assert d == pytest.approx(0.0, abs=0.01)

    def test_different_is_positive(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Hydrophobe', (), 5.0, 5.0, 5.0]]
        d = wasserstein_pharmacophore_distance(fa, fb)
        assert d > 0.0

    def test_bounded_zero_one(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Acceptor', (), 3.0, 3.0, 3.0]]
        d = wasserstein_pharmacophore_distance(fa, fb)
        assert 0.0 <= d <= 1.0

    def test_empty_both(self):
        assert wasserstein_pharmacophore_distance([], []) == 0.0

    def test_empty_one(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        assert wasserstein_pharmacophore_distance(fa, []) == 1.0
        assert wasserstein_pharmacophore_distance([], fa) == 1.0

    def test_unequal_sizes(self):
        fa = [
            ['Donor', (), 0.0, 0.0, 0.0],
            ['Acceptor', (), 1.0, 0.0, 0.0],
            ['Aromatic', (), 2.0, 0.0, 0.0],
        ]
        fb = [
            ['Donor', (), 0.1, 0.0, 0.0],
        ]
        d = wasserstein_pharmacophore_distance(fa, fb)
        assert 0.0 <= d <= 1.0

    def test_symmetric(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Acceptor', (), 2.0, 2.0, 2.0]]
        d1 = wasserstein_pharmacophore_distance(fa, fb)
        d2 = wasserstein_pharmacophore_distance(fb, fa)
        assert d1 == pytest.approx(d2, abs=0.05)

    def test_alpha_zero_spatial_only(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Donor', (), 1.0, 0.0, 0.0]]
        d = wasserstein_pharmacophore_distance(fa, fb, alpha=0.0)
        assert d > 0.0

    def test_alpha_one_type_only(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb_same = [['Donor', (), 5.0, 5.0, 5.0]]
        fb_diff = [['Hydrophobe', (), 0.0, 0.0, 0.0]]
        d_same = wasserstein_pharmacophore_distance(fa, fb_same, alpha=1.0)
        d_diff = wasserstein_pharmacophore_distance(fa, fb_diff, alpha=1.0)
        assert d_same < d_diff


class TestWassersteinSimilarity:
    """Test similarity convenience wrapper."""

    def test_identical_is_one(self):
        features = [['Donor', (), 1.0, 2.0, 3.0]]
        s = wasserstein_similarity(features, features)
        assert s == pytest.approx(1.0, abs=0.01)

    def test_inverse_of_distance(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Acceptor', (), 3.0, 3.0, 3.0]]
        d = wasserstein_pharmacophore_distance(fa, fb)
        s = wasserstein_similarity(fa, fb)
        assert s == pytest.approx(1.0 - d, abs=0.01)


@pytest.mark.skipif(not HAS_POT, reason="pot library not installed")
class TestPOTBackend:
    """Tests specific to the POT library backend."""

    def test_sinkhorn_runs(self):
        fa = [['Donor', (), 0.0, 0.0, 0.0]]
        fb = [['Donor', (), 1.0, 0.0, 0.0]]
        d = wasserstein_pharmacophore_distance(fa, fb, use_sinkhorn=True)
        assert 0.0 <= d <= 1.0

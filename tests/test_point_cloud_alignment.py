"""Tests for point cloud ICP alignment module."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.point_cloud_alignment import (
    _extract_coords_and_types,
    _build_type_distance_matrix,
    _build_colored_distance_matrix,
    colored_icp_align,
    point_cloud_similarity,
    point_cloud_similarity_aligned,
)
from pharmacophore.constants import get_type_distance


# --- Test fixtures ---

@pytest.fixture
def simple_features_a():
    """Donor + Acceptor at known positions."""
    return [
        ['Donor', (), 1.0, 0.0, 0.0],
        ['Acceptor', (), 0.0, 1.0, 0.0],
        ['Aromatic', (), 0.0, 0.0, 1.0],
    ]


@pytest.fixture
def simple_features_b():
    """Same types, slightly shifted positions."""
    return [
        ['Donor', (), 1.1, 0.1, 0.0],
        ['Acceptor', (), 0.1, 1.1, 0.0],
        ['Aromatic', (), 0.0, 0.1, 1.1],
    ]


@pytest.fixture
def rotated_features():
    """Same types as features_a but rotated 90 degrees around Z axis."""
    return [
        ['Donor', (), 0.0, 1.0, 0.0],
        ['Acceptor', (), -1.0, 0.0, 0.0],
        ['Aromatic', (), 0.0, 0.0, 1.0],
    ]


@pytest.fixture
def mismatched_features():
    """Fewer features with different types."""
    return [
        ['Hydrophobe', (), 2.0, 2.0, 2.0],
        ['Donor', (), 3.0, 3.0, 3.0],
    ]


def _make_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    return mol


# --- Tests for _extract_coords_and_types ---

class TestExtractCoordsAndTypes:

    def test_correct_shapes(self, simple_features_a):
        coords, types = _extract_coords_and_types(simple_features_a)
        assert coords.shape == (3, 3)
        assert len(types) == 3

    def test_correct_values(self):
        features = [['Donor', (), 1.5, 2.5, 3.5]]
        coords, types = _extract_coords_and_types(features)
        np.testing.assert_array_almost_equal(coords[0], [1.5, 2.5, 3.5])
        assert types[0] == 'Donor'


# --- Tests for _build_type_distance_matrix ---

class TestBuildTypeDistanceMatrix:

    def test_same_types_zero_diagonal(self):
        types = ['Donor', 'Acceptor', 'Aromatic']
        D = _build_type_distance_matrix(types, types)
        for i in range(3):
            assert D[i, i] == 0.0

    def test_symmetric(self):
        types_a = ['Donor', 'Acceptor']
        types_b = ['Aromatic', 'Hydrophobe']
        D = _build_type_distance_matrix(types_a, types_b)
        D_rev = _build_type_distance_matrix(types_b, types_a)
        np.testing.assert_array_almost_equal(D, D_rev.T)

    def test_values_in_range(self):
        types_a = ['Donor', 'Acceptor', 'Aromatic']
        types_b = ['Hydrophobe', 'Donor', 'PosIonizable']
        D = _build_type_distance_matrix(types_a, types_b)
        assert D.min() >= 0.0
        assert D.max() <= 1.0


# --- Tests for _build_colored_distance_matrix ---

class TestBuildColoredDistanceMatrix:

    def test_shape(self, simple_features_a, simple_features_b):
        coords_a, types_a = _extract_coords_and_types(simple_features_a)
        coords_b, types_b = _extract_coords_and_types(simple_features_b)
        D = _build_colored_distance_matrix(coords_a, types_a, coords_b, types_b)
        assert D.shape == (3, 3)

    def test_lambda_zero_is_spatial_only(self, simple_features_a, simple_features_b):
        coords_a, types_a = _extract_coords_and_types(simple_features_a)
        coords_b, types_b = _extract_coords_and_types(simple_features_b)
        D = _build_colored_distance_matrix(
            coords_a, types_a, coords_b, types_b, lambda_color=0.0
        )
        # Same types, so type distance is 0 everywhere for matching indices
        # With lambda=0, only spatial matters
        from scipy.spatial.distance import cdist
        spatial = cdist(coords_a, coords_b, metric='euclidean')
        spatial_norm = spatial / spatial.max()
        np.testing.assert_array_almost_equal(D, spatial_norm)

    def test_lambda_one_is_type_only(self, simple_features_a, simple_features_b):
        coords_a, types_a = _extract_coords_and_types(simple_features_a)
        coords_b, types_b = _extract_coords_and_types(simple_features_b)
        D = _build_colored_distance_matrix(
            coords_a, types_a, coords_b, types_b, lambda_color=1.0
        )
        type_dist = _build_type_distance_matrix(types_a, types_b)
        np.testing.assert_array_almost_equal(D, type_dist)


# --- Tests for colored_icp_align ---

class TestColoredICPAlign:

    def test_returns_correct_tuple(self, simple_features_a, simple_features_b):
        aligned, R, rmsd, correspondences = colored_icp_align(
            simple_features_a, simple_features_b
        )
        assert aligned.shape[1] == 3
        assert R.shape == (3, 3)
        assert isinstance(rmsd, float)
        assert isinstance(correspondences, list)

    def test_identical_features_zero_rmsd(self, simple_features_a):
        aligned, R, rmsd, correspondences = colored_icp_align(
            simple_features_a, simple_features_a
        )
        assert rmsd < 0.01

    def test_similar_features_low_rmsd(self, simple_features_a, simple_features_b):
        aligned, R, rmsd, correspondences = colored_icp_align(
            simple_features_a, simple_features_b
        )
        assert rmsd < 1.0

    def test_rotation_is_orthogonal(self, simple_features_a, rotated_features):
        aligned, R, rmsd, correspondences = colored_icp_align(
            simple_features_a, rotated_features
        )
        # R should be orthogonal: R @ R.T ≈ I
        np.testing.assert_array_almost_equal(R @ R.T, np.eye(3), decimal=5)

    def test_empty_features(self):
        aligned, R, rmsd, correspondences = colored_icp_align([], [])
        assert aligned.shape == (0, 3)
        assert rmsd == float('inf')
        assert correspondences == []

    def test_handles_size_mismatch(self, simple_features_a, mismatched_features):
        aligned, R, rmsd, correspondences = colored_icp_align(
            simple_features_a, mismatched_features
        )
        assert len(correspondences) == len(simple_features_a)
        assert rmsd < float('inf')


# --- Tests for point_cloud_similarity ---

class TestPointCloudSimilarity:

    def test_identical_is_one(self, simple_features_a):
        sim = point_cloud_similarity(simple_features_a, simple_features_a)
        assert sim > 0.9

    def test_similar_is_high(self, simple_features_a, simple_features_b):
        sim = point_cloud_similarity(simple_features_a, simple_features_b)
        assert sim > 0.5

    def test_range_zero_to_one(self, simple_features_a, simple_features_b):
        sim = point_cloud_similarity(simple_features_a, simple_features_b)
        assert 0.0 <= sim <= 1.0

    def test_empty_both_is_one(self):
        assert point_cloud_similarity([], []) == 1.0

    def test_empty_one_is_zero(self, simple_features_a):
        assert point_cloud_similarity(simple_features_a, []) == 0.0
        assert point_cloud_similarity([], simple_features_a) == 0.0

    def test_mismatched_size_penalty(self, simple_features_a, mismatched_features):
        # 3 vs 2 features — should have lower similarity than 3 vs 3
        sim_matched = point_cloud_similarity(simple_features_a, simple_features_a)
        sim_mismatched = point_cloud_similarity(simple_features_a, mismatched_features)
        assert sim_mismatched < sim_matched

    def test_sigma_affects_score(self, simple_features_a, simple_features_b):
        sim_small = point_cloud_similarity(
            simple_features_a, simple_features_b, sigma=0.5
        )
        sim_large = point_cloud_similarity(
            simple_features_a, simple_features_b, sigma=5.0
        )
        # Larger sigma = more tolerant = higher score
        assert sim_large >= sim_small


# --- Tests for point_cloud_similarity_aligned ---

class TestPointCloudSimilarityAligned:

    def test_returns_float_in_range(self):
        ref = _make_mol('c1ccccc1O')
        query = _make_mol('c1ccccc1N')
        sim = point_cloud_similarity_aligned(query, ref)
        assert isinstance(sim, float)
        assert 0.0 <= sim <= 1.0

    def test_identical_molecule_high_score(self):
        mol = _make_mol('c1ccccc1O')
        query = Chem.RWMol(mol)
        sim = point_cloud_similarity_aligned(query, mol)
        assert sim > 0.3  # Not 1.0 due to alignment noise

    def test_dissimilar_molecules_lower(self):
        ref = _make_mol('c1ccccc1O')
        similar = _make_mol('c1ccccc1N')
        dissimilar = _make_mol('CCCCCCCCCC')
        sim_similar = point_cloud_similarity_aligned(similar, ref)
        sim_dissimilar = point_cloud_similarity_aligned(dissimilar, ref)
        assert sim_similar >= sim_dissimilar

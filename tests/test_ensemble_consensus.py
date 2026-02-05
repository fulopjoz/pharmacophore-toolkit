"""Tests for ensemble consensus pharmacophore module."""

import pytest
import numpy as np
from unittest.mock import MagicMock

from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.ensemble_consensus import EnsembleConsensus


def _make_test_mol(smiles: str = "c1ccccc1O") -> Chem.Mol:
    """Create a test molecule with 3D conformer."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _make_test_mols(n: int = 3) -> list:
    """Create multiple test molecules with conformers."""
    smiles_list = [
        "c1ccccc1O",     # phenol
        "c1ccccc1N",     # aniline
        "c1ccccc1OC",    # anisole
        "c1ccccc1CC",    # ethylbenzene
        "c1ccc(O)cc1N",  # aminophenol
    ]
    mols = []
    for i in range(min(n, len(smiles_list))):
        mol = _make_test_mol(smiles_list[i])
        if mol is not None and mol.GetNumConformers() > 0:
            mols.append(mol)
    return mols


class TestEnsembleConsensusInit:
    """Test initialization."""

    def test_default_params(self):
        ec = EnsembleConsensus()
        assert ec.n_runs == 25
        assert ec.tolerance_range == (1.5, 2.5)
        assert ec.occurrence_range == (0.3, 0.7)
        assert ec.methods == ['hierarchical', 'kmeans']

    def test_custom_params(self):
        ec = EnsembleConsensus(
            n_runs=10,
            tolerance_range=(1.0, 3.0),
            methods=['hierarchical', 'grid'],
            stability_threshold=0.5
        )
        assert ec.n_runs == 10
        assert ec.stability_threshold == 0.5
        assert 'grid' in ec.methods


class TestEnsembleConsensusGeneration:
    """Test consensus generation."""

    @pytest.fixture
    def test_mols(self):
        return _make_test_mols(3)

    def test_generates_features(self, test_mols):
        if len(test_mols) < 2:
            pytest.skip("Need at least 2 molecules")
        ec = EnsembleConsensus(n_runs=5, stability_threshold=0.2)
        features = ec.generate_consensus(test_mols)
        assert isinstance(features, list)
        # Should produce at least some features from phenol/aniline/anisole
        # (all share aromatic ring + heteroatom)

    def test_feature_format(self, test_mols):
        if len(test_mols) < 2:
            pytest.skip("Need at least 2 molecules")
        ec = EnsembleConsensus(n_runs=5, stability_threshold=0.1)
        features = ec.generate_consensus(test_mols)
        for f in features:
            assert len(f) == 5, f"Feature should have 5 elements: {f}"
            assert isinstance(f[0], str), "Feature type should be string"
            assert f[1] == (), "Atom indices should be empty tuple"
            assert all(isinstance(f[i], float) for i in [2, 3, 4])

    def test_with_scores(self, test_mols):
        if len(test_mols) < 2:
            pytest.skip("Need at least 2 molecules")
        ec = EnsembleConsensus(n_runs=5, stability_threshold=0.1)
        features, scores = ec.generate_consensus_with_scores(test_mols)
        assert len(features) == len(scores)
        for s in scores:
            assert 0.0 < s <= 1.0, f"Stability score should be in (0, 1]: {s}"

    def test_deterministic(self, test_mols):
        if len(test_mols) < 2:
            pytest.skip("Need at least 2 molecules")
        ec = EnsembleConsensus(n_runs=5, random_state=42)
        f1 = ec.generate_consensus(test_mols)
        ec2 = EnsembleConsensus(n_runs=5, random_state=42)
        f2 = ec2.generate_consensus(test_mols)
        assert len(f1) == len(f2)
        for a, b in zip(f1, f2):
            assert a[0] == b[0]  # same type
            np.testing.assert_allclose(
                [a[2], a[3], a[4]],
                [b[2], b[3], b[4]],
                atol=0.01
            )

    def test_higher_stability_fewer_features(self, test_mols):
        if len(test_mols) < 2:
            pytest.skip("Need at least 2 molecules")
        ec_low = EnsembleConsensus(n_runs=10, stability_threshold=0.1)
        ec_high = EnsembleConsensus(n_runs=10, stability_threshold=0.8)
        f_low = ec_low.generate_consensus(test_mols)
        f_high = ec_high.generate_consensus(test_mols)
        assert len(f_high) <= len(f_low)


class TestVotingMechanism:
    """Test the feature voting logic."""

    def test_empty_runs(self):
        ec = EnsembleConsensus(n_runs=5)
        result = ec._vote_features([])
        assert result == []

    def test_all_empty_runs(self):
        ec = EnsembleConsensus(n_runs=3)
        result = ec._vote_features([[], [], []])
        assert result == []

    def test_perfect_agreement(self):
        ec = EnsembleConsensus(n_runs=3, stability_threshold=0.5)
        feature = ['Donor', (), 1.0, 2.0, 3.0]
        runs = [[feature], [feature], [feature]]
        result = ec._vote_features(runs)
        assert len(result) == 1
        assert result[0][0] == 'Donor'

    def test_partial_agreement(self):
        ec = EnsembleConsensus(
            n_runs=4,
            stability_threshold=0.5,
            hash_tolerance=0.5
        )
        f1 = ['Donor', (), 1.0, 2.0, 3.0]
        f2 = ['Donor', (), 1.1, 2.0, 3.0]  # close enough to match
        runs = [[f1], [f2], [], []]  # 2/4 = 0.5 stability
        result = ec._vote_features(runs)
        assert len(result) == 1  # 0.5 >= threshold

    def test_below_threshold_excluded(self):
        ec = EnsembleConsensus(
            n_runs=10,
            stability_threshold=0.5,
            hash_tolerance=0.5
        )
        feature = ['Donor', (), 1.0, 2.0, 3.0]
        runs = [[feature]] + [[] for _ in range(9)]  # 1/10 = 0.1 < 0.5
        result = ec._vote_features(runs)
        assert len(result) == 0

    def test_different_types_not_merged(self):
        ec = EnsembleConsensus(n_runs=2, stability_threshold=0.3, hash_tolerance=0.5)
        f1 = ['Donor', (), 1.0, 2.0, 3.0]
        f2 = ['Acceptor', (), 1.0, 2.0, 3.0]  # same position, different type
        runs = [[f1, f2], [f1, f2]]
        result = ec._vote_features(runs)
        types = [f[0] for f in result]
        assert 'Donor' in types
        assert 'Acceptor' in types

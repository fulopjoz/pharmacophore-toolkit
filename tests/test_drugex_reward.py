"""Tests for DrugEx reward function wrapper."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.drugex_reward import PharmacophoreReward


@pytest.fixture
def reference_mols():
    """Create reference molecules with 3D conformers."""
    smiles = ['CC(=O)O', 'CCO', 'c1ccccc1']
    mols = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol_h = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol_h, randomSeed=42)
            if mol_h.GetNumConformers() > 0:
                mols.append(mol_h)
    return mols


@pytest.fixture
def reward(reference_mols):
    return PharmacophoreReward(reference_mols, mode='pharm2d')


class TestPharmacophoreReward:
    """Tests for PharmacophoreReward."""

    def test_valid_smiles_returns_float(self, reward):
        score = reward.score('CCO')
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_invalid_smiles_returns_zero(self, reward):
        score = reward.score('INVALID_SMILES')
        assert score == 0.0

    def test_nonsense_smiles_returns_zero(self, reward):
        score = reward.score('NOT_A_MOLECULE_XYZ')
        assert score == 0.0

    def test_callable_interface(self, reward):
        """reward('CCO') should work like reward.score('CCO')."""
        score1 = reward('CCO')
        score2 = reward.score('CCO')
        assert score1 == score2

    def test_batch_scoring(self, reward):
        smiles_list = ['CCO', 'c1ccccc1', 'INVALID', 'CC(=O)O']
        scores = reward.score_batch(smiles_list)
        assert len(scores) == 4
        assert scores[2] == 0.0  # Invalid SMILES

    def test_mode_switching(self, reward):
        # Start in pharm2d
        assert reward.mode == 'pharm2d'

        # Switch to 3d
        reward.switch_mode('3d')
        assert reward.mode == '3d'

        # Switch to hybrid
        reward.switch_mode('hybrid')
        assert reward.mode == 'hybrid'

    def test_invalid_mode_raises(self, reference_mols):
        with pytest.raises(ValueError, match="mode must be"):
            PharmacophoreReward(reference_mols, mode='invalid')

    def test_invalid_switch_mode_raises(self, reward):
        with pytest.raises(ValueError, match="mode must be"):
            reward.switch_mode('invalid')

    def test_caching(self, reward):
        """Repeated SMILES should use cache."""
        _ = reward.score('CCO')
        _ = reward.score('CCO')  # Cache hit

        stats = reward.get_stats()
        assert stats['n_calls'] == 2
        assert stats['cache_hits'] >= 1

    def test_stats_tracking(self, reward):
        _ = reward.score('CCO')
        _ = reward.score('INVALID')

        stats = reward.get_stats()
        assert stats['n_calls'] == 2
        assert stats['n_valid'] >= 1
        assert stats['n_errors'] >= 1

    def test_repr(self, reward):
        r = repr(reward)
        assert 'PharmacophoreReward' in r
        assert "pharm2d" in r

    def test_3d_mode_scoring(self, reference_mols):
        reward = PharmacophoreReward(reference_mols, mode='3d', n_conformers=5)
        score = reward.score('CCO')
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_normalize(self, reference_mols):
        """Normalized scores should be in [0, 1]."""
        reward = PharmacophoreReward(
            reference_mols, mode='pharm2d', normalize=True
        )
        score = reward.score('CC(=O)O')
        assert 0.0 <= score <= 1.0

    def test_mode_switch_clears_cache(self, reward):
        _ = reward.score('CCO')
        stats1 = reward.get_stats()
        assert stats1['cache_size'] >= 1

        reward.switch_mode('3d')
        stats2 = reward.get_stats()
        assert stats2['cache_size'] == 0

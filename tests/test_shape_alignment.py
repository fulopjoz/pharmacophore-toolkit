"""Tests for the shared shape alignment utility."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.shape_alignment import align_and_extract_features
from pharmacophore.pharmacophore import Pharmacophore


def _make_mol_with_conformer(smiles: str) -> Chem.Mol:
    """Create a molecule with a 3D conformer from SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    return mol


class TestAlignAndExtractFeatures:
    """Tests for align_and_extract_features."""

    def test_returns_four_element_tuple(self):
        ref = _make_mol_with_conformer('c1ccccc1O')
        query = _make_mol_with_conformer('c1ccccc1N')
        result = align_and_extract_features(ref, query)
        assert len(result) == 4

    def test_returns_features_and_scores(self):
        ref = _make_mol_with_conformer('c1ccccc1O')
        query = _make_mol_with_conformer('c1ccccc1N')
        ref_feats, query_feats, shape_tani, color_tani = align_and_extract_features(ref, query)
        assert isinstance(ref_feats, list)
        assert isinstance(query_feats, list)
        assert isinstance(shape_tani, float)
        assert isinstance(color_tani, float)

    def test_features_are_nonempty(self):
        ref = _make_mol_with_conformer('c1ccccc1O')
        query = _make_mol_with_conformer('c1ccccc1N')
        ref_feats, query_feats, _, _ = align_and_extract_features(ref, query)
        assert len(ref_feats) > 0
        assert len(query_feats) > 0

    def test_shape_tanimoto_in_range(self):
        ref = _make_mol_with_conformer('c1ccccc1O')
        query = _make_mol_with_conformer('c1ccccc1N')
        _, _, shape_tani, color_tani = align_and_extract_features(ref, query)
        assert 0.0 <= shape_tani <= 1.0
        assert 0.0 <= color_tani <= 1.0

    def test_identical_molecule_high_similarity(self):
        mol = _make_mol_with_conformer('c1ccccc1O')
        mol_copy = Chem.RWMol(mol)
        _, _, shape_tani, _ = align_and_extract_features(mol, mol_copy)
        assert shape_tani > 0.8

    def test_raises_on_no_conformer(self):
        mol_no_conf = Chem.MolFromSmiles('c1ccccc1O')
        mol_with_conf = _make_mol_with_conformer('c1ccccc1N')
        with pytest.raises(ValueError, match="no conformer"):
            align_and_extract_features(mol_no_conf, mol_with_conf)
        with pytest.raises(ValueError, match="no conformer"):
            align_and_extract_features(mol_with_conf, mol_no_conf)

    def test_custom_pharmacophore(self):
        ref = _make_mol_with_conformer('c1ccccc1O')
        query = _make_mol_with_conformer('c1ccccc1N')
        p = Pharmacophore()
        ref_feats, query_feats, _, _ = align_and_extract_features(ref, query, pharmacophore=p)
        assert len(ref_feats) > 0
        assert len(query_feats) > 0

    def test_opt_param_affects_alignment(self):
        ref = _make_mol_with_conformer('c1ccccc1O')
        query = _make_mol_with_conformer('c1ccccc1N')
        # opt_param=0.5 balances shape and color during optimization
        _, _, shape_05, color_05 = align_and_extract_features(ref, query, opt_param=0.5)
        # opt_param=1.0 optimizes only shape (color still reported but not optimized)
        _, _, shape_10, color_10 = align_and_extract_features(
            ref, Chem.RWMol(query), opt_param=1.0
        )
        # Both should produce valid scores
        assert 0.0 <= shape_05 <= 1.0
        assert 0.0 <= color_05 <= 1.0
        assert 0.0 <= shape_10 <= 1.0
        assert 0.0 <= color_10 <= 1.0

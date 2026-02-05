"""Tests for FEATURE_TO_SMILES entries in mol_converter.

Verifies that each SMILES fragment:
1. Produces a valid RDKit molecule
2. Matches at least one SMARTS pattern for its feature type
3. Can embed a 3D conformer
"""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.constants import FEATURES


class TestFeatureSmiles:
    """Tests for FEATURE_TO_SMILES validity and SMARTS matching."""

    def test_acceptor_smiles_is_valid_molecule(self):
        """C=O should produce formaldehyde with an O atom."""
        mol = Chem.MolFromSmiles('C=O')
        assert mol is not None
        assert mol.GetNumAtoms() == 2
        # Should contain oxygen
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        assert 'O' in atoms

    def test_hydrophobe_smiles_is_valid_molecule(self):
        """C1CCCC1 should produce cyclopentane."""
        mol = Chem.MolFromSmiles('C1CCCC1')
        assert mol is not None
        assert mol.GetNumAtoms() == 5

    def test_acceptor_matches_acceptor_smarts(self):
        """C=O should match at least one Acceptor SMARTS pattern."""
        mol = Chem.MolFromSmiles('C=O')
        assert mol is not None
        matched = False
        for smarts in FEATURES['Acceptor']:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                matched = True
                break
        assert matched, "C=O does not match any Acceptor SMARTS"

    def test_hydrophobe_matches_hydrophobe_smarts(self):
        """C1CCCC1 should match at least one Hydrophobe SMARTS pattern."""
        mol = Chem.MolFromSmiles('C1CCCC1')
        assert mol is not None
        matched = False
        for smarts in FEATURES['Hydrophobe']:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                matched = True
                break
        assert matched, "C1CCCC1 does not match any Hydrophobe SMARTS"

    def test_donor_matches_donor_smarts(self):
        """[NH3] should match at least one Donor SMARTS pattern."""
        mol = Chem.MolFromSmiles('[NH3]')
        assert mol is not None
        matched = False
        for smarts in FEATURES['Donor']:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                matched = True
                break
        assert matched, "[NH3] does not match any Donor SMARTS"

    def test_aromatic_matches_aromatic_smarts(self):
        """c1ccccc1 should match at least one Aromatic SMARTS pattern."""
        mol = Chem.MolFromSmiles('c1ccccc1')
        assert mol is not None
        matched = False
        for smarts in FEATURES['Aromatic']:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                matched = True
                break
        assert matched, "c1ccccc1 does not match any Aromatic SMARTS"

    def test_all_feature_smiles_produce_valid_mols(self):
        """Every FEATURE_TO_SMILES entry should produce a valid RDKit Mol."""
        for feat_type, smiles in PharmacophoreToMol.FEATURE_TO_SMILES.items():
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None, f"{feat_type} SMILES '{smiles}' is invalid"

    def test_all_feature_smiles_embed_3d(self):
        """Every FEATURE_TO_SMILES entry should be embeddable in 3D."""
        for feat_type, smiles in PharmacophoreToMol.FEATURE_TO_SMILES.items():
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None
            mol_h = Chem.AddHs(mol)
            result = AllChem.EmbedMolecule(mol_h, randomSeed=42)
            assert result == 0, (
                f"{feat_type} SMILES '{smiles}' failed 3D embedding"
            )
            assert mol_h.GetNumConformers() > 0

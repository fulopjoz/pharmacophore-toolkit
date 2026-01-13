"""Unit tests for consensus pharmacophore implementation.

This test suite covers:
- PharmacophoreConsensus class
- PharmacophoreToMol converter
- Integration with main Pharmacophore class
- Determinism validation
- Backward compatibility
"""

import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.pharmacophore import Pharmacophore


class TestPharmacophoreConsensus(unittest.TestCase):
    """Test PharmacophoreConsensus class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.consensus = PharmacophoreConsensus(
            tolerance=2.0,
            occurrence_threshold=0.5,
            linkage='average'
        )
    
    def test_initialization_valid(self):
        """Test valid initialization parameters."""
        consensus = PharmacophoreConsensus(
            tolerance=1.5,
            occurrence_threshold=0.7,
            linkage='complete'
        )
        
        self.assertEqual(consensus.tolerance, 1.5)
        self.assertEqual(consensus.occurrence_threshold, 0.7)
        self.assertEqual(consensus.linkage, 'complete')
    
    def test_initialization_invalid_tolerance(self):
        """Test that negative tolerance raises ValueError."""
        with self.assertRaises(ValueError):
            PharmacophoreConsensus(tolerance=-1.0)
        
        with self.assertRaises(ValueError):
            PharmacophoreConsensus(tolerance=0.0)
    
    def test_initialization_invalid_threshold(self):
        """Test that invalid occurrence_threshold raises ValueError."""
        with self.assertRaises(ValueError):
            PharmacophoreConsensus(occurrence_threshold=-0.1)
        
        with self.assertRaises(ValueError):
            PharmacophoreConsensus(occurrence_threshold=1.5)
    
    def test_initialization_invalid_linkage(self):
        """Test that invalid linkage raises ValueError."""
        with self.assertRaises(ValueError):
            PharmacophoreConsensus(linkage='invalid')
    
    def test_cluster_features_single(self):
        """Test clustering with single feature."""
        coords = [np.array([1.0, 2.0, 3.0])]
        mol_indices = np.array([0])
        
        labels, cluster_to_mols = self.consensus._cluster_features(
            coordinates=coords,
            mol_indices=mol_indices
        )
        
        self.assertEqual(len(labels), 1)
        self.assertEqual(labels[0], 0)
        self.assertEqual(cluster_to_mols[0], [0])
    
    def test_cluster_features_close_points(self):
        """Test clustering of close points."""
        # Two points within tolerance
        coords = [
            np.array([1.0, 2.0, 3.0]),
            np.array([1.5, 2.1, 3.2])
        ]
        mol_indices = np.array([0, 1])
        
        labels, cluster_to_mols = self.consensus._cluster_features(
            coordinates=coords,
            mol_indices=mol_indices
        )
        
        # Should be in same cluster
        self.assertEqual(labels[0], labels[1])
    
    def test_cluster_features_far_points(self):
        """Test clustering of distant points."""
        # Two points far apart
        coords = [
            np.array([0.0, 0.0, 0.0]),
            np.array([10.0, 10.0, 10.0])
        ]
        mol_indices = np.array([0, 1])
        
        labels, cluster_to_mols = self.consensus._cluster_features(
            coordinates=coords,
            mol_indices=mol_indices
        )
        
        # Should be in different clusters
        self.assertNotEqual(labels[0], labels[1])
    
    def test_calculate_cluster_centroids(self):
        """Test centroid calculation."""
        coords = [
            np.array([0.0, 0.0, 0.0]),
            np.array([2.0, 2.0, 2.0])
        ]
        labels = np.array([0, 0])
        
        centroids = self.consensus._calculate_cluster_centroids(
            coordinates=coords,
            labels=labels
        )
        
        # Centroid should be average
        expected = np.array([1.0, 1.0, 1.0])
        np.testing.assert_array_almost_equal(centroids[0], expected)
    
    def test_filter_by_occurrence(self):
        """Test occurrence threshold filtering."""
        centroids = {
            0: np.array([1.0, 2.0, 3.0]),
            1: np.array([4.0, 5.0, 6.0])
        }
        cluster_to_mols = {
            0: [0, 1, 2],  # 3 molecules
            1: [0]         # 1 molecule
        }
        total_mols = 5
        
        # occurrence_threshold = 0.5, so need 3+ molecules
        valid = self.consensus._filter_by_occurrence(
            centroids=centroids,
            cluster_to_mols=cluster_to_mols,
            total_molecules=total_mols
        )
        
        # Only cluster 0 should pass (3/5 >= 0.5)
        self.assertEqual(len(valid), 1)
        np.testing.assert_array_almost_equal(valid[0][0], centroids[0])
        self.assertEqual(valid[0][1], 3)
    
    def test_determinism(self):
        """Test that clustering is deterministic."""
        # Create simple test molecules
        mol1 = Chem.MolFromSmiles('N')
        mol2 = Chem.MolFromSmiles('N')
        
        AllChem.EmbedMolecule(mol1, randomSeed=42)
        AllChem.EmbedMolecule(mol2, randomSeed=43)
        
        mols = [mol1, mol2]
        
        # Generate consensus twice
        features1 = self.consensus.generate_consensus(mols)
        features2 = self.consensus.generate_consensus(mols)
        
        # Should be identical
        self.assertEqual(len(features1), len(features2))
        
        for f1, f2 in zip(features1, features2):
            self.assertEqual(f1[0], f2[0])  # Same type
            np.testing.assert_array_almost_equal(
                [f1[2], f1[3], f1[4]],
                [f2[2], f2[3], f2[4]]
            )


class TestPharmacophoreToMol(unittest.TestCase):
    """Test PharmacophoreToMol converter."""
    
    def test_convert_single_feature(self):
        """Test conversion of single feature."""
        features = [['Donor', (), 1.0, 2.0, 3.0]]
        
        mol = PharmacophoreToMol.convert(features)
        
        self.assertEqual(mol.GetNumAtoms(), 1)
        self.assertEqual(mol.GetAtomWithIdx(0).GetSymbol(), 'N')
        
        # Check 3D coordinates
        conformer = mol.GetConformer()
        pos = conformer.GetAtomPosition(0)
        self.assertAlmostEqual(pos.x, 1.0)
        self.assertAlmostEqual(pos.y, 2.0)
        self.assertAlmostEqual(pos.z, 3.0)
    
    def test_convert_multiple_features(self):
        """Test conversion of multiple features."""
        features = [
            ['Donor', (), 1.0, 2.0, 3.0],
            ['Acceptor', (), 4.0, 5.0, 6.0],
            ['Aromatic', (), 7.0, 8.0, 9.0]
        ]
        
        mol = PharmacophoreToMol.convert(features)
        
        self.assertEqual(mol.GetNumAtoms(), 3)
        self.assertEqual(mol.GetAtomWithIdx(0).GetSymbol(), 'N')
        self.assertEqual(mol.GetAtomWithIdx(1).GetSymbol(), 'O')
        self.assertEqual(mol.GetAtomWithIdx(2).GetSymbol(), 'C')
    
    def test_convert_empty_features(self):
        """Test that empty features raises ValueError."""
        with self.assertRaises(ValueError):
            PharmacophoreToMol.convert([])
    
    def test_convert_invalid_feature(self):
        """Test that invalid feature raises ValueError."""
        features = [['InvalidType', (), 1.0, 2.0, 3.0]]
        
        with self.assertRaises(ValueError):
            PharmacophoreToMol.convert(features)
    
    def test_convert_with_metadata(self):
        """Test conversion with metadata."""
        features = [['Donor', (), 1.0, 2.0, 3.0]]
        
        mol = PharmacophoreToMol.convert_with_metadata(features)
        
        atom = mol.GetAtomWithIdx(0)
        self.assertEqual(atom.GetProp('FeatureType'), 'Donor')
    
    def test_validate_for_shape_alignment(self):
        """Test validation for shape alignment."""
        features = [['Donor', (), 1.0, 2.0, 3.0]]
        mol = PharmacophoreToMol.convert(features)
        
        # Should pass validation
        self.assertTrue(
            PharmacophoreToMol.validate_for_shape_alignment(mol)
        )
    
    def test_validate_none_molecule(self):
        """Test that None molecule raises ValueError."""
        with self.assertRaises(ValueError):
            PharmacophoreToMol.validate_for_shape_alignment(None)
    
    def test_feature_to_element_mapping(self):
        """Test all feature type mappings."""
        expected_mappings = {
            'Donor': 'N',
            'Acceptor': 'O',
            'Aromatic': 'C',
            'Hydrophobe': 'C',
            'LumpedHydrophobe': 'C',
            'PosIonizable': 'N'
        }
        
        for feat_type, expected_element in expected_mappings.items():
            features = [[feat_type, (), 1.0, 2.0, 3.0]]
            mol = PharmacophoreToMol.convert(features)
            
            actual_element = mol.GetAtomWithIdx(0).GetSymbol()
            self.assertEqual(actual_element, expected_element,
                           f"Failed for {feat_type}")


class TestPharmacophoreIntegration(unittest.TestCase):
    """Test integration with main Pharmacophore class."""
    
    def setUp(self):
        """Set up test molecules."""
        # Create simple molecules with 3D coordinates
        self.mol1 = Chem.MolFromSmiles('CC(=O)N')
        self.mol2 = Chem.MolFromSmiles('CC(=O)N')
        
        AllChem.EmbedMolecule(self.mol1, randomSeed=42)
        AllChem.EmbedMolecule(self.mol2, randomSeed=43)
        
        self.mols = [self.mol1, self.mol2]
    
    def test_generate_consensus_models_standard(self):
        """Test standard model set generation."""
        pharm = Pharmacophore()
        
        models = pharm.generate_consensus_models(
            mols=self.mols,
            model_set='standard'
        )
        
        # Should have 3 models
        self.assertEqual(len(models), 3)
        self.assertIn('strict', models)
        self.assertIn('moderate', models)
        self.assertIn('relaxed', models)
        
        # Each should be a tuple of (features, mol)
        for model_name, (features, mol) in models.items():
            self.assertIsInstance(features, list)
            self.assertIsInstance(mol, Chem.Mol)
            self.assertGreater(mol.GetNumAtoms(), 0)
    
    def test_generate_consensus_models_comprehensive(self):
        """Test comprehensive model set generation."""
        pharm = Pharmacophore()
        
        models = pharm.generate_consensus_models(
            mols=self.mols,
            model_set='comprehensive'
        )
        
        # Should have 5 models
        self.assertEqual(len(models), 5)
        self.assertIn('very_strict', models)
        self.assertIn('very_relaxed', models)
    
    def test_generate_consensus_models_without_mols(self):
        """Test generating only features without RDKit Mol."""
        pharm = Pharmacophore()
        
        models = pharm.generate_consensus_models(
            mols=self.mols,
            return_mols=False
        )
        
        # Should return only feature lists
        for model_name, features in models.items():
            self.assertIsInstance(features, list)
    
    def test_consensus_pharm_deprecation_warning(self):
        """Test that old method shows deprecation warning."""
        pharm = Pharmacophore()
        
        with self.assertWarns(DeprecationWarning):
            features = pharm.consensus_pharm(self.mols)
        
        # Should still work
        self.assertIsInstance(features, list)
    
    def test_determinism_across_runs(self):
        """Test that results are deterministic across runs."""
        pharm = Pharmacophore()
        
        models1 = pharm.generate_consensus_models(
            mols=self.mols,
            tolerance=2.0,
            occurrence_threshold=0.5
        )
        
        models2 = pharm.generate_consensus_models(
            mols=self.mols,
            tolerance=2.0,
            occurrence_threshold=0.5
        )
        
        # Should have same models
        self.assertEqual(set(models1.keys()), set(models2.keys()))
        
        # Compare feature counts
        for model_name in models1.keys():
            features1, mol1 = models1[model_name]
            features2, mol2 = models2[model_name]
            
            self.assertEqual(len(features1), len(features2),
                           f"Different feature counts in {model_name}")
            self.assertEqual(mol1.GetNumAtoms(), mol2.GetNumAtoms(),
                           f"Different atom counts in {model_name}")


class TestBackwardCompatibility(unittest.TestCase):
    """Test backward compatibility with old API."""
    
    def setUp(self):
        """Set up test molecules."""
        self.mol1 = Chem.MolFromSmiles('N')
        self.mol2 = Chem.MolFromSmiles('N')
        
        AllChem.EmbedMolecule(self.mol1, randomSeed=42)
        AllChem.EmbedMolecule(self.mol2, randomSeed=43)
        
        self.mols = [self.mol1, self.mol2]
    
    def test_old_method_still_works(self):
        """Test that old consensus_pharm method still works."""
        pharm = Pharmacophore()
        
        with self.assertWarns(DeprecationWarning):
            features = pharm.consensus_pharm(
                mols=self.mols,
                distance_threshold=2.0
            )
        
        self.assertIsInstance(features, list)
        
        # Check feature format
        for feature in features:
            self.assertEqual(len(feature), 5)
            self.assertIsInstance(feature[0], str)  # Type
            self.assertEqual(feature[1], ())  # Empty tuple for atom_indices


if __name__ == '__main__':
    unittest.main()

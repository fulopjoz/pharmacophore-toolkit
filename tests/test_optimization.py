"""Unit tests for consensus pharmacophore optimization.

This test suite covers:
- ConsensusOptimizer class
- optimize_consensus convenience function
- Enrichment factor calculation
- Parameter grid search
"""

import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.optimization import ConsensusOptimizer, optimize_consensus
from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol


class TestConsensusOptimizer(unittest.TestCase):
    """Test ConsensusOptimizer class."""

    def setUp(self):
        """Set up test fixtures with simple molecules."""
        self.optimizer = ConsensusOptimizer(
            n_conformers=3,
            random_state=42
        )

        # Create simple reference molecules with 3D conformers
        # These are simple molecules that can generate pharmacophore features
        smiles_list = [
            'CCO',   # Ethanol - has Donor/Acceptor
            'CCN',   # Ethylamine - has Donor
            'c1ccccc1O',  # Phenol - has Aromatic, Donor, Acceptor
        ]

        self.reference_mols = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            self.reference_mols.append(mol)

        # Simple actives (similar to references)
        active_smiles = ['CCO', 'CCCO', 'c1ccccc1N']
        self.actives = []
        for smi in active_smiles:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            self.actives.append(mol)

        # Simple decoys (different structure)
        decoy_smiles = ['CCCCCCCC', 'C1CCCCC1', 'CCCCCC']
        self.decoys = []
        for smi in decoy_smiles:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            self.decoys.append(mol)

    def test_initialization(self):
        """Test optimizer initialization."""
        optimizer = ConsensusOptimizer(
            n_conformers=10,
            linkage='complete',
            random_state=123
        )

        self.assertEqual(optimizer.n_conformers, 10)
        self.assertEqual(optimizer.linkage, 'complete')
        self.assertEqual(optimizer.random_state, 123)

    def test_enrichment_factor_perfect(self):
        """Test enrichment factor with perfect separation."""
        # All actives at top
        y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02]

        # At 30% (3 compounds), should have 100% actives
        # EF = (3/3) / (3/10) = 1.0 / 0.3 = 3.33
        ef = ConsensusOptimizer._enrichment_factor(y_true, y_scores, percentage=0.30)
        self.assertGreater(ef, 3.0)

    def test_enrichment_factor_random(self):
        """Test enrichment factor with shuffled order."""
        # Actives and decoys with similar scores (shuffled)
        y_true = [1, 0, 0, 1, 0, 0, 1, 0, 0, 0]
        y_scores = [0.55, 0.52, 0.48, 0.51, 0.49, 0.47, 0.53, 0.46, 0.45, 0.44]

        # EF at 50% should be finite and calculable
        ef = ConsensusOptimizer._enrichment_factor(y_true, y_scores, percentage=0.50)
        # Just verify it returns a valid number
        self.assertIsInstance(ef, (int, float, np.floating))

    def test_enrichment_factor_empty(self):
        """Test enrichment factor with empty lists."""
        ef = ConsensusOptimizer._enrichment_factor([], [], percentage=0.01)
        self.assertEqual(ef, 0.0)

    def test_enrichment_factor_no_actives(self):
        """Test enrichment factor with no actives."""
        y_true = [0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.6, 0.5]

        ef = ConsensusOptimizer._enrichment_factor(y_true, y_scores, percentage=0.20)
        self.assertEqual(ef, 0.0)

    def test_score_molecule_returns_tuple(self):
        """Test that score_molecule returns (shape, color, combo) tuple."""
        # Create a simple pharmacophore mol
        features = [
            ['Donor', (), 0.0, 0.0, 0.0],
            ['Acceptor', (), 2.0, 0.0, 0.0]
        ]
        pharm_mol = PharmacophoreToMol.convert(features, enable_color_features=True)

        # Score a molecule
        mol = self.actives[0]
        shape, color, combo = self.optimizer.score_molecule(mol, pharm_mol)

        self.assertIsInstance(shape, float)
        self.assertIsInstance(color, float)
        self.assertIsInstance(combo, float)
        self.assertGreaterEqual(shape, 0.0)
        self.assertGreaterEqual(color, 0.0)
        self.assertAlmostEqual(combo, shape + color, places=5)

    def test_evaluate_parameters_returns_metrics(self):
        """Test that evaluate_parameters returns expected metrics."""
        results = self.optimizer.evaluate_parameters(
            tolerance=2.0,
            occurrence_threshold=0.5,
            reference_mols=self.reference_mols,
            actives=self.actives,
            decoys=self.decoys
        )

        # Check all expected keys are present
        expected_keys = ['roc_auc', 'ef_1', 'ef_5', 'tolerance',
                        'occurrence_threshold', 'n_features']
        for key in expected_keys:
            self.assertIn(key, results)

        # Check parameter values are preserved
        self.assertEqual(results['tolerance'], 2.0)
        self.assertEqual(results['occurrence_threshold'], 0.5)

        # Check metrics are in valid ranges
        self.assertGreaterEqual(results['roc_auc'], 0.0)
        self.assertLessEqual(results['roc_auc'], 1.0)
        self.assertGreaterEqual(results['ef_1'], 0.0)
        self.assertGreaterEqual(results['ef_5'], 0.0)
        self.assertGreaterEqual(results['n_features'], 0)

    def test_evaluate_parameters_invalid_returns_empty(self):
        """Test that invalid parameters return empty results."""
        # Very strict parameters that might yield no features
        results = self.optimizer.evaluate_parameters(
            tolerance=0.01,  # Very small
            occurrence_threshold=1.0,  # Must be in all molecules
            reference_mols=self.reference_mols,
            actives=self.actives,
            decoys=self.decoys
        )

        # Should return default values
        self.assertEqual(results['roc_auc'], 0.5)
        self.assertEqual(results['tolerance'], 0.01)
        self.assertEqual(results['occurrence_threshold'], 1.0)

    def test_optimize_returns_expected_structure(self):
        """Test that optimize returns expected result structure."""
        results = self.optimizer.optimize(
            reference_mols=self.reference_mols,
            actives=self.actives,
            decoys=self.decoys,
            tolerance_range=(1.0, 2.0),
            occurrence_range=(0.3, 0.5),
            n_points=2,
            verbose=False
        )

        # Check structure
        self.assertIn('best_params', results)
        self.assertIn('best_score', results)
        self.assertIn('best_results', results)
        self.assertIn('all_results', results)
        self.assertIn('grid', results)

        # Check best_params structure
        self.assertIn('tolerance', results['best_params'])
        self.assertIn('occurrence_threshold', results['best_params'])

        # Check grid structure
        self.assertIn('tolerances', results['grid'])
        self.assertIn('occurrences', results['grid'])

        # Check all_results count (2x2 grid)
        self.assertEqual(len(results['all_results']), 4)

    def test_optimize_finds_reasonable_params(self):
        """Test that optimization finds parameters in valid range."""
        results = self.optimizer.optimize(
            reference_mols=self.reference_mols,
            actives=self.actives,
            decoys=self.decoys,
            tolerance_range=(0.5, 3.0),
            occurrence_range=(0.3, 0.7),
            n_points=3,
            verbose=False
        )

        best = results['best_params']

        # Parameters should be within specified ranges
        self.assertGreaterEqual(best['tolerance'], 0.5)
        self.assertLessEqual(best['tolerance'], 3.0)
        self.assertGreaterEqual(best['occurrence_threshold'], 0.3)
        self.assertLessEqual(best['occurrence_threshold'], 0.7)

    def test_conformer_caching(self):
        """Test that conformers are cached correctly."""
        # Create a simple pharmacophore
        features = [['Donor', (), 0.0, 0.0, 0.0]]
        pharm_mol = PharmacophoreToMol.convert(features, enable_color_features=True)

        # Score same molecule twice
        mol = self.actives[0]
        smiles = Chem.MolToSmiles(mol)

        # First call
        self.optimizer.score_molecule(mol, pharm_mol)
        self.assertIn(smiles, self.optimizer._conformer_cache)

        # Check cache has conformers
        cached = self.optimizer._conformer_cache[smiles]
        self.assertGreater(len(cached), 0)

    def test_clear_cache(self):
        """Test that clear_cache empties the conformer cache."""
        # Generate some conformers
        features = [['Donor', (), 0.0, 0.0, 0.0]]
        pharm_mol = PharmacophoreToMol.convert(features, enable_color_features=True)
        self.optimizer.score_molecule(self.actives[0], pharm_mol)

        # Cache should have entries
        self.assertGreater(len(self.optimizer._conformer_cache), 0)

        # Clear cache
        self.optimizer.clear_cache()

        # Cache should be empty
        self.assertEqual(len(self.optimizer._conformer_cache), 0)


class TestOptimizeConsensusFunction(unittest.TestCase):
    """Test optimize_consensus convenience function."""

    def setUp(self):
        """Set up minimal test molecules."""
        # Minimal setup for quick tests
        smiles = ['CCO', 'CCN']

        self.refs = []
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            self.refs.append(mol)

        # One active, one decoy
        active = Chem.MolFromSmiles('CCCO')
        active = Chem.AddHs(active)
        AllChem.EmbedMolecule(active, randomSeed=42)
        self.actives = [active]

        decoy = Chem.MolFromSmiles('CCCCCC')
        decoy = Chem.AddHs(decoy)
        AllChem.EmbedMolecule(decoy, randomSeed=42)
        self.decoys = [decoy]

    def test_function_works(self):
        """Test that convenience function executes without error."""
        results = optimize_consensus(
            reference_mols=self.refs,
            actives=self.actives,
            decoys=self.decoys,
            n_conformers=2,
            n_points=2,
            verbose=False
        )

        self.assertIn('best_params', results)
        self.assertIn('best_score', results)

    def test_function_accepts_all_parameters(self):
        """Test that function accepts all parameters."""
        results = optimize_consensus(
            reference_mols=self.refs,
            actives=self.actives,
            decoys=self.decoys,
            n_conformers=2,
            tolerance_range=(1.0, 2.0),
            occurrence_range=(0.4, 0.6),
            n_points=2,
            metric='ef_1',
            verbose=False
        )

        self.assertIsNotNone(results)


class TestEnrichmentFactorEdgeCases(unittest.TestCase):
    """Test edge cases for enrichment factor calculation."""

    def test_single_compound(self):
        """Test EF with single compound."""
        ef = ConsensusOptimizer._enrichment_factor([1], [1.0], percentage=1.0)
        self.assertGreater(ef, 0.0)

    def test_all_actives(self):
        """Test EF when all compounds are active."""
        y_true = [1, 1, 1, 1, 1]
        y_scores = [0.9, 0.8, 0.7, 0.6, 0.5]

        ef = ConsensusOptimizer._enrichment_factor(y_true, y_scores, percentage=0.20)
        # When all are actives, EF at any percentage should be 1.0
        self.assertAlmostEqual(ef, 1.0, places=2)

    def test_all_decoys(self):
        """Test EF when all compounds are decoys."""
        y_true = [0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.6, 0.5]

        ef = ConsensusOptimizer._enrichment_factor(y_true, y_scores, percentage=0.20)
        self.assertEqual(ef, 0.0)

    def test_very_small_percentage(self):
        """Test EF with very small percentage (should take at least 1)."""
        y_true = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

        # 0.1% of 10 compounds rounds up to 1
        ef = ConsensusOptimizer._enrichment_factor(y_true, y_scores, percentage=0.001)
        # Active is at top, so EF should be high
        self.assertGreater(ef, 1.0)


if __name__ == '__main__':
    unittest.main()

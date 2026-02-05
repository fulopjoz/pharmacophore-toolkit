"""Tests for DOE-based pharmacophore optimization."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Skip all tests if scikit-optimize not installed
pytest.importorskip("skopt")

from pharmacophore.doe_optimization import PharmacophoreOptimizer, optimize_pharmacophore


def create_test_molecules(n_mols: int, seed: int = 42) -> list:
    """Create simple test molecules with conformers."""
    np.random.seed(seed)
    smiles_list = [
        'CCO', 'CCCO', 'CC(C)O', 'CCCCO', 'CC(C)(C)O',
        'CCN', 'CCCN', 'CC(C)N', 'CCCCN', 'CC(C)(C)N',
        'c1ccccc1O', 'c1ccccc1N', 'c1ccc(O)cc1', 'c1ccc(N)cc1',
    ]

    mols = []
    for i in range(n_mols):
        smi = smiles_list[i % len(smiles_list)]
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=seed + i)
            if mol.GetNumConformers() > 0:
                mols.append(mol)

    return mols


def create_reference_mols(n_refs: int = 3) -> list:
    """Create reference molecules with conformers."""
    smiles_list = [
        'c1ccc(O)c(N)c1',  # Phenol with amine
        'c1ccc(N)c(O)c1',  # Aniline with hydroxyl
        'Oc1ccc(N)cc1',    # Para-aminophenol
    ]

    refs = []
    for i, smi in enumerate(smiles_list[:n_refs]):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42 + i)
            if mol.GetNumConformers() > 0:
                refs.append(mol)

    return refs


class TestPharmacophoreOptimizer:
    """Test suite for PharmacophoreOptimizer class."""

    @pytest.fixture
    def small_dataset(self):
        """Create small dataset for fast testing."""
        refs = create_reference_mols(3)
        actives = create_test_molecules(5, seed=100)
        decoys = create_test_molecules(10, seed=200)
        return refs, actives, decoys

    def test_init(self, small_dataset):
        """Test optimizer initialization."""
        refs, actives, decoys = small_dataset
        optimizer = PharmacophoreOptimizer(refs, actives, decoys, n_conformers=1)

        assert len(optimizer.reference_mols) == 3
        assert len(optimizer.actives) <= 5
        assert len(optimizer.decoys) <= 10
        assert len(optimizer.history) == 0

    def test_evaluate_single(self, small_dataset):
        """Test single parameter evaluation."""
        refs, actives, decoys = small_dataset
        optimizer = PharmacophoreOptimizer(refs, actives, decoys, n_conformers=1)

        auc = optimizer.evaluate(tolerance=2.0, occurrence=0.3, shape_weight=0.5)

        assert 0.0 <= auc <= 1.0
        assert len(optimizer.history) == 1
        assert 'auc' in optimizer.history[0]
        assert 'tolerance' in optimizer.history[0]

    def test_optimize_small(self, small_dataset):
        """Test optimization with small budget."""
        refs, actives, decoys = small_dataset
        optimizer = PharmacophoreOptimizer(refs, actives, decoys, n_conformers=1)

        result = optimizer.optimize(
            n_calls=5,
            n_random_starts=3,
            verbose=False
        )

        assert 'best_params' in result
        assert 'best_auc' in result
        assert 'history' in result
        assert result['best_params']['tolerance'] >= 0.5
        assert result['best_params']['tolerance'] <= 4.0
        assert 0.0 <= result['best_auc'] <= 1.0

    def test_parameter_ranges(self, small_dataset):
        """Test that optimization respects parameter ranges."""
        refs, actives, decoys = small_dataset
        optimizer = PharmacophoreOptimizer(refs, actives, decoys, n_conformers=1)

        result = optimizer.optimize(
            n_calls=5,
            n_random_starts=5,
            tolerance_range=(1.0, 2.0),
            occurrence_range=(0.2, 0.4),
            shape_weight_range=(0.4, 0.6),
            verbose=False
        )

        params = result['best_params']
        assert 1.0 <= params['tolerance'] <= 2.0
        assert 0.2 <= params['occurrence'] <= 0.4
        assert 0.4 <= params['shape_weight'] <= 0.6

    def test_history_accumulates(self, small_dataset):
        """Test that history accumulates across evaluations."""
        refs, actives, decoys = small_dataset
        optimizer = PharmacophoreOptimizer(refs, actives, decoys, n_conformers=1)

        optimizer.evaluate(1.0, 0.3, 0.5)
        optimizer.evaluate(2.0, 0.5, 0.6)
        optimizer.evaluate(3.0, 0.7, 0.4)

        assert len(optimizer.history) == 3

    def test_repr(self, small_dataset):
        """Test string representation."""
        refs, actives, decoys = small_dataset
        optimizer = PharmacophoreOptimizer(refs, actives, decoys, n_conformers=1)

        repr_str = repr(optimizer)
        assert 'PharmacophoreOptimizer' in repr_str
        assert 'n_refs=' in repr_str

    def test_parameter_importance(self, small_dataset):
        """Test parameter importance calculation."""
        refs, actives, decoys = small_dataset
        optimizer = PharmacophoreOptimizer(refs, actives, decoys, n_conformers=1)

        # Need at least 10 evaluations
        optimizer.optimize(n_calls=12, n_random_starts=10, verbose=False)

        importance = optimizer.get_parameter_importance()

        assert 'tolerance' in importance
        assert 'occurrence' in importance
        assert 'shape_weight' in importance
        assert sum(importance.values()) == pytest.approx(1.0, abs=0.01)


class TestConvenienceFunction:
    """Test the convenience function."""

    def test_optimize_pharmacophore(self):
        """Test optimize_pharmacophore function."""
        refs = create_reference_mols(2)
        actives = create_test_molecules(3, seed=100)
        decoys = create_test_molecules(5, seed=200)

        result = optimize_pharmacophore(
            refs, actives, decoys,
            n_calls=10,
            n_conformers=1,
            verbose=False
        )

        assert 'best_params' in result
        assert 'best_auc' in result


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_features(self):
        """Test handling when no consensus features generated."""
        refs = create_reference_mols(2)
        actives = create_test_molecules(3, seed=100)
        decoys = create_test_molecules(5, seed=200)

        optimizer = PharmacophoreOptimizer(refs, actives, decoys, n_conformers=1)

        # Very strict parameters should generate few/no features
        auc = optimizer.evaluate(tolerance=0.1, occurrence=1.0, shape_weight=0.5)

        # Should return random performance (0.5) if no features
        assert 0.0 <= auc <= 1.0

    def test_reproducibility(self):
        """Test that results are reproducible with same seed."""
        refs = create_reference_mols(2)
        actives = create_test_molecules(3, seed=100)
        decoys = create_test_molecules(5, seed=200)

        optimizer1 = PharmacophoreOptimizer(
            refs, actives, decoys, n_conformers=1, random_state=42
        )
        result1 = optimizer1.optimize(n_calls=10, n_random_starts=10, verbose=False)

        optimizer2 = PharmacophoreOptimizer(
            refs, actives, decoys, n_conformers=1, random_state=42
        )
        result2 = optimizer2.optimize(n_calls=10, n_random_starts=10, verbose=False)

        assert result1['best_auc'] == pytest.approx(result2['best_auc'], abs=0.01)

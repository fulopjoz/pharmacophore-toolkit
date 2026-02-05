"""Tests for AutoPharmacophoreOptimizer module."""

import os
import pytest
import tempfile
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.auto_optimizer import (
    AutoPharmacophoreOptimizer,
    auto_optimize_pharmacophore,
)


# Fixtures for test data
@pytest.fixture
def simple_mols():
    """Create simple test molecules."""
    smiles_list = ['CCO', 'CCCO', 'CCCCO', 'CC(C)O', 'CCC(C)O']
    mols = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        mols.append(mol)
    return mols


@pytest.fixture
def reference_mols(simple_mols):
    """Reference molecules with 3D conformers."""
    return simple_mols[:3]


@pytest.fixture
def active_mols(simple_mols):
    """Active molecules."""
    return simple_mols[:2]


@pytest.fixture
def decoy_mols(simple_mols):
    """Decoy molecules."""
    return simple_mols[2:]


@pytest.fixture
def temp_csv_files():
    """Create temporary CSV files for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create actives CSV
        actives_path = Path(tmpdir) / 'actives.csv'
        with open(actives_path, 'w') as f:
            f.write('smiles,name\n')
            f.write('CCO,ethanol\n')
            f.write('CCCO,propanol\n')

        # Create decoys CSV
        decoys_path = Path(tmpdir) / 'decoys.csv'
        with open(decoys_path, 'w') as f:
            f.write('smiles,name\n')
            f.write('CCCCO,butanol\n')
            f.write('CC(C)O,isopropanol\n')
            f.write('CCC(C)O,2-butanol\n')

        # Create reference SDF
        refs_path = Path(tmpdir) / 'refs.sdf'
        writer = Chem.SDWriter(str(refs_path))
        for smi in ['CCO', 'CCCO', 'CCCCO']:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            writer.write(mol)
        writer.close()

        yield {
            'actives': str(actives_path),
            'decoys': str(decoys_path),
            'refs': str(refs_path),
            'tmpdir': tmpdir
        }


class TestAutoPharmacophoreOptimizerInit:
    """Tests for optimizer initialization."""

    def test_init_default_params(self):
        """Test initialization with default parameters."""
        optimizer = AutoPharmacophoreOptimizer()
        assert optimizer.n_conformers == 5
        assert optimizer.random_state == 42
        assert optimizer.cache_enabled is True

    def test_init_custom_params(self):
        """Test initialization with custom parameters."""
        optimizer = AutoPharmacophoreOptimizer(
            n_conformers=10,
            random_state=123,
            cache_enabled=False
        )
        assert optimizer.n_conformers == 10
        assert optimizer.random_state == 123
        assert optimizer.cache_enabled is False

    def test_cache_disabled(self):
        """Test cache can be disabled."""
        optimizer = AutoPharmacophoreOptimizer(cache_enabled=False)
        assert optimizer.cache is None


class TestDataLoading:
    """Tests for data loading functionality."""

    def test_load_from_objects(self, reference_mols, active_mols, decoy_mols):
        """Test loading molecules from Python objects."""
        optimizer = AutoPharmacophoreOptimizer()
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        assert len(optimizer.reference_mols) == 3
        assert len(optimizer.actives) == 2
        assert len(optimizer.decoys) == 3

    def test_load_from_objects_returns_self(self, reference_mols, active_mols, decoy_mols):
        """Test load_from_objects returns self for chaining."""
        optimizer = AutoPharmacophoreOptimizer()
        result = optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)
        assert result is optimizer

    def test_load_from_objects_empty_refs_raises(self, active_mols, decoy_mols):
        """Test empty reference list raises ValueError."""
        optimizer = AutoPharmacophoreOptimizer()
        with pytest.raises(ValueError, match="reference_mols cannot be empty"):
            optimizer.load_from_objects([], active_mols, decoy_mols)

    def test_load_from_objects_empty_actives_raises(self, reference_mols, decoy_mols):
        """Test empty actives list raises ValueError."""
        optimizer = AutoPharmacophoreOptimizer()
        with pytest.raises(ValueError, match="actives cannot be empty"):
            optimizer.load_from_objects(reference_mols, [], decoy_mols)

    def test_load_from_files(self, temp_csv_files):
        """Test loading from files."""
        optimizer = AutoPharmacophoreOptimizer()
        optimizer.load_from_files(
            reference_file=temp_csv_files['refs'],
            actives_file=temp_csv_files['actives'],
            decoys_file=temp_csv_files['decoys']
        )

        assert len(optimizer.reference_mols) == 3
        assert len(optimizer.actives) == 2
        assert len(optimizer.decoys) == 3

    def test_load_from_files_missing_file(self, temp_csv_files):
        """Test missing file raises FileNotFoundError."""
        optimizer = AutoPharmacophoreOptimizer()
        with pytest.raises(FileNotFoundError):
            optimizer.load_from_files(
                reference_file='nonexistent.sdf',
                actives_file=temp_csv_files['actives'],
                decoys_file=temp_csv_files['decoys']
            )


class TestConformerGeneration:
    """Tests for conformer generation."""

    def test_conformer_caching(self, reference_mols, active_mols, decoy_mols):
        """Test conformers are cached by SMILES."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=3)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        # Access conformers twice for same molecule
        mol = active_mols[0]
        conf1 = optimizer._get_conformers(mol)
        conf2 = optimizer._get_conformers(mol)

        # Should be same cached object
        assert conf1 is conf2


class TestEvaluation:
    """Tests for parameter evaluation."""

    def test_evaluate_returns_auc(self, reference_mols, active_mols, decoy_mols):
        """Test evaluate returns a valid AUC."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        auc = optimizer.evaluate(tolerance=2.0, occurrence=0.5)

        assert 0.0 <= auc <= 1.0

    def test_evaluate_stores_history(self, reference_mols, active_mols, decoy_mols):
        """Test evaluate stores result in history."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        optimizer.evaluate(tolerance=2.0, occurrence=0.5)

        assert len(optimizer.history) == 1
        assert 'tolerance' in optimizer.history[0]
        assert 'occurrence' in optimizer.history[0]
        assert 'roc_auc' in optimizer.history[0]

    def test_evaluate_no_data_raises(self):
        """Test evaluate without data raises ValueError."""
        optimizer = AutoPharmacophoreOptimizer()
        with pytest.raises(ValueError, match="No data loaded"):
            optimizer.evaluate(2.0, 0.5)


class TestOptimization:
    """Tests for Bayesian optimization."""

    def test_optimize_minimal_budget(self, reference_mols, active_mols, decoy_mols):
        """Test optimization with minimal budget."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        result = optimizer.optimize(
            n_calls=5,
            n_random_starts=3,
            scoring_method='pharm2d',
            verbose=False
        )

        assert 'best_params' in result
        assert 'best_auc' in result
        assert 'tolerance' in result['best_params']
        assert 'occurrence' in result['best_params']
        assert 0.0 <= result['best_auc'] <= 1.0

    @pytest.mark.filterwarnings("ignore::DeprecationWarning")
    def test_optimize_shape_method(self, reference_mols, active_mols, decoy_mols):
        """Test optimization with shape scoring method (deprecated path)."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        result = optimizer.optimize(
            n_calls=5,
            n_random_starts=3,
            scoring_method='shape',
            verbose=False
        )

        assert 'shape_weight' in result['best_params']

    def test_optimize_stores_best_result(self, reference_mols, active_mols, decoy_mols):
        """Test optimization stores best_result."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        optimizer.optimize(n_calls=5, n_random_starts=3, verbose=False)

        assert optimizer.best_result is not None


class TestExport:
    """Tests for model export."""

    def test_export_json(self, reference_mols, active_mols, decoy_mols):
        """Test JSON export."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)
        optimizer.optimize(n_calls=5, n_random_starts=3, verbose=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            outputs = optimizer.export_model(tmpdir, include_pml=False)

            assert 'json' in outputs
            assert os.path.exists(outputs['json'])

            with open(outputs['json']) as f:
                data = json.load(f)
                assert 'best_params' in data
                assert 'best_auc' in data
                assert 'features' in data

    def test_export_pml(self, reference_mols, active_mols, decoy_mols):
        """Test PyMOL script export."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)
        optimizer.optimize(n_calls=5, n_random_starts=3, verbose=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            outputs = optimizer.export_model(tmpdir, include_json=False)

            assert 'pml' in outputs
            assert os.path.exists(outputs['pml'])

            with open(outputs['pml']) as f:
                content = f.read()
                assert 'pseudoatom' in content

    def test_export_no_optimization_raises(self):
        """Test export without optimization raises ValueError."""
        optimizer = AutoPharmacophoreOptimizer()
        with pytest.raises(ValueError, match="No optimization results"):
            optimizer.export_model('/tmp')


class TestComparison:
    """Tests for method comparison."""

    def test_compare_methods(self, reference_mols, active_mols, decoy_mols):
        """Test compare_methods returns DataFrame."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)
        optimizer.optimize(n_calls=5, n_random_starts=3, verbose=False)

        comparison = optimizer.compare_methods()

        assert len(comparison) == 3  # Pharm2D, Pharm3D, Shape
        assert 'method' in comparison.columns
        assert 'roc_auc' in comparison.columns


class TestConvenienceFunction:
    """Tests for auto_optimize_pharmacophore function."""

    def test_convenience_function(self, temp_csv_files):
        """Test convenience function works end-to-end."""
        result = auto_optimize_pharmacophore(
            reference_file=temp_csv_files['refs'],
            actives_file=temp_csv_files['actives'],
            decoys_file=temp_csv_files['decoys'],
            n_calls=10,  # Must be >= n_random_starts (default 10)
            scoring_method='pharm2d',
            verbose=False
        )

        assert 'best_params' in result
        assert 'best_auc' in result

    def test_convenience_function_with_export(self, temp_csv_files):
        """Test convenience function with export."""
        with tempfile.TemporaryDirectory() as output_dir:
            result = auto_optimize_pharmacophore(
                reference_file=temp_csv_files['refs'],
                actives_file=temp_csv_files['actives'],
                decoys_file=temp_csv_files['decoys'],
                n_calls=10,  # Must be >= n_random_starts (default 10)
                output_dir=output_dir,
                verbose=False
            )

            # Check output files created
            files = os.listdir(output_dir)
            assert any('json' in f for f in files)


class TestCaching:
    """Tests for caching behavior."""

    def test_cache_enabled(self, reference_mols, active_mols, decoy_mols):
        """Test caching improves performance."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2, cache_enabled=True)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        # Run two evaluations with same parameters
        optimizer.evaluate(2.0, 0.5)
        optimizer.evaluate(2.0, 0.5)

        stats = optimizer.cache.stats()
        assert stats['consensus_hits'] >= 1

    def test_clear_cache(self, reference_mols, active_mols, decoy_mols):
        """Test cache clearing."""
        optimizer = AutoPharmacophoreOptimizer(n_conformers=2)
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        optimizer.evaluate(2.0, 0.5)
        optimizer.clear_cache()

        stats = optimizer.cache.stats()
        assert stats['consensus_size'] == 0


class TestRepr:
    """Tests for string representation."""

    def test_repr_empty(self):
        """Test repr with no data."""
        optimizer = AutoPharmacophoreOptimizer()
        repr_str = repr(optimizer)
        assert 'AutoPharmacophoreOptimizer' in repr_str
        assert 'refs=0' in repr_str

    def test_repr_with_data(self, reference_mols, active_mols, decoy_mols):
        """Test repr with loaded data."""
        optimizer = AutoPharmacophoreOptimizer()
        optimizer.load_from_objects(reference_mols, active_mols, decoy_mols)

        repr_str = repr(optimizer)
        assert 'refs=3' in repr_str
        assert 'actives=2' in repr_str

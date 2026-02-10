"""Tests for Optuna-based multi-objective optimization."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    from pharmacophore.optuna_optimizer import OptunaPharmacophoreOptimizer
    HAS_OPTUNA = True
except ImportError:
    HAS_OPTUNA = False


@pytest.fixture
def sample_molecules():
    """Create sample molecules for testing."""
    smiles = [
        'CC(=O)O',  # Acetic acid
        'CCO',      # Ethanol
        'c1ccccc1', # Benzene
        'CC(C)O',   # Isopropanol
        'CCCC',     # Butane
    ]

    mols = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol_h = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol_h, randomSeed=42)
            mols.append(mol_h)

    return mols


@pytest.mark.skipif(not HAS_OPTUNA, reason="Optuna not installed")
class TestOptunaPharmacophoreOptimizer:
    """Tests for OptunaPharmacophoreOptimizer."""

    def test_initialization(self, sample_molecules):
        """Test optimizer initialization."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        assert optimizer.evaluator is not None
        assert optimizer.random_state == 42
        assert optimizer.study is None
        assert optimizer.start_time is None

    def test_repr(self, sample_molecules):
        """Test string representation."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        repr_str = repr(optimizer)
        assert 'OptunaPharmacophoreOptimizer' in repr_str
        assert 'n_refs=2' in repr_str
        assert 'n_actives=3' in repr_str
        assert 'n_decoys=3' in repr_str
        assert 'trials=0' in repr_str

    def test_gp_sampler_minimal(self, sample_molecules):
        """Test GP sampler with minimal trials."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        result = optimizer.optimize(
            sampler='gp',
            n_trials=10,
            verbose=False
        )

        # Check result structure
        assert 'sampler' in result
        assert 'n_trials' in result
        assert 'wall_time_sec' in result
        assert 'pareto_front' in result
        assert 'best_auc' in result
        assert 'best_bedroc' in result
        assert 'best_auc_params' in result
        assert 'best_bedroc_params' in result
        assert 'history' in result
        assert 'study' in result

        # Check values
        assert result['sampler'] == 'gp'
        assert result['n_trials'] == 10
        assert result['wall_time_sec'] > 0
        assert len(result['pareto_front']) > 0
        assert 0.0 <= result['best_auc'] <= 1.0
        assert 0.0 <= result['best_bedroc'] <= 1.0
        assert len(result['history']) == 10

        # Check parameter structure (reference mode: only alignment params)
        params = result['best_auc_params']
        assert 'opt_param' in params
        assert 'n_conformers' in params
        assert 'max_preiters' in params
        assert 'max_postiters' in params
        assert 'aggregation' in params

    def test_nsga2_sampler_minimal(self, sample_molecules):
        """Test NSGA-II sampler with minimal trials."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        result = optimizer.optimize(
            sampler='nsga2',
            n_trials=10,
            verbose=False
        )

        # Check result structure
        assert result['sampler'] == 'nsga2'
        assert result['n_trials'] == 10
        assert result['wall_time_sec'] > 0
        assert len(result['pareto_front']) > 0
        assert 0.0 <= result['best_auc'] <= 1.0
        assert 0.0 <= result['best_bedroc'] <= 1.0
        assert len(result['history']) == 10

    def test_invalid_sampler(self, sample_molecules):
        """Test invalid sampler raises error."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        with pytest.raises(ValueError, match="sampler must be"):
            optimizer.optimize(sampler='invalid', n_trials=10)

    def test_pareto_front_structure(self, sample_molecules):
        """Test Pareto front structure and properties."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        result = optimizer.optimize(
            sampler='gp',
            n_trials=10,
            verbose=False
        )

        pareto_front = result['pareto_front']

        # Check structure
        assert len(pareto_front) > 0

        for solution in pareto_front:
            assert 'params' in solution
            assert 'roc_auc' in solution
            assert 'bedroc' in solution
            assert 'n_features' in solution
            assert 'ef_1' in solution
            assert 'ef_5' in solution
            assert 'ef_10' in solution

            # Check objective values
            assert 0.0 <= solution['roc_auc'] <= 1.0
            assert 0.0 <= solution['bedroc'] <= 1.0

        # Verify Pareto dominance property
        # No solution should dominate another
        for i, sol1 in enumerate(pareto_front):
            for j, sol2 in enumerate(pareto_front):
                if i != j:
                    # sol1 doesn't strictly dominate sol2
                    dominates = (
                        sol1['roc_auc'] >= sol2['roc_auc'] and
                        sol1['bedroc'] >= sol2['bedroc'] and
                        (sol1['roc_auc'] > sol2['roc_auc'] or sol1['bedroc'] > sol2['bedroc'])
                    )
                    assert not dominates, "Found dominated solution in Pareto front"

    def test_get_pareto_front(self, sample_molecules):
        """Test get_pareto_front method."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        # Should raise before optimization
        with pytest.raises(ValueError, match="Run optimize"):
            optimizer.get_pareto_front()

        # Run optimization
        optimizer.optimize(sampler='gp', n_trials=10, verbose=False)

        # Now should work
        pareto = optimizer.get_pareto_front()
        assert len(pareto) > 0
        assert all('params' in sol for sol in pareto)
        assert all('roc_auc' in sol for sol in pareto)
        assert all('bedroc' in sol for sol in pareto)

    def test_get_parameter_importance(self, sample_molecules):
        """Test parameter importance estimation."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        # Should raise before optimization
        with pytest.raises(ValueError, match="Run optimize"):
            optimizer.get_parameter_importance()

        # Run optimization
        optimizer.optimize(sampler='gp', n_trials=15, verbose=False)

        # Now should work
        importance = optimizer.get_parameter_importance()

        assert isinstance(importance, dict)
        assert len(importance) > 0

        # Check parameter names (reference mode: only alignment params)
        expected_params = {'opt_param', 'n_conformers', 'max_preiters', 'max_postiters', 'aggregation'}
        assert set(importance.keys()).issubset(expected_params)

        # Check values are valid (0-1 range for normalized importance)
        for param, value in importance.items():
            assert 0.0 <= value <= 1.0

    def test_timing_information(self, sample_molecules):
        """Test timing information is tracked."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        result = optimizer.optimize(
            sampler='gp',
            n_trials=10,
            verbose=False
        )

        # Check wall time
        assert result['wall_time_sec'] > 0
        assert result['wall_time_sec'] < 600  # Should complete in < 10 min

        # Check per-trial timing
        for entry in result['history']:
            assert 'eval_time_sec' in entry
            assert entry['eval_time_sec'] > 0

    def test_reproducibility(self, sample_molecules):
        """Test reproducibility with same random_state."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        # Create two optimizers with same seed
        opt1 = OptunaPharmacophoreOptimizer(refs, actives, decoys, random_state=42)
        opt2 = OptunaPharmacophoreOptimizer(refs, actives, decoys, random_state=42)

        # Run both
        result1 = opt1.optimize(sampler='gp', n_trials=5, verbose=False)
        result2 = opt2.optimize(sampler='gp', n_trials=5, verbose=False)

        # Should produce identical results
        assert result1['best_auc'] == result2['best_auc']
        assert result1['best_bedroc'] == result2['best_bedroc']
        assert len(result1['pareto_front']) == len(result2['pareto_front'])

    def test_history_structure(self, sample_molecules):
        """Test optimization history structure."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        result = optimizer.optimize(
            sampler='gp',
            n_trials=10,
            verbose=False
        )

        history = result['history']

        assert len(history) == 10

        for entry in history:
            assert 'trial' in entry
            assert 'params' in entry
            assert 'roc_auc' in entry
            assert 'bedroc' in entry
            assert 'n_features' in entry
            assert 'eval_time_sec' in entry

            # Check types
            assert isinstance(entry['trial'], int)
            assert isinstance(entry['params'], dict)
            assert isinstance(entry['roc_auc'], (float, int))
            assert isinstance(entry['bedroc'], (float, int))
            assert isinstance(entry['n_features'], int)
            assert isinstance(entry['eval_time_sec'], (float, int))

    def test_verbose_output(self, sample_molecules, capsys):
        """Test verbose output is produced."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        result = optimizer.optimize(
            sampler='gp',
            n_trials=5,
            verbose=True
        )

        captured = capsys.readouterr()

        # Check for expected output
        assert 'Multi-Objective Pharmacophore Optimization' in captured.out
        assert 'Trial' in captured.out
        assert 'AUC=' in captured.out
        assert 'BEDROC=' in captured.out
        assert 'OPTIMIZATION COMPLETE' in captured.out

    def test_best_params_different_objectives(self, sample_molecules):
        """Test that best AUC and best BEDROC can have different params."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = OptunaPharmacophoreOptimizer(
            refs, actives, decoys, random_state=42
        )

        result = optimizer.optimize(
            sampler='nsga2',
            n_trials=20,  # More trials to explore parameter space
            verbose=False
        )

        # These MAY be the same, but check structure
        auc_params = result['best_auc_params']
        bedroc_params = result['best_bedroc_params']

        # Both should have core alignment parameters (reference mode)
        for params in [auc_params, bedroc_params]:
            assert 'opt_param' in params
            assert 'n_conformers' in params
            assert 'max_preiters' in params
            assert 'max_postiters' in params
            assert 'aggregation' in params

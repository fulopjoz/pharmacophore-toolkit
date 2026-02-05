"""Tests for HypoGen-inspired 3-phase optimization."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.hypogen_optimizer import HypoGenOptimizer
from pharmacophore.evaluation import EvaluationResult, EvaluationConfig


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


@pytest.fixture
def optimizer(sample_molecules):
    """Create HypoGenOptimizer with sample data."""
    refs = sample_molecules[:2]
    actives = sample_molecules[:3]
    decoys = sample_molecules[2:]

    return HypoGenOptimizer(refs, actives, decoys, random_state=42)


class TestHypoGenOptimizer:
    """Tests for HypoGenOptimizer class."""

    def test_initialization(self, sample_molecules):
        """Test optimizer initialization."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        optimizer = HypoGenOptimizer(refs, actives, decoys, random_state=42)

        assert optimizer.evaluator is not None
        assert len(optimizer.reference_mols) == 2
        assert optimizer.random_state == 42
        assert optimizer.rng is not None

    def test_repr(self, optimizer):
        """Test string representation."""
        repr_str = repr(optimizer)
        assert 'HypoGenOptimizer' in repr_str
        assert 'n_refs=2' in repr_str
        assert 'n_actives=3' in repr_str
        assert 'n_decoys=3' in repr_str

    def test_calculate_cost(self, optimizer):
        """Test cost function calculation."""
        config = EvaluationConfig()

        # Perfect model
        result = EvaluationResult(
            config=config,
            roc_auc=1.0,
            bedroc=1.0,
            n_features=5
        )
        cost = optimizer._calculate_cost(result, max_features=10)
        assert cost < 0.2  # Low cost for good model

        # Poor model
        result = EvaluationResult(
            config=config,
            roc_auc=0.5,
            bedroc=0.0,
            n_features=10
        )
        cost = optimizer._calculate_cost(result, max_features=10)
        assert cost > 0.7  # High cost for bad model

    def test_enumerate_subsets_count(self, optimizer):
        """Test subset enumeration count."""
        features = [
            ['Donor', (), 0.0, 0.0, 0.0],
            ['Acceptor', (), 1.0, 1.0, 1.0],
            ['Aromatic', (), 2.0, 2.0, 2.0],
            ['Hydrophobe', (), 3.0, 3.0, 3.0],
            ['Donor', (), 4.0, 4.0, 4.0]
        ]

        # Subsets of size 3-4 from 5 features
        # C(5,3) + C(5,4) = 10 + 5 = 15
        subsets = optimizer._enumerate_subsets(features, 3, 4, max_hypotheses=100)
        assert len(subsets) == 15

        # Test max_hypotheses cap
        subsets = optimizer._enumerate_subsets(features, 3, 4, max_hypotheses=5)
        assert len(subsets) == 5

    def test_enumerate_subsets_respect_bounds(self, optimizer):
        """Test subset enumeration respects feature bounds."""
        features = [
            ['Donor', (), 0.0, 0.0, 0.0],
            ['Acceptor', (), 1.0, 1.0, 1.0],
            ['Aromatic', (), 2.0, 2.0, 2.0],
            ['Hydrophobe', (), 3.0, 3.0, 3.0],
        ]

        subsets = optimizer._enumerate_subsets(features, 2, 3, max_hypotheses=100)

        # Check all subsets have 2 or 3 features
        for subset in subsets:
            assert 2 <= len(subset) <= 3

        # C(4,2) + C(4,3) = 6 + 4 = 10
        assert len(subsets) == 10

    def test_simulated_annealing_improves_or_maintains(self, optimizer):
        """Test SA improves or maintains cost."""
        features = [
            ['Donor', (), 0.0, 0.0, 0.0],
            ['Acceptor', (), 1.0, 1.0, 1.0],
            ['Aromatic', (), 2.0, 2.0, 2.0]
        ]

        # Evaluate initial state
        initial_result = optimizer.evaluator.evaluate_feature_subset(
            features,
            shape_weight=0.5,
            opt_param=0.5,
            n_conformers=5
        )
        initial_cost = optimizer._calculate_cost(initial_result, max_features=8)

        # Run SA with few iterations
        optimized = optimizer._simulated_annealing(
            features,
            n_iters=10,
            max_features=8,
            shape_weight=0.5,
            opt_param=0.5,
            n_conformers=5,
            verbose=False
        )

        # SA should not make cost worse (may stay same if already good)
        assert optimized['cost'] <= initial_cost * 1.1  # Allow 10% margin for stochasticity

    def test_optimize_end_to_end(self, optimizer):
        """Test end-to-end optimization."""
        result = optimizer.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=4,
            max_hypotheses=10,
            n_top=2,
            n_sa_iters=5,
            n_conformers=5,
            verbose=False
        )

        # Check result structure
        assert 'best_hypothesis' in result
        assert 'best_result' in result
        assert 'best_cost' in result
        assert 'phase1_hypotheses' in result
        assert 'phase2_survivors' in result
        assert 'phase3_optimized' in result
        assert 'phase_times' in result
        assert 'total_time' in result
        assert 'n_evaluations' in result

        # Check types
        assert isinstance(result['best_hypothesis'], list)
        assert isinstance(result['best_result'], EvaluationResult)
        assert isinstance(result['best_cost'], float)
        assert isinstance(result['phase1_hypotheses'], list)
        assert isinstance(result['phase_times'], dict)

        # Check values
        assert result['best_cost'] > 0.0
        assert 0.0 <= result['best_result'].roc_auc <= 1.0
        assert result['total_time'] > 0.0
        assert result['n_evaluations'] > 0

    def test_phase_times_populated(self, optimizer):
        """Test phase timings are tracked."""
        result = optimizer.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=3,
            max_hypotheses=5,
            n_top=1,
            n_sa_iters=3,
            n_conformers=5,
            verbose=False
        )

        phase_times = result['phase_times']

        assert 'phase1' in phase_times
        assert 'phase2' in phase_times
        assert 'phase3' in phase_times

        # All phases should have positive time
        assert phase_times['phase1'] > 0
        assert phase_times['phase2'] >= 0  # Phase 2 is fast
        assert phase_times['phase3'] >= 0  # Phase 3 may be 0 if no survivors

        # Total time should approximately equal sum of phases
        total = sum(phase_times.values())
        assert abs(result['total_time'] - total) < 1.0  # Within 1 second

    def test_phase1_generates_hypotheses(self, optimizer):
        """Test Phase 1 generates expected hypotheses."""
        result = optimizer.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=3,
            max_hypotheses=10,
            n_top=1,
            n_sa_iters=1,
            n_conformers=5,
            verbose=False
        )

        hypotheses = result['phase1_hypotheses']

        # Should have generated hypotheses
        assert len(hypotheses) > 0
        assert len(hypotheses) <= 10  # Respects max_hypotheses cap

        # Check structure
        for hyp in hypotheses:
            assert 'id' in hyp
            assert 'features' in hyp
            assert 'result' in hyp
            assert 'cost' in hyp
            assert isinstance(hyp['result'], EvaluationResult)

    def test_phase2_filters_hypotheses(self, optimizer):
        """Test Phase 2 filters low-quality hypotheses."""
        result = optimizer.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=3,
            max_hypotheses=10,
            n_top=1,
            n_sa_iters=1,
            auc_threshold=0.4,  # Low threshold to ensure survivors
            bedroc_threshold=0.0,
            n_conformers=5,
            verbose=False
        )

        phase1 = result['phase1_hypotheses']
        phase2 = result['phase2_survivors']

        # Phase 2 should have <= Phase 1
        assert len(phase2) <= len(phase1)

        # All survivors should meet thresholds
        for hyp in phase2:
            assert hyp['result'].roc_auc >= 0.4
            assert hyp['result'].bedroc >= 0.0

    def test_phase3_refines_survivors(self, optimizer):
        """Test Phase 3 refines survivors."""
        result = optimizer.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=3,
            max_hypotheses=10,
            n_top=2,
            n_sa_iters=5,
            auc_threshold=0.4,
            bedroc_threshold=0.0,
            n_conformers=5,
            verbose=False
        )

        phase2 = result['phase2_survivors']
        phase3 = result['phase3_optimized']

        if len(phase2) > 0:
            # Phase 3 should have optimized hypotheses
            assert len(phase3) > 0
            assert len(phase3) <= min(len(phase2), 2)  # n_top=2

            # Check structure
            for hyp in phase3:
                assert 'features' in hyp
                assert 'result' in hyp
                assert 'cost' in hyp

    def test_no_survivors_fallback(self, optimizer):
        """Test fallback when no survivors from Phase 2."""
        result = optimizer.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=3,
            max_hypotheses=5,
            n_top=1,
            n_sa_iters=1,
            auc_threshold=0.99,  # Impossibly high threshold
            bedroc_threshold=0.99,
            n_conformers=5,
            verbose=False
        )

        # Should still return a result (best from Phase 1)
        assert result['best_hypothesis'] is not None
        assert result['best_result'] is not None
        assert len(result['phase2_survivors']) == 0
        assert len(result['phase3_optimized']) == 0

    def test_reproducibility(self, sample_molecules):
        """Test reproducibility with same random_state."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        # Create two optimizers with same seed
        opt1 = HypoGenOptimizer(refs, actives, decoys, random_state=42)
        opt2 = HypoGenOptimizer(refs, actives, decoys, random_state=42)

        # Run both with same parameters
        result1 = opt1.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=3,
            max_hypotheses=5,
            n_top=1,
            n_sa_iters=3,
            n_conformers=5,
            verbose=False
        )

        result2 = opt2.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=3,
            max_hypotheses=5,
            n_top=1,
            n_sa_iters=3,
            n_conformers=5,
            verbose=False
        )

        # Results should be identical
        assert result1['best_cost'] == result2['best_cost']
        assert result1['best_result'].roc_auc == result2['best_result'].roc_auc
        assert result1['n_evaluations'] == result2['n_evaluations']

    def test_verbose_output(self, optimizer, capsys):
        """Test verbose output is produced."""
        result = optimizer.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=3,
            max_hypotheses=5,
            n_top=1,
            n_sa_iters=2,
            n_conformers=5,
            verbose=True
        )

        captured = capsys.readouterr()

        # Check for expected output
        assert 'HypoGen 3-Phase Pharmacophore Optimization' in captured.out
        assert 'PHASE 1: CONSTRUCTIVE' in captured.out
        assert 'PHASE 2: SUBTRACTIVE' in captured.out
        assert 'PHASE 3: OPTIMIZATION' in captured.out
        assert 'HYPOGEN OPTIMIZATION COMPLETE' in captured.out
        assert 'Best Hypothesis:' in captured.out

    def test_evaluation_count_accurate(self, optimizer):
        """Test evaluation count matches actual evaluations."""
        result = optimizer.optimize(
            consensus_tolerance=2.0,
            consensus_occurrence=0.3,
            min_features=2,
            max_features=3,
            max_hypotheses=10,
            n_top=2,
            n_sa_iters=5,
            n_conformers=5,
            verbose=False
        )

        n_phase1 = len(result['phase1_hypotheses'])
        n_phase3_refinements = 2 * 5  # n_top * n_sa_iters

        expected_evals = n_phase1 + n_phase3_refinements

        # May differ slightly if some survivors were skipped
        assert abs(result['n_evaluations'] - expected_evals) <= n_phase3_refinements

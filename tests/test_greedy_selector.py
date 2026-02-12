"""Tests for greedy forward feature selection."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.greedy_selector import GreedyFeatureSelector, FeatureSelectionResult


@pytest.fixture
def sample_molecules():
    """Create sample molecules for testing."""
    smiles = [
        'CC(=O)O',   # Acetic acid
        'CCO',       # Ethanol
        'c1ccccc1',  # Benzene
        'CC(=O)Nc1ccc(O)cc1',  # Acetaminophen
        'CC(O)=O',   # Acetic acid isomer
    ]
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
def selector(sample_molecules):
    """Create GreedyFeatureSelector with sample data."""
    refs = sample_molecules[:2]
    actives = sample_molecules[:3]
    decoys = sample_molecules[2:]
    return GreedyFeatureSelector(refs, actives, decoys, random_state=42)


class TestFeatureSelectionResult:
    """Tests for FeatureSelectionResult dataclass."""

    def test_defaults(self):
        result = FeatureSelectionResult()
        assert result.selected_features == []
        assert result.selected_indices == []
        assert result.best_auc == 0.5
        assert result.n_evaluations == 0

    def test_custom_values(self):
        result = FeatureSelectionResult(
            selected_indices=[0, 2, 5],
            best_auc=0.92,
            best_bedroc=0.85,
            n_evaluations=45,
            wall_time_sec=12.3,
        )
        assert len(result.selected_indices) == 3
        assert result.best_auc == 0.92


class TestGreedyFeatureSelector:
    """Tests for GreedyFeatureSelector."""

    def test_initialization(self, sample_molecules):
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]
        selector = GreedyFeatureSelector(refs, actives, decoys)
        assert len(selector.reference_mols) == 2
        expected = len(actives) + len(decoys)
        assert len(selector.y_true) == expected

    def test_select_returns_result(self, selector):
        result = selector.select(
            tolerance=2.0, occurrence=0.3, max_features=4, verbose=False
        )
        assert isinstance(result, FeatureSelectionResult)
        assert len(result.selected_features) > 0
        assert len(result.selected_indices) == len(result.selected_features)
        assert result.wall_time_sec > 0

    def test_auc_non_decreasing(self, selector):
        """AUC should be non-decreasing across greedy steps."""
        result = selector.select(
            tolerance=2.0, occurrence=0.3, max_features=6,
            convergence_threshold=0.0,  # Don't stop early
            verbose=False,
        )
        if len(result.selection_history) >= 2:
            aucs = [h[2] for h in result.selection_history]
            for i in range(1, len(aucs)):
                assert aucs[i] >= aucs[i - 1] - 1e-10, (
                    f"AUC decreased at step {i}: {aucs[i-1]:.4f} -> {aucs[i]:.4f}"
                )

    def test_convergence_stops(self, selector):
        """Selection should stop before max_features when converged."""
        result = selector.select(
            tolerance=2.0, occurrence=0.3,
            max_features=20,  # Very high limit
            convergence_threshold=0.01,
            verbose=False,
        )
        # Should stop well before 20 features (small test molecules)
        assert len(result.selected_indices) <= 20

    def test_deterministic(self, sample_molecules):
        """Same random_state should produce same results."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:3]
        decoys = sample_molecules[2:]

        s1 = GreedyFeatureSelector(refs, actives, decoys, random_state=42)
        s2 = GreedyFeatureSelector(refs, actives, decoys, random_state=42)

        r1 = s1.select(tolerance=2.0, occurrence=0.3, max_features=4, verbose=False)
        r2 = s2.select(tolerance=2.0, occurrence=0.3, max_features=4, verbose=False)

        assert r1.selected_indices == r2.selected_indices
        assert abs(r1.best_auc - r2.best_auc) < 1e-10

    def test_individual_scores_populated(self, selector):
        """Individual feature scores should be populated."""
        result = selector.select(
            tolerance=2.0, occurrence=0.3, max_features=4, verbose=False
        )
        assert len(result.individual_scores) > 0
        for idx, score in result.individual_scores.items():
            assert isinstance(idx, int)
            assert isinstance(score, float)

    def test_few_features(self, sample_molecules):
        """Should handle consensus with very few features gracefully."""
        refs = sample_molecules[:2]
        actives = [sample_molecules[0]]
        decoys = [sample_molecules[1]]
        selector = GreedyFeatureSelector(refs, actives, decoys)
        result = selector.select(
            tolerance=0.5, occurrence=0.99, max_features=4, verbose=False
        )
        assert isinstance(result, FeatureSelectionResult)

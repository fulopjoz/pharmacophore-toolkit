"""Tests for learned shape/color weighting scorer."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.evaluation import UnifiedEvaluator, EvaluationConfig
from pharmacophore.learned_scoring import LearnedScorer


@pytest.fixture
def sample_molecules():
    """Create sample molecules for testing."""
    smiles = [
        'CC(=O)O',   # Acetic acid
        'CCO',       # Ethanol
        'c1ccccc1',  # Benzene
        'CC(=O)Nc1ccc(O)cc1',  # Acetaminophen
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
def evaluator(sample_molecules):
    """Create UnifiedEvaluator with sample data."""
    refs = sample_molecules[:2]
    actives = sample_molecules[:2]
    decoys = sample_molecules[2:]
    return UnifiedEvaluator(refs, actives, decoys, random_state=42)


@pytest.fixture
def config():
    return EvaluationConfig(n_conformers=5)


class TestLearnedScorer:
    """Tests for LearnedScorer."""

    def test_fit_returns_metrics(self, evaluator, config):
        scorer = LearnedScorer(evaluator)
        result = scorer.fit(config, cv_folds=2)

        assert 'cv_auc' in result
        assert 'train_auc' in result
        assert 'coefficients' in result
        assert 0.0 <= result['cv_auc'] <= 1.0
        assert 0.0 <= result['train_auc'] <= 1.0

    def test_predict_returns_correct_length(self, evaluator, config):
        scorer = LearnedScorer(evaluator)
        scorer.fit(config, cv_folds=2)

        predictions = scorer.predict(config)
        expected_len = len(evaluator.actives) + len(evaluator.decoys)
        assert len(predictions) == expected_len

    def test_predict_in_range(self, evaluator, config):
        scorer = LearnedScorer(evaluator)
        scorer.fit(config, cv_folds=2)

        predictions = scorer.predict(config)
        assert np.all(predictions >= 0.0)
        assert np.all(predictions <= 1.0)

    def test_predict_before_fit_raises(self, evaluator, config):
        scorer = LearnedScorer(evaluator)
        with pytest.raises(RuntimeError, match="Must call fit"):
            scorer.predict(config)

    def test_feature_importance_sums_to_one(self, evaluator, config):
        scorer = LearnedScorer(evaluator)
        scorer.fit(config, cv_folds=2)

        importance = scorer.get_feature_importance()
        total = sum(imp for _, imp in importance['per_feature'])
        assert abs(total - 1.0) < 1e-6, f"Importance sums to {total}, not 1.0"

    def test_feature_importance_has_shape_color(self, evaluator, config):
        scorer = LearnedScorer(evaluator)
        scorer.fit(config, cv_folds=2)

        importance = scorer.get_feature_importance()
        assert 'shape_total' in importance
        assert 'color_total' in importance
        assert abs(importance['shape_total'] + importance['color_total'] - 1.0) < 1e-6

    def test_deterministic(self, evaluator, config):
        s1 = LearnedScorer(evaluator, random_state=42)
        s2 = LearnedScorer(evaluator, random_state=42)

        r1 = s1.fit(config, cv_folds=2)
        r2 = s2.fit(config, cv_folds=2)

        assert r1['train_auc'] == r2['train_auc']

    def test_evaluate_returns_metrics(self, evaluator, config):
        scorer = LearnedScorer(evaluator)
        scorer.fit(config, cv_folds=2)

        metrics = scorer.evaluate(config)
        assert 'roc_auc' in metrics
        assert 'bedroc' in metrics
        assert 'ef_1' in metrics

    def test_importance_before_fit_raises(self, evaluator):
        scorer = LearnedScorer(evaluator)
        with pytest.raises(RuntimeError, match="Must call fit"):
            scorer.get_feature_importance()

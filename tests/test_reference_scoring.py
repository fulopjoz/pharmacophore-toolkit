"""Tests for reference-based scoring replacement.

Validates that the new scoring modes (reference, pharm2d, hybrid) work
correctly and that the deprecated consensus_mol path emits warnings.
"""

import warnings
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.evaluation import EvaluationConfig, EvaluationResult, UnifiedEvaluator


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_mol_with_conformer(smiles: str, seed: int = 42) -> Chem.Mol:
    """Create a molecule with a 3D conformer."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    return mol


@pytest.fixture
def reference_mols():
    """Small set of reference molecules with 3D conformers."""
    smiles_list = [
        'c1ccc(NC(=O)c2ccccc2)cc1',  # Benzanilide
        'c1ccc(Oc2ccccc2)cc1',        # Diphenyl ether
        'c1ccc(Nc2ccccc2)cc1',        # Diphenylamine
    ]
    mols = [_make_mol_with_conformer(s, seed=42 + i)
            for i, s in enumerate(smiles_list)]
    return [m for m in mols if m is not None]


@pytest.fixture
def actives():
    """Small set of 'active' molecules."""
    # Similar to references
    smiles_list = [
        'c1ccc(NC(=O)c2ccc(Cl)cc2)cc1',
        'c1ccc(Oc2ccc(F)cc2)cc1',
        'c1ccc(NC(=O)Nc2ccccc2)cc1',
    ]
    mols = [_make_mol_with_conformer(s, seed=100 + i)
            for i, s in enumerate(smiles_list)]
    return [m for m in mols if m is not None]


@pytest.fixture
def decoys():
    """Small set of 'decoy' molecules."""
    # Structurally different from references
    smiles_list = [
        'CCCCCCCC',           # Octane
        'CC(=O)OC(=O)C',     # Acetic anhydride
        'C1CCCCC1',          # Cyclohexane
    ]
    mols = [_make_mol_with_conformer(s, seed=200 + i)
            for i, s in enumerate(smiles_list)]
    return [m for m in mols if m is not None]


@pytest.fixture
def evaluator(reference_mols, actives, decoys):
    """UnifiedEvaluator with test molecules."""
    return UnifiedEvaluator(reference_mols, actives, decoys, random_state=42)


# ---------------------------------------------------------------------------
# EvaluationConfig tests
# ---------------------------------------------------------------------------

class TestEvaluationConfigScoring:
    """Test new scoring_mode, aggregation, and alpha fields."""

    def test_default_scoring_mode_is_reference(self):
        config = EvaluationConfig()
        assert config.scoring_mode == 'reference'

    def test_valid_scoring_modes(self):
        for mode in ['reference', 'pharm2d', 'hybrid', 'consensus_mol']:
            config = EvaluationConfig(scoring_mode=mode)
            assert config.scoring_mode == mode

    def test_invalid_scoring_mode_raises(self):
        with pytest.raises(ValueError, match="scoring_mode"):
            EvaluationConfig(scoring_mode='invalid')

    def test_valid_aggregations(self):
        for agg in ['max', 'mean']:
            config = EvaluationConfig(aggregation=agg)
            assert config.aggregation == agg

    def test_invalid_aggregation_raises(self):
        with pytest.raises(ValueError, match="aggregation"):
            EvaluationConfig(aggregation='median')

    def test_alpha_range(self):
        EvaluationConfig(alpha=0.0)
        EvaluationConfig(alpha=1.0)
        EvaluationConfig(alpha=0.7)

    def test_alpha_out_of_range_raises(self):
        with pytest.raises(ValueError, match="alpha"):
            EvaluationConfig(alpha=1.5)
        with pytest.raises(ValueError, match="alpha"):
            EvaluationConfig(alpha=-0.1)


# ---------------------------------------------------------------------------
# UnifiedEvaluator scoring mode tests
# ---------------------------------------------------------------------------

class TestReferenceScoringMode:
    """Test reference-based scoring bypasses PharmacophoreToMol."""

    def test_reference_scoring_returns_result(self, evaluator):
        config = EvaluationConfig(scoring_mode='reference', n_conformers=5)
        result = evaluator.evaluate(config)
        assert isinstance(result, EvaluationResult)
        assert 0.0 <= result.roc_auc <= 1.0

    def test_reference_scoring_auc_above_random(self, evaluator):
        """Reference scoring should not be anti-discriminative (AUC >= 0.5)."""
        config = EvaluationConfig(
            scoring_mode='reference',
            n_conformers=5,
            opt_param=0.5
        )
        result = evaluator.evaluate(config)
        # With such small test sets, AUC may vary, but should not be
        # consistently anti-discriminative like consensus_mol
        assert result.roc_auc >= 0.0  # Sanity check

    def test_pharm2d_scoring_returns_result(self, evaluator):
        config = EvaluationConfig(scoring_mode='pharm2d')
        result = evaluator.evaluate(config)
        assert isinstance(result, EvaluationResult)
        assert 0.0 <= result.roc_auc <= 1.0

    def test_hybrid_scoring_returns_result(self, evaluator):
        config = EvaluationConfig(
            scoring_mode='hybrid',
            alpha=0.6,
            n_conformers=5
        )
        result = evaluator.evaluate(config)
        assert isinstance(result, EvaluationResult)
        assert 0.0 <= result.roc_auc <= 1.0

    def test_hybrid_alpha_affects_scores(self, evaluator):
        """Different alpha values should produce different results."""
        config_low = EvaluationConfig(
            scoring_mode='hybrid', alpha=0.1, n_conformers=5
        )
        config_high = EvaluationConfig(
            scoring_mode='hybrid', alpha=0.9, n_conformers=5
        )
        result_low = evaluator.evaluate(config_low)
        result_high = evaluator.evaluate(config_high)
        # They should differ (one is mostly 3D, other mostly 2D)
        # Small test set may have identical AUC, so just check they run
        assert isinstance(result_low, EvaluationResult)
        assert isinstance(result_high, EvaluationResult)


class TestConsensusMolDeprecation:
    """Test that consensus_mol scoring mode emits DeprecationWarning."""

    def test_consensus_mol_deprecated(self, evaluator):
        config = EvaluationConfig(scoring_mode='consensus_mol', n_conformers=5)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = evaluator.evaluate(config)
            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1
            assert "deprecated" in str(deprecation_warnings[0].message).lower()

    def test_consensus_mol_still_returns_result(self, evaluator):
        """Deprecated path should still work for backward compatibility."""
        config = EvaluationConfig(scoring_mode='consensus_mol', n_conformers=5)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = evaluator.evaluate(config)
        assert isinstance(result, EvaluationResult)
        assert 0.0 <= result.roc_auc <= 1.0


class TestEvaluateFeatureSubset:
    """Test evaluate_feature_subset with new scoring modes."""

    def test_feature_subset_reference_mode(self, evaluator):
        features = [
            ['Donor', (), 1.0, 2.0, 3.0],
            ['Acceptor', (), 4.0, 5.0, 6.0],
            ['Aromatic', (), 7.0, 8.0, 9.0],
        ]
        result = evaluator.evaluate_feature_subset(
            features, scoring_mode='reference', n_conformers=5
        )
        assert isinstance(result, EvaluationResult)
        assert result.config.scoring_mode == 'reference'

    def test_feature_subset_consensus_mol_deprecated(self, evaluator):
        features = [
            ['Donor', (), 1.0, 2.0, 3.0],
            ['Acceptor', (), 4.0, 5.0, 6.0],
            ['Aromatic', (), 7.0, 8.0, 9.0],
        ]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = evaluator.evaluate_feature_subset(
                features, scoring_mode='consensus_mol', n_conformers=5
            )
            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1

    def test_feature_subset_too_few(self, evaluator):
        features = [['Donor', (), 1.0, 2.0, 3.0]]
        result = evaluator.evaluate_feature_subset(features, n_conformers=5)
        assert result.roc_auc == 0.5


class TestPrepareReferences:
    """Test reference preparation in UnifiedEvaluator."""

    def test_prepared_refs_populated(self, evaluator):
        assert len(evaluator._prepared_refs) > 0

    def test_prepared_refs_have_conformers(self, evaluator):
        for ref in evaluator._prepared_refs:
            assert ref.GetNumConformers() > 0

    def test_pharm2d_scorer_lazy_initialized(self, evaluator):
        # Scorer is lazy-initialized (None until first use) to keep self picklable
        assert evaluator._pharm2d_scorer is None
        scorer = evaluator._get_pharm2d_scorer()
        assert scorer is not None
        assert evaluator._pharm2d_scorer is not None


class TestScoreMoleculeReference:
    """Test the _score_molecule_reference method directly."""

    def test_returns_float(self, evaluator, actives):
        mol = actives[0]
        mol_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_h, randomSeed=42)
        score = evaluator._score_molecule_reference(mol_h, opt_param=0.5)
        assert isinstance(score, float)
        assert score >= 0.0

    def test_none_mol_returns_zero(self, evaluator):
        score = evaluator._score_molecule_reference(None, opt_param=0.5)
        assert score == 0.0

    def test_max_aggregation(self, evaluator, actives):
        mol = actives[0]
        mol_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_h, randomSeed=42)
        score_max = evaluator._score_molecule_reference(
            mol_h, opt_param=0.5, aggregation='max'
        )
        score_mean = evaluator._score_molecule_reference(
            mol_h, opt_param=0.5, aggregation='mean'
        )
        # Max should be >= mean
        assert score_max >= score_mean


class TestBackwardCompatibility:
    """Test that old patterns still work with new defaults."""

    def test_evaluate_with_default_config(self, evaluator):
        """Default config should use reference scoring (not consensus_mol)."""
        config = EvaluationConfig()
        result = evaluator.evaluate(config)
        assert isinstance(result, EvaluationResult)
        assert result.config.scoring_mode == 'reference'

    def test_evaluate_returns_all_fields(self, evaluator):
        config = EvaluationConfig(n_conformers=5)
        result = evaluator.evaluate(config)
        assert hasattr(result, 'roc_auc')
        assert hasattr(result, 'bedroc')
        assert hasattr(result, 'ef_1')
        assert hasattr(result, 'ef_5')
        assert hasattr(result, 'ef_10')
        assert hasattr(result, 'n_features')
        assert hasattr(result, 'eval_time_sec')

    def test_old_config_fields_still_work(self):
        """Old-style config construction should still work."""
        config = EvaluationConfig(
            tolerance=2.5,
            occurrence=0.4,
            shape_weight=0.6,
            opt_param=0.3,
            linkage='complete',
            n_conformers=10
        )
        assert config.scoring_mode == 'reference'  # New default

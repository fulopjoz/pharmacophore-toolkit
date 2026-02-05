"""Tests for unified evaluation module."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.evaluation import (
    EvaluationConfig,
    EvaluationResult,
    UnifiedEvaluator
)


@pytest.fixture
def sample_molecules():
    """Create sample molecules for testing."""
    # Simple test molecules
    smiles = [
        'CC(=O)O',  # Acetic acid
        'CCO',      # Ethanol
        'c1ccccc1', # Benzene
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
def evaluator(sample_molecules):
    """Create UnifiedEvaluator with sample data."""
    # Use first 2 as references, split rest as actives/decoys
    refs = sample_molecules[:2]
    actives = sample_molecules[:2]
    decoys = sample_molecules[1:]

    return UnifiedEvaluator(refs, actives, decoys, random_state=42)


class TestEvaluationConfig:
    """Tests for EvaluationConfig dataclass."""

    def test_default_values(self):
        """Test default configuration values."""
        config = EvaluationConfig()

        assert config.tolerance == 2.0
        assert config.occurrence == 0.5
        assert config.shape_weight == 0.5
        assert config.opt_param == 0.5
        assert config.linkage == 'average'
        assert config.n_conformers == 25

    def test_custom_values(self):
        """Test custom configuration values."""
        config = EvaluationConfig(
            tolerance=3.0,
            occurrence=0.7,
            shape_weight=0.6,
            opt_param=0.3,
            linkage='complete',
            n_conformers=30
        )

        assert config.tolerance == 3.0
        assert config.occurrence == 0.7
        assert config.shape_weight == 0.6
        assert config.opt_param == 0.3
        assert config.linkage == 'complete'
        assert config.n_conformers == 30

    def test_tolerance_validation(self):
        """Test tolerance parameter validation."""
        # Valid range
        config = EvaluationConfig(tolerance=2.0)
        assert config.tolerance == 2.0

        # Too low
        with pytest.raises(ValueError, match="tolerance must be in"):
            EvaluationConfig(tolerance=0.1)

        # Too high
        with pytest.raises(ValueError, match="tolerance must be in"):
            EvaluationConfig(tolerance=10.0)

    def test_occurrence_validation(self):
        """Test occurrence parameter validation."""
        # Valid range
        config = EvaluationConfig(occurrence=0.5)
        assert config.occurrence == 0.5

        # Too low
        with pytest.raises(ValueError, match="occurrence must be in"):
            EvaluationConfig(occurrence=-0.1)

        # Too high
        with pytest.raises(ValueError, match="occurrence must be in"):
            EvaluationConfig(occurrence=1.5)

    def test_shape_weight_validation(self):
        """Test shape_weight parameter validation."""
        # Valid range
        config = EvaluationConfig(shape_weight=0.5)
        assert config.shape_weight == 0.5

        # Edge cases
        config = EvaluationConfig(shape_weight=0.0)
        assert config.shape_weight == 0.0

        config = EvaluationConfig(shape_weight=1.0)
        assert config.shape_weight == 1.0

        # Invalid
        with pytest.raises(ValueError, match="shape_weight must be in"):
            EvaluationConfig(shape_weight=1.5)

    def test_opt_param_validation(self):
        """Test opt_param parameter validation."""
        # Valid range
        config = EvaluationConfig(opt_param=0.5)
        assert config.opt_param == 0.5

        # Edge cases
        config = EvaluationConfig(opt_param=0.0)
        assert config.opt_param == 0.0

        config = EvaluationConfig(opt_param=1.0)
        assert config.opt_param == 1.0

        # Invalid
        with pytest.raises(ValueError, match="opt_param must be in"):
            EvaluationConfig(opt_param=-0.1)

    def test_linkage_validation(self):
        """Test linkage parameter validation."""
        # Valid options
        for linkage in ['average', 'complete', 'single', 'ward']:
            config = EvaluationConfig(linkage=linkage)
            assert config.linkage == linkage

        # Invalid
        with pytest.raises(ValueError, match="linkage must be"):
            EvaluationConfig(linkage='invalid')

    def test_n_conformers_validation(self):
        """Test n_conformers parameter validation."""
        # Valid range
        config = EvaluationConfig(n_conformers=25)
        assert config.n_conformers == 25

        # Edge cases
        config = EvaluationConfig(n_conformers=5)
        assert config.n_conformers == 5

        config = EvaluationConfig(n_conformers=50)
        assert config.n_conformers == 50

        # Too low
        with pytest.raises(ValueError, match="n_conformers must be in"):
            EvaluationConfig(n_conformers=1)

        # Too high
        with pytest.raises(ValueError, match="n_conformers must be in"):
            EvaluationConfig(n_conformers=100)


class TestEvaluationResult:
    """Tests for EvaluationResult dataclass."""

    def test_default_values(self):
        """Test default result values."""
        config = EvaluationConfig()
        result = EvaluationResult(config=config)

        assert result.config == config
        assert result.roc_auc == 0.5
        assert result.bedroc == 0.0
        assert result.ef_1 == 0.0
        assert result.ef_5 == 0.0
        assert result.ef_10 == 0.0
        assert result.n_features == 0
        assert result.eval_time_sec == 0.0

    def test_custom_values(self):
        """Test custom result values."""
        config = EvaluationConfig()
        result = EvaluationResult(
            config=config,
            roc_auc=0.85,
            bedroc=0.42,
            ef_1=8.2,
            ef_5=3.1,
            ef_10=2.5,
            n_features=12,
            eval_time_sec=15.3
        )

        assert result.roc_auc == 0.85
        assert result.bedroc == 0.42
        assert result.ef_1 == 8.2
        assert result.ef_5 == 3.1
        assert result.ef_10 == 2.5
        assert result.n_features == 12
        assert result.eval_time_sec == 15.3

    def test_summary_string(self):
        """Test summary string generation."""
        config = EvaluationConfig()
        result = EvaluationResult(
            config=config,
            roc_auc=0.85,
            bedroc=0.42,
            ef_1=8.2,
            n_features=12,
            eval_time_sec=15.3
        )

        summary = result.summary_string()
        assert 'AUC=0.8500' in summary
        assert 'BEDROC=0.4200' in summary
        assert 'EF@1%=8.2' in summary
        assert 'n_feat=12' in summary
        assert 'time=15.30s' in summary


class TestUnifiedEvaluator:
    """Tests for UnifiedEvaluator class."""

    def test_initialization(self, sample_molecules):
        """Test evaluator initialization."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:2]
        decoys = sample_molecules[1:]

        evaluator = UnifiedEvaluator(refs, actives, decoys, random_state=42)

        assert len(evaluator.reference_mols) == 2
        assert len(evaluator.actives) == 2
        assert len(evaluator.decoys) == 2
        assert len(evaluator.y_true) == 4
        assert sum(evaluator.y_true) == 2  # 2 actives
        assert evaluator.random_state == 42

    def test_repr(self, evaluator):
        """Test string representation."""
        repr_str = repr(evaluator)
        assert 'UnifiedEvaluator' in repr_str
        assert 'n_refs=2' in repr_str
        assert 'n_actives=2' in repr_str
        assert 'n_decoys=2' in repr_str

    def test_evaluate_returns_result(self, evaluator):
        """Test evaluate returns EvaluationResult."""
        config = EvaluationConfig(
            tolerance=2.0,
            occurrence=0.3,
            n_conformers=5  # Use fewer conformers for speed
        )

        result = evaluator.evaluate(config)

        assert isinstance(result, EvaluationResult)
        assert result.config == config
        assert 0.0 <= result.roc_auc <= 1.0
        assert 0.0 <= result.bedroc <= 1.0
        assert result.eval_time_sec > 0.0

    def test_degenerate_params_random_performance(self, evaluator):
        """Test degenerate consensus params give ~random AUC for consensus_mol.

        With reference scoring (default), tolerance/occurrence don't affect
        the score — queries are aligned to reference molecules directly.
        Only the deprecated consensus_mol mode routes through consensus
        feature generation, so degenerate params only hurt that path.
        """
        import warnings
        # Degenerate consensus params: very tight tolerance + high occurrence → few/no features
        config = EvaluationConfig(
            tolerance=0.5,
            occurrence=0.99,
            n_conformers=5,
            scoring_mode='consensus_mol',
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = evaluator.evaluate(config)

        # Should get random performance (AUC ~ 0.5)
        assert 0.4 <= result.roc_auc <= 0.6

    def test_timing_populated(self, evaluator):
        """Test evaluation timing is tracked."""
        config = EvaluationConfig(n_conformers=5)

        result = evaluator.evaluate(config)

        assert result.eval_time_sec > 0.0
        assert result.eval_time_sec < 60.0  # Should be reasonably fast

    def test_all_linkage_methods(self, evaluator):
        """Test all linkage methods work."""
        linkage_methods = ['average', 'complete', 'single', 'ward']

        for linkage in linkage_methods:
            config = EvaluationConfig(
                tolerance=2.0,
                occurrence=0.3,
                linkage=linkage,
                n_conformers=5
            )

            result = evaluator.evaluate(config)

            assert isinstance(result, EvaluationResult)
            assert result.config.linkage == linkage
            assert 0.0 <= result.roc_auc <= 1.0

    def test_evaluate_feature_subset(self, evaluator):
        """Test direct feature subset evaluation."""
        # Create simple features
        features = [
            ['Donor', (), 0.0, 0.0, 0.0],
            ['Acceptor', (), 1.0, 1.0, 1.0],
            ['Aromatic', (), 2.0, 2.0, 2.0]
        ]

        result = evaluator.evaluate_feature_subset(
            features,
            shape_weight=0.5,
            opt_param=0.5,
            n_conformers=5
        )

        assert isinstance(result, EvaluationResult)
        assert result.n_features == 3
        assert 0.0 <= result.roc_auc <= 1.0
        assert result.eval_time_sec > 0.0

    def test_evaluate_feature_subset_too_few(self, evaluator):
        """Test feature subset with too few features."""
        # Only 1 feature (degenerate)
        features = [['Donor', (), 0.0, 0.0, 0.0]]

        result = evaluator.evaluate_feature_subset(
            features,
            n_conformers=5
        )

        # Should return random performance
        assert result.roc_auc == 0.5
        assert result.bedroc == 0.0
        assert result.n_features == 1

    def test_reproducibility(self, sample_molecules):
        """Test reproducibility with same random_state."""
        refs = sample_molecules[:2]
        actives = sample_molecules[:2]
        decoys = sample_molecules[1:]

        # Create two evaluators with same seed
        eval1 = UnifiedEvaluator(refs, actives, decoys, random_state=42)
        eval2 = UnifiedEvaluator(refs, actives, decoys, random_state=42)

        config = EvaluationConfig(n_conformers=5)

        result1 = eval1.evaluate(config)
        result2 = eval2.evaluate(config)

        # Results should be identical
        assert result1.roc_auc == result2.roc_auc
        assert result1.bedroc == result2.bedroc
        assert result1.n_features == result2.n_features

    def test_conformer_caching(self, evaluator):
        """Test conformer generation is cached."""
        config = EvaluationConfig(n_conformers=5)

        # First evaluation - should generate conformers
        result1 = evaluator.evaluate(config)
        time1 = result1.eval_time_sec

        # Second evaluation - should use cached conformers
        result2 = evaluator.evaluate(config)
        time2 = result2.eval_time_sec

        # Second call should be faster (or at least not significantly slower)
        assert time2 <= time1 * 1.5  # Allow 50% margin for variance

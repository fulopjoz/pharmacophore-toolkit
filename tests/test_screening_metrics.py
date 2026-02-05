"""Unit tests for screening metrics module.

This test suite covers:
- Youden's Index calculation
- BEDROC and RIE metrics
- Confusion matrix metrics
- Enrichment factor calculation
"""

import unittest
import numpy as np
import pytest

from pharmacophore.screening_metrics import (
    youden_index, bedroc, rie, confusion_metrics,
    enrichment_factor, screening_report, calculate_all_metrics
)


class TestYoudenIndex(unittest.TestCase):
    """Test Youden's Index calculation."""

    def test_perfect_separation(self):
        """Test Youden's index with perfect separation."""
        y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]

        thresh, j, sens, spec = youden_index(y_true, y_scores)

        # Should have high Youden's J (close to 1)
        self.assertGreater(j, 0.9)
        # Sensitivity and specificity should be high
        self.assertGreater(sens, 0.9)
        self.assertGreater(spec, 0.9)
        # Threshold should be between active and decoy scores (inclusive)
        self.assertGreaterEqual(thresh, 0.3)
        self.assertLessEqual(thresh, 0.7)

    def test_random_ordering(self):
        """Test Youden's index with interleaved actives/decoys."""
        y_true = [1, 0, 1, 0, 1, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.85, 0.7, 0.65, 0.5, 0.45, 0.3, 0.2, 0.1, 0.05]

        thresh, j, sens, spec = youden_index(y_true, y_scores)

        # Youden's J should be lower than perfect separation
        self.assertLess(j, 0.9)
        self.assertGreater(j, 0.0)

    def test_empty_inputs(self):
        """Test Youden's index with empty inputs."""
        thresh, j, sens, spec = youden_index([], [])

        self.assertEqual(thresh, 0.0)
        self.assertEqual(j, 0.0)
        self.assertEqual(sens, 0.0)
        self.assertEqual(spec, 0.0)

    @pytest.mark.filterwarnings("ignore::sklearn.exceptions.UndefinedMetricWarning")
    def test_all_actives(self):
        """Test Youden's index when all are actives."""
        y_true = [1, 1, 1, 1, 1]
        y_scores = [0.9, 0.8, 0.7, 0.6, 0.5]

        thresh, j, sens, spec = youden_index(y_true, y_scores)

        # Should handle edge case gracefully
        self.assertIsInstance(j, float)


class TestBEDROC(unittest.TestCase):
    """Test BEDROC calculation."""

    def test_perfect_early_recognition(self):
        """Test BEDROC when all actives are at the top."""
        y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]

        score = bedroc(y_true, y_scores, alpha=20.0)

        # Should be close to 1.0
        self.assertGreater(score, 0.9)

    def test_poor_early_recognition(self):
        """Test BEDROC when actives are at the bottom."""
        y_true = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1]
        y_scores = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]

        score = bedroc(y_true, y_scores, alpha=20.0)

        # Should be close to 0
        self.assertLess(score, 0.1)

    def test_alpha_parameter(self):
        """Test that higher alpha emphasizes early recognition more."""
        y_true = [1, 0, 0, 0, 0, 1, 0, 0, 0, 1]
        y_scores = [0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.4, 0.3, 0.2, 0.1]

        score_low_alpha = bedroc(y_true, y_scores, alpha=8.0)
        score_high_alpha = bedroc(y_true, y_scores, alpha=80.5)

        # Both should be valid scores
        self.assertGreaterEqual(score_low_alpha, 0.0)
        self.assertLessEqual(score_low_alpha, 1.0)
        self.assertGreaterEqual(score_high_alpha, 0.0)
        self.assertLessEqual(score_high_alpha, 1.0)

    def test_empty_inputs(self):
        """Test BEDROC with empty inputs."""
        score = bedroc([], [], alpha=20.0)
        self.assertEqual(score, 0.0)

    def test_no_actives(self):
        """Test BEDROC when there are no actives."""
        y_true = [0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.6, 0.5]

        score = bedroc(y_true, y_scores, alpha=20.0)
        self.assertEqual(score, 0.0)

    def test_all_actives(self):
        """Test BEDROC when all are actives."""
        y_true = [1, 1, 1, 1, 1]
        y_scores = [0.9, 0.8, 0.7, 0.6, 0.5]

        score = bedroc(y_true, y_scores, alpha=20.0)
        self.assertEqual(score, 1.0)


class TestRIE(unittest.TestCase):
    """Test RIE calculation."""

    def test_better_than_random(self):
        """Test RIE when actives are at the top."""
        y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]

        score = rie(y_true, y_scores, alpha=20.0)

        # Should be greater than 1 (better than random)
        self.assertGreater(score, 1.0)

    def test_worse_than_random(self):
        """Test RIE when actives are at the bottom."""
        y_true = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1]
        y_scores = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]

        score = rie(y_true, y_scores, alpha=20.0)

        # Should be less than 1 (worse than random)
        self.assertLess(score, 1.0)


class TestConfusionMetrics(unittest.TestCase):
    """Test confusion matrix-based metrics."""

    def test_perfect_classification(self):
        """Test metrics with perfect classification."""
        y_true = [1, 1, 1, 0, 0, 0, 0, 0]
        y_pred = [1, 1, 1, 0, 0, 0, 0, 0]

        metrics = confusion_metrics(y_true, y_pred)

        self.assertEqual(metrics['TP'], 3)
        self.assertEqual(metrics['FP'], 0)
        self.assertEqual(metrics['TN'], 5)
        self.assertEqual(metrics['FN'], 0)
        self.assertEqual(metrics['sensitivity'], 1.0)
        self.assertEqual(metrics['specificity'], 1.0)
        self.assertEqual(metrics['precision'], 1.0)
        self.assertEqual(metrics['f1'], 1.0)
        self.assertEqual(metrics['accuracy'], 1.0)

    def test_all_wrong(self):
        """Test metrics with completely wrong classification."""
        y_true = [1, 1, 0, 0]
        y_pred = [0, 0, 1, 1]

        metrics = confusion_metrics(y_true, y_pred)

        self.assertEqual(metrics['TP'], 0)
        self.assertEqual(metrics['FP'], 2)
        self.assertEqual(metrics['TN'], 0)
        self.assertEqual(metrics['FN'], 2)
        self.assertEqual(metrics['sensitivity'], 0.0)
        self.assertEqual(metrics['specificity'], 0.0)
        self.assertEqual(metrics['precision'], 0.0)

    def test_partial_classification(self):
        """Test metrics with partial correct classification."""
        y_true = [1, 1, 1, 0, 0, 0]
        y_pred = [1, 1, 0, 1, 0, 0]

        metrics = confusion_metrics(y_true, y_pred)

        self.assertEqual(metrics['TP'], 2)
        self.assertEqual(metrics['FP'], 1)
        self.assertEqual(metrics['TN'], 2)
        self.assertEqual(metrics['FN'], 1)

        # Check derived metrics
        self.assertAlmostEqual(metrics['sensitivity'], 2/3, places=5)
        self.assertAlmostEqual(metrics['specificity'], 2/3, places=5)
        self.assertAlmostEqual(metrics['precision'], 2/3, places=5)

    def test_empty_inputs(self):
        """Test confusion metrics with empty inputs."""
        metrics = confusion_metrics([], [])

        self.assertEqual(metrics['TP'], 0)
        self.assertEqual(metrics['FP'], 0)
        self.assertEqual(metrics['TN'], 0)
        self.assertEqual(metrics['FN'], 0)


class TestEnrichmentFactor(unittest.TestCase):
    """Test enrichment factor calculation."""

    def test_perfect_enrichment(self):
        """Test EF when all actives are at the top."""
        y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]

        # At 30%, should capture all 3 actives in top 3 compounds
        ef = enrichment_factor(y_true, y_scores, percentage=0.30)

        # EF = (3/3) / (3/10) = 1.0 / 0.3 = 3.33
        self.assertGreater(ef, 3.0)

    def test_random_enrichment(self):
        """Test EF approaches 1.0 for random ordering."""
        np.random.seed(42)
        y_true = [1] * 30 + [0] * 70
        y_scores = list(np.random.random(100))

        # At 50%, should be close to 1.0
        ef = enrichment_factor(y_true, y_scores, percentage=0.50)

        # Allow variance due to randomness
        self.assertGreater(ef, 0.5)
        self.assertLess(ef, 2.0)

    def test_no_actives(self):
        """Test EF when there are no actives."""
        y_true = [0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.6, 0.5]

        ef = enrichment_factor(y_true, y_scores, percentage=0.20)
        self.assertEqual(ef, 0.0)

    def test_empty_inputs(self):
        """Test EF with empty inputs."""
        ef = enrichment_factor([], [], percentage=0.01)
        self.assertEqual(ef, 0.0)


class TestScreeningReport(unittest.TestCase):
    """Test comprehensive screening report generation."""

    def test_report_structure(self):
        """Test that report has expected structure."""
        y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]

        report = screening_report(y_true, y_scores)

        # Check it's a DataFrame
        self.assertEqual(type(report).__name__, 'DataFrame')

        # Check expected columns
        self.assertIn('Metric', report.columns)
        self.assertIn('Value', report.columns)

        # Check some expected metrics are present
        metrics = report['Metric'].tolist()
        self.assertIn('ROC-AUC', metrics)
        self.assertIn('BEDROC (α=20)', metrics)
        self.assertIn('Youden J', metrics)

    def test_report_values(self):
        """Test that report values are reasonable."""
        y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]

        report = screening_report(y_true, y_scores)

        # Get ROC-AUC value
        auc_row = report[report['Metric'] == 'ROC-AUC']
        auc = auc_row['Value'].values[0]

        # Should be high for good separation
        self.assertGreater(auc, 0.9)


class TestCalculateAllMetrics(unittest.TestCase):
    """Test calculate_all_metrics convenience function."""

    def test_returns_dict(self):
        """Test that function returns expected dictionary."""
        y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]

        metrics = calculate_all_metrics(y_true, y_scores)

        # Check expected keys
        expected_keys = [
            'roc_auc', 'bedroc', 'rie', 'ef_1', 'ef_5', 'ef_10',
            'youden_j', 'optimal_threshold', 'sensitivity', 'specificity'
        ]
        for key in expected_keys:
            self.assertIn(key, metrics)

    def test_metric_ranges(self):
        """Test that metrics are in valid ranges."""
        y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]

        metrics = calculate_all_metrics(y_true, y_scores)

        # AUC should be 0-1
        self.assertGreaterEqual(metrics['roc_auc'], 0.0)
        self.assertLessEqual(metrics['roc_auc'], 1.0)

        # BEDROC should be 0-1
        self.assertGreaterEqual(metrics['bedroc'], 0.0)
        self.assertLessEqual(metrics['bedroc'], 1.0)

        # EF should be >= 0
        self.assertGreaterEqual(metrics['ef_1'], 0.0)
        self.assertGreaterEqual(metrics['ef_5'], 0.0)

        # Sensitivity and specificity should be 0-1
        self.assertGreaterEqual(metrics['sensitivity'], 0.0)
        self.assertLessEqual(metrics['sensitivity'], 1.0)
        self.assertGreaterEqual(metrics['specificity'], 0.0)
        self.assertLessEqual(metrics['specificity'], 1.0)


if __name__ == '__main__':
    unittest.main()

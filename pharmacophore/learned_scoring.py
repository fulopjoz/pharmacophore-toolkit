"""Learned shape/color weighting for pharmacophore virtual screening.

Instead of fixed ``opt_param`` (equal shape/color weighting), learns the
optimal per-reference shape vs color balance from labeled data using
logistic regression over the 2*N_refs feature vector.

Literature basis: PharmRF (Kumar 2022) uses Random Forest; we use
logistic regression for interpretable coefficients and lower overfitting
risk on small datasets.

Example:
    >>> from pharmacophore.learned_scoring import LearnedScorer
    >>> from pharmacophore.evaluation import UnifiedEvaluator, EvaluationConfig
    >>> evaluator = UnifiedEvaluator(refs, actives, decoys)
    >>> scorer = LearnedScorer(evaluator)
    >>> fit_result = scorer.fit(EvaluationConfig(n_conformers=15))
    >>> print(f"CV AUC: {fit_result['cv_auc']:.4f}")
"""

from typing import List, Optional
import logging
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score

logger = logging.getLogger(__name__)


class LearnedScorer:
    """Logistic regression scorer over per-reference shape/color features.

    Builds a feature matrix ``[shape_ref0, color_ref0, shape_ref1, ...]``
    from ScoreBreakdown objects, fits a logistic regression, and predicts
    P(active) as the virtual screening score.

    Args:
        evaluator: UnifiedEvaluator instance (provides scoring + labels).
        C: Regularization strength (inverse). Default 1.0.
        random_state: Seed for reproducibility.
    """

    def __init__(self, evaluator, C: float = 1.0, random_state: int = 42):
        self.evaluator = evaluator
        self.C = C
        self.random_state = random_state
        self._model = None
        self._scaler = None
        self._fitted = False
        self._n_refs = len(evaluator._prepared_refs)

    def _breakdowns_to_matrix(self, breakdowns) -> np.ndarray:
        """Convert ScoreBreakdown list to feature matrix.

        Each row is [shape_ref0, color_ref0, shape_ref1, color_ref1, ...].
        """
        n_samples = len(breakdowns)
        n_features = 2 * self._n_refs
        X = np.zeros((n_samples, n_features))

        for i, bd in enumerate(breakdowns):
            for j in range(min(len(bd.shape_scores), self._n_refs)):
                X[i, 2 * j] = bd.shape_scores[j]
                X[i, 2 * j + 1] = bd.color_scores[j]

        return X

    def fit(self, config, cv_folds: int = 5) -> dict:
        """Fit logistic regression on per-reference shape/color features.

        Args:
            config: EvaluationConfig for scoring parameters.
            cv_folds: Number of cross-validation folds.

        Returns:
            Dict with 'cv_auc', 'train_auc', 'coefficients'.
        """
        breakdowns = self.evaluator.score_all_with_breakdown(config)
        X = self._breakdowns_to_matrix(breakdowns)
        y = self.evaluator.y_true

        # Scale features
        self._scaler = StandardScaler()
        X_scaled = self._scaler.fit_transform(X)

        # Fit model
        self._model = LogisticRegression(
            C=self.C,
            random_state=self.random_state,
            max_iter=1000,
            solver='lbfgs',
        )
        self._model.fit(X_scaled, y)
        self._fitted = True

        # Cross-validated AUC
        cv = StratifiedKFold(
            n_splits=min(cv_folds, int(y.sum())),
            shuffle=True,
            random_state=self.random_state,
        )
        cv_scores = cross_val_score(
            LogisticRegression(
                C=self.C, random_state=self.random_state,
                max_iter=1000, solver='lbfgs',
            ),
            X_scaled, y, cv=cv, scoring='roc_auc',
        )
        cv_auc = float(np.mean(cv_scores))

        # Train AUC
        train_proba = self._model.predict_proba(X_scaled)[:, 1]
        train_auc = float(roc_auc_score(y, train_proba))

        coefficients = self._model.coef_[0].tolist()

        logger.info(
            "LearnedScorer fit: CV AUC=%.4f, Train AUC=%.4f", cv_auc, train_auc
        )

        return {
            'cv_auc': cv_auc,
            'train_auc': train_auc,
            'coefficients': coefficients,
        }

    def predict(self, config) -> np.ndarray:
        """Predict P(active) for all molecules.

        Args:
            config: EvaluationConfig for scoring parameters.

        Returns:
            Array of predicted probabilities.

        Raises:
            RuntimeError: If fit() has not been called.
        """
        if not self._fitted:
            raise RuntimeError("Must call fit() before predict()")

        breakdowns = self.evaluator.score_all_with_breakdown(config)
        X = self._breakdowns_to_matrix(breakdowns)
        X_scaled = self._scaler.transform(X)
        return self._model.predict_proba(X_scaled)[:, 1]

    def evaluate(self, config) -> dict:
        """Predict and compute full screening metrics.

        Args:
            config: EvaluationConfig for scoring parameters.

        Returns:
            Dict with roc_auc, bedroc, ef_1, ef_5, ef_10.
        """
        from .screening_metrics import calculate_all_metrics

        scores = self.predict(config)
        return calculate_all_metrics(
            self.evaluator.y_true.tolist(), scores.tolist()
        )

    def get_feature_importance(self) -> dict:
        """Get per-feature importance (absolute coefficient, normalized).

        Returns:
            Dict with 'per_feature' (list of (name, importance)) and
            'shape_total', 'color_total' aggregates.

        Raises:
            RuntimeError: If fit() has not been called.
        """
        if not self._fitted:
            raise RuntimeError("Must call fit() before get_feature_importance()")

        coefs = np.abs(self._model.coef_[0])
        total = coefs.sum()
        if total < 1e-12:
            normalized = np.ones_like(coefs) / len(coefs)
        else:
            normalized = coefs / total

        per_feature = []
        shape_total = 0.0
        color_total = 0.0

        for j in range(self._n_refs):
            shape_imp = normalized[2 * j] if 2 * j < len(normalized) else 0.0
            color_imp = normalized[2 * j + 1] if 2 * j + 1 < len(normalized) else 0.0
            per_feature.append((f'shape_ref{j}', float(shape_imp)))
            per_feature.append((f'color_ref{j}', float(color_imp)))
            shape_total += shape_imp
            color_total += color_imp

        return {
            'per_feature': per_feature,
            'shape_total': float(shape_total),
            'color_total': float(color_total),
        }

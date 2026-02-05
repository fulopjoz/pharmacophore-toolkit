"""Virtual screening metrics for pharmacophore-based drug discovery.

This module provides comprehensive metrics for evaluating virtual screening
performance, including Youden's Index for optimal threshold selection and
BEDROC for early recognition assessment.

Key Metrics:
- Youden's Index: J = sensitivity + specificity - 1 (optimal threshold finder)
- BEDROC: Boltzmann-Enhanced Discrimination ROC (early recognition)
- Enrichment Factor: Actives enrichment at various percentiles
- Confusion Metrics: TP, FP, TN, FN, precision, recall, F1, MCC
"""

import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Optional, Union
from sklearn.metrics import (
    roc_curve, roc_auc_score, confusion_matrix,
    precision_score, recall_score, f1_score, matthews_corrcoef
)
from rdkit.ML.Scoring import Scoring as RDKitScoring


def youden_index(
    y_true: List[int],
    y_scores: List[float]
) -> Tuple[float, float, float, float]:
    """Find optimal threshold using Youden's J statistic.

    Youden's J = sensitivity + specificity - 1

    The optimal threshold is where J is maximized, providing the best
    trade-off between sensitivity (true positive rate) and specificity
    (true negative rate).

    Args:
        y_true: Binary labels (1=active, 0=decoy).
        y_scores: Predicted scores (higher = more likely active).

    Returns:
        Tuple of (optimal_threshold, max_youden_j, sensitivity, specificity).

    Example:
        >>> y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        >>> y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]
        >>> thresh, j, sens, spec = youden_index(y_true, y_scores)
        >>> print(f"Optimal threshold: {thresh:.2f}, J: {j:.2f}")
    """
    if len(y_true) == 0 or len(y_scores) == 0:
        return 0.0, 0.0, 0.0, 0.0

    y_true = np.array(y_true)
    y_scores = np.array(y_scores)

    # Get ROC curve points
    fpr, tpr, thresholds = roc_curve(y_true, y_scores)

    # Calculate Youden's J for each threshold
    # J = TPR - FPR = sensitivity + specificity - 1
    youden_j = tpr - fpr

    # Find optimal threshold (maximum J)
    optimal_idx = np.argmax(youden_j)
    optimal_threshold = thresholds[optimal_idx]
    max_j = youden_j[optimal_idx]

    # Calculate sensitivity and specificity at optimal threshold
    sensitivity = tpr[optimal_idx]
    specificity = 1 - fpr[optimal_idx]

    return float(optimal_threshold), float(max_j), float(sensitivity), float(specificity)


def bedroc(
    y_true: List[int],
    y_scores: List[float],
    alpha: float = 20.0
) -> float:
    """Calculate Boltzmann-Enhanced Discrimination ROC (BEDROC).

    Uses RDKit's validated implementation of BEDROC.

    BEDROC emphasizes early recognition - the ability to identify actives
    in the top portion of a ranked list. The alpha parameter controls
    the degree of early recognition emphasis:

    - alpha=20.0: Emphasizes top ~8% of ranked list (recommended for VS)
    - alpha=80.5: Emphasizes top ~2% (very stringent)
    - alpha=8.0: Emphasizes top ~20% (more permissive)

    BEDROC ranges from 0 to 1:
    - 1.0: All actives are at the very top of the ranked list
    - 0.5: Random ordering
    - 0.0: All actives are at the very bottom

    Reference: Truchon & Bayly (2007) J. Chem. Inf. Model. 47, 488-508

    Args:
        y_true: Binary labels (1=active, 0=decoy).
        y_scores: Predicted scores (higher = more likely active).
        alpha: Exponential weighting parameter (default: 20.0).

    Returns:
        BEDROC score (0-1, higher = better early recognition).

    Example:
        >>> y_true = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        >>> y_scores = [0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]
        >>> score = bedroc(y_true, y_scores, alpha=20.0)
        >>> print(f"BEDROC: {score:.3f}")  # Should be high (actives at top)
    """
    if len(y_true) == 0 or len(y_scores) == 0:
        return 0.0

    y_true = np.array(y_true)
    y_scores = np.array(y_scores)

    n_actives = int(np.sum(y_true))

    if n_actives == 0:
        return 0.0
    if n_actives == len(y_true):
        return 1.0

    # Sort by scores (descending) - RDKit expects low indices to be "better"
    sorted_indices = np.argsort(y_scores)[::-1]
    sorted_labels = y_true[sorted_indices]

    # RDKit format: list of [score, label] sorted by score descending
    # The col parameter tells which index is the label
    scores_for_rdkit = [[float(y_scores[i]), int(y_true[i])]
                        for i in sorted_indices]

    try:
        return float(RDKitScoring.CalcBEDROC(scores_for_rdkit, col=1, alpha=alpha))
    except Exception:
        return 0.0


def rie(
    y_true: List[int],
    y_scores: List[float],
    alpha: float = 20.0
) -> float:
    """Calculate Robust Initial Enhancement (RIE).

    Uses RDKit's validated implementation of RIE.

    RIE is related to BEDROC and measures early recognition ability.
    RIE > 1 indicates better than random early recognition.

    Args:
        y_true: Binary labels (1=active, 0=decoy).
        y_scores: Predicted scores (higher = more likely active).
        alpha: Exponential weighting parameter (default: 20.0).

    Returns:
        RIE score (1.0 = random, higher = better early recognition).
    """
    if len(y_true) == 0 or len(y_scores) == 0:
        return 1.0

    y_true = np.array(y_true)
    y_scores = np.array(y_scores)

    n_actives = int(np.sum(y_true))

    if n_actives == 0 or n_actives == len(y_true):
        return 1.0

    # Sort by scores (descending)
    sorted_indices = np.argsort(y_scores)[::-1]

    # RDKit format: list of [score, label] sorted by score descending
    scores_for_rdkit = [[float(y_scores[i]), int(y_true[i])]
                        for i in sorted_indices]

    try:
        return float(RDKitScoring.CalcRIE(scores_for_rdkit, col=1, alpha=alpha))
    except Exception:
        return 1.0


def confusion_metrics(
    y_true: List[int],
    y_pred: List[int]
) -> Dict[str, float]:
    """Calculate comprehensive confusion matrix-based metrics.

    Args:
        y_true: True binary labels (1=active, 0=decoy).
        y_pred: Predicted binary labels (1=active, 0=decoy).

    Returns:
        Dict with keys: TP, FP, TN, FN, sensitivity (recall), specificity,
        precision, npv, f1, mcc, accuracy, balanced_accuracy.
    """
    if len(y_true) == 0 or len(y_pred) == 0:
        return {
            'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0,
            'sensitivity': 0.0, 'specificity': 0.0,
            'precision': 0.0, 'npv': 0.0,
            'f1': 0.0, 'mcc': 0.0,
            'accuracy': 0.0, 'balanced_accuracy': 0.0
        }

    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    # Confusion matrix: [[TN, FP], [FN, TP]]
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel() if cm.size == 4 else (0, 0, 0, 0)

    # Derived metrics
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0  # TPR, Recall
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0  # TNR
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0    # PPV
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0          # Negative predictive value

    f1 = f1_score(y_true, y_pred, zero_division=0)
    mcc = matthews_corrcoef(y_true, y_pred) if len(set(y_true)) > 1 else 0.0

    accuracy = (tp + tn) / len(y_true) if len(y_true) > 0 else 0.0
    balanced_accuracy = (sensitivity + specificity) / 2

    return {
        'TP': int(tp),
        'FP': int(fp),
        'TN': int(tn),
        'FN': int(fn),
        'sensitivity': float(sensitivity),
        'specificity': float(specificity),
        'precision': float(precision),
        'npv': float(npv),
        'f1': float(f1),
        'mcc': float(mcc),
        'accuracy': float(accuracy),
        'balanced_accuracy': float(balanced_accuracy)
    }


def enrichment_factor(
    y_true: List[int],
    y_scores: List[float],
    percentage: float = 0.01
) -> float:
    """Calculate enrichment factor at a given percentage.

    EF = (actives_in_top_X% / total_actives) / (X% / 100%)

    Random selection yields EF = 1.0. Higher values indicate enrichment.

    Args:
        y_true: Binary labels (1=active, 0=decoy).
        y_scores: Predicted scores (higher = more likely active).
        percentage: Top percentage to consider (default: 0.01 = 1%).

    Returns:
        Enrichment factor value.
    """
    if len(y_true) == 0 or len(y_scores) == 0:
        return 0.0

    y_true = np.array(y_true)
    y_scores = np.array(y_scores)

    n_total = len(y_true)
    n_actives = np.sum(y_true)

    if n_actives == 0:
        return 0.0

    # Sort by score (descending)
    sorted_indices = np.argsort(y_scores)[::-1]
    sorted_labels = y_true[sorted_indices]

    # Get top percentage
    n_top = max(1, int(n_total * percentage))
    actives_in_top = np.sum(sorted_labels[:n_top])

    # Calculate EF
    expected_random = n_actives / n_total
    observed = actives_in_top / n_top

    if expected_random == 0:
        return 0.0

    return float(observed / expected_random)


def screening_report(
    y_true: List[int],
    y_scores: List[float],
    percentiles: List[float] = None
) -> pd.DataFrame:
    """Generate comprehensive screening metrics report.

    Calculates all major metrics at once and returns a DataFrame.

    Args:
        y_true: Binary labels (1=active, 0=decoy).
        y_scores: Predicted scores (higher = more likely active).
        percentiles: List of percentiles for EF (default: [0.01, 0.05, 0.10]).

    Returns:
        DataFrame with all screening metrics.
    """
    if percentiles is None:
        percentiles = [0.01, 0.05, 0.10]

    # Basic info
    n_total = len(y_true)
    n_actives = sum(y_true)
    n_decoys = n_total - n_actives

    # ROC-AUC
    try:
        auc = roc_auc_score(y_true, y_scores)
    except Exception:
        auc = 0.5

    # Youden's index
    thresh, j, sens, spec = youden_index(y_true, y_scores)

    # BEDROC
    bedroc_20 = bedroc(y_true, y_scores, alpha=20.0)
    bedroc_80 = bedroc(y_true, y_scores, alpha=80.5)

    # RIE
    rie_20 = rie(y_true, y_scores, alpha=20.0)

    # Enrichment factors
    ef_values = {f'EF@{int(p*100)}%': enrichment_factor(y_true, y_scores, p)
                 for p in percentiles}

    # Confusion metrics at optimal threshold
    y_pred = [1 if s >= thresh else 0 for s in y_scores]
    conf = confusion_metrics(y_true, y_pred)

    # Build report
    report_data = {
        'Metric': [
            'N Total', 'N Actives', 'N Decoys', 'Active Ratio',
            'ROC-AUC',
            'BEDROC (α=20)', 'BEDROC (α=80.5)', 'RIE (α=20)',
            *ef_values.keys(),
            'Youden J', 'Optimal Threshold', 'Sensitivity', 'Specificity',
            'Precision', 'F1 Score', 'MCC', 'Balanced Accuracy'
        ],
        'Value': [
            n_total, n_actives, n_decoys, n_actives / n_total if n_total > 0 else 0,
            auc,
            bedroc_20, bedroc_80, rie_20,
            *ef_values.values(),
            j, thresh, sens, spec,
            conf['precision'], conf['f1'], conf['mcc'], conf['balanced_accuracy']
        ]
    }

    return pd.DataFrame(report_data)


def calculate_all_metrics(
    y_true: List[int],
    y_scores: List[float]
) -> Dict[str, float]:
    """Calculate all screening metrics in one call.

    Convenience function that returns all metrics as a dictionary,
    suitable for benchmarking and comparison.

    Args:
        y_true: Binary labels (1=active, 0=decoy).
        y_scores: Predicted scores (higher = more likely active).

    Returns:
        Dict with all metric values.
    """
    # ROC-AUC
    try:
        auc = roc_auc_score(y_true, y_scores)
    except Exception:
        auc = 0.5

    # Youden's index
    thresh, j, sens, spec = youden_index(y_true, y_scores)

    # BEDROC and RIE
    bedroc_20 = bedroc(y_true, y_scores, alpha=20.0)
    rie_20 = rie(y_true, y_scores, alpha=20.0)

    # Enrichment factors
    ef_1 = enrichment_factor(y_true, y_scores, percentage=0.01)
    ef_5 = enrichment_factor(y_true, y_scores, percentage=0.05)
    ef_10 = enrichment_factor(y_true, y_scores, percentage=0.10)

    return {
        'roc_auc': auc,
        'bedroc': bedroc_20,
        'rie': rie_20,
        'ef_1': ef_1,
        'ef_5': ef_5,
        'ef_10': ef_10,
        'youden_j': j,
        'optimal_threshold': thresh,
        'sensitivity': sens,
        'specificity': spec
    }


if __name__ == "__main__":
    import doctest
    doctest.testmod()

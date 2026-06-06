"""Baseline scorers: equal-weight color similarity and the S3 discrimination-weighted color.

Both operate on the per-P4-type similarity-to-reference features (the Tier-1 representation).
- equal_weight: unsupervised — mean of the 6 per-type similarities (ignores train_idx).
- s3_weighted:  supervised  — logistic regression fit on train_idx (learns the donor-up /
  cation-down discrimination weights), predicts P(active) on test_idx.
"""
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from harness import register, BenchData, SEED


@register("equal_weight")
def equal_weight(data: BenchData, train_idx, test_idx) -> np.ndarray:
    return data.p4_sim_to_refs()[test_idx].mean(axis=1)


@register("s3_weighted")
def s3_weighted(data: BenchData, train_idx, test_idx) -> np.ndarray:
    X = data.p4_sim_to_refs()
    sc = StandardScaler().fit(X[train_idx])
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
    lr.fit(sc.transform(X[train_idx]), data.y[train_idx])
    return lr.predict_proba(sc.transform(X[test_idx]))[:, 1]

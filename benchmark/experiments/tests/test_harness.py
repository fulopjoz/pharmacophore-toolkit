"""Foundation tests for the comparison harness (the deep module)."""
import numpy as np
import pytest

import harness


def test_benchdata_loads_actives_and_decoys():
    d = harness.BenchData.load_default()
    assert len(d.smiles) == len(d.y)
    assert d.y.sum() >= 70 and (d.y == 0).sum() >= 400      # ~74 actives, ~500 decoys
    assert set(np.unique(d.y)) == {0, 1}


def test_featurizations_have_right_shape_and_cache():
    d = harness.BenchData.load_default()
    p4 = d.p4_counts()
    assert p4.shape == (len(d.y), 6)
    assert d.p4_counts() is p4                               # lazily cached (same object)
    assert d.morgan().shape[0] == len(d.y)


def test_register_decorator_populates_registry():
    harness.discover()
    assert {"equal_weight", "s3_weighted"} <= set(harness.REGISTRY)


def test_evaluate_returns_bounded_metrics():
    d = harness.BenchData.load_default()
    harness.discover()
    m = harness.evaluate("equal_weight", d)
    assert {"AUC", "BEDROC", "EF1%", "EF5%"} <= set(m)
    assert 0.0 <= m["AUC"] <= 1.0 and 0.0 <= m["BEDROC"] <= 1.0


def test_s3_beats_equal_weight():
    """Reproduces the Tier-1 result: discrimination-weighting > equal-weighting."""
    d = harness.BenchData.load_default()
    harness.discover()
    eq = harness.evaluate("equal_weight", d)
    s3 = harness.evaluate("s3_weighted", d)
    assert s3["AUC"] > eq["AUC"]


def test_bootstrap_delta_excludes_zero_when_s3_better():
    d = harness.BenchData.load_default()
    harness.discover()
    oof_eq = harness.evaluate_oof("equal_weight", d)
    oof_s3 = harness.evaluate_oof("s3_weighted", d)
    med, lo, hi = harness.bootstrap_delta(oof_s3, oof_eq, d.y, "AUC", n=300)
    assert med > 0 and lo > 0                                # S3 AUC reliably above equal-weight

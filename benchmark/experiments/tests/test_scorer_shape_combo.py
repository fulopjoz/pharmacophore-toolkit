"""Tests for the RDKit shape+color combo scorer (production-ROCS proxy).

TDD — these tests were written BEFORE the implementation.
"""
import numpy as np
import pytest

import harness


@pytest.fixture(scope="module")
def data():
    return harness.BenchData.load_default()


@pytest.fixture(scope="module")
def small_idx(data):
    """~20 stratified rows: 4 actives + 16 decoys."""
    rng = np.random.default_rng(0)
    active_idx = np.where(data.y == 1)[0]
    decoy_idx = np.where(data.y == 0)[0]
    chosen = np.concatenate([
        rng.choice(active_idx, size=4, replace=False),
        rng.choice(decoy_idx, size=16, replace=False),
    ])
    return np.sort(chosen)


def test_scorer_registered_after_discover(data):
    harness.discover()
    assert "shape_combo_rdkit" in harness.REGISTRY, (
        "scorer_shape_combo.py must register 'shape_combo_rdkit' via @register()"
    )


def test_scorer_output_length(data, small_idx):
    harness.discover()
    scorer = harness.REGISTRY["shape_combo_rdkit"]
    scores = scorer(data, np.array([], dtype=int), small_idx)
    assert len(scores) == len(small_idx), "Returned array must have one score per test molecule"


def test_scores_all_finite_and_non_negative(data, small_idx):
    harness.discover()
    scorer = harness.REGISTRY["shape_combo_rdkit"]
    scores = scorer(data, np.array([], dtype=int), small_idx)
    assert np.all(np.isfinite(scores)), "All scores must be finite (no NaN/Inf)"
    assert np.all(scores >= 0.0), "TanimotoCombo scores must be non-negative"


def test_scores_at_most_two(data, small_idx):
    """TanimotoCombo (shape+color) is bounded in [0, 2]."""
    harness.discover()
    scorer = harness.REGISTRY["shape_combo_rdkit"]
    scores = scorer(data, np.array([], dtype=int), small_idx)
    assert np.all(scores <= 2.0 + 1e-6), "TanimotoCombo must be <= 2.0"


def test_refs_cached_on_data(data, small_idx):
    """Calling scorer twice must reuse the cached conformers (same object)."""
    harness.discover()
    scorer = harness.REGISTRY["shape_combo_rdkit"]
    scorer(data, np.array([], dtype=int), small_idx)
    assert "_shape_refs" in data._cache, "References must be cached in data._cache['_shape_refs']"
    refs_first = data._cache["_shape_refs"]
    scorer(data, np.array([], dtype=int), small_idx)
    assert data._cache["_shape_refs"] is refs_first, "Cache must not be rebuilt on second call"

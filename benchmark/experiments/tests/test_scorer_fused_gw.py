"""TDD tests for scorer_fused_gw — Fused Gromov-Wasserstein pose-free scorer.

Tests run on a 30-molecule stratified subsample for speed (full dataset ~575 rows;
FGW O(n^3) per pair × 5 refs × 575 mols would take several minutes in a test suite).
"""
import numpy as np
import pytest

import harness


# ------------------------------------------------------------------ subsample helper

def _stratified_subsample(data, n: int = 30, seed: int = harness.SEED) -> np.ndarray:
    """Return indices for n stratified rows (roughly half active, half decoy)."""
    rng = np.random.default_rng(seed)
    active_idx = np.where(data.y == 1)[0]
    decoy_idx = np.where(data.y == 0)[0]
    n_act = max(1, n // 2)
    n_dec = n - n_act
    chosen = np.concatenate([
        rng.choice(active_idx, min(n_act, len(active_idx)), replace=False),
        rng.choice(decoy_idx, min(n_dec, len(decoy_idx)), replace=False),
    ])
    return chosen


# ------------------------------------------------------------------ registration

def test_fused_gw_registers_after_discover():
    """discover() must find and register the fused_gw scorer."""
    harness.discover()
    assert "fused_gw" in harness.REGISTRY


# ------------------------------------------------------------------ scorer contract

def test_fused_gw_output_length_matches_test_idx():
    """Score array length must equal len(test_idx)."""
    harness.discover()
    data = harness.BenchData.load_default()
    small_idx = _stratified_subsample(data, n=30)
    scores = harness.REGISTRY["fused_gw"](data, train_idx=np.array([]), test_idx=small_idx)
    assert len(scores) == len(small_idx), (
        f"Expected {len(small_idx)} scores, got {len(scores)}"
    )


def test_fused_gw_scores_all_finite():
    """All returned scores must be finite (no NaN/Inf)."""
    harness.discover()
    data = harness.BenchData.load_default()
    small_idx = _stratified_subsample(data, n=30)
    scores = harness.REGISTRY["fused_gw"](data, train_idx=np.array([]), test_idx=small_idx)
    assert np.all(np.isfinite(scores)), (
        f"Non-finite scores found: {scores[~np.isfinite(scores)]}"
    )


def test_fused_gw_scores_in_unit_interval():
    """Similarity = exp(-fgw) ∈ [0, 1]; scores must be in [0, 1]."""
    harness.discover()
    data = harness.BenchData.load_default()
    small_idx = _stratified_subsample(data, n=30)
    scores = harness.REGISTRY["fused_gw"](data, train_idx=np.array([]), test_idx=small_idx)
    assert scores.min() >= 0.0, f"Scores below 0: min={scores.min()}"
    assert scores.max() <= 1.0, f"Scores above 1: max={scores.max()}"


def test_fused_gw_scores_have_some_variance():
    """Non-trivial scorer: must not return a flat constant array."""
    harness.discover()
    data = harness.BenchData.load_default()
    small_idx = _stratified_subsample(data, n=30)
    scores = harness.REGISTRY["fused_gw"](data, train_idx=np.array([]), test_idx=small_idx)
    assert scores.std() > 1e-6, "Scorer returns constant scores — uninformative"


def test_fused_gw_ignores_train_idx():
    """Unsupervised method: scores must not change when train_idx varies."""
    harness.discover()
    data = harness.BenchData.load_default()
    small_idx = _stratified_subsample(data, n=30)
    scores_no_train = harness.REGISTRY["fused_gw"](
        data, train_idx=np.array([]), test_idx=small_idx
    )
    # Pass a different train_idx — should produce identical output
    other_train = np.array([100, 200, 300])
    scores_with_train = harness.REGISTRY["fused_gw"](
        data, train_idx=other_train, test_idx=small_idx
    )
    np.testing.assert_array_equal(
        scores_no_train, scores_with_train,
        err_msg="fused_gw is unsupervised; train_idx must not affect scores"
    )


def test_fused_gw_caches_mol_features():
    """Per-molecule features must be cached in data._cache to avoid recomputation."""
    harness.discover()
    data = harness.BenchData.load_default()
    small_idx = _stratified_subsample(data, n=10)
    # First call
    harness.REGISTRY["fused_gw"](data, train_idx=np.array([]), test_idx=small_idx)
    assert "fgw_mol_features" in data._cache, (
        "Expected 'fgw_mol_features' key in data._cache after scoring"
    )
    # Second call should hit cache (same object)
    cached = data._cache["fgw_mol_features"]
    harness.REGISTRY["fused_gw"](data, train_idx=np.array([]), test_idx=small_idx)
    assert data._cache["fgw_mol_features"] is cached, "Cache replaced on second call"

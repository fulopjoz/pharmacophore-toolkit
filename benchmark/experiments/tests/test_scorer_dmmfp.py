"""Tests for the DMMFP (Differential Multimolecule Fingerprint) scorer.

Reference: Hutter, M. C. (2022). Differential Multimolecule Fingerprints for
Virtual Screening. J. Chem. Inf. Model., 62(9), 2112–2120.
https://doi.org/10.1021/acs.jcim.2c00242
"""
import harness


def test_dmmfp_is_registered():
    """After discover(), the DMMFP scorer must appear in the registry."""
    harness.discover()
    assert "differential_mmfp" in harness.REGISTRY


def test_dmmfp_evaluate_returns_four_bounded_metrics():
    """evaluate() must return all 4 metric keys and AUC must be in [0, 1]."""
    d = harness.BenchData.load_default()
    harness.discover()
    m = harness.evaluate("differential_mmfp", d)
    assert {"AUC", "BEDROC", "EF1%", "EF5%"} <= set(m)
    assert 0.0 <= m["AUC"] <= 1.0


def test_dmmfp_auc_beats_random():
    """DMMFP is a strong discriminator — AUC must be well above 0.7."""
    d = harness.BenchData.load_default()
    harness.discover()
    m = harness.evaluate("differential_mmfp", d)
    assert m["AUC"] > 0.7, f"Expected AUC > 0.7, got {m['AUC']:.4f}"

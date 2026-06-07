"""TDD for templates.py — leakage-safe Butina template selection."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from templates import cluster_templates  # noqa: E402


def test_returns_at_most_max_templates():
    smis = ["c1ccccc1" + "C" * i for i in range(12)] + \
           ["c1ccncc1", "c1ccc2ccccc2c1", "c1ccoc1"]
    tmpl = cluster_templates(smis, sim_cutoff=0.65, max_templates=5, seed=42)
    assert 1 <= len(tmpl) <= 5
    assert all(isinstance(t, str) for t in tmpl)


def test_deterministic():
    smis = ["c1ccccc1C", "c1ccncc1", "c1ccc2ccccc2c1", "CCCCN", "c1ccoc1"]
    a = cluster_templates(smis, sim_cutoff=0.65, max_templates=8, seed=42)
    b = cluster_templates(smis, sim_cutoff=0.65, max_templates=8, seed=42)
    assert a == b


def test_distinct_scaffolds_give_multiple_templates():
    smis = ["c1ccncc1", "c1ccc2ccccc2c1", "c1ccoc1", "C1CCCCC1"]
    tmpl = cluster_templates(smis, sim_cutoff=0.65, max_templates=8, seed=42)
    assert len(tmpl) >= 3   # chemically distinct -> separate clusters


def test_empty_input_returns_empty():
    assert cluster_templates([], sim_cutoff=0.65, max_templates=8, seed=42) == []

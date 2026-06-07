"""TDD for color3d.py — per-feature-type 3D color overlap + pose decomposition."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from color3d import per_type_overlap, embed, feature_points, align_decompose, charge_points, esp_overlap  # noqa: E402
from harness import P4  # noqa: E402  ["donor","acceptor","anion","cation","hydrophobe","rings"]


def _pts(**kw):
    """Build a per-type points dict; each value an (n,3) float array."""
    return {t: np.asarray(kw.get(t, np.zeros((0, 3))), float).reshape(-1, 3) for t in P4}


def test_identical_pointsets_give_tanimoto_one_per_type():
    ref = _pts(donor=[[0, 0, 0]], hydrophobe=[[1, 1, 1], [2, 2, 2]])
    out = per_type_overlap(ref, ref, alpha=0.5)
    assert out.shape == (len(P4),)
    assert abs(out[P4.index("donor")] - 1.0) < 1e-9
    assert abs(out[P4.index("hydrophobe")] - 1.0) < 1e-9


def test_far_apart_points_give_near_zero():
    ref = _pts(donor=[[0, 0, 0]])
    qry = _pts(donor=[[50, 50, 50]])
    out = per_type_overlap(ref, qry, alpha=0.5)
    assert out[P4.index("donor")] < 1e-6


def test_type_absent_in_one_side_is_zero():
    ref = _pts(donor=[[0, 0, 0]])
    qry = _pts(cation=[[0, 0, 0]])
    out = per_type_overlap(ref, qry, alpha=0.5)
    assert out[P4.index("donor")] == 0.0
    assert out[P4.index("cation")] == 0.0


def test_overlap_is_symmetric_and_in_unit_range():
    ref = _pts(rings=[[0, 0, 0], [3, 0, 0]])
    qry = _pts(rings=[[0.5, 0, 0]])
    a = per_type_overlap(ref, qry, alpha=0.5)[P4.index("rings")]
    b = per_type_overlap(qry, ref, alpha=0.5)[P4.index("rings")]
    assert abs(a - b) < 1e-9
    assert 0.0 <= a <= 1.0


def test_embed_and_feature_points():
    m = embed("c1ccccc1CCN", nconf=2, seed=42)
    assert m is not None and m.GetNumConformers() >= 1
    pts = feature_points(m, m.GetConformer(0).GetId())
    # phenethylamine: an aromatic ring and a donor/cation amine
    assert pts["rings"].shape[0] >= 1
    assert pts["donor"].shape[0] >= 1


def test_align_decompose_self_is_high_color():
    m = embed("c1ccc(CCN)cc1O", nconf=3, seed=42)
    vec = align_decompose(m, m, alpha=0.5)   # align a molecule to itself
    assert vec.shape == (6,)
    present = vec[vec > 0]
    assert present.size >= 1 and present.max() > 0.5


def test_align_decompose_returns_zero_vector_on_bad_input():
    assert np.allclose(align_decompose(None, None, alpha=0.5), 0.0)


def test_esp_overlap_identical_is_one():
    xyz = np.array([[0,0,0],[1.5,0,0]], float); q = np.array([0.3,-0.3], float)
    assert abs(esp_overlap(xyz, q, xyz, q, alpha=0.3) - 1.0) < 1e-9

def test_esp_overlap_sign_flip_is_minus_one():
    xyz = np.array([[0,0,0]], float)
    assert abs(esp_overlap(xyz, np.array([0.5]), xyz, np.array([-0.5]), alpha=0.3) + 1.0) < 1e-9

def test_esp_overlap_far_apart_near_zero():
    assert abs(esp_overlap(np.array([[0,0,0]]), np.array([0.5]),
                           np.array([[50,0,0]]), np.array([0.5]), alpha=0.3)) < 1e-6

def test_charge_points_returns_coords_and_charges():
    m = embed("c1ccccc1O", nconf=1, seed=42)
    xyz, q = charge_points(m, m.GetConformer(0).GetId())
    assert xyz.shape[0] == q.shape[0] >= 1 and np.isfinite(q).all()

def test_align_decompose_with_esp_has_seven_columns():
    m = embed("c1ccc(CCN)cc1O", nconf=2, seed=42)
    vec = align_decompose(m, m, alpha=0.5, with_esp=True)
    assert vec.shape == (7,) and np.isfinite(vec).all()

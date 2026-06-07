"""TDD for prism_core: leakage-safe templates + per-(template×type) feature matrix."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from harness import BenchData, P4
from prism_core import make_templates, feature_matrix, prism_features


def _data():
    act = ["c1ccc(CCN)cc1O", "c1ccc(CCN)cc1", "c1ccc2ccccc2c1CN"]
    dec = ["CCCCCCN", "CCCCCCO", "CCCCCCCC"]
    ref = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                       "tutorials", "data", "CCR2_reference_ligands.sdf"))
    return BenchData.from_lists(act, dec, ref)


def test_make_templates_from_train_actives():
    tmpls = make_templates(["c1ccc(CCN)cc1O", "c1ccc(CCN)cc1"], nconf=2, max_templates=4)
    assert len(tmpls) >= 1
    assert all(m.GetNumConformers() >= 1 for m in tmpls)


def test_feature_matrix_color_only_shape():
    data = _data()
    tmpls = make_templates(["c1ccc(CCN)cc1O"], nconf=2, max_templates=4)
    X = feature_matrix(data, tmpls, with_esp=False, nconf=2, njobs=1)
    assert X.shape == (len(data.smiles), len(tmpls) * len(P4))   # K*6
    assert np.all(np.isfinite(X))


def test_feature_matrix_with_esp_is_seven_wide():
    """Regression guard: with_esp=True must yield K*7 columns (6 color + 1 ESP),
    so an accidental with_esp drop can't silently degrade prism_esp to non-ESP."""
    data = _data()
    tmpls = make_templates(["c1ccc(CCN)cc1O"], nconf=2, max_templates=4)
    X = feature_matrix(data, tmpls, with_esp=True, nconf=2, njobs=1)
    assert X.shape == (len(data.smiles), len(tmpls) * (len(P4) + 1))   # K*7
    assert np.all(np.isfinite(X))


def test_prism_features_caches_per_split_and_with_esp(monkeypatch):
    """Same data+split+with_esp -> reuse the SAME matrix (prism & prism_fixed share it,
    guaranteeing an identical-feature ablation); a different split rebuilds."""
    import prism_core
    monkeypatch.setattr(prism_core, "NCONF", 2)
    data = _data()
    tr = np.array([0, 1, 3, 4])
    X1, _ = prism_core.prism_features(data, tr, with_esp=False)
    X2, _ = prism_core.prism_features(data, tr, with_esp=False)
    assert X1 is X2                                   # cache hit -> identical object (shared by prism & prism_fixed)
    Xesp, _ = prism_core.prism_features(data, tr, with_esp=True)  # esp is a separate cache entry
    assert Xesp is not X1 and Xesp.shape[1] > X1.shape[1]   # K*7 > K*6
    X3, _ = prism_core.prism_features(data, np.array([1, 2, 4, 5]), with_esp=False)
    assert X3 is not X1                               # different split -> rebuild

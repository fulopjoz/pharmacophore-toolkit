"""TDD for prism_core: leakage-safe templates + per-(template×type) feature matrix."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from harness import BenchData, P4
from prism_core import make_templates, feature_matrix


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

# Bake-off 3D Additions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add three transferable 3D comparators to the bake-off — `s3_3d_fixed` (fixed-weight ablation of `s3_3d`), `s3_3d_esp` (adds an electrostatic-similarity channel), and (deferred/optional) a ROSHAMBO engine scorer — so we can tell whether `s3_3d`'s win comes from the learned weighting, whether electrostatics adds signal, and whether a different engine helps.

**Architecture:** Refactor `scorer_s3_3d.py`'s template-build + feature-matrix into a shared `s3_3d_core.py` (DRY); `s3_3d_fixed` reuses the same `(K×6)` per-type color features but aggregates them with equal fixed weights instead of a logistic; `s3_3d_esp` extends `color3d.align_decompose` to also emit a per-template Carbó electrostatic-similarity scalar (RDKit Gasteiger charges, same pose), giving `(K×7)` features for the logistic. ROSHAMBO is an isolated-venv engine scorer added last.

**Tech Stack:** RDKit (`rdShapeAlign`, `ComputeGasteigerCharges`, ETKDGv3, feature factory), scikit-learn (LogisticRegression, StandardScaler), numpy, joblib; the existing `benchmark/experiments/{harness,color3d,templates,scorer_s3_3d}.py`.

**Process guardrails (user-mandated):** every scorer is TDD'd; after implementation there is a **front-node smoke test before any PBS run**, then **two independent code reviews via parallel agents**, fixes applied, and only then the full PBS run.

---

### Task 1: Extract shared core (DRY) without changing `s3_3d` behavior

**Files:**
- Create: `benchmark/experiments/s3_3d_core.py`
- Modify: `benchmark/experiments/scorer_s3_3d.py`
- Test: `benchmark/experiments/tests/test_s3_3d_core.py`

- [ ] **Step 1: Write the failing test**

```python
# benchmark/experiments/tests/test_s3_3d_core.py
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from harness import BenchData, P4
from s3_3d_core import make_templates, feature_matrix


def _data():
    act = ["c1ccc(CCN)cc1O", "c1ccc(CCN)cc1", "c1ccc2ccccc2c1CN"]
    dec = ["CCCCCCN", "CCCCCCO", "CCCCCCCC"]
    ref = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                       "tutorials", "data", "CCR2_reference_ligands.sdf"))
    return BenchData.from_lists(act, dec, ref)


def test_make_templates_from_train_actives():
    data = _data()
    tmpls = make_templates(["c1ccc(CCN)cc1O", "c1ccc(CCN)cc1"], nconf=2, max_templates=4)
    assert len(tmpls) >= 1
    assert all(m.GetNumConformers() >= 1 for m in tmpls)


def test_feature_matrix_shape_color_only():
    data = _data()
    tmpls = make_templates(["c1ccc(CCN)cc1O"], nconf=2, max_templates=4)
    X = feature_matrix(data, tmpls, with_esp=False, nconf=2, njobs=1)
    assert X.shape == (len(data.smiles), len(tmpls) * len(P4))   # K*6
    assert np.all(np.isfinite(X))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=2 ~/miniconda3/bin/python -m pytest tests/test_s3_3d_core.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 's3_3d_core'`.

- [ ] **Step 3: Write minimal implementation**

```python
# benchmark/experiments/s3_3d_core.py
"""Shared core for the s3_3d family: leakage-safe templates + the per-(template×type)
3D-color feature matrix (optionally with an electrostatic-similarity column per template).
Factored out so s3_3d / s3_3d_fixed / s3_3d_esp share one tested feature builder (DRY)."""
from __future__ import annotations

import os
import sys

import numpy as np
from joblib import Parallel, delayed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import BenchData, SEED, P4
from templates import cluster_templates
from color3d import embed, align_decompose


def make_templates(train_act_smiles, nconf: int, max_templates: int = 8, seed: int = SEED):
    """Cluster TRAIN actives -> embedded template mols (leakage-safe)."""
    smis = cluster_templates(train_act_smiles, sim_cutoff=0.65,
                             max_templates=max_templates, seed=seed)
    return [m for m in (embed(s, max(1, nconf // 2), seed) for s in smis) if m]


def _row(smi, tmpl_mols, with_esp, nconf, alpha):
    q = embed(smi, nconf, SEED)
    width = len(P4) + (1 if with_esp else 0)
    if q is None:
        return np.zeros(len(tmpl_mols) * width, dtype=float)
    return np.concatenate([align_decompose(q, t, alpha, with_esp=with_esp) for t in tmpl_mols])


def feature_matrix(data: BenchData, tmpl_mols, with_esp: bool, nconf: int,
                   njobs: int, alpha: float = 0.5) -> np.ndarray:
    """(n_mols, K*(6 or 7)) per-(template×type) color (+esp) overlap matrix."""
    rows = Parallel(n_jobs=njobs, backend="loky")(
        delayed(_row)(data.smiles[i], tmpl_mols, with_esp, nconf, alpha)
        for i in range(len(data.smiles)))
    return np.vstack(rows)
```

Then refactor `scorer_s3_3d.py` to use the core (delete its local `_row`; keep behavior identical):

```python
# benchmark/experiments/scorer_s3_3d.py  (replace body after the imports/env block)
from s3_3d_core import make_templates, feature_matrix  # noqa: E402

@register("s3_3d")
def s3_3d(data: BenchData, train_idx, test_idx) -> np.ndarray:
    tr = np.asarray(train_idx)
    train_act = [data.smiles[i] for i in tr if data.y[i] == 1]
    tmpl_mols = make_templates(train_act, nconf=_NCONF, max_templates=_MAX_TEMPLATES, seed=SEED)
    if not tmpl_mols:
        return np.zeros(len(test_idx), dtype=float)
    X = feature_matrix(data, tmpl_mols, with_esp=False, nconf=_NCONF, njobs=_NJOBS, alpha=_ALPHA)
    sc = StandardScaler().fit(X[tr])
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
    lr.fit(sc.transform(X[tr]), data.y[tr])
    return lr.predict_proba(sc.transform(X[np.asarray(test_idx)]))[:, 1]
```
(Keep the existing `_NCONF/_NJOBS/_MAX_TEMPLATES/_ALPHA` env constants and imports of `register, BenchData, SEED, StandardScaler, LogisticRegression, np` at the top of `scorer_s3_3d.py`.)

- [ ] **Step 4: Run tests (core + existing s3_3d regression)**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 ~/miniconda3/bin/python -m pytest tests/test_s3_3d_core.py tests/test_scorer_s3_3d.py -q`
Expected: PASS (3 passed) — `s3_3d` still registers and scores aligned length.

- [ ] **Step 5: Commit**

```bash
git add benchmark/experiments/s3_3d_core.py benchmark/experiments/scorer_s3_3d.py benchmark/experiments/tests/test_s3_3d_core.py
git commit -m "refactor(benchmark): extract s3_3d_core (templates + feature matrix), DRY"
```

---

### Task 2: `s3_3d_fixed` — fixed-weight ablation

**Files:**
- Create: `benchmark/experiments/scorer_s3_3d_fixed.py`
- Test: `benchmark/experiments/tests/test_scorer_s3_3d_fixed.py`

- [ ] **Step 1: Write the failing test**

```python
# benchmark/experiments/tests/test_scorer_s3_3d_fixed.py
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from harness import BenchData, discover, REGISTRY


def test_fixed_registers_runs_and_differs_from_learned():
    discover()
    assert "s3_3d_fixed" in REGISTRY
    act = ["c1ccc(CCN)cc1O", "c1ccc(CCN)cc1", "c1ccc2ccccc2c1CN"]
    dec = ["CCCCCCN", "CCCCCCO", "CCCCCCCC"]
    ref = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                       "tutorials", "data", "CCR2_reference_ligands.sdf"))
    data = BenchData.from_lists(act, dec, ref)
    tr, te = np.array([0, 1, 3, 4]), np.array([2, 5])
    s_fixed = REGISTRY["s3_3d_fixed"](data, tr, te)
    assert s_fixed.shape == (2,) and np.all(np.isfinite(s_fixed))
    # fixed (equal-weight mean) must not be identical to the learned logistic
    s_learned = REGISTRY["s3_3d"](data, tr, te)
    assert not np.allclose(s_fixed, s_learned)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 ~/miniconda3/bin/python -m pytest tests/test_scorer_s3_3d_fixed.py -q`
Expected: FAIL with `assert 's3_3d_fixed' in REGISTRY`.

- [ ] **Step 3: Write minimal implementation**

```python
# benchmark/experiments/scorer_s3_3d_fixed.py
"""s3_3d_fixed — the SHAFTS-style fixed-weight ablation of s3_3d.

IDENTICAL per-(template×type) 3D-color features as s3_3d, but aggregated with EQUAL
FIXED weights (mean) instead of a learned logistic. Isolates the single variable
'learned vs fixed weighting' so we can attribute s3_3d's enrichment to the
discrimination weighting rather than the 3D features. UNSUPERVISED given templates
(templates still come from TRAIN actives -> leakage-safe). Honest label: this is a
fixed-weight ablation of OUR features, not literal SHAFTS (not a Python lib)."""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData, SEED
from s3_3d_core import make_templates, feature_matrix

_NCONF = int(os.environ.get("BAKEOFF_NCONF", "5"))
_NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))
_MAX_TEMPLATES = int(os.environ.get("BAKEOFF_MAX_TEMPLATES", "8"))
_ALPHA = float(os.environ.get("BAKEOFF_COLOR_ALPHA", "0.5"))


@register("s3_3d_fixed")
def s3_3d_fixed(data: BenchData, train_idx, test_idx) -> np.ndarray:
    tr = np.asarray(train_idx)
    train_act = [data.smiles[i] for i in tr if data.y[i] == 1]
    tmpl_mols = make_templates(train_act, nconf=_NCONF, max_templates=_MAX_TEMPLATES, seed=SEED)
    if not tmpl_mols:
        return np.zeros(len(test_idx), dtype=float)
    X = feature_matrix(data, tmpl_mols, with_esp=False, nconf=_NCONF, njobs=_NJOBS, alpha=_ALPHA)
    return X[np.asarray(test_idx)].mean(axis=1)   # equal fixed weights
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 ~/miniconda3/bin/python -m pytest tests/test_scorer_s3_3d_fixed.py -q`
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
git add benchmark/experiments/scorer_s3_3d_fixed.py benchmark/experiments/tests/test_scorer_s3_3d_fixed.py
git commit -m "feat(benchmark): s3_3d_fixed — fixed-weight ablation (isolates learned weighting)"
```

---

### Task 3: Electrostatic-similarity overlap in `color3d`

**Files:**
- Modify: `benchmark/experiments/color3d.py`
- Test: `benchmark/experiments/tests/test_color3d.py`

- [ ] **Step 1: Write the failing test** (append to `tests/test_color3d.py`)

```python
from color3d import charge_points, esp_overlap   # add to imports at top of file

def test_esp_overlap_identical_is_one():
    xyz = np.array([[0, 0, 0], [1.5, 0, 0]], float)
    q = np.array([0.3, -0.3], float)
    assert abs(esp_overlap(xyz, q, xyz, q, alpha=0.3) - 1.0) < 1e-9

def test_esp_overlap_sign_flip_is_minus_one():
    xyz = np.array([[0, 0, 0]], float)
    a = esp_overlap(xyz, np.array([0.5]), xyz, np.array([-0.5]), alpha=0.3)
    assert abs(a + 1.0) < 1e-9          # opposite charge, same place -> Carbo = -1

def test_esp_overlap_far_apart_near_zero():
    a = esp_overlap(np.array([[0, 0, 0]]), np.array([0.5]),
                    np.array([[50, 0, 0]]), np.array([0.5]), alpha=0.3)
    assert abs(a) < 1e-6

def test_charge_points_returns_coords_and_charges():
    m = embed("c1ccccc1O", nconf=1, seed=42)
    xyz, q = charge_points(m, m.GetConformer(0).GetId())
    assert xyz.shape[0] == q.shape[0] >= 1
    assert np.isfinite(q).all()

def test_align_decompose_with_esp_has_seven_columns():
    m = embed("c1ccc(CCN)cc1O", nconf=2, seed=42)
    vec = align_decompose(m, m, alpha=0.5, with_esp=True)
    assert vec.shape == (7,)             # 6 color types + 1 esp scalar
    assert np.isfinite(vec).all()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd benchmark/experiments && ~/miniconda3/bin/python -m pytest tests/test_color3d.py -q`
Expected: FAIL with `ImportError: cannot import name 'charge_points'`.

- [ ] **Step 3: Write minimal implementation** (append to `color3d.py`; and extend `align_decompose` signature)

```python
def charge_points(mol, conf_id: int):
    """Return (coords (n,3), gasteiger_charges (n,)) for the given conformer."""
    if not mol.GetAtomWithIdx(0).HasProp("_GasteigerCharge"):
        AllChem.ComputeGasteigerCharges(mol)
    conf = mol.GetConformer(conf_id)
    xyz, q = [], []
    for a in mol.GetAtoms():
        c = a.GetDoubleProp("_GasteigerCharge")
        if not np.isfinite(c):
            c = 0.0
        p = conf.GetAtomPosition(a.GetIdx())
        xyz.append([p.x, p.y, p.z]); q.append(c)
    return np.asarray(xyz, float), np.asarray(q, float)


def esp_overlap(A_xyz, A_q, B_xyz, B_q, alpha: float = 0.3) -> float:
    """Carbo electrostatic-similarity index in [-1,1] (same-sign charges that coincide
    score high). O(A,B)=sum_ij q_i q_j exp(-alpha|a_i-b_j|^2); C=O_AB/sqrt(O_AA*O_BB)."""
    o_ab = _q_gauss(A_xyz, A_q, B_xyz, B_q, alpha)
    o_aa = _q_gauss(A_xyz, A_q, A_xyz, A_q, alpha)
    o_bb = _q_gauss(B_xyz, B_q, B_xyz, B_q, alpha)
    denom = (o_aa * o_bb) ** 0.5
    return float(o_ab / denom) if denom > 1e-12 else 0.0


def _q_gauss(A_xyz, A_q, B_xyz, B_q, alpha: float) -> float:
    if len(A_xyz) == 0 or len(B_xyz) == 0:
        return 0.0
    d2 = ((A_xyz[:, None, :] - B_xyz[None, :, :]) ** 2).sum(-1)
    return float((np.outer(A_q, B_q) * np.exp(-alpha * d2)).sum())
```

Extend `align_decompose` to optionally append the ESP scalar from the SAME best pose:

```python
def align_decompose(query_mol, template_mol, alpha: float = 0.5, with_esp: bool = False) -> np.ndarray:
    width = len(P4) + (1 if with_esp else 0)
    if query_mol is None or template_mol is None:
        return np.zeros(width, dtype=float)
    if query_mol.GetNumConformers() == 0 or template_mol.GetNumConformers() == 0:
        return np.zeros(width, dtype=float)
    ref_conf = template_mol.GetConformer(0).GetId()
    ref_pts = feature_points(template_mol, ref_conf)
    ref_cxyz, ref_cq = charge_points(template_mol, ref_conf) if with_esp else (None, None)
    best, best_sum = np.zeros(width, dtype=float), -1.0
    for qc in query_mol.GetConformers():
        probe = Chem.Mol(query_mol)
        try:
            rdShapeAlign.AlignMol(template_mol, probe, refConfId=ref_conf,
                                  probeConfId=qc.GetId(), useColors=True, opt_param=0.5)
        except (RuntimeError, ValueError):
            continue
        vec = per_type_overlap(ref_pts, feature_points(probe, qc.GetId()), alpha)
        if with_esp:
            q_cxyz, q_cq = charge_points(probe, qc.GetId())
            vec = np.append(vec, esp_overlap(ref_cxyz, ref_cq, q_cxyz, q_cq, alpha=0.3))
        if vec[:len(P4)].sum() > best_sum:           # choose pose by COLOR overlap (consistent with s3_3d)
            best, best_sum = vec, vec[:len(P4)].sum()
    return best
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd benchmark/experiments && ~/miniconda3/bin/python -m pytest tests/test_color3d.py -q`
Expected: PASS (12 passed — 7 prior + 5 new). The prior `align_decompose` tests still pass because `with_esp` defaults to False.

- [ ] **Step 5: Commit**

```bash
git add benchmark/experiments/color3d.py benchmark/experiments/tests/test_color3d.py
git commit -m "feat(benchmark): Carbo electrostatic-similarity overlap + with_esp decomposition"
```

---

### Task 4: `s3_3d_esp` scorer

**Files:**
- Create: `benchmark/experiments/scorer_s3_3d_esp.py`
- Test: `benchmark/experiments/tests/test_scorer_s3_3d_esp.py`

- [ ] **Step 1: Write the failing test**

```python
# benchmark/experiments/tests/test_scorer_s3_3d_esp.py
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from harness import BenchData, discover, REGISTRY


def test_esp_registers_and_scores_aligned_length():
    discover()
    assert "s3_3d_esp" in REGISTRY
    act = ["c1ccc(CCN)cc1O", "c1ccc(CCN)cc1", "c1ccc2ccccc2c1CN"]
    dec = ["CCCCCCN", "CCCCCCO", "CCCCCCCC"]
    ref = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                       "tutorials", "data", "CCR2_reference_ligands.sdf"))
    data = BenchData.from_lists(act, dec, ref)
    tr, te = np.array([0, 1, 3, 4]), np.array([2, 5])
    s = REGISTRY["s3_3d_esp"](data, tr, te)
    assert s.shape == (2,) and np.all(np.isfinite(s))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 ~/miniconda3/bin/python -m pytest tests/test_scorer_s3_3d_esp.py -q`
Expected: FAIL with `assert 's3_3d_esp' in REGISTRY`.

- [ ] **Step 3: Write minimal implementation**

```python
# benchmark/experiments/scorer_s3_3d_esp.py
"""s3_3d_esp — s3_3d plus an electrostatic-similarity channel.

Same templates + per-(template×type) color features as s3_3d, but each template also
contributes a Carbo electrostatic-similarity scalar (RDKit Gasteiger charges, computed
on the same shape+color-optimized pose) -> (K×7) features -> logistic. Tests whether
electrostatics adds discrimination beyond shape+color. Honest scope: an ESP-similarity
FEATURE on a color-optimized pose, not field-based alignment (ShaEP/EON)."""
from __future__ import annotations

import os
import sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData, SEED
from s3_3d_core import make_templates, feature_matrix

_NCONF = int(os.environ.get("BAKEOFF_NCONF", "5"))
_NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))
_MAX_TEMPLATES = int(os.environ.get("BAKEOFF_MAX_TEMPLATES", "8"))
_ALPHA = float(os.environ.get("BAKEOFF_COLOR_ALPHA", "0.5"))


@register("s3_3d_esp")
def s3_3d_esp(data: BenchData, train_idx, test_idx) -> np.ndarray:
    tr = np.asarray(train_idx)
    train_act = [data.smiles[i] for i in tr if data.y[i] == 1]
    tmpl_mols = make_templates(train_act, nconf=_NCONF, max_templates=_MAX_TEMPLATES, seed=SEED)
    if not tmpl_mols:
        return np.zeros(len(test_idx), dtype=float)
    X = feature_matrix(data, tmpl_mols, with_esp=True, nconf=_NCONF, njobs=_NJOBS, alpha=_ALPHA)
    sc = StandardScaler().fit(X[tr])
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
    lr.fit(sc.transform(X[tr]), data.y[tr])
    return lr.predict_proba(sc.transform(X[np.asarray(test_idx)]))[:, 1]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 ~/miniconda3/bin/python -m pytest tests/test_scorer_s3_3d_esp.py -q`
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
git add benchmark/experiments/scorer_s3_3d_esp.py benchmark/experiments/tests/test_scorer_s3_3d_esp.py
git commit -m "feat(benchmark): s3_3d_esp — adds Carbo electrostatic-similarity channel"
```

---

### Task 5: Wire into the driver + run the full new-scorer test suite + front-node smoke

**Files:**
- Modify: `benchmark/experiments/bakeoff.py` (`SLOW_METHODS`)

- [ ] **Step 1: Add both to `SLOW_METHODS`**

Change:
```python
SLOW_METHODS = ["shape_combo_rdkit", "rdshape_ensemble", "learned_scorer", "s3_3d"]
```
to:
```python
SLOW_METHODS = ["shape_combo_rdkit", "rdshape_ensemble", "learned_scorer",
                "s3_3d", "s3_3d_fixed", "s3_3d_esp"]
```

- [ ] **Step 2: Full new-scorer test suite (must be green before any run)**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 ~/miniconda3/bin/python -m pytest tests/test_color3d.py tests/test_templates.py tests/test_s3_3d_core.py tests/test_scorer_s3_3d.py tests/test_scorer_s3_3d_fixed.py tests/test_scorer_s3_3d_esp.py -q`
Expected: PASS (all green).

- [ ] **Step 3: Front-node smoke (code test only) — the user-mandated "test before the run"**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 BAKEOFF_NJOBS=4 ~/miniconda3/bin/python bakeoff.py --methods s3_3d,s3_3d_fixed,s3_3d_esp --datasets ccr2_project --bootstrap 200`
Expected: three BEDROC/AUC lines on `ccr2_project` (`test 19a/127d`), finite, no crash; `results/BAKEOFF.md` written.

- [ ] **Step 4: Commit**

```bash
git add benchmark/experiments/bakeoff.py
git commit -m "feat(benchmark): include s3_3d_fixed + s3_3d_esp in the bake-off"
```

---

### Task 6: Two independent code reviews (user-mandated, parallel agents)

**Files:** none (review gate before the PBS run)

- [ ] **Step 1: Dispatch TWO code-review agents in parallel (one message, two Agent tool calls)** over the working-tree diff of `color3d.py`, `s3_3d_core.py`, `scorer_s3_3d.py`, `scorer_s3_3d_fixed.py`, `scorer_s3_3d_esp.py`. Reviewer A focuses on **correctness** (per-type/ESP overlap math, pose selection, leakage via train-only templates, feature-matrix alignment to `data.smiles`, NaN/empty guards). Reviewer B focuses on **scientific validity** (is the fixed ablation a fair control? is the ESP overlap defined correctly as Carbó same-sign similarity? overfitting from extra columns? determinism). Each returns a JSON list of {file, line, severity, issue}.
- [ ] **Step 2:** Consolidate, fix every CONFIRMED critical/high finding, re-run the Task 5 Step 2 test suite (must stay green), commit fixes.
- [ ] **Step 3 (second review pass):** Re-dispatch one code-review agent on the *post-fix* diff to confirm the fixes are correct and introduced no regressions. This satisfies "review the code at least two times." Commit any further fixes.

```bash
git add -A && git commit -m "fix(benchmark): address two-pass code review of s3_3d_fixed/esp"
```

---

### Task 7: Full PBS run (only after Tasks 5–6 pass)

**Files:**
- Create: `ccr2_gen/scripts/62_bakeoff_3d_additions.sh` (DrugEx repo) — copy `61_s3_3d_bakeoff.sh`, change `#PBS -N bake3d_add`, log name, and rely on the default full method set (now 8 methods incl. fixed+esp). Keep `ncpus=46`, `/scratch` TMPDIR, `NJOBS=44`, `NCONF=15`, `MAX_TEMPLATES=8`.

- [ ] **Step 1:** Write `62_bakeoff_3d_additions.sh`; `bash -n` it.
- [ ] **Step 2:** Submit `/opt/pbs/bin/qsub ccr2_gen/scripts/62_bakeoff_3d_additions.sh`; record job id; background-monitor `qstat` until done.
- [ ] **Step 3:** Read `benchmark/results/BAKEOFF.md`. **Decision gates:** (a) `s3_3d` vs `s3_3d_fixed` on `created_CCR2`/`ccr2_mubd` — if `s3_3d` > `s3_3d_fixed` with non-overlapping CIs, the **learned weighting** (not the features) drives the win; (b) `s3_3d_esp` vs `s3_3d` on MUBD — if `+` with non-overlapping CIs, electrostatics adds real signal; if only on property-matched, it's bias.
- [ ] **Step 4:** Update `benchmark/results/S3_3D_VERDICT.md` with the ablation + ESP outcomes; commit results + PBS script.

```bash
git -C /home/fulopj/drx-run/DrugEx add ccr2_gen/scripts/62_bakeoff_3d_additions.sh
git -C /home/fulopj/drx-run/DrugEx commit -m "feat(ccr2): PBS run for s3_3d_fixed + s3_3d_esp bake-off"
git add benchmark/results/BAKEOFF.md benchmark/results/bakeoff_master.tsv benchmark/results/S3_3D_VERDICT.md
git commit -m "results(benchmark): ablation (fixed) + electrostatics (esp) outcomes"
```

---

### Task 8 (DEFERRED / OPTIONAL): ROSHAMBO engine scorer

> Per the critical assessment (`docs/research/TRANSFERABLE_3D_METHODS.md` §3C): engine swap = confounded, predictably collapses unweighted on MUBD, and needs an isolated install. Do ONLY after Tasks 1–7 and only if explicitly prioritized.

- [ ] **Step 1:** Create an isolated venv on a PBS node (NOT the shared conda env): `python -m venv /scratch/$USER/roshambo_env && source .../bin/activate && pip install roshambo` (or from source per its repo); verify `import roshambo` + GPU availability there.
- [ ] **Step 2:** Write `benchmark/experiments/scorer_roshambo.py` mirroring `scorer_rdshape_ensemble.py` (unsupervised: align test mols to the curated references, return combo score), guarded so it returns zeros + logs a clear skip if `roshambo` is unimportable (keeps the bake-off runnable without it).
- [ ] **Step 3:** TDD a smoke test (skipped via `pytest.importorskip("roshambo")`); front-node smoke; two-pass code review; then a PBS run on the venv. Label results as an **engine comparison**, not an ablation.

---

## Self-Review

**Spec coverage:** ✅ `s3_3d_fixed` weighting ablation (Task 2) · ✅ electrostatics channel (Tasks 3–4) · ✅ ROSHAMBO (Task 8, deferred per critical assessment) · ✅ DRY core (Task 1) · ✅ test-before-run (Task 5 Step 3) · ✅ two code reviews via parallel agents (Task 6) · ✅ PBS only after tests+reviews (Task 7).

**Placeholder scan:** none — every code step has runnable code; constants are env-defaulted; ESP `alpha=0.3`, color `alpha=0.5` are concrete.

**Type consistency:** `make_templates(train_act_smiles, nconf, max_templates, seed) -> list[mol]`; `feature_matrix(data, tmpl_mols, with_esp, nconf, njobs, alpha) -> (n, K*(6|7))`; `align_decompose(query, template, alpha, with_esp) -> length-6|7`; `charge_points -> (xyz, q)`; `esp_overlap(A_xyz,A_q,B_xyz,B_q,alpha) -> float`. Signatures are consistent across Tasks 1–4. `align_decompose`'s new `with_esp` defaults False, preserving the Task-3-of-the-prior-plan callers.

**Risk note:** `s3_3d_esp` adds K extra columns (overfitting risk on small actives) — mitigated by L2 + the MUBD held-out gate (Task 7 decision gate (b)); trust ESP only if it lifts the *unbiased* set with non-overlapping CIs.

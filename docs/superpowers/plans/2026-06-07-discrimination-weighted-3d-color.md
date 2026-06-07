# Discrimination-Weighted 3D Color (S3-in-3D) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a bake-off scorer that makes 3D ROCS *color* discriminative by learning per-(feature-type × cluster-template) weights from actives-vs-decoys, and prove on the unbiased benchmark suite whether it beats unweighted 3D color (`rdshape_ensemble`/`shape_combo`) and 2D `s3_weighted`.

**Architecture:** Cluster the *train* actives into K templates (leakage-safe Butina). For each molecule, align it to each template with `rdShapeAlign` (shape+color-optimized pose, `opt_param=0.5`), then compute the **6 per-feature-type Gaussian color overlaps** ourselves from `feat.GetPos()` (rdShapeAlign only returns the aggregate color scalar). That yields a `(K×6)` feature matrix; a logistic regression (the S3 step) learns the discrimination weights `w[t,f]` and predicts P(active). Validated by re-scoring on the committed CCR2 datasets through the existing `bakeoff.py` (held-out scaffold split + bootstrap + Friedman).

**Tech Stack:** RDKit (`rdShapeAlign`, `ChemicalFeatures`, `rdMolDescriptors`, ETKDGv3, `ML.Cluster.Butina`), scikit-learn (`LogisticRegression`, `StandardScaler`), numpy, joblib; the existing `benchmark/experiments/harness.py` plugin registry.

---

## File structure (what each file owns)

- `benchmark/experiments/color3d.py` — **pure geometry + alignment** (no ML, no dataset): per-type feature points, per-type Gaussian Tanimoto overlap, 3D embed, align-and-decompose to a 6-vector. The novel, unit-testable core.
- `benchmark/experiments/templates.py` — **leakage-safe template selection**: Butina cluster a SMILES list → ≤K representative SMILES.
- `benchmark/experiments/scorer_s3_3d.py` — **bake-off adapter**: train-derived templates → `(K×6)` features for all rows → logistic fit-on-train / predict-test. Registers `s3_3d`.
- `benchmark/experiments/tests/test_color3d.py` — TDD for the overlap geometry + decomposition.
- `benchmark/experiments/tests/test_templates.py` — TDD for clustering.
- `benchmark/experiments/bakeoff.py` — **modify**: add `s3_3d` to `SLOW_METHODS`.
- `ccr2_gen/scripts/61_s3_3d_bakeoff.sh` (DrugEx repo) — **modify/create**: PBS run including `s3_3d`.

Shared constants reused from `harness.py`: `P4` (the 6 type names), `_FACTORY` (BaseFeatures fdef), `_FAM2P4` (RDKit family → P4 type), `SEED`.

---

### Task 1: Per-feature-type Gaussian color overlap (the core math)

**Files:**
- Create: `benchmark/experiments/color3d.py`
- Test: `benchmark/experiments/tests/test_color3d.py`

- [ ] **Step 1: Write the failing test**

```python
# benchmark/experiments/tests/test_color3d.py
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from color3d import per_type_overlap
from harness import P4   # ["donor","acceptor","anion","cation","hydrophobe","rings"]


def _pts(**kw):
    """Build a per-type points dict; each value an (n,3) float array."""
    return {t: np.asarray(kw.get(t, np.zeros((0, 3))), float).reshape(-1, 3) for t in P4}


def test_identical_pointsets_give_tanimoto_one_per_type():
    ref = _pts(donor=[[0, 0, 0]], hydrophobe=[[1, 1, 1], [2, 2, 2]])
    out = per_type_overlap(ref, ref, alpha=0.5)
    assert out.shape == (len(P4),)
    assert out[P4.index("donor")] == np.float64(1.0) or abs(out[P4.index("donor")] - 1.0) < 1e-9
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
    assert out[P4.index("cation")] == 0.0  # absent in ref

def test_overlap_is_symmetric_and_in_unit_range():
    ref = _pts(rings=[[0, 0, 0], [3, 0, 0]])
    qry = _pts(rings=[[0.5, 0, 0]])
    a = per_type_overlap(ref, qry, alpha=0.5)[P4.index("rings")]
    b = per_type_overlap(qry, ref, alpha=0.5)[P4.index("rings")]
    assert abs(a - b) < 1e-9
    assert 0.0 <= a <= 1.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd benchmark/experiments && ~/miniconda3/bin/python -m pytest tests/test_color3d.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'color3d'`.

- [ ] **Step 3: Write minimal implementation**

```python
# benchmark/experiments/color3d.py
"""Per-feature-type 3D color overlap + pose decomposition for the S3-in-3D scorer.

rdShapeAlign returns only an AGGREGATE color scalar; to weight color by feature
type we recompute the overlap ourselves from feature centroids (feat.GetPos()).
Gaussian volume overlap per type, Tanimoto-normalised, in [0,1]."""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import P4  # noqa: E402


def per_type_overlap(ref_pts: dict, qry_pts: dict, alpha: float = 0.5) -> np.ndarray:
    """Gaussian Tanimoto overlap per P4 feature type.

    O(A,B) = sum_ij exp(-alpha * |a_i - b_j|^2); Tanimoto = O_AB / (O_AA + O_BB - O_AB).
    Returns a length-6 array (order = P4); 0.0 for any type absent on either side."""
    out = np.zeros(len(P4), dtype=float)
    for k, t in enumerate(P4):
        A = np.asarray(ref_pts.get(t, np.zeros((0, 3))), float).reshape(-1, 3)
        B = np.asarray(qry_pts.get(t, np.zeros((0, 3))), float).reshape(-1, 3)
        if len(A) == 0 or len(B) == 0:
            continue
        o_ab = _gauss(A, B, alpha)
        o_aa = _gauss(A, A, alpha)
        o_bb = _gauss(B, B, alpha)
        denom = o_aa + o_bb - o_ab
        out[k] = float(o_ab / denom) if denom > 1e-12 else 0.0
    return out


def _gauss(A: np.ndarray, B: np.ndarray, alpha: float) -> float:
    d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)   # pairwise squared dist
    return float(np.exp(-alpha * d2).sum())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd benchmark/experiments && ~/miniconda3/bin/python -m pytest tests/test_color3d.py -q`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add benchmark/experiments/color3d.py benchmark/experiments/tests/test_color3d.py
git commit -m "feat(benchmark): per-feature-type 3D Gaussian color overlap (S3-in-3D core)"
```

---

### Task 2: Feature points + embed + align-and-decompose

**Files:**
- Modify: `benchmark/experiments/color3d.py`
- Test: `benchmark/experiments/tests/test_color3d.py`

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_color3d.py
from color3d import embed, feature_points, align_decompose

def test_embed_and_feature_points():
    m = embed("c1ccccc1CCN", nconf=2, seed=42)
    assert m is not None and m.GetNumConformers() >= 1
    pts = feature_points(m, m.GetConformer(0).GetId())
    # aniline-ethyl-benzene has a ring and a donor/cation amine
    assert pts["rings"].shape[0] >= 1
    assert pts["donor"].shape[0] >= 1

def test_align_decompose_self_is_high_color():
    m = embed("c1ccc(CCN)cc1O", nconf=3, seed=42)
    vec = align_decompose(m, m, alpha=0.5)   # align a molecule to itself
    assert vec.shape == (6,)
    # self-alignment should make the present types overlap strongly
    present = vec[vec > 0]
    assert present.size >= 1 and present.max() > 0.5

def test_align_decompose_returns_zero_vector_on_bad_input():
    assert np.allclose(align_decompose(None, None, alpha=0.5), 0.0)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd benchmark/experiments && ~/miniconda3/bin/python -m pytest tests/test_color3d.py -q`
Expected: FAIL with `ImportError: cannot import name 'embed'`.

- [ ] **Step 3: Write minimal implementation** (append to `color3d.py`)

```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeAlign

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import _FACTORY, _FAM2P4  # noqa: E402


def embed(smiles: str, nconf: int = 5, seed: int = 42):
    """Return an AddHs'd mol with up to `nconf` ETKDGv3 conformers, or None."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.numThreads = 1
    AllChem.EmbedMultipleConfs(mol, numConfs=nconf, params=p)
    return mol if mol.GetNumConformers() > 0 else None


def feature_points(mol, conf_id: int) -> dict:
    """Map each P4 type -> (n,3) array of feature centroids for the given conformer."""
    acc = {t: [] for t in P4}
    for f in _FACTORY.GetFeaturesForMol(mol):
        t = _FAM2P4.get(f.GetFamily())
        if t is None:
            continue
        pos = f.GetPos(conf_id)
        acc[t].append([pos.x, pos.y, pos.z])
    return {t: np.asarray(v, float).reshape(-1, 3) for t, v in acc.items()}


def align_decompose(query_mol, template_mol, alpha: float = 0.5) -> np.ndarray:
    """Align query to template (shape+color pose) and return the 6 per-type color overlaps.

    Template uses its first conformer; query is tried over all its conformers and the
    pose with the largest total per-type overlap is kept. rdShapeAlign mutates the probe
    copy's coordinates to the aligned pose, which feature_points() then reads."""
    if query_mol is None or template_mol is None:
        return np.zeros(len(P4), dtype=float)
    if query_mol.GetNumConformers() == 0 or template_mol.GetNumConformers() == 0:
        return np.zeros(len(P4), dtype=float)
    ref_conf = template_mol.GetConformer(0).GetId()
    ref_pts = feature_points(template_mol, ref_conf)
    best, best_sum = np.zeros(len(P4), dtype=float), -1.0
    for qc in query_mol.GetConformers():
        probe = Chem.Mol(query_mol)
        try:
            rdShapeAlign.AlignMol(template_mol, probe, refConfId=ref_conf,
                                  probeConfId=qc.GetId(), useColors=True, opt_param=0.5)
        except (RuntimeError, ValueError):
            continue
        qry_pts = feature_points(probe, qc.GetId())
        vec = per_type_overlap(ref_pts, qry_pts, alpha)
        if vec.sum() > best_sum:
            best, best_sum = vec, vec.sum()
    return best
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd benchmark/experiments && ~/miniconda3/bin/python -m pytest tests/test_color3d.py -q`
Expected: PASS (7 passed).

- [ ] **Step 5: Commit**

```bash
git add benchmark/experiments/color3d.py benchmark/experiments/tests/test_color3d.py
git commit -m "feat(benchmark): 3D embed + align-and-decompose to per-type color vector"
```

---

### Task 3: Leakage-safe template clustering

**Files:**
- Create: `benchmark/experiments/templates.py`
- Test: `benchmark/experiments/tests/test_templates.py`

- [ ] **Step 1: Write the failing test**

```python
# benchmark/experiments/tests/test_templates.py
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from templates import cluster_templates


def test_returns_at_most_max_templates():
    # 12 near-duplicate benzenes + 3 distinct rings -> few clusters
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd benchmark/experiments && ~/miniconda3/bin/python -m pytest tests/test_templates.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'templates'`.

- [ ] **Step 3: Write minimal implementation**

```python
# benchmark/experiments/templates.py
"""Leakage-safe template selection: Butina-cluster a SMILES list, return one
representative per cluster (largest clusters first), capped at `max_templates`.

Used to derive 3D-color query templates from the TRAIN actives only — so the
multi-active query is never built from molecules in the held-out test set."""
from __future__ import annotations

from typing import List

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def cluster_templates(smiles: List[str], sim_cutoff: float = 0.65,
                      max_templates: int = 8, seed: int = 42) -> List[str]:
    """Return up to `max_templates` representative SMILES (one per Butina cluster).

    Clusters merge at Tanimoto >= `sim_cutoff` (Butina distance cutoff = 1 - cutoff).
    Representative = first member of the cluster (Butina lists the centroid first).
    Largest clusters are taken first so the templates cover the most actives."""
    mols, keep = [], []
    for s in smiles:
        m = Chem.MolFromSmiles(str(s))
        if m is not None:
            mols.append(m)
            keep.append(str(s))
    if not mols:
        return []
    fps = [_MORGAN.GetFingerprint(m) for m in mols]

    dists = []
    for i in range(1, len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend(1.0 - s for s in sims)
    clusters = Butina.ClusterData(dists, len(fps), 1.0 - sim_cutoff, isDistData=True)

    clusters = sorted(clusters, key=len, reverse=True)   # largest coverage first
    return [keep[c[0]] for c in clusters[:max_templates]]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd benchmark/experiments && ~/miniconda3/bin/python -m pytest tests/test_templates.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add benchmark/experiments/templates.py benchmark/experiments/tests/test_templates.py
git commit -m "feat(benchmark): leakage-safe Butina template selection from train actives"
```

---

### Task 4: The `s3_3d` bake-off scorer

**Files:**
- Create: `benchmark/experiments/scorer_s3_3d.py`

- [ ] **Step 1: Write the failing test** (a smoke check that the scorer registers + runs end-to-end on a tiny set)

```python
# benchmark/experiments/tests/test_scorer_s3_3d.py
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from harness import BenchData, discover, REGISTRY


def test_s3_3d_registers_and_scores_aligned_length():
    discover()
    assert "s3_3d" in REGISTRY
    # tiny CCR2-like set: a few ring actives, a few chain decoys
    act = ["c1ccc(CCN)cc1O", "c1ccc(CCN)cc1", "c1ccc2ccccc2c1CN"]
    dec = ["CCCCCCN", "CCCCCCO", "CCCCCCCC"]
    ref = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                       "tutorials", "data", "CCR2_reference_ligands.sdf"))
    data = BenchData.from_lists(act, dec, ref)
    train_idx = np.array([0, 1, 3, 4]); test_idx = np.array([2, 5])
    scores = REGISTRY["s3_3d"](data, train_idx, test_idx)
    assert scores.shape == (2,)
    assert np.all(np.isfinite(scores))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 ~/miniconda3/bin/python -m pytest tests/test_scorer_s3_3d.py -q`
Expected: FAIL with `assert 's3_3d' in REGISTRY` (scorer file not yet created).

- [ ] **Step 3: Write minimal implementation**

```python
# benchmark/experiments/scorer_s3_3d.py
"""Discrimination-weighted 3D color (S3-in-3D), feature-type x template.

Per held-out split:
  1. cluster the TRAIN actives -> K templates (leakage-safe), embed + annotate;
  2. for every molecule, align to each template and decompose into 6 per-type
     color overlaps -> (K*6) feature matrix;
  3. logistic regression (StandardScaler + LogisticRegression) fit on train rows
     learns the per-(type x template) discrimination weights and predicts test.

The logistic IS the S3 weighting; the 3D-overlap features make it 3D. Compare on
the benchmark vs rdshape_ensemble (unweighted 3D color) and s3_weighted (2D)."""
from __future__ import annotations

import os
import sys

import numpy as np
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import register, BenchData, SEED, P4  # noqa: E402
from templates import cluster_templates  # noqa: E402
from color3d import embed, align_decompose  # noqa: E402

_NCONF = int(os.environ.get("BAKEOFF_NCONF", "5"))
_NJOBS = int(os.environ.get("BAKEOFF_NJOBS", "1"))
_MAX_TEMPLATES = int(os.environ.get("BAKEOFF_MAX_TEMPLATES", "8"))
_ALPHA = float(os.environ.get("BAKEOFF_COLOR_ALPHA", "0.5"))


def _row(smi, tmpl_mols):
    """(K*6) feature row for one molecule: per-type color overlaps vs each template."""
    q = embed(smi, _NCONF, SEED)
    if q is None:
        return np.zeros(len(tmpl_mols) * len(P4), dtype=float)
    return np.concatenate([align_decompose(q, t, _ALPHA) for t in tmpl_mols])


@register("s3_3d")
def s3_3d(data: BenchData, train_idx, test_idx) -> np.ndarray:
    tr = np.asarray(train_idx)
    train_act = [data.smiles[i] for i in tr if data.y[i] == 1]
    templates = cluster_templates(train_act, sim_cutoff=0.65,
                                  max_templates=_MAX_TEMPLATES, seed=SEED)
    tmpl_mols = [m for m in (embed(s, max(1, _NCONF // 2), SEED) for s in templates) if m]
    if not tmpl_mols:
        return np.zeros(len(test_idx), dtype=float)

    rows = Parallel(n_jobs=_NJOBS, backend="loky")(
        delayed(_row)(data.smiles[i], tmpl_mols) for i in range(len(data.smiles)))
    X = np.vstack(rows)

    sc = StandardScaler().fit(X[tr])
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
    lr.fit(sc.transform(X[tr]), data.y[tr])
    return lr.predict_proba(sc.transform(X[np.asarray(test_idx)]))[:, 1]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 ~/miniconda3/bin/python -m pytest tests/test_scorer_s3_3d.py -q`
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
git add benchmark/experiments/scorer_s3_3d.py benchmark/experiments/tests/test_scorer_s3_3d.py
git commit -m "feat(benchmark): s3_3d scorer — discrimination-weighted 3D color (type x template)"
```

---

### Task 5: Wire into the bake-off driver + front-node smoke

**Files:**
- Modify: `benchmark/experiments/bakeoff.py:~46` (`SLOW_METHODS`)

- [ ] **Step 1: Add `s3_3d` to the slow (3D) method set**

Change `SLOW_METHODS` from:
```python
SLOW_METHODS = ["shape_combo_rdkit", "rdshape_ensemble", "learned_scorer"]
```
to:
```python
SLOW_METHODS = ["shape_combo_rdkit", "rdshape_ensemble", "learned_scorer", "s3_3d"]
```

- [ ] **Step 2: Front-node smoke on the smallest dataset (code test only)**

Run: `cd benchmark/experiments && BAKEOFF_NCONF=3 BAKEOFF_NJOBS=2 ~/miniconda3/bin/python bakeoff.py --methods s3_3d --datasets ccr2_project --bootstrap 200`
Expected: prints a BEDROC/AUC line for `s3_3d` on `ccr2_project` with `test 19a/127d`, writes `results/BAKEOFF.md`. No crash; finite numbers.

- [ ] **Step 3: Commit**

```bash
git add benchmark/experiments/bakeoff.py
git commit -m "feat(benchmark): include s3_3d in the bake-off method set"
```

---

### Task 6: Full PBS validation run + interpretation

**Files:**
- Create: `ccr2_gen/scripts/61_s3_3d_bakeoff.sh` (DrugEx repo) — copy `60_pharmacophore_bakeoff.sh`, set `BAKEOFF_MAX_TEMPLATES=8`, and run the full method set so `s3_3d` is compared head-to-head with `rdshape_ensemble`/`shape_combo`/`s3_weighted` on all 3 CCR2 datasets.

- [ ] **Step 1:** Duplicate `60_pharmacophore_bakeoff.sh` → `61_s3_3d_bakeoff.sh`; add `export BAKEOFF_MAX_TEMPLATES=8`; change `#PBS -N s3_3d_bake`. Keep `ncpus=46`, `/scratch` TMPDIR, `BAKEOFF_NJOBS=44`, `BAKEOFF_NCONF=15`.
- [ ] **Step 2:** `bash -n ccr2_gen/scripts/61_s3_3d_bakeoff.sh` (syntax check).
- [ ] **Step 3:** Submit: `/opt/pbs/bin/qsub ccr2_gen/scripts/61_s3_3d_bakeoff.sh`; record job id.
- [ ] **Step 4:** When finished, read `benchmark/results/BAKEOFF.md`. **Decision gate:** does `s3_3d` beat `rdshape_ensemble`/`shape_combo` (unweighted 3D color) on **ccr2_mubd** (the max-unbiased set where unweighted 3D collapsed to ≈0), and does it match/beat 2D `s3_weighted`? If yes → discrimination-weighting rescues 3D color → port `w[t,f]` into the DrugEx `rocs_*` color path as the production lever. If no → log the negative result; the 2D `s3` weighting does not transfer to 3D overlap features on this target.
- [ ] **Step 5:** Commit the PBS script and the resulting `BAKEOFF.md`/`bakeoff_master.tsv`.

```bash
git -C /home/fulopj/drx-run/DrugEx add ccr2_gen/scripts/61_s3_3d_bakeoff.sh
git -C /home/fulopj/drx-run/DrugEx commit -m "feat(ccr2): PBS run for s3_3d (discrimination-weighted 3D color) bake-off"
git add benchmark/results/BAKEOFF.md benchmark/results/bakeoff_master.tsv
git commit -m "results(benchmark): s3_3d vs unweighted 3D color on CCR2 suite"
```

---

## Self-Review

**Spec coverage:** ✅ feature-type × template weighting (Task 4 logistic over `K×6` columns) · ✅ leakage-safe templates from train actives (Task 3, used in Task 4) · ✅ per-type 3D overlap not available from rdShapeAlign → computed from `feat.GetPos()` (Tasks 1–2) · ✅ rescoring-first validation on the benchmark suite (Tasks 5–6) · ✅ head-to-head vs unweighted 3D + 2D s3 (Task 6 decision gate).

**Placeholder scan:** none — every code step has runnable code; alpha/nconf/templates are concrete env-defaulted constants.

**Type consistency:** `per_type_overlap`/`align_decompose` return length-6 arrays (order `P4`); `feature_points` returns `dict[str,(n,3)]`; `cluster_templates` returns `list[str]`; scorer returns `(len(test_idx),)` — consistent across Tasks 1→4. `embed(smiles, nconf, seed)` signature identical in Tasks 2 and 4.

**Known cost / risk:** building `X` aligns `n_mols × K × nconf` poses — heavy on `created_CCR2` (≈3938 × 8 × 15); mitigated by `joblib` (`BAKEOFF_NJOBS`) and template `nconf//2`. PBS only for the full run; front node uses `--datasets ccr2_project --methods s3_3d` with `NCONF=3`.

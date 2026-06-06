"""
Comparison harness — the deep module.

Simple interface: write a `scorer_<name>.py` that calls `@register("name")` on a
function `score(data: BenchData, train_idx, test_idx) -> np.ndarray` (scores for the
test rows). The harness hides everything else: data loading + lazy featurization,
held-out 5-fold CV, RDKit early-recognition metrics, and bootstrap CIs.

Contract (design-by-contract): a scorer receives a BenchData and two index arrays.
SUPERVISED scorers fit on `train_idx` and predict `test_idx`; UNSUPERVISED scorers
(similarity-to-reference) ignore `train_idx` and just score `test_idx`. One signature,
no special cases.
"""
from __future__ import annotations

import importlib
import os
from dataclasses import dataclass, field
from glob import glob
from typing import Callable, Dict, Optional

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures, rdFingerprintGenerator
from rdkit.ML.Scoring.Scoring import CalcBEDROC, CalcEnrichment
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "..", "tutorials", "data")  # ../../tutorials/data
SEED = 42
P4 = ["donor", "acceptor", "anion", "cation", "hydrophobe", "rings"]
_FAM2P4 = {"Donor": "donor", "Acceptor": "acceptor", "NegIonizable": "anion",
           "PosIonizable": "cation", "Hydrophobe": "hydrophobe",
           "LumpedHydrophobe": "hydrophobe", "Aromatic": "rings"}
_FACTORY = ChemicalFeatures.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


# --------------------------------------------------------------------------- data
@dataclass
class BenchData:
    """Held CCR2 actives+decoys+references with lazily-cached, shared featurizations.

    Centralizing featurization here means every scorer sees identical features
    (no per-scorer drift) and each expensive representation is computed once.
    """
    smiles: list
    y: np.ndarray
    ref_sdf: str
    _cache: dict = field(default_factory=dict)

    @classmethod
    def load_default(cls) -> "BenchData":
        return cls.load(
            os.path.join(DATA, "actives_ccr2_N75.csv"),
            os.path.join(DATA, "decoys_ccr2_N500.csv"),
            os.path.join(DATA, "CCR2_reference_ligands.sdf"),
        )

    @classmethod
    def load(cls, actives_csv, decoys_csv, ref_sdf) -> "BenchData":
        import pandas as pd

        def _smiles(path):
            df = pd.read_csv(path)
            col = next(c for c in df.columns if c.lower() == "smiles")
            out = [s for s in df[col] if Chem.MolFromSmiles(str(s)) is not None]
            return out

        a, d = _smiles(actives_csv), _smiles(decoys_csv)
        return cls(smiles=a + d, y=np.r_[np.ones(len(a)), np.zeros(len(d))], ref_sdf=ref_sdf)

    # -- lazy featurizations (each computed once, then cached) --
    def mols(self):
        return self._memo("mols", lambda: [Chem.MolFromSmiles(s) for s in self.smiles])

    def p4_counts(self) -> np.ndarray:
        return self._memo("p4", lambda: np.array([self._p4(m) for m in self.mols()], float))

    def morgan(self) -> np.ndarray:
        return self._memo("morgan", lambda: np.array(
            [np.frombuffer(_MORGAN.GetFingerprint(m).ToBitString().encode(), "u1") - 48
             for m in self.mols()], np.uint8))

    def ref_p4(self) -> np.ndarray:
        return self._memo("refp4", lambda: np.array(
            [self._p4(m) for m in Chem.SDMolSupplier(self.ref_sdf, removeHs=False) if m], float))

    def p4_sim_to_refs(self) -> np.ndarray:
        """Per-P4-type best-matching-reference agreement in (0,1] — the Tier-1 feature."""
        return self._memo("p4sim", self._p4_sim_to_refs)

    # -- helpers --
    def _memo(self, key, fn):
        if key not in self._cache:
            self._cache[key] = fn()
        return self._cache[key]

    @staticmethod
    def _p4(mol) -> np.ndarray:
        c = {k: 0 for k in P4}
        for f in _FACTORY.GetFeaturesForMol(mol):
            t = _FAM2P4.get(f.GetFamily())
            if t:
                c[t] += 1
        return np.array([c[k] for k in P4], float)

    def _p4_sim_to_refs(self) -> np.ndarray:
        X, refs = self.p4_counts(), self.ref_p4()
        # Scale from the REFERENCES only (the fixed query side) — never from the
        # actives/decoys being ranked — so no test-fold statistics leak into the
        # features under held-out CV.
        scale = np.maximum(1.0, refs.max(axis=0))
        sims = np.zeros((len(X), len(P4)))
        for j in range(len(P4)):
            d = np.abs(X[:, [j]] - refs[:, j][None, :]) / scale[j]
            sims[:, j] = np.exp(-d).max(axis=1)
        return sims


# ----------------------------------------------------------------------- registry
Scorer = Callable[[BenchData, np.ndarray, np.ndarray], np.ndarray]
REGISTRY: Dict[str, Scorer] = {}


def register(name: str):
    def deco(fn: Scorer) -> Scorer:
        REGISTRY[name] = fn
        return fn
    return deco


def discover() -> None:
    """Import every scorer_*.py in this directory so each self-registers (plugin pattern)."""
    import sys
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    for path in sorted(glob(os.path.join(HERE, "scorer_*.py"))):
        importlib.import_module(os.path.splitext(os.path.basename(path))[0])


# ------------------------------------------------------------------------ metrics
def _ranked(y, s):
    order = np.argsort(-np.asarray(s, float))
    return [[int(y[i])] for i in order]


def metrics(y, scores) -> dict:
    r = _ranked(y, scores)
    ef1, ef5 = CalcEnrichment(r, 0, [0.01, 0.05])
    return {"AUC": roc_auc_score(y, scores), "EF1%": ef1, "EF5%": ef5,
            "BEDROC": CalcBEDROC(r, 0, 20.0)}


# ----------------------------------------------------------------------- evaluate
def evaluate_oof(name: str, data: BenchData, n_splits: int = 5, seed: int = SEED) -> np.ndarray:
    """5-fold held-out (out-of-fold) scores for a registered scorer."""
    if name not in REGISTRY:
        discover()
    scorer = REGISTRY[name]
    y = data.y
    oof = np.full(len(y), np.nan)
    skf = StratifiedKFold(n_splits, shuffle=True, random_state=seed)
    for tr, te in skf.split(np.zeros(len(y)), y):
        oof[te] = scorer(data, tr, te)
    if np.isnan(oof).any():                                 # post-condition: every row scored
        raise RuntimeError(f"scorer '{name}' left {int(np.isnan(oof).sum())} rows unscored")
    return oof


def evaluate(name: str, data: BenchData, **kw) -> dict:
    return metrics(data.y, evaluate_oof(name, data, **kw))


def bootstrap_delta(oof_a, oof_b, y, metric: str, n: int = 1000, seed: int = SEED):
    """Median + 95% CI of metric(a) - metric(b) by bootstrap resampling of molecules."""
    rng = np.random.default_rng(seed)
    diffs = []
    y = np.asarray(y)
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if y[idx].sum() < 3 or (y[idx] == 0).sum() < 3:
            continue
        diffs.append(metrics(y[idx], oof_a[idx])[metric] - metrics(y[idx], oof_b[idx])[metric])
    if not diffs:                                            # every resample degenerate (tiny data)
        return float("nan"), float("nan"), float("nan")
    return float(np.median(diffs)), float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))

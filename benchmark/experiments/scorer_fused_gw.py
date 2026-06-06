"""Fused Gromov-Wasserstein (FGW) pose-free molecular similarity scorer.

Reference: Vayer et al., "Optimal Transport for structured data with application on
graphs", arXiv:1811.02834 (2018). Implemented via POT >= 0.9 (ot.gromov.fused_gromov_wasserstein2).

Method overview
---------------
Pose-free similarity between a query molecule and each of the 5 CCR2 reference ligands
using FGW optimal transport over topological (graph) distances + atom-feature matrices:

  1. Build C (n×n) = topological atom-atom distance matrix from RDKit GetDistanceMatrix,
     normalised to [0,1] by its own maximum (heavy atoms only — fast, no 3-D conformers).
  2. Build feat (n×10) = per-atom features:
       - one-hot of element in {C,N,O,F,S,Cl,other}  (7 dims)
       - is_aromatic                                   (1 dim)
       - is_H-bond_donor  (N or O with ≥1 H)          (1 dim)
       - is_H-bond_acceptor  (N or O)                  (1 dim)
  3. FGW distance: fgw = ot.gromov.fused_gromov_wasserstein2(
         M=ot.dist(feat_q, feat_ref),  C1=C_q, C2=C_ref,
         p=uniform(n_q), q=uniform(n_ref), alpha=0.5)
     alpha=0.5 balances structural (graph) and feature (colour) contributions equally.
  4. Similarity to one reference  = exp(-fgw)  ∈ (0,1].
  5. Score(query) = max similarity over the 5 references.
  6. Each FGW call is wrapped in try/except → 0.0 on failure (degenerate molecule).

Caching
-------
Per-molecule (C, feat) tuples are stored once in ``data._cache["fgw_mol_features"]``
(a dict keyed by molecule index) so that different folds never recompute the same mol.

Speed note
----------
FGW is O(n^3) per pair via conditional-gradient iterations, where n = heavy-atom count
(median ~23 in this CCR2 dataset, max ~45 for references). Expect ~0.3–1 s/molecule on a
single CPU. Full 575-molecule evaluate() takes ~10–20 min; subsample to ~150 for a ~5 min
estimate. Tests use n=30 and finish in < 30 s.
"""
from __future__ import annotations

import numpy as np
import ot
import ot.gromov
from rdkit import Chem

from harness import BenchData, register

# --------------------------------------------------------------------------- atom features

_ATOM_IDX: dict = {"C": 0, "N": 1, "O": 2, "F": 3, "S": 4, "Cl": 5}  # 6 = "other"
_N_FEAT = 10  # 7 one-hot elements + aromatic + hbd + hba


def _atom_features(mol: Chem.Mol) -> np.ndarray:
    """Return (n_atoms, 10) float feature matrix for heavy-atom-only mol."""
    rows = []
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        feat = [0.0] * _N_FEAT
        feat[_ATOM_IDX.get(sym, 6)] = 1.0                   # element one-hot
        feat[7] = float(atom.GetIsAromatic())                # aromatic
        feat[8] = float(sym in ("N", "O") and atom.GetTotalNumHs() > 0)  # HBD
        feat[9] = float(sym in ("N", "O"))                   # HBA
        rows.append(feat)
    return np.array(rows, dtype=np.float64)


def _mol_repr(mol: Chem.Mol) -> tuple[np.ndarray, np.ndarray] | None:
    """Return (C_normalised, feat_matrix) for *mol* or None if degenerate."""
    try:
        mol = Chem.RemoveHs(mol)
        n = mol.GetNumAtoms()
        if n < 2:
            return None
        C = Chem.GetDistanceMatrix(mol).astype(np.float64)
        mx = C.max()
        if mx > 0:
            C = C / mx
        feat = _atom_features(mol)
        return C, feat
    except Exception:
        return None


# --------------------------------------------------------------------------- reference loading

def _load_refs(ref_sdf: str) -> list[tuple[np.ndarray, np.ndarray]]:
    """Load the 5 CCR2 reference ligands; return list of (C, feat) pairs."""
    results = []
    sup = Chem.SDMolSupplier(ref_sdf, removeHs=False)
    for mol in sup:
        if mol is None:
            continue
        r = _mol_repr(mol)
        if r is not None:
            results.append(r)
    return results


# --------------------------------------------------------------------------- FGW similarity

def _fgw_sim(C1: np.ndarray, feat1: np.ndarray,
             C2: np.ndarray, feat2: np.ndarray,
             alpha: float = 0.5) -> float:
    """FGW distance between two molecules → similarity in (0,1] or 0.0 on failure."""
    try:
        n1, n2 = len(feat1), len(feat2)
        p = np.ones(n1, dtype=np.float64) / n1
        q = np.ones(n2, dtype=np.float64) / n2
        M = ot.dist(feat1, feat2, metric="sqeuclidean")
        # normalise M to [0,1] to match C's [0,1] scale, so alpha genuinely
        # balances feature vs structure (else sqeuclidean spans [0,~5] and the
        # feature term silently dominates the structural Gromov term).
        mmax = M.max()
        if mmax > 0:
            M = M / mmax
        fgw = float(ot.gromov.fused_gromov_wasserstein2(
            M, C1, C2, p, q, loss_fun="square_loss", alpha=alpha,
        ))
        # fgw >= 0 normally → sim ∈ (0,1]; clip guards solver oscillation (fgw<0 → sim>1).
        return float(min(1.0, np.exp(-fgw)))
    except Exception:
        return 0.0


# --------------------------------------------------------------------------- scorer

@register("fused_gw")
def fused_gw(data: BenchData, train_idx: np.ndarray, test_idx: np.ndarray) -> np.ndarray:
    """Pose-free FGW similarity scorer (UNSUPERVISED — train_idx is ignored).

    Score(mol) = max_{ref ∈ 5 CCR2 refs} exp( -FGW(mol, ref) )

    FGW balances structural topology (graph distances, alpha=0.5) with
    atom-type features (element, aromaticity, HBD/HBA, alpha=0.5).
    """
    # -- load references (cached once per data object) --
    if "fgw_refs" not in data._cache:
        data._cache["fgw_refs"] = _load_refs(data.ref_sdf)
    refs = data._cache["fgw_refs"]

    if not refs:
        return np.zeros(len(test_idx), dtype=float)

    # -- build/retrieve per-molecule feature cache --
    if "fgw_mol_features" not in data._cache:
        data._cache["fgw_mol_features"] = {}
    mol_cache: dict = data._cache["fgw_mol_features"]

    mols = data.mols()
    scores = np.zeros(len(test_idx), dtype=np.float64)

    for out_i, mol_i in enumerate(test_idx):
        mol_i = int(mol_i)
        if mol_i not in mol_cache:
            mol = mols[mol_i]
            mol_cache[mol_i] = _mol_repr(mol) if mol is not None else None

        repr_q = mol_cache[mol_i]
        if repr_q is None:
            scores[out_i] = 0.0
            continue

        C_q, feat_q = repr_q
        best = 0.0
        for C_r, feat_r in refs:
            s = _fgw_sim(C_q, feat_q, C_r, feat_r, alpha=0.5)
            if s > best:
                best = s
        scores[out_i] = best

    return scores

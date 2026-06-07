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
    """Sum of pairwise Gaussian volume overlaps between two point sets."""
    d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)   # pairwise squared distance
    return float(np.exp(-alpha * d2).sum())


# --------------------------------------------------------------------- 3D embed + decompose
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem, rdShapeAlign  # noqa: E402

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

    Template uses its first conformer; the query is tried over all its conformers and the
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

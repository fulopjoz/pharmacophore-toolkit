"""RDKit shape+color combo scorer — production-ROCS proxy.

For each query molecule:
  1. Embed a single 3D conformer (ETKDGv3, seed=SEED); score 0.0 if embedding fails.
  2. Load the 5 CCR2 reference ligands from data.ref_sdf ONCE (cached on data._cache).
  3. Align probe to each reference via rdShapeAlign.AlignMol(ref, probe, opt_param=0.5)
     which optimises for shape+color jointly and returns (shape_tanimoto, color_tanimoto).
  4. Score = max over 5 refs of (shape + color)  ≈  TanimotoCombo ∈ [0, 2].

Embedded per-query conformers are cached in data._cache['_shape_conf_<idx>'] so that
repeated folds do not re-embed the same molecule.

UNSUPERVISED — train_idx is ignored.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeAlign

from harness import register, BenchData, SEED


def _load_refs(data: BenchData):
    """Load reference mols from data.ref_sdf, cached once on data._cache['_shape_refs']."""
    if "_shape_refs" not in data._cache:
        data._cache["_shape_refs"] = [
            m for m in Chem.SDMolSupplier(data.ref_sdf, removeHs=False) if m is not None
        ]
    return data._cache["_shape_refs"]


def _embed_mol(smi: str, seed: int) -> Chem.Mol | None:
    """Return a 3D mol (with Hs) or None if embedding fails."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) == -1:
        return None
    return mol


def _shape_combo(ref: Chem.Mol, probe: Chem.Mol) -> float:
    """TanimotoCombo = shape + color for one ref–probe pair.

    opt_param=0.5 balances shape and colour during the optimisation, so both
    components of the returned tuple are non-zero.
    """
    shape, color = rdShapeAlign.AlignMol(ref, probe, -1, -1, True, 0.5)
    return float(shape + color)


@register("shape_combo_rdkit")
def shape_combo_rdkit(data: BenchData, train_idx, test_idx) -> np.ndarray:
    """RDKit TanimotoCombo (shape+color) vs 5 CCR2 refs — unsupervised; train_idx ignored."""
    refs = _load_refs(data)
    scores = np.zeros(len(test_idx), dtype=float)

    for out_i, mol_i in enumerate(test_idx):
        cache_key = f"_shape_conf_{mol_i}"
        if cache_key not in data._cache:
            data._cache[cache_key] = _embed_mol(data.smiles[mol_i], SEED)
        probe = data._cache[cache_key]

        if probe is None:
            scores[out_i] = 0.0
            continue

        best = max(_shape_combo(ref, Chem.RWMol(probe)) for ref in refs)
        scores[out_i] = best

    return scores

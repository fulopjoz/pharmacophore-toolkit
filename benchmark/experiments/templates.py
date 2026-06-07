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
    Largest clusters are taken first so the templates cover the most actives.
    `seed` is accepted for API stability; Butina on a fixed distance matrix is
    already deterministic."""
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

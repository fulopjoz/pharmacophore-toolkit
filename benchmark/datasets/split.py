"""
benchmark/datasets/split.py
============================
Dataset splitting utilities for the pharmacophore-toolkit benchmark suite.

WHY SCAFFOLD SPLIT?
-------------------
Random splits allow the same chemical scaffold to appear in both training and
test sets, so a model can memorise scaffold-activity relationships rather than
learning transferable pharmacophore features.  Scaffold-disjoint ("group-
disjoint") splits ensure that every scaffold seen at test time is *new* to the
model, giving a cleaner estimate of generalisation to novel chemotypes.

This matters for unbiased optimizer selection: if two optimizers score
identically on a random-split benchmark but differently on a scaffold-split one,
the scaffold-split result reveals which optimizer truly generalises.

References
----------
Bemis, G. W. & Murcko, M. A. (1996).  The properties of known drugs. 1.
    Molecular frameworks.  J. Med. Chem., 39(15), 2887–2893.
Sheridan, R. P. (2013).  Time-split cross-validation as a method for estimating
    the goodness of prospective prediction.  J. Chem. Inf. Model., 53(4),
    783–790.
Hu, Y. & Bajorath, J. (2019).  Introducing a new compound generation benchmark
    that bridges scaffold- and property-based evaluation.  J. Chem. Inf. Model.,
    59(11), 4724–4734.
"""

from __future__ import annotations

import random
from collections import defaultdict
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _generic_scaffold(smiles: str) -> str:
    """Return the Bemis-Murcko *generic* scaffold SMILES for *smiles*.

    Generic scaffolds replace all atoms with carbon and all bonds with single
    bonds, making the grouping insensitive to heteroatom substitution.  Molecules
    that fail to parse or have no ring system are assigned the empty string as
    their scaffold (they form their own singleton group).

    Uses RDKit ``MurckoScaffold.GetScaffoldForMol`` followed by
    ``MurckoScaffold.MakeScaffoldGeneric``.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles  # keep invalid SMILES in their own singleton group
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    generic_mol  = MurckoScaffold.MakeScaffoldGeneric(scaffold_mol)
    return Chem.MolToSmiles(generic_mol)


def _scaffold_groups(smiles_list: List[str]) -> dict[str, List[int]]:
    """Map generic-scaffold SMILES -> list of indices into *smiles_list*."""
    groups: dict[str, List[int]] = defaultdict(list)
    for idx, smi in enumerate(smiles_list):
        key = _generic_scaffold(smi)
        groups[key].append(idx)
    return dict(groups)


def _split_groups_by_fraction(
    groups: dict[str, List[int]],
    test_frac: float,
    rng: random.Random,
) -> tuple[List[int], List[int]]:
    """Shuffle scaffold groups and assign whole groups to train / test.

    Groups are processed in random order; each group goes entirely into either
    the train or the test partition so that *no scaffold appears in both splits*.
    The split target is *test_frac* of the total number of molecules.

    Strategy: greedy — once the test set has reached *test_frac* of total
    molecules, remaining groups go to train.  This is deterministic given *rng*.
    """
    total = sum(len(v) for v in groups.values())
    target_test = total * test_frac

    group_keys = list(groups.keys())
    rng.shuffle(group_keys)

    train_idx: List[int] = []
    test_idx:  List[int] = []

    for key in group_keys:
        idxs = groups[key]
        # Assign whole group to test until target is met, then to train.
        if len(test_idx) < target_test:
            test_idx.extend(idxs)
        else:
            train_idx.extend(idxs)

    return train_idx, test_idx


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def scaffold_split(
    actives_smiles: List[str],
    decoys_smiles: List[str],
    test_frac: float = 0.25,
    seed: int = 42,
) -> dict:
    """Partition actives and decoys into train / test by Bemis-Murcko generic scaffold.

    **Guarantee**: no generic scaffold appears in both the train and test
    partitions of the same list (actives or decoys are split independently, so
    a scaffold shared between actives and decoys is handled separately for each).

    Parameters
    ----------
    actives_smiles:
        List of SMILES strings for active compounds.
    decoys_smiles:
        List of SMILES strings for decoy compounds.
    test_frac:
        Approximate fraction of molecules to place in the test set.  The actual
        fraction may differ because whole scaffold groups are moved together.
    seed:
        Random seed for reproducibility.

    Returns
    -------
    dict with keys:
        ``train_actives``, ``test_actives`` – SMILES lists (actives)
        ``train_decoys``,  ``test_decoys``  – SMILES lists (decoys)

    Notes
    -----
    WHY SCAFFOLD SPLIT MATTERS FOR OPTIMIZER SELECTION:
    A scaffold-disjoint split prevents data leakage where a model exploits
    memorised scaffold-activity patterns rather than truly generalising to new
    chemical space.  Optimizer benchmarks built on random splits systematically
    over-estimate performance on novel chemotypes (Sheridan, 2013).  Using
    generic scaffolds (Bemis & Murcko, 1996) further removes sensitivity to
    heteroatom decoration, making the grouping more conservative.
    """
    rng = random.Random(seed)

    act_groups = _scaffold_groups(actives_smiles)
    dec_groups = _scaffold_groups(decoys_smiles)

    act_train_idx, act_test_idx = _split_groups_by_fraction(act_groups, test_frac, rng)
    dec_train_idx, dec_test_idx = _split_groups_by_fraction(dec_groups, test_frac, rng)

    return {
        "train_actives": [actives_smiles[i] for i in act_train_idx],
        "test_actives":  [actives_smiles[i] for i in act_test_idx],
        "train_decoys":  [decoys_smiles[i]  for i in dec_train_idx],
        "test_decoys":   [decoys_smiles[i]  for i in dec_test_idx],
    }


def random_split(
    actives_smiles: List[str],
    decoys_smiles: List[str],
    test_frac: float = 0.25,
    seed: int = 42,
) -> dict:
    """Partition actives and decoys into train / test by simple random shuffle.

    Unlike :func:`scaffold_split`, this function makes **no scaffold-disjoint
    guarantee** — the same scaffold can appear in both train and test sets.  It
    is provided as a baseline for comparison to quantify the optimism introduced
    by ignoring scaffold structure.

    Parameters
    ----------
    actives_smiles:
        List of SMILES strings for active compounds.
    decoys_smiles:
        List of SMILES strings for decoy compounds.
    test_frac:
        Fraction of molecules to place in the test set.
    seed:
        Random seed for reproducibility.

    Returns
    -------
    dict with keys:
        ``train_actives``, ``test_actives``, ``train_decoys``, ``test_decoys``
    """
    rng = random.Random(seed)

    def _split(smiles_list: List[str]) -> tuple[List[str], List[str]]:
        shuffled = list(smiles_list)
        rng.shuffle(shuffled)
        n_test = max(1, round(len(shuffled) * test_frac)) if shuffled else 0
        return shuffled[n_test:], shuffled[:n_test]

    train_act, test_act = _split(actives_smiles)
    train_dec, test_dec = _split(decoys_smiles)

    return {
        "train_actives": train_act,
        "test_actives":  test_act,
        "train_decoys":  train_dec,
        "test_decoys":   test_dec,
    }

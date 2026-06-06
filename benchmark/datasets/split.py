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
    bonds, making the grouping insensitive to heteroatom substitution.

    Molecules with **no ring system** (acyclic) have an empty generic scaffold.
    Returning the empty string for all of them would collapse every chemically
    unrelated acyclic molecule into a single huge group that then moves as one
    indivisible block across the split.  Instead, each acyclic molecule is keyed
    by its own canonical SMILES so it forms its own singleton group (identical
    acyclic molecules still share a group).  Unparseable SMILES are likewise
    given their own singleton key.

    Uses RDKit ``MurckoScaffold.GetScaffoldForMol`` followed by
    ``MurckoScaffold.MakeScaffoldGeneric``.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f"__invalid__::{smiles}"          # own singleton group
    try:
        scaffold_mol   = MurckoScaffold.GetScaffoldForMol(mol)
        generic_mol    = MurckoScaffold.MakeScaffoldGeneric(scaffold_mol)
        generic_smiles = Chem.MolToSmiles(generic_mol)
    except Exception:
        return f"__error__::{Chem.MolToSmiles(mol)}"  # own singleton group
    if not generic_smiles:                        # acyclic -> own singleton group
        return f"__acyclic__::{Chem.MolToSmiles(mol)}"
    return generic_smiles


def _assign_test_scaffolds(
    group_sizes: dict[str, int],
    test_frac: float,
    rng: random.Random,
) -> set:
    """Choose whole scaffold groups for the TEST partition.

    Each scaffold group goes entirely to either train or test, so no scaffold
    crosses the split boundary.  The target test size is *test_frac* of the
    total molecule count.

    Strategy — **fill train largest-first, overflow to test** (the DeepChem
    ``ScaffoldSplitter`` polarity).  Groups are considered from largest to
    smallest (ties broken by a seeded shuffle for reproducibility) and added to
    *train* while doing so keeps train at or below its target size
    ``(1 - test_frac) * total``; every group that does not fit goes to test.

    Two properties motivate this design:

    * **No starvation / no silent overshoot.**  The naive ``len(test) < target``
      check assigned groups to test in random order, so a stream of small groups
      could fill test to just under the target and then a single large group
      arriving late would blow far past it -- in the worst case sweeping every
      remaining group into test and leaving *train* empty.  Filling train first,
      largest-first, removes that ordering dependence: whenever there are at
      least two scaffold groups, the group(s) that overflow the train target
      guarantee a non-empty test.

    * **Right generalisation semantics.**  Common (large) scaffolds land in
      train and the rarer (small) scaffolds become the test set, so the
      benchmark measures generalisation to *novel* chemotypes rather than
      re-scoring well-represented ones.

    Returns the set of scaffold keys assigned to the test partition.
    """
    total = sum(group_sizes.values())
    train_target = total * (1.0 - test_frac)

    keys = list(group_sizes.keys())
    rng.shuffle(keys)                                  # reproducible tie-break
    keys.sort(key=lambda k: group_sizes[k], reverse=True)

    train_keys: set = set()
    n_train = 0
    for key in keys:
        size = group_sizes[key]
        if n_train + size <= train_target:
            train_keys.add(key)
            n_train += size

    return {k for k in keys if k not in train_keys}


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
    partitions — *and this holds across the two lists jointly*.  Actives and
    decoys are pooled before grouping, so a scaffold shared by an active and a
    decoy is assigned to the same side for both (it never sits in, e.g.,
    train-actives and test-decoys at once).  This removes a subtle confound:
    with independent per-list splits the same scaffold could be learned as
    "active" in training and then appear as a test decoy, creating
    label-conflicting test cases that add noise to the enrichment estimate.

    If the dataset has too few distinct scaffolds in either class to place at
    least one molecule of each class in each partition, a ``ValueError`` is
    raised rather than silently returning a degenerate (empty-partition) split.

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
    smiles = list(actives_smiles) + list(decoys_smiles)
    labels = [1] * len(actives_smiles) + [0] * len(decoys_smiles)
    train_idx, test_idx = scaffold_split_indices(smiles, labels, test_frac, seed)

    return {
        "train_actives": [smiles[i] for i in train_idx if labels[i] == 1],
        "test_actives":  [smiles[i] for i in test_idx  if labels[i] == 1],
        "train_decoys":  [smiles[i] for i in train_idx if labels[i] == 0],
        "test_decoys":   [smiles[i] for i in test_idx  if labels[i] == 0],
    }


def scaffold_split_indices(
    smiles: List[str],
    labels: List[int],
    test_frac: float = 0.25,
    seed: int = 42,
    stratify: bool = True,
) -> tuple[List[int], List[int]]:
    """Index-level Bemis-Murcko scaffold split (the shared core of :func:`scaffold_split`).

    This is the form the benchmark harness consumes: it returns ``(train_idx,
    test_idx)`` into the *combined* ``smiles`` list so a scorer can fit on the
    train rows of a feature matrix and predict the test rows.

    Molecules are grouped by generic scaffold and whole groups are assigned
    train/test with the train-fill, largest-first strategy of
    :func:`_assign_test_scaffolds`.  No scaffold crosses the partition boundary.

    Parameters
    ----------
    smiles:
        SMILES for every molecule (any class ordering).
    labels:
        Parallel 0/1 labels (1 = active, 0 = decoy).  Used to guarantee each
        class appears in both partitions and (when ``stratify``) to balance the
        test fraction per class.
    test_frac, seed:
        As in :func:`scaffold_split`.
    stratify:
        If True (default), aim for ``test_frac`` of *each* class in the test set.
        Scaffolds unique to one class are split to hit that class's target;
        scaffolds shared by both classes are decided once (as a unit) so they
        still never cross the partition boundary.  This avoids the pooled split's
        failure mode where one class is starved in test (e.g. 3 of 74 actives) —
        which makes per-class enrichment metrics like BEDROC unreliable.
        If False, scaffolds are pooled across classes and split by total size.

    Returns
    -------
    ``(train_idx, test_idx)`` : two lists of integer indices into ``smiles``.

    Raises
    ------
    ValueError
        If either partition is empty, or either class is absent from either
        partition — a degenerate split whose enrichment metrics would be
        undefined.  Failing loudly beats returning a silently broken split.
    """
    if len(smiles) != len(labels):
        raise ValueError(f"smiles ({len(smiles)}) and labels ({len(labels)}) length mismatch")

    rng = random.Random(seed)
    scaffolds = [_generic_scaffold(s) for s in smiles]

    if stratify:
        act_count: dict[str, int] = defaultdict(int)
        dec_count: dict[str, int] = defaultdict(int)
        for key, lab in zip(scaffolds, labels):
            (act_count if lab == 1 else dec_count)[key] += 1
        all_keys = set(act_count) | set(dec_count)
        shared = {k for k in all_keys if act_count[k] > 0 and dec_count[k] > 0}
        act_only = {k: act_count[k] for k in all_keys if k not in shared and act_count[k] > 0}
        dec_only = {k: dec_count[k] for k in all_keys if k not in shared and dec_count[k] > 0}

        test_keys: set = set()
        # Shared scaffolds decided once (as whole-class-spanning units) so a
        # scaffold never lands train for actives + test for decoys.
        if shared:
            test_keys |= _assign_test_scaffolds({k: act_count[k] + dec_count[k] for k in shared},
                                                test_frac, rng)
        # Each class's own scaffolds split to hit that class's test_frac.
        if act_only:
            test_keys |= _assign_test_scaffolds(act_only, test_frac, rng)
        if dec_only:
            test_keys |= _assign_test_scaffolds(dec_only, test_frac, rng)
    else:
        group_sizes: dict[str, int] = defaultdict(int)
        for key in scaffolds:
            group_sizes[key] += 1
        test_keys = _assign_test_scaffolds(dict(group_sizes), test_frac, rng)

    train_idx = [i for i, key in enumerate(scaffolds) if key not in test_keys]
    test_idx  = [i for i, key in enumerate(scaffolds) if key in test_keys]

    def _has_both_classes(idx: List[int]) -> bool:
        seen = {labels[i] for i in idx}
        return 0 in seen and 1 in seen

    if not train_idx or not test_idx or not _has_both_classes(train_idx) or not _has_both_classes(test_idx):
        n_tr_act = sum(labels[i] == 1 for i in train_idx)
        n_te_act = sum(labels[i] == 1 for i in test_idx)
        n_tr_dec = sum(labels[i] == 0 for i in train_idx)
        n_te_dec = sum(labels[i] == 0 for i in test_idx)
        raise ValueError(
            "scaffold_split produced a degenerate partition "
            f"(train_actives={n_tr_act}, test_actives={n_te_act}, "
            f"train_decoys={n_tr_dec}, test_decoys={n_te_dec}). "
            "The dataset likely has too few distinct scaffolds in one class to "
            "split at this test_frac; use a larger/more diverse set, adjust "
            "test_frac, or fall back to random_split."
        )

    return train_idx, test_idx


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

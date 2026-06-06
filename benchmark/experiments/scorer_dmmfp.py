"""Differential Multimolecule Fingerprint (DMMFP) scorer for virtual screening.

Implements the method from Hutter (2022): for each Morgan fingerprint bit, compute
a differential weight as the mean frequency in actives minus the mean frequency in
decoys over the training set.  Test molecules are scored by the dot product of their
bit vectors with these weights — bits enriched in actives push the score up; bits
enriched in inactives push it down.  No external dependencies beyond numpy and RDKit
(both already required by the harness).

Reference: Hutter, M. C. (2022). Differential Multimolecule Fingerprints.
J. Chem. Inf. Model., 62(9), 2112–2120. https://doi.org/10.1021/acs.jcim.2c00242
"""
import numpy as np

from harness import register, BenchData


@register("differential_mmfp")
def differential_mmfp(data: BenchData, train_idx, test_idx) -> np.ndarray:
    """Score test molecules by the DMMFP dot-product rule (Hutter 2022).

    Parameters
    ----------
    data : BenchData
        Benchmark dataset providing shared Morgan fingerprint matrix via
        ``data.morgan()`` (shape N × 2048, dtype uint8).
    train_idx : array-like of int
        Indices of the training split used to derive differential weights.
    test_idx : array-like of int
        Indices of molecules to score; returned in the same order.

    Returns
    -------
    np.ndarray, shape (len(test_idx),)
        Scalar scores — higher means more active-like.
    """
    morgan = data.morgan().astype(np.float32)  # N x 2048

    train_idx = np.asarray(train_idx)
    test_idx = np.asarray(test_idx)

    y_train = data.y[train_idx]
    active_mask = y_train == 1
    decoy_mask = y_train == 0

    active_fps = morgan[train_idx[active_mask]]   # n_actives x 2048
    decoy_fps = morgan[train_idx[decoy_mask]]      # n_decoys  x 2048

    # Per-bit differential weight: freq(actives) − freq(inactives)
    w = active_fps.mean(axis=0) - decoy_fps.mean(axis=0)  # shape (2048,)

    # Score = dot product of test bit vector with weights
    return morgan[test_idx] @ w  # shape (len(test_idx),)

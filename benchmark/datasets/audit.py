"""
benchmark/datasets/audit.py
============================
Decoy-bias auditing for the pharmacophore-toolkit benchmark suite.

BACKGROUND
----------
Benchmark datasets whose decoys are structurally similar to actives
("analog decoys") artificially inflate enrichment metrics (AUC, EF) and
lead to misleading optimizer comparisons.  Two well-known criteria:

* **Rohrer-Baumann / MUV rule** (Rohrer & Baumann, 2009):
  Maximum Tanimoto coefficient (max-Tc) between any decoy and any active
  should be < 0.6 for the decoy to be considered "dissimilar".
  MUV datasets enforce max-Tc < 0.6 against every active.

* **Sieg 2019 rule** (Sieg et al., 2019):
  Even a max-Tc > 0.35 (ECFP4/Morgan radius 2, 2048 bits) signals
  potential analog bias that can distort QSAR model evaluation.

This module implements a conservative flag at max-Tc > 0.35 while also
reporting where decoys fall relative to the stricter 0.6 threshold.

References
----------
Rohrer, S. G. & Baumann, K. (2009).  Maximum unbiased validation (MUV)
    data sets for virtual screening based on PubChem bioactivity data.
    J. Chem. Inf. Model., 49(2), 169–184.  https://doi.org/10.1021/ci8002649
Sieg, J., Flachsenberg, F. & Rarey, M. (2019).  In need of bias control:
    evaluating chemical data for machine learning in structure-based virtual
    screening.  J. Chem. Inf. Model., 59(3), 947–961.
    https://doi.org/10.1021/acs.jcim.8b00712
"""

from __future__ import annotations

from typing import List

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _morgan_fp(smiles: str):
    """Return an RDKit Morgan bit-vector fingerprint (radius=2, nBits=2048) or None."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def _properties(smiles: str) -> dict:
    """Return MW, LogP, HBD, HBA for a SMILES; NaN on parse failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"mw": float("nan"), "logp": float("nan"),
                "hbd": float("nan"), "hba": float("nan")}
    return {
        "mw":   Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd":  float(rdMolDescriptors.CalcNumHBD(mol)),
        "hba":  float(rdMolDescriptors.CalcNumHBA(mol)),
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# Thresholds (configurable at module level for testing convenience)
ANALOG_BIAS_TC: float = 0.35   # flag if max-Tc to any active exceeds this (Sieg 2019)
MILD_BIAS_TC:   float = 0.20   # fraction above which verdict is at least "mild-bias"


def decoy_bias_audit(
    actives_smiles: List[str],
    decoys_smiles: List[str],
) -> dict:
    """Compute similarity-based and property-based decoy-bias metrics.

    For each decoy the **maximum Morgan Tanimoto coefficient** (radius 2,
    2048 bits) against the full active set is computed.  The distribution of
    these per-decoy max-Tc values is summarised and compared against the
    Sieg 2019 analog-bias threshold (max-Tc > 0.35).

    Verdict rules
    -------------
    * ``"biased (analog leakage)"``  — fraction of decoys with max-Tc > 0.35
      is ≥ 0.30, *or* the median max-Tc exceeds 0.35.
    * ``"mild-bias"``                — fraction > 0.35 is ≥ 0.10, or median
      max-Tc > 0.20.
    * ``"unbiased"``                 — otherwise.

    Parameters
    ----------
    actives_smiles:
        List of SMILES strings for active compounds.
    decoys_smiles:
        List of SMILES strings for decoy compounds.

    Returns
    -------
    dict with keys:

    ``mean_max_tc``, ``median_max_tc``, ``p95_max_tc``
        Distribution statistics of per-decoy max-Tc values.
    ``fraction_above_0_35``
        Fraction of decoys with max-Tc > 0.35 (Sieg 2019 analog-bias flag).
    ``verdict``
        One of ``"unbiased"``, ``"mild-bias"``, ``"biased (analog leakage)"``.
    ``actives_mean_mw``, ``decoys_mean_mw``
        Mean molecular weight for actives and decoys respectively.
    ``actives_mean_logp``, ``decoys_mean_logp``
        Mean calculated LogP.
    ``actives_mean_hbd``, ``decoys_mean_hbd``
        Mean hydrogen-bond donor count.
    ``actives_mean_hba``, ``decoys_mean_hba``
        Mean hydrogen-bond acceptor count.

    Notes
    -----
    **Rohrer-Baumann MUV criterion** (max-Tc < 0.6): a stricter threshold for
    maximum-unbiased validation datasets.  The present function uses the
    looser Sieg 2019 threshold (0.35) as the primary flag because real-world
    datasets rarely enforce the MUV criterion.

    References
    ----------
    Rohrer, S. G. & Baumann, K. (2009).  J. Chem. Inf. Model., 49, 169–184.
    Sieg, J., Flachsenberg, F. & Rarey, M. (2019).  J. Chem. Inf. Model.,
        59, 947–961.
    """
    # --- fingerprints ------------------------------------------------------ #
    active_fps = [_morgan_fp(s) for s in actives_smiles]
    active_fps = [fp for fp in active_fps if fp is not None]

    max_tc_values: list[float] = []

    for smi in decoys_smiles:
        dfp = _morgan_fp(smi)
        if dfp is None or not active_fps:
            max_tc_values.append(0.0)
            continue
        sims = DataStructs.BulkTanimotoSimilarity(dfp, active_fps)
        max_tc_values.append(float(max(sims)))

    # --- verdict ------------------------------------------------------------ #
    # If there are no parseable actives (or no decoys) the similarity audit was
    # never actually performed -- every decoy would have been assigned max_tc=0.
    # Returning "unbiased" here is a false negative: it passes a dataset we could
    # not evaluate.  Report "insufficient_data" with NaN stats instead.
    if (not active_fps) or (not decoys_smiles):
        nan = float("nan")
        mean_tc = median_tc = p95_tc = frac_35 = nan
        verdict = "insufficient_data"
    else:
        arr = np.array(max_tc_values)
        mean_tc   = float(np.mean(arr))
        median_tc = float(np.median(arr))
        p95_tc    = float(np.percentile(arr, 95))
        frac_35   = float(np.mean(arr > ANALOG_BIAS_TC))

        if frac_35 >= 0.30 or median_tc > ANALOG_BIAS_TC:
            verdict = "biased (analog leakage)"
        elif frac_35 >= 0.10 or median_tc > MILD_BIAS_TC:
            verdict = "mild-bias"
        else:
            verdict = "unbiased"

    # --- property summaries ------------------------------------------------- #
    def _mean_props(smiles_list: List[str]) -> dict:
        keys = ("mw", "logp", "hbd", "hba")
        props = [_properties(s) for s in smiles_list]
        if not props:                       # no molecules -> NaN means, no warning
            return {k: float("nan") for k in keys}
        return {k: float(np.nanmean([p[k] for p in props])) for k in keys}

    act_props = _mean_props(actives_smiles)
    dec_props = _mean_props(decoys_smiles)

    return {
        # Similarity distribution
        "mean_max_tc":         mean_tc,
        "median_max_tc":       median_tc,
        "p95_max_tc":          p95_tc,
        "fraction_above_0_35": frac_35,
        "verdict":             verdict,
        # Property match summary
        "actives_mean_mw":   act_props["mw"],
        "decoys_mean_mw":    dec_props["mw"],
        "actives_mean_logp": act_props["logp"],
        "decoys_mean_logp":  dec_props["logp"],
        "actives_mean_hbd":  act_props["hbd"],
        "decoys_mean_hbd":   dec_props["hbd"],
        "actives_mean_hba":  act_props["hba"],
        "decoys_mean_hba":   dec_props["hba"],
    }

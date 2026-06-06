"""Load the existing DrugEx CCR2 project actives/decoys into the uniform benchmark format.

Data files (bundled alongside this module):
  actives_ccr2_N75.csv  — col: SMILES  (74 valid actives after parsing)
  decoys_ccr2_N500.csv  — col: Smiles  (500 property-matched decoys)

Usage
-----
    from benchmark.datasets.ccr2_project.load import load
    actives, decoys, meta = load()

Returns
-------
actives_smiles : list[str]   — 74 canonical SMILES (CCR2 actives, pChEMBL >= 6)
decoys_smiles  : list[str]   — ~498-500 canonical SMILES (property-matched decoys)
meta           : dict        — contents of meta.json, augmented with runtime counts

Bias note
---------
Decoys are DUD-E-style property-matched: bias_suspect=True in meta.
Do NOT use for unbiased enrichment benchmarking without caveat.
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).parent.resolve()

_ACTIVES_CSV = HERE / "actives_ccr2_N75.csv"
_DECOYS_CSV = HERE / "decoys_ccr2_N500.csv"
_META_JSON = HERE / "meta.json"

# Expected approximate counts (used in validation)
_EXPECTED_N_ACT = 74
_EXPECTED_N_DEC_MIN = 490  # allow small parse failures
_EXPECTED_N_DEC_MAX = 502


def _csv_to_smiles(csv_path: Path, smiles_col: str) -> list[str]:
    """Parse a CSV and return canonical SMILES for valid molecules."""
    try:
        from rdkit import Chem
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
        _rdkit_available = True
    except ImportError:
        _rdkit_available = False

    import csv

    smiles_list = []
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            raw = row.get(smiles_col, "").strip()
            if not raw:
                continue
            if _rdkit_available:
                mol = Chem.MolFromSmiles(raw)
                if mol is None:
                    continue
                smi = Chem.MolToSmiles(mol)
            else:
                smi = raw  # fallback: return as-is (no validation)
            if smi:
                smiles_list.append(smi)
    return smiles_list


def load() -> tuple[list[str], list[str], dict]:
    """Return (actives_smiles, decoys_smiles, meta).

    Raises
    ------
    FileNotFoundError
        If the CSV files are not found next to this script.
    ValueError
        If the parsed counts deviate significantly from expected values.
    """
    for p in (_ACTIVES_CSV, _DECOYS_CSV, _META_JSON):
        if not p.exists():
            raise FileNotFoundError(f"Required file not found: {p}")

    actives_smiles = _csv_to_smiles(_ACTIVES_CSV, smiles_col="SMILES")
    decoys_smiles = _csv_to_smiles(_DECOYS_CSV, smiles_col="Smiles")

    # Validate counts
    if len(actives_smiles) != _EXPECTED_N_ACT:
        raise ValueError(
            f"Expected {_EXPECTED_N_ACT} valid actives, got {len(actives_smiles)}. "
            f"Check {_ACTIVES_CSV}."
        )
    if not (_EXPECTED_N_DEC_MIN <= len(decoys_smiles) <= _EXPECTED_N_DEC_MAX):
        raise ValueError(
            f"Expected {_EXPECTED_N_DEC_MIN}-{_EXPECTED_N_DEC_MAX} valid decoys, "
            f"got {len(decoys_smiles)}. Check {_DECOYS_CSV}."
        )

    with open(_META_JSON) as fh:
        meta = json.load(fh)

    meta["n_actives_loaded"] = len(actives_smiles)
    meta["n_decoys_loaded"] = len(decoys_smiles)

    return actives_smiles, decoys_smiles, meta


if __name__ == "__main__":
    actives, decoys, meta = load()
    print(f"Loaded {len(actives)} actives and {len(decoys)} decoys")
    print(f"Dataset: {meta['name']}")
    print(f"Bias suspect: {meta['bias_suspect']}")
    print(f"First active: {actives[0]}")
    print(f"First decoy:  {decoys[0]}")

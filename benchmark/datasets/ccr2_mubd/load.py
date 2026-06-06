"""Load the CCR2 MUBD benchmark (Xia et al. 2018, DOI:10.1021/acs.jcim.8b00004).

Usage
-----
    from benchmark.datasets.ccr2_mubd.load import load
    actives, decoys, meta = load()

If the SDF files are not present, ``fetch.py`` is invoked automatically.

Returns
-------
actives_smiles : list[str]   — 60 canonical SMILES for CCR2 actives (pIC50 >= 6)
decoys_smiles  : list[str]   — 2340 canonical SMILES for max-unbiased decoys
meta           : dict        — contents of meta.json
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).parent.resolve()

_ACTIVES_SDF = HERE / "ccr2_actives.sdf"
_DECOYS_SDF = HERE / "ccr2_decoys.sdf"
_META_JSON = HERE / "meta.json"


def _ensure_rdkit() -> None:
    if importlib.util.find_spec("rdkit") is None:
        raise ImportError("RDKit is required: pip install rdkit-pypi")


def _sdf_to_smiles(sdf_path: Path) -> list[str]:
    """Read an SDF file and return canonical SMILES for valid molecules."""
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=True)
    smiles = []
    for mol in suppl:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        if smi:
            smiles.append(smi)
    return smiles


def _auto_fetch() -> None:
    """Run fetch.py to download the SDF files."""
    fetch_script = HERE / "fetch.py"
    if not fetch_script.exists():
        raise FileNotFoundError(f"fetch.py not found at {fetch_script}")

    spec = importlib.util.spec_from_file_location("ccr2_mubd_fetch", fetch_script)
    fetch_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fetch_mod)
    fetch_mod.fetch(HERE)


def load(auto_fetch: bool = True) -> tuple[list[str], list[str], dict]:
    """Return (actives_smiles, decoys_smiles, meta).

    Parameters
    ----------
    auto_fetch : bool
        If True (default) and SDF files are missing, fetch.py is called
        automatically. Set to False to raise FileNotFoundError instead.
    """
    _ensure_rdkit()

    # ------------------------------------------------------------------
    # Ensure SDF files are present
    # ------------------------------------------------------------------
    if not _ACTIVES_SDF.exists() or not _DECOYS_SDF.exists():
        if not auto_fetch:
            raise FileNotFoundError(
                f"SDF files not found in {HERE}. Run fetch.py first, or call load(auto_fetch=True)."
            )
        print("[ccr2_mubd] SDF files not found — running fetch.py ...")
        _auto_fetch()

    # ------------------------------------------------------------------
    # Parse SDF → SMILES
    # ------------------------------------------------------------------
    actives_smiles = _sdf_to_smiles(_ACTIVES_SDF)
    decoys_smiles = _sdf_to_smiles(_DECOYS_SDF)

    if not actives_smiles:
        raise RuntimeError(f"No valid actives parsed from {_ACTIVES_SDF}")
    if not decoys_smiles:
        raise RuntimeError(f"No valid decoys parsed from {_DECOYS_SDF}")

    # ------------------------------------------------------------------
    # Load metadata
    # ------------------------------------------------------------------
    with open(_META_JSON) as fh:
        meta = json.load(fh)

    # Patch with runtime counts in case they differ from JSON
    meta["n_actives_loaded"] = len(actives_smiles)
    meta["n_decoys_loaded"] = len(decoys_smiles)

    return actives_smiles, decoys_smiles, meta


if __name__ == "__main__":
    actives, decoys, meta = load()
    print(f"Loaded {len(actives)} actives and {len(decoys)} decoys")
    print(f"Dataset: {meta['name']}")
    print(f"DOI: {meta['doi']}")
    print(f"First active: {actives[0]}")
    print(f"First decoy:  {decoys[0]}")

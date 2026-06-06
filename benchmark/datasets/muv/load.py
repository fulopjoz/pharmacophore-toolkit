"""
Load MUV (Maximum Unbiased Validation) benchmark data into the uniform benchmark format.

Returns (actives_smiles, decoys_smiles, meta) where:
  actives_smiles – list[str]  SMILES of confirmed actives (label == 1)
  decoys_smiles  – list[str]  SMILES of spatially-unbiased decoys (label == 0)
  meta           – dict matching benchmark/datasets/SPEC.md meta.json schema

Source: Rohrer SG & Baumann K (2009). J. Chem. Inf. Model. 49(2):169-184.
DOI:    10.1021/ci8002649
Data:   https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/muv.csv.gz

CSV format:
  Columns: MUV-466, MUV-548, ..., MUV-859, mol_id, smiles
  Values:  '1' = active, '0' = decoy/inactive, '' = unscreened

Decoy selection (max-unbiased):
  MUV decoys are selected from PubChem confirmatory assay negatives using an
  algorithm that minimises artificial enrichment by removing nearest-neighbour
  bias, MW/cLogP clustering artefacts, and topological redundancy.  Each assay
  has ~30 actives and ~15 000 decoys, giving a 1:500 ratio.
"""

import csv
import gzip
import json
from pathlib import Path

RAW_FILE = Path(__file__).parent / "raw" / "muv.csv.gz"
META_DIR = Path(__file__).parent

DOI        = "10.1021/ci8002649"
SOURCE_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/muv.csv.gz"
SOURCE     = "MUV"

MUV_TARGETS = [
    "MUV-466", "MUV-548", "MUV-600", "MUV-644", "MUV-652",
    "MUV-689", "MUV-692", "MUV-712", "MUV-713", "MUV-733",
    "MUV-737", "MUV-810", "MUV-832", "MUV-846", "MUV-852",
    "MUV-858", "MUV-859",
]

# Human-readable target descriptions (from Rohrer & Baumann 2009, Table 1)
_TARGET_DESCRIPTIONS = {
    "MUV-466": "SF1 (Steroidogenic factor 1, NR5A1)",
    "MUV-548": "PKA (cAMP-dependent protein kinase catalytic subunit, PRKACA)",
    "MUV-600": "THR (Thyroid hormone receptor, THRB)",
    "MUV-644": "Rab1 GTPase activating protein TBC1D20",
    "MUV-652": "HIV RT (Reverse transcriptase, p51/p66 heterodimer)",
    "MUV-689": "ER (Estrogen receptor alpha, ESR1)",
    "MUV-692": "Ras and Rab interactor 1 (RIN1)",
    "MUV-712": "HSP90α (Heat shock protein 90α)",
    "MUV-713": "Androgen receptor (AR, NR3C4)",
    "MUV-733": "GR (Glucocorticoid receptor, NR3C1)",
    "MUV-737": "Akt (Protein kinase B, PKB, AKT1)",
    "MUV-810": "FAK (Focal adhesion kinase, PTK2)",
    "MUV-832": "Caspase-3 (CASP3)",
    "MUV-846": "FXa (Coagulation factor Xa, F10)",
    "MUV-852": "Urokinase (uPA, PLAU)",
    "MUV-858": "D4 receptor (Dopamine receptor D4, DRD4)",
    "MUV-859": "SF1 inhibition (second assay, NR5A1)",
}


def _load_csv() -> list[dict]:
    """Read and return all rows from the gzipped CSV (skips blank rows)."""
    if not RAW_FILE.exists():
        raise FileNotFoundError(
            f"MUV raw file not found: {RAW_FILE}\n"
            "Run:  python benchmark/datasets/muv/fetch.py"
        )
    rows = []
    with gzip.open(RAW_FILE, "rt", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            smiles = (row.get("smiles") or "").strip()
            if smiles:  # skip blank rows (first row in this CSV is empty)
                rows.append(row)
    return rows


def load(target: str) -> tuple[list[str], list[str], dict]:
    """
    Load MUV actives + decoys for a single assay.

    Parameters
    ----------
    target : str
        MUV assay ID, e.g. ``"MUV-466"`` or ``"466"`` (prefix inferred).

    Returns
    -------
    (actives_smiles, decoys_smiles, meta)
        actives_smiles – SMILES strings for confirmed actives (label == 1)
        decoys_smiles  – SMILES strings for spatially-unbiased decoys (label == 0)
        meta           – dict matching SPEC.md schema; also written as meta_<target>.json
    """
    # Normalise target id
    if not target.startswith("MUV-"):
        target = f"MUV-{target}"
    if target not in MUV_TARGETS:
        raise ValueError(
            f"Unknown MUV target: {target!r}. "
            f"Valid targets: {MUV_TARGETS}"
        )

    rows = _load_csv()

    actives: list[str] = []
    decoys:  list[str] = []
    for row in rows:
        v   = (row.get(target) or "").strip()
        smi = (row.get("smiles") or "").strip()
        if not smi:
            continue
        if v == "1":
            actives.append(smi)
        elif v == "0":
            decoys.append(smi)
        # blank == unscreened → omit

    meta = {
        "target":      target,
        "target_name": _TARGET_DESCRIPTIONS.get(target, target),
        "source":      SOURCE,
        "doi":         DOI,
        "decoy_type":  "max-unbiased",
        "n_act":       len(actives),
        "n_dec":       len(decoys),
        "url":         SOURCE_URL,
        "notes": (
            "MUV decoys are spatially-unbiased PubChem confirmatory-assay negatives. "
            "Selection minimises nearest-neighbour, MW/cLogP clustering and topological "
            "artefacts (Rohrer & Baumann 2009, doi:10.1021/ci8002649). "
            "~30 actives : ~15 000 decoys per assay (1:500 ratio)."
        ),
    }

    # Write meta alongside this loader
    safe_id = target.replace("-", "_")
    meta_path = META_DIR / f"meta_{safe_id}.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    return actives, decoys, meta


def available_targets() -> list[str]:
    """Return the list of all 17 MUV target IDs (data present when raw file exists)."""
    if not RAW_FILE.exists():
        return []
    return list(MUV_TARGETS)


if __name__ == "__main__":
    import sys

    targets = available_targets()
    if not targets:
        print("MUV raw file not found.  Run: python fetch.py", file=sys.stderr)
        sys.exit(1)

    for t in targets:
        act, dec, meta = load(t)
        print(f"{t} ({meta['target_name']}): "
              f"{meta['n_act']} actives, {meta['n_dec']} decoys")
        if act:
            print(f"  sample active: {act[0]}")
        if dec:
            print(f"  sample decoy:  {dec[0]}")

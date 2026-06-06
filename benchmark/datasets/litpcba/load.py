"""
Load LIT-PCBA benchmark data into the uniform benchmark format.

Returns (actives_smiles, decoys_smiles, meta) where:
  actives_smiles  – list[str]  SMILES of confirmed active compounds
  decoys_smiles   – list[str]  SMILES of confirmed inactive compounds
                               (experimental inactives from PubChem dose-response assays)
  meta            – dict matching benchmark/datasets/SPEC.md meta.json schema

File format (raw/*.smi):
  Two-column whitespace-separated: <SMILES>  <PubChem_CID>
  No header line.
"""

import json
from pathlib import Path

RAW_DIR = Path(__file__).parent / "raw"
META_DIR = Path(__file__).parent

DOI = "10.1021/acs.jcim.0c00155"
SOURCE_URL = "http://drugdesign.unistra.fr/LIT-PCBA/"

# Reference counts from the original publication (Tran-Nguyen et al. 2020)
# Used to populate meta when raw files are absent.
_PUBLISHED_COUNTS = {
    "ADRB2":    {"n_act": 17,  "n_dec": 312483},
    "ALDH1":    {"n_act": 35,  "n_dec": 304073},
    "ESR1_ago": {"n_act": 354, "n_dec": 14975},
    "ESR1_ant": {"n_act": 1045,"n_dec": 13512},
    "FEN1":     {"n_act": 269, "n_dec": 65691},
    "GBA":      {"n_act": 127, "n_dec": 96534},
    "IDH1":     {"n_act": 104, "n_dec": 18285},
    "KAT2A":    {"n_act": 193, "n_dec": 48555},
    "MAPK1":    {"n_act": 2472,"n_dec": 33048},
    "MTORC1":   {"n_act": 686, "n_dec": 28000},
    "OPRK1":    {"n_act": 24,  "n_dec": 269816},
    "PKM2":     {"n_act": 645, "n_dec": 38936},
    "PPARG":    {"n_act": 294, "n_dec": 8042},
    "TP53":     {"n_act": 90,  "n_dec": 18700},
    "VDR":      {"n_act": 391, "n_dec": 48744},
}


def _read_smi(path: Path) -> list[str]:
    """Read a two-column *.smi file and return the first column (SMILES)."""
    smiles = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            smiles.append(line.split()[0])
    return smiles


def load(target: str) -> tuple[list[str], list[str], dict]:
    """
    Load LIT-PCBA actives + inactives for a single target.

    Parameters
    ----------
    target : str
        Target identifier, e.g. "ADRB2" or "OPRK1".
        Raw files must exist under raw/<target>/ (run fetch.py first).

    Returns
    -------
    (actives_smiles, inactives_smiles, meta)
    """
    act_path = RAW_DIR / target / "actives.smi"
    inact_path = RAW_DIR / target / "inactives.smi"

    if not act_path.exists() or not inact_path.exists():
        raise FileNotFoundError(
            f"Raw files for target '{target}' not found under {RAW_DIR}/{target}/. "
            "Run: python fetch.py --targets " + target
        )

    actives = _read_smi(act_path)
    inactives = _read_smi(inact_path)

    meta = {
        "target": target,
        "source": "LIT-PCBA",
        "doi": DOI,
        "decoy_type": "experimental-inactive",
        "n_act": len(actives),
        "n_dec": len(inactives),
        "url": SOURCE_URL,
        "notes": (
            "Inactives are experimentally confirmed from PubChem dose-response assays; "
            "actives and inactives are within similar molecular property ranges (AVE-unbiased). "
            "SMILES column 1, PubChem CID column 2."
        ),
    }

    # Persist meta.json alongside this loader
    meta_path = META_DIR / f"meta_{target}.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    return actives, inactives, meta


def available_targets() -> list[str]:
    """Return target IDs for which raw files are present."""
    return sorted(
        d.name for d in RAW_DIR.iterdir()
        if d.is_dir()
        and (d / "actives.smi").exists()
        and (d / "inactives.smi").exists()
    )


if __name__ == "__main__":
    import sys
    targets = available_targets()
    if not targets:
        print("No targets available. Run: python fetch.py", file=sys.stderr)
        sys.exit(1)
    for t in targets:
        act, inact, meta = load(t)
        print(f"{t}: {meta['n_act']} actives, {meta['n_dec']} inactives")
        print(f"  sample active:   {act[0]}")
        print(f"  sample inactive: {inact[0]}")
        print(f"  meta written:    meta_{t}.json")

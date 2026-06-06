"""
Load DUD-E benchmark data into the uniform benchmark format.

Returns (actives_smiles, decoys_smiles, meta) where:
  actives_smiles – list[str]  SMILES of confirmed active compounds (from ChEMBL)
  decoys_smiles  – list[str]  SMILES of property-matched synthetic decoys (from ZINC)
  meta           – dict matching benchmark/datasets/SPEC.md meta.json schema

File format (raw/<target>/actives_final.ism):
  Three-column: <SMILES>  <ZINC_ID>  <CHEMBL_ID>   (actives)
  Two-column:   <SMILES>  <ZINC_ID>                 (decoys)
  No header line.

⚠️  BIAS WARNING: DUD-E decoys are property-matched synthetic compounds, NOT
experimental inactives. The dataset is known to have analogue bias and decoy
bias. Use only for bias-control diagnostics alongside LIT-PCBA.
"""

import json
from pathlib import Path

RAW_DIR = Path(__file__).parent / "raw"
META_DIR = Path(__file__).parent

DOI = "10.1021/jm300687e"
SOURCE_URL = "https://dude.docking.org/"

# Human-readable target names for meta
_TARGET_NAMES = {
    "aa2ar":  "Adenosine A2a receptor",
    "adrb1":  "Beta-1 adrenergic receptor",
    "adrb2":  "Beta-2 adrenergic receptor",
    "cxcr4":  "C-X-C chemokine receptor type 4",
    "drd3":   "Dopamine D3 receptor",
}


def _read_ism(path: Path) -> list[str]:
    """Read a DUD-E ISM file and return first column (SMILES)."""
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
    Load DUD-E actives + decoys for a single target.

    Parameters
    ----------
    target : str
        Target identifier (lowercase), e.g. "adrb2" or "aa2ar".
        Raw files must exist under raw/<target>/ (run fetch.py first).

    Returns
    -------
    (actives_smiles, decoys_smiles, meta)
    """
    act_path = RAW_DIR / target / "actives_final.ism"
    dec_path = RAW_DIR / target / "decoys_final.ism"

    if not act_path.exists() or not dec_path.exists():
        raise FileNotFoundError(
            f"Raw files for target '{target}' not found under {RAW_DIR}/{target}/. "
            "Run: python fetch.py --targets " + target
        )

    actives = _read_ism(act_path)
    decoys  = _read_ism(dec_path)

    meta = {
        "target": target,
        "target_name": _TARGET_NAMES.get(target, target.upper()),
        "source": "DUD-E",
        "doi": DOI,
        "decoy_type": "property-matched",
        "bias_warning": (
            "DUD-E decoys are synthetic property-matched compounds from ZINC, "
            "NOT experimental inactives. Known to exhibit analogue bias and decoy bias. "
            "See: doi:10.1371/journal.pone.0220113"
        ),
        "n_act": len(actives),
        "n_dec": len(decoys),
        "url": SOURCE_URL,
        "notes": (
            "Actives from ChEMBL (clustered by Bemis-Murcko scaffold). "
            "50 property-matched decoys per active, topology-dissimilar, "
            "net-charge matched. SMILES column 1, ZINC_ID column 2."
        ),
    }

    # Persist meta.json alongside this loader
    meta_path = META_DIR / f"meta_{target}.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    return actives, decoys, meta


def available_targets() -> list[str]:
    """Return target IDs for which raw files are present."""
    return sorted(
        d.name for d in RAW_DIR.iterdir()
        if d.is_dir()
        and (d / "actives_final.ism").exists()
        and (d / "decoys_final.ism").exists()
    )


if __name__ == "__main__":
    import sys
    targets = available_targets()
    if not targets:
        print("No targets available. Run: python fetch.py", file=sys.stderr)
        sys.exit(1)
    for t in targets:
        act, dec, meta = load(t)
        print(f"{t} ({meta['target_name']}): {meta['n_act']} actives, {meta['n_dec']} decoys")
        print(f"  sample active: {act[0]}")
        print(f"  sample decoy:  {dec[0]}")
        print(f"  meta written:  meta_{t}.json")

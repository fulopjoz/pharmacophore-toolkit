"""
load.py — Uniform loader for ChEMBL-derived actives + MUV-style unbiased decoys.

Conforms to benchmark/datasets/SPEC.md interface:

    load(target) -> (actives_smiles, decoys_smiles, meta)

Available targets: CCR2, CCR5, CXCR4.

Run ``python build.py`` first to generate the CSV files.
"""

from __future__ import annotations

import json
from pathlib import Path

_HERE = Path(__file__).parent


def load(target: str) -> tuple[list[str], list[str], dict]:
    """
    Load ChEMBL-derived actives and MUV-style unbiased decoys for a GPCR target.

    Parameters
    ----------
    target : str
        One of "CCR2", "CCR5", "CXCR4".

    Returns
    -------
    actives_smiles : list[str]
        Canonical SMILES of active compounds (pChEMBL >= 7, i.e. IC50/Ki <= 100 nM).
    decoys_smiles : list[str]
        Canonical SMILES of property-matched, similarity-unbiased decoy compounds.
    meta : dict
        Metadata dict conforming to SPEC.md meta.json schema.

    Raises
    ------
    FileNotFoundError
        If ``actives.csv`` or ``decoys.csv`` are absent for the requested target.
        Run ``python build.py --targets <target>`` to generate them.
    """
    target = target.upper()
    tgt_dir = _HERE / target

    act_path = tgt_dir / "actives.csv"
    dec_path = tgt_dir / "decoys.csv"
    meta_path = tgt_dir / "meta.json"

    if not act_path.exists() or not dec_path.exists():
        raise FileNotFoundError(
            f"Dataset files for target '{target}' not found under {tgt_dir}. "
            "Run: python build.py --targets " + target
        )

    actives = _read_csv(act_path)
    decoys = _read_csv(dec_path)

    if meta_path.exists():
        with open(meta_path) as f:
            meta = json.load(f)
    else:
        meta = {
            "target": target,
            "source": "ChEMBL REST API",
            "decoy_type": "created-muv-unbiased",
            "n_act": len(actives),
            "n_dec": len(decoys),
        }

    return actives, decoys, meta


def _read_csv(path: Path) -> list[str]:
    """Read a single-column CSV file with header 'smiles'; return SMILES list."""
    lines = path.read_text().splitlines()
    if not lines:
        return []
    # Skip header
    if lines[0].strip().lower() == "smiles":
        lines = lines[1:]
    return [line.strip() for line in lines if line.strip()]


def available_targets() -> list[str]:
    """Return target IDs for which CSV files are present."""
    targets = []
    for d in _HERE.iterdir():
        if d.is_dir() and (d / "actives.csv").exists() and (d / "decoys.csv").exists():
            targets.append(d.name)
    return sorted(targets)


if __name__ == "__main__":
    import sys

    targets = available_targets()
    if not targets:
        print("No datasets found. Run: python build.py", file=sys.stderr)
        sys.exit(1)

    for t in targets:
        act, dec, meta = load(t)
        print(f"{t}: {len(act)} actives, {len(dec)} decoys")
        if act:
            print(f"  sample active: {act[0][:60]}...")
        if dec:
            print(f"  sample decoy:  {dec[0][:60]}...")

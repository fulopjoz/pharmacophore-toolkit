"""
TDD-lite: verify load.py returns non-empty lists with valid SMILES.
Run with: python test_load.py
"""

import sys
from pathlib import Path

# Allow running from repo root or from this directory
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def validate_smiles_rdkit(smiles_list: list[str], label: str, n_check: int = 50) -> int:
    """Return count of invalid SMILES (0 = all good). Prints failures."""
    try:
        from rdkit import Chem
    except ImportError:
        print("  WARNING: RDKit not available, skipping SMILES validation")
        return 0

    invalid = 0
    for smi in smiles_list[:n_check]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  INVALID SMILES in {label}: {smi!r}")
            invalid += 1
    return invalid


def test_litpcba_load():
    from load import load, available_targets

    targets = available_targets()
    assert targets, "No LIT-PCBA targets found — run fetch.py first"
    print(f"LIT-PCBA targets available: {targets}")

    for target in targets:
        print(f"\n  Testing target: {target}")
        actives, inactives, meta = load(target)

        assert isinstance(actives, list),   f"actives must be list, got {type(actives)}"
        assert isinstance(inactives, list), f"inactives must be list, got {type(inactives)}"
        assert len(actives) > 0,   f"{target}: actives list is empty"
        assert len(inactives) > 0, f"{target}: inactives list is empty"

        assert meta["source"] == "LIT-PCBA"
        assert meta["decoy_type"] == "experimental-inactive"
        assert meta["n_act"] == len(actives)
        assert meta["n_dec"] == len(inactives)
        assert "doi" in meta

        bad_act = validate_smiles_rdkit(actives, f"{target}/actives")
        bad_inact = validate_smiles_rdkit(inactives, f"{target}/inactives")
        assert bad_act == 0,   f"{target}: {bad_act} invalid SMILES in actives"
        assert bad_inact == 0, f"{target}: {bad_inact} invalid SMILES in inactives (first 50)"

        print(f"  PASS  {target}: {len(actives)} actives, {len(inactives)} inactives")

    print("\nAll LIT-PCBA load tests PASSED")


if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent)
    test_litpcba_load()

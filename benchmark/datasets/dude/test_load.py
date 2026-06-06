"""
TDD-lite: verify load.py returns non-empty lists with valid SMILES.
Run with: python test_load.py
"""

import sys
from pathlib import Path

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


def test_dude_load():
    from load import load, available_targets

    targets = available_targets()
    assert targets, "No DUD-E targets found — run fetch.py first"
    print(f"DUD-E targets available: {targets}")

    for target in targets:
        print(f"\n  Testing target: {target}")
        actives, decoys, meta = load(target)

        assert isinstance(actives, list), f"actives must be list, got {type(actives)}"
        assert isinstance(decoys, list),  f"decoys must be list, got {type(decoys)}"
        assert len(actives) > 0, f"{target}: actives list is empty"
        assert len(decoys) > 0,  f"{target}: decoys list is empty"

        assert meta["source"] == "DUD-E"
        assert meta["decoy_type"] == "property-matched"
        assert meta["n_act"] == len(actives)
        assert meta["n_dec"] == len(decoys)
        assert "doi" in meta
        assert "bias_warning" in meta, "DUD-E meta must include bias_warning"

        bad_act = validate_smiles_rdkit(actives, f"{target}/actives")
        bad_dec = validate_smiles_rdkit(decoys, f"{target}/decoys")
        assert bad_act == 0, f"{target}: {bad_act} invalid SMILES in actives"
        assert bad_dec == 0, f"{target}: {bad_dec} invalid SMILES in decoys (first 50)"

        print(f"  PASS  {target}: {len(actives)} actives, {len(decoys)} decoys")

    print("\nAll DUD-E load tests PASSED")


if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent)
    test_dude_load()

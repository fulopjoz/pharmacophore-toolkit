"""
TDD-lite: verify load.py returns ~30 active SMILES (all RDKit-valid) + many decoys.

Run with:
    python test_load.py
or:
    cd benchmark/datasets/muv && python test_load.py
"""

import sys
from pathlib import Path

# Allow running from repo root
sys.path.insert(0, str(Path(__file__).parent))


def validate_smiles_rdkit(smiles_list: list[str], label: str, n_check: int = 50) -> int:
    """Return count of invalid SMILES (0 = all good). Prints each failure."""
    try:
        from rdkit import Chem
    except ImportError:
        print("  WARNING: RDKit not available — skipping SMILES validation")
        return 0

    invalid = 0
    for smi in smiles_list[:n_check]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  INVALID SMILES in {label}: {smi!r}")
            invalid += 1
    return invalid


def test_muv_load_single_target():
    """
    Load one MUV target and assert:
      - ~30 actives (between 20 and 40 inclusive)
      - at least 10 000 decoys
      - all active SMILES are RDKit-parseable
      - first 50 decoy SMILES are RDKit-parseable
      - meta fields match SPEC.md schema
    """
    from load import load, available_targets

    targets = available_targets()
    assert targets, (
        "No MUV targets available — run fetch.py first:\n"
        "  python benchmark/datasets/muv/fetch.py"
    )

    # Test with first available target (MUV-466 if file is present)
    target = targets[0]
    print(f"Testing target: {target}")
    actives, decoys, meta = load(target)

    # Type checks
    assert isinstance(actives, list), f"actives must be list, got {type(actives)}"
    assert isinstance(decoys,  list), f"decoys must be list, got {type(decoys)}"

    # Count checks
    assert 20 <= len(actives) <= 40, (
        f"{target}: expected ~30 actives (20-40), got {len(actives)}"
    )
    assert len(decoys) >= 10_000, (
        f"{target}: expected >=10 000 decoys, got {len(decoys)}"
    )

    # Meta schema checks (SPEC.md)
    for field in ("target", "source", "doi", "decoy_type", "n_act", "n_dec", "url"):
        assert field in meta, f"meta missing required field '{field}'"
    assert meta["source"]     == "MUV",           f"meta source wrong: {meta['source']!r}"
    assert meta["decoy_type"] == "max-unbiased",  f"meta decoy_type wrong: {meta['decoy_type']!r}"
    assert meta["doi"]        == "10.1021/ci8002649"
    assert meta["n_act"]      == len(actives),    "meta n_act mismatch"
    assert meta["n_dec"]      == len(decoys),     "meta n_dec mismatch"

    # SMILES validity
    bad_act = validate_smiles_rdkit(actives,       f"{target}/actives", n_check=len(actives))
    bad_dec = validate_smiles_rdkit(decoys[:50],   f"{target}/decoys",  n_check=50)
    assert bad_act == 0, f"{target}: {bad_act} invalid SMILES in actives"
    assert bad_dec == 0, f"{target}: {bad_dec} invalid SMILES in decoys (first 50)"

    print(f"  PASS  {target}: {len(actives)} actives, {len(decoys)} decoys — all SMILES valid")


def test_muv_available_targets():
    """All 17 MUV targets must be reported as available once the CSV is downloaded."""
    from load import available_targets, MUV_TARGETS

    targets = available_targets()
    if not targets:
        print("  SKIP  test_muv_available_targets — raw file not downloaded yet")
        return

    assert len(targets) == 17, f"Expected 17 MUV targets, got {len(targets)}: {targets}"
    assert set(targets) == set(MUV_TARGETS), "Target list mismatch"
    print(f"  PASS  available_targets() returns all 17 MUV targets")


if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent)

    print("=" * 60)
    print("MUV benchmark load tests")
    print("=" * 60)

    test_muv_available_targets()
    test_muv_load_single_target()

    print("\nAll MUV load tests PASSED")

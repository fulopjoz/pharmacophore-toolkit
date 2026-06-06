"""TDD-lite tests for benchmark/datasets/ccr2_project/load.py.

Run with:
    python -m pytest benchmark/datasets/ccr2_project/test_load.py -v
or:
    python benchmark/datasets/ccr2_project/test_load.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running as a standalone script without pytest on sys.path
HERE = Path(__file__).parent.resolve()
REPO_ROOT = HERE.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark.datasets.ccr2_project.load import load  # noqa: E402


def _validate_smiles_list(smiles_list: list[str], label: str) -> None:
    """Ensure every entry is a non-empty string that RDKit can parse."""
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")

    assert isinstance(smiles_list, list), f"{label}: expected list, got {type(smiles_list)}"
    for i, smi in enumerate(smiles_list):
        assert isinstance(smi, str) and smi, f"{label}[{i}]: empty or non-string SMILES"
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"{label}[{i}]: invalid SMILES: {smi!r}"


def test_load_returns_three_items() -> None:
    result = load()
    assert len(result) == 3, "load() must return a 3-tuple (actives, decoys, meta)"


def test_actives_count() -> None:
    actives, _, _ = load()
    assert len(actives) == 74, f"Expected 74 valid actives, got {len(actives)}"


def test_decoys_count() -> None:
    _, decoys, _ = load()
    assert 490 <= len(decoys) <= 502, (
        f"Expected 490-502 valid decoys (nominal 500), got {len(decoys)}"
    )


def test_actives_are_valid_smiles() -> None:
    actives, _, _ = load()
    _validate_smiles_list(actives, "actives")


def test_decoys_are_valid_smiles() -> None:
    _, decoys, _ = load()
    _validate_smiles_list(decoys, "decoys")


def test_meta_fields() -> None:
    _, _, meta = load()
    assert isinstance(meta, dict), "meta must be a dict"
    assert meta.get("decoy_type") == "property-matched"
    assert meta.get("bias_suspect") is True
    assert meta.get("n_actives") == 74
    assert meta.get("n_decoys") == 500


def test_meta_runtime_counts_present() -> None:
    _, _, meta = load()
    assert "n_actives_loaded" in meta
    assert "n_decoys_loaded" in meta
    assert meta["n_actives_loaded"] == 74
    assert 490 <= meta["n_decoys_loaded"] <= 502


def test_no_duplicate_smiles_in_actives() -> None:
    actives, _, _ = load()
    assert len(set(actives)) == len(actives), (
        f"Duplicate SMILES found in actives ({len(actives) - len(set(actives))} duplicates)"
    )


# ---------------------------------------------------------------------------
# Standalone runner (no pytest required)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    tests = [
        test_load_returns_three_items,
        test_actives_count,
        test_decoys_count,
        test_actives_are_valid_smiles,
        test_decoys_are_valid_smiles,
        test_meta_fields,
        test_meta_runtime_counts_present,
        test_no_duplicate_smiles_in_actives,
    ]

    passed = 0
    failed = 0
    for fn in tests:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
            passed += 1
        except Exception as exc:  # noqa: BLE001
            print(f"  FAIL  {fn.__name__}: {exc}")
            failed += 1

    print(f"\n{passed}/{passed + failed} tests passed.")
    sys.exit(0 if failed == 0 else 1)

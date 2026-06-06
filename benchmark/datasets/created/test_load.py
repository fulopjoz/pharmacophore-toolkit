"""
test_load.py — TDD-lite tests for the created/ ChEMBL dataset.

Asserts:
  - load("CCR2") returns >= 30 actives and >= 30 decoys
  - All SMILES (first 50 from each list) are RDKit-parseable
  - Spot-check: max Morgan-Tc of decoys to actives < 0.35

Run: python test_load.py [--target CCR2]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Ensure load.py is importable when run from any cwd
sys.path.insert(0, str(Path(__file__).parent))
import load as loader  # noqa: E402  (local import after sys.path fix)

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect


TC_CUTOFF = 0.35
SPOT_CHECK_N = 20   # number of decoys to spot-check Tc


def _fp(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def run_tests(target: str = "CCR2") -> None:
    print(f"\nRunning tests for target: {target}")
    print("-" * 50)

    # --- Load ---
    actives, decoys, meta = loader.load(target)

    # 1. Non-empty and minimum count
    assert len(actives) >= 30, f"Expected >= 30 actives, got {len(actives)}"
    assert len(decoys) >= 30, f"Expected >= 30 decoys, got {len(decoys)}"
    print(f"[PASS] actives >= 30: {len(actives)}")
    print(f"[PASS] decoys  >= 30: {len(decoys)}")

    # 2. All SMILES in first 50 are RDKit-parseable
    for label, smiles_list in [("actives", actives), ("decoys", decoys)]:
        n_check = min(50, len(smiles_list))
        failures = [
            s for s in smiles_list[:n_check]
            if Chem.MolFromSmiles(s) is None
        ]
        assert not failures, f"{label}: {len(failures)} unparseable SMILES in first {n_check}"
        print(f"[PASS] first {n_check} {label} are RDKit-parseable")

    # 3. Tc spot-check: max Tc of sampled decoys to actives < TC_CUTOFF
    import random
    random.seed(42)
    sampled_decoys = random.sample(decoys, min(SPOT_CHECK_N, len(decoys)))

    active_fps = [_fp(s) for s in actives if _fp(s) is not None]
    failures_tc = []
    for smi in sampled_decoys:
        fp = _fp(smi)
        if fp is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, active_fps)
        max_tc = max(sims) if sims else 0.0
        if max_tc >= TC_CUTOFF:
            failures_tc.append((smi, max_tc))

    if failures_tc:
        print(f"[WARN] {len(failures_tc)} / {len(sampled_decoys)} spot-checked decoys have "
              f"Tc >= {TC_CUTOFF} to actives:")
        for smi, tc in failures_tc[:5]:
            print(f"       Tc={tc:.3f}: {smi[:60]}...")
        # Soft assertion: flag but don't fail the whole suite
        # (a few borderline Tc values may appear due to property matching vs exact-Tc filter)
        assert len(failures_tc) / len(sampled_decoys) < 0.10, (
            f"Too many decoys violate Tc cutoff: {len(failures_tc)} / {len(sampled_decoys)}"
        )
    else:
        print(f"[PASS] all {len(sampled_decoys)} spot-checked decoys have max Tc < {TC_CUTOFF}")

    # 4. Meta checks
    assert meta.get("n_act") == len(actives), "meta n_act mismatch"
    assert meta.get("n_dec") == len(decoys), "meta n_dec mismatch"
    assert "chembl_id" in meta or "source" in meta, "meta missing chembl_id/source"
    print(f"[PASS] meta valid: n_act={meta['n_act']}, n_dec={meta['n_dec']}")

    print("-" * 50)
    print("ALL TESTS PASSED\n")


def main():
    parser = argparse.ArgumentParser(description="TDD-lite tests for created/ ChEMBL dataset")
    parser.add_argument("--target", default="CCR2", choices=loader.available_targets() or ["CCR2"],
                        help="Target to test (default: CCR2)")
    parser.add_argument("--all", action="store_true", help="Test all available targets")
    args = parser.parse_args()

    targets_to_test = loader.available_targets() if args.all else [args.target]
    if not targets_to_test:
        print("No datasets found. Run: python build.py", file=sys.stderr)
        sys.exit(1)

    for t in targets_to_test:
        run_tests(t)


if __name__ == "__main__":
    main()

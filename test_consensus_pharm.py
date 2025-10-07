#!/usr/bin/env python
"""
Test script for consensus pharmacophore functionality.

This script tests the consensus_pharm method with various scenarios:
1. Multiple conformations of the same molecule
2. Similar molecules with alignment
3. Different distance thresholds
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from pharmacophore import Pharmacophore
import numpy as np


def test_same_molecule_conformations():
    """Test consensus with multiple conformations of the same molecule."""
    print("=" * 60)
    print("Test 1: Same molecule, different conformations")
    print("=" * 60)
    
    smiles = "CCO"  # Ethanol
    mols = []
    
    # Generate 3 different conformations
    for i in range(3):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42 + i)
        mols.append(mol)
    
    pharm = Pharmacophore()
    consensus = pharm.consensus_pharm(mols, distance_threshold=2.0)
    
    print(f"\nNumber of molecules: {len(mols)}")
    print(f"Number of consensus features: {len(consensus)}")
    print("\nConsensus features:")
    for f in consensus:
        print(f"  {f[0]}: ({f[2]:.3f}, {f[3]:.3f}, {f[4]:.3f})")
    
    assert len(consensus) > 0, "Should have at least one consensus feature"
    print("\n✓ Test passed!")
    return True


def test_aligned_molecules():
    """Test consensus with aligned similar molecules."""
    print("\n" + "=" * 60)
    print("Test 2: Aligned similar molecules")
    print("=" * 60)
    
    # Using molecules from the tutorial
    smiles_list = [
        "NCCc1c[nH]c2ccc(O)cc12",  # Serotonin
        "CC(C)NCC(O)c1ccc(O)c(O)c1",  # Dobutamine
    ]
    
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    mols = [Chem.AddHs(m) for m in mols]
    
    # Generate 3D conformations
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    for m in mols:
        AllChem.EmbedMolecule(m, ps)
    
    # Align molecules
    mol_aligned = rdMolAlign.GetO3A(mols[1], mols[0])
    rmsd = mol_aligned.Align()
    
    print(f"\nAlignment RMSD: {rmsd:.3f} Å")
    
    # Calculate consensus
    pharm = Pharmacophore()
    consensus = pharm.consensus_pharm(mols, distance_threshold=2.0)
    
    print(f"Number of molecules: {len(mols)}")
    print(f"Number of consensus features: {len(consensus)}")
    
    # Group by type
    by_type = {}
    for f in consensus:
        feat_type = f[0]
        if feat_type not in by_type:
            by_type[feat_type] = []
        by_type[feat_type].append(f)
    
    print("\nConsensus features by type:")
    for feat_type, features in sorted(by_type.items()):
        print(f"  {feat_type}: {len(features)} features")
        for f in features:
            print(f"    ({f[2]:.3f}, {f[3]:.3f}, {f[4]:.3f})")
    
    assert len(consensus) > 0, "Should have consensus features"
    print("\n✓ Test passed!")
    return True


def test_distance_threshold():
    """Test different distance thresholds."""
    print("\n" + "=" * 60)
    print("Test 3: Different distance thresholds")
    print("=" * 60)
    
    # Create slightly shifted versions of the same molecule
    smiles = "c1ccc(O)cc1"  # Phenol
    mols = []
    
    for i in range(3):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42 + i)
        mols.append(mol)
    
    pharm = Pharmacophore()
    
    # Test different thresholds
    thresholds = [1.0, 2.0, 3.0]
    results = []
    
    for threshold in thresholds:
        consensus = pharm.consensus_pharm(mols, distance_threshold=threshold)
        results.append(len(consensus))
        print(f"\nThreshold {threshold:.1f} Å: {len(consensus)} consensus features")
        for f in consensus:
            print(f"  {f[0]}: ({f[2]:.3f}, {f[3]:.3f}, {f[4]:.3f})")
    
    # With higher threshold, we should generally have fewer or equal features
    # (features are merged more aggressively)
    print(f"\nFeature counts: {results}")
    print("Note: Higher thresholds may produce fewer features as they merge more aggressively")
    
    assert all(count > 0 for count in results), "All thresholds should produce features"
    print("\n✓ Test passed!")
    return True


def test_with_rdkit_features():
    """Test consensus with RDKit default features."""
    print("\n" + "=" * 60)
    print("Test 4: Using RDKit default features")
    print("=" * 60)
    
    smiles_list = ["CCO", "CCCO"]
    mols = []
    
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        mols.append(mol)
    
    pharm = Pharmacophore()
    consensus = pharm.consensus_pharm(mols, distance_threshold=2.0, features='rdkit')
    
    print(f"\nNumber of consensus features (RDKit): {len(consensus)}")
    print("\nConsensus features:")
    for f in consensus:
        print(f"  {f[0]}: ({f[2]:.3f}, {f[3]:.3f}, {f[4]:.3f})")
    
    assert len(consensus) > 0, "Should have consensus features with RDKit features"
    print("\n✓ Test passed!")
    return True


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("CONSENSUS PHARMACOPHORE TEST SUITE")
    print("=" * 60)
    
    tests = [
        test_same_molecule_conformations,
        test_aligned_molecules,
        test_distance_threshold,
        test_with_rdkit_features,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"\n✗ Test failed with error: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return failed == 0


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)

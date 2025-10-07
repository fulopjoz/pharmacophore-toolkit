#!/usr/bin/env python
"""
Visual demonstration of consensus pharmacophore generation.

This script generates a consensus pharmacophore from neurotransmitter molecules
and saves the output for visualization.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from pharmacophore import Pharmacophore

def main():
    print("Consensus Pharmacophore Demonstration")
    print("=" * 60)
    
    # Define neurotransmitter molecules
    molecules = {
        "Serotonin": "NCCc1c[nH]c2ccc(O)cc12",
        "Dopamine": "NCCc1ccc(O)c(O)c1",
        "Norepinephrine": "NCC(O)c1ccc(O)c(O)c1"
    }
    
    print("\nMolecules:")
    for name, smiles in molecules.items():
        print(f"  {name}: {smiles}")
    
    # Create and prepare molecules
    mols = []
    for name, smiles in molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        mols.append(mol)
    
    # Generate 3D conformations
    print("\nGenerating 3D conformations...")
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    for mol in mols:
        AllChem.EmbedMolecule(mol, ps)
    print("  ✓ 3D conformations generated")
    
    # Align molecules
    print("\nAligning molecules to Serotonin...")
    for i in range(1, len(mols)):
        alignment = rdMolAlign.GetO3A(mols[i], mols[0])
        rmsd = alignment.Align()
        print(f"  {list(molecules.keys())[i]} RMSD: {rmsd:.3f} Å")
    
    # Calculate individual pharmacophores
    print("\nIndividual pharmacophore features:")
    pharm = Pharmacophore()
    for i, (name, mol) in enumerate(zip(molecules.keys(), mols)):
        features = pharm.calc_pharm(mol)
        print(f"  {name}: {len(features)} features")
        
        # Count by type
        type_counts = {}
        for f in features:
            type_counts[f[0]] = type_counts.get(f[0], 0) + 1
        for feat_type, count in sorted(type_counts.items()):
            print(f"    {feat_type}: {count}")
    
    # Generate consensus pharmacophore
    print("\nGenerating consensus pharmacophore...")
    consensus = pharm.consensus_pharm(mols, distance_threshold=2.0)
    print(f"  ✓ Consensus features: {len(consensus)}")
    
    # Analyze consensus features
    print("\nConsensus pharmacophore by type:")
    type_counts = {}
    for f in consensus:
        type_counts[f[0]] = type_counts.get(f[0], 0) + 1
    
    for feat_type, count in sorted(type_counts.items()):
        print(f"  {feat_type}: {count} feature(s)")
    
    print("\nConsensus feature details:")
    for i, f in enumerate(consensus, 1):
        print(f"  {i}. {f[0]:12} at ({f[2]:6.2f}, {f[3]:6.2f}, {f[4]:6.2f})")
    
    # Save for PyMOL visualization
    output_file = 'consensus_demo.pml'
    pharm.output_features(
        feature_list=consensus,
        savepath=output_file,
        sphere_size=0.8,
        transparency=0.3
    )
    print(f"\n✓ PyMOL visualization script saved to: {output_file}")
    
    # Test different thresholds
    print("\nEffect of distance threshold:")
    for threshold in [1.0, 2.0, 3.0]:
        temp_consensus = pharm.consensus_pharm(mols, distance_threshold=threshold)
        print(f"  {threshold:.1f} Å threshold: {len(temp_consensus)} features")
    
    print("\n" + "=" * 60)
    print("Demonstration complete!")
    print("\nKey takeaways:")
    print("  - Consensus pharmacophore identifies common features")
    print("  - Spatial clustering groups nearby features")
    print("  - Centroid calculation produces representative positions")
    print("  - Distance threshold controls clustering sensitivity")
    print("=" * 60)

if __name__ == "__main__":
    main()

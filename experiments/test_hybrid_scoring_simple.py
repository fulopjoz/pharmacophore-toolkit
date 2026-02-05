"""
Hybrid 2D+3D Scoring Validation

Tests hybrid scoring with alpha parameter sweep.
Expected: alpha=0.6 achieves +10-15% improvement over pure 3D.

Literature: Sanders et al. 2012, Moshawih et al. 2024
"""

import sys, time, warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import roc_auc_score  # FIXED: Added missing import

sys.path.insert(0, '.')
from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.hybrid_scoring import HybridScorer
from pharmacophore.screening_metrics import calculate_all_metrics

print("="*80)
print("HYBRID 2D+3D SCORING - VALIDATION")
print("="*80)

# Load CCR2
print("\n📦 Loading dataset...")
refs = [m for m in Chem.SDMolSupplier('tutorials/data/CCR2_reference_ligands.sdf') if m]

def load_molecules_safe(csv_path):
    """Safely load molecules from CSV."""
    df = pd.read_csv(csv_path)
    for col in df.columns:
        col_upper = str(col).upper().strip()
        if 'SMILES' in col_upper and 'MURCKO' not in col_upper:
            smiles_list = df[col].dropna().tolist()
            mols = []
            for smi in smiles_list:
                mol = Chem.MolFromSmiles(str(smi))
                if mol:
                    mols.append(mol)
            return mols
    raise ValueError(f"No SMILES column in {csv_path}")

actives = load_molecules_safe('tutorials/data/actives_ccr2_N75.csv')
decoys = load_molecules_safe('tutorials/data/decoys_ccr2_N500.csv')
print(f"   ✓ {len(refs)} refs, {len(actives)} actives, {len(decoys)} decoys")

# Generate consensus (Phase 5 baseline params)
print("\n🔬 Generating consensus pharmacophore...")
consensus = PharmacophoreConsensus(tolerance=1.5, occurrence_threshold=0.3, linkage='complete')
features = consensus.generate_consensus(refs)
consensus_mol = PharmacophoreToMol.convert(features, enable_color_features=True)
print(f"   ✓ {len(features)} features")

# Prepare molecules with conformers
print("\n⚙️  Preparing test molecules (generating conformers)...")
all_mols = actives + decoys
y_true = [1]*len(actives) + [0]*len(decoys)

prepared_mols = []
valid_y_true = []

for i, mol in enumerate(all_mols):
    try:
        mol_h = Chem.AddHs(mol)
        conf_ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=5, randomSeed=42)
        if not conf_ids:
            AllChem.EmbedMolecule(mol_h, randomSeed=42)
        if mol_h.GetNumConformers() > 0:
            prepared_mols.append(mol_h)
            valid_y_true.append(y_true[i])
    except:
        pass

y_true = valid_y_true  # Use only valid labels

print(f"   ✓ {len(prepared_mols)} molecules prepared ({len(prepared_mols)/len(all_mols)*100:.1f}% success)")

# Test alpha values
alphas = [0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
results = []

print(f"\n{'='*80}")
print(f"ALPHA PARAMETER SWEEP (N={len(alphas)} values)")
print(f"{'='*80}")
print(f"{'Alpha':<10} {'2D Weight':<12} {'3D Weight':<12} {'ROC-AUC':<10} {'Time(s)':<10} {'Improvement'}")
print("-" * 90)

baseline_auc = None

for alpha in alphas:
    print(f"\n[Alpha={alpha:.1f}] Running...", end=" ", flush=True)
    
    # Create scorer
    scorer = HybridScorer(
        reference_mols=refs,
        consensus_mol=consensus_mol,
        alpha=alpha,
        fingerprint_type='morgan',
        max_align_iters=50
    )
    
    # Score molecules
    t0 = time.time()
    scores_full = scorer.score_batch(prepared_mols)  # No slicing needed
    elapsed = time.time() - t0
    
    # Extract hybrid scores from list of dicts
    scores = [s['hybrid'] for s in scores_full]
    
    # Calculate metrics
    auc = roc_auc_score(y_true, scores)
    
    # Determine improvement
    if baseline_auc is None:
        baseline_auc = auc
        improvement_str = "baseline"
    else:
        improvement = (auc - baseline_auc) / baseline_auc * 100
        improvement_str = f"{improvement:+.2f}%"
    
    # Alpha interpretation
    if alpha == 0.0:
        mode_desc = "Pure 3D"
    elif alpha == 1.0:
        mode_desc = "Pure 2D"
    else:
        mode_desc = f"{alpha*100:.0f}%2D+{(1-alpha)*100:.0f}%3D"
    
    print(f"\r{alpha:<10.1f} {alpha:<12.2f} {1-alpha:<12.2f} {auc:<10.4f} {elapsed:<10.1f} {improvement_str:<12} ({mode_desc})")
    
    results.append({
        'alpha': alpha,
        '2d_weight': alpha,
        '3d_weight': 1-alpha,
        'auc': auc,
        'time_sec': elapsed,
        'mode': mode_desc
    })

# Analysis
results_df = pd.DataFrame(results)
best = results_df.loc[results_df['auc'].idxmax()]
pure_3d = results_df[results_df['alpha'] == 0.0].iloc[0]
pure_2d = results_df[results_df['alpha'] == 1.0].iloc[0]

print(f"\n{'='*80}")
print(f"RESULTS SUMMARY")
print(f"{'='*80}")
print(f"Pure 3D (α=0.0):         {pure_3d['auc']:.4f}")
print(f"Pure 2D (α=1.0):         {pure_2d['auc']:.4f}")
print(f"Best Hybrid (α={best['alpha']:.1f}): {best['auc']:.4f}")
print(f"")
print(f"Improvement over 3D:     {(best['auc']-pure_3d['auc'])/pure_3d['auc']*100:+.2f}%")
print(f"Improvement over 2D:     {(best['auc']-pure_2d['auc'])/pure_2d['auc']*100:+.2f}%")

# Literature comparison
print(f"\n📊 Literature Comparison:")
print(f"   Expected hybrid gain: +10-15% (Sanders et al. 2012)")
hybrid_gain = (best['auc'] - pure_3d['auc']) / pure_3d['auc'] * 100
status = "✅ Within range" if 10 <= hybrid_gain <= 15 else "⚠️ Outside range"
print(f"   Observed hybrid gain: {hybrid_gain:+.2f}%")
print(f"   Status: {status}")

# Save results
output_path = 'docs/research/experiments/results/hybrid_scoring_alpha_sweep.csv'
results_df.to_csv(output_path, index=False)
print(f"\n💾 Results saved to {output_path}")

# Recommendation
print(f"\n{'='*80}")
print(f"RECOMMENDATION")
print(f"{'='*80}")
print(f"Optimal α = {best['alpha']:.1f} ({best['mode']})")
print(f"ROC-AUC = {best['auc']:.4f} ({(best['auc']-baseline_auc)/baseline_auc*100:+.2f}% vs baseline)")
print(f"Runtime = {best['time_sec']:.1f}s")
print(f"\nUse HybridScorer(alpha={best['alpha']}) for production screening.")
print(f"{'='*80}\n")

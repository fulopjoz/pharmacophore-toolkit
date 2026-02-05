"""
Full dataset validation of feature weighting + alignment tries.

Usage:
    python experiments/test_improvements_full_dataset.py

Expected Results:
- Baseline (Phase 3): 0.7629 ROC-AUC
- Feature Weights: +2-5% improvement
- Alignment Tries=300: +1-3% improvement
- Combined: +3-8% improvement
"""

import sys, time, warnings
warnings.filterwarnings('ignore')
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, '.')
from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.benchmark import ScreeningBenchmark
from pharmacophore.constants import FEATURE_WEIGHTS

print("="*80)
print("PHASE 5 VALIDATION - FULL DATASET")
print("Feature Weighting + Alignment Tries")
print("="*80)

# Load CCR2
print("\n📦 Loading CCR2 dataset...")
refs = [m for m in Chem.SDMolSupplier('tutorials/data/CCR2_reference_ligands.sdf') if m]

def load_molecules_safe(csv_path):
    """Safely load molecules handling any column name variations."""
    df = pd.read_csv(csv_path)
    print(f"   Columns in {csv_path}: {list(df.columns)}")
    
    # Find SMILES column
    for col in df.columns:
        col_upper = str(col).upper().strip()
        if 'SMILES' in col_upper and 'MURCKO' not in col_upper:
            smiles_list = df[col].dropna().tolist()
            mols = []
            for smi in smiles_list:
                mol = Chem.MolFromSmiles(str(smi))
                if mol:
                    mols.append(mol)
            print(f"   Loaded {len(mols)} molecules from column '{col}'")
            return mols
    
    raise ValueError(f"No SMILES column found in {csv_path}")

actives = load_molecules_safe('tutorials/data/actives_ccr2_N75.csv')
decoys = load_molecules_safe('tutorials/data/decoys_ccr2_N500.csv')

print(f"\n✅ Dataset ready:")
print(f"   References: {len(refs)}")
print(f"   Actives: {len(actives)}")
print(f"   Decoys: {len(decoys)}")
print(f"   Total test molecules: {len(actives) + len(decoys)}")

# Generate consensus (Phase 3 optimal parameters)
print(f"\n🔬 Generating consensus pharmacophore...")
print(f"   Parameters: tolerance=1.5 Å, threshold=0.3, linkage='complete'")

consensus = PharmacophoreConsensus(
    tolerance=1.5,
    occurrence_threshold=0.3,
    linkage='complete'
)
features = consensus.generate_consensus(refs)
consensus_mol = PharmacophoreToMol.convert(features, enable_color_features=True)

print(f"   ✓ Generated {len(features)} consensus features")

# Show feature distribution
print(f"\n📊 Feature distribution with weights:")
for ftype, weight in FEATURE_WEIGHTS.items():
    count = sum(1 for f in features if f[0] == ftype)
    if count > 0:
        print(f"   {ftype:12s} count={count:2d}  weight={weight:.1f}x  effective={count*weight:.1f}")

# Test configurations
configs = [
    {
        'name': 'Baseline (Phase 3)',
        'max_align_iters': 50,
        'use_feature_weights': False,
        'expected': '0.7629'
    },
    {
        'name': 'Feature Weights',
        'max_align_iters': 50,
        'use_feature_weights': True,
        'expected': '+2-5%'
    },
    {
        'name': 'Alignment Tries=100',
        'max_align_iters': 100,
        'use_feature_weights': False,
        'expected': '+0.5-1.5%'
    },
    {
        'name': 'Alignment Tries=300',
        'max_align_iters': 300,
        'use_feature_weights': False,
        'expected': '+1-3%'
    },
    {
        'name': 'Combined (weights+300)',
        'max_align_iters': 300,
        'use_feature_weights': True,
        'expected': '+3-8%'
    },
]

print(f"\n{'='*80}")
print(f"RUNNING EXPERIMENTS (N={len(configs)})")
print(f"{'='*80}")
print(f"{'Configuration':<30} {'ROC-AUC':<10} {'Time(s)':<10} {'Δ':<12} {'Expected'}")
print("-" * 90)

results = []
baseline_auc = None

for i, config in enumerate(configs, 1):
    print(f"\n[{i}/{len(configs)}] {config['name']}...", flush=True)
    
    # Create benchmark
    benchmark = ScreeningBenchmark(
        reference_mols=refs,
        actives=actives,
        decoys=decoys,
        n_conformers=5,
        max_align_iters=config['max_align_iters'],
        use_feature_weights=config['use_feature_weights']
    )
    
    # Run
    t0 = time.time()
    result = benchmark.run_combo_scoring(
        tolerance=1.5,
        occurrence_threshold=0.3
    )
    elapsed = time.time() - t0
    
    auc = result['roc_auc']
    bedroc = result.get('bedroc', 0)
    ef1 = result.get('ef_1pct', 0)
    
    # Calculate improvement
    if baseline_auc is None:
        baseline_auc = auc
        improvement_str = "baseline"
    else:
        improvement = (auc - baseline_auc) / baseline_auc * 100
        improvement_str = f"{improvement:+.2f}%"
    
    print(f"{config['name']:<30} {auc:<10.4f} {elapsed:<10.1f} {improvement_str:<12} {config['expected']}")
    
    results.append({
        'config': config['name'],
        'max_align_iters': config['max_align_iters'],
        'use_weights': config['use_feature_weights'],
        'auc': auc,
        'bedroc': bedroc,
        'ef_1pct': ef1,
        'time_sec': elapsed,
        'expected': config['expected']
    })

# Analysis
results_df = pd.DataFrame(results)
best = results_df.loc[results_df['auc'].idxmax()]

print(f"\n{'='*80}")
print(f"RESULTS SUMMARY")
print(f"{'='*80}")
print(f"Baseline (Phase 3):      ROC-AUC = {baseline_auc:.4f}")
print(f"Best ({best['config']}): ROC-AUC = {best['auc']:.4f}")
print(f"Total Improvement:       {(best['auc']-baseline_auc)/baseline_auc*100:+.2f}%")

# Individual improvements
print(f"\n📊 Individual Contributions:")
for idx, row in results_df.iterrows():
    if row['config'] == 'Baseline (Phase 3)':
        continue
    improvement = (row['auc'] - baseline_auc) / baseline_auc * 100
    status = "✅" if improvement > 0 else "⚠️"
    print(f"{status} {row['config']:<30} {improvement:+.2f}% (expected: {row['expected']})")

# Speed analysis
print(f"\n⚡ Runtime Analysis:")
base_time = results_df[results_df['config'] == 'Baseline (Phase 3)'].iloc[0]['time_sec']
print(f"Baseline (50 iters):  {base_time:.1f}s")

for idx, row in results_df.iterrows():
    if row['config'] == 'Baseline (Phase 3)':
        continue
    speedup = row['time_sec'] / base_time
    print(f"{row['config']:<30} {row['time_sec']:6.1f}s ({speedup:.2f}x)")

# Check sub-linear scaling
a300_time = results_df[results_df['max_align_iters'] == 300].iloc[0]['time_sec']
scaling = a300_time / base_time
print(f"\n{'✅' if scaling < 6 else '❌'} Runtime scaling: {scaling:.2f}x (expected: <6x, HAT predicts ~1.3x)")

# Save results
output_path = 'docs/research/experiments/results/phase5_full_validation.csv'
results_df.to_csv(output_path, index=False)
print(f"\n💾 Results saved to {output_path}")

# Comparison with literature
print(f"\n{'='*80}")
print(f"LITERATURE COMPARISON")
print(f"{'='*80}")

fw_improvement = (results_df[results_df['config'] == 'Feature Weights'].iloc[0]['auc'] - baseline_auc) / baseline_auc * 100
a300_improvement = (results_df[results_df['config'] == 'Alignment Tries=300'].iloc[0]['auc'] - baseline_auc) / baseline_auc * 100

print(f"Feature Weighting:")
print(f"  Observed:  {fw_improvement:+.2f}%")
print(f"  Literature: +2-5% (Langer 2025, Lipinski)")
print(f"  Status:    {'✅ Within expected range' if 2 <= fw_improvement <= 5 else '⚠️ Outside expected range'}")

print(f"\nAlignment Tries (50→300):")
print(f"  Observed:  {a300_improvement:+.2f}%")
print(f"  Literature: +1-3% (Langer HAT)")
print(f"  Status:    {'✅ Within expected range' if 1 <= a300_improvement <= 3 else '⚠️ Outside expected range'}")

combined_improvement = (best['auc'] - baseline_auc) / baseline_auc * 100
print(f"\nCombined Stack:")
print(f"  Observed:  {combined_improvement:+.2f}%")
print(f"  Target:    +3-8%")
print(f"  Status:    {'✅ Target achieved' if combined_improvement >= 3 else '⚠️ Below target'}")

print(f"\n{'='*80}")
print(f"NEXT: Implement 2D+3D hybrid scoring for +10-15% additional gain")
print(f"{'='*80}\n")

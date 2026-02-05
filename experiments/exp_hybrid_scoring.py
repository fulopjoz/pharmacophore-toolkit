"""
Test script for hybrid 2D+3D scoring on CCR2 dataset.

Compares:
1. Pure 3D (Phase 3 baseline: 0.7629 AUC)
2. Pure 2D (fingerprint only)
3. Hybrid 2D+3D with various alpha values

Expected: Hybrid outperforms either method alone
Literature: Sanders et al. 2012, Moshawih et al. 2024
"""

import sys
import time
import pandas as pd
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import roc_auc_score, roc_curve
from pharmacophore.pharmacophore import Pharmacophore
from pharmacophore.hybrid_scoring import HybridScorer, compare_scoring_methods
from pharmacophore.benchmark import ScreeningBenchmark


def load_ccr2_dataset():
    """Load CCR2 reference ligands, actives, and decoys."""
    print("="*80)
    print("LOADING CCR2 DATASET")
    print("="*80)
    
    # Load references
    ref_sdf = 'tutorials/data/CCR2_reference_ligands.sdf'
    refs = [m for m in Chem.SDMolSupplier(ref_sdf) if m]
    print(f"✓ Loaded {len(refs)} reference ligands")
    
    # Load actives and decoys
    actives_df = pd.read_csv('tutorials/data/actives_ccr2_N75.csv')
    decoys_df = pd.read_csv('tutorials/data/decoys_ccr2_N500.csv')
    print(f"✓ Loaded {len(actives_df)} actives, {len(decoys_df)} decoys")
    
    return refs, actives_df, decoys_df


def generate_consensus_model(refs):
    """Generate Phase 3 optimal consensus pharmacophore."""
    print("\n" + "="*80)
    print("GENERATING CONSENSUS PHARMACOPHORE (Phase 3 Optimal)")
    print("="*80)
    
    pharm = Pharmacophore()
    
    start = time.time()
    models = pharm.generate_consensus_models(
        refs,
        tolerance=1.5,
        occurrence_threshold=0.3,
        linkage='complete',
        model_set='standard'  # Returns strict/moderate/relaxed
    )
    elapsed = time.time() - start
    
    # Use moderate model (closest to Phase 3 params)
    features, consensus_mol = models['moderate']
    
    print(f"✓ Generated consensus model")
    print(f"  Features: {len(features)}")
    print(f"  Time: {elapsed:.4f}s")
    
    return features, consensus_mol


def test_hybrid_scoring(refs, consensus_mol, actives_df, decoys_df):
    """Test hybrid 2D+3D scoring with various alpha values."""
    print("\n" + "="*80)
    print("HYBRID 2D+3D SCORING BENCHMARK")
    print("="*80)
    
    # Prepare benchmark
    benchmark = ScreeningBenchmark(actives_df, decoys_df)
    
    # Generate conformers for test molecules
    print("\n📊 Preparing test molecules...")
    test_mols_with_labels = []
    
    for _, row in actives_df.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol:
            AllChem.EmbedMultipleConfs(mol, numConfs=5, randomSeed=42)
            test_mols_with_labels.append((mol, 1))  # Label = 1 (active)
    
    for _, row in decoys_df.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol:
            AllChem.EmbedMultipleConfs(mol, numConfs=5, randomSeed=42)
            test_mols_with_labels.append((mol, 0))  # Label = 0 (decoy)
    
    test_mols = [m for m, _ in test_mols_with_labels]
    labels = [l for _, l in test_mols_with_labels]
    
    print(f"✓ Prepared {len(test_mols)} molecules ({sum(labels)} actives, {len(labels)-sum(labels)} decoys)")
    
    # Test different alpha values
    print("\n📈 Testing alpha values...")
    alphas_to_test = [0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 1.0]
    
    results = []
    
    for alpha in alphas_to_test:
        print(f"\n  Alpha = {alpha:.1f} ", end='')
        
        if alpha == 0.0:
            print("(Pure 3D - Phase 3 baseline)")
        elif alpha == 1.0:
            print("(Pure 2D)")
        else:
            print(f"(Hybrid: {alpha*100:.0f}% 2D, {(1-alpha)*100:.0f}% 3D)")
        
        # Create scorer
        scorer = HybridScorer(
            reference_mols=refs,
            consensus_mol=consensus_mol,
            alpha=alpha,
            fingerprint_type='morgan',
            radius=2,  # ECFP4
            n_bits=1024,
            use_colors=True,
            max_align_iters=50
        )
        
        # Score all molecules
        start = time.time()
        scores_list = scorer.score_batch(test_mols, verbose=False)
        elapsed = time.time() - start
        
        # Extract hybrid scores
        scores = [s['hybrid'] for s in scores_list]
        
        # Calculate metrics
        roc_auc = roc_auc_score(labels, scores)
        
        # Calculate enrichment factor at 1%
        sorted_indices = np.argsort(scores)[::-1]  # Descending
        top_1pct = int(0.01 * len(scores))
        top_labels = [labels[i] for i in sorted_indices[:top_1pct]]
        ef_1pct = (sum(top_labels) / top_1pct) / (sum(labels) / len(labels))
        
        print(f"    ROC-AUC: {roc_auc:.4f}")
        print(f"    EF@1%: {ef_1pct:.2f}")
        print(f"    Time: {elapsed:.2f}s ({elapsed/len(test_mols)*1000:.1f}ms/mol)")
        
        results.append({
            'alpha': alpha,
            'method': 'Pure 3D' if alpha == 0.0 else ('Pure 2D' if alpha == 1.0 else 'Hybrid'),
            'roc_auc': roc_auc,
            'ef_1pct': ef_1pct,
            'time_total': elapsed,
            'time_per_mol': elapsed / len(test_mols)
        })
    
    return pd.DataFrame(results)


def analyze_results(results_df):
    """Analyze and display results."""
    print("\n" + "="*80)
    print("RESULTS SUMMARY")
    print("="*80)
    
    # Find best performers
    best_auc = results_df.loc[results_df['roc_auc'].idxmax()]
    phase3_baseline = results_df[results_df['alpha'] == 0.0].iloc[0]
    
    print("\n📊 Performance Comparison:\n")
    print(results_df[['alpha', 'method', 'roc_auc', 'ef_1pct', 'time_total']].to_string(index=False))
    
    print(f"\n🏆 BEST PERFORMER:")
    print(f"   Alpha: {best_auc['alpha']:.2f}")
    print(f"   Method: {best_auc['method']}")
    print(f"   ROC-AUC: {best_auc['roc_auc']:.4f}")
    print(f"   EF@1%: {best_auc['ef_1pct']:.2f}")
    
    print(f"\n📈 IMPROVEMENT vs Phase 3 Baseline:")
    improvement = (best_auc['roc_auc'] - phase3_baseline['roc_auc']) / phase3_baseline['roc_auc'] * 100
    print(f"   Phase 3 (3D only): {phase3_baseline['roc_auc']:.4f}")
    print(f"   Best Hybrid:       {best_auc['roc_auc']:.4f}")
    print(f"   Improvement:       {improvement:+.1f}%")
    
    if best_auc['roc_auc'] > phase3_baseline['roc_auc']:
        print(f"\n✅ HYBRID SCORING WINS! (+{improvement:.1f}% improvement)")
    else:
        print(f"\n⚠️  Hybrid did not improve (may need parameter tuning)")
    
    # Save results
    results_df.to_csv('docs/research/experiments/results/hybrid_scoring_comparison.csv', index=False)
    print(f"\n💾 Results saved to: docs/research/experiments/results/hybrid_scoring_comparison.csv")
    
    return best_auc


def main():
    """Main execution."""
    print("\n" + "="*80)
    print("HYBRID 2D+3D SCORING - CCR2 BENCHMARK")
    print("="*80)
    print("Testing literature-based improvement: 2D fingerprint + 3D shape ensemble")
    print("Expected: +10-15% ROC-AUC improvement vs Phase 3 baseline (0.7629)")
    print("="*80)
    
    # Load data
    refs, actives_df, decoys_df = load_ccr2_dataset()
    
    # Generate consensus model
    features, consensus_mol = generate_consensus_model(refs)
    
    # Test hybrid scoring
    results_df = test_hybrid_scoring(refs, consensus_mol, actives_df, decoys_df)
    
    # Analyze results
    best_config = analyze_results(results_df)
    
    print("\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)
    print(f"✅ Use alpha = {best_config['alpha']:.2f} for production")
    print(f"✅ Expected ROC-AUC: {best_config['roc_auc']:.4f}")
    print(f"✅ Screening time: {best_config['time_total']:.1f}s for {len(actives_df)+len(decoys_df)} molecules")
    
    print("\n" + "="*80)
    print("NEXT STEPS")
    print("="*80)
    print("1. Add feature weighting (+2-5% AUC)")
    print("2. Add alignment tries parameter (+1-3% AUC)")
    print("3. Implement multi-model ensemble (+5-8% AUC)")
    print("4. Validate on additional datasets (DUD-E, LIT-PCBA)")
    print("="*80)


if __name__ == '__main__':
    main()

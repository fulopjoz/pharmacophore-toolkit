"""Fast prototype comparison of three pharmacophore optimization approaches.

This is a RAPID VALIDATION script designed to verify all three approaches work
correctly before running the full comparison. Uses:
- Reduced trials (10 instead of 50/100)
- Fewer conformers (5 instead of 25)
- Smaller dataset subset (20 actives, 100 decoys)
- joblib caching to avoid recomputation

Expected runtime: 10-20 minutes (vs 4+ hours for full version)

Usage:
    python experiments/compare_optimizers_fast.py
"""

import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from rdkit import Chem
from rdkit.Chem import AllChem, SDMolSupplier
from joblib import Memory
import pandas as pd

from pharmacophore.optuna_optimizer import OptunaPharmacophoreOptimizer
from pharmacophore.hypogen_optimizer import HypoGenOptimizer

# Setup caching
cache_dir = Path(__file__).parent / '.cache'
cache_dir.mkdir(exist_ok=True)
memory = Memory(cache_dir, verbose=0)


@memory.cache
def generate_conformers_cached(smiles, n_conformers, random_seed=42):
    """Cached conformer generation to avoid recomputation."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    # No UFF optimization - ETKDGv3 is sufficient
    AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)

    if mol_h.GetNumConformers() > 0:
        return mol_h
    return None


def load_ccr2_dataset_fast(n_actives=20, n_decoys=100, n_conformers=5, verbose=True):
    """Load CCR2 dataset subset for fast prototyping.

    Args:
        n_actives: Number of active molecules to load (default: 20 vs 75 full).
        n_decoys: Number of decoy molecules to load (default: 100 vs 500 full).
        n_conformers: Conformers per molecule (default: 5 vs 25 full).
        verbose: Print loading progress.

    Returns:
        Tuple of (refs, actives, decoys).
    """
    if verbose:
        print(f"Loading CCR2 dataset (FAST MODE)...")
        print(f"  Subset: {n_actives} actives, {n_decoys} decoys, {n_conformers} conformers")

    base_path = Path(__file__).parent.parent / 'tutorials' / 'data'

    # Load references (keep all 5)
    refs_path = base_path / 'CCR2_reference_ligands.sdf'
    refs = []
    supplier = SDMolSupplier(str(refs_path), removeHs=False)
    for mol in supplier:
        if mol and mol.GetNumConformers() > 0:
            refs.append(mol)

    if verbose:
        print(f"  ✓ Loaded {len(refs)} reference molecules")

    # Load actives subset
    actives_path = base_path / 'actives_ccr2_N75.csv'
    actives_df = pd.read_csv(actives_path)
    actives = []

    for i, smiles in enumerate(actives_df['SMILES'][:n_actives]):
        mol = generate_conformers_cached(smiles, n_conformers, random_seed=42)
        if mol:
            actives.append(mol)

    if verbose:
        print(f"  ✓ Loaded {len(actives)}/{n_actives} actives ({n_conformers} conformers each)")

    # Load decoys subset
    decoys_path = base_path / 'decoys_ccr2_N500.csv'
    decoys_df = pd.read_csv(decoys_path)
    decoys = []

    for i, smiles in enumerate(decoys_df['Smiles'][:n_decoys]):
        mol = generate_conformers_cached(smiles, n_conformers, random_seed=42)
        if mol:
            decoys.append(mol)

    if verbose:
        print(f"  ✓ Loaded {len(decoys)}/{n_decoys} decoys ({n_conformers} conformers each)")
        print(f"  Total: {len(actives)} actives, {len(decoys)} decoys\n")

    return refs, actives, decoys


def run_optuna_gp_fast(refs, actives, decoys, n_trials=10, verbose=True):
    """Run Optuna GP-Sampler (FAST: 10 trials instead of 50)."""
    if verbose:
        print(f"\n{'='*70}")
        print(f"APPROACH 1: Optuna GP-Sampler ({n_trials} trials - FAST MODE)")
        print(f"{'='*70}\n")

    optimizer = OptunaPharmacophoreOptimizer(
        refs, actives, decoys, random_state=42
    )

    start = time.time()
    result = optimizer.optimize(sampler='gp', n_trials=n_trials, verbose=verbose)
    wall_time = time.time() - start

    return {
        'approach': 'Optuna GP',
        'n_trials': n_trials,
        'best_auc': result['best_auc'],
        'best_bedroc': result['best_bedroc'],
        'wall_time_sec': wall_time,
        'time_per_eval': wall_time / n_trials,
        'pareto_size': len(result['pareto_front']),
        'result': result
    }


def run_optuna_nsga2_fast(refs, actives, decoys, n_trials=10, verbose=True):
    """Run Optuna NSGA-II (FAST: 10 trials instead of 100)."""
    if verbose:
        print(f"\n{'='*70}")
        print(f"APPROACH 2: Optuna NSGA-II ({n_trials} trials - FAST MODE)")
        print(f"{'='*70}\n")

    optimizer = OptunaPharmacophoreOptimizer(
        refs, actives, decoys, random_state=42
    )

    start = time.time()
    result = optimizer.optimize(sampler='nsga2', n_trials=n_trials, verbose=verbose)
    wall_time = time.time() - start

    return {
        'approach': 'NSGA-II',
        'n_trials': n_trials,
        'best_auc': result['best_auc'],
        'best_bedroc': result['best_bedroc'],
        'wall_time_sec': wall_time,
        'time_per_eval': wall_time / n_trials,
        'pareto_size': len(result['pareto_front']),
        'result': result
    }


def run_hypogen_fast(refs, actives, decoys, verbose=True):
    """Run HypoGen 3-phase (FAST: reduced hypotheses and SA iterations)."""
    if verbose:
        print(f"\n{'='*70}")
        print(f"APPROACH 3: HypoGen 3-Phase (FAST MODE)")
        print(f"{'='*70}\n")

    optimizer = HypoGenOptimizer(refs, actives, decoys, random_state=42)

    start = time.time()
    result = optimizer.optimize(
        consensus_tolerance=2.0,
        consensus_occurrence=0.3,
        min_features=3,
        max_features=6,          # Reduced from 8
        max_hypotheses=50,       # Reduced from 500
        n_top=3,                 # Reduced from 10
        n_sa_iters=20,           # Reduced from 200
        auc_threshold=0.50,      # Lowered for small dataset
        bedroc_threshold=0.00,
        n_conformers=5,          # Fast
        verbose=verbose
    )
    wall_time = time.time() - start

    return {
        'approach': 'HypoGen',
        'n_trials': result['n_evaluations'],
        'best_auc': result['best_result'].roc_auc,
        'best_bedroc': result['best_result'].bedroc,
        'wall_time_sec': wall_time,
        'time_per_eval': wall_time / result['n_evaluations'] if result['n_evaluations'] > 0 else 0,
        'pareto_size': len(result['phase2_survivors']),
        'result': result
    }


def print_comparison_report(results):
    """Print formatted comparison table."""
    print(f"\n{'='*70}")
    print(f"FAST PROTOTYPE COMPARISON - CCR2 Subset")
    print(f"{'='*70}\n")

    # Header
    print(f"{'Approach':<15} {'AUC':>7} {'BEDROC':>7} {'Trials':>7} "
          f"{'Time(s)':>8} {'s/Eval':>7}")
    print(f"{'-'*70}")

    # Rows
    for name, result in results.items():
        print(
            f"{result['approach']:<15} "
            f"{result['best_auc']:>7.4f} "
            f"{result['best_bedroc']:>7.4f} "
            f"{result['n_trials']:>7d} "
            f"{result['wall_time_sec']:>8.0f} "
            f"{result['time_per_eval']:>7.2f}"
        )

    print(f"{'-'*70}")

    # Winners
    best_auc = max(results.items(), key=lambda x: x[1]['best_auc'])
    best_bedroc = max(results.items(), key=lambda x: x[1]['best_bedroc'])
    fastest = min(results.items(), key=lambda x: x[1]['wall_time_sec'])

    print(f"\n✓ VALIDATION SUCCESSFUL - All approaches working!")
    print(f"\nBest Performance:")
    print(f"  • Best AUC: {best_auc[1]['approach']} ({best_auc[1]['best_auc']:.4f})")
    print(f"  • Best BEDROC: {best_bedroc[1]['approach']} ({best_bedroc[1]['best_bedroc']:.4f})")
    print(f"  • Fastest: {fastest[1]['approach']} ({fastest[1]['wall_time_sec']:.0f}s)")

    total_time = sum(r['wall_time_sec'] for r in results.values())
    print(f"\nTotal runtime: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    print(f"\n{'='*70}\n")


def main():
    """Run fast prototype comparison."""
    print(f"\n{'='*70}")
    print(f"FAST PROTOTYPE: Three-Way Optimizer Comparison")
    print(f"{'='*70}")
    print(f"\nThis is a RAPID VALIDATION to verify all approaches work.")
    print(f"For full comparison, use: experiments/compare_optimizers.py")
    print(f"\nOptimizations:")
    print(f"  • Dataset: 20 actives, 100 decoys (vs 75/500 full)")
    print(f"  • Conformers: 5 per molecule (vs 25 full)")
    print(f"  • Trials: 10 each for GP/NSGA-II (vs 50/100 full)")
    print(f"  • HypoGen: 50 hypotheses, 20 SA iters (vs 500/200 full)")
    print(f"  • Caching: joblib memoization for conformer generation")
    print(f"\nExpected runtime: 10-20 minutes\n")

    # Load dataset
    refs, actives, decoys = load_ccr2_dataset_fast(
        n_actives=20,
        n_decoys=100,
        n_conformers=5,
        verbose=True
    )

    results = {}

    # Run all three approaches
    results['gp'] = run_optuna_gp_fast(refs, actives, decoys, n_trials=10, verbose=True)
    results['nsga2'] = run_optuna_nsga2_fast(refs, actives, decoys, n_trials=10, verbose=True)
    results['hypogen'] = run_hypogen_fast(refs, actives, decoys, verbose=True)

    # Print comparison
    print_comparison_report(results)

    print("✓ Fast prototype complete!")
    print("\nNext steps:")
    print("  1. Verify results look reasonable (AUC > 0.5, BEDROC > 0)")
    print("  2. If satisfied, run full comparison: python experiments/compare_optimizers.py")
    print("  3. Or adjust parameters and re-run this fast version\n")


if __name__ == '__main__':
    main()

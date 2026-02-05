"""End-to-end comparison of three pharmacophore optimization approaches.

Compares:
1. Optuna GP-Sampler: Gaussian Process Bayesian optimization (50 trials)
2. NSGA-II: Evolutionary multi-objective optimization (100 trials)
3. HypoGen: Constructive-subtractive-optimization 3-phase approach

Dataset: CCR2 receptor ligands (5 references, 75 actives, 500 decoys)

Usage:
    python experiments/compare_optimizers.py

Output:
    - Console: Formatted comparison report
    - File: experiments/results/optimizer_comparison_TIMESTAMP.csv
"""

import os
import sys
import time
import pandas as pd
from datetime import datetime
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rdkit import Chem
from rdkit.Chem import AllChem, SDMolSupplier

from pharmacophore.optuna_optimizer import OptunaPharmacophoreOptimizer
from pharmacophore.hypogen_optimizer import HypoGenOptimizer


def load_ccr2_dataset(n_conformers=25, verbose=True):
    """Load CCR2 dataset with references, actives, and decoys.

    Args:
        n_conformers: Number of conformers to generate for actives/decoys.
        verbose: Print loading progress.

    Returns:
        Tuple of (refs, actives, decoys) as RDKit molecules with conformers.
    """
    if verbose:
        print("Loading CCR2 dataset...")

    base_path = Path(__file__).parent.parent / 'tutorials' / 'data'

    # Load references from SDF
    refs_path = base_path / 'CCR2_reference_ligands.sdf'
    refs = []
    supplier = SDMolSupplier(str(refs_path), removeHs=False)
    for mol in supplier:
        if mol and mol.GetNumConformers() > 0:
            refs.append(mol)

    if verbose:
        print(f"  Loaded {len(refs)} reference molecules")

    # Load actives from CSV
    actives_path = base_path / 'actives_ccr2_N75.csv'
    actives_df = pd.read_csv(actives_path)
    actives = []

    for smiles in actives_df['SMILES']:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol_h = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)
            if mol_h.GetNumConformers() > 0:
                actives.append(mol_h)

    if verbose:
        print(f"  Loaded {len(actives)} active molecules ({n_conformers} conformers each)")

    # Load decoys from CSV
    decoys_path = base_path / 'decoys_ccr2_N500.csv'
    decoys_df = pd.read_csv(decoys_path)
    decoys = []

    for smiles in decoys_df['Smiles']:  # Note: capital 'S' for decoys
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol_h = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)
            if mol_h.GetNumConformers() > 0:
                decoys.append(mol_h)

    if verbose:
        print(f"  Loaded {len(decoys)} decoy molecules ({n_conformers} conformers each)")
        print(f"  Total dataset: {len(actives)} actives, {len(decoys)} decoys\n")

    return refs, actives, decoys


def run_optuna_gp(refs, actives, decoys, n_trials=50, verbose=True):
    """Run Optuna GP-Sampler optimization (Approach 1).

    Args:
        refs: Reference molecules.
        actives: Active molecules.
        decoys: Decoy molecules.
        n_trials: Number of optimization trials.
        verbose: Print progress.

    Returns:
        Result dictionary.
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"APPROACH 1: Optuna GP-Sampler ({n_trials} trials)")
        print(f"{'='*70}\n")

    optimizer = OptunaPharmacophoreOptimizer(
        refs, actives, decoys, random_state=42
    )

    result = optimizer.optimize(
        sampler='gp',
        n_trials=n_trials,
        verbose=verbose
    )

    return {
        'approach': 'Optuna GP',
        'sampler': 'gp',
        'n_trials': n_trials,
        'best_auc': result['best_auc'],
        'best_bedroc': result['bedroc'],
        'ef_1': result['pareto_front'][0]['ef_1'] if result['pareto_front'] else 0.0,
        'ef_5': result['pareto_front'][0]['ef_5'] if result['pareto_front'] else 0.0,
        'ef_10': result['pareto_front'][0]['ef_10'] if result['pareto_front'] else 0.0,
        'wall_time_sec': result['wall_time_sec'],
        'n_evaluations': n_trials,
        'time_per_eval': result['wall_time_sec'] / n_trials,
        'pareto_size': len(result['pareto_front']),
        'best_params': result['best_auc_params'],
        'raw_result': result
    }


def run_optuna_nsga2(refs, actives, decoys, n_trials=100, verbose=True):
    """Run Optuna NSGA-II optimization (Approach 2).

    Args:
        refs: Reference molecules.
        actives: Active molecules.
        decoys: Decoy molecules.
        n_trials: Number of optimization trials.
        verbose: Print progress.

    Returns:
        Result dictionary.
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"APPROACH 2: Optuna NSGA-II ({n_trials} trials)")
        print(f"{'='*70}\n")

    optimizer = OptunaPharmacophoreOptimizer(
        refs, actives, decoys, random_state=42
    )

    result = optimizer.optimize(
        sampler='nsga2',
        n_trials=n_trials,
        verbose=verbose
    )

    return {
        'approach': 'NSGA-II',
        'sampler': 'nsga2',
        'n_trials': n_trials,
        'best_auc': result['best_auc'],
        'best_bedroc': result['best_bedroc'],
        'ef_1': result['pareto_front'][0]['ef_1'] if result['pareto_front'] else 0.0,
        'ef_5': result['pareto_front'][0]['ef_5'] if result['pareto_front'] else 0.0,
        'ef_10': result['pareto_front'][0]['ef_10'] if result['pareto_front'] else 0.0,
        'wall_time_sec': result['wall_time_sec'],
        'n_evaluations': n_trials,
        'time_per_eval': result['wall_time_sec'] / n_trials,
        'pareto_size': len(result['pareto_front']),
        'best_params': result['best_auc_params'],
        'raw_result': result
    }


def run_hypogen(refs, actives, decoys, verbose=True):
    """Run HypoGen 3-phase optimization (Approach 3).

    Args:
        refs: Reference molecules.
        actives: Active molecules.
        decoys: Decoy molecules.
        verbose: Print progress.

    Returns:
        Result dictionary.
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"APPROACH 3: HypoGen 3-Phase")
        print(f"{'='*70}\n")

    optimizer = HypoGenOptimizer(
        refs, actives, decoys, random_state=42
    )

    result = optimizer.optimize(
        consensus_tolerance=2.0,
        consensus_occurrence=0.3,
        min_features=3,
        max_features=8,
        max_hypotheses=500,
        n_top=10,
        n_sa_iters=200,
        auc_threshold=0.55,
        bedroc_threshold=0.05,
        n_conformers=25,
        verbose=verbose
    )

    return {
        'approach': 'HypoGen',
        'sampler': '3-phase',
        'n_trials': result['n_evaluations'],
        'best_auc': result['best_result'].roc_auc,
        'best_bedroc': result['best_result'].bedroc,
        'ef_1': result['best_result'].ef_1,
        'ef_5': result['best_result'].ef_5,
        'ef_10': result['best_result'].ef_10,
        'wall_time_sec': result['total_time'],
        'n_evaluations': result['n_evaluations'],
        'time_per_eval': result['total_time'] / result['n_evaluations'],
        'pareto_size': len(result['phase2_survivors']),
        'best_params': {
            'n_features': result['best_result'].n_features,
            'cost': result['best_cost']
        },
        'raw_result': result
    }


def print_comparison_report(results):
    """Print formatted comparison table.

    Args:
        results: Dictionary mapping approach names to result dicts.
    """
    print(f"\n{'='*80}")
    print(f"OPTIMIZER COMPARISON REPORT - CCR2 Dataset")
    print(f"{'='*80}\n")

    # Header
    print(f"{'Approach':<15} {'AUC':>7} {'BEDROC':>7} {'EF@1%':>7} {'EF@5%':>7} "
          f"{'Wall(s)':>8} {'#Evals':>7} {'s/Eval':>7}")
    print(f"{'-'*80}")

    # Rows
    for name, result in results.items():
        print(
            f"{result['approach']:<15} "
            f"{result['best_auc']:>7.4f} "
            f"{result['best_bedroc']:>7.4f} "
            f"{result['ef_1']:>7.1f} "
            f"{result['ef_5']:>7.1f} "
            f"{result['wall_time_sec']:>8.0f} "
            f"{result['n_evaluations']:>7d} "
            f"{result['time_per_eval']:>7.2f}"
        )

    print(f"{'-'*80}")

    # Winners
    best_auc = max(results.items(), key=lambda x: x[1]['best_auc'])
    best_bedroc = max(results.items(), key=lambda x: x[1]['best_bedroc'])
    fastest = min(results.items(), key=lambda x: x[1]['wall_time_sec'])
    most_efficient = min(results.items(), key=lambda x: x[1]['time_per_eval'])

    print(f"\nWINNERS:")
    print(f"  Best AUC: {best_auc[1]['approach']} ({best_auc[1]['best_auc']:.4f})")
    print(f"  Best BEDROC: {best_bedroc[1]['approach']} ({best_bedroc[1]['best_bedroc']:.4f})")
    print(f"  Fastest: {fastest[1]['approach']} ({fastest[1]['wall_time_sec']:.0f}s)")
    print(f"  Most Efficient: {most_efficient[1]['approach']} ({most_efficient[1]['time_per_eval']:.2f}s/eval)")
    print(f"\n{'='*80}\n")


def save_results(results, output_dir='experiments/results'):
    """Save results to CSV file.

    Args:
        results: Dictionary of results.
        output_dir: Output directory for results.
    """
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)

    # Generate timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_path = Path(output_dir) / f'optimizer_comparison_{timestamp}.csv'

    # Build DataFrame
    rows = []
    for name, result in results.items():
        row = {
            'approach': result['approach'],
            'sampler': result['sampler'],
            'n_trials': result['n_trials'],
            'best_auc': result['best_auc'],
            'best_bedroc': result['best_bedroc'],
            'ef_1': result['ef_1'],
            'ef_5': result['ef_5'],
            'ef_10': result['ef_10'],
            'wall_time_sec': result['wall_time_sec'],
            'n_evaluations': result['n_evaluations'],
            'time_per_eval': result['time_per_eval'],
            'pareto_size': result['pareto_size']
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)

    print(f"Results saved to: {output_path}")


def main():
    """Run full optimizer comparison."""
    print(f"\n{'='*80}")
    print(f"THREE-WAY OPTIMIZER COMPARISON")
    print(f"{'='*80}\n")

    # Load dataset
    refs, actives, decoys = load_ccr2_dataset(n_conformers=25, verbose=True)

    results = {}

    # Run Approach 1: Optuna GP (50 trials)
    results['gp'] = run_optuna_gp(refs, actives, decoys, n_trials=50, verbose=True)

    # Run Approach 2: NSGA-II (100 trials)
    results['nsga2'] = run_optuna_nsga2(refs, actives, decoys, n_trials=100, verbose=True)

    # Run Approach 3: HypoGen
    results['hypogen'] = run_hypogen(refs, actives, decoys, verbose=True)

    # Print comparison report
    print_comparison_report(results)

    # Save results
    save_results(results)

    # Final summary
    print("\nComparison complete!")
    print(f"Total wall time: {sum(r['wall_time_sec'] for r in results.values()):.1f}s")
    print(f"Total evaluations: {sum(r['n_evaluations'] for r in results.values())}")


if __name__ == '__main__':
    main()

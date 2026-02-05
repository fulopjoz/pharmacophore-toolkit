#!/usr/bin/env python
"""
Phase 4: Bayesian Optimization for Fine-Tuning DBSCAN Consensus Parameters

Uses Gaussian Process regression with Expected Improvement acquisition
to sequentially explore parameter space and find global optimum.

Advantages over grid search:
- Sample-efficient: ~30 evaluations vs 100+ for grid
- Explores intelligently based on past results
- Balances exploration vs exploitation
- Provides uncertainty estimates

Search Space (narrowed from Phase 3):
- tolerance: [1.3, 1.7] Å
- threshold: [0.25, 0.35]
- min_samples: [1, 4] (DBSCAN density parameter)

Target: ROC-AUC > 0.80
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from datetime import datetime
import time
import warnings
from typing import Dict, List, Any, Tuple

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from skopt import gp_minimize
from skopt.space import Real, Integer
from skopt.utils import use_named_args
from skopt import dump, load
from skopt.plots import plot_convergence, plot_objective, plot_evaluations
import matplotlib.pyplot as plt

from rdkit import Chem
from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.screening_metrics import calculate_all_metrics
from pharmacophore.clustering_algorithms import cluster_dbscan

from experiments.parameter_sweep import (
    load_ccr2_dataset, score_molecule
)
from experiments.experiment_logger import ExperimentLogger


class BayesianOptimizedConsensus:
    """DBSCAN consensus with Bayesian-optimized parameters."""
    
    def __init__(
        self,
        tolerance: float,
        occurrence_threshold: float,
        min_samples: int = 2
    ):
        """Initialize with DBSCAN-specific parameters.
        
        Args:
            tolerance: DBSCAN eps parameter (Å)
            occurrence_threshold: Minimum fraction of molecules
            min_samples: DBSCAN min_samples parameter
        """
        self.tolerance = tolerance
        self.occurrence_threshold = occurrence_threshold
        self.min_samples = min_samples
    
    def generate_consensus(
        self,
        mols: List[Chem.Mol]
    ) -> List[List]:
        """Generate consensus using DBSCAN clustering.
        
        Args:
            mols: List of aligned RDKit molecules
        
        Returns:
            List of consensus features [type, (), x, y, z]
        """
        from pharmacophore.pharmacophore import Pharmacophore
        
        # Extract features from each molecule
        all_features_by_type = {}
        
        for mol in mols:
            pharm = Pharmacophore()
            features = pharm.calc_pharm(mol)
            
            for feat in features:
                feat_type = feat[0]
                coords = np.array(feat[2:5])
                
                if feat_type not in all_features_by_type:
                    all_features_by_type[feat_type] = []
                
                all_features_by_type[feat_type].append(coords)
        
        # Cluster each feature type with DBSCAN
        consensus_features = []
        
        for feat_type, coords_list in all_features_by_type.items():
            if not coords_list:
                continue
            
            coords_array = np.array(coords_list)
            
            # Apply DBSCAN clustering
            centroids = cluster_dbscan(
                coords=coords_array,
                tolerance=self.tolerance,
                occurrence_threshold=self.occurrence_threshold,
                n_molecules=len(mols)
            )
            
            # Convert centroids to feature format
            for centroid in centroids:
                consensus_features.append([
                    feat_type,
                    (),  # Empty tuple (no atom indices for consensus)
                    float(centroid[0]),
                    float(centroid[1]),
                    float(centroid[2])
                ])
        
        return consensus_features


def evaluate_parameters(
    tolerance: float,
    threshold: float,
    min_samples: int,
    verbose: bool = False
) -> Dict[str, Any]:
    """Evaluate consensus parameters on CCR2 dataset.
    
    Args:
        tolerance: DBSCAN eps (Å)
        threshold: Occurrence threshold
        min_samples: DBSCAN min_samples
        verbose: Print progress
    
    Returns:
        Dict with performance metrics
    """
    min_samples = int(round(min_samples))  # Ensure integer
    
    if verbose:
        print(f"\nEvaluating: tolerance={tolerance:.4f}, threshold={threshold:.4f}, min_samples={min_samples}")
    
    # Load dataset
    ref_mols, actives, decoys = load_ccr2_dataset()
    
    try:
        # Generate consensus
        start_time = time.time()
        consensus = BayesianOptimizedConsensus(
            tolerance=tolerance,
            occurrence_threshold=threshold,
            min_samples=min_samples
        )
        features = consensus.generate_consensus(ref_mols)
        consensus_time = time.time() - start_time
        
        if verbose:
            print(f"  Features: {len(features)}, Time: {consensus_time:.4f}s")
        
        if not features:
            return {
                'roc_auc': 0.5,
                'ef_1': 0.0,
                'bedroc': 0.0,
                'n_features': 0,
                'consensus_time': consensus_time,
                'screening_time': 0,
                'error': 'No features'
            }
        
        # Convert to mol
        pharm_mol = PharmacophoreToMol.convert(features, enable_color_features=True)
        
        # Screen molecules
        start_screen = time.time()
        all_mols = actives + decoys
        y_true = [1] * len(actives) + [0] * len(decoys)
        
        scores = []
        for mol in all_mols:
            score = score_molecule(mol, pharm_mol, n_conformers=5)
            scores.append(score)
        
        screening_time = time.time() - start_screen
        
        # Calculate metrics
        metrics = calculate_all_metrics(y_true, scores)
        
        if verbose:
            print(f"  ROC-AUC: {metrics['roc_auc']:.4f}, EF@1%: {metrics['ef_1']:.2f}")
        
        return {
            **metrics,
            'n_features': len(features),
            'consensus_time': consensus_time,
            'screening_time': screening_time,
            'total_time': consensus_time + screening_time
        }
        
    except Exception as e:
        if verbose:
            print(f"  ERROR: {str(e)}")
        
        return {
            'roc_auc': 0.5,
            'ef_1': 0.0,
            'bedroc': 0.0,
            'n_features': 0,
            'consensus_time': 0,
            'screening_time': 0,
            'error': str(e)
        }


def run_bayesian_optimization(
    n_calls: int = 30,
    n_initial_points: int = 5,
    random_state: int = 42,
    output_dir: Path = None
) -> Tuple[Any, pd.DataFrame]:
    """Run Bayesian optimization for DBSCAN consensus parameters.
    
    Args:
        n_calls: Total number of evaluations
        n_initial_points: Number of random initial points
        random_state: Random seed for reproducibility
        output_dir: Directory to save results
    
    Returns:
        Tuple of (optimization result, results dataframe)
    """
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "docs/research/experiments/results"
    
    # Initialize logger
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = ExperimentLogger(
        experiment_name=f"bayesian_optimization_{timestamp}",
        output_dir=output_dir
    )
    
    print(f"\n{'='*80}")
    print(f"PHASE 4: BAYESIAN OPTIMIZATION")
    print(f"{'='*80}")
    print(f"Algorithm: DBSCAN consensus")
    print(f"Total evaluations: {n_calls}")
    print(f"Initial random points: {n_initial_points}")
    print(f"Target: ROC-AUC > 0.80")
    print(f"{'='*80}\n")
    
    # Define search space
    space = [
        Real(1.3, 1.7, name='tolerance'),
        Real(0.25, 0.35, name='threshold'),
        Integer(1, 4, name='min_samples')
    ]
    
    # Store results for logging
    results_list = []
    iteration = [0]  # Use list to allow modification in nested function
    
    # Define objective function
    @use_named_args(space)
    def objective(tolerance, threshold, min_samples):
        """Objective function to minimize (negative AUC)."""
        iteration[0] += 1
        
        print(f"\n[Iteration {iteration[0]}/{n_calls}]")
        
        metrics = evaluate_parameters(
            tolerance=tolerance,
            threshold=threshold,
            min_samples=min_samples,
            verbose=True
        )
        
        # Log result
        params = {
            'tolerance': tolerance,
            'threshold': threshold,
            'min_samples': min_samples
        }
        
        logger.log_run(
            run_id=iteration[0],
            parameters=params,
            metrics=metrics,
            metadata={'phase': 'Phase4_BayesianOptimization'}
        )
        
        results_list.append({
            'iteration': iteration[0],
            'tolerance': tolerance,
            'threshold': threshold,
            'min_samples': min_samples,
            **metrics
        })
        
        # Return negative AUC (we minimize, so negate to maximize)
        return -metrics['roc_auc']
    
    # Run optimization
    print("\nStarting Bayesian optimization...")
    
    result = gp_minimize(
        func=objective,
        dimensions=space,
        n_calls=n_calls,
        n_initial_points=n_initial_points,
        acq_func='EI',  # Expected Improvement
        random_state=random_state,
        verbose=False
    )
    
    # Save results
    logger.save_csv()
    logger.save_json()
    logger.generate_summary_report()
    
    df_results = pd.DataFrame(results_list)
    
    # Find best configuration
    best_idx = df_results['roc_auc'].idxmax()
    best = df_results.loc[best_idx]
    
    print(f"\n{'='*80}")
    print(f"BAYESIAN OPTIMIZATION COMPLETE")
    print(f"{'='*80}")
    print(f"\nBest Configuration:")
    print(f"  Tolerance: {best['tolerance']:.4f} Å")
    print(f"  Threshold: {best['threshold']:.4f}")
    print(f"  Min Samples: {int(best['min_samples'])}")
    print(f"\nPerformance:")
    print(f"  ROC-AUC: {best['roc_auc']:.4f}")
    print(f"  EF@1%: {best['ef_1']:.2f}")
    print(f"  BEDROC: {best['bedroc']:.4f}")
    print(f"  Features: {int(best['n_features'])}")
    
    # Compare to Phase 3
    phase3_auc = 0.7629
    improvement = (best['roc_auc'] - phase3_auc) / phase3_auc * 100
    print(f"\nImprovement over Phase 3:")
    print(f"  Phase 3 DBSCAN: {phase3_auc:.4f}")
    print(f"  Phase 4 Optimized: {best['roc_auc']:.4f}")
    print(f"  Change: {improvement:+.2f}%")
    
    if best['roc_auc'] >= 0.80:
        print(f"\n🎉 TARGET ACHIEVED: ROC-AUC ≥ 0.80!")
    
    print(f"{'='*80}\n")
    
    # Generate plots
    print("Generating optimization plots...")
    plot_optimization_results(result, df_results, output_dir / "plots")
    
    # Save optimization object
    dump(result, str(output_dir / f"bayesian_opt_result_{timestamp}.pkl"))
    
    return result, df_results


def plot_optimization_results(
    result: Any,
    df_results: pd.DataFrame,
    output_dir: Path
):
    """Generate visualization plots for Bayesian optimization.
    
    Args:
        result: skopt optimization result
        df_results: DataFrame with results
        output_dir: Directory to save plots
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Convergence plot
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_convergence(result, ax=ax)
    ax.set_title('Bayesian Optimization Convergence', fontsize=14, fontweight='bold')
    ax.set_xlabel('Number of Evaluations', fontsize=12)
    ax.set_ylabel('Best ROC-AUC Found', fontsize=12)
    # Invert y-axis since we minimized negative AUC
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_yticklabels([f"{abs(float(t.get_text())):.3f}" for t in ax.get_yticklabels()])
    plt.tight_layout()
    plt.savefig(output_dir / "bayesian_convergence.png", dpi=300)
    plt.close()
    
    # 2. Parameter progression
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    axes[0].plot(df_results['iteration'], df_results['tolerance'], 'o-')
    axes[0].set_ylabel('Tolerance (Å)', fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].axhline(y=1.5, color='r', linestyle='--', alpha=0.5, label='Phase 3 value')
    axes[0].legend()
    
    axes[1].plot(df_results['iteration'], df_results['threshold'], 'o-')
    axes[1].set_ylabel('Threshold', fontsize=12, fontweight='bold')
    axes[1].grid(True, alpha=0.3)
    axes[1].axhline(y=0.3, color='r', linestyle='--', alpha=0.5, label='Phase 3 value')
    axes[1].legend()
    
    axes[2].plot(df_results['iteration'], df_results['min_samples'], 'o-')
    axes[2].set_ylabel('Min Samples', fontsize=12, fontweight='bold')
    axes[2].set_xlabel('Iteration', fontsize=12, fontweight='bold')
    axes[2].grid(True, alpha=0.3)
    
    plt.suptitle('Parameter Exploration Over Iterations', fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(output_dir / "parameter_progression.png", dpi=300)
    plt.close()
    
    # 3. ROC-AUC progression
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(df_results['iteration'], df_results['roc_auc'], 'o-', label='Observed AUC')
    
    # Running best
    running_best = df_results['roc_auc'].cummax()
    ax.plot(df_results['iteration'], running_best, 'r-', linewidth=2, label='Best So Far')
    
    ax.axhline(y=0.7629, color='g', linestyle='--', alpha=0.7, label='Phase 3 DBSCAN')
    ax.axhline(y=0.80, color='orange', linestyle='--', alpha=0.7, label='Target (0.80)')
    
    ax.set_xlabel('Iteration', fontsize=12, fontweight='bold')
    ax.set_ylabel('ROC-AUC', fontsize=12, fontweight='bold')
    ax.set_title('Optimization Progress', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "optimization_progress.png", dpi=300)
    plt.close()
    
    print(f"✓ Plots saved to: {output_dir}")


def run_confirmation_experiments(
    tolerance: float,
    threshold: float,
    min_samples: int,
    n_replicates: int = 5
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Run confirmation experiments at optimal parameters.
    
    Args:
        tolerance: Optimal tolerance
        threshold: Optimal threshold
        min_samples: Optimal min_samples
        n_replicates: Number of replicate runs
    
    Returns:
        Tuple of (results dataframe, statistics dict)
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = ExperimentLogger(
        experiment_name=f"bayesian_confirmation_{timestamp}",
        output_dir=Path(__file__).parent.parent / "docs/research/experiments/results"
    )
    
    print(f"\n{'='*80}")
    print(f"CONFIRMATION EXPERIMENTS")
    print(f"{'='*80}")
    print(f"Parameters: tolerance={tolerance:.4f}, threshold={threshold:.4f}, min_samples={min_samples}")
    print(f"Replicates: {n_replicates}")
    print(f"{'='*80}\n")
    
    results = []
    
    for rep in range(1, n_replicates + 1):
        print(f"[{rep}/{n_replicates}] Running confirmation...")
        
        metrics = evaluate_parameters(
            tolerance=tolerance,
            threshold=threshold,
            min_samples=min_samples,
            verbose=True
        )
        
        metrics['replicate'] = rep
        metrics['tolerance'] = tolerance
        metrics['threshold'] = threshold
        metrics['min_samples'] = min_samples
        
        logger.log_run(
            run_id=rep,
            parameters={'tolerance': tolerance, 'threshold': threshold, 'min_samples': min_samples},
            metrics=metrics,
            metadata={'replicate': rep, 'phase': 'Phase4_Confirmation'}
        )
        
        results.append(metrics)
    
    logger.save_csv()
    logger.save_json()
    logger.generate_summary_report()
    
    df_confirm = pd.DataFrame(results)
    
    # Statistics
    mean_auc = df_confirm['roc_auc'].mean()
    std_auc = df_confirm['roc_auc'].std()
    ci_95 = 1.96 * std_auc / np.sqrt(n_replicates)
    
    print(f"\n{'='*80}")
    print(f"CONFIRMATION RESULTS")
    print(f"{'='*80}")
    print(f"Mean ROC-AUC: {mean_auc:.4f} ± {std_auc:.4f}")
    print(f"95% CI: [{mean_auc - ci_95:.4f}, {mean_auc + ci_95:.4f}]")
    print(f"Range: {df_confirm['roc_auc'].min():.4f} - {df_confirm['roc_auc'].max():.4f}")
    print(f"{'='*80}\n")
    
    stats = {
        'mean_auc': mean_auc,
        'std_auc': std_auc,
        'ci_95': ci_95,
        'min_auc': df_confirm['roc_auc'].min(),
        'max_auc': df_confirm['roc_auc'].max()
    }
    
    return df_confirm, stats


def main():
    """Main execution pipeline for Phase 4."""
    
    # Run Bayesian optimization
    result, df_results = run_bayesian_optimization(
        n_calls=30,
        n_initial_points=5,
        random_state=42
    )
    
    # Get best parameters
    best_idx = df_results['roc_auc'].idxmax()
    best = df_results.loc[best_idx]
    
    # Run confirmation experiments
    df_confirm, stats = run_confirmation_experiments(
        tolerance=best['tolerance'],
        threshold=best['threshold'],
        min_samples=int(best['min_samples']),
        n_replicates=5
    )
    
    return result, df_results, df_confirm, stats


if __name__ == "__main__":
    result, df_results, df_confirm, stats = main()

#!/usr/bin/env python
"""
Phase 3: Algorithm Comparison for Consensus Pharmacophore Generation

Compare different clustering algorithms:
1. Hierarchical (current default)
2. K-means
3. DBSCAN
4. Grid-based binning

Test at optimal parameters from Phase 2:
- Tolerance: 1.5 Å
- Threshold: 0.3
- Dataset: CCR2 (5 refs, 74 actives, 499 decoys)

Metrics:
- Accuracy: ROC-AUC, EF@1%, BEDROC
- Speed: Consensus time, screening time, total time
- Features: Number of consensus features generated
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from datetime import datetime
import time
from typing import Dict, List, Any, Tuple
import warnings

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol

from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.clustering_algorithms import (
    cluster_features_generic, ALGORITHM_INFO
)
from pharmacophore.screening_metrics import calculate_all_metrics
from experiments.parameter_sweep import (
    load_ccr2_dataset, score_molecule
)
from experiments.experiment_logger import ExperimentLogger


class PharmacophoreConsensusWithAlgorithm(PharmacophoreConsensus):
    """Extended consensus class that supports multiple clustering algorithms."""
    
    def __init__(
        self,
        tolerance: float = 2.0,
        occurrence_threshold: float = 0.5,
        linkage: str = 'average',
        clustering_method: str = 'hierarchical'
    ):
        """Initialize with clustering algorithm selection.
        
        Args:
            tolerance: Distance threshold (Å)
            occurrence_threshold: Minimum fraction (0.0-1.0)
            linkage: For hierarchical only
            clustering_method: 'hierarchical', 'kmeans', 'dbscan', 'grid'
        """
        super().__init__(tolerance, occurrence_threshold, linkage)
        self.clustering_method = clustering_method
    
    def _cluster_features(
        self,
        coordinates: List[np.ndarray],
        mol_indices: np.ndarray
    ) -> Tuple[np.ndarray, Dict[int, List[int]]]:
        """Cluster features using selected algorithm.
        
        Override parent method to use different clustering algorithms.
        """
        if not coordinates:
            return np.array([]), {}
        
        # Convert to numpy array
        coords_array = np.array(coordinates)
        n_molecules = len(set(mol_indices))
        
        # Use alternative clustering if not hierarchical
        if self.clustering_method != 'hierarchical':
            from pharmacophore.clustering_algorithms import cluster_features_generic
            
            # Get centroids from alternative algorithm
            centroids = cluster_features_generic(
                coords=coords_array,
                tolerance=self.tolerance,
                occurrence_threshold=self.occurrence_threshold,
                n_molecules=n_molecules,
                method=self.clustering_method,
                linkage=self.linkage
            )
            
            if not centroids:
                return np.array([]), {}
            
            # Map features to nearest centroid (pseudo-labeling)
            from scipy.spatial.distance import cdist
            
            centroid_array = np.array(centroids)
            distances = cdist(coords_array, centroid_array, metric='euclidean')
            labels = np.argmin(distances, axis=1)
            
            # Map clusters to molecule indices
            cluster_to_mols = {}
            for cluster_id, mol_idx in zip(labels, mol_indices):
                cluster_id = int(cluster_id)
                mol_idx = int(mol_idx)
                
                if cluster_id not in cluster_to_mols:
                    cluster_to_mols[cluster_id] = []
                
                if mol_idx not in cluster_to_mols[cluster_id]:
                    cluster_to_mols[cluster_id].append(mol_idx)
            
            return labels, cluster_to_mols
        
        # Use default hierarchical clustering
        return super()._cluster_features(coordinates, mol_indices)


def evaluate_algorithm(
    algorithm: str,
    tolerance: float = 1.5,
    threshold: float = 0.3,
    linkage: str = 'complete',
    verbose: bool = True
) -> Dict[str, Any]:
    """Evaluate a clustering algorithm on CCR2 dataset.
    
    Args:
        algorithm: 'hierarchical', 'kmeans', 'dbscan', 'grid'
        tolerance: Spatial tolerance (Å)
        threshold: Occurrence threshold
        linkage: Linkage method (for hierarchical only)
        verbose: Print progress
    
    Returns:
        Dict with metrics and timing information
    """
    if verbose:
        print(f"\n{'='*80}")
        print(f"Algorithm: {algorithm.upper()}")
        print(f"{'='*80}")
    
    # Load dataset
    ref_mols, actives, decoys = load_ccr2_dataset()
    
    # Generate consensus with timing
    start_time = time.time()
    
    consensus = PharmacophoreConsensusWithAlgorithm(
        tolerance=tolerance,
        occurrence_threshold=threshold,
        linkage=linkage,
        clustering_method=algorithm
    )
    
    try:
        features = consensus.generate_consensus(ref_mols)
        consensus_time = time.time() - start_time
        
        if verbose:
            print(f"  Consensus generation: {consensus_time:.4f} s")
            print(f"  Features generated: {len(features)}")
        
        if not features:
            return {
                'algorithm': algorithm,
                'n_features': 0,
                'consensus_time': consensus_time,
                'screening_time': 0,
                'total_time': consensus_time,
                'roc_auc': 0.5,
                'ef_1': 0.0,
                'bedroc': 0.0,
                'error': 'No features generated'
            }
        
        # Convert to mol
        pharm_mol = PharmacophoreToMol.convert(
            features, enable_color_features=True
        )
        
        if verbose:
            print(f"  Pharmacophore mol created: {pharm_mol.GetNumAtoms()} atoms")
        
        # Score molecules with timing
        start_screening = time.time()
        all_mols = actives + decoys
        y_true = [1] * len(actives) + [0] * len(decoys)
        
        scores = []
        for mol in all_mols:
            score = score_molecule(mol, pharm_mol, n_conformers=5)
            scores.append(score)
        
        screening_time = time.time() - start_screening
        total_time = consensus_time + screening_time
        
        if verbose:
            print(f"  Screening time: {screening_time:.4f} s")
            print(f"  Total time: {total_time:.4f} s")
        
        # Calculate metrics
        metrics = calculate_all_metrics(y_true, scores)
        
        if verbose:
            print(f"  ROC-AUC: {metrics['roc_auc']:.4f}")
            print(f"  EF@1%: {metrics['ef_1']:.2f}")
            print(f"  BEDROC: {metrics['bedroc']:.4f}")
        
        return {
            'algorithm': algorithm,
            'n_features': len(features),
            'consensus_time': consensus_time,
            'screening_time': screening_time,
            'total_time': total_time,
            'time_per_mol_ms': screening_time / len(all_mols) * 1000,
            **metrics
        }
        
    except Exception as e:
        if verbose:
            print(f"  ERROR: {str(e)}")
        
        return {
            'algorithm': algorithm,
            'n_features': 0,
            'consensus_time': 0,
            'screening_time': 0,
            'total_time': 0,
            'roc_auc': 0.5,
            'ef_1': 0.0,
            'bedroc': 0.0,
            'error': str(e)
        }


def run_algorithm_comparison(
    algorithms: List[str] = None,
    n_replicates: int = 3,
    tolerance: float = 1.5,
    threshold: float = 0.3,
    linkage: str = 'complete'
) -> pd.DataFrame:
    """Run comparative benchmark of clustering algorithms.
    
    Args:
        algorithms: List of algorithms to test (default: all)
        n_replicates: Number of replicate runs per algorithm
        tolerance: Spatial tolerance (Å)
        threshold: Occurrence threshold
        linkage: Linkage method for hierarchical
    
    Returns:
        DataFrame with comparison results
    """
    if algorithms is None:
        algorithms = ['hierarchical', 'kmeans', 'grid']

    # Initialize logger
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = ExperimentLogger(
        experiment_name=f"algorithm_comparison_{timestamp}",
        output_dir=Path(__file__).parent.parent / "docs/research/experiments/results"
    )
    
    print(f"\n{'='*80}")
    print(f"PHASE 3: ALGORITHM COMPARISON")
    print(f"{'='*80}")
    print(f"Algorithms: {', '.join(algorithms)}")
    print(f"Replicates: {n_replicates}")
    print(f"Parameters: tolerance={tolerance}Å, threshold={threshold}, linkage={linkage}")
    print(f"{'='*80}\n")
    
    results = []
    run_id = 1
    
    for algorithm in algorithms:
        print(f"\n{'='*80}")
        print(f"Testing Algorithm: {ALGORITHM_INFO[algorithm]['name']}")
        print(f"Complexity: {ALGORITHM_INFO[algorithm]['complexity']}")
        print(f"{'='*80}")
        
        for rep in range(1, n_replicates + 1):
            print(f"\n[Replicate {rep}/{n_replicates}]")
            
            metrics = evaluate_algorithm(
                algorithm=algorithm,
                tolerance=tolerance,
                threshold=threshold,
                linkage=linkage,
                verbose=True
            )
            
            metrics['replicate'] = rep
            metrics['tolerance'] = tolerance
            metrics['threshold'] = threshold
            metrics['linkage'] = linkage
            
            # Log result
            logger.log_run(
                run_id=run_id,
                parameters={
                    'algorithm': algorithm,
                    'tolerance': tolerance,
                    'threshold': threshold,
                    'linkage': linkage
                },
                metrics=metrics,
                metadata={'replicate': rep, 'phase': 'Phase3_AlgorithmComparison'}
            )
            
            results.append(metrics)
            run_id += 1
    
    # Save results
    logger.save_csv()
    logger.save_json()
    logger.generate_summary_report()
    
    df_results = pd.DataFrame(results)
    
    # Generate comparison summary
    print(f"\n{'='*80}")
    print(f"ALGORITHM COMPARISON SUMMARY")
    print(f"{'='*80}\n")
    
    # Group by algorithm and calculate statistics
    summary = df_results.groupby('algorithm').agg({
        'roc_auc': ['mean', 'std'],
        'ef_1': ['mean', 'std'],
        'bedroc': ['mean', 'std'],
        'n_features': ['mean', 'std'],
        'consensus_time': ['mean', 'std'],
        'screening_time': ['mean', 'std'],
        'total_time': ['mean', 'std']
    }).round(4)
    
    print("Performance Comparison:")
    print(summary)
    
    # Rank algorithms
    print(f"\n{'='*80}")
    print("RANKINGS")
    print(f"{'='*80}\n")
    
    algorithm_means = df_results.groupby('algorithm').mean()
    
    print("By Accuracy (ROC-AUC):")
    ranked_auc = algorithm_means.sort_values('roc_auc', ascending=False)
    for i, (alg, row) in enumerate(ranked_auc.iterrows(), 1):
        print(f"  {i}. {alg:12s}: {row['roc_auc']:.4f}")
    
    print("\nBy Speed (Total Time):")
    ranked_speed = algorithm_means.sort_values('total_time', ascending=True)
    for i, (alg, row) in enumerate(ranked_speed.iterrows(), 1):
        print(f"  {i}. {alg:12s}: {row['total_time']:.4f} s")
    
    print("\nBy Efficiency (AUC per second):")
    algorithm_means['efficiency'] = algorithm_means['roc_auc'] / algorithm_means['total_time']
    ranked_eff = algorithm_means.sort_values('efficiency', ascending=False)
    for i, (alg, row) in enumerate(ranked_eff.iterrows(), 1):
        print(f"  {i}. {alg:12s}: {row['efficiency']:.4f}")
    
    # Save detailed comparison
    comparison_path = logger.results_dir / f"PHASE3_ALGORITHM_COMPARISON_{datetime.now().strftime('%Y%m%d')}.md"
    with open(comparison_path, 'w') as f:
        f.write("# Phase 3: Algorithm Comparison Results\n\n")
        f.write(f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"**Parameters**: tolerance={tolerance}Å, threshold={threshold}, linkage={linkage}\n")
        f.write(f"**Replicates**: {n_replicates} per algorithm\n\n")
        
        f.write("## Summary Statistics\n\n")
        f.write(summary.to_markdown())
        f.write("\n\n")
        
        f.write("## Rankings\n\n")
        f.write("### By Accuracy (ROC-AUC)\n")
        for i, (alg, row) in enumerate(ranked_auc.iterrows(), 1):
            f.write(f"{i}. **{alg}**: {row['roc_auc']:.4f}\n")
        
        f.write("\n### By Speed (Total Time)\n")
        for i, (alg, row) in enumerate(ranked_speed.iterrows(), 1):
            f.write(f"{i}. **{alg}**: {row['total_time']:.4f} s\n")
        
        f.write("\n### By Efficiency (AUC/second)\n")
        for i, (alg, row) in enumerate(ranked_eff.iterrows(), 1):
            f.write(f"{i}. **{alg}**: {row['efficiency']:.4f}\n")
        
        f.write("\n## Algorithm Information\n\n")
        for alg in algorithms:
            info = ALGORITHM_INFO[alg]
            f.write(f"### {info['name']}\n\n")
            f.write(f"- **Complexity**: {info['complexity']}\n")
            f.write(f"- **Advantages**: {', '.join(info['advantages'])}\n")
            f.write(f"- **Disadvantages**: {', '.join(info['disadvantages'])}\n")
            f.write(f"- **Best for**: {info['best_for']}\n\n")
    
    print(f"\n{'='*80}")
    print(f"✓ Algorithm comparison complete!")
    print(f"✓ Results saved to: {logger.results_dir}")
    print(f"✓ Summary saved to: {comparison_path}")
    print(f"{'='*80}\n")
    
    return df_results


def main():
    """Main execution pipeline for Phase 3."""
    
    # Run comparison at optimal parameters from Phase 2
    df_results = run_algorithm_comparison(
        algorithms=['hierarchical', 'kmeans', 'grid'],
        n_replicates=3,
        tolerance=1.5,
        threshold=0.3,
        linkage='complete'
    )
    
    return df_results


if __name__ == "__main__":
    results = main()

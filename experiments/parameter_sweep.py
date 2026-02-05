#!/usr/bin/env python3
"""Parameter sweep experiments for consensus pharmacophore optimization.

This script runs systematic experiments to find optimal parameters for
consensus pharmacophore generation based on screening performance.

Experiments:
1. Tolerance Sweep: Test spatial tolerance 0.8-2.5 Å
2. Threshold Sweep: Test occurrence threshold 0.3-1.0
3. Linkage Comparison: Compare average/complete/single/ward
4. Combined Optimization: Grid search over best parameters

Usage:
    python experiments/parameter_sweep.py --experiment tolerance
    python experiments/parameter_sweep.py --experiment threshold
    python experiments/parameter_sweep.py --experiment linkage
    python experiments/parameter_sweep.py --experiment all
"""

import argparse
import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Any
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

# Import pharmacophore toolkit
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.screening_metrics import calculate_all_metrics
from rdkit.Chem.rdShapeAlign import AlignMol


def load_ccr2_dataset() -> Tuple[List[Chem.Mol], List[Chem.Mol], List[Chem.Mol]]:
    """Load CCR2 reference, actives, and decoys from tutorial data."""
    data_dir = Path(__file__).parent.parent / "tutorials" / "data"
    
    # Load reference ligands
    ref_path = data_dir / "CCR2_reference_ligands.sdf"
    ref_mols = [m for m in Chem.SDMolSupplier(str(ref_path)) if m]
    
    # Load actives
    actives_path = data_dir / "actives_ccr2_N75.csv"
    actives_df = pd.read_csv(actives_path)
    actives = []
    for smi in actives_df['SMILES']:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            if mol.GetNumConformers() > 0:
                actives.append(mol)
    
    # Load decoys
    decoys_path = data_dir / "decoys_ccr2_N500.csv"
    decoys_df = pd.read_csv(decoys_path)
    decoys = []
    for smi in decoys_df['Smiles']:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            if mol.GetNumConformers() > 0:
                decoys.append(mol)
    
    return ref_mols, actives, decoys


def score_molecule(mol: Chem.Mol, pharm_mol: Chem.Mol, n_conformers: int = 5) -> float:
    """Score a molecule against pharmacophore model."""
    try:
        mol_h = Chem.AddHs(mol)
        conf_ids = AllChem.EmbedMultipleConfs(
            mol_h, numConfs=n_conformers, randomSeed=42
        )
        
        if not conf_ids:
            AllChem.EmbedMolecule(mol_h, randomSeed=42)
            conf_ids = [0] if mol_h.GetNumConformers() > 0 else []
        
        best_score = 0.0
        for conf_id in conf_ids:
            try:
                conf_mol = Chem.Mol(mol_h)
                shape, color = AlignMol(
                    ref=pharm_mol,
                    probe=conf_mol,
                    useColors=True,
                    opt_param=0.5
                )
                combo = shape + color
                if combo > best_score:
                    best_score = combo
            except Exception:
                continue
        
        return best_score
    except Exception:
        return 0.0


def evaluate_model(
    ref_mols: List[Chem.Mol],
    actives: List[Chem.Mol],
    decoys: List[Chem.Mol],
    tolerance: float,
    threshold: float,
    linkage: str = 'average'
) -> Dict[str, Any]:
    """Evaluate a consensus model with given parameters."""
    
    # Generate consensus
    start_time = time.time()
    consensus = PharmacophoreConsensus(
        tolerance=tolerance,
        occurrence_threshold=threshold,
        linkage=linkage
    )
    features = consensus.generate_consensus(ref_mols)
    consensus_time = time.time() - start_time
    
    if not features:
        return {
            'n_features': 0,
            'consensus_time': consensus_time,
            'screening_time': 0,
            'roc_auc': 0.5,
            'ef_1': 0.0,
            'bedroc': 0.0,
            'error': 'No features generated'
        }
    
    # Convert to mol
    pharm_mol = PharmacophoreToMol.convert(
        features, enable_color_features=True
    )
    
    # Score molecules
    start_time = time.time()
    all_mols = actives + decoys
    y_true = [1] * len(actives) + [0] * len(decoys)
    
    scores = []
    for mol in all_mols:
        score = score_molecule(mol, pharm_mol)
        scores.append(score)
    
    screening_time = time.time() - start_time
    
    # Calculate metrics
    metrics = calculate_all_metrics(y_true, scores)
    
    return {
        'n_features': len(features),
        'consensus_time': consensus_time,
        'screening_time': screening_time,
        'time_per_mol_ms': screening_time / len(all_mols) * 1000,
        **metrics
    }


def run_tolerance_sweep(
    ref_mols: List[Chem.Mol],
    actives: List[Chem.Mol],
    decoys: List[Chem.Mol],
    threshold: float = 0.5,
    linkage: str = 'average'
) -> pd.DataFrame:
    """Sweep tolerance parameter and record results."""
    
    tolerances = [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5]
    results = []
    
    for tol in tolerances:
        print(f"  Testing tolerance={tol}...", end=" ", flush=True)
        result = evaluate_model(
            ref_mols, actives, decoys, tol, threshold, linkage
        )
        result['tolerance'] = tol
        result['threshold'] = threshold
        result['linkage'] = linkage
        results.append(result)
        print(f"AUC={result['roc_auc']:.3f}, EF@1%={result['ef_1']:.1f}, "
              f"features={result['n_features']}")
    
    return pd.DataFrame(results)


def run_threshold_sweep(
    ref_mols: List[Chem.Mol],
    actives: List[Chem.Mol],
    decoys: List[Chem.Mol],
    tolerance: float = 1.5,
    linkage: str = 'average'
) -> pd.DataFrame:
    """Sweep occurrence threshold and record results."""
    
    thresholds = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    results = []
    
    for thresh in thresholds:
        print(f"  Testing threshold={thresh}...", end=" ", flush=True)
        result = evaluate_model(
            ref_mols, actives, decoys, tolerance, thresh, linkage
        )
        result['tolerance'] = tolerance
        result['threshold'] = thresh
        result['linkage'] = linkage
        results.append(result)
        print(f"AUC={result['roc_auc']:.3f}, EF@1%={result['ef_1']:.1f}, "
              f"features={result['n_features']}")
    
    return pd.DataFrame(results)


def run_linkage_comparison(
    ref_mols: List[Chem.Mol],
    actives: List[Chem.Mol],
    decoys: List[Chem.Mol],
    tolerance: float = 1.5,
    threshold: float = 0.5
) -> pd.DataFrame:
    """Compare different linkage methods."""
    
    linkages = ['average', 'complete', 'single', 'ward']
    results = []
    
    for link in linkages:
        print(f"  Testing linkage={link}...", end=" ", flush=True)
        result = evaluate_model(
            ref_mols, actives, decoys, tolerance, threshold, link
        )
        result['tolerance'] = tolerance
        result['threshold'] = threshold
        result['linkage'] = link
        results.append(result)
        print(f"AUC={result['roc_auc']:.3f}, EF@1%={result['ef_1']:.1f}, "
              f"features={result['n_features']}")
    
    return pd.DataFrame(results)


def run_grid_search(
    ref_mols: List[Chem.Mol],
    actives: List[Chem.Mol],
    decoys: List[Chem.Mol]
) -> pd.DataFrame:
    """Grid search over tolerance and threshold."""
    
    tolerances = [1.0, 1.3, 1.5, 1.8, 2.0]
    thresholds = [0.4, 0.5, 0.6, 0.7, 0.8]
    results = []
    
    total = len(tolerances) * len(thresholds)
    count = 0
    
    for tol in tolerances:
        for thresh in thresholds:
            count += 1
            print(f"  [{count}/{total}] tol={tol}, thresh={thresh}...", 
                  end=" ", flush=True)
            result = evaluate_model(
                ref_mols, actives, decoys, tol, thresh, 'average'
            )
            result['tolerance'] = tol
            result['threshold'] = thresh
            result['linkage'] = 'average'
            results.append(result)
            print(f"AUC={result['roc_auc']:.3f}")
    
    return pd.DataFrame(results)


def save_results(df: pd.DataFrame, experiment_name: str, output_dir: Path):
    """Save experiment results to files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # CSV
    csv_path = output_dir / f"{experiment_name}_{timestamp}.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")
    
    # Summary
    summary_path = output_dir / f"{experiment_name}_{timestamp}_summary.md"
    
    with open(summary_path, 'w') as f:
        f.write(f"# {experiment_name.replace('_', ' ').title()} Results\n\n")
        f.write(f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        
        # Best result
        best_idx = df['roc_auc'].idxmax()
        best = df.loc[best_idx]
        f.write("## Best Configuration\n\n")
        f.write(f"- **ROC-AUC**: {best['roc_auc']:.3f}\n")
        f.write(f"- **EF@1%**: {best['ef_1']:.1f}\n")
        f.write(f"- **BEDROC**: {best['bedroc']:.3f}\n")
        f.write(f"- **Features**: {best['n_features']}\n")
        if 'tolerance' in best:
            f.write(f"- **Tolerance**: {best['tolerance']} Å\n")
        if 'threshold' in best:
            f.write(f"- **Threshold**: {best['threshold']}\n")
        if 'linkage' in best:
            f.write(f"- **Linkage**: {best['linkage']}\n")
        
        f.write("\n## Full Results\n\n")
        f.write(df.to_markdown(index=False))
    
    print(f"Summary saved to: {summary_path}")
    
    return csv_path, summary_path


def main():
    parser = argparse.ArgumentParser(
        description="Parameter sweep experiments for consensus pharmacophore"
    )
    parser.add_argument(
        '--experiment',
        choices=['tolerance', 'threshold', 'linkage', 'grid', 'all'],
        default='all',
        help='Which experiment to run'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path(__file__).parent.parent / "docs" / "research" / "experiments",
        help='Output directory for results'
    )
    args = parser.parse_args()
    
    print("=" * 60)
    print("Consensus Pharmacophore Parameter Sweep Experiments")
    print("=" * 60)
    
    # Load data
    print("\nLoading CCR2 dataset...")
    try:
        ref_mols, actives, decoys = load_ccr2_dataset()
        print(f"  References: {len(ref_mols)}")
        print(f"  Actives: {len(actives)}")
        print(f"  Decoys: {len(decoys)}")
    except FileNotFoundError as e:
        print(f"Error: Could not load dataset: {e}")
        print("Make sure CCR2 data exists in tutorials/data/")
        return 1
    
    experiments_to_run = []
    if args.experiment == 'all':
        experiments_to_run = ['tolerance', 'threshold', 'linkage', 'grid']
    else:
        experiments_to_run = [args.experiment]
    
    for exp in experiments_to_run:
        print(f"\n{'=' * 60}")
        print(f"Running: {exp.upper()} experiment")
        print("=" * 60)
        
        if exp == 'tolerance':
            df = run_tolerance_sweep(ref_mols, actives, decoys)
        elif exp == 'threshold':
            df = run_threshold_sweep(ref_mols, actives, decoys)
        elif exp == 'linkage':
            df = run_linkage_comparison(ref_mols, actives, decoys)
        elif exp == 'grid':
            df = run_grid_search(ref_mols, actives, decoys)
        
        save_results(df, f"{exp}_sweep", args.output_dir)
    
    print("\n" + "=" * 60)
    print("All experiments complete!")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    exit(main())

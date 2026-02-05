"""CLI script for running combinatorial pharmacophore optimization.

Usage:
    python experiments/combinatorial_search.py \\
        --refs tutorials/data/CCR2_reference_ligands.sdf \\
        --actives tutorials/data/actives_ccr2_N75.csv \\
        --decoys tutorials/data/decoys_ccr2_N500.csv \\
        --n-trials 30 \\
        --output experiments/results/combinatorial/
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from pharmacophore.combinatorial_optimizer import (
    CombinatorialPharmacophoreOptimizer,
    CombinatorialResult,
)


def load_sdf(path: str):
    """Load molecules from SDF file."""
    supplier = Chem.SDMolSupplier(path, removeHs=False)
    return [mol for mol in supplier if mol is not None]


def load_csv_smiles(path: str, smiles_col: str = 'SMILES'):
    """Load molecules from CSV with SMILES column."""
    df = pd.read_csv(path)
    # Try common SMILES column names
    for col in [smiles_col, 'smiles', 'Smiles', 'SMILES', 'canonical_smiles']:
        if col in df.columns:
            mols = []
            for smi in df[col]:
                mol = Chem.MolFromSmiles(str(smi))
                if mol is not None:
                    mols.append(mol)
            return mols
    raise ValueError(f"No SMILES column found in {path}. Columns: {list(df.columns)}")


def save_results(result: CombinatorialResult, output_dir: Path):
    """Save optimization results to files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Best config as JSON
    config_path = output_dir / 'best_config.json'
    serializable_config = {}
    for k, v in result.best_config.items():
        if isinstance(v, (tuple, list)):
            serializable_config[k] = list(v)
        elif isinstance(v, np.integer):
            serializable_config[k] = int(v)
        elif isinstance(v, np.floating):
            serializable_config[k] = float(v)
        else:
            serializable_config[k] = v

    with open(config_path, 'w') as f:
        json.dump({
            'best_config': serializable_config,
            'best_metrics': result.best_metrics,
            'total_evaluations': result.total_evaluations,
            'wall_time_sec': result.wall_time_sec,
        }, f, indent=2)

    # Top-k configs as CSV
    if result.top_k_configs:
        rows = []
        for entry in result.top_k_configs:
            row = {'roc_auc': entry['roc_auc']}
            if entry.get('config'):
                for k, v in entry['config'].items():
                    if isinstance(v, (tuple, list)):
                        row[k] = str(v)
                    else:
                        row[k] = v
            rows.append(row)
        df = pd.DataFrame(rows)
        df.to_csv(output_dir / 'top_configs.csv', index=False)

    # Feature importance
    if result.feature_importance:
        imp_path = output_dir / 'feature_importance.json'
        with open(imp_path, 'w') as f:
            json.dump(result.feature_importance, f, indent=2)

    print(f"\nResults saved to {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description='Combinatorial pharmacophore optimization'
    )
    parser.add_argument(
        '--refs', required=True,
        help='Reference molecules SDF file'
    )
    parser.add_argument(
        '--actives', required=True,
        help='Active molecules (SDF or CSV with SMILES)'
    )
    parser.add_argument(
        '--decoys', required=True,
        help='Decoy molecules (SDF or CSV with SMILES)'
    )
    parser.add_argument(
        '--n-trials', type=int, default=30,
        help='Optuna trials per discrete combo (default: 30)'
    )
    parser.add_argument(
        '--min-subset-size', type=int, default=3,
        help='Minimum reference subset size (default: 3)'
    )
    parser.add_argument(
        '--n-conformers', type=int, default=10,
        help='Fixed conformer count for 3D scoring (default: 10)'
    )
    parser.add_argument(
        '--scoring-modes', type=str, nargs='+',
        default=None,
        help='Scoring modes to test (default: reference hybrid)'
    )
    parser.add_argument(
        '--no-minimize', action='store_true',
        help='Skip MMFF minimization option (faster)'
    )
    parser.add_argument(
        '--no-multi-fidelity', action='store_true',
        help='Disable multi-fidelity pruning'
    )
    parser.add_argument(
        '--output', type=str,
        default='experiments/results/combinatorial',
        help='Output directory'
    )
    parser.add_argument(
        '--seed', type=int, default=42,
        help='Random seed (default: 42)'
    )

    args = parser.parse_args()

    # Load molecules
    print("Loading molecules...")
    refs = load_sdf(args.refs)
    print(f"  References: {len(refs)}")

    actives_path = args.actives
    if actives_path.endswith('.sdf'):
        actives = load_sdf(actives_path)
    else:
        actives = load_csv_smiles(actives_path)
    print(f"  Actives: {len(actives)}")

    decoys_path = args.decoys
    if decoys_path.endswith('.sdf'):
        decoys = load_sdf(decoys_path)
    else:
        decoys = load_csv_smiles(decoys_path)
    print(f"  Decoys: {len(decoys)}")

    # Run optimization
    optimizer = CombinatorialPharmacophoreOptimizer(
        reference_mols=refs,
        actives=actives,
        decoys=decoys,
        random_state=args.seed
    )

    minimize_options = [False] if args.no_minimize else [False, True]

    result = optimizer.optimize(
        n_trials=args.n_trials,
        n_conformers=args.n_conformers,
        multi_fidelity=not args.no_multi_fidelity,
        min_subset_size=args.min_subset_size,
        scoring_modes=args.scoring_modes,
        minimize_options=minimize_options,
        verbose=True
    )

    # Save results
    save_results(result, Path(args.output))


if __name__ == '__main__':
    main()

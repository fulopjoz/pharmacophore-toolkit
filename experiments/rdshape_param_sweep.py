#!/usr/bin/env python
"""
rdShapeAlign Parameter Sweep for Consensus Pharmacophore Optimization

This script systematically tests rdShapeAlign parameters to find optimal settings
for discriminating actives from decoys using consensus pharmacophore models.

Usage:
    python experiments/rdshape_param_sweep.py

Output:
    - experiments/results/rdshape_param_sweep_results.csv
    - Console summary of best parameters

Key Parameters Tested:
    - opt_param: Balance between shape (1.0) and color (0.0) optimization
    - max_preiters: Phase 1 iterations on all starting poses
    - max_postiters: Phase 2 iterations on best poses
    - tolerance: Consensus clustering distance threshold
    - occurrence_threshold: Minimum feature frequency

Reference:
    See docs/RDSHAPEALIGN_OPTIMIZATION_PLAN.md for full context
"""

import sys
import os
import time
import itertools
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol, PrepareConformer
from sklearn.metrics import roc_auc_score

from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.screening_metrics import bedroc, enrichment_factor


# =============================================================================
# Configuration
# =============================================================================

# Parameter grid for sweep
PARAM_GRID = {
    # rdShapeAlign parameters
    'opt_param': [0.0, 0.25, 0.5, 0.75, 1.0],      # 0=color only, 1=shape only
    'max_preiters': [10, 30, 50],                   # Phase 1 iterations
    'max_postiters': [30, 50, 100],                 # Phase 2 iterations (NEW!)

    # Consensus parameters
    'tolerance': [1.5, 2.0, 2.5],                   # Angstroms
    'occurrence_threshold': [0.4, 0.5, 0.6],        # Feature frequency
}

# Reduced grid for quick testing
QUICK_GRID = {
    'opt_param': [0.0, 0.5, 1.0],
    'max_preiters': [10, 50],
    'max_postiters': [30, 100],
    'tolerance': [2.0],
    'occurrence_threshold': [0.5],
}

# Data paths
DATA_DIR = project_root / 'tutorials' / 'data'
REFERENCE_FILE = DATA_DIR / 'CCR2_reference_ligands.sdf'
ACTIVES_FILE = DATA_DIR / 'actives_ccr2_N75.csv'
DECOYS_FILE = DATA_DIR / 'decoys_ccr2_N500.csv'

# Output
RESULTS_DIR = project_root / 'experiments' / 'results'
RESULTS_FILE = RESULTS_DIR / 'rdshape_param_sweep_results.csv'

# Scoring settings
N_CONFORMERS = 10
RANDOM_SEED = 42


# =============================================================================
# Data Loading
# =============================================================================

def load_reference_mols(filepath: Path) -> list:
    """Load reference molecules from SDF file."""
    supplier = Chem.SDMolSupplier(str(filepath))
    mols = [mol for mol in supplier if mol is not None]
    print(f"Loaded {len(mols)} reference molecules")
    return mols


def load_mols_from_csv(filepath: Path, smiles_col: str = 'smiles') -> list:
    """Load molecules from CSV with SMILES column."""
    df = pd.read_csv(filepath)

    # Try different column names
    for col in [smiles_col, 'SMILES', 'Smiles', 'canonical_smiles']:
        if col in df.columns:
            smiles_list = df[col].tolist()
            break
    else:
        raise ValueError(f"No SMILES column found in {filepath}")

    mols = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append(mol)

    print(f"Loaded {len(mols)} molecules from {filepath.name}")
    return mols


def generate_conformers(mol: Chem.Mol, n_conformers: int = 10) -> list:
    """Generate conformers for a molecule."""
    mol_h = Chem.AddHs(mol)

    try:
        conf_ids = AllChem.EmbedMultipleConfs(
            mol_h,
            numConfs=n_conformers,
            randomSeed=RANDOM_SEED,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True
        )

        if not conf_ids:
            AllChem.EmbedMolecule(mol_h, randomSeed=RANDOM_SEED)

        if mol_h.GetNumConformers() > 0:
            return [mol_h]
    except Exception:
        pass

    return []


# =============================================================================
# Scoring Functions
# =============================================================================

def score_molecule(
    probe_mol: Chem.Mol,
    ref_mol: Chem.Mol,
    opt_param: float = 0.5,
    max_preiters: int = 50,
    max_postiters: int = 100,
    use_colors: bool = True,
    n_conformers: int = 10
) -> tuple:
    """
    Score a molecule against a reference using rdShapeAlign.

    Returns:
        Tuple of (best_shape, best_color, best_combo)
    """
    conformers = generate_conformers(probe_mol, n_conformers)

    if not conformers:
        return 0.0, 0.0, 0.0

    best_shape = 0.0
    best_color = 0.0
    best_combo = 0.0

    for conf_mol in conformers:
        for conf_id in range(conf_mol.GetNumConformers()):
            try:
                shape, color = AlignMol(
                    ref=ref_mol,
                    probe=conf_mol,
                    probeConfId=conf_id,
                    useColors=use_colors,
                    opt_param=opt_param,
                    max_preiters=max_preiters,
                    max_postiters=max_postiters
                )

                combo = shape + color

                if combo > best_combo:
                    best_combo = combo
                    best_shape = shape
                    best_color = color

            except Exception as e:
                continue

    return best_shape, best_color, best_combo


def evaluate_params(
    params: dict,
    reference_mols: list,
    actives: list,
    decoys: list,
    verbose: bool = False
) -> dict:
    """
    Evaluate a parameter combination.

    Returns dict with metrics and timing.
    """
    start_time = time.time()

    # Generate consensus pharmacophore
    consensus = PharmacophoreConsensus(
        tolerance=params['tolerance'],
        occurrence_threshold=params['occurrence_threshold']
    )

    try:
        features = consensus.generate_consensus(reference_mols)
    except Exception as e:
        return {**params, 'error': str(e), 'roc_auc': 0.5}

    if not features or len(features) < 2:
        return {**params, 'error': 'Too few features', 'roc_auc': 0.5, 'n_features': len(features) if features else 0}

    # Convert to mol
    try:
        pharm_mol = PharmacophoreToMol.convert(
            features,
            name='Consensus',
            enable_color_features=True
        )
    except Exception as e:
        return {**params, 'error': f'Mol conversion: {e}', 'roc_auc': 0.5}

    # Score all molecules
    active_scores = []
    for mol in actives:
        _, _, combo = score_molecule(
            mol, pharm_mol,
            opt_param=params['opt_param'],
            max_preiters=params['max_preiters'],
            max_postiters=params['max_postiters'],
            use_colors=True,
            n_conformers=N_CONFORMERS
        )
        active_scores.append(combo)

    decoy_scores = []
    for mol in decoys:
        _, _, combo = score_molecule(
            mol, pharm_mol,
            opt_param=params['opt_param'],
            max_preiters=params['max_preiters'],
            max_postiters=params['max_postiters'],
            use_colors=True,
            n_conformers=N_CONFORMERS
        )
        decoy_scores.append(combo)

    elapsed = time.time() - start_time

    # Calculate metrics
    y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
    y_scores = active_scores + decoy_scores

    try:
        auc = roc_auc_score(y_true, y_scores)
    except Exception:
        auc = 0.5

    try:
        bedroc_score = bedroc(y_true, y_scores, alpha=20.0)
    except Exception:
        bedroc_score = 0.0

    try:
        ef1 = enrichment_factor(y_true, y_scores, percentage=0.01)
        ef5 = enrichment_factor(y_true, y_scores, percentage=0.05)
    except Exception:
        ef1, ef5 = 0.0, 0.0

    # Score statistics
    active_mean = np.mean(active_scores)
    decoy_mean = np.mean(decoy_scores)
    separation = active_mean - decoy_mean

    return {
        **params,
        'roc_auc': auc,
        'bedroc': bedroc_score,
        'ef_1': ef1,
        'ef_5': ef5,
        'active_mean': active_mean,
        'decoy_mean': decoy_mean,
        'separation': separation,
        'n_features': len(features),
        'time_sec': elapsed,
        'time_per_mol_ms': elapsed / (len(actives) + len(decoys)) * 1000,
        'error': None
    }


# =============================================================================
# Main Sweep
# =============================================================================

def run_parameter_sweep(
    reference_mols: list,
    actives: list,
    decoys: list,
    param_grid: dict = None,
    quick: bool = False,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Run full parameter sweep.

    Args:
        reference_mols: Reference molecules for consensus
        actives: Active compounds
        decoys: Decoy compounds
        param_grid: Parameter grid (uses default if None)
        quick: Use reduced grid for quick testing
        verbose: Print progress

    Returns:
        DataFrame with all results
    """
    if param_grid is None:
        param_grid = QUICK_GRID if quick else PARAM_GRID

    # Generate all parameter combinations
    param_names = list(param_grid.keys())
    param_values = list(param_grid.values())
    combinations = list(itertools.product(*param_values))

    n_total = len(combinations)
    print(f"\n{'='*60}")
    print(f"rdShapeAlign Parameter Sweep")
    print(f"{'='*60}")
    print(f"Parameter combinations: {n_total}")
    print(f"Reference molecules: {len(reference_mols)}")
    print(f"Actives: {len(actives)}")
    print(f"Decoys: {len(decoys)}")
    print(f"Conformers per mol: {N_CONFORMERS}")
    print(f"{'='*60}\n")

    results = []
    best_auc = 0.0
    best_params = None

    for i, values in enumerate(combinations, 1):
        params = dict(zip(param_names, values))

        if verbose:
            print(f"[{i}/{n_total}] opt={params['opt_param']:.2f}, "
                  f"pre={params['max_preiters']}, post={params['max_postiters']}, "
                  f"tol={params['tolerance']}, occ={params['occurrence_threshold']}",
                  end=' ... ')

        result = evaluate_params(params, reference_mols, actives, decoys)
        results.append(result)

        if verbose:
            if result.get('error'):
                print(f"ERROR: {result['error']}")
            else:
                print(f"AUC={result['roc_auc']:.3f}, "
                      f"BEDROC={result['bedroc']:.3f}, "
                      f"EF@1%={result['ef_1']:.1f}, "
                      f"n_feat={result['n_features']}, "
                      f"time={result['time_sec']:.1f}s")

        if result['roc_auc'] > best_auc:
            best_auc = result['roc_auc']
            best_params = params.copy()

    # Create DataFrame
    df = pd.DataFrame(results)
    df = df.sort_values('roc_auc', ascending=False).reset_index(drop=True)

    # Print summary
    print(f"\n{'='*60}")
    print(f"RESULTS SUMMARY")
    print(f"{'='*60}")
    print(f"\nBest AUC: {best_auc:.4f}")
    print(f"Best parameters:")
    for k, v in best_params.items():
        print(f"  {k}: {v}")

    print(f"\nTop 5 configurations:")
    print(df[['opt_param', 'max_preiters', 'max_postiters', 'tolerance',
              'occurrence_threshold', 'roc_auc', 'bedroc', 'ef_1', 'n_features']].head())

    return df


def test_fragment_color_features():
    """
    Verify that our pharmacophore fragments are recognized by rdShapeAlign.

    This tests whether color scoring works with disconnected fragments.
    """
    print("\n" + "="*60)
    print("Testing Fragment Color Feature Recognition")
    print("="*60)

    # Test fragments
    fragments = {
        'Donor (NH3)': '[NH3]',
        'Acceptor (O)': '[O]',
        'Aromatic (benzene)': 'c1ccccc1',
        'Hydrophobe (CH4)': '[CH4]',
    }

    for name, smiles in fragments.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  {name}: FAILED to create mol")
            continue

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)

        if mol.GetNumConformers() == 0:
            print(f"  {name}: FAILED to embed")
            continue

        # Self-alignment should give shape=1.0, color depends on recognition
        try:
            shape, color = AlignMol(
                ref=mol,
                probe=mol,
                useColors=True,
                opt_param=0.5
            )
            status = "OK" if color > 0.5 else "LOW COLOR"
            print(f"  {name}: shape={shape:.3f}, color={color:.3f} [{status}]")
        except Exception as e:
            print(f"  {name}: ERROR - {e}")

    print()


# =============================================================================
# Main
# =============================================================================

def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description='rdShapeAlign Parameter Sweep')
    parser.add_argument('--quick', action='store_true',
                        help='Use reduced parameter grid for quick testing')
    parser.add_argument('--test-fragments', action='store_true',
                        help='Only test fragment color recognition')
    args = parser.parse_args()

    # Test fragment recognition
    if args.test_fragments:
        test_fragment_color_features()
        return

    # Always run fragment test first
    test_fragment_color_features()

    # Load data
    print("Loading data...")
    reference_mols = load_reference_mols(REFERENCE_FILE)
    actives = load_mols_from_csv(ACTIVES_FILE)
    decoys = load_mols_from_csv(DECOYS_FILE)

    # Run sweep
    results_df = run_parameter_sweep(
        reference_mols=reference_mols,
        actives=actives,
        decoys=decoys,
        quick=args.quick,
        verbose=True
    )

    # Save results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(RESULTS_FILE, index=False)
    print(f"\nResults saved to: {RESULTS_FILE}")

    # Summary statistics
    print("\n" + "="*60)
    print("PARAMETER IMPORTANCE ANALYSIS")
    print("="*60)

    for param in ['opt_param', 'max_preiters', 'max_postiters', 'tolerance', 'occurrence_threshold']:
        grouped = results_df.groupby(param)['roc_auc'].agg(['mean', 'std', 'max'])
        print(f"\n{param}:")
        print(grouped.to_string())


if __name__ == '__main__':
    main()

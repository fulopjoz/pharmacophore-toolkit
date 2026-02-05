"""Ablation study: measure the contribution of each improvement individually.

Tests 8 configurations on the CCR2 dataset:
  - RefEnsemble: raw shape similarity baseline (no consensus model)
  - Configs that test consensus/scoring improvements use fixed parameters
  - Configs that test BO improvements use ComboPharmacophoreOptimizer

Configurations:
  0. RefEnsemble  — direct reference ensemble scoring (no consensus, no clustering)
  1. Baseline     — no improvements (spatial off, EI, no retest, standard 3D, agglom)
  2. +spatial     — SPATIAL_SCALE per-type tolerance on
  3. +gp_hedge    — adaptive acquisition function (via BO)
  4. +retest      — post-hoc retest policy (via BO)
  5. +ensemble_clust — stability-aware ensemble clustering
  6. +ot_scoring  — cascading multi-level RRF scoring
  7. +all         — all improvements combined (via BO)

Usage:
    python experiments/ablation_study.py
    python experiments/ablation_study.py --n-conformers 10 --bo-budget 15
    python experiments/ablation_study.py --run-bo   # include BO configs (slow)
"""

import os
import sys
import time
import json
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, SDMolSupplier
from joblib import Memory, Parallel, delayed
from sklearn.metrics import roc_auc_score

from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.evaluation import UnifiedEvaluator, EvaluationConfig
from pharmacophore.screening_metrics import calculate_all_metrics
from pharmacophore.rdshape_optimizer import ReferenceEnsembleScorer

# Optional dependencies
try:
    from pharmacophore.ensemble_consensus import EnsembleConsensus
    HAS_ENSEMBLE = True
except ImportError:
    HAS_ENSEMBLE = False

try:
    from pharmacophore.combo_optimizer import ComboPharmacophoreOptimizer
    HAS_COMBO = True
except ImportError:
    HAS_COMBO = False

# Setup
N_JOBS = -1
cache_dir = Path(__file__).parent / '.cache'
cache_dir.mkdir(exist_ok=True)
memory = Memory(cache_dir, verbose=0)


@memory.cache
def generate_conformers_cached(smiles, n_conformers, random_seed=42, minimize=False):
    """Cached conformer generation (ETKDGv3, optional MMFF94s minimization)."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.numThreads = 0
    params.pruneRmsThresh = 0.5
    AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)
    if minimize and mol_h.GetNumConformers() > 0:
        AllChem.MMFFOptimizeMoleculeConfs(mol_h, numThreads=0, maxIters=200)
    if mol_h.GetNumConformers() > 0:
        return mol_h
    return None


def load_ccr2_data(n_conformers=5, minimize=False, verbose=True):
    """Load CCR2 dataset with conformers."""
    base = Path(__file__).parent.parent / 'tutorials' / 'data'

    # References from SDF
    refs = []
    ref_smiles = []
    supplier = SDMolSupplier(str(base / 'CCR2_reference_ligands.sdf'), removeHs=False)
    for mol in supplier:
        if mol and mol.GetNumConformers() > 0:
            refs.append(mol)
            ref_smiles.append(Chem.MolToSmiles(mol))

    # Actives
    actives_df = pd.read_csv(base / 'actives_ccr2_N75.csv')
    active_smiles = actives_df['SMILES'].tolist()

    # Decoys
    decoys_df = pd.read_csv(base / 'decoys_ccr2_N500.csv')
    decoy_smiles = decoys_df['Smiles'].tolist()

    if verbose:
        print(f"  References: {len(refs)}")
        print(f"  Actives: {len(active_smiles)}, generating {n_conformers} conformers...")

    actives = Parallel(n_jobs=N_JOBS, backend='loky')(
        delayed(generate_conformers_cached)(smi, n_conformers, minimize=minimize)
        for smi in active_smiles
    )
    actives = [m for m in actives if m is not None]

    if verbose:
        print(f"  Actives loaded: {len(actives)}")
        print(f"  Decoys: {len(decoy_smiles)}, generating {n_conformers} conformers...")

    decoys = Parallel(n_jobs=N_JOBS, backend='loky')(
        delayed(generate_conformers_cached)(smi, n_conformers, minimize=minimize)
        for smi in decoy_smiles
    )
    decoys = [m for m in decoys if m is not None]

    if verbose:
        print(f"  Decoys loaded: {len(decoys)}")

    return refs, ref_smiles, actives, decoys, active_smiles, decoy_smiles


# ---------------------------------------------------------------------------
# Evaluation helpers
# ---------------------------------------------------------------------------

def evaluate_ref_ensemble(refs, actives, decoys, *, n_conformers=5,
                          aggregation='max', opt_param=0.5, seed=42,
                          n_jobs=1, minimize=False):
    """Evaluate using direct reference ensemble scoring (no consensus model).

    Scores each query molecule against the reference molecules directly using
    ReferenceEnsembleScorer. This is the raw shape similarity baseline —
    no clustering, no pharmacophore construction.
    """
    scorer = ReferenceEnsembleScorer(
        reference_mols=refs, opt_param=opt_param,
        n_conformers=n_conformers, aggregation=aggregation,
        random_seed=seed, n_jobs=n_jobs, minimize=minimize,
    )
    active_scores = scorer.score_batch(actives, verbose=False)
    decoy_scores = scorer.score_batch(decoys, verbose=False)
    y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
    y_scores = active_scores + decoy_scores
    metrics = calculate_all_metrics(y_true, y_scores)
    return {
        'roc_auc': metrics['roc_auc'],
        'bedroc': metrics['bedroc'],
        'ef_1': metrics['ef_1'],
        'n_features': len(refs),
    }


def evaluate_fixed_params(refs, actives, decoys, *,
                          use_spatial_scaling=True,
                          use_ensemble_clustering=False,
                          use_ensemble_scoring=False,
                          minimize=False,
                          n_conformers=10, seed=42, n_jobs=1):
    """Evaluate with fixed parameters (no BO).

    Used for configs that test consensus/scoring improvements where we want
    to isolate the effect without confounding with BO search randomness.
    """
    evaluator = UnifiedEvaluator(refs, actives, decoys, seed, n_jobs=n_jobs)

    # --- Build consensus features ---
    if use_ensemble_clustering and HAS_ENSEMBLE:
        ec = EnsembleConsensus(
            n_runs=25,
            tolerance_range=(1.5, 2.5),
            occurrence_range=(0.3, 0.7),
            stability_threshold=0.3,
            random_state=seed,
        )
        features, stability = ec.generate_consensus_with_scores(refs)
    else:
        consensus = PharmacophoreConsensus(
            tolerance=2.0,
            occurrence_threshold=0.5,
            use_spatial_scaling=use_spatial_scaling,
        )
        features = consensus.generate_consensus(refs)

    if len(features) < 2:
        return {'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0, 'n_features': len(features)}

    # --- Score ---
    if use_ensemble_scoring:
        config = EvaluationConfig(
            tolerance=2.0, occurrence=0.5, shape_weight=0.5,
            opt_param=0.5, n_conformers=n_conformers,
            minimize=minimize,
        )
        result = evaluator.evaluate_ensemble(config, quick_mode=False)
    else:
        result = evaluator.evaluate_feature_subset(
            features, shape_weight=0.5, opt_param=0.5,
            n_conformers=n_conformers, minimize=minimize,
        )

    return {
        'roc_auc': result.roc_auc,
        'bedroc': result.bedroc,
        'ef_1': result.ef_1 if hasattr(result, 'ef_1') else 0.0,
        'n_features': result.n_features,
    }


def evaluate_with_bo(ref_smiles, active_smiles, decoy_smiles, *,
                     acq_func="EI", enable_retest=False,
                     use_spatial_scaling=True,
                     use_ensemble_scoring=False,
                     n_calls=15, seed=42):
    """Evaluate using ComboPharmacophoreOptimizer for BO ablation.

    Used for configs that test acq_func and retest improvements.
    """
    if not HAS_COMBO:
        return {'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0,
                'n_features': 0, 'note': 'scikit-optimize not installed'}

    optimizer = ComboPharmacophoreOptimizer(
        verbose=False, random_state=seed,
    )
    optimizer.load_from_smiles(
        reference_smiles=ref_smiles,
        active_smiles=active_smiles,
        decoy_smiles=decoy_smiles,
    )

    result = optimizer.optimize(
        n_calls=n_calls,
        n_random_starts=min(5, max(1, n_calls - 2)),
        acq_func=acq_func,
        enable_retest=enable_retest,
    )

    metrics = result.get('best_metrics', {})
    return {
        'roc_auc': result.get('best_auc', 0.5),
        'bedroc': metrics.get('bedroc', 0.0),
        'ef_1': metrics.get('ef_1', 0.0),
        'n_features': len(result.get('best_features', [])),
        'best_params': result.get('best_params', {}),
    }


# ---------------------------------------------------------------------------
# Ablation configurations
# ---------------------------------------------------------------------------

ABLATION_CONFIGS = [
    {
        'name': 'RefEnsemble (max)',
        'method': 'ref_ensemble',
    },
    {
        'name': 'Baseline',
        'method': 'fixed',
        'use_spatial_scaling': False,
        'use_ensemble_clustering': False,
        'use_ensemble_scoring': False,
    },
    {
        'name': '+spatial',
        'method': 'fixed',
        'use_spatial_scaling': True,
        'use_ensemble_clustering': False,
        'use_ensemble_scoring': False,
    },
    {
        'name': '+gp_hedge',
        'method': 'bo',
        'acq_func': 'gp_hedge',
        'enable_retest': False,
    },
    {
        'name': '+retest',
        'method': 'bo',
        'acq_func': 'EI',
        'enable_retest': True,
    },
    {
        'name': '+ensemble_clust',
        'method': 'fixed',
        'use_spatial_scaling': False,
        'use_ensemble_clustering': True,
        'use_ensemble_scoring': False,
    },
    {
        'name': '+ot_scoring',
        'method': 'fixed',
        'use_spatial_scaling': False,
        'use_ensemble_clustering': False,
        'use_ensemble_scoring': True,
    },
    {
        'name': '+minimize',
        'method': 'fixed',
        'use_spatial_scaling': False,
        'use_ensemble_clustering': False,
        'use_ensemble_scoring': False,
        'minimize': True,
    },
    {
        'name': '+all (BO)',
        'method': 'bo',
        'acq_func': 'gp_hedge',
        'enable_retest': True,
        'use_spatial_scaling': True,
        'use_ensemble_scoring': True,
    },
]


def run_single_config(config, refs, ref_smiles, actives, decoys,
                      active_smiles, decoy_smiles, n_conformers, bo_budget, seed):
    """Run a single ablation configuration and return results dict."""
    name = config['name']
    start = time.time()

    if config['method'] == 'ref_ensemble':
        metrics = evaluate_ref_ensemble(
            refs, actives, decoys,
            n_conformers=n_conformers, seed=seed,
            n_jobs=-1,
        )
    elif config['method'] == 'fixed':
        metrics = evaluate_fixed_params(
            refs, actives, decoys,
            use_spatial_scaling=config.get('use_spatial_scaling', False),
            use_ensemble_clustering=config.get('use_ensemble_clustering', False),
            use_ensemble_scoring=config.get('use_ensemble_scoring', False),
            minimize=config.get('minimize', False),
            n_conformers=n_conformers,
            seed=seed, n_jobs=-1,
        )
    else:  # 'bo'
        metrics = evaluate_with_bo(
            ref_smiles, active_smiles, decoy_smiles,
            acq_func=config.get('acq_func', 'EI'),
            enable_retest=config.get('enable_retest', False),
            use_spatial_scaling=config.get('use_spatial_scaling', False),
            use_ensemble_scoring=config.get('use_ensemble_scoring', False),
            n_calls=bo_budget,
            seed=seed,
        )

    elapsed = time.time() - start

    return {
        'config': name,
        'method': config['method'],
        'roc_auc': metrics.get('roc_auc', 0.5),
        'bedroc': metrics.get('bedroc', 0.0),
        'ef_1': metrics.get('ef_1', 0.0),
        'n_features': metrics.get('n_features', 0),
        'wall_time_sec': elapsed,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Ablation study for pharmacophore improvements")
    parser.add_argument('--n-conformers', type=int, default=5)
    parser.add_argument('--bo-budget', type=int, default=15,
                        help='BO evaluation budget for optimizer configs')
    parser.add_argument('--run-bo', action='store_true', default=False,
                        help='Include BO configs (slow; skipped by default)')
    parser.add_argument('--minimize', action='store_true', default=False,
                        help='Enable MMFF94s energy minimization after ETKDGv3')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print(f"\n{'='*75}")
    print(f"ABLATION STUDY: Individual Improvement Contributions")
    print(f"{'='*75}")
    print(f"  Dataset: CCR2 (74 actives, 500 decoys)")
    print(f"  Conformers: {args.n_conformers}")
    print(f"  BO budget: {args.bo_budget} calls (for optimizer configs)")
    print(f"  Run BO configs: {args.run_bo}")
    print(f"  Minimize (MMFF94s): {args.minimize}")
    print(f"  Seed: {args.seed}")
    print(f"  CPU cores: {os.cpu_count()}")
    print(f"{'='*75}\n")

    # Load data
    print("Loading CCR2 dataset...")
    refs, ref_smiles, actives, decoys, active_smi, decoy_smi = load_ccr2_data(
        n_conformers=args.n_conformers, minimize=args.minimize, verbose=True,
    )
    print()

    # Run all configs sequentially (each prints progress)
    results = []
    for i, config in enumerate(ABLATION_CONFIGS, 1):
        name = config['name']
        method = config['method']
        method_label = {'ref_ensemble': 'ref', 'fixed': 'fixed', 'bo': 'BO'}[method]

        # Skip BO configs unless --run-bo
        if method == 'bo' and not args.run_bo:
            print(f"[{i}/{len(ABLATION_CONFIGS)}] Skipping {name} (use --run-bo)")
            results.append({
                'config': name, 'method': method,
                'roc_auc': float('nan'), 'bedroc': float('nan'),
                'ef_1': float('nan'), 'n_features': 0,
                'wall_time_sec': 0.0, 'skipped': True,
            })
            continue

        print(f"[{i}/{len(ABLATION_CONFIGS)}] Running {name} ({method_label})...")

        try:
            r = run_single_config(
                config, refs, ref_smiles, actives, decoys,
                active_smi, decoy_smi,
                n_conformers=args.n_conformers,
                bo_budget=args.bo_budget,
                seed=args.seed,
            )
            results.append(r)
            print(f"         AUC={r['roc_auc']:.4f}  BEDROC={r['bedroc']:.4f}  "
                  f"EF@1%={r['ef_1']:.2f}  features={r['n_features']}  "
                  f"time={r['wall_time_sec']:.1f}s")
        except Exception as e:
            print(f"         FAILED: {e}")
            results.append({
                'config': name, 'method': config['method'],
                'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0,
                'n_features': 0, 'wall_time_sec': 0.0,
                'error': str(e),
            })

    # === Print results table ===
    print(f"\n{'='*75}")
    print(f"ABLATION RESULTS")
    print(f"{'='*75}\n")

    # Delta column uses Baseline row, not first row
    baseline_auc = next(
        (r['roc_auc'] for r in results if r['config'] == 'Baseline'), 0.5
    )

    print(f"{'Config':<20} {'Method':>6} {'AUC':>7} {'delta':>7} "
          f"{'BEDROC':>8} {'EF@1%':>7} {'Feat':>5} {'Time(s)':>8}")
    print(f"{'-'*75}")

    for r in results:
        is_skipped = r.get('skipped', False)
        is_ref = r['method'] == 'ref_ensemble'
        is_baseline = r['config'] == 'Baseline'

        if is_skipped:
            delta_str = " skip"
        elif is_ref:
            delta_str = "  ref"
        elif is_baseline:
            delta_str = "  ---"
        else:
            delta = r['roc_auc'] - baseline_auc
            delta_str = f"{delta:+.4f}"

        method_label = {'ref_ensemble': 'ref', 'fixed': 'fixed', 'bo': 'BO'}[r['method']]

        if is_skipped:
            print(f"{r['config']:<20} {method_label:>6} {'—':>7} {delta_str:>7} "
                  f"{'—':>8} {'—':>7} {r['n_features']:>5} "
                  f"{r['wall_time_sec']:>8.1f}")
        else:
            print(f"{r['config']:<20} {method_label:>6} {r['roc_auc']:>7.4f} {delta_str:>7} "
                  f"{r['bedroc']:>8.4f} {r['ef_1']:>7.2f} {r['n_features']:>5} "
                  f"{r['wall_time_sec']:>8.1f}")

    print(f"{'-'*75}")

    # Footer for skipped configs
    n_skipped = sum(1 for r in results if r.get('skipped', False))
    if n_skipped > 0:
        print(f"\n  {n_skipped} BO config(s) skipped. Run with --run-bo to include.")

    # Best improvement (excluding ref_ensemble and skipped)
    improvements = [
        (r['config'], r['roc_auc'] - baseline_auc)
        for r in results
        if r['config'] != 'Baseline'
        and r['method'] != 'ref_ensemble'
        and not r.get('skipped', False)
        and r['roc_auc'] > baseline_auc
    ]
    if improvements:
        best_name, best_delta = max(improvements, key=lambda x: x[1])
        print(f"\n  Largest single improvement: {best_name} (+{best_delta:.4f} AUC)")

    # Show +all result if not skipped
    all_result = next((r for r in results if r['config'] == '+all (BO)'), None)
    if all_result and not all_result.get('skipped', False):
        print(f"  Combined (+all):           AUC={all_result['roc_auc']:.4f} "
              f"(delta={all_result['roc_auc'] - baseline_auc:+.4f})")

    # Save results
    output_dir = Path(__file__).parent / 'results'
    output_dir.mkdir(exist_ok=True)

    json_file = output_dir / f'ablation_{timestamp}.json'
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  JSON saved: {json_file}")

    csv_file = output_dir / f'ablation_{timestamp}.csv'
    pd.DataFrame(results).to_csv(csv_file, index=False)
    print(f"  CSV saved:  {csv_file}")


if __name__ == '__main__':
    main()

"""Full-dataset comparison of ALL pharmacophore scoring approaches.

Compares 7 approaches on the complete CCR2 dataset (74 actives, 500 decoys):
  1. Optuna GP-Sampler     (3D rdShapeAlign, Bayesian optimization)
  2. Optuna NSGA-II        (3D rdShapeAlign, evolutionary multi-objective)
  3. HypoGen 3-Phase       (3D rdShapeAlign, constructive-subtractive)
  4. Pharm2D Fingerprint    (2D pharmacophore fingerprints)
  5. Ensemble Consensus     (stability-aware clustering + 3D scoring)
  6. Strategy Selector      (tournament across clustering strategies)
  7. Ensemble Scoring       (cascading multi-level scoring with RRF)

Usage:
    python experiments/compare_all_approaches.py
    python experiments/compare_all_approaches.py --n-conformers 10 --n-trials 15
"""

import os
import sys
import time
import json
import itertools
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, SDMolSupplier
from joblib import Memory, Parallel, delayed
from sklearn.metrics import roc_auc_score

from pharmacophore.evaluation import UnifiedEvaluator, EvaluationConfig
from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.pharm2d_scoring import Pharm2DScorer
from pharmacophore.screening_metrics import calculate_all_metrics

try:
    import optuna
    from optuna.samplers import NSGAIISampler
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    HAS_OPTUNA = True
except ImportError:
    HAS_OPTUNA = False

# Setup
N_JOBS = -1
cache_dir = Path(__file__).parent / '.cache'
cache_dir.mkdir(exist_ok=True)
memory = Memory(cache_dir, verbose=0)


@memory.cache
def generate_conformers_cached(smiles, n_conformers, random_seed=42):
    """Cached conformer generation (ETKDGv3 only, no UFF)."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)
    if mol_h.GetNumConformers() > 0:
        return mol_h
    return None


def load_full_ccr2(n_conformers=5, verbose=True):
    """Load complete CCR2 dataset with parallel conformer generation."""
    base = Path(__file__).parent.parent / 'tutorials' / 'data'

    # References
    refs = []
    supplier = SDMolSupplier(str(base / 'CCR2_reference_ligands.sdf'), removeHs=False)
    for mol in supplier:
        if mol and mol.GetNumConformers() > 0:
            refs.append(mol)

    # Actives (ALL)
    actives_df = pd.read_csv(base / 'actives_ccr2_N75.csv')
    active_smiles = actives_df['SMILES'].tolist()

    # Decoys (ALL)
    decoys_df = pd.read_csv(base / 'decoys_ccr2_N500.csv')
    decoy_smiles = decoys_df['Smiles'].tolist()

    if verbose:
        print(f"  References: {len(refs)}")
        print(f"  Actives: {len(active_smiles)} SMILES, generating {n_conformers} conformers each...")

    actives = Parallel(n_jobs=N_JOBS, backend='loky')(
        delayed(generate_conformers_cached)(smi, n_conformers) for smi in active_smiles
    )
    actives = [m for m in actives if m is not None]

    if verbose:
        print(f"  Actives loaded: {len(actives)} (of {len(active_smiles)})")
        print(f"  Decoys: {len(decoy_smiles)} SMILES, generating {n_conformers} conformers each...")

    decoys = Parallel(n_jobs=N_JOBS, backend='loky')(
        delayed(generate_conformers_cached)(smi, n_conformers) for smi in decoy_smiles
    )
    decoys = [m for m in decoys if m is not None]

    if verbose:
        print(f"  Decoys loaded: {len(decoys)} (of {len(decoy_smiles)})")
        print(f"  Total: {len(actives)} actives + {len(decoys)} decoys = {len(actives) + len(decoys)} molecules\n")

    return refs, actives, decoys, active_smiles, decoy_smiles


def run_optuna_approach(refs, actives, decoys, sampler_name, n_trials=10, seed=42):
    """Run Optuna-based optimization (GP or NSGA-II)."""
    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)

    if sampler_name == 'gp':
        try:
            from optuna.samplers import GPSampler
            sampler = GPSampler(seed=seed)
            label = 'Optuna GP'
        except Exception:
            sampler = optuna.samplers.TPESampler(seed=seed)
            label = 'Optuna TPE'
    else:
        sampler = NSGAIISampler(seed=seed)
        label = 'Optuna NSGA-II'

    study = optuna.create_study(directions=['maximize', 'maximize'], sampler=sampler)

    def objective(trial):
        config = EvaluationConfig(
            tolerance=trial.suggest_float('tolerance', 0.5, 4.0),
            occurrence=trial.suggest_float('occurrence', 0.1, 1.0),
            shape_weight=trial.suggest_float('shape_weight', 0.3, 0.9),
            opt_param=trial.suggest_float('opt_param', 0.0, 1.0),
            linkage=trial.suggest_categorical('linkage', ['average', 'complete', 'single', 'ward']),
            n_conformers=trial.suggest_int('n_conformers', 5, 20),
        )
        result = evaluator.evaluate(config)
        return result.roc_auc, result.bedroc

    start = time.time()
    study.optimize(objective, n_trials=n_trials, n_jobs=N_JOBS, show_progress_bar=True)
    elapsed = time.time() - start

    completed = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    if not completed:
        return {'approach': label, 'best_auc': 0.5, 'best_bedroc': 0.0,
                'wall_time_sec': elapsed, 'n_evals': 0, 'best_params': {}}

    best = max(completed, key=lambda t: t.values[0])
    return {
        'approach': label,
        'best_auc': best.values[0],
        'best_bedroc': best.values[1],
        'wall_time_sec': elapsed,
        'n_evals': len(completed),
        'best_params': best.params,
    }


def run_hypogen(refs, actives, decoys, seed=42):
    """Run HypoGen 3-phase with parallel Phase 1."""
    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)

    start = time.time()

    # Phase 1: Generate consensus and enumerate subsets
    consensus = PharmacophoreConsensus(tolerance=2.0, occurrence_threshold=0.3)
    features = consensus.generate_consensus(refs)

    subsets = []
    for n_feat in range(3, min(7, len(features) + 1)):
        combos = list(itertools.combinations(features, n_feat))
        subsets.extend([list(c) for c in combos[:20]])
        if len(subsets) >= 50:
            break

    print(f"    Phase 1: Evaluating {len(subsets)} hypotheses...")

    def eval_subset(subset):
        result = evaluator.evaluate_feature_subset(subset, 0.5, 0.5, 5)
        cost = 0.5 * (1 - result.roc_auc) + 0.3 * (1 - result.bedroc) + 0.2 * (len(subset) / 8)
        return {'features': subset, 'roc_auc': result.roc_auc, 'bedroc': result.bedroc, 'cost': cost}

    hypotheses = Parallel(n_jobs=N_JOBS, backend='loky')(
        delayed(eval_subset)(s) for s in subsets
    )

    # Phase 2: Filter
    survivors = [h for h in hypotheses if h['roc_auc'] >= 0.50]
    print(f"    Phase 2: {len(survivors)}/{len(hypotheses)} hypotheses survived")

    # Phase 3: Select best
    best = min(survivors, key=lambda h: h['cost']) if survivors else min(hypotheses, key=lambda h: h['cost'])
    elapsed = time.time() - start

    return {
        'approach': 'HypoGen 3-Phase',
        'best_auc': best['roc_auc'],
        'best_bedroc': best['bedroc'],
        'wall_time_sec': elapsed,
        'n_evals': len(hypotheses),
        'best_params': {'n_features': len(best['features']), 'cost': best['cost']},
    }


def run_pharm2d(refs, actives, decoys):
    """Run 2D pharmacophore fingerprint scoring."""
    start = time.time()

    scorer = Pharm2DScorer(refs)

    active_scores = scorer.score_all(actives)
    decoy_scores = scorer.score_all(decoys)

    y_true = np.concatenate([np.ones(len(active_scores)), np.zeros(len(decoy_scores))])
    y_scores = np.concatenate([active_scores, decoy_scores])

    # Filter out NaN
    valid = ~np.isnan(y_scores)
    y_true_valid = y_true[valid]
    y_scores_valid = y_scores[valid]

    auc = roc_auc_score(y_true_valid, y_scores_valid)

    # Calculate BEDROC using our metrics
    metrics = calculate_all_metrics(y_true_valid, y_scores_valid)

    elapsed = time.time() - start

    return {
        'approach': 'Pharm2D (2D)',
        'best_auc': auc,
        'best_bedroc': metrics.get('bedroc', 0.0),
        'wall_time_sec': elapsed,
        'n_evals': 1,
        'best_params': {'method': 'Gobbi_Pharm2D', 'similarity': 'Tanimoto'},
        'ef_1': metrics.get('ef_1', 0.0),
        'ef_5': metrics.get('ef_5', 0.0),
    }


def run_ensemble_consensus(refs, actives, decoys, seed=42):
    """Run ensemble consensus clustering + 3D shape scoring."""
    from pharmacophore.ensemble_consensus import EnsembleConsensus

    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)

    start = time.time()

    ec = EnsembleConsensus(
        n_runs=25,
        tolerance_range=(1.5, 2.5),
        occurrence_range=(0.3, 0.7),
        stability_threshold=0.3,
        random_state=seed,
    )
    features, stability = ec.generate_consensus_with_scores(refs)

    if len(features) < 2:
        elapsed = time.time() - start
        return {
            'approach': 'Ensemble Consensus',
            'best_auc': 0.5,
            'best_bedroc': 0.0,
            'wall_time_sec': elapsed,
            'n_evals': 1,
            'best_params': {'n_features': len(features)},
        }

    result = evaluator.evaluate_feature_subset(features, 0.5, 0.5, 10)
    elapsed = time.time() - start

    return {
        'approach': 'Ensemble Consensus',
        'best_auc': result.roc_auc,
        'best_bedroc': result.bedroc,
        'wall_time_sec': elapsed,
        'n_evals': 1,
        'best_params': {
            'n_features': len(features),
            'n_stable': sum(1 for s in stability if s >= 0.5),
            'avg_stability': float(np.mean(stability)) if stability else 0.0,
        },
    }


def run_strategy_selector(refs, actives, decoys, seed=42):
    """Run StrategySelector tournament across clustering strategies."""
    from pharmacophore.auto_strategy import StrategySelector

    start = time.time()

    selector = StrategySelector(
        refs, actives, decoys, random_state=seed, verbose=True
    )
    best = selector.select_best()

    elapsed = time.time() - start

    return {
        'approach': f'StrategySelector ({best.strategy_name})',
        'best_auc': best.eval_result.roc_auc,
        'best_bedroc': best.eval_result.bedroc,
        'wall_time_sec': elapsed,
        'n_evals': best.eval_result.n_features,
        'best_params': {
            'strategy': best.strategy_name,
            'n_features': best.eval_result.n_features,
            'sdbw': float(best.sdbw) if best.sdbw != float('inf') else None,
        },
    }


def run_ensemble_scoring(refs, actives, decoys, seed=42):
    """Run multi-level ensemble scoring (cascading 3D + OT with RRF)."""
    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)

    start = time.time()

    config = EvaluationConfig(
        tolerance=2.0, occurrence=0.5, shape_weight=0.5,
        opt_param=0.5, n_conformers=10,
    )
    result = evaluator.evaluate_ensemble(config, quick_mode=False)

    elapsed = time.time() - start

    return {
        'approach': 'Ensemble Scoring (RRF)',
        'best_auc': result.roc_auc,
        'best_bedroc': result.bedroc,
        'wall_time_sec': elapsed,
        'n_evals': 1,
        'best_params': {
            'scoring': 'reciprocal_rank_fusion',
            'n_features': result.n_features,
        },
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Full-dataset comparison of all approaches")
    parser.add_argument('--n-conformers', type=int, default=5, help='Conformers per molecule')
    parser.add_argument('--n-trials', type=int, default=10, help='Optuna trials per approach')
    args = parser.parse_args()

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print(f"\n{'='*70}")
    print(f"FULL-DATASET COMPARISON: All Pharmacophore Approaches")
    print(f"{'='*70}")
    print(f"  Dataset: Complete CCR2 (74 actives, 500 decoys)")
    print(f"  Conformers: {args.n_conformers} per molecule (ETKDGv3, no UFF)")
    print(f"  Optuna trials: {args.n_trials} per approach")
    print(f"  CPU cores: {os.cpu_count()}")
    print(f"  Timestamp: {timestamp}")
    print(f"{'='*70}\n")

    total_start = time.time()

    # Load data
    print("Loading CCR2 dataset...")
    refs, actives, decoys, active_smi, decoy_smi = load_full_ccr2(
        n_conformers=args.n_conformers, verbose=True
    )

    results = {}

    # --- Approach 1: Optuna GP ---
    print(f"{'='*70}")
    print(f"APPROACH 1: Optuna GP-Sampler ({args.n_trials} trials, 3D rdShapeAlign)")
    print(f"{'='*70}")
    results['gp'] = run_optuna_approach(refs, actives, decoys, 'gp', args.n_trials)
    print(f"  -> AUC: {results['gp']['best_auc']:.4f}, "
          f"BEDROC: {results['gp']['best_bedroc']:.4f}, "
          f"Time: {results['gp']['wall_time_sec']:.1f}s\n")

    # --- Approach 2: Optuna NSGA-II ---
    print(f"{'='*70}")
    print(f"APPROACH 2: Optuna NSGA-II ({args.n_trials} trials, 3D rdShapeAlign)")
    print(f"{'='*70}")
    results['nsga2'] = run_optuna_approach(refs, actives, decoys, 'nsga2', args.n_trials)
    print(f"  -> AUC: {results['nsga2']['best_auc']:.4f}, "
          f"BEDROC: {results['nsga2']['best_bedroc']:.4f}, "
          f"Time: {results['nsga2']['wall_time_sec']:.1f}s\n")

    # --- Approach 3: HypoGen ---
    print(f"{'='*70}")
    print(f"APPROACH 3: HypoGen 3-Phase (3D rdShapeAlign)")
    print(f"{'='*70}")
    results['hypogen'] = run_hypogen(refs, actives, decoys)
    print(f"  -> AUC: {results['hypogen']['best_auc']:.4f}, "
          f"BEDROC: {results['hypogen']['best_bedroc']:.4f}, "
          f"Time: {results['hypogen']['wall_time_sec']:.1f}s\n")

    # --- Approach 4: Pharm2D (2D fingerprints) ---
    print(f"{'='*70}")
    print(f"APPROACH 4: Pharm2D Fingerprints (2D, no alignment)")
    print(f"{'='*70}")
    results['pharm2d'] = run_pharm2d(refs, actives, decoys)
    print(f"  -> AUC: {results['pharm2d']['best_auc']:.4f}, "
          f"BEDROC: {results['pharm2d']['best_bedroc']:.4f}, "
          f"Time: {results['pharm2d']['wall_time_sec']:.1f}s\n")

    # --- Approach 5: Ensemble Consensus ---
    print(f"{'='*70}")
    print(f"APPROACH 5: Ensemble Consensus (stability-aware clustering + 3D)")
    print(f"{'='*70}")
    results['ensemble_consensus'] = run_ensemble_consensus(refs, actives, decoys)
    print(f"  -> AUC: {results['ensemble_consensus']['best_auc']:.4f}, "
          f"BEDROC: {results['ensemble_consensus']['best_bedroc']:.4f}, "
          f"Time: {results['ensemble_consensus']['wall_time_sec']:.1f}s\n")

    # --- Approach 6: Strategy Selector ---
    print(f"{'='*70}")
    print(f"APPROACH 6: Strategy Selector (tournament across strategies)")
    print(f"{'='*70}")
    results['strategy_selector'] = run_strategy_selector(refs, actives, decoys)
    print(f"  -> AUC: {results['strategy_selector']['best_auc']:.4f}, "
          f"BEDROC: {results['strategy_selector']['best_bedroc']:.4f}, "
          f"Time: {results['strategy_selector']['wall_time_sec']:.1f}s\n")

    # --- Approach 7: Ensemble Scoring (RRF) ---
    print(f"{'='*70}")
    print(f"APPROACH 7: Ensemble Scoring (cascading 3D + OT with RRF)")
    print(f"{'='*70}")
    results['ensemble_scoring'] = run_ensemble_scoring(refs, actives, decoys)
    print(f"  -> AUC: {results['ensemble_scoring']['best_auc']:.4f}, "
          f"BEDROC: {results['ensemble_scoring']['best_bedroc']:.4f}, "
          f"Time: {results['ensemble_scoring']['wall_time_sec']:.1f}s\n")

    total_elapsed = time.time() - total_start

    # === Print results ===
    print(f"\n{'='*70}")
    print(f"COMPARISON RESULTS — Complete CCR2 Dataset")
    print(f"{'='*70}\n")

    all_keys = ['gp', 'nsga2', 'hypogen', 'pharm2d',
                'ensemble_consensus', 'strategy_selector', 'ensemble_scoring']

    print(f"{'Approach':<30} {'Type':>5} {'AUC':>7} {'BEDROC':>8} {'Evals':>6} {'Time(s)':>8}")
    print(f"{'-'*75}")
    for key in all_keys:
        r = results[key]
        atype = '2D' if '2D' in r['approach'] else '3D'
        print(f"{r['approach']:<30} {atype:>5} {r['best_auc']:>7.4f} {r['best_bedroc']:>8.4f} "
              f"{r['n_evals']:>6} {r['wall_time_sec']:>8.1f}")
    print(f"{'-'*75}")

    best_auc_key = max(results, key=lambda k: results[k]['best_auc'])
    best_bedroc_key = max(results, key=lambda k: results[k]['best_bedroc'])
    fastest_key = min(results, key=lambda k: results[k]['wall_time_sec'])

    print(f"\n  Best AUC:    {results[best_auc_key]['approach']} ({results[best_auc_key]['best_auc']:.4f})")
    print(f"  Best BEDROC: {results[best_bedroc_key]['approach']} ({results[best_bedroc_key]['best_bedroc']:.4f})")
    print(f"  Fastest:     {results[fastest_key]['approach']} ({results[fastest_key]['wall_time_sec']:.1f}s)")
    print(f"  Total time:  {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")
    print(f"{'='*70}\n")

    # Save results
    output_dir = Path(__file__).parent / 'results'
    output_dir.mkdir(exist_ok=True)

    results_clean = {}
    for k, v in results.items():
        results_clean[k] = {key: val for key, val in v.items()
                            if key != 'study'}

    output_file = output_dir / f'comparison_{timestamp}.json'
    with open(output_file, 'w') as f:
        json.dump(results_clean, f, indent=2, default=str)
    print(f"Results saved: {output_file}")

    # Also save CSV summary
    csv_file = output_dir / f'comparison_{timestamp}.csv'
    rows = []
    for k, r in results.items():
        rows.append({
            'approach': r['approach'],
            'type': '2D' if '2D' in r['approach'] else '3D',
            'roc_auc': r['best_auc'],
            'bedroc': r['best_bedroc'],
            'n_evals': r['n_evals'],
            'wall_time_sec': r['wall_time_sec'],
        })
    pd.DataFrame(rows).to_csv(csv_file, index=False)
    print(f"CSV saved: {csv_file}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python
"""Phase 2 validation experiment: LORO CV, greedy selection, comparison.

Tests:
1. Leave-One-Reference-Out CV: For each of 5 refs, train on 4, measure
   AUC/BEDROC across baseline (max), Z-score, and learned weighting.
2. Greedy feature selection: Run on full CCR2, report selected features
   and AUC improvement curve.
3. Comparison table: All methods side-by-side.

Expected runtime: ~30 min (5 LORO folds x 3 methods + greedy selection).

Usage:
    python experiments/exp_validation_phases.py
    python experiments/exp_validation_phases.py --quick  # Skip greedy
"""

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pharmacophore.evaluation import UnifiedEvaluator, EvaluationConfig
from pharmacophore.learned_scoring import LearnedScorer
from pharmacophore.greedy_selector import GreedyFeatureSelector
from pharmacophore.screening_metrics import calculate_all_metrics

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
logger = logging.getLogger(__name__)

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'tutorials', 'data')


def load_ccr2_dataset():
    """Load CCR2 benchmark dataset."""
    refs_path = os.path.join(DATA_DIR, 'CCR2_reference_ligands.sdf')
    actives_path = os.path.join(DATA_DIR, 'actives_ccr2_N75.csv')
    decoys_path = os.path.join(DATA_DIR, 'decoys_ccr2_N500.csv')

    # Load references
    supplier = Chem.SDMolSupplier(refs_path, removeHs=False)
    ref_mols = [mol for mol in supplier if mol is not None]
    logger.info("Loaded %d reference molecules", len(ref_mols))

    # Load actives
    actives_df = pd.read_csv(actives_path)
    smiles_col = 'SMILES' if 'SMILES' in actives_df.columns else 'Smiles'
    actives = []
    for smi in actives_df[smiles_col]:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            actives.append(mol)
    logger.info("Loaded %d actives", len(actives))

    # Load decoys
    decoys_df = pd.read_csv(decoys_path)
    smiles_col = 'Smiles' if 'Smiles' in decoys_df.columns else 'SMILES'
    decoys = []
    for smi in decoys_df[smiles_col]:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            decoys.append(mol)
    logger.info("Loaded %d decoys", len(decoys))

    return ref_mols, actives, decoys


def run_loro_cv(ref_mols, actives, decoys, config):
    """Leave-One-Reference-Out cross-validation.

    For each reference, hold it out, train/evaluate with remaining refs.
    """
    n_refs = len(ref_mols)
    results = {
        'baseline_max': [],
        'zscore': [],
        'learned': [],
    }

    for hold_out_idx in range(n_refs):
        logger.info("LORO fold %d/%d (holding out ref %d)",
                     hold_out_idx + 1, n_refs, hold_out_idx)

        # Train refs = all except hold-out
        train_refs = [m for i, m in enumerate(ref_mols) if i != hold_out_idx]

        # Create evaluator with training refs
        evaluator = UnifiedEvaluator(train_refs, actives, decoys, random_state=42)

        # Baseline: max aggregation
        config_max = EvaluationConfig(
            n_conformers=config.n_conformers,
            opt_param=config.opt_param,
            aggregation='max',
            max_preiters=config.max_preiters,
            max_postiters=config.max_postiters,
        )
        result_max = evaluator.evaluate(config_max)
        results['baseline_max'].append({
            'fold': hold_out_idx,
            'roc_auc': result_max.roc_auc,
            'bedroc': result_max.bedroc,
            'ef_1': result_max.ef_1,
        })
        logger.info("  Baseline (max): AUC=%.4f BEDROC=%.4f",
                     result_max.roc_auc, result_max.bedroc)

        # Z-score: compute params then evaluate
        config_zscore = EvaluationConfig(
            n_conformers=config.n_conformers,
            opt_param=config.opt_param,
            aggregation='zscore',
            max_preiters=config.max_preiters,
            max_postiters=config.max_postiters,
        )
        evaluator.compute_zscore_params(config_zscore)
        result_zscore = evaluator.evaluate(config_zscore)
        results['zscore'].append({
            'fold': hold_out_idx,
            'roc_auc': result_zscore.roc_auc,
            'bedroc': result_zscore.bedroc,
            'ef_1': result_zscore.ef_1,
        })
        logger.info("  Z-score: AUC=%.4f BEDROC=%.4f",
                     result_zscore.roc_auc, result_zscore.bedroc)

        # Learned weighting
        learned = LearnedScorer(evaluator, random_state=42)
        fit_result = learned.fit(config_max, cv_folds=3)
        learned_metrics = learned.evaluate(config_max)
        results['learned'].append({
            'fold': hold_out_idx,
            'roc_auc': learned_metrics['roc_auc'],
            'bedroc': learned_metrics['bedroc'],
            'ef_1': learned_metrics['ef_1'],
            'cv_auc': fit_result['cv_auc'],
        })
        logger.info("  Learned: AUC=%.4f BEDROC=%.4f (CV AUC=%.4f)",
                     learned_metrics['roc_auc'], learned_metrics['bedroc'],
                     fit_result['cv_auc'])

    return results


def run_greedy_selection(ref_mols, actives, decoys):
    """Run greedy feature selection on full dataset."""
    logger.info("Running greedy feature selection...")
    selector = GreedyFeatureSelector(ref_mols, actives, decoys, random_state=42)
    result = selector.select(
        tolerance=2.0, occurrence=0.3,
        max_features=10, convergence_threshold=0.001,
        verbose=True,
    )

    logger.info("Selected %d features in %.1fs (%d evaluations)",
                len(result.selected_features), result.wall_time_sec,
                result.n_evaluations)
    logger.info("Final AUC=%.4f BEDROC=%.4f",
                result.best_auc, result.best_bedroc)

    return {
        'selected_indices': result.selected_indices,
        'selected_features': [
            [f[0], list(f[1]), f[2], f[3], f[4]]
            for f in result.selected_features
        ],
        'selection_history': [
            {'step': h[0], 'feature_idx': h[1], 'auc': h[2], 'bedroc': h[3]}
            for h in result.selection_history
        ],
        'best_auc': result.best_auc,
        'best_bedroc': result.best_bedroc,
        'n_evaluations': result.n_evaluations,
        'wall_time_sec': result.wall_time_sec,
    }


def print_comparison_table(loro_results):
    """Print comparison table of LORO results."""
    print("\n" + "=" * 70)
    print("LEAVE-ONE-REFERENCE-OUT CROSS-VALIDATION RESULTS")
    print("=" * 70)

    for method_name, folds in loro_results.items():
        aucs = [f['roc_auc'] for f in folds]
        bedrocs = [f['bedroc'] for f in folds]
        ef1s = [f['ef_1'] for f in folds]
        print(f"\n{method_name}:")
        print(f"  AUC:    {np.mean(aucs):.4f} +/- {np.std(aucs):.4f}")
        print(f"  BEDROC: {np.mean(bedrocs):.4f} +/- {np.std(bedrocs):.4f}")
        print(f"  EF@1%:  {np.mean(ef1s):.1f} +/- {np.std(ef1s):.1f}")

    print("\n" + "-" * 70)
    print("Per-fold AUC:")
    header = f"{'Fold':<6}"
    for method_name in loro_results:
        header += f"{method_name:<18}"
    print(header)

    n_folds = len(next(iter(loro_results.values())))
    for fold_idx in range(n_folds):
        row = f"{fold_idx:<6}"
        for method_name, folds in loro_results.items():
            row += f"{folds[fold_idx]['roc_auc']:<18.4f}"
        print(row)


def main():
    parser = argparse.ArgumentParser(description='Phase 2 validation experiment')
    parser.add_argument('--quick', action='store_true',
                        help='Skip greedy feature selection')
    parser.add_argument('--n-conformers', type=int, default=15,
                        help='Conformers per molecule (default: 15)')
    args = parser.parse_args()

    start_time = time.time()

    # Load dataset
    ref_mols, actives, decoys = load_ccr2_dataset()

    config = EvaluationConfig(
        n_conformers=args.n_conformers,
        opt_param=0.5,
        max_preiters=10,
        max_postiters=30,
    )

    # Run LORO CV
    logger.info("Starting LORO CV with %d references", len(ref_mols))
    loro_results = run_loro_cv(ref_mols, actives, decoys, config)

    # Run greedy selection (unless --quick)
    greedy_result = None
    if not args.quick:
        greedy_result = run_greedy_selection(ref_mols, actives, decoys)

    # Print results
    print_comparison_table(loro_results)

    if greedy_result:
        print("\n" + "=" * 70)
        print("GREEDY FEATURE SELECTION")
        print("=" * 70)
        print(f"Selected {len(greedy_result['selected_indices'])} features")
        print(f"Final AUC: {greedy_result['best_auc']:.4f}")
        print(f"Final BEDROC: {greedy_result['best_bedroc']:.4f}")
        print(f"Time: {greedy_result['wall_time_sec']:.1f}s")
        print("\nSelection curve:")
        for h in greedy_result['selection_history']:
            print(f"  Step {h['step']}: +feature {h['feature_idx']}, "
                  f"AUC={h['auc']:.4f}")

    total_time = time.time() - start_time
    logger.info("Total experiment time: %.1f min", total_time / 60)

    # Save results
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = os.path.join(
        os.path.dirname(__file__), 'results'
    )
    os.makedirs(output_dir, exist_ok=True)

    output = {
        'timestamp': timestamp,
        'config': {
            'n_conformers': args.n_conformers,
            'opt_param': config.opt_param,
        },
        'loro_cv': loro_results,
        'greedy_selection': greedy_result,
        'total_time_sec': total_time,
    }

    output_path = os.path.join(output_dir, f'validation_phases_{timestamp}.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    logger.info("Results saved to %s", output_path)


if __name__ == '__main__':
    main()

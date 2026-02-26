#!/usr/bin/env python
"""Validate fixed scoring methods against broken baselines on CCR2 dataset.

Runs side-by-side comparison:
  1. Hungarian (broken, unaligned) vs Hungarian (fixed, aligned)
  2. OT (broken, unaligned) vs OT (fixed, aligned)
  3. RRF ensemble (with improved two-tier gate)
  4. ICP Point Cloud (new)
  5. Reference 3D baseline (for comparison)

Usage:
    python experiments/exp_validate_scoring_fixes.py
    python experiments/exp_validate_scoring_fixes.py --n-conformers 15
"""

import os
import sys
import time
import json
import argparse
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.screening_metrics import calculate_all_metrics
from pharmacophore.pharmacophore import Pharmacophore
from pharmacophore.evaluation import UnifiedEvaluator, EvaluationConfig


def load_ccr2_dataset():
    """Load CCR2 reference, active, and decoy molecules."""
    data_dir = Path(__file__).parent.parent / 'tutorials' / 'data'

    # Load references from SDF (already have 3D coords)
    ref_path = data_dir / 'CCR2_reference_ligands.sdf'
    supplier = Chem.SDMolSupplier(str(ref_path), removeHs=False)
    refs = [m for m in supplier if m is not None and m.GetNumConformers() > 0]
    print(f"  References: {len(refs)}")

    # Load actives (SMILES column is 'SMILES')
    actives_df = pd.read_csv(data_dir / 'actives_ccr2_N75.csv')
    actives = []
    for smi in actives_df['SMILES']:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            if mol.GetNumConformers() > 0:
                actives.append(mol)
    print(f"  Actives:    {len(actives)}")

    # Load decoys (SMILES column is 'Smiles' — capital S)
    decoys_df = pd.read_csv(data_dir / 'decoys_ccr2_N500.csv')
    decoys = []
    for smi in decoys_df['Smiles']:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            if mol.GetNumConformers() > 0:
                decoys.append(mol)
    print(f"  Decoys:     {len(decoys)}")

    return refs, actives, decoys


def run_broken_hungarian(refs, actives, decoys):
    """Run the broken (unaligned) Hungarian matching."""
    from pharmacophore.hungarian_matching import pharmacophore_distance

    p = Pharmacophore()
    ref_features = [p.calc_pharm(mol=r) for r in refs]
    ref_features = [f for f in ref_features if f]

    def score_mol(mol):
        try:
            feats = p.calc_pharm(mol=mol)
            if not feats:
                return 0.0
            best = 0.0
            for rf in ref_features:
                dist = pharmacophore_distance(rf, feats)
                sim = 1.0 / (1.0 + dist)
                best = max(best, sim)
            return best
        except Exception:
            return 0.0

    active_scores = [score_mol(m) for m in actives]
    decoy_scores = [score_mol(m) for m in decoys]
    return active_scores, decoy_scores


def run_fixed_hungarian(refs, actives, decoys):
    """Run the fixed (aligned) Hungarian matching."""
    from pharmacophore.hungarian_matching import pharmacophore_similarity_aligned

    def score_mol(mol):
        try:
            best = 0.0
            for ref in refs:
                sim = pharmacophore_similarity_aligned(mol, ref, alpha=0.3)
                best = max(best, sim)
            return best
        except Exception:
            return 0.0

    active_scores = [score_mol(m) for m in actives]
    decoy_scores = [score_mol(m) for m in decoys]
    return active_scores, decoy_scores


def run_broken_ot(refs, actives, decoys):
    """Run the broken (unaligned) OT scoring."""
    from pharmacophore.ot_scoring import wasserstein_similarity

    p = Pharmacophore()
    ref_features = [p.calc_pharm(mol=r) for r in refs]
    ref_features = [f for f in ref_features if f]

    def score_mol(mol):
        try:
            feats = p.calc_pharm(mol=mol)
            if not feats:
                return 0.0
            best = 0.0
            for rf in ref_features:
                sim = wasserstein_similarity(rf, feats)
                best = max(best, sim)
            return best
        except Exception:
            return 0.0

    active_scores = [score_mol(m) for m in actives]
    decoy_scores = [score_mol(m) for m in decoys]
    return active_scores, decoy_scores


def run_fixed_ot(refs, actives, decoys):
    """Run the fixed (aligned) OT scoring."""
    from pharmacophore.ot_scoring import wasserstein_similarity_aligned

    def score_mol(mol):
        try:
            best = 0.0
            for ref in refs:
                sim = wasserstein_similarity_aligned(mol, ref, blend_alpha=0.3)
                best = max(best, sim)
            return best
        except Exception:
            return 0.0

    active_scores = [score_mol(m) for m in actives]
    decoy_scores = [score_mol(m) for m in decoys]
    return active_scores, decoy_scores


def run_icp(refs, actives, decoys):
    """Run the new ICP point cloud scoring."""
    from pharmacophore.point_cloud_alignment import point_cloud_similarity_aligned

    def score_mol(mol):
        try:
            best = 0.0
            for ref in refs:
                sim = point_cloud_similarity_aligned(
                    mol, ref, alpha=0.3, lambda_color=0.5, sigma=2.0
                )
                best = max(best, sim)
            return best
        except Exception:
            return 0.0

    active_scores = [score_mol(m) for m in actives]
    decoy_scores = [score_mol(m) for m in decoys]
    return active_scores, decoy_scores


def run_rrf_ensemble(refs, actives, decoys, n_conformers=15):
    """Run the RRF ensemble with improved gate."""
    evaluator = UnifiedEvaluator(refs, actives, decoys, random_state=42)
    config = EvaluationConfig(
        tolerance=2.0, occurrence=0.5,
        n_conformers=n_conformers, opt_param=0.5,
    )
    result = evaluator.evaluate_ensemble(config, quick_mode=False)
    return result


def run_3d_baseline(refs, actives, decoys, n_conformers=15):
    """Run the reference 3D baseline for comparison."""
    evaluator = UnifiedEvaluator(refs, actives, decoys, random_state=42)
    config = EvaluationConfig(
        tolerance=2.0, occurrence=0.5,
        n_conformers=n_conformers, opt_param=0.5,
    )
    result = evaluator.evaluate(config)
    return result


def main():
    parser = argparse.ArgumentParser(description='Validate fixed scoring methods')
    parser.add_argument('--n-conformers', type=int, default=15)
    args = parser.parse_args()

    print("=" * 70)
    print("  SCORING FIX VALIDATION — CCR2 Dataset")
    print("=" * 70)

    print("\nLoading CCR2 dataset...")
    refs, actives, decoys = load_ccr2_dataset()

    results = {}

    # 1. 3D Baseline
    print("\n[1/7] Reference 3D baseline...")
    t = time.time()
    baseline_result = run_3d_baseline(refs, actives, decoys, args.n_conformers)
    elapsed = time.time() - t
    results['3D Baseline'] = {
        'roc_auc': baseline_result.roc_auc,
        'bedroc': baseline_result.bedroc,
        'ef_1': baseline_result.ef_1,
        'ef_5': baseline_result.ef_5,
        'time_sec': elapsed,
    }
    print(f"  AUC={baseline_result.roc_auc:.4f}  BEDROC={baseline_result.bedroc:.4f}  "
          f"EF@1%={baseline_result.ef_1:.1f}  ({elapsed:.1f}s)")

    # 2-3. Hungarian
    print("\n[2/7] Hungarian (broken, unaligned)...")
    t = time.time()
    a_scores, d_scores = run_broken_hungarian(refs, actives, decoys)
    elapsed = time.time() - t
    y_true = [1] * len(a_scores) + [0] * len(d_scores)
    metrics = calculate_all_metrics(y_true, a_scores + d_scores)
    results['Hungarian (broken)'] = {**metrics, 'time_sec': elapsed}
    print(f"  AUC={metrics['roc_auc']:.4f}  BEDROC={metrics['bedroc']:.4f}  ({elapsed:.1f}s)")

    print("\n[3/7] Hungarian (fixed, aligned)...")
    t = time.time()
    a_scores, d_scores = run_fixed_hungarian(refs, actives, decoys)
    elapsed = time.time() - t
    y_true = [1] * len(a_scores) + [0] * len(d_scores)
    metrics = calculate_all_metrics(y_true, a_scores + d_scores)
    results['Hungarian (fixed)'] = {**metrics, 'time_sec': elapsed}
    print(f"  AUC={metrics['roc_auc']:.4f}  BEDROC={metrics['bedroc']:.4f}  ({elapsed:.1f}s)")

    # 4-5. OT
    print("\n[4/7] OT (broken, unaligned)...")
    t = time.time()
    a_scores, d_scores = run_broken_ot(refs, actives, decoys)
    elapsed = time.time() - t
    y_true = [1] * len(a_scores) + [0] * len(d_scores)
    metrics = calculate_all_metrics(y_true, a_scores + d_scores)
    results['OT (broken)'] = {**metrics, 'time_sec': elapsed}
    print(f"  AUC={metrics['roc_auc']:.4f}  BEDROC={metrics['bedroc']:.4f}  ({elapsed:.1f}s)")

    print("\n[5/7] OT (fixed, aligned)...")
    t = time.time()
    a_scores, d_scores = run_fixed_ot(refs, actives, decoys)
    elapsed = time.time() - t
    y_true = [1] * len(a_scores) + [0] * len(d_scores)
    metrics = calculate_all_metrics(y_true, a_scores + d_scores)
    results['OT (fixed)'] = {**metrics, 'time_sec': elapsed}
    print(f"  AUC={metrics['roc_auc']:.4f}  BEDROC={metrics['bedroc']:.4f}  ({elapsed:.1f}s)")

    # 6. RRF ensemble
    print("\n[6/7] RRF Ensemble (improved gate)...")
    t = time.time()
    rrf_result = run_rrf_ensemble(refs, actives, decoys, args.n_conformers)
    elapsed = time.time() - t
    results['RRF Ensemble (gated)'] = {
        'roc_auc': rrf_result.roc_auc,
        'bedroc': rrf_result.bedroc,
        'ef_1': rrf_result.ef_1,
        'ef_5': rrf_result.ef_5,
        'time_sec': elapsed,
    }
    print(f"  AUC={rrf_result.roc_auc:.4f}  BEDROC={rrf_result.bedroc:.4f}  ({elapsed:.1f}s)")

    # 7. ICP
    print("\n[7/7] ICP Point Cloud (new)...")
    t = time.time()
    a_scores, d_scores = run_icp(refs, actives, decoys)
    elapsed = time.time() - t
    y_true = [1] * len(a_scores) + [0] * len(d_scores)
    metrics = calculate_all_metrics(y_true, a_scores + d_scores)
    results['ICP Point Cloud'] = {**metrics, 'time_sec': elapsed}
    print(f"  AUC={metrics['roc_auc']:.4f}  BEDROC={metrics['bedroc']:.4f}  ({elapsed:.1f}s)")

    # Summary table
    print("\n" + "=" * 70)
    print("  RESULTS SUMMARY")
    print("=" * 70)
    print(f"\n{'Approach':<30} {'AUC':>8} {'BEDROC':>8} {'EF@1%':>8} {'Time':>8}")
    print("-" * 70)
    for name, r in results.items():
        auc = r.get('roc_auc', 0.5)
        bedroc = r.get('bedroc', 0.0)
        ef1 = r.get('ef_1', 0.0)
        t = r.get('time_sec', 0.0)
        print(f"{name:<30} {auc:>8.4f} {bedroc:>8.4f} {ef1:>8.1f} {t:>7.1f}s")

    # Improvement summary
    print("\n" + "-" * 70)
    print("  IMPROVEMENTS")
    print("-" * 70)

    h_broken = results.get('Hungarian (broken)', {}).get('roc_auc', 0.5)
    h_fixed = results.get('Hungarian (fixed)', {}).get('roc_auc', 0.5)
    print(f"  Hungarian: {h_broken:.4f} -> {h_fixed:.4f} ({h_fixed-h_broken:+.4f})")

    ot_broken = results.get('OT (broken)', {}).get('roc_auc', 0.5)
    ot_fixed = results.get('OT (fixed)', {}).get('roc_auc', 0.5)
    print(f"  OT:        {ot_broken:.4f} -> {ot_fixed:.4f} ({ot_fixed-ot_broken:+.4f})")

    rrf_auc = results.get('RRF Ensemble (gated)', {}).get('roc_auc', 0.5)
    baseline_auc = results.get('3D Baseline', {}).get('roc_auc', 0.5)
    print(f"  RRF:       {rrf_auc:.4f} (baseline: {baseline_auc:.4f})")

    icp_auc = results.get('ICP Point Cloud', {}).get('roc_auc', 0.5)
    print(f"  ICP:       {icp_auc:.4f} (new)")

    # Save results
    results_dir = Path(__file__).parent / 'results'
    results_dir.mkdir(exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_path = results_dir / f'scoring_fixes_{timestamp}.json'
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {output_path}")


if __name__ == '__main__':
    main()

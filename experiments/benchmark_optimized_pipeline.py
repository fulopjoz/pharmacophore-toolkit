#!/usr/bin/env python
"""
Ablation benchmark: measure the impact of each pipeline optimization.

Tests the optimized 3D shape matching pipeline against the baseline on CCR2.
Each optimization is tested incrementally to measure individual contribution.

Configurations tested:
  1. baseline_old    — AlignMol, no PrepareConformer, 25 conformers (old defaults)
  2. baseline_new    — AlignShapes + PrepareConformer, 10 conformers (new defaults)
  3. early_stop      — + early termination at threshold 1.8
  4. fast_iters      — + reduced iterations (5/15)
  5. fewer_confs     — baseline_new with only 5 conformers
  6. ref_weighted    — + reference weighting by per-ref AUC

Usage:
    python experiments/benchmark_optimized_pipeline.py
"""

import sys
from pathlib import Path
import time
import json
from datetime import datetime
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol, PrepareConformer, AlignShapes
from rdkit.Chem.rdmolfiles import SDMolSupplier
from sklearn.metrics import roc_auc_score
from joblib import Parallel, delayed, Memory
import pandas as pd

from pharmacophore.screening_metrics import calculate_all_metrics


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
CACHE_DIR = Path(__file__).parent / '.cache'
CACHE_DIR.mkdir(exist_ok=True)
memory = Memory(CACHE_DIR, verbose=0)

RESULTS_DIR = Path(__file__).parent / 'results'
RESULTS_DIR.mkdir(exist_ok=True)

N_JOBS = -1


@dataclass
class BenchResult:
    name: str
    roc_auc: float
    bedroc: float
    ef_1: float
    ef_5: float
    wall_time_sec: float
    n_align_calls: int = 0
    params: Optional[Dict] = None


# ---------------------------------------------------------------------------
# Data loading (cached)
# ---------------------------------------------------------------------------
@memory.cache
def generate_conformers_cached(smiles, n_conformers, random_seed=42):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.pruneRmsThresh = 0.5
    AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)
    if mol_h.GetNumConformers() > 0:
        return mol_h
    return None


def load_dataset(n_conformers=10, verbose=True):
    base = Path(__file__).parent.parent / 'tutorials' / 'data'

    refs = []
    supplier = SDMolSupplier(str(base / 'CCR2_reference_ligands.sdf'), removeHs=False)
    for mol in supplier:
        if mol and mol.GetNumConformers() > 0:
            refs.append(mol)

    actives_df = pd.read_csv(base / 'actives_ccr2_N75.csv')
    active_smiles = actives_df['SMILES'].tolist()
    decoys_df = pd.read_csv(base / 'decoys_ccr2_N500.csv')
    decoy_smiles = decoys_df['Smiles'].tolist()

    if verbose:
        print(f"  Loading {len(active_smiles)} actives + {len(decoy_smiles)} decoys "
              f"with {n_conformers} conformers each...")

    actives = Parallel(n_jobs=N_JOBS, backend='loky')(
        delayed(generate_conformers_cached)(smi, n_conformers) for smi in active_smiles
    )
    actives = [m for m in actives if m is not None]

    decoys = Parallel(n_jobs=N_JOBS, backend='loky')(
        delayed(generate_conformers_cached)(smi, n_conformers) for smi in decoy_smiles
    )
    decoys = [m for m in decoys if m is not None]

    if verbose:
        print(f"  Loaded: {len(refs)} refs, {len(actives)} actives, {len(decoys)} decoys\n")

    return refs, actives, decoys


# ---------------------------------------------------------------------------
# Prepare references
# ---------------------------------------------------------------------------
def prepare_refs(refs):
    """Prepare reference molecules with AddHs and pre-computed shapes."""
    prepared = []
    shapes = []
    for mol in refs:
        mol_h = Chem.AddHs(mol)
        if mol_h.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol_h, randomSeed=42)
        if mol_h.GetNumConformers() > 0:
            prepared.append(mol_h)
            try:
                s = PrepareConformer(mol_h, confId=0)
                shapes.append(s)
            except Exception:
                shapes.append(None)
    return prepared, shapes


# ---------------------------------------------------------------------------
# Scoring functions (different optimization levels)
# ---------------------------------------------------------------------------
def score_baseline_old(mol, refs, opt_param=0.5):
    """Old baseline: AlignMol with default RDKit params, no PrepareConformer."""
    n_calls = 0
    ref_scores = []
    for ref in refs:
        best = 0.0
        for cid in range(mol.GetNumConformers()):
            try:
                s, c = AlignMol(
                    ref=ref, probe=mol, probeConfId=cid,
                    useColors=True, opt_param=opt_param
                )
                best = max(best, s + c)
                n_calls += 1
            except Exception:
                continue
        ref_scores.append(best)
    return max(ref_scores) if ref_scores else 0.0, n_calls


def score_alignshapes(mol, refs, ref_shapes, opt_param=0.5,
                      max_preiters=10, max_postiters=30,
                      early_stop=2.0):
    """Optimized: PrepareConformer + AlignShapes with configurable params."""
    n_calls = 0
    n_confs = mol.GetNumConformers()

    # Pre-compute probe shapes once
    probe_shapes = []
    for cid in range(n_confs):
        try:
            ps = PrepareConformer(mol, confId=cid)
            probe_shapes.append(ps)
        except Exception:
            probe_shapes.append(None)

    ref_scores = []
    for ref_idx, ref in enumerate(refs):
        ref_shape = ref_shapes[ref_idx] if ref_idx < len(ref_shapes) else None
        best = 0.0

        for cid in range(n_confs):
            ps = probe_shapes[cid]
            try:
                if ref_shape is not None and ps is not None:
                    s, c, _ = AlignShapes(
                        ref_shape, ps,
                        opt_param=opt_param,
                        max_preiters=max_preiters,
                        max_postiters=max_postiters
                    )
                else:
                    s, c = AlignMol(
                        ref=ref, probe=mol, probeConfId=cid,
                        useColors=True, opt_param=opt_param,
                        max_preiters=max_preiters,
                        max_postiters=max_postiters
                    )
                best = max(best, s + c)
                n_calls += 1
            except Exception:
                continue

            if best >= early_stop:
                break

        ref_scores.append(best)
    return max(ref_scores) if ref_scores else 0.0, n_calls


def score_ref_weighted(mol, refs, ref_shapes, ref_weights, opt_param=0.5,
                       max_preiters=10, max_postiters=30, early_stop=1.8):
    """Weighted aggregation: score * ref_weight instead of simple max."""
    n_calls = 0
    n_confs = mol.GetNumConformers()

    probe_shapes = []
    for cid in range(n_confs):
        try:
            ps = PrepareConformer(mol, confId=cid)
            probe_shapes.append(ps)
        except Exception:
            probe_shapes.append(None)

    ref_scores = []
    for ref_idx, ref in enumerate(refs):
        ref_shape = ref_shapes[ref_idx] if ref_idx < len(ref_shapes) else None
        best = 0.0

        for cid in range(n_confs):
            ps = probe_shapes[cid]
            try:
                if ref_shape is not None and ps is not None:
                    s, c, _ = AlignShapes(
                        ref_shape, ps,
                        opt_param=opt_param,
                        max_preiters=max_preiters,
                        max_postiters=max_postiters
                    )
                else:
                    s, c = AlignMol(
                        ref=ref, probe=mol, probeConfId=cid,
                        useColors=True, opt_param=opt_param
                    )
                best = max(best, s + c)
                n_calls += 1
            except Exception:
                continue
            if best >= early_stop:
                break
        ref_scores.append(best)

    # Weighted aggregation
    weighted = sum(s * w for s, w in zip(ref_scores, ref_weights))
    return weighted, n_calls


# ---------------------------------------------------------------------------
# Run a single configuration
# ---------------------------------------------------------------------------
def run_config(name, all_mols, y_true, score_fn, **kwargs):
    """Run a scoring configuration and measure time + metrics."""
    params_info = kwargs.pop('params_info', {})
    print(f"  [{name}] Scoring {len(all_mols)} molecules...", end=' ', flush=True)
    t0 = time.perf_counter()

    total_calls = 0
    scores = []
    for mol in all_mols:
        score, n_calls = score_fn(mol, **kwargs)
        scores.append(score)
        total_calls += n_calls

    wall = time.perf_counter() - t0
    metrics = calculate_all_metrics(y_true, scores)

    print(f"AUC={metrics['roc_auc']:.4f}, BEDROC={metrics['bedroc']:.4f}, "
          f"time={wall:.1f}s, calls={total_calls}")

    return BenchResult(
        name=name,
        roc_auc=metrics['roc_auc'],
        bedroc=metrics['bedroc'],
        ef_1=metrics['ef_1'],
        ef_5=metrics['ef_5'],
        wall_time_sec=round(wall, 2),
        n_align_calls=total_calls,
        params=kwargs.get('params_info', {})
    )


# ---------------------------------------------------------------------------
# Compute reference weights
# ---------------------------------------------------------------------------
def compute_ref_weights(refs, ref_shapes, all_mols, y_true, opt_param=0.5):
    """Compute per-reference AUC and derive weights."""
    print("\n  Computing per-reference weights...")
    ref_aucs = []
    for ref_idx in range(len(refs)):
        ref_shape = ref_shapes[ref_idx]
        scores = []
        for mol in all_mols:
            best = 0.0
            for cid in range(mol.GetNumConformers()):
                try:
                    ps = PrepareConformer(mol, confId=cid)
                    if ref_shape is not None and ps is not None:
                        s, c, _ = AlignShapes(ref_shape, ps, opt_param=opt_param)
                    else:
                        s, c = AlignMol(
                            ref=refs[ref_idx], probe=mol, probeConfId=cid,
                            useColors=True, opt_param=opt_param
                        )
                    best = max(best, s + c)
                except Exception:
                    continue
            scores.append(best)
        auc = roc_auc_score(y_true, scores)
        ref_aucs.append(auc)
        print(f"    Ref {ref_idx}: AUC={auc:.4f}")

    # Weight by (AUC - 0.5), zero out anti-discriminative refs
    raw_weights = [max(0, a - 0.5) for a in ref_aucs]
    total = sum(raw_weights)
    if total > 0:
        weights = [w / total for w in raw_weights]
    else:
        weights = [1.0 / len(refs)] * len(refs)

    print(f"    Weights: {[f'{w:.3f}' for w in weights]}\n")
    return weights, ref_aucs


# ---------------------------------------------------------------------------
# Main benchmark
# ---------------------------------------------------------------------------
def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"\n{'='*70}")
    print(f"  ABLATION BENCHMARK: Optimized 3D Shape Pipeline")
    print(f"  {timestamp}")
    print(f"{'='*70}\n")

    # Load data with 25 conformers (max we'll test)
    print("Loading dataset with 25 conformers (for baseline_old)...")
    refs_25, actives_25, decoys_25 = load_dataset(n_conformers=25)
    all_mols_25 = actives_25 + decoys_25
    y_true = [1] * len(actives_25) + [0] * len(decoys_25)

    # Also load with 10 conformers
    print("Loading dataset with 10 conformers (for optimized configs)...")
    refs_10, actives_10, decoys_10 = load_dataset(n_conformers=10)
    all_mols_10 = actives_10 + decoys_10

    # And 5 conformers
    print("Loading dataset with 5 conformers (for few-confs test)...")
    _, actives_5, decoys_5 = load_dataset(n_conformers=5)
    all_mols_5 = actives_5 + decoys_5

    # Prepare references
    prep_refs, ref_shapes = prepare_refs(refs_25)

    results = []

    print(f"\n{'='*70}")
    print(f"  Running ablation configs...")
    print(f"{'='*70}\n")

    # Config 1: baseline_old — AlignMol, 25 conformers, no tricks
    r = run_config(
        "1_baseline_old",
        all_mols_25, y_true,
        score_baseline_old,
        refs=prep_refs, opt_param=0.5,
        params_info={'n_conformers': 25, 'method': 'AlignMol', 'iters': 'default(10/30)'}
    )
    results.append(r)

    # Config 2: baseline_new — AlignShapes + PrepareConformer, 10 conformers
    r = run_config(
        "2_alignshapes_10conf",
        all_mols_10, y_true,
        score_alignshapes,
        refs=prep_refs, ref_shapes=ref_shapes,
        opt_param=0.5, max_preiters=10, max_postiters=30, early_stop=2.0,
        params_info={'n_conformers': 10, 'method': 'AlignShapes', 'iters': '10/30'}
    )
    results.append(r)

    # Config 3: + early termination at 1.8
    r = run_config(
        "3_early_stop_1.8",
        all_mols_10, y_true,
        score_alignshapes,
        refs=prep_refs, ref_shapes=ref_shapes,
        opt_param=0.5, max_preiters=10, max_postiters=30, early_stop=1.8,
        params_info={'n_conformers': 10, 'method': 'AlignShapes', 'iters': '10/30', 'early_stop': 1.8}
    )
    results.append(r)

    # Config 4: fast iterations (5/15)
    r = run_config(
        "4_fast_iters_5_15",
        all_mols_10, y_true,
        score_alignshapes,
        refs=prep_refs, ref_shapes=ref_shapes,
        opt_param=0.5, max_preiters=5, max_postiters=15, early_stop=1.8,
        params_info={'n_conformers': 10, 'method': 'AlignShapes', 'iters': '5/15', 'early_stop': 1.8}
    )
    results.append(r)

    # Config 5: only 5 conformers (extreme speed)
    r = run_config(
        "5_5_conformers",
        all_mols_5, y_true,
        score_alignshapes,
        refs=prep_refs, ref_shapes=ref_shapes,
        opt_param=0.5, max_preiters=10, max_postiters=30, early_stop=1.8,
        params_info={'n_conformers': 5, 'method': 'AlignShapes', 'iters': '10/30', 'early_stop': 1.8}
    )
    results.append(r)

    # Config 6: reference weighting
    ref_weights, ref_aucs = compute_ref_weights(
        prep_refs, ref_shapes, all_mols_10, y_true
    )

    r = run_config(
        "6_ref_weighted",
        all_mols_10, y_true,
        score_ref_weighted,
        refs=prep_refs, ref_shapes=ref_shapes,
        ref_weights=ref_weights,
        opt_param=0.5, max_preiters=10, max_postiters=30, early_stop=1.8,
        params_info={'n_conformers': 10, 'method': 'AlignShapes+weighted', 'iters': '10/30',
                     'early_stop': 1.8, 'ref_weights': [round(w, 3) for w in ref_weights],
                     'ref_aucs': [round(a, 4) for a in ref_aucs]}
    )
    results.append(r)

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    print(f"\n{'='*70}")
    print(f"  RESULTS SUMMARY")
    print(f"{'='*70}\n")

    baseline_time = results[0].wall_time_sec
    baseline_auc = results[0].roc_auc

    print(f"{'Config':<30s} {'AUC':>7s} {'BEDROC':>8s} {'EF@1%':>7s} "
          f"{'Time':>7s} {'Speedup':>8s} {'Calls':>8s} {'dAUC':>7s}")
    print('-' * 84)

    for r in results:
        speedup = baseline_time / r.wall_time_sec if r.wall_time_sec > 0 else 0
        d_auc = r.roc_auc - baseline_auc
        print(f"{r.name:<30s} {r.roc_auc:>7.4f} {r.bedroc:>8.4f} {r.ef_1:>7.2f} "
              f"{r.wall_time_sec:>6.1f}s {speedup:>7.1f}x {r.n_align_calls:>8d} "
              f"{d_auc:>+7.4f}")

    # Save results
    out = {
        'timestamp': timestamp,
        'dataset': 'CCR2',
        'results': [asdict(r) for r in results],
        'ref_weights': [round(w, 4) for w in ref_weights],
        'ref_aucs': [round(a, 4) for a in ref_aucs]
    }
    out_path = RESULTS_DIR / f'ablation_pipeline_{timestamp}.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nResults saved to: {out_path}")


if __name__ == '__main__':
    main()

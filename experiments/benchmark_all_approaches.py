#!/usr/bin/env python
"""Comprehensive benchmark of ALL 17 pharmacophore approaches on CCR2 dataset.

Compares clustering, scoring, optimization, and meta approaches:

Clustering Methods (5):
  1. Baseline Agglomerative     - Standard hierarchical clustering
  2. Spatial-Scaled             - Per-feature-type spatial scaling
  3. Ensemble Consensus         - Stability-aware multi-run voting
  4. K-Means Clustering         - K-means with estimated cluster count
  5. Grid-Based Binning         - Ultra-fast grid binning (O(N))

Scoring Methods (6):
  6. Pharm2D Fingerprints       - 2D pharmacophore fingerprints
  7. Reference Ensemble (3D)    - rdShapeAlign to reference molecules
  8. Hybrid 2D+3D              - Weighted combination (best AUC)
  9. Hungarian Matching         - Aligned feature assignment distance
  10. Optimal Transport         - Aligned Wasserstein pharmacophore distance
  17. Point Cloud ICP           - Colored ICP alignment scoring

Optimization Methods (4):
  11. Optuna GP (Bayesian)      - Gaussian Process sampler
  12. Optuna NSGA-II            - Evolutionary multi-objective
  13. HypoGen 3-Phase           - Constructive-subtractive-refinement
  14. Multi-Fidelity BO         - Combinatorial optimizer

Meta Methods (2):
  15. Strategy Selector         - Tournament across strategies
  16. Ensemble RRF Scoring      - Cascading multi-level with RRF

Usage:
    python experiments/benchmark_all_approaches.py
    python experiments/benchmark_all_approaches.py --n-conformers 10 --n-trials 30
    python experiments/benchmark_all_approaches.py --skip-slow  # skip optimizers
"""

import os
import sys
import time
import json
import argparse
import itertools
import traceback
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Any

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, SDMolSupplier
from joblib import Memory, Parallel, delayed
from sklearn.metrics import roc_auc_score

from pharmacophore.evaluation import UnifiedEvaluator, EvaluationConfig, EvaluationResult
from pharmacophore.consensus import PharmacophoreConsensus
from pharmacophore.mol_converter import PharmacophoreToMol
from pharmacophore.pharm2d_scoring import Pharm2DScorer
from pharmacophore.screening_metrics import calculate_all_metrics

try:
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    HAS_OPTUNA = True
except ImportError:
    HAS_OPTUNA = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
N_JOBS = -1
CACHE_DIR = Path(__file__).parent / '.cache'
CACHE_DIR.mkdir(exist_ok=True)
memory = Memory(CACHE_DIR, verbose=0)

RESULTS_DIR = Path(__file__).parent / 'results'
RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR = RESULTS_DIR / 'figures'
FIGURES_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Approach result container
# ---------------------------------------------------------------------------
@dataclass
class ApproachResult:
    """Unified result container for all approaches."""
    id: int
    approach: str
    category: str  # clustering, scoring, optimization, meta
    roc_auc: float = 0.5
    bedroc: float = 0.0
    ef_1: float = 0.0
    ef_5: float = 0.0
    n_features: int = 0
    wall_time_sec: float = 0.0
    best_params: Optional[Dict] = None
    error: Optional[str] = None
    pros: str = ''
    cons: str = ''
    # Clustering quality metrics (approaches 1-5 only)
    silhouette: float = float('nan')
    sdbw: float = float('nan')
    n_clusters: int = 0

    def to_dict(self) -> dict:
        d = asdict(self)
        if d['best_params'] is None:
            d['best_params'] = {}
        # JSON can't serialize nan — use None instead
        for key in ('silhouette', 'sdbw'):
            if isinstance(d[key], float) and np.isnan(d[key]):
                d[key] = None
        return d


# ---------------------------------------------------------------------------
# Data loading (cached)
# ---------------------------------------------------------------------------
@memory.cache
def generate_conformers_cached(smiles, n_conformers, random_seed=42):
    """Cached conformer generation (ETKDGv3)."""
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


def load_ccr2_dataset(n_conformers=10, verbose=True):
    """Load complete CCR2 dataset with parallel conformer generation."""
    base = Path(__file__).parent.parent / 'tutorials' / 'data'

    # References (already have 3D coords in SDF)
    refs = []
    supplier = SDMolSupplier(str(base / 'CCR2_reference_ligands.sdf'), removeHs=False)
    for mol in supplier:
        if mol and mol.GetNumConformers() > 0:
            refs.append(mol)

    # Actives
    actives_df = pd.read_csv(base / 'actives_ccr2_N75.csv')
    active_smiles = actives_df['SMILES'].tolist()

    # Decoys
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
        n_total = len(actives) + len(decoys)
        print(f"  Total: {len(actives)} actives + {len(decoys)} decoys = {n_total} molecules\n")

    return refs, actives, decoys


# ---------------------------------------------------------------------------
# Helper: clustering quality from consensus metadata
# ---------------------------------------------------------------------------
def _compute_clustering_quality(metadata: dict) -> dict:
    """Compute silhouette and S_Dbw from consensus metadata.

    Args:
        metadata: Dict from generate_consensus(return_metadata=True).
            Keys are feature types, values have 'coordinates' and 'labels'.

    Returns:
        Dict with 'silhouette', 'sdbw', 'n_clusters'.
    """
    from pharmacophore.evaluation import compute_sdbw
    from sklearn.metrics import silhouette_score

    sdbw_val = compute_sdbw(metadata)

    # Aggregate coordinates and labels across feature types for silhouette
    all_coords = []
    all_labels = []
    label_offset = 0
    for feat_type, info in metadata.items():
        coords = info['coordinates']
        labels = info['labels']
        if len(coords) == 0:
            continue
        all_coords.append(coords)
        all_labels.append(labels + label_offset)
        label_offset += labels.max() + 1 if len(labels) > 0 else 0

    if not all_coords:
        return {'silhouette': float('nan'), 'sdbw': sdbw_val, 'n_clusters': 0}

    coords = np.vstack(all_coords)
    labels = np.concatenate(all_labels)
    n_clusters = len(set(labels))

    if n_clusters < 2 or len(coords) < n_clusters + 1:
        return {'silhouette': float('nan'), 'sdbw': sdbw_val, 'n_clusters': n_clusters}

    try:
        sil = silhouette_score(coords, labels)
    except Exception:
        sil = float('nan')

    return {'silhouette': sil, 'sdbw': sdbw_val, 'n_clusters': n_clusters}


# ---------------------------------------------------------------------------
# Helper: generate consensus features + evaluate via UnifiedEvaluator
# ---------------------------------------------------------------------------
def _consensus_and_evaluate(
    refs, actives, decoys, seed,
    tolerance=2.0, occurrence=0.5, linkage='average',
    use_spatial_scaling=False, scoring_mode='reference',
    n_conformers=10, opt_param=0.5, aggregation='max',
    alpha=0.6, clustering_method='hierarchical',
    return_quality=False,
):
    """Generate consensus with given params and evaluate.

    Args:
        return_quality: If True, return a third element with clustering
            quality dict: {'silhouette', 'sdbw', 'n_clusters'}.
    """
    consensus = PharmacophoreConsensus(
        tolerance=tolerance,
        occurrence_threshold=occurrence,
        linkage=linkage,
        use_spatial_scaling=use_spatial_scaling,
    )

    # Override clustering method if not hierarchical
    if clustering_method != 'hierarchical':
        from pharmacophore.clustering_algorithms import cluster_features_with_labels

        class _AltConsensus(PharmacophoreConsensus):
            """Consensus with pluggable clustering algorithm.

            Uses native labels from each clustering method instead of
            post-hoc cdist+argmin reassignment (which erases k-means/grid
            differences and makes all methods produce identical results).
            """
            def __init__(self, tol, occ, link, method):
                super().__init__(tol, occ, link)
                self._method = method

            def _cluster_features(self, coordinates, mol_indices, feat_type=None):
                if not coordinates:
                    return np.array([]), {}
                coords_array = np.array(coordinates)
                n_molecules = len(set(int(m) for m in mol_indices))

                # filter_by_occurrence=False: return ALL cluster labels
                # (contiguous, >= 0).  The parent class handles occurrence
                # filtering via _filter_by_occurrence — we must not pre-filter
                # or the parent's _calculate_cluster_centroids will crash on -1.
                labels, _centroids = cluster_features_with_labels(
                    coords=coords_array,
                    tolerance=self.tolerance,
                    occurrence_threshold=self.occurrence_threshold,
                    n_molecules=n_molecules,
                    method=self._method,
                    linkage=self.linkage,
                    filter_by_occurrence=False,
                )
                if len(_centroids) == 0:
                    return np.array([]), {}

                cluster_to_mols = {}
                for cluster_id, mol_idx in zip(labels, mol_indices):
                    cid, mid = int(cluster_id), int(mol_idx)
                    if cid not in cluster_to_mols:
                        cluster_to_mols[cid] = []
                    if mid not in cluster_to_mols[cid]:
                        cluster_to_mols[cid].append(mid)
                return labels, cluster_to_mols

        consensus = _AltConsensus(tolerance, occurrence, linkage, clustering_method)

    if return_quality:
        features, metadata = consensus.generate_consensus(refs, return_metadata=True)
    else:
        features = consensus.generate_consensus(refs)
        metadata = None

    quality = None
    if return_quality and metadata:
        quality = _compute_clustering_quality(metadata)

    if len(features) < 2:
        if return_quality:
            return features, None, quality
        return features, None

    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)
    config = EvaluationConfig(
        tolerance=tolerance,
        occurrence=occurrence,
        shape_weight=0.5,
        opt_param=opt_param,
        linkage=linkage,
        n_conformers=n_conformers,
        scoring_mode=scoring_mode,
        aggregation=aggregation,
        alpha=alpha,
    )
    result = evaluator.evaluate(config)
    if return_quality:
        return features, result, quality
    return features, result


# ---------------------------------------------------------------------------
# APPROACH RUNNERS (1-16)
# ---------------------------------------------------------------------------

def run_01_baseline_agglomerative(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """1. Baseline Agglomerative Clustering."""
    start = time.time()
    features, result, quality = _consensus_and_evaluate(
        refs, actives, decoys, seed,
        tolerance=2.0, occurrence=0.5, linkage='average',
        scoring_mode='reference', n_conformers=n_conformers,
        return_quality=True,
    )
    elapsed = time.time() - start
    q = quality or {}

    if result is None:
        return ApproachResult(
            id=1, approach='Baseline Agglomerative', category='clustering',
            n_features=len(features), wall_time_sec=elapsed,
            error='<2 features',
            pros='Fast, deterministic, simple',
            cons='Sensitive to tolerance/occurrence parameters',
            silhouette=q.get('silhouette', float('nan')),
            sdbw=q.get('sdbw', float('nan')),
            n_clusters=q.get('n_clusters', 0),
        )

    return ApproachResult(
        id=1, approach='Baseline Agglomerative', category='clustering',
        roc_auc=result.roc_auc, bedroc=result.bedroc,
        ef_1=result.ef_1, ef_5=result.ef_5,
        n_features=result.n_features, wall_time_sec=elapsed,
        best_params={'tolerance': 2.0, 'occurrence': 0.5, 'linkage': 'average'},
        pros='Fast, deterministic, simple',
        cons='Sensitive to tolerance/occurrence parameters',
        silhouette=q.get('silhouette', float('nan')),
        sdbw=q.get('sdbw', float('nan')),
        n_clusters=q.get('n_clusters', 0),
    )


def run_02_spatial_scaled(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """2. Spatial-Scaled Clustering."""
    start = time.time()
    features, result, quality = _consensus_and_evaluate(
        refs, actives, decoys, seed,
        tolerance=2.0, occurrence=0.5, linkage='average',
        use_spatial_scaling=True,
        scoring_mode='reference', n_conformers=n_conformers,
        return_quality=True,
    )
    elapsed = time.time() - start
    q = quality or {}

    if result is None:
        return ApproachResult(
            id=2, approach='Spatial-Scaled', category='clustering',
            n_features=len(features), wall_time_sec=elapsed,
            error='<2 features',
            pros='Physics-based per-type tolerances',
            cons='Minimal improvement on small CCR2 dataset',
            silhouette=q.get('silhouette', float('nan')),
            sdbw=q.get('sdbw', float('nan')),
            n_clusters=q.get('n_clusters', 0),
        )

    return ApproachResult(
        id=2, approach='Spatial-Scaled', category='clustering',
        roc_auc=result.roc_auc, bedroc=result.bedroc,
        ef_1=result.ef_1, ef_5=result.ef_5,
        n_features=result.n_features, wall_time_sec=elapsed,
        best_params={'tolerance': 2.0, 'occurrence': 0.5, 'use_spatial_scaling': True},
        pros='Physics-based per-type tolerances',
        cons='Minimal improvement on small CCR2 dataset',
        silhouette=q.get('silhouette', float('nan')),
        sdbw=q.get('sdbw', float('nan')),
        n_clusters=q.get('n_clusters', 0),
    )


def run_03_ensemble_consensus(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """3. Ensemble Consensus (stability-aware multi-run)."""
    from pharmacophore.ensemble_consensus import EnsembleConsensus

    start = time.time()

    ec = EnsembleConsensus(
        n_runs=25,
        tolerance_range=(1.5, 2.5),
        occurrence_range=(0.3, 0.7),
        stability_threshold=0.3,
        random_state=seed,
    )
    features, stability = ec.generate_consensus_with_scores(refs)

    # Also generate standard consensus for quality metrics
    quality_consensus = PharmacophoreConsensus(tolerance=2.0, occurrence_threshold=0.5)
    _, meta = quality_consensus.generate_consensus(refs, return_metadata=True)
    q = _compute_clustering_quality(meta) if meta else {}

    if len(features) < 2:
        return ApproachResult(
            id=3, approach='Ensemble Consensus', category='clustering',
            n_features=len(features), wall_time_sec=time.time() - start,
            error='<2 features',
            pros='Robust, stable features identified',
            cons='Slower (25 internal runs)',
            silhouette=q.get('silhouette', float('nan')),
            sdbw=q.get('sdbw', float('nan')),
            n_clusters=q.get('n_clusters', 0),
        )

    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)
    result = evaluator.evaluate_feature_subset(features, 0.5, 0.5, n_conformers)
    elapsed = time.time() - start

    return ApproachResult(
        id=3, approach='Ensemble Consensus', category='clustering',
        roc_auc=result.roc_auc, bedroc=result.bedroc,
        ef_1=result.ef_1, ef_5=result.ef_5,
        n_features=result.n_features, wall_time_sec=elapsed,
        best_params={
            'n_runs': 25, 'stability_threshold': 0.3,
            'n_stable_features': sum(1 for s in stability if s >= 0.5),
            'avg_stability': float(np.mean(stability)) if stability else 0.0,
        },
        pros='Robust, stable features identified',
        cons='Slower (25 internal runs)',
        silhouette=q.get('silhouette', float('nan')),
        sdbw=q.get('sdbw', float('nan')),
        n_clusters=q.get('n_clusters', 0),
    )


def run_04_kmeans(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """4. K-Means Clustering."""
    start = time.time()
    features, result, quality = _consensus_and_evaluate(
        refs, actives, decoys, seed,
        tolerance=2.0, occurrence=0.5,
        clustering_method='kmeans',
        scoring_mode='reference', n_conformers=n_conformers,
        return_quality=True,
    )
    elapsed = time.time() - start
    q = quality or {}

    if result is None:
        return ApproachResult(
            id=4, approach='K-Means Clustering', category='clustering',
            n_features=len(features), wall_time_sec=elapsed,
            error='<2 features',
            pros='Fast for large datasets',
            cons='Assumes spherical clusters, K estimation needed',
            silhouette=q.get('silhouette', float('nan')),
            sdbw=q.get('sdbw', float('nan')),
            n_clusters=q.get('n_clusters', 0),
        )

    return ApproachResult(
        id=4, approach='K-Means Clustering', category='clustering',
        roc_auc=result.roc_auc, bedroc=result.bedroc,
        ef_1=result.ef_1, ef_5=result.ef_5,
        n_features=result.n_features, wall_time_sec=elapsed,
        best_params={'method': 'kmeans', 'tolerance': 2.0, 'occurrence': 0.5},
        pros='Fast for large datasets',
        cons='Assumes spherical clusters, K estimation needed',
        silhouette=q.get('silhouette', float('nan')),
        sdbw=q.get('sdbw', float('nan')),
        n_clusters=q.get('n_clusters', 0),
    )


def run_05_grid_binning(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """5. Grid-Based Binning."""
    start = time.time()
    features, result, quality = _consensus_and_evaluate(
        refs, actives, decoys, seed,
        tolerance=2.0, occurrence=0.5,
        clustering_method='grid',
        scoring_mode='reference', n_conformers=n_conformers,
        return_quality=True,
    )
    elapsed = time.time() - start
    q = quality or {}

    if result is None:
        return ApproachResult(
            id=5, approach='Grid-Based Binning', category='clustering',
            n_features=len(features), wall_time_sec=elapsed,
            error='<2 features',
            pros='O(N) complexity, fastest clustering',
            cons='Grid alignment artifacts, axis-dependent',
            silhouette=q.get('silhouette', float('nan')),
            sdbw=q.get('sdbw', float('nan')),
            n_clusters=q.get('n_clusters', 0),
        )

    return ApproachResult(
        id=5, approach='Grid-Based Binning', category='clustering',
        roc_auc=result.roc_auc, bedroc=result.bedroc,
        ef_1=result.ef_1, ef_5=result.ef_5,
        n_features=result.n_features, wall_time_sec=elapsed,
        best_params={'method': 'grid', 'tolerance': 2.0, 'occurrence': 0.5},
        pros='O(N) complexity, fastest clustering',
        cons='Grid alignment artifacts, axis-dependent',
        silhouette=q.get('silhouette', float('nan')),
        sdbw=q.get('sdbw', float('nan')),
        n_clusters=q.get('n_clusters', 0),
    )


def run_06_pharm2d(refs, actives, decoys, **_kwargs):
    """6. Pharm2D Fingerprint Scoring (2D only)."""
    start = time.time()

    scorer = Pharm2DScorer(refs)
    active_scores = scorer.score_all(actives)
    decoy_scores = scorer.score_all(decoys)

    y_true = np.concatenate([np.ones(len(active_scores)), np.zeros(len(decoy_scores))])
    y_scores = np.concatenate([active_scores, decoy_scores])

    valid = ~np.isnan(y_scores)
    metrics = calculate_all_metrics(y_true[valid].tolist(), y_scores[valid].tolist())
    elapsed = time.time() - start

    return ApproachResult(
        id=6, approach='Pharm2D Fingerprints', category='scoring',
        roc_auc=metrics.get('roc_auc', 0.5),
        bedroc=metrics.get('bedroc', 0.0),
        ef_1=metrics.get('ef_1', 0.0),
        ef_5=metrics.get('ef_5', 0.0),
        n_features=0,  # 2D — no spatial features
        wall_time_sec=elapsed,
        best_params={'factory': 'Gobbi_Pharm2D', 'similarity': 'Tanimoto'},
        pros='Fastest, no alignment needed, no 3D conformers required',
        cons='No scaffold hopping, sensitive to scaffold bias in dataset',
    )


def run_07_reference_ensemble(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """7. Reference Ensemble (3D rdShapeAlign)."""
    from pharmacophore.rdshape_optimizer import ReferenceEnsembleScorer

    start = time.time()

    scorer = ReferenceEnsembleScorer(
        reference_mols=refs,
        opt_param=0.5,
        n_conformers=n_conformers,
        aggregation='max',
        random_seed=seed,
    )

    active_scores = scorer.score_batch(actives, verbose=False)
    decoy_scores = scorer.score_batch(decoys, verbose=False)

    y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
    y_scores = list(active_scores) + list(decoy_scores)

    metrics = calculate_all_metrics(y_true, y_scores)
    elapsed = time.time() - start

    return ApproachResult(
        id=7, approach='Reference Ensemble (3D)', category='scoring',
        roc_auc=metrics.get('roc_auc', 0.5),
        bedroc=metrics.get('bedroc', 0.0),
        ef_1=metrics.get('ef_1', 0.0),
        ef_5=metrics.get('ef_5', 0.0),
        n_features=len(refs),  # number of reference mols
        wall_time_sec=elapsed,
        best_params={'opt_param': 0.5, 'aggregation': 'max', 'n_conformers': n_conformers},
        pros='Interpretable 3D shape+color matching, best pure-3D approach',
        cons='Needs aligned reference molecules, slower than 2D',
    )


def run_08_hybrid_2d3d(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """8. Hybrid 2D+3D Scoring."""
    start = time.time()

    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)
    config = EvaluationConfig(
        tolerance=2.0, occurrence=0.5,
        scoring_mode='hybrid',
        alpha=0.7,  # 70% 2D, 30% 3D — empirically best on CCR2
        n_conformers=n_conformers,
        opt_param=0.5,
    )
    result = evaluator.evaluate(config)
    elapsed = time.time() - start

    return ApproachResult(
        id=8, approach='Hybrid 2D+3D', category='scoring',
        roc_auc=result.roc_auc, bedroc=result.bedroc,
        ef_1=result.ef_1, ef_5=result.ef_5,
        n_features=result.n_features, wall_time_sec=elapsed,
        best_params={'alpha': 0.7, 'scoring_mode': 'hybrid'},
        pros='Best overall AUC, combines 2D scaffold with 3D shape info',
        cons='Requires both 2D fingerprints and 3D conformers',
    )


def run_09_hungarian_matching(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """9. Hungarian Matching (aligned) as a scoring function.

    Aligns each query molecule to each reference via RDKit shape alignment,
    then computes Hungarian pharmacophore distance on aligned features.
    Takes max similarity across references (ensemble scoring).
    """
    from pharmacophore.hungarian_matching import pharmacophore_similarity_aligned

    start = time.time()

    # Generate conformers for query molecules
    all_mols = actives + decoys
    prepared = []
    for mol in all_mols:
        if mol is None:
            prepared.append(None)
            continue
        m = Chem.AddHs(mol)
        if m.GetNumConformers() == 0:
            AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
        prepared.append(m)

    def score_mol(mol):
        if mol is None or mol.GetNumConformers() == 0:
            return 0.0
        try:
            best_sim = 0.0
            for ref in refs:
                sim = pharmacophore_similarity_aligned(
                    mol, ref, alpha=0.3,
                )
                if sim > best_sim:
                    best_sim = sim
            return best_sim
        except Exception:
            return 0.0

    scores = [score_mol(m) for m in prepared]
    active_scores = scores[:len(actives)]
    decoy_scores = scores[len(actives):]

    y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
    y_scores = active_scores + decoy_scores

    metrics = calculate_all_metrics(y_true, y_scores)
    elapsed = time.time() - start

    return ApproachResult(
        id=9, approach='Hungarian Matching (aligned)', category='scoring',
        roc_auc=metrics.get('roc_auc', 0.5),
        bedroc=metrics.get('bedroc', 0.0),
        ef_1=metrics.get('ef_1', 0.0),
        ef_5=metrics.get('ef_5', 0.0),
        n_features=len(refs), wall_time_sec=elapsed,
        best_params={'alpha': 0.3, 'aggregation': 'max'},
        pros='Aligned feature comparison, optimal assignment, metric-space distance',
        cons='Slow (AlignMol per query-ref pair), needs 3D conformers',
    )


def run_10_optimal_transport(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """10. Optimal Transport (aligned) / Wasserstein Distance scoring.

    Aligns each query molecule to each reference via RDKit shape alignment,
    then computes Wasserstein pharmacophore distance on aligned features.
    Takes max similarity across references.
    """
    from pharmacophore.ot_scoring import wasserstein_similarity_aligned

    start = time.time()

    # Generate conformers for query molecules
    all_mols = actives + decoys
    prepared = []
    for mol in all_mols:
        if mol is None:
            prepared.append(None)
            continue
        m = Chem.AddHs(mol)
        if m.GetNumConformers() == 0:
            AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
        prepared.append(m)

    def score_mol(mol):
        if mol is None or mol.GetNumConformers() == 0:
            return 0.0
        try:
            best_sim = 0.0
            for ref in refs:
                sim = wasserstein_similarity_aligned(
                    mol, ref, blend_alpha=0.3,
                )
                if sim > best_sim:
                    best_sim = sim
            return best_sim
        except Exception:
            return 0.0

    scores = [score_mol(m) for m in prepared]
    active_scores = scores[:len(actives)]
    decoy_scores = scores[len(actives):]

    y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
    y_scores = active_scores + decoy_scores

    metrics = calculate_all_metrics(y_true, y_scores)
    elapsed = time.time() - start

    return ApproachResult(
        id=10, approach='Optimal Transport (aligned)', category='scoring',
        roc_auc=metrics.get('roc_auc', 0.5),
        bedroc=metrics.get('bedroc', 0.0),
        ef_1=metrics.get('ef_1', 0.0),
        ef_5=metrics.get('ef_5', 0.0),
        n_features=len(refs), wall_time_sec=elapsed,
        best_params={'blend_alpha': 0.3, 'ot_alpha': 0.5, 'aggregation': 'max'},
        pros='Aligned feature comparison, theoretically sound OT, handles partial matching',
        cons='Slow (AlignMol per query-ref pair), needs 3D conformers',
    )


def run_11_optuna_gp(refs, actives, decoys, seed=42, n_trials=30, **_kwargs):
    """11. Optuna GP (Bayesian Optimization)."""
    if not HAS_OPTUNA:
        return ApproachResult(
            id=11, approach='Optuna GP (Bayesian)', category='optimization',
            error='optuna not installed',
            pros='Sample efficient, Gaussian Process surrogate',
            cons='O(n^3) scaling with trials, needs optuna',
        )

    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)

    try:
        from optuna.samplers import GPSampler
        sampler = GPSampler(seed=seed)
        label = 'Optuna GP (Bayesian)'
    except Exception:
        sampler = optuna.samplers.TPESampler(seed=seed)
        label = 'Optuna TPE (fallback)'

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
    study.optimize(objective, n_trials=n_trials, show_progress_bar=True)
    elapsed = time.time() - start

    completed = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    if not completed:
        return ApproachResult(
            id=11, approach=label, category='optimization',
            wall_time_sec=elapsed, error='No completed trials',
            pros='Sample efficient, Gaussian Process surrogate',
            cons='O(n^3) scaling with trials',
        )

    best = max(completed, key=lambda t: t.values[0])
    return ApproachResult(
        id=11, approach=label, category='optimization',
        roc_auc=best.values[0], bedroc=best.values[1],
        ef_1=0.0, ef_5=0.0,  # Not tracked per-trial
        n_features=0, wall_time_sec=elapsed,
        best_params=best.params,
        pros='Sample efficient, Gaussian Process surrogate',
        cons='O(n^3) scaling with trials',
    )


def run_12_optuna_nsga2(refs, actives, decoys, seed=42, n_trials=30, **_kwargs):
    """12. Optuna NSGA-II (Multi-objective evolutionary)."""
    if not HAS_OPTUNA:
        return ApproachResult(
            id=12, approach='Optuna NSGA-II', category='optimization',
            error='optuna not installed',
            pros='Multi-objective Pareto front, population-based',
            cons='Needs many trials for convergence',
        )

    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)
    sampler = optuna.samplers.NSGAIISampler(seed=seed)
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
    study.optimize(objective, n_trials=n_trials, show_progress_bar=True)
    elapsed = time.time() - start

    completed = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    if not completed:
        return ApproachResult(
            id=12, approach='Optuna NSGA-II', category='optimization',
            wall_time_sec=elapsed, error='No completed trials',
            pros='Multi-objective Pareto front, population-based',
            cons='Needs many trials for convergence',
        )

    best = max(completed, key=lambda t: t.values[0])
    pareto = optuna.visualization.is_available()  # just check availability

    return ApproachResult(
        id=12, approach='Optuna NSGA-II', category='optimization',
        roc_auc=best.values[0], bedroc=best.values[1],
        ef_1=0.0, ef_5=0.0,
        n_features=0, wall_time_sec=elapsed,
        best_params=best.params,
        pros='Multi-objective Pareto front, population-based',
        cons='Needs many trials for convergence',
    )


def run_13_hypogen(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """13. HypoGen 3-Phase (constructive-subtractive-refinement)."""
    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)

    start = time.time()

    # Phase 1: Generate consensus and enumerate feature subsets
    consensus = PharmacophoreConsensus(tolerance=2.0, occurrence_threshold=0.3)
    features = consensus.generate_consensus(refs)

    subsets = []
    for n_feat in range(3, min(8, len(features) + 1)):
        combos = list(itertools.combinations(features, n_feat))
        subsets.extend([list(c) for c in combos[:20]])
        if len(subsets) >= 100:
            break

    print(f"    Phase 1: Evaluating {len(subsets)} hypotheses...")

    def eval_subset(subset):
        result = evaluator.evaluate_feature_subset(subset, 0.5, 0.5, min(n_conformers, 5))
        cost = 0.5 * (1 - result.roc_auc) + 0.3 * (1 - result.bedroc) + 0.2 * (len(subset) / 8)
        return {
            'features': subset,
            'roc_auc': result.roc_auc,
            'bedroc': result.bedroc,
            'ef_1': result.ef_1,
            'ef_5': result.ef_5,
            'n_features': len(subset),
            'cost': cost,
        }

    # Use threading backend: eval_subset captures `evaluator` which contains
    # ShapeInput objects (_ref_shapes) that can't be pickled by loky.
    # Threading avoids pickling and is GIL-safe since RDKit scoring
    # releases the GIL in C++ code.
    hypotheses = Parallel(n_jobs=N_JOBS, backend='threading')(
        delayed(eval_subset)(s) for s in subsets
    )

    # Phase 2: Filter survivors
    survivors = [h for h in hypotheses if h['roc_auc'] >= 0.50]
    print(f"    Phase 2: {len(survivors)}/{len(hypotheses)} hypotheses survived")

    # Phase 3: Select best by cost
    best = min(survivors, key=lambda h: h['cost']) if survivors else min(hypotheses, key=lambda h: h['cost'])
    elapsed = time.time() - start

    return ApproachResult(
        id=13, approach='HypoGen 3-Phase', category='optimization',
        roc_auc=best['roc_auc'], bedroc=best['bedroc'],
        ef_1=best.get('ef_1', 0.0), ef_5=best.get('ef_5', 0.0),
        n_features=best['n_features'], wall_time_sec=elapsed,
        best_params={'n_features': best['n_features'], 'cost': best['cost'],
                     'n_hypotheses_evaluated': len(hypotheses)},
        pros='Interpretable feature selection, fast exhaustive enumeration',
        cons='SA refinement is stochastic, limited to enumerable subset space',
    )


def run_14_multifidelity_bo(refs, actives, decoys, seed=42, n_trials=30, **_kwargs):
    """14. Multi-Fidelity Bayesian Optimization (combinatorial optimizer).

    Hardened against segfaults: the CombinatorialPharmacophoreOptimizer
    creates one UnifiedEvaluator per discrete combo (128 total), each
    with PrepareConformer objects. Accumulated resource leaks cause
    RDKit C++ memory corruption around combo ~50. Mitigation:
    - Limit parallelism (max 4 jobs)
    - gc.collect() is called internally between combos
    - Catch SystemError/segfault-survivors and return partial results
    """
    import gc
    from pharmacophore.combinatorial_optimizer import CombinatorialPharmacophoreOptimizer

    start = time.time()

    optimizer = CombinatorialPharmacophoreOptimizer(
        reference_mols=refs,
        actives=actives,
        decoys=decoys,
        random_state=seed,
    )

    # n_trials here means "per discrete combo" — cap at 3 to keep runtime sane
    # The optimizer has 128 discrete combos; multi-fidelity filtering prunes most
    mf_trials = min(n_trials, 3)

    try:
        result = optimizer.optimize(
            n_trials=mf_trials, multi_fidelity=True, verbose=True,
        )
        gc.collect()
    except (SystemError, OSError, MemoryError) as e:
        elapsed = time.time() - start
        print(f"    MF-BO crashed ({type(e).__name__}: {e}), returning partial results")
        # Try to recover partial results from optimizer internals
        if hasattr(optimizer, '_best_result') and optimizer._best_result is not None:
            result = optimizer._best_result
        else:
            return ApproachResult(
                id=14, approach='Multi-Fidelity BO', category='optimization',
                wall_time_sec=elapsed,
                error=f'Crashed: {type(e).__name__}: {e}',
                pros='Explores cheap 2D first, then expensive 3D — efficient',
                cons='Resource leaks cause crashes on large combo spaces',
            )

    elapsed = time.time() - start
    best_metrics = result.best_metrics

    return ApproachResult(
        id=14, approach='Multi-Fidelity BO', category='optimization',
        roc_auc=best_metrics.get('roc_auc', 0.5),
        bedroc=best_metrics.get('bedroc', 0.0),
        ef_1=best_metrics.get('ef_1', 0.0),
        ef_5=best_metrics.get('ef_5', 0.0),
        n_features=0, wall_time_sec=elapsed,
        best_params=result.best_config if isinstance(result.best_config, dict) else {},
        pros='Explores cheap 2D first, then expensive 3D — efficient',
        cons='Complex setup, many degrees of freedom',
    )


def run_15_strategy_selector(refs, actives, decoys, seed=42, **_kwargs):
    """15. Strategy Selector (tournament)."""
    from pharmacophore.auto_strategy import StrategySelector

    start = time.time()

    selector = StrategySelector(
        refs, actives, decoys, random_state=seed, verbose=True
    )
    best = selector.select_best()
    elapsed = time.time() - start

    return ApproachResult(
        id=15, approach=f'Strategy Selector ({best.strategy_name})', category='meta',
        roc_auc=best.eval_result.roc_auc,
        bedroc=best.eval_result.bedroc,
        ef_1=best.eval_result.ef_1,
        ef_5=best.eval_result.ef_5,
        n_features=best.eval_result.n_features, wall_time_sec=elapsed,
        best_params={
            'strategy': best.strategy_name,
            'sdbw': float(best.sdbw) if best.sdbw != float('inf') else None,
        },
        pros='Automatic strategy selection, no manual tuning',
        cons='Limited to predefined config grid (9 configs)',
    )


def run_16_ensemble_rrf(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """16. Ensemble RRF Scoring (cascading multi-level with rank fusion)."""
    evaluator = UnifiedEvaluator(refs, actives, decoys, seed)

    start = time.time()

    config = EvaluationConfig(
        tolerance=2.0, occurrence=0.5, shape_weight=0.5,
        opt_param=0.5, n_conformers=n_conformers,
    )
    result = evaluator.evaluate_ensemble(config, quick_mode=False)
    elapsed = time.time() - start

    return ApproachResult(
        id=16, approach='Ensemble RRF Scoring', category='meta',
        roc_auc=result.roc_auc, bedroc=result.bedroc,
        ef_1=result.ef_1, ef_5=result.ef_5,
        n_features=result.n_features, wall_time_sec=elapsed,
        best_params={'rrf_constant': 60, 'scoring': 'reference_3d+ot_fusion'},
        pros='Multi-signal fusion, robust via rank combination',
        cons='Slower (runs both 3D + OT), overkill for simple cases',
    )


def run_17_point_cloud_icp(refs, actives, decoys, seed=42, n_conformers=10, **_kwargs):
    """17. Point Cloud ICP alignment scoring.

    Aligns each query molecule to each reference via RDKit shape alignment,
    then runs colored ICP (combining spatial + type distance) on the aligned
    pharmacophore features. Takes max similarity across references.

    Based on Zhou, Griffith & Gaeta (BMC Bioinformatics 2014).
    """
    from pharmacophore.point_cloud_alignment import point_cloud_similarity_aligned

    start = time.time()

    # Generate conformers for query molecules
    all_mols = actives + decoys
    prepared = []
    for mol in all_mols:
        if mol is None:
            prepared.append(None)
            continue
        m = Chem.AddHs(mol)
        if m.GetNumConformers() == 0:
            AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
        prepared.append(m)

    def score_mol(mol):
        if mol is None or mol.GetNumConformers() == 0:
            return 0.0
        try:
            best_sim = 0.0
            for ref in refs:
                sim = point_cloud_similarity_aligned(
                    mol, ref, alpha=0.3, lambda_color=0.5, sigma=2.0,
                )
                if sim > best_sim:
                    best_sim = sim
            return best_sim
        except Exception:
            return 0.0

    scores = [score_mol(m) for m in prepared]
    active_scores = scores[:len(actives)]
    decoy_scores = scores[len(actives):]

    y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
    y_scores = active_scores + decoy_scores

    metrics = calculate_all_metrics(y_true, y_scores)
    elapsed = time.time() - start

    return ApproachResult(
        id=17, approach='Point Cloud ICP', category='scoring',
        roc_auc=metrics.get('roc_auc', 0.5),
        bedroc=metrics.get('bedroc', 0.0),
        ef_1=metrics.get('ef_1', 0.0),
        ef_5=metrics.get('ef_5', 0.0),
        n_features=len(refs), wall_time_sec=elapsed,
        best_params={'alpha': 0.3, 'lambda_color': 0.5, 'sigma': 2.0},
        pros='Combines spatial + type alignment, handles partial overlap, literature-validated',
        cons='Slow (AlignMol + ICP per pair), sensitive to lambda_color',
    )


# ---------------------------------------------------------------------------
# Analysis & visualization
# ---------------------------------------------------------------------------

def generate_rankings(results: List[ApproachResult]) -> Dict[str, List]:
    """Generate multiple rankings from results."""
    valid = [r for r in results if r.error is None]

    rankings = {}

    # By AUC
    by_auc = sorted(valid, key=lambda r: r.roc_auc, reverse=True)
    rankings['by_auc'] = [{'rank': i+1, 'approach': r.approach, 'roc_auc': r.roc_auc}
                          for i, r in enumerate(by_auc)]

    # By BEDROC
    by_bedroc = sorted(valid, key=lambda r: r.bedroc, reverse=True)
    rankings['by_bedroc'] = [{'rank': i+1, 'approach': r.approach, 'bedroc': r.bedroc}
                             for i, r in enumerate(by_bedroc)]

    # By EF@1%
    by_ef1 = sorted(valid, key=lambda r: r.ef_1, reverse=True)
    rankings['by_ef1'] = [{'rank': i+1, 'approach': r.approach, 'ef_1': r.ef_1}
                          for i, r in enumerate(by_ef1)]

    # By speed (fastest first)
    by_speed = sorted(valid, key=lambda r: r.wall_time_sec)
    rankings['by_speed'] = [{'rank': i+1, 'approach': r.approach, 'wall_time_sec': r.wall_time_sec}
                            for i, r in enumerate(by_speed)]

    # By efficiency (AUC / time)
    for r in valid:
        r._efficiency = r.roc_auc / max(r.wall_time_sec, 0.01)
    by_efficiency = sorted(valid, key=lambda r: r._efficiency, reverse=True)
    rankings['by_efficiency'] = [
        {'rank': i+1, 'approach': r.approach,
         'auc_per_sec': r._efficiency}
        for i, r in enumerate(by_efficiency)
    ]

    return rankings


def generate_plots(results: List[ApproachResult], timestamp: str):
    """Generate visualization plots."""
    if not HAS_MATPLOTLIB:
        print("  matplotlib not available, skipping plots")
        return

    valid = [r for r in results if r.error is None]
    if not valid:
        return

    names = [r.approach for r in valid]
    aucs = [r.roc_auc for r in valid]
    bedrocs = [r.bedroc for r in valid]
    times = [r.wall_time_sec for r in valid]
    categories = [r.category for r in valid]

    category_colors = {
        'clustering': '#4285f4',
        'scoring': '#ea4335',
        'optimization': '#fbbc04',
        'meta': '#34a853',
    }
    colors = [category_colors.get(c, '#999999') for c in categories]

    # --- Plot 1: AUC + BEDROC bar chart ---
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # Sort by AUC for the bar chart
    sorted_idx = np.argsort(aucs)
    sorted_names = [names[i] for i in sorted_idx]
    sorted_aucs = [aucs[i] for i in sorted_idx]
    sorted_bedrocs = [bedrocs[i] for i in sorted_idx]
    sorted_colors = [colors[i] for i in sorted_idx]

    axes[0].barh(range(len(sorted_names)), sorted_aucs, color=sorted_colors, alpha=0.8)
    axes[0].set_yticks(range(len(sorted_names)))
    axes[0].set_yticklabels(sorted_names, fontsize=9)
    axes[0].set_xlabel('ROC-AUC')
    axes[0].set_title('ROC-AUC by Approach')
    axes[0].axvline(x=0.5, color='gray', linestyle='--', alpha=0.5, label='Random')
    axes[0].set_xlim(0, 1.05)

    # Sort by BEDROC
    sorted_idx_b = np.argsort(bedrocs)
    sorted_names_b = [names[i] for i in sorted_idx_b]
    sorted_bedrocs_b = [bedrocs[i] for i in sorted_idx_b]
    sorted_colors_b = [colors[i] for i in sorted_idx_b]

    axes[1].barh(range(len(sorted_names_b)), sorted_bedrocs_b, color=sorted_colors_b, alpha=0.8)
    axes[1].set_yticks(range(len(sorted_names_b)))
    axes[1].set_yticklabels(sorted_names_b, fontsize=9)
    axes[1].set_xlabel('BEDROC (alpha=20)')
    axes[1].set_title('BEDROC by Approach')
    axes[1].set_xlim(0, 1.05)

    # Legend
    from matplotlib.patches import Patch
    legend_handles = [Patch(color=c, label=cat) for cat, c in category_colors.items()]
    axes[0].legend(handles=legend_handles, loc='lower right', fontsize=8)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / f'benchmark_bars_{timestamp}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {FIGURES_DIR / f'benchmark_bars_{timestamp}.png'}")

    # --- Plot 2: AUC vs Speed (Pareto front) ---
    fig, ax = plt.subplots(figsize=(12, 8))

    for i, r in enumerate(valid):
        ax.scatter(r.wall_time_sec, r.roc_auc, c=category_colors.get(r.category, '#999'),
                   s=100, zorder=5, alpha=0.8)
        ax.annotate(r.approach, (r.wall_time_sec, r.roc_auc),
                    textcoords='offset points', xytext=(5, 5), fontsize=7)

    # Draw Pareto front
    pareto_points = []
    sorted_by_time = sorted(valid, key=lambda r: r.wall_time_sec)
    best_auc_so_far = -1
    for r in sorted_by_time:
        if r.roc_auc > best_auc_so_far:
            pareto_points.append(r)
            best_auc_so_far = r.roc_auc

    if len(pareto_points) > 1:
        pareto_times = [r.wall_time_sec for r in pareto_points]
        pareto_aucs = [r.roc_auc for r in pareto_points]
        ax.plot(pareto_times, pareto_aucs, 'k--', alpha=0.4, linewidth=1, label='Pareto front')

    ax.set_xlabel('Wall Time (seconds)')
    ax.set_ylabel('ROC-AUC')
    ax.set_title('AUC vs Speed Trade-off')
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.3)
    ax.legend(handles=legend_handles, loc='lower right', fontsize=8)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / f'benchmark_pareto_{timestamp}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {FIGURES_DIR / f'benchmark_pareto_{timestamp}.png'}")

    # --- Plot 3: Heatmap of all metrics ---
    fig, ax = plt.subplots(figsize=(14, 8))

    metric_names = ['ROC-AUC', 'BEDROC', 'EF@1%', 'EF@5%']
    data = np.array([
        [r.roc_auc for r in valid],
        [r.bedroc for r in valid],
        [r.ef_1 for r in valid],
        [r.ef_5 for r in valid],
    ])

    # Normalize each metric to [0, 1] for heatmap
    data_norm = np.zeros_like(data)
    for i in range(data.shape[0]):
        row_min, row_max = data[i].min(), data[i].max()
        if row_max > row_min:
            data_norm[i] = (data[i] - row_min) / (row_max - row_min)
        else:
            data_norm[i] = 0.5

    im = ax.imshow(data_norm, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)

    ax.set_xticks(range(len(valid)))
    ax.set_xticklabels([r.approach for r in valid], rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(metric_names)))
    ax.set_yticklabels(metric_names)

    # Annotate with actual values
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            fmt = f'{val:.3f}' if val < 10 else f'{val:.1f}'
            ax.text(j, i, fmt, ha='center', va='center', fontsize=7,
                    color='black' if data_norm[i, j] > 0.3 else 'white')

    plt.colorbar(im, ax=ax, label='Normalized Score (0=worst, 1=best)')
    ax.set_title('Metric Heatmap: All Approaches')

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / f'benchmark_heatmap_{timestamp}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {FIGURES_DIR / f'benchmark_heatmap_{timestamp}.png'}")


def generate_report(results: List[ApproachResult], rankings: Dict,
                    timestamp: str, total_time: float):
    """Generate markdown comparison report."""
    report_path = RESULTS_DIR / f'APPROACH_COMPARISON_REPORT_{timestamp}.md'

    valid = [r for r in results if r.error is None]
    failed = [r for r in results if r.error is not None]

    with open(report_path, 'w') as f:
        f.write('# Comprehensive Pharmacophore Approach Comparison Report\n\n')
        f.write(f'**Date**: {datetime.now().strftime("%Y-%m-%d %H:%M")}\n')
        f.write(f'**Dataset**: CCR2 (5 refs, ~74 actives, ~499 decoys)\n')
        f.write(f'**Total runtime**: {total_time:.1f}s ({total_time/60:.1f} min)\n')
        f.write(f'**Approaches evaluated**: {len(valid)}/{len(results)}\n\n')

        # Summary table
        f.write('## Summary Table\n\n')
        f.write(f'| {"#":>2} | {"Approach":<30} | {"Cat":>8} | {"AUC":>7} | '
                f'{"BEDROC":>7} | {"EF@1%":>7} | {"EF@5%":>7} | '
                f'{"Feats":>5} | {"Time(s)":>8} |\n')
        f.write(f'|{"---":>4}|{"-"*32}|{"-"*10}|{"-"*9}|{"-"*9}|'
                f'{"-"*9}|{"-"*9}|{"-"*7}|{"-"*10}|\n')

        for r in sorted(results, key=lambda x: x.id):
            auc_str = f'{r.roc_auc:.4f}' if r.error is None else 'FAIL'
            bed_str = f'{r.bedroc:.4f}' if r.error is None else '-'
            ef1_str = f'{r.ef_1:.2f}' if r.error is None else '-'
            ef5_str = f'{r.ef_5:.2f}' if r.error is None else '-'
            f.write(f'| {r.id:>2} | {r.approach:<30} | {r.category:>8} | '
                    f'{auc_str:>7} | {bed_str:>7} | {ef1_str:>7} | '
                    f'{ef5_str:>7} | {r.n_features:>5} | '
                    f'{r.wall_time_sec:>8.1f} |\n')

        f.write('\n> **Note**: Approaches 1-5 share identical VS metrics '
                '(AUC, BEDROC, EF) because `scoring_mode=\'reference\'` '
                'scores queries against reference molecules directly — '
                'consensus features are not used for scoring. The '
                'Clustering Quality section below differentiates them.\n')

        # Clustering quality comparison
        clustering_results = [r for r in results
                              if r.category == 'clustering' and r.error is None]
        if clustering_results:
            f.write('\n## Clustering Quality Comparison\n\n')
            f.write('These metrics measure how well each clustering method '
                    'partitions pharmacophore features, independent of VS scoring.\n\n')
            f.write(f'| {"#":>2} | {"Approach":<30} | {"Silhouette":>10} | '
                    f'{"S_Dbw":>8} | {"Clusters":>8} | {"Features":>8} |\n')
            f.write(f'|{"---":>4}|{"-"*32}|{"-"*12}|{"-"*10}|{"-"*10}|{"-"*10}|\n')

            for r in sorted(clustering_results, key=lambda x: x.id):
                sil_str = f'{r.silhouette:.4f}' if not np.isnan(r.silhouette) else '-'
                sdbw_str = f'{r.sdbw:.4f}' if not np.isnan(r.sdbw) else '-'
                f.write(f'| {r.id:>2} | {r.approach:<30} | '
                        f'{sil_str:>10} | {sdbw_str:>8} | '
                        f'{r.n_clusters:>8} | {r.n_features:>8} |\n')

            f.write('\n*Silhouette: higher is better ([-1, 1]). '
                    'S_Dbw: lower is better (compactness + separation).*\n')

        # Rankings
        f.write('\n## Rankings\n\n')

        f.write('### By ROC-AUC (Overall Discrimination)\n\n')
        for item in rankings['by_auc']:
            f.write(f'{item["rank"]}. **{item["approach"]}** — {item["roc_auc"]:.4f}\n')

        f.write('\n### By BEDROC (Early Enrichment)\n\n')
        for item in rankings['by_bedroc']:
            f.write(f'{item["rank"]}. **{item["approach"]}** — {item["bedroc"]:.4f}\n')

        f.write('\n### By Speed (Fastest First)\n\n')
        for item in rankings['by_speed']:
            f.write(f'{item["rank"]}. **{item["approach"]}** — {item["wall_time_sec"]:.1f}s\n')

        f.write('\n### By Efficiency (AUC / Time)\n\n')
        for item in rankings['by_efficiency']:
            f.write(f'{item["rank"]}. **{item["approach"]}** — {item["auc_per_sec"]:.4f} AUC/s\n')

        # Pros & Cons
        f.write('\n## Approach Details (Pros & Cons)\n\n')
        f.write(f'| {"Approach":<30} | {"Pros":<55} | {"Cons":<55} |\n')
        f.write(f'|{"-"*32}|{"-"*57}|{"-"*57}|\n')
        for r in sorted(results, key=lambda x: x.id):
            f.write(f'| {r.approach:<30} | {r.pros:<55} | {r.cons:<55} |\n')

        # Category analysis
        f.write('\n## Analysis by Category\n\n')
        for cat in ['clustering', 'scoring', 'optimization', 'meta']:
            cat_results = [r for r in valid if r.category == cat]
            if not cat_results:
                continue
            f.write(f'### {cat.title()} Methods\n\n')
            best = max(cat_results, key=lambda r: r.roc_auc)
            fastest = min(cat_results, key=lambda r: r.wall_time_sec)
            f.write(f'- **Best AUC**: {best.approach} ({best.roc_auc:.4f})\n')
            f.write(f'- **Fastest**: {fastest.approach} ({fastest.wall_time_sec:.1f}s)\n')
            avg_auc = np.mean([r.roc_auc for r in cat_results])
            f.write(f'- **Average AUC**: {avg_auc:.4f}\n\n')

        # Failed approaches
        if failed:
            f.write('\n## Failed Approaches\n\n')
            for r in failed:
                f.write(f'- **{r.approach}**: {r.error}\n')

        # Recommendations
        f.write('\n## Recommendations\n\n')
        if valid:
            top = max(valid, key=lambda r: r.roc_auc)
            fast = min(valid, key=lambda r: r.wall_time_sec)
            f.write(f'1. **Best overall**: {top.approach} (AUC={top.roc_auc:.4f})\n')
            f.write(f'2. **Fastest screening**: {fast.approach} '
                    f'({fast.wall_time_sec:.1f}s, AUC={fast.roc_auc:.4f})\n')

            # Best per category
            for cat in ['clustering', 'scoring', 'optimization', 'meta']:
                cat_valid = [r for r in valid if r.category == cat]
                if cat_valid:
                    best_cat = max(cat_valid, key=lambda r: r.roc_auc)
                    f.write(f'3. **Best {cat}**: {best_cat.approach} '
                            f'(AUC={best_cat.roc_auc:.4f})\n')

    print(f"  Report saved: {report_path}")
    return report_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive benchmark of all 16 pharmacophore approaches")
    parser.add_argument('--n-conformers', type=int, default=10,
                        help='Conformers per molecule (default: 10)')
    parser.add_argument('--n-trials', type=int, default=30,
                        help='Optuna trials per optimizer approach (default: 30)')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--skip-slow', action='store_true',
                        help='Skip optimizer approaches (11-14)')
    parser.add_argument('--only', type=str, default=None,
                        help='Comma-separated approach IDs to run (e.g. "1,6,8")')
    args = parser.parse_args()

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print(f"\n{'='*70}")
    print(f"COMPREHENSIVE BENCHMARK: All Pharmacophore Approaches")
    print(f"{'='*70}")
    print(f"  Dataset: CCR2 (5 refs, ~74 actives, ~499 decoys)")
    print(f"  Conformers: {args.n_conformers} per molecule (ETKDGv3)")
    print(f"  Optuna trials: {args.n_trials} per optimizer")
    print(f"  Seed: {args.seed}")
    print(f"  CPU cores: {os.cpu_count()}")
    print(f"  Skip slow: {args.skip_slow}")
    print(f"  Timestamp: {timestamp}")
    print(f"{'='*70}\n")

    total_start = time.time()

    # --- Load data ---
    print("Loading CCR2 dataset...")
    refs, actives, decoys = load_ccr2_dataset(
        n_conformers=args.n_conformers, verbose=True)

    # --- Define all 16 approach runners ---
    all_approaches = {
        1:  ('Baseline Agglomerative', run_01_baseline_agglomerative),
        2:  ('Spatial-Scaled', run_02_spatial_scaled),
        3:  ('Ensemble Consensus', run_03_ensemble_consensus),
        4:  ('K-Means Clustering', run_04_kmeans),
        5:  ('Grid-Based Binning', run_05_grid_binning),
        6:  ('Pharm2D Fingerprints', run_06_pharm2d),
        7:  ('Reference Ensemble (3D)', run_07_reference_ensemble),
        8:  ('Hybrid 2D+3D', run_08_hybrid_2d3d),
        9:  ('Hungarian Matching', run_09_hungarian_matching),
        10: ('Optimal Transport', run_10_optimal_transport),
        11: ('Optuna GP (Bayesian)', run_11_optuna_gp),
        12: ('Optuna NSGA-II', run_12_optuna_nsga2),
        13: ('HypoGen 3-Phase', run_13_hypogen),
        14: ('Multi-Fidelity BO', run_14_multifidelity_bo),
        15: ('Strategy Selector', run_15_strategy_selector),
        16: ('Ensemble RRF Scoring', run_16_ensemble_rrf),
        17: ('Point Cloud ICP', run_17_point_cloud_icp),
    }

    # Filter approaches
    if args.only:
        selected_ids = [int(x.strip()) for x in args.only.split(',')]
    elif args.skip_slow:
        selected_ids = [i for i in range(1, 17) if i not in (11, 12, 14)]
    else:
        selected_ids = list(range(1, 18))

    # --- Run all selected approaches ---
    results = []
    for approach_id in selected_ids:
        name, runner = all_approaches[approach_id]
        print(f"\n{'='*70}")
        print(f"APPROACH {approach_id}: {name}")
        print(f"{'='*70}")

        try:
            result = runner(
                refs, actives, decoys,
                seed=args.seed,
                n_conformers=args.n_conformers,
                n_trials=args.n_trials,
            )
            results.append(result)

            if result.error:
                print(f"  -> ERROR: {result.error}")
            else:
                print(f"  -> AUC: {result.roc_auc:.4f}, "
                      f"BEDROC: {result.bedroc:.4f}, "
                      f"EF@1%: {result.ef_1:.2f}, "
                      f"Features: {result.n_features}, "
                      f"Time: {result.wall_time_sec:.1f}s")

        except Exception as e:
            print(f"  -> EXCEPTION: {e}")
            traceback.print_exc()
            results.append(ApproachResult(
                id=approach_id, approach=name,
                category='unknown', error=str(e),
            ))

    total_elapsed = time.time() - total_start

    # --- Print summary ---
    print(f"\n{'='*70}")
    print(f"BENCHMARK RESULTS — {len(results)} approaches")
    print(f"{'='*70}\n")

    print(f"{'#':>2} {'Approach':<32} {'Cat':>10} {'AUC':>7} "
          f"{'BEDROC':>8} {'EF@1%':>7} {'EF@5%':>7} {'Feats':>5} {'Time(s)':>8}")
    print(f"{'-'*100}")

    for r in sorted(results, key=lambda x: x.id):
        if r.error:
            print(f"{r.id:>2} {r.approach:<32} {r.category:>10}   FAIL   "
                  f"{'':>8} {'':>7} {'':>7} {'':>5} {r.wall_time_sec:>8.1f}")
        else:
            print(f"{r.id:>2} {r.approach:<32} {r.category:>10} {r.roc_auc:>7.4f} "
                  f"{r.bedroc:>8.4f} {r.ef_1:>7.2f} {r.ef_5:>7.2f} "
                  f"{r.n_features:>5} {r.wall_time_sec:>8.1f}")

    print(f"{'-'*100}")

    # Winners
    valid = [r for r in results if r.error is None]
    if valid:
        best_auc = max(valid, key=lambda r: r.roc_auc)
        best_bedroc = max(valid, key=lambda r: r.bedroc)
        fastest = min(valid, key=lambda r: r.wall_time_sec)

        print(f"\n  Best AUC:    {best_auc.approach} ({best_auc.roc_auc:.4f})")
        print(f"  Best BEDROC: {best_bedroc.approach} ({best_bedroc.bedroc:.4f})")
        print(f"  Fastest:     {fastest.approach} ({fastest.wall_time_sec:.1f}s)")

    print(f"  Total time:  {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")

    # --- Generate rankings ---
    rankings = generate_rankings(results)

    # --- Save results ---
    results_dicts = [r.to_dict() for r in results]

    output = {
        'timestamp': timestamp,
        'n_conformers': args.n_conformers,
        'n_trials': args.n_trials,
        'seed': args.seed,
        'total_time_sec': total_elapsed,
        'results': results_dicts,
        'rankings': rankings,
    }

    json_path = RESULTS_DIR / f'benchmark_all_{timestamp}.json'
    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n  JSON saved: {json_path}")

    # CSV summary
    csv_path = RESULTS_DIR / f'benchmark_all_{timestamp}.csv'
    rows = []
    for r in results:
        row = {
            'id': r.id,
            'approach': r.approach,
            'category': r.category,
            'roc_auc': r.roc_auc if r.error is None else None,
            'bedroc': r.bedroc if r.error is None else None,
            'ef_1': r.ef_1 if r.error is None else None,
            'ef_5': r.ef_5 if r.error is None else None,
            'n_features': r.n_features,
            'wall_time_sec': r.wall_time_sec,
            'error': r.error,
        }
        if r.category == 'clustering':
            row['silhouette'] = r.silhouette if not np.isnan(r.silhouette) else None
            row['sdbw'] = r.sdbw if not np.isnan(r.sdbw) else None
            row['n_clusters'] = r.n_clusters
        rows.append(row)
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    print(f"  CSV saved: {csv_path}")

    # --- Generate plots ---
    print("\nGenerating visualizations...")
    generate_plots(results, timestamp)

    # --- Generate report ---
    print("\nGenerating comparison report...")
    generate_report(results, rankings, timestamp, total_elapsed)

    print(f"\n{'='*70}")
    print(f"BENCHMARK COMPLETE")
    print(f"{'='*70}\n")


if __name__ == '__main__':
    main()

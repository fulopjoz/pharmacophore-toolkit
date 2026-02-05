"""Combinatorial pharmacophore optimizer with multi-fidelity scoring.

This module implements a two-tier search over the 9 degrees of freedom
for pharmacophore-based virtual screening:

Continuous (Bayesian optimization via Optuna):
    - tolerance (0.5-4.0 Å)
    - occurrence (0.1-1.0)
    - opt_param (0.0-1.0)
    - alpha (0.3-0.9, hybrid mode only)

Fixed (set once, cached for speed):
    - n_conformers (default 10)

Discrete (enumerated):
    - reference_subset (C(N,k) subsets)
    - aggregation ('max', 'mean')
    - minimize (True, False)
    - scoring_mode ('reference', 'hybrid')

Architecture:
    1. Enumerate discrete combinations
    2. For each, run Optuna BO over continuous dimensions
    3. Multi-fidelity pruning: pharm2d screen → reference → full hybrid

Example:
    >>> from pharmacophore.combinatorial_optimizer import CombinatorialPharmacophoreOptimizer
    >>> optimizer = CombinatorialPharmacophoreOptimizer(ref_mols, actives, decoys)
    >>> result = optimizer.optimize(n_trials=30)
    >>> print(f"Best AUC: {result.best_metrics['roc_auc']:.4f}")
"""

import itertools
import logging
import time
import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .evaluation import EvaluationConfig, EvaluationResult, UnifiedEvaluator

logger = logging.getLogger(__name__)

try:
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    HAS_OPTUNA = True
except ImportError:
    HAS_OPTUNA = False


@dataclass
class CombinatorialResult:
    """Results from combinatorial pharmacophore optimization.

    Attributes:
        best_config: All 9 degrees of freedom for the best configuration.
        best_metrics: AUC, BEDROC, EF@1% for the best configuration.
        pareto_front: List of (config, AUC, wall_time) triples on the
            Pareto front of quality vs speed.
        feature_importance: Relative importance of each DoF (from fANOVA
            if available, else from variance analysis).
        top_k_configs: Top 10 configurations for inspection.
        total_evaluations: Total number of configurations evaluated.
        wall_time_sec: Total wall-clock time.
    """
    best_config: Dict[str, Any]
    best_metrics: Dict[str, float]
    pareto_front: List[Tuple[Dict, float, float]] = field(default_factory=list)
    feature_importance: Dict[str, float] = field(default_factory=dict)
    top_k_configs: List[Dict[str, Any]] = field(default_factory=list)
    total_evaluations: int = 0
    wall_time_sec: float = 0.0


def _generate_reference_subsets(
    n_refs: int,
    min_k: int = 3
) -> List[Tuple[int, ...]]:
    """Generate all reference molecule subsets of size >= min_k.

    Args:
        n_refs: Total number of reference molecules.
        min_k: Minimum subset size (default: 3).

    Returns:
        List of tuples of reference indices.
    """
    subsets = []
    for k in range(min_k, n_refs + 1):
        subsets.extend(itertools.combinations(range(n_refs), k))
    return subsets


class CombinatorialPharmacophoreOptimizer:
    """Two-tier combinatorial + Bayesian optimizer for pharmacophore screening.

    Enumerates discrete parameter combinations and runs Optuna Bayesian
    optimization for continuous parameters within each combination.

    Supports multi-fidelity pruning: screens with fast pharm2d first,
    promotes top fraction to reference scoring, and top fraction of those
    to full hybrid scoring.

    Args:
        reference_mols: Aligned reference molecules.
        actives: Active compounds.
        decoys: Decoy compounds.
        random_state: Random seed for reproducibility.
        n_jobs: Parallel jobs for scoring.

    Example:
        >>> optimizer = CombinatorialPharmacophoreOptimizer(refs, actives, decoys)
        >>> result = optimizer.optimize(n_trials=30)
        >>> print(result.best_config)
    """

    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol],
        random_state: int = 42,
        n_jobs: int = 1
    ):
        if not HAS_OPTUNA:
            raise ImportError(
                "optuna required for CombinatorialPharmacophoreOptimizer. "
                "Install with: pip install optuna"
            )

        self.reference_mols = reference_mols
        self.actives = actives
        self.decoys = decoys
        self.random_state = random_state
        self.n_jobs = n_jobs

        # Full evaluator with all references
        self.evaluator = UnifiedEvaluator(
            reference_mols, actives, decoys,
            random_state=random_state, n_jobs=n_jobs
        )

        self._all_results: List[Dict[str, Any]] = []
        self._evaluator_cache: Dict[Tuple[int, ...], UnifiedEvaluator] = {}

    def _create_subset_evaluator(
        self,
        ref_indices: Tuple[int, ...]
    ) -> UnifiedEvaluator:
        """Create evaluator for a reference subset (cached per subset)."""
        if ref_indices in self._evaluator_cache:
            return self._evaluator_cache[ref_indices]
        subset_mols = [self.reference_mols[i] for i in ref_indices]
        evaluator = UnifiedEvaluator(
            subset_mols, self.actives, self.decoys,
            random_state=self.random_state, n_jobs=self.n_jobs
        )
        self._evaluator_cache[ref_indices] = evaluator
        return evaluator

    def _evaluate_discrete_combo(
        self,
        ref_indices: Tuple[int, ...],
        aggregation: str,
        minimize: bool,
        scoring_mode: str,
        n_trials: int = 20,
        n_conformers: int = 10,
        multi_fidelity: bool = True,
        min_auc_threshold: float = 0.6
    ) -> Dict[str, Any]:
        """Evaluate one discrete combination with BO over continuous params.

        Args:
            ref_indices: Reference molecule indices.
            aggregation: 'max' or 'mean'.
            minimize: Whether to MMFF-minimize conformers.
            scoring_mode: 'reference' or 'hybrid'.
            n_trials: Optuna trials for continuous optimization.
            n_conformers: Fixed conformer count (not optimized to avoid
                regenerating 3D coords every trial).
            multi_fidelity: Use pharm2d pre-screening.
            min_auc_threshold: Minimum pharm2d AUC to proceed.

        Returns:
            Dict with best config and metrics for this combo.
        """
        evaluator = self._create_subset_evaluator(ref_indices)

        # Multi-fidelity pre-screen: quick pharm2d check
        if multi_fidelity:
            quick_config = EvaluationConfig(
                scoring_mode='pharm2d',
                aggregation=aggregation,
            )
            quick_result = evaluator.evaluate(quick_config)
            if quick_result.roc_auc < min_auc_threshold:
                return {
                    'ref_indices': ref_indices,
                    'aggregation': aggregation,
                    'minimize': minimize,
                    'scoring_mode': scoring_mode,
                    'best_auc': quick_result.roc_auc,
                    'pruned': True,
                    'best_config': None,
                    'best_result': quick_result
                }

        # Optuna study for continuous parameters
        # n_conformers is fixed to avoid costly conformer regeneration
        # per trial — conformers are cached by (n_conformers, minimize).
        # For reference/hybrid modes, tolerance/occurrence don't affect
        # scoring (they control consensus generation only), so we only
        # optimize the parameters that actually influence the score.
        def objective(trial):
            opt_param = trial.suggest_float('opt_param', 0.0, 1.0)

            config_kwargs = {
                'opt_param': opt_param,
                'n_conformers': n_conformers,
                'minimize': minimize,
                'scoring_mode': scoring_mode,
                'aggregation': aggregation,
            }

            if scoring_mode == 'hybrid':
                alpha = trial.suggest_float('alpha', 0.3, 0.9)
                config_kwargs['alpha'] = alpha

            if scoring_mode == 'consensus_mol':
                # Only consensus_mol uses tolerance/occurrence for scoring
                config_kwargs['tolerance'] = trial.suggest_float(
                    'tolerance', 0.5, 4.0)
                config_kwargs['occurrence'] = trial.suggest_float(
                    'occurrence', 0.1, 1.0)

            config = EvaluationConfig(**config_kwargs)
            result = evaluator.evaluate(config)
            return result.roc_auc

        study = optuna.create_study(
            direction='maximize',
            sampler=optuna.samplers.TPESampler(seed=self.random_state)
        )
        study.optimize(objective, n_trials=n_trials, show_progress_bar=False)

        best_trial = study.best_trial
        best_config = {
            'ref_indices': ref_indices,
            'aggregation': aggregation,
            'minimize': minimize,
            'scoring_mode': scoring_mode,
            'n_conformers': n_conformers,
            **best_trial.params
        }

        # Re-evaluate best to get full metrics
        config_kwargs = {
            k: v for k, v in best_trial.params.items()
            if k in EvaluationConfig.__dataclass_fields__
        }
        config_kwargs.update({
            'minimize': minimize,
            'scoring_mode': scoring_mode,
            'aggregation': aggregation,
            'n_conformers': n_conformers,
        })
        best_eval_config = EvaluationConfig(**config_kwargs)
        best_result = evaluator.evaluate(best_eval_config)

        return {
            'ref_indices': ref_indices,
            'aggregation': aggregation,
            'minimize': minimize,
            'scoring_mode': scoring_mode,
            'best_auc': best_result.roc_auc,
            'pruned': False,
            'best_config': best_config,
            'best_result': best_result,
            'study': study
        }

    def optimize(
        self,
        n_trials: int = 30,
        n_conformers: int = 10,
        multi_fidelity: bool = True,
        min_subset_size: int = 3,
        scoring_modes: Optional[List[str]] = None,
        aggregations: Optional[List[str]] = None,
        minimize_options: Optional[List[bool]] = None,
        min_auc_threshold: float = 0.6,
        verbose: bool = True
    ) -> CombinatorialResult:
        """Run combinatorial optimization across all degrees of freedom.

        Args:
            n_trials: Optuna trials per discrete combination.
            n_conformers: Fixed number of conformers for 3D scoring.
                Kept constant across trials to allow conformer caching
                (varying it would regenerate 3D coords per trial).
            multi_fidelity: Use pharm2d pre-screening for pruning.
            min_subset_size: Minimum reference subset size.
            scoring_modes: Scoring modes to try (default: ['reference', 'hybrid']).
            aggregations: Aggregation methods (default: ['max', 'mean']).
            minimize_options: Minimize options (default: [False, True]).
            min_auc_threshold: Minimum pharm2d AUC for multi-fidelity.
            verbose: Print progress.

        Returns:
            CombinatorialResult with best configuration and analysis.
        """
        start_time = time.time()

        if scoring_modes is None:
            scoring_modes = ['reference', 'hybrid']
        if aggregations is None:
            aggregations = ['max', 'mean']
        if minimize_options is None:
            minimize_options = [False, True]

        n_refs = len(self.reference_mols)
        ref_subsets = _generate_reference_subsets(n_refs, min_k=min_subset_size)

        # Build discrete combos
        discrete_combos = list(itertools.product(
            ref_subsets, aggregations, minimize_options, scoring_modes
        ))
        total_combos = len(discrete_combos)

        if verbose:
            print(f"Combinatorial search space:")
            print(f"  Reference subsets: {len(ref_subsets)}")
            print(f"  Aggregations: {aggregations}")
            print(f"  Minimize options: {minimize_options}")
            print(f"  Scoring modes: {scoring_modes}")
            print(f"  n_conformers (fixed): {n_conformers}")
            print(f"  Total discrete combos: {total_combos}")
            print(f"  Optuna trials per combo: {n_trials}")
            print(f"  Max total evaluations: {total_combos * n_trials}")
            print()

        all_results = []
        best_auc = 0.0
        best_result_info = None

        for i, (refs, agg, mini, mode) in enumerate(discrete_combos):
            if verbose:
                ref_str = f"refs={list(refs)}"
                print(
                    f"[{i+1}/{total_combos}] {ref_str}, "
                    f"agg={agg}, min={mini}, mode={mode}",
                    end=" ... ",
                    flush=True
                )

            try:
                result = self._evaluate_discrete_combo(
                    ref_indices=refs,
                    aggregation=agg,
                    minimize=mini,
                    scoring_mode=mode,
                    n_trials=n_trials,
                    n_conformers=n_conformers,
                    multi_fidelity=multi_fidelity,
                    min_auc_threshold=min_auc_threshold
                )
            except Exception as e:
                logger.warning("Combo %d failed: %s", i, e)
                if verbose:
                    print(f"FAILED: {e}")
                continue

            all_results.append(result)

            if verbose:
                status = "PRUNED" if result['pruned'] else f"AUC={result['best_auc']:.4f}"
                print(status)

            if result['best_auc'] > best_auc:
                best_auc = result['best_auc']
                best_result_info = result

        wall_time = time.time() - start_time

        # Build output
        if best_result_info is None or best_result_info['best_config'] is None:
            return CombinatorialResult(
                best_config={},
                best_metrics={'roc_auc': 0.5},
                total_evaluations=sum(
                    1 for r in all_results if not r['pruned']
                ) * n_trials,
                wall_time_sec=wall_time
            )

        best_eval = best_result_info['best_result']
        best_metrics = {
            'roc_auc': best_eval.roc_auc,
            'bedroc': best_eval.bedroc,
            'ef_1': best_eval.ef_1,
            'ef_5': best_eval.ef_5,
            'ef_10': best_eval.ef_10,
            'n_features': best_eval.n_features,
        }

        # Top-k configs
        non_pruned = [r for r in all_results if not r['pruned']]
        non_pruned.sort(key=lambda r: r['best_auc'], reverse=True)
        top_k = []
        for r in non_pruned[:10]:
            top_k.append({
                'config': r['best_config'],
                'roc_auc': r['best_auc'],
            })

        # Feature importance (variance-based)
        feature_importance = self._compute_feature_importance(non_pruned)

        # Pareto front (AUC vs eval time)
        pareto_front = self._compute_pareto_front(non_pruned)

        total_evals = sum(
            n_trials for r in all_results if not r['pruned']
        ) + sum(1 for r in all_results if r['pruned'])

        if verbose:
            print(f"\n{'='*60}")
            print(f"OPTIMIZATION COMPLETE")
            print(f"{'='*60}")
            print(f"Best AUC: {best_metrics['roc_auc']:.4f}")
            print(f"Best BEDROC: {best_metrics['bedroc']:.4f}")
            print(f"Best config: {best_result_info['best_config']}")
            print(f"Total evaluations: {total_evals}")
            print(f"Wall time: {wall_time:.1f}s")
            if feature_importance:
                print(f"\nFeature importance:")
                for name, imp in sorted(
                    feature_importance.items(), key=lambda x: -x[1]
                ):
                    print(f"  {name}: {imp:.3f}")

        return CombinatorialResult(
            best_config=best_result_info['best_config'],
            best_metrics=best_metrics,
            pareto_front=pareto_front,
            feature_importance=feature_importance,
            top_k_configs=top_k,
            total_evaluations=total_evals,
            wall_time_sec=wall_time
        )

    def _compute_feature_importance(
        self,
        results: List[Dict]
    ) -> Dict[str, float]:
        """Estimate feature importance from variance in AUC across DoFs."""
        if len(results) < 3:
            return {}

        importance = {}

        # Discrete DoFs: variance of AUC grouped by each
        for dof_name in ['aggregation', 'scoring_mode', 'minimize']:
            groups = {}
            for r in results:
                val = r.get(dof_name, r.get('best_config', {}).get(dof_name))
                if val is not None:
                    groups.setdefault(str(val), []).append(r['best_auc'])
            if len(groups) >= 2:
                group_means = [np.mean(v) for v in groups.values()]
                importance[dof_name] = float(np.var(group_means))

        # Reference subset: grouped by size
        size_groups = {}
        for r in results:
            size = len(r.get('ref_indices', ()))
            size_groups.setdefault(size, []).append(r['best_auc'])
        if len(size_groups) >= 2:
            size_means = [np.mean(v) for v in size_groups.values()]
            importance['ref_subset_size'] = float(np.var(size_means))

        # Normalize
        total = sum(importance.values())
        if total > 0:
            importance = {k: v / total for k, v in importance.items()}

        return importance

    def _compute_pareto_front(
        self,
        results: List[Dict]
    ) -> List[Tuple[Dict, float, float]]:
        """Find Pareto front of AUC vs evaluation time."""
        if not results:
            return []

        points = []
        for r in results:
            if r.get('best_result') is not None:
                auc = r['best_auc']
                eval_time = r['best_result'].eval_time_sec
                points.append((r.get('best_config', {}), auc, eval_time))

        if not points:
            return []

        # Sort by AUC descending
        points.sort(key=lambda x: -x[1])

        # Filter Pareto-optimal (no point dominates on both AUC and time)
        pareto = [points[0]]
        min_time = points[0][2]

        for config, auc, t in points[1:]:
            if t < min_time:
                pareto.append((config, auc, t))
                min_time = t

        return pareto

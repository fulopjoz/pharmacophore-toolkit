"""Optuna-based multi-objective pharmacophore optimization.

This module implements two optimization approaches using the Optuna framework:

1. Gaussian Process (GP-Sampler): Bayesian optimization with GP surrogate model.
   Recommended for smaller budgets (~50 trials). Reaches optimum efficiently.

2. NSGA-II (Evolutionary): Multi-objective genetic algorithm with non-dominated sorting.
   Recommended for larger budgets (~100 trials). Explores Pareto front broadly.

Both approaches maximize ROC-AUC and BEDROC simultaneously, producing a Pareto
front of non-dominated solutions.

Example:
    >>> from pharmacophore.optuna_optimizer import OptunaPharmacophoreOptimizer
    >>> optimizer = OptunaPharmacophoreOptimizer(ref_mols, actives, decoys)
    >>> result = optimizer.optimize(sampler='gp', n_trials=50)
    >>> print(f"Best AUC: {result['best_auc']:.4f}")
    >>> print(f"Pareto front size: {len(result['pareto_front'])}")
"""

from typing import List, Dict, Any, Optional
import time
import numpy as np
from rdkit import Chem

try:
    import optuna
    from optuna.samplers import GPSampler, NSGAIISampler, TPESampler
    HAS_OPTUNA = True
except ImportError:
    HAS_OPTUNA = False

from .evaluation import UnifiedEvaluator, EvaluationConfig


class OptunaPharmacophoreOptimizer:
    """Multi-objective pharmacophore optimization using Optuna.

    Supports two sampling strategies:
    - 'gp': Gaussian Process (GPSampler) for efficient optimization
    - 'nsga2': NSGA-II evolutionary algorithm for broad Pareto exploration

    Optimizes 6 parameters simultaneously:
    - tolerance (float): 0.5-4.0 Å
    - occurrence (float): 0.1-1.0
    - shape_weight (float): 0.3-0.9
    - opt_param (float): 0.0-1.0
    - linkage (categorical): average, complete, single, ward
    - n_conformers (int): 5-50

    Objectives:
    - Maximize ROC-AUC (overall discrimination)
    - Maximize BEDROC (early recognition)

    Attributes:
        evaluator: UnifiedEvaluator instance for scoring configurations.
        study: Optuna study object (created during optimization).
        start_time: Optimization start timestamp.

    Example:
        >>> optimizer = OptunaPharmacophoreOptimizer(refs, actives, decoys)
        >>> result = optimizer.optimize(sampler='gp', n_trials=50, verbose=True)
        >>> print(f"Found {len(result['pareto_front'])} Pareto-optimal solutions")
    """

    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol],
        random_state: int = 42
    ):
        """Initialize optimizer with molecules.

        Args:
            reference_mols: Aligned reference molecules for consensus.
            actives: Active compounds (will generate conformers if needed).
            decoys: Decoy compounds (will generate conformers if needed).
            random_state: Random seed for reproducibility.

        Raises:
            ImportError: If Optuna is not installed.
        """
        if not HAS_OPTUNA:
            raise ImportError(
                "Optuna required. Install with: pip install optuna"
            )

        self.evaluator = UnifiedEvaluator(
            reference_mols, actives, decoys, random_state=random_state
        )
        self.random_state = random_state
        self.study: Optional[optuna.Study] = None
        self.start_time: Optional[float] = None

    def _objective(self, trial: optuna.Trial) -> tuple:
        """Optuna objective function.

        Args:
            trial: Optuna trial object for suggesting parameters.

        Returns:
            Tuple of (roc_auc, bedroc) for multi-objective optimization.
        """
        # Suggest parameter values
        tolerance = trial.suggest_float('tolerance', 0.5, 4.0)
        occurrence = trial.suggest_float('occurrence', 0.1, 1.0)
        shape_weight = trial.suggest_float('shape_weight', 0.3, 0.9)
        opt_param = trial.suggest_float('opt_param', 0.0, 1.0)
        linkage = trial.suggest_categorical('linkage', ['average', 'complete', 'single', 'ward'])
        n_conformers = trial.suggest_int('n_conformers', 5, 50)

        # Create configuration
        config = EvaluationConfig(
            tolerance=tolerance,
            occurrence=occurrence,
            shape_weight=shape_weight,
            opt_param=opt_param,
            linkage=linkage,
            n_conformers=n_conformers
        )

        # Evaluate
        result = self.evaluator.evaluate(config)

        # Store additional info in trial
        trial.set_user_attr('n_features', result.n_features)
        trial.set_user_attr('ef_1', result.ef_1)
        trial.set_user_attr('ef_5', result.ef_5)
        trial.set_user_attr('ef_10', result.ef_10)
        trial.set_user_attr('eval_time_sec', result.eval_time_sec)

        # Return objectives (to maximize)
        return result.roc_auc, result.bedroc

    def optimize(
        self,
        sampler: str = 'gp',
        n_trials: int = 50,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """Run multi-objective optimization.

        Args:
            sampler: Sampling strategy ('gp' or 'nsga2').
            n_trials: Number of optimization trials.
            verbose: Print progress during optimization.

        Returns:
            Dictionary with:
            - sampler: Sampler type used
            - n_trials: Number of trials run
            - wall_time_sec: Total optimization time (seconds)
            - pareto_front: List of Pareto-optimal solutions
            - best_auc: Best ROC-AUC achieved
            - best_bedroc: Best BEDROC achieved
            - best_auc_params: Parameters for best AUC
            - best_bedroc_params: Parameters for best BEDROC
            - history: All trial results
            - study: Optuna study object

        Example:
            >>> result = optimizer.optimize(sampler='gp', n_trials=50)
            >>> print(f"Wall time: {result['wall_time_sec']:.1f}s")
            >>> print(f"Pareto solutions: {len(result['pareto_front'])}")
        """
        if sampler not in ['gp', 'nsga2']:
            raise ValueError(f"sampler must be 'gp' or 'nsga2', got {sampler}")

        # Create sampler
        if sampler == 'gp':
            try:
                import warnings as _w
                with _w.catch_warnings():
                    _w.simplefilter("ignore", category=FutureWarning)
                    _w.filterwarnings("ignore", message=".*experimental.*")
                    sampler_obj = GPSampler(seed=self.random_state)
                sampler_name = 'Gaussian Process'
            except (ImportError, AttributeError):
                # Fallback if GPSampler unavailable (Optuna < 3.4)
                if verbose:
                    print("Warning: GPSampler unavailable, using TPESampler as fallback")
                sampler_obj = TPESampler(seed=self.random_state)
                sampler_name = 'Tree-structured Parzen Estimator (fallback)'
        else:  # nsga2
            sampler_obj = NSGAIISampler(seed=self.random_state)
            sampler_name = 'NSGA-II Evolutionary'

        if verbose:
            print(f"{'='*60}")
            print(f"Multi-Objective Pharmacophore Optimization")
            print(f"{'='*60}")
            print(f"Sampler: {sampler_name}")
            print(f"Trials: {n_trials}")
            print(f"Objectives: Maximize ROC-AUC + BEDROC")
            print(f"Parameters: 6 (tolerance, occurrence, shape_weight, opt_param, linkage, n_conformers)")
            print(f"{'='*60}\n")

        # Create study
        self.study = optuna.create_study(
            directions=['maximize', 'maximize'],  # (AUC, BEDROC)
            sampler=sampler_obj
        )

        # Run optimization
        self.start_time = time.time()

        if verbose:
            # Custom callback for progress
            def print_progress(study, trial):
                elapsed = time.time() - self.start_time
                auc, bedroc = trial.values
                n_feat = trial.user_attrs.get('n_features', 0)
                print(
                    f"Trial {trial.number:3d}/{n_trials}: "
                    f"AUC={auc:.4f}, BEDROC={bedroc:.4f}, "
                    f"n_feat={n_feat}, elapsed={elapsed:.1f}s"
                )

            self.study.optimize(
                self._objective,
                n_trials=n_trials,
                callbacks=[print_progress],
                show_progress_bar=False
            )
        else:
            self.study.optimize(self._objective, n_trials=n_trials)

        wall_time = time.time() - self.start_time

        # Extract Pareto front
        pareto_trials = self.study.best_trials
        pareto_front = []

        for trial in pareto_trials:
            auc, bedroc = trial.values
            pareto_front.append({
                'params': trial.params,
                'roc_auc': auc,
                'bedroc': bedroc,
                'n_features': trial.user_attrs.get('n_features', 0),
                'ef_1': trial.user_attrs.get('ef_1', 0.0),
                'ef_5': trial.user_attrs.get('ef_5', 0.0),
                'ef_10': trial.user_attrs.get('ef_10', 0.0)
            })

        # Find best individual objectives
        all_trials = self.study.trials
        best_auc_trial = max(all_trials, key=lambda t: t.values[0])
        best_bedroc_trial = max(all_trials, key=lambda t: t.values[1])

        best_auc = best_auc_trial.values[0]
        best_bedroc = best_bedroc_trial.values[1]
        best_auc_params = best_auc_trial.params
        best_bedroc_params = best_bedroc_trial.params

        # Build history
        history = []
        for trial in all_trials:
            if trial.state == optuna.trial.TrialState.COMPLETE:
                auc, bedroc = trial.values
                history.append({
                    'trial': trial.number,
                    'params': trial.params,
                    'roc_auc': auc,
                    'bedroc': bedroc,
                    'n_features': trial.user_attrs.get('n_features', 0),
                    'eval_time_sec': trial.user_attrs.get('eval_time_sec', 0.0)
                })

        if verbose:
            print(f"\n{'='*60}")
            print(f"OPTIMIZATION COMPLETE")
            print(f"{'='*60}")
            print(f"Wall time: {wall_time:.2f} seconds")
            print(f"Time per trial: {wall_time/n_trials:.2f} seconds")
            print(f"Pareto front size: {len(pareto_front)} solutions")
            print(f"\nBest ROC-AUC: {best_auc:.4f}")
            print(f"  Parameters: tolerance={best_auc_params['tolerance']:.2f}, "
                  f"occurrence={best_auc_params['occurrence']:.2f}, "
                  f"shape_weight={best_auc_params['shape_weight']:.2f}")
            print(f"\nBest BEDROC: {best_bedroc:.4f}")
            print(f"  Parameters: tolerance={best_bedroc_params['tolerance']:.2f}, "
                  f"occurrence={best_bedroc_params['occurrence']:.2f}, "
                  f"shape_weight={best_bedroc_params['shape_weight']:.2f}")
            print(f"{'='*60}")

        return {
            'sampler': sampler,
            'n_trials': n_trials,
            'wall_time_sec': wall_time,
            'pareto_front': pareto_front,
            'best_auc': best_auc,
            'best_bedroc': best_bedroc,
            'best_auc_params': best_auc_params,
            'best_bedroc_params': best_bedroc_params,
            'history': history,
            'study': self.study
        }

    def get_parameter_importance(self) -> Dict[str, float]:
        """Estimate parameter importance using fANOVA.

        Requires optimization to have been run first.

        Returns:
            Dictionary mapping parameter names to importance scores (0-1).
            Higher values indicate more important parameters.

        Raises:
            ValueError: If optimize() hasn't been called yet.

        Example:
            >>> result = optimizer.optimize(sampler='gp', n_trials=50)
            >>> importance = optimizer.get_parameter_importance()
            >>> print(f"Most important: {max(importance, key=importance.get)}")
        """
        if self.study is None:
            raise ValueError("Run optimize() first before getting parameter importance")

        try:
            # Use Optuna's built-in importance evaluator
            # Target first objective (ROC-AUC)
            evaluator = optuna.importance.FanovaImportanceEvaluator()
            importance = optuna.importance.get_param_importances(
                self.study,
                evaluator=evaluator,
                target=lambda t: t.values[0]  # ROC-AUC
            )
            return importance
        except Exception:
            # Fallback to simple variance-based importance
            return self._simple_parameter_importance()

    def _simple_parameter_importance(self) -> Dict[str, float]:
        """Simple variance-based parameter importance (fallback)."""
        if self.study is None:
            return {}

        completed_trials = [t for t in self.study.trials if t.state == optuna.trial.TrialState.COMPLETE]

        if len(completed_trials) < 10:
            return {}

        # Extract parameter values and AUC scores
        param_names = ['tolerance', 'occurrence', 'shape_weight', 'opt_param', 'linkage', 'n_conformers']
        importance = {}

        for param in param_names:
            if param == 'linkage':
                # Categorical - use group variance
                groups = {}
                for trial in completed_trials:
                    val = trial.params[param]
                    auc = trial.values[0]
                    if val not in groups:
                        groups[val] = []
                    groups[val].append(auc)

                # Calculate between-group variance
                group_means = [np.mean(g) for g in groups.values()]
                importance[param] = float(np.var(group_means)) if len(group_means) > 1 else 0.0
            else:
                # Continuous - use correlation
                param_values = [trial.params[param] for trial in completed_trials]
                auc_values = [trial.values[0] for trial in completed_trials]
                importance[param] = abs(float(np.corrcoef(param_values, auc_values)[0, 1]))

        # Normalize
        total = sum(importance.values())
        if total > 0:
            importance = {k: v/total for k, v in importance.items()}

        return importance

    def get_pareto_front(self) -> List[Dict]:
        """Get Pareto-optimal solutions.

        Requires optimization to have been run first.

        Returns:
            List of dictionaries with params and objective values.

        Raises:
            ValueError: If optimize() hasn't been called yet.
        """
        if self.study is None:
            raise ValueError("Run optimize() first before getting Pareto front")

        pareto_trials = self.study.best_trials
        pareto_front = []

        for trial in pareto_trials:
            auc, bedroc = trial.values
            pareto_front.append({
                'params': trial.params,
                'roc_auc': auc,
                'bedroc': bedroc,
                'n_features': trial.user_attrs.get('n_features', 0)
            })

        return pareto_front

    def __repr__(self) -> str:
        n_refs = len(self.evaluator.reference_mols)
        n_actives = len(self.evaluator.actives)
        n_decoys = len(self.evaluator.decoys)
        n_trials = len(self.study.trials) if self.study else 0

        return (
            f"OptunaPharmacophoreOptimizer("
            f"n_refs={n_refs}, n_actives={n_actives}, "
            f"n_decoys={n_decoys}, trials={n_trials})"
        )

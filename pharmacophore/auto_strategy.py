"""Unified strategy selection for pharmacophore model optimization.

Automatically selects and evaluates the best combination of clustering
strategy and scoring approach for a given dataset. Runs a lightweight
tournament across strategies using multi-fidelity evaluation.

This module integrates:
- PharmacophoreConsensus (agglomerative hierarchical clustering)
- EnsembleConsensus (stability-aware ensemble voting)
- Hungarian matching and Wasserstein distance scoring
- Multi-level cascading evaluation

Named StrategySelector to avoid collision with the existing
AutoPharmacophoreOptimizer in auto_optimizer.py, which optimizes
parameters within a single strategy.
"""

import time
import numpy as np
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from rdkit import Chem

from .consensus import PharmacophoreConsensus
from .ensemble_consensus import EnsembleConsensus
from .evaluation import (
    UnifiedEvaluator, EvaluationConfig, EvaluationResult, compute_sdbw
)


@dataclass
class StrategyResult:
    """Result from a single strategy evaluation.

    Attributes:
        strategy_name: Name of the clustering strategy.
        eval_result: Full evaluation metrics.
        sdbw: S_Dbw cluster validation score (lower = better).
        features: Consensus pharmacophore features.
        stability_scores: Per-feature stability (ensemble only).
    """
    strategy_name: str
    eval_result: EvaluationResult
    sdbw: float = float('inf')
    features: Optional[List[List]] = None
    stability_scores: Optional[List[float]] = None


class StrategySelector:
    """Automatically select the best pharmacophore optimization strategy.

    Runs a lightweight tournament across clustering strategies and scoring
    approaches. Each strategy is evaluated on the provided dataset, and the
    best one is selected based on ROC-AUC.

    Available clustering strategies:
        - 'agglomerative': Standard hierarchical clustering (fast, deterministic)
        - 'ensemble': Stability-aware ensemble voting (robust, slower)

    Available scoring modes:
        - 'standard': AlignMol shape+color scoring
        - 'ensemble': Cascading multi-level scoring with OT

    Example:
        >>> selector = StrategySelector(ref_mols, actives, decoys)
        >>> best = selector.select_best()
        >>> print(f"Best strategy: {best.strategy_name}, AUC: {best.eval_result.roc_auc:.4f}")
    """

    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        actives: Optional[List[Chem.Mol]] = None,
        decoys: Optional[List[Chem.Mol]] = None,
        random_state: int = 42,
        verbose: bool = True
    ):
        """Initialize strategy selector.

        Args:
            reference_mols: Aligned reference molecules.
            actives: Active compounds (required for evaluation).
            decoys: Decoy compounds (required for evaluation).
            random_state: Random seed for reproducibility.
            verbose: Print progress information.
        """
        self.reference_mols = reference_mols
        self.actives = actives
        self.decoys = decoys
        self.random_state = random_state
        self.verbose = verbose

        self._evaluator = None
        if actives and decoys:
            self._evaluator = UnifiedEvaluator(
                reference_mols, actives, decoys, random_state=random_state
            )

    def select_best(
        self,
        configs: Optional[List[EvaluationConfig]] = None,
        strategies: Optional[List[str]] = None
    ) -> StrategyResult:
        """Run tournament and return the best strategy.

        Args:
            configs: Configurations to try. If None, uses a default grid.
            strategies: Which strategies to evaluate. Default: all available.

        Returns:
            StrategyResult for the best strategy.
        """
        if self._evaluator is None:
            raise ValueError("Actives and decoys required for strategy selection")

        if configs is None:
            configs = self._default_configs()

        if strategies is None:
            strategies = ['agglomerative', 'ensemble']

        results = []
        for strategy in strategies:
            for config in configs:
                result = self._evaluate_strategy(strategy, config)
                results.append(result)

                if self.verbose:
                    print(
                        f"  {result.strategy_name:20s} | "
                        f"tol={config.tolerance:.1f} occ={config.occurrence:.1f} | "
                        f"AUC={result.eval_result.roc_auc:.4f} "
                        f"BEDROC={result.eval_result.bedroc:.4f} "
                        f"n_feat={result.eval_result.n_features}"
                    )

        # Select best by ROC-AUC
        best = max(results, key=lambda r: r.eval_result.roc_auc)

        if self.verbose:
            print(f"\nBest: {best.strategy_name} "
                  f"(AUC={best.eval_result.roc_auc:.4f})")

        return best

    def _evaluate_strategy(
        self,
        strategy: str,
        config: EvaluationConfig
    ) -> StrategyResult:
        """Evaluate a single strategy with a given config."""
        start_time = time.time()

        if strategy == 'agglomerative':
            consensus = PharmacophoreConsensus(
                tolerance=config.tolerance,
                occurrence_threshold=config.occurrence,
                linkage=config.linkage
            )
            features, metadata = consensus.generate_consensus(
                self.reference_mols, return_metadata=True
            )
            sdbw = compute_sdbw(metadata)
            stability = None

        elif strategy == 'ensemble':
            ec = EnsembleConsensus(
                n_runs=15,  # Reduced for speed in tournament
                tolerance_range=(config.tolerance * 0.75, config.tolerance * 1.25),
                occurrence_range=(max(0.1, config.occurrence - 0.15),
                                  min(1.0, config.occurrence + 0.15)),
                stability_threshold=0.3,
                random_state=self.random_state
            )
            features, stability = ec.generate_consensus_with_scores(
                self.reference_mols
            )
            sdbw = float('inf')  # Not available for ensemble

        else:
            raise ValueError(f"Unknown strategy: {strategy}")

        # Evaluate features
        if len(features) >= 2:
            eval_result = self._evaluator.evaluate_feature_subset(
                features,
                shape_weight=config.shape_weight,
                opt_param=config.opt_param,
                n_conformers=config.n_conformers
            )
        else:
            eval_result = EvaluationResult(
                config=config, roc_auc=0.5, bedroc=0.0,
                ef_1=0.0, ef_5=0.0, ef_10=0.0,
                n_features=len(features),
                eval_time_sec=time.time() - start_time
            )

        return StrategyResult(
            strategy_name=strategy,
            eval_result=eval_result,
            sdbw=sdbw,
            features=features,
            stability_scores=stability
        )

    def _default_configs(self) -> List[EvaluationConfig]:
        """Default configuration grid for tournament."""
        configs = []
        for tol in [1.5, 2.0, 2.5]:
            for occ in [0.3, 0.5, 0.7]:
                configs.append(EvaluationConfig(
                    tolerance=tol,
                    occurrence=occ,
                    shape_weight=0.5,
                    opt_param=0.5,
                    n_conformers=10  # Fast for tournament
                ))
        return configs

"""HypoGen-inspired pharmacophore optimization algorithm.

This module implements a 3-phase constructive-subtractive-optimization approach
inspired by Catalyst/Discovery Studio's HypoGen algorithm:

Phase 1 (Constructive): Enumerate feature subsets from consensus pharmacophore
Phase 2 (Subtractive): Filter out hypotheses with poor selectivity
Phase 3 (Optimization): Refine feature positions using simulated annealing

Unlike black-box optimizers (GP, NSGA-II), HypoGen uses pharmacophore-specific
heuristics and domain knowledge.

Example:
    >>> from pharmacophore.hypogen_optimizer import HypoGenOptimizer
    >>> optimizer = HypoGenOptimizer(ref_mols, actives, decoys)
    >>> result = optimizer.optimize(
    ...     consensus_tolerance=2.0,
    ...     consensus_occurrence=0.3,
    ...     min_features=3,
    ...     max_features=8
    ... )
    >>> print(f"Best hypothesis cost: {result['best_cost']:.4f}")
"""

from typing import List, Dict, Any, Tuple
import time
import itertools
import numpy as np
from rdkit import Chem

from .evaluation import UnifiedEvaluator, EvaluationResult
from .consensus import PharmacophoreConsensus


class HypoGenOptimizer:
    """HypoGen-inspired 3-phase pharmacophore optimization.

    Phases:
    1. Constructive: Generate candidate hypotheses from feature subsets
    2. Subtractive: Remove hypotheses that score decoys highly
    3. Optimization: Refine feature positions with simulated annealing

    Cost function:
        cost = 0.5*(1-AUC) + 0.3*(1-BEDROC) + 0.2*complexity
        where complexity = n_features / max_features

    Attributes:
        evaluator: UnifiedEvaluator for scoring hypotheses.
        random_state: Random seed for reproducibility.

    Example:
        >>> optimizer = HypoGenOptimizer(refs, actives, decoys)
        >>> result = optimizer.optimize()
        >>> print(f"Found {len(result['survivors'])} viable hypotheses")
    """

    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol],
        random_state: int = 42
    ):
        """Initialize HypoGen optimizer.

        Args:
            reference_mols: Aligned reference molecules for consensus.
            actives: Active compounds (will generate conformers if needed).
            decoys: Decoy compounds (will generate conformers if needed).
            random_state: Random seed for reproducibility.
        """
        self.evaluator = UnifiedEvaluator(
            reference_mols, actives, decoys, random_state=random_state
        )
        self.reference_mols = reference_mols
        self.random_state = random_state
        self.rng = np.random.RandomState(random_state)

    def _calculate_cost(
        self,
        result: EvaluationResult,
        max_features: int
    ) -> float:
        """Calculate HypoGen-style cost function.

        Lower cost is better. Penalizes poor discrimination and complexity.

        Args:
            result: Evaluation result with metrics.
            max_features: Maximum number of features for normalization.

        Returns:
            Cost value (lower is better).
        """
        # Discrimination penalties (inverted metrics)
        auc_penalty = 1.0 - result.roc_auc
        bedroc_penalty = 1.0 - result.bedroc

        # Complexity penalty (prefer fewer features)
        complexity = result.n_features / max_features if max_features > 0 else 0.0

        # Weighted cost
        cost = (
            0.5 * auc_penalty +
            0.3 * bedroc_penalty +
            0.2 * complexity
        )

        return cost

    def _enumerate_subsets(
        self,
        features: List,
        min_features: int,
        max_features: int,
        max_hypotheses: int = 500
    ) -> List[List]:
        """Enumerate feature subsets for Phase 1.

        Args:
            features: Full consensus feature list.
            min_features: Minimum features per hypothesis.
            max_features: Maximum features per hypothesis.
            max_hypotheses: Cap on number of hypotheses to generate.

        Returns:
            List of feature subsets (hypotheses).
        """
        subsets = []

        # Generate all combinations from min to max features
        for n_feat in range(min_features, min(max_features + 1, len(features) + 1)):
            combinations = list(itertools.combinations(features, n_feat))

            # Add to subsets
            for combo in combinations:
                subsets.append(list(combo))

                # Cap at max_hypotheses
                if len(subsets) >= max_hypotheses:
                    return subsets

        return subsets

    def _phase1_constructive(
        self,
        tolerance: float,
        occurrence: float,
        min_features: int,
        max_features: int,
        max_hypotheses: int,
        n_conformers: int,
        shape_weight: float,
        opt_param: float,
        verbose: bool
    ) -> Tuple[List[Dict], List, float]:
        """Phase 1: Generate candidate hypotheses from consensus.

        Args:
            tolerance: Consensus clustering tolerance.
            occurrence: Consensus occurrence threshold.
            min_features: Minimum features per hypothesis.
            max_features: Maximum features per hypothesis.
            max_hypotheses: Cap on hypotheses to generate.
            n_conformers: Conformers per molecule.
            shape_weight: Shape vs color weight.
            opt_param: AlignMol optimization parameter.
            verbose: Print progress.

        Returns:
            Tuple of (scored_hypotheses, full_consensus_features, phase_time).
        """
        phase_start = time.time()

        if verbose:
            print(f"\n{'='*60}")
            print(f"PHASE 1: CONSTRUCTIVE")
            print(f"{'='*60}")

        # Generate full consensus
        consensus = PharmacophoreConsensus(
            tolerance=tolerance,
            occurrence_threshold=occurrence
        )
        features = consensus.generate_consensus(self.reference_mols)

        if verbose:
            print(f"Consensus features: {len(features)}")
            print(f"Enumerating subsets ({min_features}-{max_features} features)...")

        # Enumerate subsets
        subsets = self._enumerate_subsets(
            features, min_features, max_features, max_hypotheses
        )

        if verbose:
            print(f"Generated {len(subsets)} candidate hypotheses")
            print(f"Evaluating hypotheses...")

        # Score each hypothesis
        hypotheses = []
        for i, subset in enumerate(subsets):
            result = self.evaluator.evaluate_feature_subset(
                subset,
                shape_weight=shape_weight,
                opt_param=opt_param,
                n_conformers=n_conformers
            )

            cost = self._calculate_cost(result, max_features)

            hypotheses.append({
                'id': i,
                'features': subset,
                'result': result,
                'cost': cost
            })

            if verbose and (i + 1) % 50 == 0:
                print(f"  Evaluated {i+1}/{len(subsets)} hypotheses...")

        phase_time = time.time() - phase_start

        if verbose:
            print(f"Phase 1 complete: {phase_time:.2f}s")
            print(f"Best cost: {min(h['cost'] for h in hypotheses):.4f}")
            print(f"Best AUC: {max(h['result'].roc_auc for h in hypotheses):.4f}")

        return hypotheses, features, phase_time

    def _phase2_subtractive(
        self,
        hypotheses: List[Dict],
        auc_threshold: float,
        bedroc_threshold: float,
        verbose: bool
    ) -> Tuple[List[Dict], float]:
        """Phase 2: Remove hypotheses with poor selectivity.

        Args:
            hypotheses: Scored hypotheses from Phase 1.
            auc_threshold: Minimum AUC to survive.
            bedroc_threshold: Minimum BEDROC to survive.
            verbose: Print progress.

        Returns:
            Tuple of (surviving_hypotheses, phase_time).
        """
        phase_start = time.time()

        if verbose:
            print(f"\n{'='*60}")
            print(f"PHASE 2: SUBTRACTIVE")
            print(f"{'='*60}")
            print(f"Filtering {len(hypotheses)} hypotheses...")
            print(f"Thresholds: AUC >= {auc_threshold}, BEDROC >= {bedroc_threshold}")

        # Filter by thresholds
        survivors = [
            h for h in hypotheses
            if h['result'].roc_auc >= auc_threshold and h['result'].bedroc >= bedroc_threshold
        ]

        phase_time = time.time() - phase_start

        if verbose:
            print(f"Survivors: {len(survivors)}/{len(hypotheses)}")
            print(f"Phase 2 complete: {phase_time:.2f}s")

        return survivors, phase_time

    def _phase3_optimization(
        self,
        survivors: List[Dict],
        n_top: int,
        n_sa_iters: int,
        max_features: int,
        shape_weight: float,
        opt_param: float,
        n_conformers: int,
        verbose: bool
    ) -> Tuple[Dict, List[Dict], float]:
        """Phase 3: Refine top hypotheses with simulated annealing.

        Args:
            survivors: Surviving hypotheses from Phase 2.
            n_top: Number of top hypotheses to optimize.
            n_sa_iters: SA iterations per hypothesis.
            max_features: Maximum features for cost calculation.
            shape_weight: Shape vs color weight.
            opt_param: AlignMol optimization parameter.
            n_conformers: Conformers per molecule.
            verbose: Print progress.

        Returns:
            Tuple of (best_hypothesis, optimized_hypotheses, phase_time).
        """
        phase_start = time.time()

        if verbose:
            print(f"\n{'='*60}")
            print(f"PHASE 3: OPTIMIZATION")
            print(f"{'='*60}")
            print(f"Refining top {n_top} hypotheses with simulated annealing...")

        # Sort by cost
        survivors_sorted = sorted(survivors, key=lambda h: h['cost'])
        top_hypotheses = survivors_sorted[:n_top]

        optimized = []

        for i, hypothesis in enumerate(top_hypotheses):
            if verbose:
                print(f"\nOptimizing hypothesis {i+1}/{len(top_hypotheses)}...")
                print(f"  Initial cost: {hypothesis['cost']:.4f}, AUC: {hypothesis['result'].roc_auc:.4f}")

            # Run simulated annealing
            optimized_hyp = self._simulated_annealing(
                hypothesis['features'],
                n_iters=n_sa_iters,
                max_features=max_features,
                shape_weight=shape_weight,
                opt_param=opt_param,
                n_conformers=n_conformers,
                verbose=verbose
            )

            optimized.append(optimized_hyp)

            if verbose:
                print(f"  Final cost: {optimized_hyp['cost']:.4f}, AUC: {optimized_hyp['result'].roc_auc:.4f}")

        # Find best overall
        best_hypothesis = min(optimized, key=lambda h: h['cost'])

        phase_time = time.time() - phase_start

        if verbose:
            print(f"\nPhase 3 complete: {phase_time:.2f}s")
            print(f"Best hypothesis cost: {best_hypothesis['cost']:.4f}")
            print(f"Best AUC: {best_hypothesis['result'].roc_auc:.4f}")
            print(f"Best BEDROC: {best_hypothesis['result'].bedroc:.4f}")

        return best_hypothesis, optimized, phase_time

    def _simulated_annealing(
        self,
        features: List,
        n_iters: int,
        max_features: int,
        shape_weight: float,
        opt_param: float,
        n_conformers: int,
        verbose: bool = False
    ) -> Dict:
        """Simulated annealing to refine feature positions.

        Args:
            features: Initial feature list.
            n_iters: Number of SA iterations.
            max_features: Maximum features for cost calculation.
            shape_weight: Shape vs color weight.
            opt_param: AlignMol optimization parameter.
            n_conformers: Conformers per molecule.
            verbose: Print SA progress.

        Returns:
            Optimized hypothesis dictionary.
        """
        # Evaluate initial state
        current_result = self.evaluator.evaluate_feature_subset(
            features, shape_weight, opt_param, n_conformers
        )
        current_cost = self._calculate_cost(current_result, max_features)

        best_features = [f[:] for f in features]  # Deep copy
        best_result = current_result
        best_cost = current_cost

        # SA parameters
        T = 1.0  # Initial temperature
        alpha = 0.95  # Cooling rate

        for iteration in range(n_iters):
            # Perturb features
            perturbed_features = []
            for feat in features:
                feat_type, indices, x, y, z = feat
                # Add Gaussian noise to coordinates
                dx, dy, dz = self.rng.normal(0, 0.5, 3)
                new_feat = [feat_type, indices, x + dx, y + dy, z + dz]
                perturbed_features.append(new_feat)

            # Evaluate perturbed state
            try:
                new_result = self.evaluator.evaluate_feature_subset(
                    perturbed_features, shape_weight, opt_param, n_conformers
                )
                new_cost = self._calculate_cost(new_result, max_features)

                # Metropolis acceptance criterion
                delta_cost = new_cost - current_cost
                if delta_cost < 0 or self.rng.rand() < np.exp(-delta_cost / T):
                    # Accept new state
                    features = perturbed_features
                    current_result = new_result
                    current_cost = new_cost

                    # Update best if improved
                    if new_cost < best_cost:
                        best_features = [f[:] for f in features]
                        best_result = new_result
                        best_cost = new_cost

            except Exception:
                # Skip invalid perturbations
                pass

            # Cool down
            T *= alpha

        return {
            'features': best_features,
            'result': best_result,
            'cost': best_cost
        }

    def optimize(
        self,
        consensus_tolerance: float = 2.0,
        consensus_occurrence: float = 0.3,
        min_features: int = 3,
        max_features: int = 8,
        max_hypotheses: int = 500,
        n_top: int = 10,
        n_sa_iters: int = 200,
        auc_threshold: float = 0.55,
        bedroc_threshold: float = 0.05,
        shape_weight: float = 0.5,
        opt_param: float = 0.5,
        n_conformers: int = 25,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """Run HypoGen 3-phase optimization.

        Args:
            consensus_tolerance: Tolerance for initial consensus generation.
            consensus_occurrence: Occurrence threshold for consensus.
            min_features: Minimum features per hypothesis.
            max_features: Maximum features per hypothesis.
            max_hypotheses: Cap on Phase 1 hypotheses.
            n_top: Number of top hypotheses to optimize in Phase 3.
            n_sa_iters: SA iterations per hypothesis.
            auc_threshold: Minimum AUC for Phase 2 survival.
            bedroc_threshold: Minimum BEDROC for Phase 2 survival.
            shape_weight: Shape vs color weight (0-1).
            opt_param: AlignMol optimization parameter (0-1).
            n_conformers: Conformers per query molecule.
            verbose: Print detailed progress.

        Returns:
            Dictionary with:
            - best_hypothesis: Best hypothesis found
            - best_result: Evaluation result for best hypothesis
            - best_cost: Cost of best hypothesis
            - phase1_hypotheses: All Phase 1 hypotheses
            - phase2_survivors: Surviving hypotheses after Phase 2
            - phase3_optimized: Refined hypotheses from Phase 3
            - phase_times: Dictionary of timings per phase
            - total_time: Total optimization time
            - n_evaluations: Total number of evaluations

        Example:
            >>> result = optimizer.optimize(verbose=True)
            >>> print(f"Best AUC: {result['best_result'].roc_auc:.4f}")
            >>> print(f"Total time: {result['total_time']:.1f}s")
        """
        opt_start = time.time()

        if verbose:
            print(f"{'='*60}")
            print(f"HypoGen 3-Phase Pharmacophore Optimization")
            print(f"{'='*60}")
            print(f"Parameters:")
            print(f"  Consensus: tolerance={consensus_tolerance}, occurrence={consensus_occurrence}")
            print(f"  Features: {min_features}-{max_features}")
            print(f"  Max hypotheses: {max_hypotheses}")
            print(f"  Top-N for refinement: {n_top}")
            print(f"  SA iterations: {n_sa_iters}")

        # Phase 1: Constructive
        phase1_hypotheses, full_features, phase1_time = self._phase1_constructive(
            consensus_tolerance, consensus_occurrence,
            min_features, max_features, max_hypotheses,
            n_conformers, shape_weight, opt_param, verbose
        )

        # Phase 2: Subtractive
        phase2_survivors, phase2_time = self._phase2_subtractive(
            phase1_hypotheses, auc_threshold, bedroc_threshold, verbose
        )

        # Phase 3: Optimization
        if len(phase2_survivors) > 0:
            best_hypothesis, phase3_optimized, phase3_time = self._phase3_optimization(
                phase2_survivors, n_top, n_sa_iters, max_features,
                shape_weight, opt_param, n_conformers, verbose
            )
        else:
            # No survivors - return best from Phase 1
            if verbose:
                print("\nWarning: No survivors from Phase 2. Returning best Phase 1 hypothesis.")
            best_hypothesis = min(phase1_hypotheses, key=lambda h: h['cost'])
            phase3_optimized = []
            phase3_time = 0.0

        total_time = time.time() - opt_start

        # Count total evaluations
        n_evaluations = len(phase1_hypotheses) + n_top * n_sa_iters

        if verbose:
            print(f"\n{'='*60}")
            print(f"HYPOGEN OPTIMIZATION COMPLETE")
            print(f"{'='*60}")
            print(f"Total time: {total_time:.2f}s")
            print(f"Total evaluations: {n_evaluations}")
            print(f"Phase times: P1={phase1_time:.1f}s, P2={phase2_time:.1f}s, P3={phase3_time:.1f}s")
            print(f"\nBest Hypothesis:")
            print(f"  Cost: {best_hypothesis['cost']:.4f}")
            print(f"  ROC-AUC: {best_hypothesis['result'].roc_auc:.4f}")
            print(f"  BEDROC: {best_hypothesis['result'].bedroc:.4f}")
            print(f"  Features: {best_hypothesis['result'].n_features}")
            print(f"  EF@1%: {best_hypothesis['result'].ef_1:.2f}")
            print(f"{'='*60}")

        return {
            'best_hypothesis': best_hypothesis['features'],
            'best_result': best_hypothesis['result'],
            'best_cost': best_hypothesis['cost'],
            'phase1_hypotheses': phase1_hypotheses,
            'phase2_survivors': phase2_survivors,
            'phase3_optimized': phase3_optimized,
            'phase_times': {
                'phase1': phase1_time,
                'phase2': phase2_time,
                'phase3': phase3_time
            },
            'total_time': total_time,
            'n_evaluations': n_evaluations
        }

    def __repr__(self) -> str:
        n_refs = len(self.evaluator.reference_mols)
        n_actives = len(self.evaluator.actives)
        n_decoys = len(self.evaluator.decoys)

        return (
            f"HypoGenOptimizer("
            f"n_refs={n_refs}, n_actives={n_actives}, "
            f"n_decoys={n_decoys})"
        )

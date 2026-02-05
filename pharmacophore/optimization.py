"""Consensus pharmacophore parameter optimization.

This module provides automatic parameter optimization for consensus
pharmacophore generation, finding optimal tolerance and occurrence_threshold
values that maximize separation between actives and decoys.

Uses grid search over parameter space with enrichment factor and ROC-AUC
as evaluation metrics.
"""

import warnings

import numpy as np
from typing import List, Dict, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol
from sklearn.metrics import roc_auc_score

from .consensus import PharmacophoreConsensus
from .mol_converter import PharmacophoreToMol
from .screening_metrics import bedroc, youden_index, enrichment_factor


class ConsensusOptimizer:
    """Optimize consensus pharmacophore parameters for active/decoy separation.

    This class performs grid search over tolerance and occurrence_threshold
    parameters to find optimal consensus pharmacophore models that maximize
    separation between active compounds and decoys.

    The optimization uses shape-based scoring with pharmacophore color features
    via RDKit's rdShapeAlign.AlignMol function.

    Attributes:
        n_conformers: Number of conformers to generate per molecule.
        linkage: Clustering linkage method for consensus generation.
        random_state: Random seed for conformer generation.

    Example:
        >>> from pharmacophore.optimization import ConsensusOptimizer
        >>>
        >>> optimizer = ConsensusOptimizer(n_conformers=10)
        >>> results = optimizer.optimize(
        ...     reference_mols=ref_mols,
        ...     actives=active_mols,
        ...     decoys=decoy_mols,
        ...     tolerance_range=(0.5, 4.0),
        ...     occurrence_range=(0.3, 0.9),
        ...     n_points=5
        ... )
        >>> print(f"Best params: {results['best_params']}")
    """

    def __init__(
        self,
        n_conformers: int = 10,
        linkage: str = 'average',
        random_state: int = 42,
        scoring_mode: str = 'reference'
    ):
        """Initialize optimizer.

        Args:
            n_conformers: Number of conformers per molecule (default: 10).
            linkage: Clustering linkage method (default: 'average').
            random_state: Random seed for reproducibility (default: 42).
            scoring_mode: Scoring strategy - 'reference' (recommended) or
                'consensus_mol' (legacy PharmacophoreToMol path).
        """
        self.n_conformers = n_conformers
        self.linkage = linkage
        self.random_state = random_state
        self.scoring_mode = scoring_mode
        self._prepared_refs: List[Chem.Mol] = []
        self._conformer_cache: Dict[str, List[Chem.Mol]] = {}

    def _prepare_references(self, reference_mols: List[Chem.Mol]) -> List[Chem.Mol]:
        """Prepare reference molecules for reference-based scoring."""
        prepared = []
        for mol in reference_mols:
            if mol is None:
                continue
            mol_h = Chem.AddHs(mol)
            if mol_h.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol_h, randomSeed=self.random_state)
            if mol_h.GetNumConformers() > 0:
                prepared.append(mol_h)
        return prepared

    def _score_molecule_reference(
        self,
        mol: Chem.Mol,
        reference_mols: List[Chem.Mol]
    ) -> Tuple[float, float, float]:
        """Score a molecule against reference molecules directly.

        Returns:
            Tuple of (shape_score, color_score, combo_score).
        """
        conformers = self._get_conformers(mol)
        if not conformers:
            return 0.0, 0.0, 0.0

        best_combo = 0.0
        best_shape = 0.0
        best_color = 0.0

        for ref in reference_mols:
            for conf_mol in conformers:
                try:
                    shape, color = AlignMol(
                        ref=ref,
                        probe=conf_mol,
                        useColors=True,
                        opt_param=0.5
                    )
                    combo = shape + color
                    if combo > best_combo:
                        best_combo = combo
                        best_shape = shape
                        best_color = color
                except Exception:
                    continue

        return best_shape, best_color, best_combo

    def score_molecule(
        self,
        mol: Chem.Mol,
        pharmacophore_mol: Chem.Mol,
        use_colors: bool = True
    ) -> Tuple[float, float, float]:
        """Score a molecule against a pharmacophore using shape alignment.

        Generates multiple conformers and returns the best score across all
        conformers. Uses rdShapeAlign.AlignMol with color features enabled.

        Args:
            mol: Query molecule to score.
            pharmacophore_mol: Reference pharmacophore as RDKit Mol.
            use_colors: Enable pharmacophore color scoring (default: True).

        Returns:
            Tuple of (shape_score, color_score, combo_score).
            Combo score = shape + color, range 0-2.

        Raises:
            ValueError: If molecule cannot be processed.
        """
        # Get or generate conformers
        conformers = self._get_conformers(mol)

        if not conformers:
            return 0.0, 0.0, 0.0

        best_combo = 0.0
        best_shape = 0.0
        best_color = 0.0

        for conf_mol in conformers:
            try:
                shape, color = AlignMol(
                    ref=pharmacophore_mol,
                    probe=conf_mol,
                    useColors=use_colors,
                    opt_param=0.5
                )
                combo = shape + color

                if combo > best_combo:
                    best_combo = combo
                    best_shape = shape
                    best_color = color
            except Exception:
                continue

        return best_shape, best_color, best_combo

    def _get_conformers(self, mol: Chem.Mol) -> List[Chem.Mol]:
        """Get or generate conformers for a molecule.

        Uses caching to avoid regenerating conformers for the same molecule.

        Args:
            mol: Input molecule.

        Returns:
            List of conformer molecules.
        """
        # Use SMILES as cache key
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            smiles = None

        if smiles and smiles in self._conformer_cache:
            return self._conformer_cache[smiles]

        # Generate conformers
        conformers = []
        try:
            mol_h = Chem.AddHs(mol)

            # Generate multiple conformers
            conf_ids = AllChem.EmbedMultipleConfs(
                mol_h,
                numConfs=self.n_conformers,
                randomSeed=self.random_state,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True
            )

            if not conf_ids:
                # Try simpler embedding
                AllChem.EmbedMolecule(mol_h, randomSeed=self.random_state)
                if mol_h.GetNumConformers() > 0:
                    conformers.append(mol_h)
            else:
                # Create separate mol for each conformer
                for conf_id in conf_ids:
                    conf_mol = Chem.Mol(mol_h)
                    # Keep only this conformer
                    conf_mol.RemoveAllConformers()
                    conf_mol.AddConformer(
                        mol_h.GetConformer(conf_id),
                        assignId=True
                    )
                    conformers.append(conf_mol)
        except Exception:
            # If conformer generation fails, try using existing conformer
            if mol.GetNumConformers() > 0:
                conformers.append(mol)

        # Cache results
        if smiles:
            self._conformer_cache[smiles] = conformers

        return conformers

    def evaluate_parameters(
        self,
        tolerance: float,
        occurrence_threshold: float,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol]
    ) -> Dict[str, float]:
        """Evaluate a parameter set using enrichment metrics.

        Generates consensus pharmacophore with given parameters, converts to
        mol, scores all actives and decoys, and calculates metrics.

        Args:
            tolerance: Distance threshold for feature clustering (Angstroms).
            occurrence_threshold: Minimum fraction of molecules for features.
            reference_mols: Aligned reference molecules for consensus.
            actives: Active compounds to score.
            decoys: Decoy compounds to score.

        Returns:
            Dict with keys: 'roc_auc', 'ef_1', 'ef_5', 'tolerance',
            'occurrence_threshold', 'n_features'.
        """
        # Generate consensus pharmacophore
        consensus = PharmacophoreConsensus(
            tolerance=tolerance,
            occurrence_threshold=occurrence_threshold,
            linkage=self.linkage
        )

        try:
            features = consensus.generate_consensus(reference_mols)
        except Exception:
            return self._empty_results(tolerance, occurrence_threshold)

        if not features:
            return self._empty_results(tolerance, occurrence_threshold)

        if self.scoring_mode == 'reference':
            # Reference-based scoring (recommended)
            if not self._prepared_refs:
                self._prepared_refs = self._prepare_references(reference_mols)

            # Score actives
            active_scores = []
            for mol in actives:
                _, _, combo = self._score_molecule_reference(mol, self._prepared_refs)
                active_scores.append(combo)

            # Score decoys
            decoy_scores = []
            for mol in decoys:
                _, _, combo = self._score_molecule_reference(mol, self._prepared_refs)
                decoy_scores.append(combo)
        else:
            # Legacy PharmacophoreToMol path (deprecated)
            warnings.warn(
                "scoring_mode='consensus_mol' uses anti-discriminative "
                "PharmacophoreToMol scoring. Use scoring_mode='reference' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            try:
                pharm_mol = PharmacophoreToMol.convert(
                    features,
                    name='Consensus',
                    enable_color_features=True
                )
            except Exception:
                return self._empty_results(tolerance, occurrence_threshold)

            # Score actives
            active_scores = []
            for mol in actives:
                _, _, combo = self.score_molecule(mol, pharm_mol)
                active_scores.append(combo)

            # Score decoys
            decoy_scores = []
            for mol in decoys:
                _, _, combo = self.score_molecule(mol, pharm_mol)
                decoy_scores.append(combo)

        # Calculate metrics
        y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
        y_scores = active_scores + decoy_scores

        # ROC-AUC
        try:
            auc = roc_auc_score(y_true, y_scores)
        except Exception:
            auc = 0.5

        # Enrichment factors
        ef_1 = enrichment_factor(y_true, y_scores, percentage=0.01)
        ef_5 = enrichment_factor(y_true, y_scores, percentage=0.05)
        ef_10 = enrichment_factor(y_true, y_scores, percentage=0.10)

        # BEDROC for early recognition
        bedroc_score = bedroc(y_true, y_scores, alpha=20.0)

        # Youden's index for optimal threshold
        thresh, youden_j, sensitivity, specificity = youden_index(y_true, y_scores)

        return {
            'roc_auc': auc,
            'bedroc': bedroc_score,
            'ef_1': ef_1,
            'ef_5': ef_5,
            'ef_10': ef_10,
            'youden_j': youden_j,
            'optimal_threshold': thresh,
            'sensitivity': sensitivity,
            'specificity': specificity,
            'tolerance': tolerance,
            'occurrence_threshold': occurrence_threshold,
            'n_features': len(features)
        }

    def _empty_results(
        self,
        tolerance: float,
        occurrence_threshold: float
    ) -> Dict[str, float]:
        """Return empty results for failed parameter evaluation."""
        return {
            'roc_auc': 0.5,
            'bedroc': 0.0,
            'ef_1': 0.0,
            'ef_5': 0.0,
            'ef_10': 0.0,
            'youden_j': 0.0,
            'optimal_threshold': 0.0,
            'sensitivity': 0.0,
            'specificity': 0.0,
            'tolerance': tolerance,
            'occurrence_threshold': occurrence_threshold,
            'n_features': 0
        }

    def optimize(
        self,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol],
        tolerance_range: Tuple[float, float] = (0.5, 4.0),
        occurrence_range: Tuple[float, float] = (0.3, 0.9),
        n_points: int = 5,
        metric: str = 'roc_auc',
        verbose: bool = True
    ) -> Dict:
        """Optimize parameters using grid search.

        Performs grid search over tolerance and occurrence_threshold
        parameters to find optimal values that maximize the chosen metric.

        Args:
            reference_mols: Aligned reference molecules for consensus.
            actives: Active compounds for validation.
            decoys: Decoy compounds for validation.
            tolerance_range: (min, max) for tolerance (default: 0.5-4.0 A).
            occurrence_range: (min, max) for occurrence (default: 0.3-0.9).
            n_points: Grid points per dimension (default: 5).
            metric: Metric to optimize ('roc_auc', 'ef_1', 'ef_5').
            verbose: Print progress (default: True).

        Returns:
            Dict with keys:
                'best_params': Dict with optimal tolerance and occurrence.
                'best_score': Best metric value achieved.
                'best_results': Full results dict for best params.
                'all_results': List of all evaluation results.
                'grid': Dict with tolerance and occurrence arrays.

        Example:
            >>> results = optimizer.optimize(refs, actives, decoys)
            >>> best = results['best_params']
            >>> print(f"Optimal: tol={best['tolerance']:.2f}, "
            ...       f"occ={best['occurrence_threshold']:.2f}")
        """
        # Generate parameter grid
        tolerances = np.linspace(
            tolerance_range[0],
            tolerance_range[1],
            n_points
        )
        occurrences = np.linspace(
            occurrence_range[0],
            occurrence_range[1],
            n_points
        )

        all_results = []
        best_score = -float('inf')
        best_results = None

        total = len(tolerances) * len(occurrences)
        current = 0

        for tol in tolerances:
            for occ in occurrences:
                current += 1

                if verbose:
                    print(f"[{current}/{total}] tol={tol:.2f}, occ={occ:.2f}",
                          end='')

                results = self.evaluate_parameters(
                    tolerance=tol,
                    occurrence_threshold=occ,
                    reference_mols=reference_mols,
                    actives=actives,
                    decoys=decoys
                )

                all_results.append(results)

                score = results[metric]
                if verbose:
                    print(f" -> {metric}={score:.3f}, "
                          f"BEDROC={results['bedroc']:.3f}, "
                          f"n_features={results['n_features']}")

                if score > best_score:
                    best_score = score
                    best_results = results

        return {
            'best_params': {
                'tolerance': best_results['tolerance'],
                'occurrence_threshold': best_results['occurrence_threshold']
            },
            'best_score': best_score,
            'best_results': best_results,
            'all_results': all_results,
            'grid': {
                'tolerances': tolerances.tolist(),
                'occurrences': occurrences.tolist()
            }
        }

    @staticmethod
    def _enrichment_factor(
        y_true: List[int],
        y_scores: List[float],
        percentage: float = 0.01
    ) -> float:
        """Calculate enrichment factor at a given percentage.

        This is a wrapper for backward compatibility.
        Use pharmacophore.screening_metrics.enrichment_factor for new code.

        Args:
            y_true: Binary labels (1=active, 0=decoy).
            y_scores: Predicted scores (higher = more likely active).
            percentage: Top percentage to consider (default: 0.01 = 1%).

        Returns:
            Enrichment factor value.
        """
        return enrichment_factor(y_true, y_scores, percentage)

    def clear_cache(self):
        """Clear the conformer cache to free memory."""
        self._conformer_cache.clear()


def optimize_consensus(
    reference_mols: List[Chem.Mol],
    actives: List[Chem.Mol],
    decoys: List[Chem.Mol],
    n_conformers: int = 10,
    tolerance_range: Tuple[float, float] = (0.5, 4.0),
    occurrence_range: Tuple[float, float] = (0.3, 0.9),
    n_points: int = 5,
    metric: str = 'roc_auc',
    verbose: bool = True,
    scoring_mode: str = 'reference'
) -> Dict:
    """Convenience function for parameter optimization.

    Creates a ConsensusOptimizer and runs optimization in one call.

    Args:
        reference_mols: Aligned reference molecules for consensus.
        actives: Active compounds for validation.
        decoys: Decoy compounds for validation.
        n_conformers: Conformers per molecule (default: 10).
        tolerance_range: (min, max) for tolerance (default: 0.5-4.0 A).
        occurrence_range: (min, max) for occurrence (default: 0.3-0.9).
        n_points: Grid points per dimension (default: 5).
        metric: Metric to optimize ('roc_auc', 'ef_1', 'ef_5').
        verbose: Print progress (default: True).
        scoring_mode: Scoring strategy - 'reference' (recommended) or
            'consensus_mol' (legacy). Default: 'reference'.

    Returns:
        Optimization results dict (see ConsensusOptimizer.optimize).

    Example:
        >>> from pharmacophore.optimization import optimize_consensus
        >>> results = optimize_consensus(refs, actives, decoys)
        >>> print(f"Best AUC: {results['best_score']:.3f}")
    """
    optimizer = ConsensusOptimizer(
        n_conformers=n_conformers,
        scoring_mode=scoring_mode
    )
    return optimizer.optimize(
        reference_mols=reference_mols,
        actives=actives,
        decoys=decoys,
        tolerance_range=tolerance_range,
        occurrence_range=occurrence_range,
        n_points=n_points,
        metric=metric,
        verbose=verbose
    )


if __name__ == "__main__":
    import doctest
    doctest.testmod()

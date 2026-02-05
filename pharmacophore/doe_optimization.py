"""Design of Experiments (DOE) optimization for pharmacophore models.

This module provides Bayesian optimization for finding optimal consensus
pharmacophore parameters. It uses Gaussian Process surrogate models with
Expected Improvement acquisition to efficiently explore the parameter space.

Key advantages over grid search:
- Typically finds optimum with 50% fewer evaluations
- Builds response surface model for parameter sensitivity insights
- Handles continuous parameter spaces naturally

Example:
    >>> from pharmacophore.doe_optimization import PharmacophoreOptimizer
    >>> optimizer = PharmacophoreOptimizer(ref_mols, actives, decoys)
    >>> result = optimizer.optimize(n_calls=50)
    >>> print(f"Best AUC: {result['best_auc']:.4f}")
    >>> print(f"Best params: {result['best_params']}")
"""

from typing import List, Dict, Optional, Tuple, Any, Callable
import warnings
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol
from sklearn.metrics import roc_auc_score

try:
    from skopt import gp_minimize
    from skopt.space import Real
    from skopt.plots import plot_convergence, plot_objective
    HAS_SKOPT = True
except ImportError:
    HAS_SKOPT = False

from .consensus import PharmacophoreConsensus
from .mol_converter import PharmacophoreToMol


class PharmacophoreOptimizer:
    """Bayesian optimization for consensus pharmacophore parameters.

    Uses Gaussian Process with Expected Improvement to find optimal
    tolerance, occurrence_threshold, and shape_weight parameters.

    Attributes:
        reference_mols: List of aligned reference molecules.
        actives: List of active molecules with conformers.
        decoys: List of decoy molecules with conformers.
        history: List of (params, score) tuples from optimization.

    Example:
        >>> optimizer = PharmacophoreOptimizer(refs, actives, decoys)
        >>> result = optimizer.optimize(n_calls=50)
        >>> optimizer.plot_convergence()  # Visualize optimization
    """

    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol],
        n_conformers: int = 3,
        random_state: int = 42,
        scoring_mode: str = 'reference'
    ):
        """Initialize optimizer with molecules.

        Args:
            reference_mols: Aligned reference molecules for consensus.
            actives: Active compounds (SMILES-based, will generate conformers).
            decoys: Decoy compounds (SMILES-based, will generate conformers).
            n_conformers: Number of conformers per query molecule.
            random_state: Random seed for reproducibility.
        """
        if not HAS_SKOPT:
            raise ImportError(
                "scikit-optimize required. Install with: pip install scikit-optimize"
            )

        self.reference_mols = reference_mols
        self.n_conformers = n_conformers
        self.random_state = random_state
        self.scoring_mode = scoring_mode
        self._prepared_refs = []
        self.history: List[Dict[str, Any]] = []
        self._result = None

        # Prepare molecules with conformers
        self.actives = self._prepare_molecules(actives)
        self.decoys = self._prepare_molecules(decoys)

        # Ground truth labels
        self.y_true = np.array(
            [1] * len(self.actives) + [0] * len(self.decoys)
        )

    def _prepare_molecules(self, mols: List[Chem.Mol]) -> List[Chem.Mol]:
        """Generate conformers for molecules if needed."""
        prepared = []
        for mol in mols:
            if mol is None:
                continue

            # Check if already has conformers
            if mol.GetNumConformers() >= self.n_conformers:
                prepared.append(mol)
                continue

            # Generate conformers
            mol_h = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = self.random_state
            AllChem.EmbedMultipleConfs(
                mol_h, numConfs=self.n_conformers, params=params
            )

            if mol_h.GetNumConformers() > 0:
                prepared.append(mol_h)

        return prepared

    def _score_molecule(
        self,
        mol: Chem.Mol,
        pharm_mol: Chem.Mol,
        shape_weight: float
    ) -> float:
        """Score a single molecule against pharmacophore."""
        if mol is None or pharm_mol is None:
            return 0.0

        color_weight = 1.0 - shape_weight
        best_score = 0.0

        for conf_id in range(mol.GetNumConformers()):
            try:
                shape, color = AlignMol(
                    ref=pharm_mol,
                    probe=mol,
                    probeConfId=conf_id,
                    useColors=True,
                    opt_param=0.5
                )
                score = shape_weight * shape + color_weight * color
                best_score = max(best_score, score)
            except Exception:
                continue

        return best_score

    def _score_molecule_reference(
        self,
        mol: Chem.Mol,
        shape_weight: float
    ) -> float:
        """Score a single molecule against reference ensemble."""
        if mol is None:
            return 0.0

        color_weight = 1.0 - shape_weight
        best_score = 0.0

        for ref in self._prepared_refs:
            for conf_id in range(mol.GetNumConformers()):
                try:
                    shape, color = AlignMol(
                        ref=ref,
                        probe=mol,
                        probeConfId=conf_id,
                        useColors=True,
                        opt_param=0.5
                    )
                    score = shape_weight * shape + color_weight * color
                    best_score = max(best_score, score)
                except Exception:
                    continue

        return best_score

    def evaluate(
        self,
        tolerance: float,
        occurrence: float,
        shape_weight: float
    ) -> float:
        """Evaluate a parameter configuration.

        Args:
            tolerance: Consensus clustering tolerance (Angstroms).
            occurrence: Minimum feature occurrence fraction.
            shape_weight: Weight for shape vs color (0-1).

        Returns:
            ROC-AUC score for this configuration.
        """
        # Generate consensus
        consensus = PharmacophoreConsensus(
            tolerance=tolerance,
            occurrence_threshold=occurrence
        )
        features = consensus.generate_consensus(self.reference_mols)

        if len(features) < 2:
            return 0.5  # Random performance

        if self.scoring_mode == 'reference':
            # Reference-based scoring (recommended)
            if not self._prepared_refs:
                self._prepared_refs = []
                for mol in self.reference_mols:
                    if mol is None:
                        continue
                    mol_h = Chem.AddHs(mol)
                    if mol_h.GetNumConformers() == 0:
                        AllChem.EmbedMolecule(mol_h, randomSeed=self.random_state)
                    if mol_h.GetNumConformers() > 0:
                        self._prepared_refs.append(mol_h)

            active_scores = [
                self._score_molecule_reference(m, shape_weight)
                for m in self.actives
            ]
            decoy_scores = [
                self._score_molecule_reference(m, shape_weight)
                for m in self.decoys
            ]
        else:
            # Legacy PharmacophoreToMol path (deprecated)
            warnings.warn(
                "scoring_mode='consensus_mol' uses anti-discriminative "
                "PharmacophoreToMol scoring. Use scoring_mode='reference' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            pharm_mol = PharmacophoreToMol.convert(
                features,
                name='Consensus',
                enable_color_features=True
            )

            if pharm_mol is None:
                return 0.5

            active_scores = [
                self._score_molecule(m, pharm_mol, shape_weight)
                for m in self.actives
            ]
            decoy_scores = [
                self._score_molecule(m, pharm_mol, shape_weight)
                for m in self.decoys
            ]

        y_scores = np.array(active_scores + decoy_scores)

        # Calculate AUC
        try:
            auc = roc_auc_score(self.y_true, y_scores)
        except ValueError:
            auc = 0.5

        # Store in history
        self.history.append({
            'tolerance': tolerance,
            'occurrence': occurrence,
            'shape_weight': shape_weight,
            'auc': auc,
            'n_features': len(features)
        })

        return auc

    def optimize(
        self,
        n_calls: int = 50,
        n_random_starts: int = 10,
        tolerance_range: Tuple[float, float] = (0.5, 4.0),
        occurrence_range: Tuple[float, float] = (0.1, 1.0),
        shape_weight_range: Tuple[float, float] = (0.3, 0.9),
        verbose: bool = True
    ) -> Dict[str, Any]:
        """Run Bayesian optimization to find optimal parameters.

        Args:
            n_calls: Total number of evaluations (budget).
            n_random_starts: Initial random exploration points.
            tolerance_range: Search range for tolerance parameter.
            occurrence_range: Search range for occurrence parameter.
            shape_weight_range: Search range for shape_weight parameter.
            verbose: Print progress during optimization.

        Returns:
            Dictionary with:
            - best_params: Optimal parameter values
            - best_auc: AUC at optimal parameters
            - history: Full optimization history
            - result: scikit-optimize result object
        """
        # Define parameter space
        space = [
            Real(tolerance_range[0], tolerance_range[1], name='tolerance'),
            Real(occurrence_range[0], occurrence_range[1], name='occurrence'),
            Real(shape_weight_range[0], shape_weight_range[1], name='shape_weight'),
        ]

        # Objective function (minimize negative AUC)
        def objective(params):
            tol, occ, sw = params
            auc = self.evaluate(tol, occ, sw)
            if verbose:
                print(
                    f"  Eval {len(self.history):3d}: "
                    f"tol={tol:.2f}, occ={occ:.2f}, sw={sw:.2f} → AUC={auc:.4f}"
                )
            return -auc  # Minimize negative = maximize positive

        if verbose:
            print(f"Starting Bayesian Optimization ({n_calls} evaluations)")
            print(f"Parameter space:")
            print(f"  tolerance: {tolerance_range}")
            print(f"  occurrence: {occurrence_range}")
            print(f"  shape_weight: {shape_weight_range}")
            print()

        # Run optimization
        self._result = gp_minimize(
            objective,
            space,
            n_calls=n_calls,
            n_initial_points=n_random_starts,
            acq_func='EI',  # Expected Improvement
            random_state=self.random_state,
            verbose=False
        )

        # Extract results
        best_params = {
            'tolerance': self._result.x[0],
            'occurrence': self._result.x[1],
            'shape_weight': self._result.x[2]
        }
        best_auc = -self._result.fun

        if verbose:
            print(f"\n{'='*50}")
            print(f"OPTIMIZATION COMPLETE")
            print(f"{'='*50}")
            print(f"Best AUC: {best_auc:.4f}")
            print(f"Best parameters:")
            print(f"  tolerance: {best_params['tolerance']:.3f} Å")
            print(f"  occurrence: {best_params['occurrence']:.3f}")
            print(f"  shape_weight: {best_params['shape_weight']:.3f}")

        return {
            'best_params': best_params,
            'best_auc': best_auc,
            'history': self.history,
            'result': self._result
        }

    def plot_convergence(self, save_path: Optional[str] = None):
        """Plot optimization convergence curve.

        Args:
            save_path: If provided, save figure to this path.
        """
        if self._result is None:
            raise ValueError("Run optimize() first")

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))
        plot_convergence(self._result, ax=ax)
        ax.set_ylabel('Negative AUC (minimize)')
        ax.set_title('Bayesian Optimization Convergence')

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')

        return fig

    def plot_response_surface(self, save_path: Optional[str] = None):
        """Plot response surface showing parameter effects.

        Args:
            save_path: If provided, save figure to this path.
        """
        if self._result is None:
            raise ValueError("Run optimize() first")

        import matplotlib.pyplot as plt

        fig = plot_objective(
            self._result,
            dimensions=['tolerance', 'occurrence', 'shape_weight'],
            n_points=20
        )
        fig.suptitle('Parameter Response Surface', y=1.02)

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')

        return fig

    def get_parameter_importance(self) -> Dict[str, float]:
        """Estimate parameter importance from optimization history.

        Returns:
            Dictionary mapping parameter names to importance scores.
        """
        if len(self.history) < 10:
            raise ValueError("Need at least 10 evaluations for importance analysis")

        import pandas as pd
        from sklearn.ensemble import RandomForestRegressor

        # Convert history to DataFrame
        df = pd.DataFrame(self.history)
        X = df[['tolerance', 'occurrence', 'shape_weight']].values
        y = df['auc'].values

        # Fit random forest for feature importance
        rf = RandomForestRegressor(n_estimators=100, random_state=self.random_state)
        rf.fit(X, y)

        importance = {
            'tolerance': rf.feature_importances_[0],
            'occurrence': rf.feature_importances_[1],
            'shape_weight': rf.feature_importances_[2]
        }

        return importance

    def __repr__(self) -> str:
        return (
            f"PharmacophoreOptimizer("
            f"n_refs={len(self.reference_mols)}, "
            f"n_actives={len(self.actives)}, "
            f"n_decoys={len(self.decoys)}, "
            f"evaluations={len(self.history)})"
        )


def optimize_pharmacophore(
    reference_mols: List[Chem.Mol],
    actives: List[Chem.Mol],
    decoys: List[Chem.Mol],
    n_calls: int = 50,
    n_conformers: int = 3,
    verbose: bool = True
) -> Dict[str, Any]:
    """Convenience function for pharmacophore parameter optimization.

    Args:
        reference_mols: Aligned reference molecules.
        actives: Active compounds.
        decoys: Decoy compounds.
        n_calls: Number of optimization evaluations.
        n_conformers: Conformers per query molecule.
        verbose: Print progress.

    Returns:
        Optimization results with best_params and best_auc.

    Example:
        >>> result = optimize_pharmacophore(refs, actives, decoys)
        >>> print(f"Optimal tolerance: {result['best_params']['tolerance']:.2f}")
    """
    optimizer = PharmacophoreOptimizer(
        reference_mols, actives, decoys,
        n_conformers=n_conformers
    )
    return optimizer.optimize(n_calls=n_calls, verbose=verbose)

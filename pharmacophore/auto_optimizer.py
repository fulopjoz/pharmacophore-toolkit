"""End-to-end automated pharmacophore model optimization.

This module provides a complete pipeline from input files (reference structures,
actives, decoys) to an optimal pharmacophore model with minimal manual intervention.
Uses Bayesian optimization with multi-level caching for efficient parameter search.

Key features:
- Automatic data loading from SDF and CSV files
- Multiple scoring methods: Pharm2D (best), Pharm3D, Shape-based
- Consensus feature caching for accelerated evaluation
- Early stopping to save computational budget
- Model export to PyMOL (.pml) and JSON formats

Example:
    >>> from pharmacophore import AutoPharmacophoreOptimizer
    >>> optimizer = AutoPharmacophoreOptimizer(n_conformers=5)
    >>> optimizer.load_from_files(
    ...     reference_file='refs.sdf',
    ...     actives_file='actives.csv',
    ...     decoys_file='decoys.csv'
    ... )
    >>> result = optimizer.optimize(n_calls=50, scoring_method='pharm2d')
    >>> optimizer.export_model('output/')
"""

import json
import os
import warnings
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any, Union
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol
from sklearn.metrics import roc_auc_score

try:
    from skopt import gp_minimize
    from skopt.space import Real
    from skopt.callbacks import DeltaYStopper
    from skopt.plots import plot_convergence, plot_objective
    HAS_SKOPT = True
except ImportError:
    HAS_SKOPT = False

from .consensus import PharmacophoreConsensus
from .mol_converter import PharmacophoreToMol
from .caching import (
    PharmacophoreCache,
    ConsensusCacheKey,
    LRUCache,
    hash_features,
    hash_reference_mols,
    create_cache_key,
)
from .pharm2d_scoring import Pharm2DScorer
from .pharm3d_scoring import Pharm3DScorer
from .screening_metrics import (
    bedroc,
    enrichment_factor,
    youden_index,
    calculate_all_metrics,
)


class AutoPharmacophoreOptimizer:
    """End-to-end automated pharmacophore model optimization.

    Provides a complete pipeline from input files to optimal pharmacophore model:
    1. Load reference structures, actives, and decoys from files or objects
    2. Preprocess molecules (conformer generation with ETKDGv3)
    3. Run Bayesian optimization to find optimal parameters
    4. Export final model with visualization

    Attributes:
        n_conformers: Number of conformers per query molecule.
        random_state: Random seed for reproducibility.
        cache: PharmacophoreCache for accelerated evaluation.
        reference_mols: Aligned reference molecules.
        actives: Active molecules with conformers.
        decoys: Decoy molecules with conformers.
        history: List of evaluation results.
        best_result: Best optimization result.

    Example:
        >>> optimizer = AutoPharmacophoreOptimizer(n_conformers=5)
        >>> optimizer.load_from_files(
        ...     reference_file='tutorials/data/CCR2_reference_ligands.sdf',
        ...     actives_file='tutorials/data/actives_ccr2_N75.csv',
        ...     decoys_file='tutorials/data/decoys_ccr2_N500.csv'
        ... )
        >>> result = optimizer.optimize(n_calls=50, scoring_method='pharm2d')
        >>> print(f"Best AUC: {result['best_auc']:.4f}")
    """

    def __init__(
        self,
        n_conformers: int = 5,
        random_state: int = 42,
        cache_enabled: bool = True,
        cache_max_size: int = 100
    ):
        """Initialize auto-optimizer.

        Args:
            n_conformers: Number of conformers per query molecule (default: 5).
            random_state: Random seed for reproducibility (default: 42).
            cache_enabled: Enable consensus/pharm_mol caching (default: True).
            cache_max_size: Maximum cache entries per level (default: 100).
        """
        if not HAS_SKOPT:
            raise ImportError(
                "scikit-optimize required. Install with: pip install scikit-optimize"
            )

        self.n_conformers = n_conformers
        self.random_state = random_state
        self.cache_enabled = cache_enabled

        # Initialize caches
        self.cache = PharmacophoreCache(max_size=cache_max_size) if cache_enabled else None
        self._conformer_cache = LRUCache(max_size=200)

        # Data placeholders
        self.reference_mols: Optional[List[Chem.Mol]] = None
        self.actives: Optional[List[Chem.Mol]] = None
        self.decoys: Optional[List[Chem.Mol]] = None
        self._ref_hash: Optional[str] = None

        # Prepare reference molecules for reference-based scoring
        self._prepared_refs = []

        # Results
        self.history: List[Dict[str, Any]] = []
        self.best_result: Optional[Dict[str, Any]] = None
        self._skopt_result = None

        # Scoring state
        self._scoring_method: str = 'pharm2d'
        self._pharm2d_scorer: Optional[Pharm2DScorer] = None

    # -------------------------------------------------------------------------
    # Data Loading
    # -------------------------------------------------------------------------

    def load_from_files(
        self,
        reference_file: str,
        actives_file: str,
        decoys_file: str,
        smiles_column: Optional[str] = None
    ) -> 'AutoPharmacophoreOptimizer':
        """Load molecules from files.

        Args:
            reference_file: Path to SDF file with aligned reference molecules.
            actives_file: Path to CSV file with active SMILES.
            decoys_file: Path to CSV file with decoy SMILES.
            smiles_column: Column name for SMILES (auto-detected if None).

        Returns:
            self for method chaining.

        Raises:
            FileNotFoundError: If any file doesn't exist.
            ValueError: If files cannot be parsed.
        """
        # Validate file paths
        reference_path = Path(reference_file).resolve()
        actives_path = Path(actives_file).resolve()
        decoys_path = Path(decoys_file).resolve()
        for p, label in [(reference_path, "Reference"), (actives_path, "Actives"), (decoys_path, "Decoys")]:
            if not p.exists():
                raise FileNotFoundError(f"{label} file not found: {p}")

        # Load references from SDF
        self.reference_mols = self._load_sdf(str(reference_path))
        if not self.reference_mols:
            raise ValueError(f"No molecules loaded from {reference_path}")

        # Load actives and decoys from CSV
        self.actives = self._load_csv_smiles(str(actives_path), smiles_column)
        if not self.actives:
            raise ValueError(f"No molecules loaded from {actives_path}")

        self.decoys = self._load_csv_smiles(str(decoys_path), smiles_column)
        if not self.decoys:
            raise ValueError(f"No molecules loaded from {decoys_path}")

        # Compute reference hash for caching
        self._ref_hash = hash_reference_mols(self.reference_mols)

        # Initialize Pharm2D scorer with references
        self._pharm2d_scorer = Pharm2DScorer(self.reference_mols)

        return self

    def load_from_objects(
        self,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol]
    ) -> 'AutoPharmacophoreOptimizer':
        """Load molecules from Python objects.

        Args:
            reference_mols: Aligned reference molecules (must have 3D conformers).
            actives: Active compounds.
            decoys: Decoy compounds.

        Returns:
            self for method chaining.

        Raises:
            ValueError: If any list is empty or contains None.
        """
        if not reference_mols:
            raise ValueError("reference_mols cannot be empty")
        if not actives:
            raise ValueError("actives cannot be empty")
        if not decoys:
            raise ValueError("decoys cannot be empty")

        self.reference_mols = [m for m in reference_mols if m is not None]
        self.actives = [m for m in actives if m is not None]
        self.decoys = [m for m in decoys if m is not None]

        # Compute reference hash for caching
        self._ref_hash = hash_reference_mols(self.reference_mols)

        # Initialize Pharm2D scorer with references
        self._pharm2d_scorer = Pharm2DScorer(self.reference_mols)

        return self

    def _load_sdf(self, filepath: str) -> List[Chem.Mol]:
        """Load molecules from SDF file."""
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"SDF file not found: {filepath}")

        suppl = Chem.SDMolSupplier(filepath, removeHs=False)
        mols = [mol for mol in suppl if mol is not None]
        return mols

    def _load_csv_smiles(
        self,
        filepath: str,
        smiles_column: Optional[str] = None
    ) -> List[Chem.Mol]:
        """Load molecules from CSV file with SMILES.

        Auto-detects SMILES column if not specified.
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"CSV file not found: {filepath}")

        df = pd.read_csv(filepath)

        # Auto-detect SMILES column
        if smiles_column is None:
            candidates = ['smiles', 'SMILES', 'Smiles', 'smi', 'canonical_smiles']
            for col in candidates:
                if col in df.columns:
                    smiles_column = col
                    break

            if smiles_column is None:
                # Try first column
                smiles_column = df.columns[0]

        if smiles_column not in df.columns:
            raise ValueError(f"SMILES column '{smiles_column}' not found in {filepath}")

        mols = []
        for smiles in df[smiles_column]:
            try:
                mol = Chem.MolFromSmiles(str(smiles))
                if mol is not None:
                    mols.append(mol)
            except Exception:
                continue

        return mols

    # -------------------------------------------------------------------------
    # Conformer Generation
    # -------------------------------------------------------------------------

    def _get_conformers(self, mol: Chem.Mol) -> List[Chem.Mol]:
        """Get or generate conformers for a molecule (cached)."""
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            smiles = None

        if smiles:
            cached = self._conformer_cache.get(smiles)
            if cached is not None:
                return cached

        conformers = []
        try:
            mol_h = Chem.AddHs(mol)

            # Use ETKDGv3 for conformer generation (no force field needed)
            params = AllChem.ETKDGv3()
            params.randomSeed = self.random_state

            conf_ids = AllChem.EmbedMultipleConfs(
                mol_h, numConfs=self.n_conformers, params=params
            )

            if not conf_ids:
                # Fallback: try single conformer
                AllChem.EmbedMolecule(mol_h, params)
                if mol_h.GetNumConformers() > 0:
                    conformers.append(mol_h)
            else:
                for conf_id in conf_ids:
                    conf_mol = Chem.Mol(mol_h)
                    conf_mol.RemoveAllConformers()
                    conf_mol.AddConformer(mol_h.GetConformer(conf_id), assignId=True)
                    conformers.append(conf_mol)

        except Exception:
            if mol.GetNumConformers() > 0:
                conformers.append(mol)

        if smiles:
            self._conformer_cache.set(smiles, conformers)

        return conformers

    # -------------------------------------------------------------------------
    # Scoring Methods
    # -------------------------------------------------------------------------

    def _score_pharm2d(self, features: List) -> Tuple[float, Dict[str, float]]:
        """Score using Pharm2D fingerprints (best method).

        Returns:
            Tuple of (roc_auc, metrics_dict).
        """
        if self._pharm2d_scorer is None:
            raise ValueError("Pharm2D scorer not initialized. Load data first.")

        # Score all molecules
        active_scores = self._pharm2d_scorer.score_all(self.actives)
        decoy_scores = self._pharm2d_scorer.score_all(self.decoys)

        y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
        y_scores = list(active_scores) + list(decoy_scores)

        # Calculate metrics
        try:
            auc = roc_auc_score(y_true, y_scores)
        except ValueError:
            auc = 0.5

        metrics = {
            'roc_auc': auc,
            'bedroc': bedroc(y_true, y_scores, alpha=20.0),
            'ef_1': enrichment_factor(y_true, y_scores, percentage=0.01),
            'ef_5': enrichment_factor(y_true, y_scores, percentage=0.05),
        }

        return auc, metrics

    def _score_pharm3d(self, features: List) -> Tuple[float, Dict[str, float]]:
        """Score using Pharm3D distance-based method.

        Returns:
            Tuple of (roc_auc, metrics_dict).
        """
        scorer = Pharm3DScorer(features, distance_tolerance=1.5)

        active_scores = []
        for mol in self.actives:
            score = scorer.score(mol)
            active_scores.append(score)

        decoy_scores = []
        for mol in self.decoys:
            score = scorer.score(mol)
            decoy_scores.append(score)

        y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
        y_scores = active_scores + decoy_scores

        try:
            auc = roc_auc_score(y_true, y_scores)
        except ValueError:
            auc = 0.5

        metrics = {
            'roc_auc': auc,
            'bedroc': bedroc(y_true, y_scores, alpha=20.0),
            'ef_1': enrichment_factor(y_true, y_scores, percentage=0.01),
            'ef_5': enrichment_factor(y_true, y_scores, percentage=0.05),
        }

        return auc, metrics

    def _score_shape(
        self,
        features: List,
        shape_weight: float
    ) -> Tuple[float, Dict[str, float]]:
        """Score using shape alignment with pharmacophore colors.

        Returns:
            Tuple of (roc_auc, metrics_dict).
        """
        # Convert features to pharmacophore mol
        feature_hash = hash_features(features)

        if self.cache_enabled:
            pharm_mol = self.cache.get_pharm_mol(feature_hash)
        else:
            pharm_mol = None

        if pharm_mol is None:
            try:
                pharm_mol = PharmacophoreToMol.convert(
                    features, name='Consensus', enable_color_features=True
                )
            except Exception:
                return 0.5, {'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0, 'ef_5': 0.0}

            if self.cache_enabled and pharm_mol is not None:
                self.cache.set_pharm_mol(feature_hash, pharm_mol)

        if pharm_mol is None:
            return 0.5, {'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0, 'ef_5': 0.0}

        # Score molecules
        color_weight = 1.0 - shape_weight

        active_scores = []
        for mol in self.actives:
            score = self._score_molecule_shape(mol, pharm_mol, shape_weight, color_weight)
            active_scores.append(score)

        decoy_scores = []
        for mol in self.decoys:
            score = self._score_molecule_shape(mol, pharm_mol, shape_weight, color_weight)
            decoy_scores.append(score)

        y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
        y_scores = active_scores + decoy_scores

        try:
            auc = roc_auc_score(y_true, y_scores)
        except ValueError:
            auc = 0.5

        metrics = {
            'roc_auc': auc,
            'bedroc': bedroc(y_true, y_scores, alpha=20.0),
            'ef_1': enrichment_factor(y_true, y_scores, percentage=0.01),
            'ef_5': enrichment_factor(y_true, y_scores, percentage=0.05),
        }

        return auc, metrics

    def _score_molecule_shape(
        self,
        mol: Chem.Mol,
        pharm_mol: Chem.Mol,
        shape_weight: float,
        color_weight: float
    ) -> float:
        """Score single molecule using shape alignment."""
        conformers = self._get_conformers(mol)
        if not conformers:
            return 0.0

        best_score = 0.0
        for conf_mol in conformers:
            try:
                shape, color = AlignMol(
                    ref=pharm_mol,
                    probe=conf_mol,
                    useColors=True,
                    opt_param=0.5
                )
                score = shape_weight * shape + color_weight * color
                best_score = max(best_score, score)
            except Exception:
                continue

        return best_score

    def _score_reference(
        self,
        features: List,
        shape_weight: float
    ) -> Tuple[float, Dict[str, float]]:
        """Score using reference-based alignment (bypasses PharmacophoreToMol).

        Aligns query molecules to reference ligands directly, avoiding
        the anti-discriminative PharmacophoreToMol conversion.

        Returns:
            Tuple of (roc_auc, metrics_dict).
        """
        # Prepare references lazily
        if not hasattr(self, '_prepared_refs') or not self._prepared_refs:
            self._prepared_refs = []
            for mol in self.reference_mols:
                if mol is None:
                    continue
                mol_h = Chem.AddHs(mol)
                if mol_h.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol_h, randomSeed=self.random_state)
                if mol_h.GetNumConformers() > 0:
                    self._prepared_refs.append(mol_h)

        if not self._prepared_refs:
            return 0.5, {'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0, 'ef_5': 0.0}

        # Score molecules against reference ensemble
        active_scores = []
        for mol in self.actives:
            score = self._score_molecule_against_refs(mol, shape_weight)
            active_scores.append(score)

        decoy_scores = []
        for mol in self.decoys:
            score = self._score_molecule_against_refs(mol, shape_weight)
            decoy_scores.append(score)

        y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
        y_scores = active_scores + decoy_scores

        try:
            auc = roc_auc_score(y_true, y_scores)
        except ValueError:
            auc = 0.5

        metrics = {
            'roc_auc': auc,
            'bedroc': bedroc(y_true, y_scores, alpha=20.0),
            'ef_1': enrichment_factor(y_true, y_scores, percentage=0.01),
            'ef_5': enrichment_factor(y_true, y_scores, percentage=0.05),
        }

        return auc, metrics

    def _score_molecule_against_refs(
        self,
        mol: Chem.Mol,
        shape_weight: float
    ) -> float:
        """Score single molecule against reference ensemble."""
        conformers = self._get_conformers(mol)
        if not conformers:
            return 0.0

        color_weight = 1.0 - shape_weight
        best_score = 0.0

        for ref in self._prepared_refs:
            for conf_mol in conformers:
                try:
                    shape, color = AlignMol(
                        ref=ref,
                        probe=conf_mol,
                        useColors=True,
                        opt_param=0.5
                    )
                    score = shape_weight * shape + color_weight * color
                    best_score = max(best_score, score)
                except Exception:
                    continue

        return best_score

    # -------------------------------------------------------------------------
    # Evaluation
    # -------------------------------------------------------------------------

    def evaluate(
        self,
        tolerance: float,
        occurrence: float,
        shape_weight: float = 0.5
    ) -> float:
        """Evaluate a parameter configuration (with caching).

        Args:
            tolerance: Consensus clustering tolerance (Angstroms).
            occurrence: Minimum feature occurrence fraction.
            shape_weight: Weight for shape vs color (only for 'shape' method).

        Returns:
            ROC-AUC score for this configuration.
        """
        if self.reference_mols is None:
            raise ValueError("No data loaded. Call load_from_files() first.")

        # Check consensus cache
        cache_key = ConsensusCacheKey(
            tolerance=tolerance,
            occurrence=occurrence,
            linkage='average',
            ref_hash=self._ref_hash
        )

        features = None
        if self.cache_enabled:
            features = self.cache.get_consensus(cache_key)

        if features is None:
            # Generate consensus
            consensus = PharmacophoreConsensus(
                tolerance=tolerance,
                occurrence_threshold=occurrence
            )
            try:
                features = consensus.generate_consensus(self.reference_mols)
            except Exception:
                features = []

            # Cache result
            if self.cache_enabled:
                self.cache.set_consensus(cache_key, features)

        if len(features) < 2:
            return 0.5  # Random performance for degenerate cases

        # Score based on method
        if self._scoring_method == 'pharm2d':
            auc, metrics = self._score_pharm2d(features)
        elif self._scoring_method == 'pharm3d':
            auc, metrics = self._score_pharm3d(features)
        elif self._scoring_method == 'reference':
            auc, metrics = self._score_reference(features, shape_weight)
        elif self._scoring_method == 'shape':
            warnings.warn(
                "scoring_method='shape' uses anti-discriminative PharmacophoreToMol "
                "scoring (AUC < 0.5). Use 'reference' or 'pharm2d' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            auc, metrics = self._score_shape(features, shape_weight)
        else:
            auc, metrics = self._score_shape(features, shape_weight)

        # Store in history
        self.history.append({
            'tolerance': tolerance,
            'occurrence': occurrence,
            'shape_weight': shape_weight,
            'n_features': len(features),
            **metrics
        })

        return auc

    # -------------------------------------------------------------------------
    # Optimization
    # -------------------------------------------------------------------------

    def optimize(
        self,
        n_calls: int = 50,
        n_random_starts: int = 10,
        scoring_method: str = 'pharm2d',
        metric: str = 'roc_auc',
        tolerance_range: Tuple[float, float] = (0.5, 4.0),
        occurrence_range: Tuple[float, float] = (0.1, 1.0),
        shape_weight_range: Tuple[float, float] = (0.3, 0.9),
        early_stopping: bool = True,
        early_stopping_delta: float = 0.001,
        early_stopping_n: int = 10,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """Run Bayesian optimization to find optimal parameters.

        Args:
            n_calls: Total number of evaluations (budget).
            n_random_starts: Initial random exploration points.
            scoring_method: 'pharm2d' (best), 'reference' (recommended), 'pharm3d', or 'shape' (deprecated).
            metric: Metric to optimize ('roc_auc', 'bedroc', 'ef_1').
            tolerance_range: Search range for tolerance parameter.
            occurrence_range: Search range for occurrence parameter.
            shape_weight_range: Search range for shape_weight (shape method only).
            early_stopping: Enable early stopping if improvement stalls.
            early_stopping_delta: Minimum improvement threshold.
            early_stopping_n: Consecutive iterations without improvement.
            verbose: Print progress during optimization.

        Returns:
            Dictionary with:
            - best_params: Optimal parameter values
            - best_auc: AUC at optimal parameters
            - history: Full optimization history
            - cache_stats: Cache statistics (if enabled)
        """
        if self.reference_mols is None:
            raise ValueError("No data loaded. Call load_from_files() first.")

        self._scoring_method = scoring_method
        self.history = []  # Reset history

        # Define parameter space
        if scoring_method == 'shape':
            space = [
                Real(tolerance_range[0], tolerance_range[1], name='tolerance'),
                Real(occurrence_range[0], occurrence_range[1], name='occurrence'),
                Real(shape_weight_range[0], shape_weight_range[1], name='shape_weight'),
            ]
        else:
            # Pharm2D and Pharm3D only need tolerance and occurrence
            space = [
                Real(tolerance_range[0], tolerance_range[1], name='tolerance'),
                Real(occurrence_range[0], occurrence_range[1], name='occurrence'),
            ]

        # Objective function
        def objective(params):
            if scoring_method == 'shape':
                tol, occ, sw = params
            else:
                tol, occ = params
                sw = 0.5  # Not used

            auc = self.evaluate(tol, occ, sw)

            if verbose:
                n_eval = len(self.history)
                print(
                    f"  [{n_eval:3d}/{n_calls}] tol={tol:.2f}, occ={occ:.2f}"
                    + (f", sw={sw:.2f}" if scoring_method == 'shape' else "")
                    + f" -> AUC={auc:.4f}"
                )

            return -auc  # Minimize negative = maximize positive

        # Callbacks
        callbacks = []
        if early_stopping:
            callbacks.append(
                DeltaYStopper(delta=early_stopping_delta, n_best=early_stopping_n)
            )

        if verbose:
            print(f"\n{'='*60}")
            print(f"STARTING BAYESIAN OPTIMIZATION")
            print(f"{'='*60}")
            print(f"Scoring method: {scoring_method}")
            print(f"Budget: {n_calls} evaluations")
            print(f"Dataset: {len(self.actives)} actives, {len(self.decoys)} decoys")
            print(f"References: {len(self.reference_mols)} molecules")
            print(f"Parameter space:")
            print(f"  tolerance: {tolerance_range}")
            print(f"  occurrence: {occurrence_range}")
            if scoring_method == 'shape':
                print(f"  shape_weight: {shape_weight_range}")
            print()

        # Run optimization
        self._skopt_result = gp_minimize(
            objective,
            space,
            n_calls=n_calls,
            n_initial_points=n_random_starts,
            acq_func='EI',
            random_state=self.random_state,
            verbose=False,
            callback=callbacks if callbacks else None
        )

        # Extract results
        if scoring_method == 'shape':
            best_params = {
                'tolerance': self._skopt_result.x[0],
                'occurrence': self._skopt_result.x[1],
                'shape_weight': self._skopt_result.x[2]
            }
        else:
            best_params = {
                'tolerance': self._skopt_result.x[0],
                'occurrence': self._skopt_result.x[1]
            }

        best_auc = -self._skopt_result.fun

        self.best_result = {
            'best_params': best_params,
            'best_auc': best_auc,
            'scoring_method': scoring_method,
            'n_evaluations': len(self.history),
            'history': self.history,
            'cache_stats': self.cache.stats() if self.cache_enabled else None
        }

        if verbose:
            print(f"\n{'='*60}")
            print(f"OPTIMIZATION COMPLETE")
            print(f"{'='*60}")
            print(f"Best AUC: {best_auc:.4f}")
            print(f"Best parameters:")
            for k, v in best_params.items():
                print(f"  {k}: {v:.4f}")
            if self.cache_enabled:
                stats = self.cache.stats()
                print(f"\nCache statistics:")
                print(f"  Consensus hit rate: {stats['consensus_hit_rate']:.1%}")
                print(f"  Evaluations: {len(self.history)}/{n_calls}")

        return self.best_result

    # -------------------------------------------------------------------------
    # Export
    # -------------------------------------------------------------------------

    def export_model(
        self,
        output_dir: str,
        include_pml: bool = True,
        include_json: bool = True,
        model_name: str = 'optimal_pharmacophore'
    ) -> Dict[str, str]:
        """Export optimal pharmacophore model.

        Args:
            output_dir: Directory to save outputs.
            include_pml: Generate PyMOL script (.pml).
            include_json: Generate parameters JSON.
            model_name: Base name for output files.

        Returns:
            Dict mapping output type to file path.
        """
        if self.best_result is None:
            raise ValueError("No optimization results. Run optimize() first.")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        outputs = {}

        # Generate final consensus with best parameters
        params = self.best_result['best_params']
        consensus = PharmacophoreConsensus(
            tolerance=params['tolerance'],
            occurrence_threshold=params['occurrence']
        )
        features = consensus.generate_consensus(self.reference_mols)

        # Save JSON
        if include_json:
            json_path = output_dir / f"{model_name}_params.json"
            result_data = {
                'best_params': params,
                'best_auc': self.best_result['best_auc'],
                'scoring_method': self.best_result['scoring_method'],
                'n_features': len(features),
                'features': [
                    {'type': f[0], 'x': f[2], 'y': f[3], 'z': f[4]}
                    for f in features
                ],
                'n_evaluations': self.best_result['n_evaluations'],
                'dataset': {
                    'n_references': len(self.reference_mols),
                    'n_actives': len(self.actives),
                    'n_decoys': len(self.decoys)
                }
            }
            with open(json_path, 'w') as f:
                json.dump(result_data, f, indent=2)
            outputs['json'] = str(json_path)

        # Save PML (PyMOL script)
        if include_pml:
            pml_path = output_dir / f"{model_name}.pml"
            self._generate_pml(features, pml_path, model_name)
            outputs['pml'] = str(pml_path)

        return outputs

    def _generate_pml(self, features: List, filepath: Path, name: str):
        """Generate PyMOL script for pharmacophore visualization."""
        from .constants import FEATURE_COLORS

        lines = [
            f"# PyMOL script for pharmacophore: {name}",
            f"# Generated by AutoPharmacophoreOptimizer",
            f"# Best AUC: {self.best_result['best_auc']:.4f}",
            "",
            "# Create pseudoatoms for pharmacophore features",
        ]

        for i, feat in enumerate(features):
            feat_type = feat[0]
            x, y, z = feat[2], feat[3], feat[4]
            color = FEATURE_COLORS.get(feat_type, 'gray')

            lines.append(
                f"pseudoatom {name}_{feat_type}_{i}, pos=[{x:.3f}, {y:.3f}, {z:.3f}]"
            )
            lines.append(f"color {color}, {name}_{feat_type}_{i}")
            lines.append(f"show spheres, {name}_{feat_type}_{i}")
            lines.append(f"set sphere_scale, 1.0, {name}_{feat_type}_{i}")
            lines.append("")

        # Group all features
        feature_names = [f"{name}_{f[0]}_{i}" for i, f in enumerate(features)]
        lines.append(f"group {name}, {' '.join(feature_names)}")

        with open(filepath, 'w') as f:
            f.write('\n'.join(lines))

    # -------------------------------------------------------------------------
    # Comparison
    # -------------------------------------------------------------------------

    def compare_methods(self) -> pd.DataFrame:
        """Compare all scoring methods on optimal parameters.

        Returns:
            DataFrame with AUC and other metrics for each method.
        """
        if self.best_result is None:
            raise ValueError("No optimization results. Run optimize() first.")

        params = self.best_result['best_params']

        # Generate consensus with optimal parameters
        consensus = PharmacophoreConsensus(
            tolerance=params['tolerance'],
            occurrence_threshold=params['occurrence']
        )
        features = consensus.generate_consensus(self.reference_mols)

        results = []

        # Pharm2D
        auc, metrics = self._score_pharm2d(features)
        results.append({
            'method': 'Pharm2D',
            'roc_auc': auc,
            **metrics
        })

        # Pharm3D
        auc, metrics = self._score_pharm3d(features)
        results.append({
            'method': 'Pharm3D',
            'roc_auc': auc,
            **metrics
        })

        # Shape (use shape_weight from params or default)
        sw = params.get('shape_weight', 0.5)
        auc, metrics = self._score_shape(features, sw)
        results.append({
            'method': f'Shape (sw={sw:.2f})',
            'roc_auc': auc,
            **metrics
        })

        return pd.DataFrame(results).sort_values('roc_auc', ascending=False)

    # -------------------------------------------------------------------------
    # Visualization
    # -------------------------------------------------------------------------

    def plot_convergence(self, save_path: Optional[str] = None):
        """Plot optimization convergence curve.

        Args:
            save_path: If provided, save figure to this path.

        Returns:
            matplotlib Figure object.
        """
        if self._skopt_result is None:
            raise ValueError("Run optimize() first")

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))
        plot_convergence(self._skopt_result, ax=ax)
        ax.set_ylabel('Negative AUC (minimize)')
        ax.set_title('Bayesian Optimization Convergence')

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')

        return fig

    def plot_response_surface(self, save_path: Optional[str] = None):
        """Plot response surface showing parameter effects.

        Args:
            save_path: If provided, save figure to this path.

        Returns:
            matplotlib Figure object.
        """
        if self._skopt_result is None:
            raise ValueError("Run optimize() first")

        import matplotlib.pyplot as plt

        dimensions = ['tolerance', 'occurrence']
        if self._scoring_method == 'shape':
            dimensions.append('shape_weight')

        fig = plot_objective(self._skopt_result, dimensions=dimensions, n_points=20)
        fig.suptitle('Parameter Response Surface', y=1.02)

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')

        return fig

    # -------------------------------------------------------------------------
    # Utilities
    # -------------------------------------------------------------------------

    def clear_cache(self):
        """Clear all caches."""
        if self.cache:
            self.cache.clear()
        self._conformer_cache.clear()

    def __repr__(self) -> str:
        n_refs = len(self.reference_mols) if self.reference_mols else 0
        n_act = len(self.actives) if self.actives else 0
        n_dec = len(self.decoys) if self.decoys else 0
        return (
            f"AutoPharmacophoreOptimizer("
            f"refs={n_refs}, actives={n_act}, decoys={n_dec}, "
            f"evaluations={len(self.history)})"
        )


# -----------------------------------------------------------------------------
# Convenience Function
# -----------------------------------------------------------------------------

def auto_optimize_pharmacophore(
    reference_file: str,
    actives_file: str,
    decoys_file: str,
    n_calls: int = 50,
    n_random_starts: int = 10,
    scoring_method: str = 'pharm2d',
    output_dir: Optional[str] = None,
    verbose: bool = True
) -> Dict[str, Any]:
    """Convenience function for end-to-end pharmacophore optimization.

    Args:
        reference_file: Path to SDF file with aligned reference molecules.
        actives_file: Path to CSV file with active SMILES.
        decoys_file: Path to CSV file with decoy SMILES.
        n_calls: Number of optimization evaluations.
        n_random_starts: Initial random exploration points.
        scoring_method: 'pharm2d' (best), 'pharm3d', or 'shape'.
        output_dir: If provided, export model to this directory.
        verbose: Print progress.

    Returns:
        Optimization results with best_params and best_auc.

    Example:
        >>> result = auto_optimize_pharmacophore(
        ...     'refs.sdf', 'actives.csv', 'decoys.csv',
        ...     n_calls=50, output_dir='output/'
        ... )
        >>> print(f"Best AUC: {result['best_auc']:.4f}")
    """
    optimizer = AutoPharmacophoreOptimizer(n_conformers=5)

    optimizer.load_from_files(
        reference_file=reference_file,
        actives_file=actives_file,
        decoys_file=decoys_file
    )

    # Ensure n_random_starts doesn't exceed n_calls
    n_random_starts = min(n_random_starts, n_calls)

    result = optimizer.optimize(
        n_calls=n_calls,
        n_random_starts=n_random_starts,
        scoring_method=scoring_method,
        verbose=verbose
    )

    if output_dir:
        optimizer.export_model(output_dir)

    return result

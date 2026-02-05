"""Unified evaluation module for pharmacophore optimization.

This module provides a shared evaluation framework used by all three optimization
approaches (Optuna GP, NSGA-II, HypoGen). It ensures fair comparisons by using
identical evaluation logic across all optimizers.

Scoring modes:
    - ``'reference'`` (default): Align queries to reference molecules directly.
      This is the recommended mode as it avoids the anti-discriminative behavior
      of PharmacophoreToMol-based scoring. AUC > 0.85 on CCR2 dataset.
    - ``'pharm2d'``: 2D pharmacophore fingerprint scoring. Fastest mode,
      AUC ~0.91 on CCR2 dataset.
    - ``'hybrid'``: Weighted combination of 2D (pharm2d) + 3D (reference)
      scoring. Best overall discrimination, AUC ~0.95 at alpha=0.7.
    - ``'consensus_mol'`` (deprecated): Legacy PharmacophoreToMol → AlignMol
      path. Anti-discriminative (AUC < 0.5). Kept for backward compatibility.

Example:
    >>> from pharmacophore.evaluation import UnifiedEvaluator, EvaluationConfig
    >>> evaluator = UnifiedEvaluator(ref_mols, actives, decoys)
    >>> config = EvaluationConfig(tolerance=2.0, occurrence=0.5)
    >>> result = evaluator.evaluate(config)
    >>> print(f"AUC: {result.roc_auc:.4f}, BEDROC: {result.bedroc:.4f}")
"""

from dataclasses import dataclass, field
from typing import List, Optional
import logging
import time
import warnings
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol
from joblib import Parallel, delayed

from .consensus import PharmacophoreConsensus
from .mol_converter import PharmacophoreToMol
from .screening_metrics import calculate_all_metrics

logger = logging.getLogger(__name__)

# Reciprocal Rank Fusion constant (Cormack, Clarke & Buettcher, SIGIR 2009).
# Higher k reduces influence of top-ranked items.
RRF_CONSTANT = 60


_VALID_SCORING_MODES = {'reference', 'pharm2d', 'hybrid', 'consensus_mol'}
_VALID_AGGREGATIONS = {'max', 'mean'}


@dataclass
class EvaluationConfig:
    """Configuration for pharmacophore evaluation.

    Attributes:
        tolerance: Consensus clustering spatial tolerance (Angstroms).
        occurrence: Minimum feature occurrence fraction (0.0-1.0).
        shape_weight: Weight for shape vs color scoring (0.0-1.0).
        opt_param: AlignMol optimization parameter (0.0=color-only, 1.0=shape-only).
        linkage: Hierarchical clustering linkage method.
        n_conformers: Number of conformers per query molecule.
        minimize: Whether to MMFF-minimize conformers.
        scoring_mode: Scoring approach to use. One of:
            ``'reference'`` — align to reference molecules (recommended),
            ``'pharm2d'`` — 2D pharmacophore fingerprints,
            ``'hybrid'`` — weighted 2D + 3D combination,
            ``'consensus_mol'`` — legacy PharmacophoreToMol path (deprecated).
        aggregation: How to aggregate scores across references.
            ``'max'`` — best match (default), ``'mean'`` — average.
        alpha: Weight for 2D component in hybrid mode (0.0-1.0).
            Only used when ``scoring_mode='hybrid'``.
    """
    tolerance: float = 2.0
    occurrence: float = 0.5
    shape_weight: float = 0.5
    opt_param: float = 0.5
    linkage: str = 'average'
    n_conformers: int = 25  # Literature-backed optimal (JCIM 2023)
    minimize: bool = False
    scoring_mode: str = 'reference'
    aggregation: str = 'max'
    alpha: float = 0.6

    def __post_init__(self):
        """Validate parameter ranges."""
        if not 0.5 <= self.tolerance <= 4.0:
            raise ValueError(f"tolerance must be in [0.5, 4.0], got {self.tolerance}")
        if not 0.0 <= self.occurrence <= 1.0:
            raise ValueError(f"occurrence must be in [0.0, 1.0], got {self.occurrence}")
        if not 0.0 <= self.shape_weight <= 1.0:
            raise ValueError(f"shape_weight must be in [0.0, 1.0], got {self.shape_weight}")
        if not 0.0 <= self.opt_param <= 1.0:
            raise ValueError(f"opt_param must be in [0.0, 1.0], got {self.opt_param}")
        if self.linkage not in ['average', 'complete', 'single', 'ward']:
            raise ValueError(f"linkage must be average/complete/single/ward, got {self.linkage}")
        if not 5 <= self.n_conformers <= 50:
            raise ValueError(f"n_conformers must be in [5, 50], got {self.n_conformers}")
        if self.scoring_mode not in _VALID_SCORING_MODES:
            raise ValueError(
                f"scoring_mode must be one of {_VALID_SCORING_MODES}, "
                f"got '{self.scoring_mode}'"
            )
        if self.aggregation not in _VALID_AGGREGATIONS:
            raise ValueError(
                f"aggregation must be one of {_VALID_AGGREGATIONS}, "
                f"got '{self.aggregation}'"
            )
        if not 0.0 <= self.alpha <= 1.0:
            raise ValueError(f"alpha must be in [0.0, 1.0], got {self.alpha}")


@dataclass
class EvaluationResult:
    """Results from pharmacophore evaluation.

    Attributes:
        config: Configuration used for this evaluation.
        roc_auc: ROC-AUC score (0.5=random, 1.0=perfect).
        bedroc: BEDROC score for early recognition (α=20).
        ef_1: Enrichment factor at 1%.
        ef_5: Enrichment factor at 5%.
        ef_10: Enrichment factor at 10%.
        n_features: Number of pharmacophore features in consensus.
        eval_time_sec: Wall-clock time for evaluation (seconds).
    """
    config: EvaluationConfig
    roc_auc: float = 0.5
    bedroc: float = 0.0
    ef_1: float = 0.0
    ef_5: float = 0.0
    ef_10: float = 0.0
    n_features: int = 0
    eval_time_sec: float = 0.0

    def summary_string(self) -> str:
        """Generate one-line summary of result."""
        return (
            f"AUC={self.roc_auc:.4f}, BEDROC={self.bedroc:.4f}, "
            f"EF@1%={self.ef_1:.1f}, n_feat={self.n_features}, "
            f"time={self.eval_time_sec:.2f}s"
        )


def compute_sdbw(metadata: dict) -> float:
    """Compute S_Dbw cluster validation index from consensus metadata.

    S_Dbw(c) = Scat(c) + Dens_bw(c), where:
    - Scat: average cluster scatter / total scatter (intra-cluster compactness)
    - Dens_bw: inter-cluster density at midpoints (cluster separation)

    Lower S_Dbw = better clustering quality.

    Source: Braun & Fayne 2022 (MoPBS), Halkidi et al. 2002.

    Args:
        metadata: Dict from PharmacophoreConsensus.generate_consensus(return_metadata=True).
            Keys are feature types, values have 'coordinates' and 'labels'.

    Returns:
        S_Dbw score (lower = better). Returns float('inf') if computation fails.
    """
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
        return float('inf')

    coords = np.vstack(all_coords)
    labels = np.concatenate(all_labels)
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    if n_clusters < 2 or len(coords) < 2:
        return float('inf')

    # Total dataset variance
    total_var = np.var(coords, axis=0)
    total_std = np.sqrt(np.sum(total_var))

    if total_std < 1e-12:
        return float('inf')

    # Scat: average intra-cluster scatter / total scatter
    cluster_stds = []
    centroids = {}
    for label in unique_labels:
        mask = labels == label
        cluster_coords = coords[mask]
        centroids[label] = np.mean(cluster_coords, axis=0)
        cluster_var = np.var(cluster_coords, axis=0)
        cluster_stds.append(np.sqrt(np.sum(cluster_var)))

    scat = np.mean(cluster_stds) / total_std

    # Dens_bw: inter-cluster density at midpoints
    stdev_avg = np.mean(cluster_stds)
    if stdev_avg < 1e-12:
        stdev_avg = 1e-12

    dens_bw_sum = 0.0
    n_pairs = 0
    for i, li in enumerate(unique_labels):
        for j, lj in enumerate(unique_labels):
            if i >= j:
                continue
            midpoint = (centroids[li] + centroids[lj]) / 2.0

            # Count points within stdev_avg of midpoint
            dists = np.linalg.norm(coords - midpoint, axis=1)
            n_mid = np.sum(dists <= stdev_avg)

            # Count points in each cluster within stdev_avg of their centroid
            dists_i = np.linalg.norm(
                coords[labels == li] - centroids[li], axis=1)
            n_i = max(np.sum(dists_i <= stdev_avg), 1)

            dists_j = np.linalg.norm(
                coords[labels == lj] - centroids[lj], axis=1)
            n_j = max(np.sum(dists_j <= stdev_avg), 1)

            dens_bw_sum += n_mid / max(n_i, n_j)
            n_pairs += 1

    dens_bw = dens_bw_sum / max(n_pairs, 1)

    return scat + dens_bw


class UnifiedEvaluator:
    """Unified evaluation framework for all optimization approaches.

    Supports multiple scoring modes via ``EvaluationConfig.scoring_mode``:

    - ``'reference'`` (default, recommended): Aligns query molecules to the
      actual reference ligands, bypassing PharmacophoreToMol. AUC > 0.85.
    - ``'pharm2d'``: Uses Pharm2D fingerprint Tanimoto similarity. AUC ~0.91.
    - ``'hybrid'``: Weighted 2D + 3D combination. AUC ~0.95 at alpha=0.7.
    - ``'consensus_mol'`` (deprecated): Legacy PharmacophoreToMol → AlignMol.
      Anti-discriminative (AUC < 0.5).

    Note:
        ``evaluate()`` is not thread-safe for concurrent calls. Use separate
        instances for concurrent use. Parallel scoring within a single call
        is safe (loky uses separate processes).

    Attributes:
        reference_mols: Aligned reference molecules for consensus generation.
        actives: Active compounds with conformers.
        decoys: Decoy compounds with conformers.
        y_true: Ground truth labels (1=active, 0=decoy).
        random_state: Random seed for reproducibility.

    Example:
        >>> evaluator = UnifiedEvaluator(ref_mols, actives, decoys)
        >>> config = EvaluationConfig(tolerance=2.0, occurrence=0.3)
        >>> result = evaluator.evaluate(config)
    """

    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol],
        random_state: int = 42,
        n_jobs: int = 1
    ):
        """Initialize evaluator with molecules.

        Args:
            reference_mols: Aligned reference molecules for consensus.
            actives: Active compounds (will generate conformers if needed).
            decoys: Decoy compounds (will generate conformers if needed).
            random_state: Random seed for conformer generation.
            n_jobs: Number of parallel jobs for scoring (-1 = all cores).
        """
        self.reference_mols = reference_mols
        self.random_state = random_state
        self.n_jobs = n_jobs

        # Store prepared molecules (conformers will be generated on demand)
        self.actives = actives
        self.decoys = decoys

        # Ground truth labels
        self.y_true = np.array(
            [1] * len(actives) + [0] * len(decoys)
        )

        # Cache for prepared molecules with conformers
        self._prepared_mols_cache = {}

        # Prepare reference molecules for reference-based scoring
        self._prepared_refs = self._prepare_references(reference_mols)

        # Pharm2D scorer: lazy-initialized on first use to keep self picklable
        # (MolChemicalFeatureFactory cannot be pickled by loky workers)
        self._pharm2d_scorer = None

    def _prepare_references(self, mols: List[Chem.Mol]) -> List[Chem.Mol]:
        """Prepare reference molecules with 3D conformers for alignment.

        Args:
            mols: Reference molecules (may or may not have conformers).

        Returns:
            List of reference molecules with at least one conformer.
        """
        prepared = []
        for mol in mols:
            if mol is None:
                continue
            mol_h = Chem.AddHs(mol)
            if mol_h.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol_h, randomSeed=self.random_state)
            if mol_h.GetNumConformers() > 0:
                prepared.append(mol_h)
        return prepared

    def _get_pharm2d_scorer(self):
        """Lazy-initialize Pharm2DScorer on first use.

        Kept out of __init__ because MolChemicalFeatureFactory is not
        picklable, which breaks joblib loky parallelization of bound
        methods on this class.
        """
        if self._pharm2d_scorer is None:
            from .pharm2d_scoring import Pharm2DScorer
            self._pharm2d_scorer = Pharm2DScorer(self.reference_mols)
        return self._pharm2d_scorer

    def _prepare_molecules(
        self,
        mols: List[Chem.Mol],
        n_conformers: int,
        minimize: bool = False
    ) -> List[Chem.Mol]:
        """Generate conformers for molecules if needed.

        Uses caching to avoid regenerating conformers for the same configuration.
        """
        if mols:
            try:
                mol_id = Chem.MolToSmiles(mols[0], canonical=True)
            except Exception:
                mol_id = id(mols[0])
        else:
            mol_id = 0
        cache_key = (mol_id, len(mols), n_conformers, minimize)
        if cache_key in self._prepared_mols_cache:
            return self._prepared_mols_cache[cache_key]

        prepared = []
        for mol in mols:
            if mol is None:
                continue

            # Reuse existing conformers if not minimizing and enough exist
            if not minimize and mol.GetNumConformers() >= n_conformers:
                prepared.append(mol)
                continue

            # Generate conformers
            # numThreads=1 to avoid RDKit C++ thread pool heap corruption
            # when loky workers later pickle these molecules for scoring.
            mol_h = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = self.random_state
            params.numThreads = 1
            params.pruneRmsThresh = 0.5
            AllChem.EmbedMultipleConfs(
                mol_h, numConfs=n_conformers, params=params
            )
            if minimize and mol_h.GetNumConformers() > 0:
                AllChem.MMFFOptimizeMoleculeConfs(mol_h, numThreads=1, maxIters=200)

            if mol_h.GetNumConformers() > 0:
                prepared.append(mol_h)

        self._prepared_mols_cache[cache_key] = prepared
        return prepared

    def _score_molecule(
        self,
        mol: Chem.Mol,
        pharm_mol: Chem.Mol,
        shape_weight: float,
        opt_param: float
    ) -> float:
        """Score a single molecule against pharmacophore.

        Args:
            mol: Query molecule with conformers.
            pharm_mol: Pharmacophore mol from PharmacophoreToMol.convert().
            shape_weight: Weight for shape vs color (0-1).
            opt_param: AlignMol optimization parameter (0-1).

        Returns:
            Best score across all conformers.
        """
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
                    opt_param=opt_param
                )
                score = shape_weight * shape + color_weight * color
                best_score = max(best_score, score)
            except Exception as e:
                logger.debug("AlignMol failed for conformer %d: %s", conf_id, e)
                continue

        return best_score

    def _score_molecule_reference(
        self,
        mol: Chem.Mol,
        opt_param: float,
        aggregation: str = 'max'
    ) -> float:
        """Score a molecule against the reference ensemble.

        Aligns query molecule conformers to each prepared reference molecule
        and returns the aggregated combo Tanimoto score. This bypasses
        PharmacophoreToMol and aligns to real (connected) molecules.

        Args:
            mol: Query molecule with conformers.
            opt_param: AlignMol optimization parameter (0.0-1.0).
            aggregation: ``'max'`` for best reference, ``'mean'`` for average.

        Returns:
            Aggregated combo Tanimoto score (0.0-2.0).
        """
        if mol is None:
            return 0.0

        ref_scores = []
        for ref in self._prepared_refs:
            best_combo = 0.0
            for conf_id in range(mol.GetNumConformers()):
                try:
                    shape, color = AlignMol(
                        ref=ref,
                        probe=mol,
                        probeConfId=conf_id,
                        useColors=True,
                        opt_param=opt_param
                    )
                    best_combo = max(best_combo, shape + color)
                except Exception:
                    continue
            ref_scores.append(best_combo)

        if not ref_scores:
            return 0.0

        if aggregation == 'mean':
            return float(np.mean(ref_scores))
        return float(max(ref_scores))

    def _score_all_molecules(
        self,
        config: EvaluationConfig,
        prepared_actives: List[Chem.Mol],
        prepared_decoys: List[Chem.Mol],
    ) -> np.ndarray:
        """Score all molecules using the configured scoring mode.

        Dispatches to the appropriate scoring method based on
        ``config.scoring_mode``.

        Args:
            config: Evaluation configuration.
            prepared_actives: Active molecules with conformers.
            prepared_decoys: Decoy molecules with conformers.

        Returns:
            Array of scores for all molecules (actives then decoys).
        """
        if config.scoring_mode == 'consensus_mol':
            warnings.warn(
                "scoring_mode='consensus_mol' is deprecated and produces "
                "anti-discriminative results (AUC < 0.5). "
                "Use scoring_mode='reference' instead.",
                DeprecationWarning,
                stacklevel=3
            )
            # Legacy PharmacophoreToMol path — caller must pass pharm_mol
            # via the evaluate() method directly; this branch is handled
            # inline in evaluate() for backward compatibility.
            return None

        elif config.scoring_mode == 'reference':
            all_mols = prepared_actives + prepared_decoys
            if self.n_jobs == 1:
                scores = [
                    self._score_molecule_reference(m, config.opt_param, config.aggregation)
                    for m in all_mols
                ]
            else:
                scores = Parallel(n_jobs=self.n_jobs, backend='loky')(
                    delayed(self._score_molecule_reference)(
                        m, config.opt_param, config.aggregation
                    )
                    for m in all_mols
                )
            return np.array(scores)

        elif config.scoring_mode == 'pharm2d':
            all_mols = self.actives + self.decoys  # No conformers needed
            scores = [self._get_pharm2d_scorer().score(m) for m in all_mols]
            return np.array(scores)

        elif config.scoring_mode == 'hybrid':
            # 2D component (no conformers needed)
            scores_2d = np.array([
                self._get_pharm2d_scorer().score(m)
                for m in (self.actives + self.decoys)
            ])
            # 3D component (reference-based, needs conformers)
            all_mols = prepared_actives + prepared_decoys
            if self.n_jobs == 1:
                raw_3d = [
                    self._score_molecule_reference(m, config.opt_param, config.aggregation)
                    for m in all_mols
                ]
            else:
                raw_3d = Parallel(n_jobs=self.n_jobs, backend='loky')(
                    delayed(self._score_molecule_reference)(
                        m, config.opt_param, config.aggregation
                    )
                    for m in all_mols
                )
            # Normalize [0,2] → [0,1]
            scores_3d = np.array(raw_3d) / 2.0
            return config.alpha * scores_2d + (1.0 - config.alpha) * scores_3d

        else:
            raise ValueError(f"Unknown scoring_mode: {config.scoring_mode}")

    def evaluate(self, config: EvaluationConfig) -> EvaluationResult:
        """Evaluate a parameter configuration.

        Pipeline:
        1. Generate consensus pharmacophore from references
        2. Prepare query molecules with conformers
        3. Score all molecules using ``config.scoring_mode``
        4. Calculate all screening metrics

        For ``scoring_mode='reference'`` (default), queries are aligned
        directly to reference molecules—bypassing PharmacophoreToMol.

        Args:
            config: Configuration to evaluate.

        Returns:
            EvaluationResult with all metrics populated.
        """
        start_time = time.time()

        # Step 1: Generate consensus features
        # For reference/pharm2d/hybrid modes, consensus is only used for
        # feature count reporting. Skip it if not needed for scoring
        # (consensus_mol mode still requires it for PharmacophoreToMol).
        if config.scoring_mode == 'consensus_mol':
            consensus = PharmacophoreConsensus(
                tolerance=config.tolerance,
                occurrence_threshold=config.occurrence,
                linkage=config.linkage
            )
            features = consensus.generate_consensus(self.reference_mols)

            if len(features) < 2:
                return EvaluationResult(
                    config=config,
                    roc_auc=0.5,
                    bedroc=0.0,
                    ef_1=0.0,
                    ef_5=0.0,
                    ef_10=0.0,
                    n_features=len(features),
                    eval_time_sec=time.time() - start_time
                )
        else:
            # For reference-based modes, use cached consensus or skip
            # Generate lazily only for n_features reporting
            if not hasattr(self, '_cached_consensus_features'):
                consensus = PharmacophoreConsensus(
                    tolerance=config.tolerance,
                    occurrence_threshold=config.occurrence,
                    linkage=config.linkage
                )
                self._cached_consensus_features = consensus.generate_consensus(
                    self.reference_mols
                )
            features = self._cached_consensus_features

        # Step 2: Prepare molecules with conformers
        # Conformers are cached by (n_conformers, minimize) so repeated
        # calls with the same config skip regeneration.
        prepared_actives = self._prepare_molecules(
            self.actives, config.n_conformers, minimize=config.minimize)
        prepared_decoys = self._prepare_molecules(
            self.decoys, config.n_conformers, minimize=config.minimize)

        # Step 3: Score using the configured mode
        if config.scoring_mode == 'consensus_mol':
            # Legacy path (deprecated)
            warnings.warn(
                "scoring_mode='consensus_mol' is deprecated and produces "
                "anti-discriminative results (AUC < 0.5). "
                "Use scoring_mode='reference' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            pharm_mol = PharmacophoreToMol.convert(
                features, name='Consensus', enable_color_features=True
            )
            if pharm_mol is None:
                return EvaluationResult(
                    config=config,
                    roc_auc=0.5, bedroc=0.0,
                    ef_1=0.0, ef_5=0.0, ef_10=0.0,
                    n_features=len(features),
                    eval_time_sec=time.time() - start_time
                )
            if self.n_jobs == 1:
                active_scores = [
                    self._score_molecule(m, pharm_mol, config.shape_weight, config.opt_param)
                    for m in prepared_actives
                ]
                decoy_scores = [
                    self._score_molecule(m, pharm_mol, config.shape_weight, config.opt_param)
                    for m in prepared_decoys
                ]
            else:
                active_scores = Parallel(n_jobs=self.n_jobs, backend='loky')(
                    delayed(self._score_molecule)(m, pharm_mol, config.shape_weight, config.opt_param)
                    for m in prepared_actives
                )
                decoy_scores = Parallel(n_jobs=self.n_jobs, backend='loky')(
                    delayed(self._score_molecule)(m, pharm_mol, config.shape_weight, config.opt_param)
                    for m in prepared_decoys
                )
            y_scores = np.array(active_scores + decoy_scores)
        else:
            y_scores = self._score_all_molecules(config, prepared_actives, prepared_decoys)

        # Step 4: Calculate all metrics
        try:
            metrics = calculate_all_metrics(self.y_true.tolist(), y_scores.tolist())
        except Exception as e:
            logger.warning("calculate_all_metrics failed: %s", e)
            metrics = {
                'roc_auc': 0.5,
                'bedroc': 0.0,
                'ef_1': 0.0,
                'ef_5': 0.0,
                'ef_10': 0.0
            }

        eval_time = time.time() - start_time

        return EvaluationResult(
            config=config,
            roc_auc=metrics['roc_auc'],
            bedroc=metrics['bedroc'],
            ef_1=metrics['ef_1'],
            ef_5=metrics['ef_5'],
            ef_10=metrics['ef_10'],
            n_features=len(features),
            eval_time_sec=eval_time
        )

    def evaluate_feature_subset(
        self,
        features: List,
        shape_weight: float = 0.5,
        opt_param: float = 0.5,
        n_conformers: int = 25,
        minimize: bool = False,
        scoring_mode: str = 'reference'
    ) -> EvaluationResult:
        """Evaluate a specific feature list directly (for HypoGen).

        This bypasses consensus generation and evaluates a pre-defined
        feature set. Used by HypoGen's constructive-subtractive phase.

        Args:
            features: List of pharmacophore features.
            shape_weight: Weight for shape vs color (0-1).
            opt_param: AlignMol optimization parameter (0-1).
            n_conformers: Number of conformers per query molecule.
            minimize: Whether to MMFF-minimize conformers.
            scoring_mode: Scoring approach (default ``'reference'``).

        Returns:
            EvaluationResult with all metrics.
        """
        start_time = time.time()

        # Create config for tracking
        config = EvaluationConfig(
            tolerance=2.0,  # Not applicable for direct feature evaluation
            occurrence=0.5,  # Not applicable
            shape_weight=shape_weight,
            opt_param=opt_param,
            n_conformers=n_conformers,
            scoring_mode=scoring_mode
        )

        # Degenerate case
        if len(features) < 2:
            return EvaluationResult(
                config=config,
                roc_auc=0.5,
                bedroc=0.0,
                ef_1=0.0,
                ef_5=0.0,
                ef_10=0.0,
                n_features=len(features),
                eval_time_sec=time.time() - start_time
            )

        # Prepare molecules
        prepared_actives = self._prepare_molecules(
            self.actives, n_conformers, minimize=minimize)
        prepared_decoys = self._prepare_molecules(
            self.decoys, n_conformers, minimize=minimize)

        # Score molecules based on scoring_mode
        if scoring_mode == 'consensus_mol':
            warnings.warn(
                "scoring_mode='consensus_mol' is deprecated and produces "
                "anti-discriminative results. Use 'reference' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            pharm_mol = PharmacophoreToMol.convert(
                features, name='Direct', enable_color_features=True
            )
            if pharm_mol is None:
                return EvaluationResult(
                    config=config,
                    roc_auc=0.5, bedroc=0.0,
                    ef_1=0.0, ef_5=0.0, ef_10=0.0,
                    n_features=len(features),
                    eval_time_sec=time.time() - start_time
                )
            if self.n_jobs == 1:
                active_scores = [
                    self._score_molecule(m, pharm_mol, shape_weight, opt_param)
                    for m in prepared_actives
                ]
                decoy_scores = [
                    self._score_molecule(m, pharm_mol, shape_weight, opt_param)
                    for m in prepared_decoys
                ]
            else:
                active_scores = Parallel(n_jobs=self.n_jobs, backend='loky')(
                    delayed(self._score_molecule)(m, pharm_mol, shape_weight, opt_param)
                    for m in prepared_actives
                )
                decoy_scores = Parallel(n_jobs=self.n_jobs, backend='loky')(
                    delayed(self._score_molecule)(m, pharm_mol, shape_weight, opt_param)
                    for m in prepared_decoys
                )
            y_scores = np.array(active_scores + decoy_scores)
        else:
            y_scores = self._score_all_molecules(config, prepared_actives, prepared_decoys)

        # Calculate metrics
        try:
            metrics = calculate_all_metrics(self.y_true.tolist(), y_scores.tolist())
        except Exception as e:
            logger.warning("calculate_all_metrics failed in evaluate_feature_subset: %s", e)
            metrics = {
                'roc_auc': 0.5,
                'bedroc': 0.0,
                'ef_1': 0.0,
                'ef_5': 0.0,
                'ef_10': 0.0
            }

        eval_time = time.time() - start_time

        return EvaluationResult(
            config=config,
            roc_auc=metrics['roc_auc'],
            bedroc=metrics['bedroc'],
            ef_1=metrics['ef_1'],
            ef_5=metrics['ef_5'],
            ef_10=metrics['ef_10'],
            n_features=len(features),
            eval_time_sec=eval_time
        )

    def evaluate_ensemble(
        self,
        config: EvaluationConfig,
        quick_mode: bool = True,
        min_auc_threshold: float = 0.60,
        rrf_constant: int = RRF_CONSTANT
    ) -> EvaluationResult:
        """Evaluate with cascading multi-level scoring ensemble.

        Uses reference-based 3D scoring as the primary signal. If quick_mode
        is enabled, checks whether the 3D AUC exceeds the threshold before
        running the more expensive OT-based scoring for rank fusion.

        This implements the multi-fidelity insight from McDonald et al. 2025:
        use cheap evaluation for exploration, expensive for exploitation.

        Args:
            config: Configuration to evaluate.
            quick_mode: If True, skip OT scoring when 3D score
                is below threshold.
            min_auc_threshold: Minimum 3D AUC to proceed to OT
                scoring (only used when quick_mode=True).
            rrf_constant: Reciprocal rank fusion constant (default: 60).

        Returns:
            EvaluationResult with ensemble metrics.
        """
        start_time = time.time()

        # Step 1: Generate consensus (for feature count + OT scoring)
        consensus = PharmacophoreConsensus(
            tolerance=config.tolerance,
            occurrence_threshold=config.occurrence,
            linkage=config.linkage
        )
        features = consensus.generate_consensus(self.reference_mols)

        if len(features) < 2:
            return EvaluationResult(
                config=config,
                roc_auc=0.5, bedroc=0.0,
                ef_1=0.0, ef_5=0.0, ef_10=0.0,
                n_features=len(features),
                eval_time_sec=time.time() - start_time
            )

        # Step 2: Prepare molecules with conformers
        prepared_actives = self._prepare_molecules(
            self.actives, config.n_conformers, minimize=config.minimize)
        prepared_decoys = self._prepare_molecules(
            self.decoys, config.n_conformers, minimize=config.minimize)

        # Step 3: Reference-based 3D scoring (replaces PharmacophoreToMol path)
        all_prepared = prepared_actives + prepared_decoys
        if self.n_jobs == 1:
            scores_3d = [
                self._score_molecule_reference(m, config.opt_param, config.aggregation)
                for m in all_prepared
            ]
        else:
            scores_3d = Parallel(n_jobs=self.n_jobs, backend='loky')(
                delayed(self._score_molecule_reference)(
                    m, config.opt_param, config.aggregation
                )
                for m in all_prepared
            )
        y_scores_3d = np.array(scores_3d)

        # Step 4: Quick mode check — is this config worth scoring further?
        try:
            from sklearn.metrics import roc_auc_score
            quick_auc = roc_auc_score(self.y_true, y_scores_3d)
        except Exception as e:
            logger.debug("Quick AUC calculation failed: %s", e)
            quick_auc = 0.5

        if quick_mode and quick_auc < min_auc_threshold:
            metrics = calculate_all_metrics(self.y_true.tolist(), y_scores_3d.tolist())
            return EvaluationResult(
                config=config,
                roc_auc=metrics['roc_auc'],
                bedroc=metrics['bedroc'],
                ef_1=metrics['ef_1'],
                ef_5=metrics['ef_5'],
                ef_10=metrics['ef_10'],
                n_features=len(features),
                eval_time_sec=time.time() - start_time
            )

        # Step 5: Also try OT-based scoring for rank fusion
        try:
            from .ot_scoring import wasserstein_similarity

            # Score each molecule by Wasserstein similarity to consensus
            from pharmacophore.pharmacophore import Pharmacophore
            p = Pharmacophore()

            ot_scores = []
            for m in all_prepared:
                try:
                    mol_feats = p.calc_pharm(mol=m)
                    sim = wasserstein_similarity(features, mol_feats)
                except Exception as e:
                    logger.debug("OT scoring failed: %s", e)
                    sim = 0.0
                ot_scores.append(sim)

            y_scores_ot = np.array(ot_scores)

            # Reciprocal rank fusion: combine 3D and OT rankings
            ranks_3d = _rankdata_descending(y_scores_3d)
            ranks_ot = _rankdata_descending(y_scores_ot)
            k = rrf_constant
            rrf_scores = 1.0 / (k + ranks_3d) + 1.0 / (k + ranks_ot)
            y_scores_final = rrf_scores

        except (ImportError, Exception):
            y_scores_final = y_scores_3d

        # Final metrics
        try:
            metrics = calculate_all_metrics(self.y_true.tolist(), y_scores_final.tolist())
        except Exception as e:
            logger.warning("Final metrics calculation failed: %s", e)
            metrics = {
                'roc_auc': 0.5, 'bedroc': 0.0,
                'ef_1': 0.0, 'ef_5': 0.0, 'ef_10': 0.0
            }

        return EvaluationResult(
            config=config,
            roc_auc=metrics['roc_auc'],
            bedroc=metrics['bedroc'],
            ef_1=metrics['ef_1'],
            ef_5=metrics['ef_5'],
            ef_10=metrics['ef_10'],
            n_features=len(features),
            eval_time_sec=time.time() - start_time
        )

    def __repr__(self) -> str:
        return (
            f"UnifiedEvaluator("
            f"n_refs={len(self.reference_mols)}, "
            f"n_actives={len(self.actives)}, "
            f"n_decoys={len(self.decoys)})"
        )


def _rankdata_descending(scores: np.ndarray) -> np.ndarray:
    """Rank scores in descending order (highest score = rank 1)."""
    temp = np.argsort(-scores)
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(1, len(scores) + 1)
    return ranks

# TECHNICAL DEBT: 1237 lines, 18 methods, 8 responsibilities.
# Consider decomposing into: DataLoader, ConformerGenerator,
# ReferenceAligner, Scorer, BayesianOptimizer, ModelExporter.

"""Combo Pharmacophore Optimizer using rdShapeAlign combo Tanimoto.

Builds an optimal consensus pharmacophore model from reference SMILES
(no binding poses) and optimizes parameters via Bayesian optimization
to maximize ROC AUC for discriminating actives vs decoys.

Pipeline:
    1. Load 5 reference SMILES (discard any 3D poses)
    2. Generate conformers de novo
    3. Align references using rdShapeAlign (template = hyperparameter)
    4. Build consensus pharmacophore
    5. Convert to RDKit Mol with color features
    6. Score actives/decoys using combo Tanimoto (shape + color)
    7. Optimize via Bayesian GP to maximize ROC AUC

Example:
    >>> from pharmacophore.combo_optimizer import ComboPharmacophoreOptimizer
    >>> opt = ComboPharmacophoreOptimizer()
    >>> opt.load_from_files(
    ...     'tutorials/data/CCR2_reference_ligands.sdf',
    ...     'tutorials/data/actives_ccr2_N75.csv',
    ...     'tutorials/data/decoys_ccr2_N500.csv'
    ... )
    >>> result = opt.optimize(n_calls=50)
    >>> print(f"Best AUC: {result['best_auc']:.4f}")
    >>> opt.export_model('output/combo_optimal')
"""

import json
import time
import logging
import warnings
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeAlign
from sklearn.metrics import roc_auc_score
from tqdm import tqdm

from .caching import LRUCache
from .consensus import PharmacophoreConsensus
from .mol_converter import PharmacophoreToMol
from .screening_metrics import calculate_all_metrics
from .constants import FEATURE_COLORS

logger = logging.getLogger(__name__)

# Optional: scikit-optimize
try:
    from skopt import gp_minimize
    from skopt.space import Categorical, Integer, Real
    from skopt.callbacks import DeltaYStopper

    HAS_SKOPT = True
except ImportError:
    HAS_SKOPT = False


class ComboPharmacophoreOptimizer:
    """Bayesian optimizer for consensus pharmacophore using combo Tanimoto.

    Uses rdShapeAlign (shape + color Tanimoto, range 0-2) to score molecules
    against a consensus pharmacophore model. Optimizes consensus parameters
    and scoring settings via Gaussian Process Bayesian optimization.

    The 5 reference molecules are loaded from SMILES only (3D poses discarded)
    to avoid information leakage from known binding poses.

    Args:
        random_state: Random seed for reproducibility.
        verbose: Print progress during optimization.
        scoring_mode: Scoring strategy for evaluation. Options:
            ``'reference'`` (default) -- score query molecules directly
            against the aligned reference molecules, bypassing the
            PharmacophoreToMol conversion step.
            ``'consensus_mol'`` (legacy) -- convert the consensus
            pharmacophore to an RDKit Mol and score against that.
            Emits a DeprecationWarning.
    """

    def __init__(
        self,
        random_state: int = 42,
        verbose: bool = True,
        scoring_mode: str = "reference",
    ):
        if not HAS_SKOPT:
            raise ImportError(
                "scikit-optimize required. Install with: pip install scikit-optimize"
            )

        self.random_state = random_state
        self.verbose = verbose
        self.scoring_mode = scoring_mode
        self.n_ref_conformers = 10

        # Data
        self._ref_smiles: Optional[List[str]] = None
        self._active_smiles: Optional[List[str]] = None
        self._decoy_smiles: Optional[List[str]] = None

        # Prepared reference molecules for reference-based scoring
        self._prepared_refs: List[Chem.Mol] = []

        # Caches (4 levels, bounded via LRUCache)
        self._ref_conformers = LRUCache(max_size=50)      # L1: SMILES -> Mol
        self._query_conformers = LRUCache(max_size=600)    # L2: (SMILES,n) -> Mol
        self._aligned_refs = LRUCache(max_size=20)         # L3: template_idx -> mols
        self._consensus_cache = LRUCache(max_size=200)     # L4

        # Results
        self.history: List[Dict[str, Any]] = []
        self.best_result: Optional[Dict[str, Any]] = None

    # -------------------------------------------------------------------------
    # Data Loading
    # -------------------------------------------------------------------------

    def load_from_files(
        self,
        reference_file: str,
        actives_file: str,
        decoys_file: str,
        smiles_column: Optional[str] = None,
    ) -> "ComboPharmacophoreOptimizer":
        """Load reference SDF, actives CSV, and decoys CSV.

        Extracts SMILES from the SDF file (3D coordinates are discarded).
        Auto-detects SMILES column in CSV files.

        Args:
            reference_file: Path to SDF with reference molecules.
            actives_file: Path to CSV with active SMILES.
            decoys_file: Path to CSV with decoy SMILES.
            smiles_column: SMILES column name (auto-detect if None).

        Returns:
            self for method chaining.
        """
        # Validate file paths
        reference_path = Path(reference_file).resolve()
        actives_path = Path(actives_file).resolve()
        decoys_path = Path(decoys_file).resolve()
        for p, label in [(reference_path, "Reference"), (actives_path, "Actives"), (decoys_path, "Decoys")]:
            if not p.exists():
                raise FileNotFoundError(f"{label} file not found: {p}")

        # Extract SMILES from SDF (discard 3D)
        ref_smiles = []
        supplier = Chem.SDMolSupplier(str(reference_path), removeHs=True)
        for mol in supplier:
            if mol is not None:
                smi = Chem.MolToSmiles(mol)
                ref_smiles.append(smi)

        if not ref_smiles:
            raise ValueError(f"No valid molecules in {reference_file}")

        # Load actives CSV
        active_smiles = self._load_smiles_from_csv(actives_file, smiles_column)
        # Load decoys CSV
        decoy_smiles = self._load_smiles_from_csv(decoys_file, smiles_column)

        return self.load_from_smiles(ref_smiles, active_smiles, decoy_smiles)

    def load_from_smiles(
        self,
        reference_smiles: List[str],
        active_smiles: List[str],
        decoy_smiles: List[str],
    ) -> "ComboPharmacophoreOptimizer":
        """Load molecules from SMILES lists.

        Args:
            reference_smiles: Reference SMILES (typically 5).
            active_smiles: Active compound SMILES.
            decoy_smiles: Decoy compound SMILES.

        Returns:
            self for method chaining.
        """
        if len(reference_smiles) < 2:
            raise ValueError(
                f"Need at least 2 reference SMILES, got {len(reference_smiles)}"
            )
        if not active_smiles:
            raise ValueError("Active SMILES list is empty")
        if not decoy_smiles:
            raise ValueError("Decoy SMILES list is empty")

        self._ref_smiles = list(reference_smiles)
        self._active_smiles = list(active_smiles)
        self._decoy_smiles = list(decoy_smiles)

        # Clear caches (new data)
        self.clear_cache()

        # Pre-generate reference conformers
        self._generate_ref_conformers()

        if self.verbose:
            print(f"Loaded: {len(self._ref_smiles)} refs, "
                  f"{len(self._active_smiles)} actives, "
                  f"{len(self._decoy_smiles)} decoys")

        return self

    @staticmethod
    def _load_smiles_from_csv(
        filepath: str, smiles_column: Optional[str] = None
    ) -> List[str]:
        """Load SMILES from CSV, auto-detecting column name."""
        csv_path = Path(filepath).resolve()
        if not csv_path.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_path}")
        df = pd.read_csv(csv_path)
        if smiles_column and smiles_column in df.columns:
            col = smiles_column
        else:
            candidates = ["SMILES", "smiles", "Smiles", "smi", "canonical_smiles"]
            col = None
            for c in candidates:
                if c in df.columns:
                    col = c
                    break
            if col is None:
                raise ValueError(
                    f"No SMILES column found in {filepath}. "
                    f"Columns: {list(df.columns)}"
                )
        smiles = df[col].dropna().astype(str).tolist()
        # Validate SMILES
        valid = [s for s in smiles if Chem.MolFromSmiles(s) is not None]
        if len(valid) < len(smiles):
            logger.warning(
                f"Skipped {len(smiles) - len(valid)} invalid SMILES in {filepath}"
            )
        return valid

    # -------------------------------------------------------------------------
    # Conformer Generation
    # -------------------------------------------------------------------------

    def _generate_conformers(
        self,
        smiles: str,
        n_conformers: int,
        cache: LRUCache,
    ) -> Optional[Chem.Mol]:
        """Generate 3D conformers for a SMILES string (cached)."""
        cache_key = f"{smiles}_{n_conformers}"
        cached = cache.get(cache_key)
        if cached is not None:
            return cached

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol_h = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = self.random_state
        AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)

        if mol_h.GetNumConformers() == 0:
            # Fallback: single conformer
            AllChem.EmbedMolecule(mol_h, randomSeed=self.random_state)

        if mol_h.GetNumConformers() == 0:
            logger.warning(f"Failed to generate conformers for: {smiles}")
            return None

        cache.set(cache_key, mol_h)
        return mol_h

    def _generate_ref_conformers(self, n_conformers: Optional[int] = None) -> None:
        """Pre-generate conformers for all reference molecules."""
        if n_conformers is None:
            n_conformers = self.n_ref_conformers
        for smi in self._ref_smiles:
            self._generate_conformers(smi, n_conformers, self._ref_conformers)

    # -------------------------------------------------------------------------
    # Reference Alignment
    # -------------------------------------------------------------------------

    def _align_references_to_template(
        self, template_idx: int
    ) -> List[Chem.Mol]:
        """Align reference molecules to a chosen template using rdShapeAlign.

        For each non-template reference, tries all conformer pairs against
        the template and keeps the best combo Tanimoto alignment.

        Args:
            template_idx: Index of the template reference molecule.

        Returns:
            List of aligned RDKit Mol objects (all in template's frame).
        """
        cached_aligned = self._aligned_refs.get(template_idx)
        if cached_aligned is not None:
            return cached_aligned

        n_refs = len(self._ref_smiles)
        template_smi = self._ref_smiles[template_idx]
        template_mol = self._generate_conformers(
            template_smi, self.n_ref_conformers, self._ref_conformers
        )

        if template_mol is None:
            raise ValueError(
                f"Cannot generate conformers for template ref {template_idx}"
            )

        aligned = [None] * n_refs

        # Template stays as-is (use first conformer)
        template_copy = Chem.Mol(template_mol)
        aligned[template_idx] = template_copy

        for i in range(n_refs):
            if i == template_idx:
                continue

            probe_smi = self._ref_smiles[i]
            probe_mol = self._generate_conformers(
                probe_smi, self.n_ref_conformers, self._ref_conformers
            )
            if probe_mol is None:
                logger.warning(f"Skipping ref {i}: no conformers")
                continue

            best_combo = -1.0
            best_probe = None

            # Try all conformer pairs
            n_template_confs = template_mol.GetNumConformers()
            n_probe_confs = probe_mol.GetNumConformers()

            early_stop = False
            for ref_cid in range(n_template_confs):
                if early_stop:
                    break
                for probe_cid in range(n_probe_confs):
                    # Deep copy probe (AlignMol mutates in-place)
                    probe_copy = Chem.Mol(probe_mol)
                    try:
                        shape, color = rdShapeAlign.AlignMol(
                            ref=template_mol,
                            probe=probe_copy,
                            refConfId=ref_cid,
                            probeConfId=probe_cid,
                            useColors=True,
                            opt_param=0.5,
                            max_preiters=50,
                            max_postiters=100,
                        )
                        combo = shape + color
                        if combo > best_combo:
                            best_combo = combo
                            best_probe = probe_copy
                        if best_combo > 1.5:  # Near-optimal; combo Tanimoto max is 2.0
                            early_stop = True
                            break
                    except Exception as e:
                        logger.debug(f"Alignment failed ref {template_idx} "
                                     f"conf {ref_cid} -> probe {i} "
                                     f"conf {probe_cid}: {e}")
                        continue

            if best_probe is not None:
                aligned[i] = best_probe
            else:
                logger.warning(
                    f"All alignments failed for ref {i} to template {template_idx}"
                )

        # Filter out None entries
        valid_aligned = [m for m in aligned if m is not None]
        if len(valid_aligned) < 2:
            raise ValueError(
                f"Only {len(valid_aligned)} refs aligned successfully "
                f"(need at least 2)"
            )

        self._aligned_refs.set(template_idx, valid_aligned)
        return valid_aligned

    # -------------------------------------------------------------------------
    # Consensus Building
    # -------------------------------------------------------------------------

    def _build_consensus(
        self,
        template_idx: int,
        tolerance: float,
        occurrence_threshold: float,
    ) -> Tuple[List, Optional[Chem.Mol]]:
        """Build consensus pharmacophore from aligned references.

        Returns:
            Tuple of (consensus_features, pharm_mol).
            pharm_mol is None if fewer than 2 features found.
        """
        cache_key = f"{template_idx}_{tolerance:.4f}_{occurrence_threshold:.4f}"
        cached_result = self._consensus_cache.get(cache_key)
        if cached_result is not None:
            return cached_result

        aligned_mols = self._align_references_to_template(template_idx)

        consensus_gen = PharmacophoreConsensus(
            tolerance=tolerance,
            occurrence_threshold=occurrence_threshold,
        )

        try:
            features = consensus_gen.generate_consensus(aligned_mols)
        except Exception as e:
            logger.warning(f"Consensus generation failed: {e}")
            features = []

        if len(features) < 2:
            result = (features, None)
        else:
            try:
                pharm_mol = PharmacophoreToMol.convert(
                    features, enable_color_features=True
                )
                result = (features, pharm_mol)
            except Exception as e:
                logger.warning(f"Mol conversion failed: {e}")
                result = (features, None)

        self._consensus_cache.set(cache_key, result)
        return result

    # -------------------------------------------------------------------------
    # Scoring
    # -------------------------------------------------------------------------

    def _score_molecule(
        self,
        smiles: str,
        pharm_mol: Chem.Mol,
        opt_param: float,
        n_conformers: int,
    ) -> float:
        """Score a single molecule against the consensus pharmacophore Mol.

        Returns best combo Tanimoto score (0-2), or 0.0 on failure.
        """
        query_mol = self._generate_conformers(
            smiles, n_conformers, self._query_conformers
        )
        if query_mol is None or query_mol.GetNumConformers() == 0:
            return 0.0

        best_combo = 0.0

        for conf_id in range(query_mol.GetNumConformers()):
            try:
                shape, color = rdShapeAlign.AlignMol(
                    ref=pharm_mol,
                    probe=query_mol,
                    probeConfId=conf_id,
                    useColors=True,
                    opt_param=opt_param,
                    max_preiters=50,
                    max_postiters=100,
                )
                combo = shape + color
                best_combo = max(best_combo, combo)
            except Exception:
                continue

        return best_combo

    def _score_molecule_reference(
        self,
        mol: Chem.Mol,
        shape_weight: float = 0.5,
        color_weight: float = 0.5
    ) -> float:
        """Score a single molecule against prepared reference molecules."""
        if mol is None or not self._prepared_refs:
            return 0.0

        best_score = 0.0
        for ref in self._prepared_refs:
            for conf_id in range(mol.GetNumConformers()):
                try:
                    shape, color = rdShapeAlign.AlignMol(
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

    def _score_all_reference(
        self,
        smiles_list: List[str],
        n_conformers: int,
        shape_weight: float = 0.5,
        color_weight: float = 0.5,
        desc: str = "Scoring",
    ) -> List[float]:
        """Score a list of molecules against prepared references.

        Returns list of best combo scores per molecule.
        """
        scores = []
        iterator = smiles_list
        if self.verbose:
            iterator = tqdm(smiles_list, desc=desc, leave=False)

        for smi in iterator:
            query_mol = self._generate_conformers(
                smi, n_conformers, self._query_conformers
            )
            score = self._score_molecule_reference(
                query_mol, shape_weight, color_weight
            )
            scores.append(score)

        return scores

    def _score_all(
        self,
        smiles_list: List[str],
        pharm_mol: Chem.Mol,
        opt_param: float,
        n_conformers: int,
        desc: str = "Scoring",
    ) -> List[float]:
        """Score a list of molecules. Returns list of combo scores."""
        scores = []
        iterator = smiles_list
        if self.verbose:
            iterator = tqdm(smiles_list, desc=desc, leave=False)

        for smi in iterator:
            score = self._score_molecule(smi, pharm_mol, opt_param, n_conformers)
            scores.append(score)

        return scores

    # -------------------------------------------------------------------------
    # Evaluation
    # -------------------------------------------------------------------------

    def evaluate(
        self,
        template_idx: int,
        tolerance: float,
        occurrence_threshold: float,
        opt_param: float,
        n_conformers: int,
    ) -> float:
        """Evaluate a single parameter configuration.

        Builds consensus, scores all actives and decoys, computes ROC AUC.
        Scoring method depends on ``self.scoring_mode``:

        * ``'reference'`` -- score query molecules directly against the
          aligned reference molecules stored in ``self._prepared_refs``.
        * ``'consensus_mol'`` (legacy) -- convert the consensus
          pharmacophore to an RDKit Mol via :class:`PharmacophoreToMol`
          and score against that.  Emits a :class:`DeprecationWarning`.

        Args:
            template_idx: Template reference index.
            tolerance: Consensus clustering tolerance (Angstroms).
            occurrence_threshold: Min feature occurrence fraction.
            opt_param: Shape/color scoring balance (0=color, 1=shape).
            n_conformers: Conformers per query molecule.

        Returns:
            ROC AUC score (0.5 = random).
        """
        if self._ref_smiles is None:
            raise ValueError("No data loaded. Call load_from_files() or "
                             "load_from_smiles() first.")

        # Always build consensus (needed for feature count / export)
        features, pharm_mol = self._build_consensus(
            template_idx, tolerance, occurrence_threshold
        )

        # -----------------------------------------------------------------
        # Reference-based scoring path (new default)
        # -----------------------------------------------------------------
        if self.scoring_mode == "reference":
            # Populate _prepared_refs from the aligned reference molecules
            aligned_mols = self._align_references_to_template(template_idx)
            self._prepared_refs = [
                m for m in aligned_mols if m is not None
            ]

            if not self._prepared_refs:
                result = {
                    "template_idx": template_idx,
                    "tolerance": tolerance,
                    "occurrence_threshold": occurrence_threshold,
                    "opt_param": opt_param,
                    "n_conformers": n_conformers,
                    "n_features": len(features),
                    "auc": 0.5,
                    "scoring_mode": self.scoring_mode,
                }
                self.history.append(result)
                return 0.5

            shape_weight = 1.0 - opt_param
            color_weight = opt_param

            active_scores = self._score_all_reference(
                self._active_smiles, n_conformers,
                shape_weight=shape_weight,
                color_weight=color_weight,
                desc="Actives (ref)",
            )
            decoy_scores = self._score_all_reference(
                self._decoy_smiles, n_conformers,
                shape_weight=shape_weight,
                color_weight=color_weight,
                desc="Decoys (ref)",
            )

        # -----------------------------------------------------------------
        # Legacy consensus-mol scoring path
        # -----------------------------------------------------------------
        else:
            warnings.warn(
                "scoring_mode='consensus_mol' is deprecated and will be "
                "removed in a future release. Use scoring_mode='reference' "
                "for improved accuracy.",
                DeprecationWarning,
                stacklevel=2,
            )

            if pharm_mol is None:
                result = {
                    "template_idx": template_idx,
                    "tolerance": tolerance,
                    "occurrence_threshold": occurrence_threshold,
                    "opt_param": opt_param,
                    "n_conformers": n_conformers,
                    "n_features": len(features),
                    "auc": 0.5,
                    "scoring_mode": self.scoring_mode,
                }
                self.history.append(result)
                return 0.5

            active_scores = self._score_all(
                self._active_smiles, pharm_mol, opt_param, n_conformers,
                desc="Actives"
            )
            decoy_scores = self._score_all(
                self._decoy_smiles, pharm_mol, opt_param, n_conformers,
                desc="Decoys"
            )

        y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
        y_scores = active_scores + decoy_scores

        try:
            auc = roc_auc_score(y_true, y_scores)
        except Exception:
            auc = 0.5

        result = {
            "template_idx": template_idx,
            "tolerance": tolerance,
            "occurrence_threshold": occurrence_threshold,
            "opt_param": opt_param,
            "n_conformers": n_conformers,
            "n_features": len(features),
            "auc": auc,
            "scoring_mode": self.scoring_mode,
        }
        self.history.append(result)

        return auc

    # -------------------------------------------------------------------------
    # Optimization
    # -------------------------------------------------------------------------

    def optimize(
        self,
        n_calls: int = 50,
        n_random_starts: int = 15,
        tolerance_range: Tuple[float, float] = (0.5, 4.0),
        occurrence_range: Tuple[float, float] = (0.1, 1.0),
        opt_param_range: Tuple[float, float] = (0.0, 1.0),
        n_conformers_range: Tuple[int, int] = (1, 10),
        early_stopping: bool = True,
        early_stopping_delta: float = 0.001,
        early_stopping_n: int = 10,
        acq_func: str = "gp_hedge",
        enable_retest: bool = True,
    ) -> Dict[str, Any]:
        """Run Bayesian optimization to find optimal parameters.

        Args:
            n_calls: Total number of evaluations.
            n_random_starts: Initial random explorations.
            tolerance_range: Search range for tolerance.
            occurrence_range: Search range for occurrence.
            opt_param_range: Search range for opt_param.
            n_conformers_range: Search range for n_conformers.
            early_stopping: Enable early stopping.
            early_stopping_delta: Min improvement threshold.
            early_stopping_n: No-improvement iterations before stopping.
            acq_func: Acquisition function for GP optimizer (default:
                "gp_hedge"). Use "EI" for ablation baseline.
            enable_retest: Run post-hoc GP retest policy (default: True).
                Set False for ablation baseline.

        Returns:
            Dict with best_params, best_auc, best_metrics,
            best_features, n_evaluations, history, elapsed_sec.
        """
        if self._ref_smiles is None:
            raise ValueError("No data loaded.")

        # Auto-adjust n_random_starts if it exceeds n_calls
        if n_random_starts >= n_calls:
            n_random_starts = max(1, n_calls - 1)

        n_refs = len(self._ref_smiles)
        space = [
            Categorical(list(range(n_refs)), name="template_idx"),
            Real(tolerance_range[0], tolerance_range[1], name="tolerance"),
            Real(occurrence_range[0], occurrence_range[1], name="occurrence_threshold"),
            Real(opt_param_range[0], opt_param_range[1], name="opt_param"),
            Integer(n_conformers_range[0], n_conformers_range[1], name="n_conformers"),
        ]

        eval_count = [0]
        best_auc_so_far = [0.0]

        def objective(params):
            template_idx, tol, occ, opt_p, n_conf = params
            eval_count[0] += 1

            auc = self.evaluate(
                int(template_idx), float(tol), float(occ),
                float(opt_p), int(n_conf)
            )

            if auc > best_auc_so_far[0]:
                best_auc_so_far[0] = auc

            if self.verbose:
                print(
                    f"  [{eval_count[0]:3d}/{n_calls}] "
                    f"tmpl={int(template_idx)} tol={tol:.2f} "
                    f"occ={occ:.2f} opt={opt_p:.2f} "
                    f"nc={int(n_conf)} -> AUC={auc:.4f} "
                    f"(best={best_auc_so_far[0]:.4f})"
                )

            return -auc  # minimize negative = maximize

        callbacks = []
        if early_stopping:
            callbacks.append(
                DeltaYStopper(delta=early_stopping_delta, n_best=early_stopping_n)
            )

        if self.verbose:
            print("=" * 70)
            print("Combo Pharmacophore Optimizer")
            print("=" * 70)
            print(f"References:  {n_refs}")
            print(f"Actives:     {len(self._active_smiles)}")
            print(f"Decoys:      {len(self._decoy_smiles)}")
            print(f"Evaluations: {n_calls} ({n_random_starts} random starts)")
            print(f"Scoring:     rdShapeAlign combo Tanimoto (shape + color)")
            print("=" * 70)

        start_time = time.time()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            skopt_result = gp_minimize(
                objective,
                space,
                n_calls=n_calls,
                n_random_starts=n_random_starts,
                acq_func=acq_func,
                random_state=self.random_state,
                callback=callbacks if callbacks else None,
            )

        elapsed = time.time() - start_time

        # Extract best parameters
        best_params = {
            "template_idx": int(skopt_result.x[0]),
            "tolerance": float(skopt_result.x[1]),
            "occurrence_threshold": float(skopt_result.x[2]),
            "opt_param": float(skopt_result.x[3]),
            "n_conformers": int(skopt_result.x[4]),
        }
        best_auc = -skopt_result.fun

        # Post-hoc retest: identify high-variance predictions and re-evaluate
        # with 2x conformers (Bellamy 2022 retest policy)
        if enable_retest:
            try:
                gp_model = skopt_result.models[-1] if skopt_result.models else None
                if gp_model is not None and len(skopt_result.x_iters) >= 5:
                    X_all = np.array(skopt_result.x_iters)
                    y_all = np.array(skopt_result.func_vals)
                    pred_mean, pred_std = gp_model.predict(X_all, return_std=True)

                    # Find points where GP prediction disagrees with measurement
                    residuals = np.abs(pred_mean - y_all)
                    high_var_mask = residuals > pred_std

                    if np.any(high_var_mask):
                        # Re-evaluate top candidates with higher conformers
                        retest_indices = np.where(high_var_mask)[0]
                        # Only retest top 3 (sorted by best predicted score)
                        retest_sorted = retest_indices[
                            np.argsort(pred_mean[retest_indices])[:3]
                        ]

                        for idx in retest_sorted:
                            x_point = X_all[idx]
                            retest_n_conf = min(int(x_point[4]) * 2, 50)
                            try:
                                retest_params = {
                                    "template_idx": int(x_point[0]),
                                    "tolerance": float(x_point[1]),
                                    "occurrence_threshold": float(x_point[2]),
                                    "opt_param": float(x_point[3]),
                                    "n_conformers": retest_n_conf,
                                }
                                retest_features, retest_mol = self._build_consensus(
                                    retest_params["template_idx"],
                                    retest_params["tolerance"],
                                    retest_params["occurrence_threshold"],
                                )
                                if retest_mol is not None:
                                    retest_active = self._score_all(
                                        self._active_smiles, retest_mol,
                                        retest_params["opt_param"],
                                        retest_params["n_conformers"],
                                    )
                                    retest_decoy = self._score_all(
                                        self._decoy_smiles, retest_mol,
                                        retest_params["opt_param"],
                                        retest_params["n_conformers"],
                                    )
                                    from sklearn.metrics import roc_auc_score
                                    y_t = [1] * len(retest_active) + [0] * len(retest_decoy)
                                    retest_auc = roc_auc_score(y_t, retest_active + retest_decoy)

                                    if retest_auc > best_auc:
                                        best_auc = retest_auc
                                        best_params = retest_params
                                        if self.verbose:
                                            print(f"  [Retest] Improved AUC: {retest_auc:.4f}")
                            except Exception:
                                continue
            except Exception:
                pass  # Retest is optional, don't fail the optimization

        # Re-evaluate best to get full metrics
        features, pharm_mol = self._build_consensus(
            best_params["template_idx"],
            best_params["tolerance"],
            best_params["occurrence_threshold"],
        )

        best_metrics = {}
        if pharm_mol is not None:
            active_scores = self._score_all(
                self._active_smiles, pharm_mol,
                best_params["opt_param"], best_params["n_conformers"],
                desc="Final actives"
            )
            decoy_scores = self._score_all(
                self._decoy_smiles, pharm_mol,
                best_params["opt_param"], best_params["n_conformers"],
                desc="Final decoys"
            )
            y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
            y_scores = active_scores + decoy_scores
            best_metrics = calculate_all_metrics(y_true, y_scores)

        self.best_result = {
            "best_params": best_params,
            "best_auc": best_auc,
            "best_metrics": best_metrics,
            "best_features": features,
            "n_evaluations": eval_count[0],
            "history": list(self.history),
            "elapsed_sec": elapsed,
        }

        if self.verbose:
            print("\n" + "=" * 70)
            print("OPTIMIZATION COMPLETE")
            print("=" * 70)
            print(f"Best AUC:      {best_auc:.4f}")
            print(f"Best params:   {best_params}")
            print(f"N features:    {len(features)}")
            print(f"Evaluations:   {eval_count[0]}")
            print(f"Elapsed:       {elapsed:.1f}s")
            if best_metrics:
                print(f"BEDROC(20):    {best_metrics.get('bedroc', 'N/A')}")
                print(f"EF(1%):        {best_metrics.get('ef_1', 'N/A')}")
            print("=" * 70)

        return self.best_result

    # -------------------------------------------------------------------------
    # Multi-Fidelity Optimization
    # -------------------------------------------------------------------------

    def optimize_multifidelity(
        self,
        n_calls: int = 50,
        n_random_starts: int = 15,
        tolerance_range: Tuple[float, float] = (0.5, 4.0),
        occurrence_range: Tuple[float, float] = (0.1, 1.0),
        opt_param_range: Tuple[float, float] = (0.0, 1.0),
        n_conformers_final: int = 10,
        early_stopping: bool = True,
        early_stopping_delta: float = 0.001,
        early_stopping_n: int = 10,
        explore_fraction: float = 0.7,
        refine_top_k: int = 5,
        acq_func: str = "gp_hedge",
    ) -> Dict[str, Any]:
        """Multi-fidelity Bayesian optimization with cascading evaluation.

        Two-stage approach (McDonald et al. 2025):
        1. Exploration: BO with n_conformers=1 (fast, covers more space)
        2. Refinement: Re-evaluate top-K configurations with full conformers

        Reduces wall time by 3-5x versus single-fidelity optimization while
        typically matching or exceeding final AUC.

        Args:
            n_calls: BO evaluation budget for the exploration stage.
            n_random_starts: Initial random explorations.
            tolerance_range: Search range for tolerance.
            occurrence_range: Search range for occurrence.
            opt_param_range: Search range for opt_param.
            n_conformers_final: Conformers per query in refinement stage.
            early_stopping: Enable early stopping in exploration.
            early_stopping_delta: Min improvement threshold.
            early_stopping_n: No-improvement iterations before stopping.
            explore_fraction: Fraction of budget for exploration (default: 0.7).
            refine_top_k: Top configs to refine in stage 2 (default: 5).
            acq_func: Acquisition function for BO (default: "gp_hedge").
                Use "EI" for ablation baseline.

        Returns:
            Dict with best_params, best_auc, best_metrics, best_features,
            n_evaluations, history, elapsed_sec, and stage_info with
            per-stage details.
        """
        if self._ref_smiles is None:
            raise ValueError("No data loaded.")

        # Auto-adjust budget
        explore_budget = max(3, int(n_calls * explore_fraction))
        explore_random = min(n_random_starts, max(1, explore_budget - 1))

        n_refs = len(self._ref_smiles)
        space = [
            Categorical(list(range(n_refs)), name="template_idx"),
            Real(tolerance_range[0], tolerance_range[1], name="tolerance"),
            Real(occurrence_range[0], occurrence_range[1],
                 name="occurrence_threshold"),
            Real(opt_param_range[0], opt_param_range[1], name="opt_param"),
        ]

        eval_count = [0]
        explore_results = []

        def explore_objective(params):
            template_idx, tol, occ, opt_p = params
            eval_count[0] += 1

            auc = self.evaluate(
                int(template_idx), float(tol), float(occ),
                float(opt_p), 1  # n_conformers=1 for speed
            )

            explore_results.append({
                "template_idx": int(template_idx),
                "tolerance": float(tol),
                "occurrence_threshold": float(occ),
                "opt_param": float(opt_p),
                "auc_low_fidelity": auc,
            })

            if self.verbose:
                print(
                    f"  [Explore {eval_count[0]:3d}/{explore_budget}] "
                    f"tmpl={int(template_idx)} tol={tol:.2f} "
                    f"occ={occ:.2f} opt={opt_p:.2f} "
                    f"-> AUC={auc:.4f} (1-conf)"
                )

            return -auc

        callbacks = []
        if early_stopping:
            callbacks.append(
                DeltaYStopper(
                    delta=early_stopping_delta, n_best=early_stopping_n
                )
            )

        if self.verbose:
            print("=" * 70)
            print("Multi-Fidelity Combo Pharmacophore Optimizer")
            print("=" * 70)
            print(f"References:      {n_refs}")
            print(f"Actives:         {len(self._active_smiles)}")
            print(f"Decoys:          {len(self._decoy_smiles)}")
            print(f"Stage 1 budget:  {explore_budget} evals (1 conformer)")
            print(f"Stage 2 refine:  top-{refine_top_k} "
                  f"({n_conformers_final} conformers)")
            print("=" * 70)
            print("\n--- Stage 1: Exploration ---")

        start_time = time.time()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gp_minimize(
                explore_objective,
                space,
                n_calls=explore_budget,
                n_random_starts=explore_random,
                acq_func=acq_func,
                random_state=self.random_state,
                callback=callbacks if callbacks else None,
            )

        stage1_time = time.time() - start_time

        # --- Stage 2: Refinement ---
        if self.verbose:
            print(f"\n--- Stage 2: Refinement (top-{refine_top_k}) ---")

        # Deduplicate and take top-K by low-fidelity AUC
        explore_results.sort(
            key=lambda r: r["auc_low_fidelity"], reverse=True
        )
        seen = set()
        top_configs = []
        for r in explore_results:
            key = (r["template_idx"],
                   round(r["tolerance"], 2),
                   round(r["occurrence_threshold"], 2),
                   round(r["opt_param"], 2))
            if key not in seen:
                seen.add(key)
                top_configs.append(r)
            if len(top_configs) >= refine_top_k:
                break

        best_auc = 0.0
        best_params = None
        best_metrics = {}
        best_features = []

        for i, config in enumerate(top_configs):
            auc = self.evaluate(
                config["template_idx"],
                config["tolerance"],
                config["occurrence_threshold"],
                config["opt_param"],
                n_conformers_final,
            )

            if self.verbose:
                delta = auc - config["auc_low_fidelity"]
                print(
                    f"  [Refine {i + 1}/{len(top_configs)}] "
                    f"tmpl={config['template_idx']} "
                    f"tol={config['tolerance']:.2f} "
                    f"occ={config['occurrence_threshold']:.2f} "
                    f"opt={config['opt_param']:.2f} "
                    f"-> AUC={auc:.4f} (delta={delta:+.4f})"
                )

            if auc > best_auc:
                best_auc = auc
                best_params = {
                    "template_idx": config["template_idx"],
                    "tolerance": config["tolerance"],
                    "occurrence_threshold": config["occurrence_threshold"],
                    "opt_param": config["opt_param"],
                    "n_conformers": n_conformers_final,
                }

        # Re-evaluate best to get full metrics
        if best_params is not None:
            features, pharm_mol = self._build_consensus(
                best_params["template_idx"],
                best_params["tolerance"],
                best_params["occurrence_threshold"],
            )
            best_features = features

            if pharm_mol is not None:
                active_scores = self._score_all(
                    self._active_smiles, pharm_mol,
                    best_params["opt_param"],
                    best_params["n_conformers"],
                    desc="Final actives"
                )
                decoy_scores = self._score_all(
                    self._decoy_smiles, pharm_mol,
                    best_params["opt_param"],
                    best_params["n_conformers"],
                    desc="Final decoys"
                )
                y_true = ([1] * len(active_scores)
                          + [0] * len(decoy_scores))
                y_scores = active_scores + decoy_scores
                best_metrics = calculate_all_metrics(y_true, y_scores)
        else:
            best_params = {
                "template_idx": 0,
                "tolerance": 2.0,
                "occurrence_threshold": 0.5,
                "opt_param": 0.5,
                "n_conformers": n_conformers_final,
            }

        elapsed = time.time() - start_time

        self.best_result = {
            "best_params": best_params,
            "best_auc": best_auc,
            "best_metrics": best_metrics,
            "best_features": best_features,
            "n_evaluations": eval_count[0] + len(top_configs),
            "history": list(self.history),
            "elapsed_sec": elapsed,
            "stage_info": {
                "explore_evals": eval_count[0],
                "explore_time_sec": stage1_time,
                "refine_configs": len(top_configs),
                "refine_time_sec": elapsed - stage1_time,
            },
        }

        if self.verbose:
            print("\n" + "=" * 70)
            print("MULTI-FIDELITY OPTIMIZATION COMPLETE")
            print("=" * 70)
            print(f"Best AUC:        {best_auc:.4f}")
            print(f"Best params:     {best_params}")
            print(f"N features:      {len(best_features)}")
            print(f"Explore evals:   {eval_count[0]}")
            print(f"Refine evals:    {len(top_configs)}")
            print(f"Total time:      {elapsed:.1f}s")
            print(f"  Stage 1:       {stage1_time:.1f}s")
            print(f"  Stage 2:       {elapsed - stage1_time:.1f}s")
            if best_metrics:
                print(f"BEDROC(20):      "
                      f"{best_metrics.get('bedroc', 'N/A')}")
                print(f"EF(1%):          "
                      f"{best_metrics.get('ef_1', 'N/A')}")
            print("=" * 70)

        return self.best_result

    # -------------------------------------------------------------------------
    # Export
    # -------------------------------------------------------------------------

    def export_model(
        self,
        output_dir: str,
        model_name: str = "combo_optimal",
    ) -> Dict[str, str]:
        """Export optimal model to PyMOL .pml and JSON.

        Args:
            output_dir: Directory for output files.
            model_name: Base name for files.

        Returns:
            Dict mapping format to file path.
        """
        if self.best_result is None:
            raise ValueError("No optimization result. Call optimize() first.")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        outputs = {}

        features = self.best_result["best_features"]
        best_params = self.best_result["best_params"]
        best_auc = self.best_result["best_auc"]

        # JSON export
        json_path = output_path / f"{model_name}.json"
        json_data = {
            "best_params": best_params,
            "best_auc": best_auc,
            "best_metrics": {
                k: float(v) if isinstance(v, (np.floating, float)) else v
                for k, v in self.best_result.get("best_metrics", {}).items()
            },
            "n_features": len(features),
            "features": [
                {"type": f[0], "x": float(f[2]), "y": float(f[3]), "z": float(f[4])}
                for f in features
            ],
            "dataset": {
                "n_references": len(self._ref_smiles),
                "n_actives": len(self._active_smiles),
                "n_decoys": len(self._decoy_smiles),
            },
            "reference_smiles": self._ref_smiles,
            "n_evaluations": self.best_result["n_evaluations"],
            "elapsed_sec": self.best_result["elapsed_sec"],
        }
        with open(json_path, "w") as f:
            json.dump(json_data, f, indent=2)
        outputs["json"] = str(json_path)

        # PML export
        if features:
            pml_path = output_path / f"{model_name}.pml"
            self._generate_pml(features, pml_path, model_name, best_auc)
            outputs["pml"] = str(pml_path)

        if self.verbose:
            print(f"Exported to: {outputs}")

        return outputs

    def _generate_pml(
        self, features: List, filepath: Path, name: str, auc: float
    ) -> None:
        """Generate PyMOL script for pharmacophore visualization."""
        lines = [
            f"# PyMOL script for pharmacophore: {name}",
            f"# Generated by ComboPharmacophoreOptimizer",
            f"# Best AUC: {auc:.4f}",
            "",
            "# Create pseudoatoms for pharmacophore features",
        ]

        for i, feat in enumerate(features):
            feat_type = feat[0]
            x, y, z = feat[2], feat[3], feat[4]
            color = FEATURE_COLORS.get(feat_type, (0.5, 0.5, 0.5))

            lines.append(
                f"pseudoatom {name}_{feat_type}_{i}, "
                f"pos=[{x:.3f}, {y:.3f}, {z:.3f}]"
            )
            lines.append(f"color {color}, {name}_{feat_type}_{i}")
            lines.append(f"show spheres, {name}_{feat_type}_{i}")
            lines.append(f"set sphere_scale, 1.0, {name}_{feat_type}_{i}")
            lines.append("")

        # Group all features
        feature_names = [f"{name}_{f[0]}_{i}" for i, f in enumerate(features)]
        lines.append(f"group {name}, {' '.join(feature_names)}")

        with open(filepath, "w") as f:
            f.write("\n".join(lines))

    # -------------------------------------------------------------------------
    # Utilities
    # -------------------------------------------------------------------------

    def clear_cache(self) -> None:
        """Clear all in-memory caches."""
        self._ref_conformers.clear()
        self._query_conformers.clear()
        self._aligned_refs.clear()
        self._consensus_cache.clear()

    def get_cache_stats(self) -> Dict[str, int]:
        """Return cache sizes for each level."""
        return {
            "L1_ref_conformers": len(self._ref_conformers),
            "L2_query_conformers": len(self._query_conformers),
            "L3_aligned_refs": len(self._aligned_refs),
            "L4_consensus": len(self._consensus_cache),
        }

    def __repr__(self) -> str:
        n_refs = len(self._ref_smiles) if self._ref_smiles else 0
        n_act = len(self._active_smiles) if self._active_smiles else 0
        n_dec = len(self._decoy_smiles) if self._decoy_smiles else 0
        best = f", best_auc={self.best_result['best_auc']:.4f}" if self.best_result else ""
        return (
            f"ComboPharmacophoreOptimizer("
            f"refs={n_refs}, actives={n_act}, decoys={n_dec}{best})"
        )


# =============================================================================
# CLI
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Combo Pharmacophore Optimizer - "
                    "Find optimal consensus model using rdShapeAlign combo Tanimoto"
    )
    parser.add_argument(
        "--refs",
        default="tutorials/data/CCR2_reference_ligands.sdf",
        help="Path to reference molecules SDF",
    )
    parser.add_argument(
        "--actives",
        default="tutorials/data/actives_ccr2_N75.csv",
        help="Path to actives CSV",
    )
    parser.add_argument(
        "--decoys",
        default="tutorials/data/decoys_ccr2_N500.csv",
        help="Path to decoys CSV",
    )
    parser.add_argument(
        "--n-calls", type=int, default=50,
        help="Number of Bayesian optimization evaluations (default: 50)",
    )
    parser.add_argument(
        "--n-random-starts", type=int, default=15,
        help="Number of random initial evaluations (default: 15)",
    )
    parser.add_argument(
        "--output-dir",
        default="output/combo_optimal",
        help="Output directory for model files",
    )
    parser.add_argument(
        "--no-early-stop", action="store_true",
        help="Disable early stopping",
    )
    args = parser.parse_args()

    optimizer = ComboPharmacophoreOptimizer()
    optimizer.load_from_files(args.refs, args.actives, args.decoys)

    result = optimizer.optimize(
        n_calls=args.n_calls,
        n_random_starts=args.n_random_starts,
        early_stopping=not args.no_early_stop,
    )

    outputs = optimizer.export_model(args.output_dir)

    print("\nHow to run:")
    print("  cd /home/dodo/Documents/projects/pharmacophore-toolkit")
    print("  conda activate pharmacophore-toolkit")
    print("  python -m pharmacophore.combo_optimizer \\")
    print("    --refs tutorials/data/CCR2_reference_ligands.sdf \\")
    print("    --actives tutorials/data/actives_ccr2_N75.csv \\")
    print("    --decoys tutorials/data/decoys_ccr2_N500.csv \\")
    print("    --n-calls 50 \\")
    print(f"    --output-dir {args.output_dir}")

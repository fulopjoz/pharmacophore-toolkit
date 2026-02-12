"""DrugEx reward function wrapper for pharmacophore-guided generation.

Provides a SMILES-in/float-out interface compatible with DrugEx and
other reinforcement learning molecular generators. Supports multi-fidelity
mode switching (pharm2d → 3d → hybrid) for curriculum training.

Performance estimates:
- Pharm2D mode: ~15 ms/mol → 100K molecules in 25 min (RL-compatible)
- 3D mode: ~40 ms/mol → 100K molecules in 67 min (fine-tuning only)

Example:
    >>> from pharmacophore.drugex_reward import PharmacophoreReward
    >>> reward = PharmacophoreReward(reference_mols, mode='pharm2d')
    >>> score = reward('CCO')  # Callable interface
    >>> scores = reward.score_batch(['CCO', 'c1ccccc1', 'INVALID'])
"""

from typing import List, Optional
import logging
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

logger = logging.getLogger(__name__)

_VALID_MODES = {'pharm2d', '3d', 'hybrid'}


class PharmacophoreReward:
    """SMILES-in/float-out reward for RL-based molecular generation.

    Args:
        reference_mols: Aligned reference molecules (RDKit Mol objects).
        mode: Scoring mode — 'pharm2d' (fast), '3d' (accurate), or 'hybrid'.
        n_conformers: Conformers per molecule for 3D scoring.
        normalize: If True, normalize scores to [0, 1].
        cache_size: Max SMILES to cache (LRU).
    """

    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        mode: str = 'pharm2d',
        n_conformers: int = 5,
        normalize: bool = True,
        cache_size: int = 10000,
    ):
        if mode not in _VALID_MODES:
            raise ValueError(
                f"mode must be one of {_VALID_MODES}, got '{mode}'"
            )

        self.reference_mols = reference_mols
        self.mode = mode
        self.n_conformers = n_conformers
        self.normalize = normalize
        self._cache_size = cache_size

        # Stats tracking
        self._n_calls = 0
        self._n_valid = 0
        self._n_errors = 0
        self._n_cache_hits = 0

        # Lazy-initialized scorers
        self._pharm2d_scorer = None
        self._ref_ensemble_scorer = None

        # LRU cache for scoring
        self._score_cached = lru_cache(maxsize=cache_size)(self._score_uncached)

    def _get_pharm2d_scorer(self):
        """Lazy-initialize Pharm2DScorer."""
        if self._pharm2d_scorer is None:
            from .pharm2d_scoring import Pharm2DScorer
            self._pharm2d_scorer = Pharm2DScorer(self.reference_mols)
        return self._pharm2d_scorer

    def _get_ref_ensemble_scorer(self):
        """Lazy-initialize ReferenceEnsembleScorer."""
        if self._ref_ensemble_scorer is None:
            from .rdshape_optimizer import ReferenceEnsembleScorer
            self._ref_ensemble_scorer = ReferenceEnsembleScorer(
                self.reference_mols, n_conformers=self.n_conformers
            )
        return self._ref_ensemble_scorer

    def _score_uncached(self, smiles: str) -> float:
        """Score a single SMILES string (uncached).

        Returns:
            Score in [0, 1] if normalize=True, raw score otherwise.
            0.0 for invalid SMILES.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            self._n_errors += 1
            return 0.0

        self._n_valid += 1

        try:
            if self.mode == 'pharm2d':
                raw = self._get_pharm2d_scorer().score(mol)
                # Pharm2D Tanimoto is already in [0, 1]
                return float(raw)

            elif self.mode == '3d':
                # Need 3D conformer
                mol_h = Chem.AddHs(mol)
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                params.numThreads = 1
                AllChem.EmbedMultipleConfs(
                    mol_h, numConfs=self.n_conformers, params=params
                )
                if mol_h.GetNumConformers() == 0:
                    return 0.0

                scorer = self._get_ref_ensemble_scorer()
                raw = scorer.score(mol_h)
                if self.normalize:
                    return float(min(raw / 2.0, 1.0))  # [0,2] → [0,1]
                return float(raw)

            elif self.mode == 'hybrid':
                # Pharm2D component
                score_2d = self._get_pharm2d_scorer().score(mol)

                # 3D component
                mol_h = Chem.AddHs(mol)
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                params.numThreads = 1
                AllChem.EmbedMultipleConfs(
                    mol_h, numConfs=self.n_conformers, params=params
                )
                if mol_h.GetNumConformers() == 0:
                    return float(score_2d)  # Fall back to 2D only

                scorer = self._get_ref_ensemble_scorer()
                raw_3d = scorer.score(mol_h)
                score_3d = min(raw_3d / 2.0, 1.0)

                # 60/40 blend (same as evaluation.py default alpha=0.6)
                return float(0.6 * score_2d + 0.4 * score_3d)

        except Exception as e:
            logger.debug("Scoring error for %s: %s", smiles, e)
            self._n_errors += 1
            return 0.0

    def score(self, smiles: str) -> float:
        """Score a single SMILES string.

        Args:
            smiles: SMILES string.

        Returns:
            Score in [0, 1]. 0.0 for invalid SMILES.
        """
        self._n_calls += 1
        # Check cache info to track hits
        info_before = self._score_cached.cache_info()
        result = self._score_cached(smiles)
        info_after = self._score_cached.cache_info()
        if info_after.hits > info_before.hits:
            self._n_cache_hits += 1
        return result

    def score_batch(self, smiles_list: List[str]) -> List[float]:
        """Score a batch of SMILES strings.

        Args:
            smiles_list: List of SMILES strings.

        Returns:
            List of scores.
        """
        return [self.score(smi) for smi in smiles_list]

    def switch_mode(self, new_mode: str) -> None:
        """Switch scoring mode (multi-fidelity curriculum).

        Clears the score cache since scores change between modes.

        Args:
            new_mode: One of 'pharm2d', '3d', 'hybrid'.

        Raises:
            ValueError: If new_mode is invalid.
        """
        if new_mode not in _VALID_MODES:
            raise ValueError(
                f"mode must be one of {_VALID_MODES}, got '{new_mode}'"
            )
        if new_mode != self.mode:
            self.mode = new_mode
            self._score_cached.cache_clear()
            logger.info("Switched to mode '%s', cache cleared", new_mode)

    def get_stats(self) -> dict:
        """Get scoring statistics.

        Returns:
            Dict with n_calls, n_valid, n_errors, cache_size,
            cache_hit_rate.
        """
        cache_info = self._score_cached.cache_info()
        return {
            'n_calls': self._n_calls,
            'n_valid': self._n_valid,
            'n_errors': self._n_errors,
            'cache_size': cache_info.currsize,
            'cache_hits': self._n_cache_hits,
            'cache_hit_rate': (
                self._n_cache_hits / self._n_calls
                if self._n_calls > 0 else 0.0
            ),
        }

    def __call__(self, smiles: str) -> float:
        """Callable interface: ``reward('CCO')``."""
        return self.score(smiles)

    def __repr__(self) -> str:
        return (
            f"PharmacophoreReward(mode='{self.mode}', "
            f"n_refs={len(self.reference_mols)}, "
            f"n_calls={self._n_calls})"
        )

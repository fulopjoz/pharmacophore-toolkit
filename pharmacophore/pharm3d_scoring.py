"""Pharm3D Distance-based scoring for pharmacophore virtual screening.

This module provides 3D pharmacophore distance scoring using consensus features.
Unlike shape alignment (which fails with disconnected fragments), this approach
scores molecules by how well their feature pairwise distances match the
consensus pharmacophore pairwise distances.

Key advantages:
- Coordinate-system independent (no alignment needed)
- Captures spatial constraints between features
- Works with consensus pharmacophore features
- Achieves AUC ~0.84 vs ~0.47 for shape-only methods

Example:
    >>> from pharmacophore.pharm3d_scoring import Pharm3DScorer
    >>> from pharmacophore import PharmacophoreConsensus
    >>>
    >>> # Generate consensus
    >>> consensus = PharmacophoreConsensus(tolerance=1.0, occurrence_threshold=0.5)
    >>> features = consensus.generate_consensus(reference_mols)
    >>>
    >>> # Create scorer
    >>> scorer = Pharm3DScorer(features, distance_tolerance=1.5)
    >>>
    >>> # Score molecules
    >>> score = scorer.score(query_mol)
"""

from typing import List, Optional, Tuple, Dict
from collections import defaultdict
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .pharmacophore import Pharmacophore


class Pharm3DScorer:
    """Score molecules by matching pairwise feature distances to consensus.

    This scorer extracts pharmacophore features from query molecules and
    compares their pairwise distances to the expected distances from the
    consensus pharmacophore. Features that match at similar distances
    receive high scores.

    The scoring function uses a Gaussian:
        score = exp(-diff²/(2σ²))
    where diff is the difference between actual and expected distance,
    and σ is the distance_tolerance parameter.

    Attributes:
        consensus_features: List of consensus feature definitions.
        distance_tolerance: Gaussian width for distance matching (Angstroms).
        feature_pairs: List of (type1, type2, distance) tuples.
        required_types: Set of feature types in the consensus.

    Example:
        >>> features = [
        ...     ['Acceptor', (), 1.0, 2.0, 3.0],
        ...     ['Donor', (), 4.0, 5.0, 6.0],
        ... ]
        >>> scorer = Pharm3DScorer(features, distance_tolerance=1.5)
        >>> print(f"Feature pairs: {len(scorer.feature_pairs)}")
        Feature pairs: 1
    """

    def __init__(
        self,
        consensus_features: List[List],
        distance_tolerance: float = 1.5
    ):
        """Initialize scorer with consensus features.

        Args:
            consensus_features: List of consensus features, where each feature
                is [type, atom_indices, x, y, z]. Typically from
                PharmacophoreConsensus.generate_consensus().
            distance_tolerance: Gaussian width for distance matching in
                Angstroms. Smaller values are stricter. Default: 1.5 Å.

        Raises:
            ValueError: If consensus_features is empty or distance_tolerance <= 0.
        """
        if not consensus_features:
            raise ValueError("consensus_features cannot be empty")
        if distance_tolerance <= 0:
            raise ValueError(f"distance_tolerance must be > 0, got {distance_tolerance}")

        self.consensus_features = consensus_features
        self.distance_tolerance = distance_tolerance
        self.feature_pairs: List[Tuple[str, str, float]] = []
        self.required_types: set = set()

        if len(consensus_features) < 2:
            # Cannot compute pairwise distances with fewer than 2 features
            return

        # Extract positions and types
        positions = []
        types = []
        for f in consensus_features:
            positions.append([f[2], f[3], f[4]])
            types.append(f[0])
            self.required_types.add(f[0])

        positions = np.array(positions)
        n_features = len(positions)

        # Compute all pairwise distances
        for i in range(n_features):
            for j in range(i + 1, n_features):
                dist = np.linalg.norm(positions[i] - positions[j])
                pair_type = tuple(sorted([types[i], types[j]]))
                self.feature_pairs.append((pair_type[0], pair_type[1], dist))

    def score(self, mol: Chem.Mol) -> float:
        """Score a molecule against the consensus pharmacophore.

        For each consensus feature pair, finds the best matching feature
        pair in the query molecule and scores by distance similarity.

        Args:
            mol: RDKit Mol object with 3D conformer.

        Returns:
            Score between 0.0 and 1.0. Higher = better match.
            Returns 0.0 if molecule is None, has no conformer, or
            has insufficient pharmacophore features.

        Example:
            >>> mol = Chem.MolFromSmiles('CCO')
            >>> AllChem.EmbedMolecule(mol)
            >>> score = scorer.score(mol)
            >>> print(f"Score: {score:.3f}")
        """
        if mol is None or len(self.feature_pairs) == 0:
            return 0.0

        # Extract pharmacophore features from query molecule
        pharm = Pharmacophore()
        try:
            query_features = pharm.calc_pharm(mol)
        except Exception:
            return 0.0

        if not query_features or len(query_features) < 2:
            return 0.0

        # Organize query features by type
        query_by_type: Dict[str, List[np.ndarray]] = defaultdict(list)
        for f in query_features:
            query_by_type[f[0]].append(np.array([f[2], f[3], f[4]]))

        # Score each consensus pair
        pair_scores = []

        for type1, type2, expected_dist in self.feature_pairs:
            if type1 not in query_by_type or type2 not in query_by_type:
                # Missing feature type - score 0
                pair_scores.append(0.0)
                continue

            if type1 == type2:
                # Same type: pairs within same list
                positions = query_by_type[type1]
                if len(positions) < 2:
                    pair_scores.append(0.0)
                    continue

                min_diff = float('inf')
                for i, p1 in enumerate(positions):
                    for j, p2 in enumerate(positions):
                        if i >= j:
                            continue
                        d = np.linalg.norm(p1 - p2)
                        diff = abs(d - expected_dist)
                        min_diff = min(min_diff, diff)
            else:
                # Different types: pairs across lists
                min_diff = float('inf')
                for p1 in query_by_type[type1]:
                    for p2 in query_by_type[type2]:
                        d = np.linalg.norm(p1 - p2)
                        diff = abs(d - expected_dist)
                        min_diff = min(min_diff, diff)

            if min_diff == float('inf'):
                pair_scores.append(0.0)
            else:
                # Gaussian scoring
                score = np.exp(-min_diff**2 / (2 * self.distance_tolerance**2))
                pair_scores.append(score)

        return float(np.mean(pair_scores)) if pair_scores else 0.0

    def score_with_conformers(
        self,
        mol: Chem.Mol,
        n_conformers: int = 10,
        random_seed: int = 42
    ) -> float:
        """Score molecule using multiple conformers.

        Generates conformers for the molecule and returns the best score
        across all conformers.

        Args:
            mol: RDKit Mol object (can be 2D, will generate 3D conformers).
            n_conformers: Number of conformers to generate. Default: 10.
            random_seed: Random seed for conformer generation. Default: 42.

        Returns:
            Best score across all conformers.

        Example:
            >>> mol = Chem.MolFromSmiles('CCO')
            >>> score = scorer.score_with_conformers(mol, n_conformers=5)
        """
        if mol is None:
            return 0.0

        mol_h = Chem.AddHs(mol)

        # Generate conformers
        params = AllChem.ETKDGv3()
        params.randomSeed = random_seed
        AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)

        if mol_h.GetNumConformers() == 0:
            # Fallback: try single conformer
            AllChem.EmbedMolecule(mol_h, randomSeed=random_seed)

        if mol_h.GetNumConformers() == 0:
            return 0.0

        # Score each conformer
        best_score = 0.0
        for conf_id in range(mol_h.GetNumConformers()):
            conf_mol = Chem.Mol(mol_h)
            conf_mol.RemoveAllConformers()
            conf_mol.AddConformer(mol_h.GetConformer(conf_id), assignId=True)
            score = self.score(conf_mol)
            best_score = max(best_score, score)

        return best_score

    def score_all(
        self,
        mols: List[Chem.Mol],
        n_conformers: int = 10,
        random_seed: int = 42
    ) -> np.ndarray:
        """Score multiple molecules.

        Args:
            mols: List of RDKit Mol objects.
            n_conformers: Number of conformers per molecule. Default: 10.
            random_seed: Random seed for conformer generation. Default: 42.

        Returns:
            numpy array of scores, shape (len(mols),).

        Example:
            >>> scores = scorer.score_all(query_mols, n_conformers=5)
            >>> print(f"Mean score: {scores.mean():.3f}")
        """
        scores = np.zeros(len(mols))
        for i, mol in enumerate(mols):
            scores[i] = self.score_with_conformers(mol, n_conformers, random_seed)
        return scores

    def get_feature_pair_summary(self) -> Dict[Tuple[str, str], List[float]]:
        """Get summary of feature pair distances.

        Returns:
            Dictionary mapping (type1, type2) to list of distances.

        Example:
            >>> summary = scorer.get_feature_pair_summary()
            >>> for pair, dists in summary.items():
            ...     print(f"{pair}: {min(dists):.1f}-{max(dists):.1f} Å")
        """
        summary: Dict[Tuple[str, str], List[float]] = defaultdict(list)
        for type1, type2, dist in self.feature_pairs:
            summary[(type1, type2)].append(dist)
        return dict(summary)

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"Pharm3DScorer(n_features={len(self.consensus_features)}, "
            f"n_pairs={len(self.feature_pairs)}, "
            f"distance_tolerance={self.distance_tolerance})"
        )


def score_molecules_pharm3d(
    reference_mols: List[Chem.Mol],
    query_mols: List[Chem.Mol],
    tolerance: float = 1.0,
    occurrence_threshold: float = 0.5,
    distance_tolerance: float = 1.5,
    n_conformers: int = 10
) -> np.ndarray:
    """Convenience function for Pharm3D distance scoring.

    Generates consensus from reference molecules and scores query molecules.

    Args:
        reference_mols: Aligned reference molecules for consensus.
        query_mols: Molecules to score.
        tolerance: Consensus clustering tolerance in Angstroms. Default: 1.0.
        occurrence_threshold: Minimum fraction of references for feature.
            Default: 0.5.
        distance_tolerance: Gaussian width for distance matching. Default: 1.5.
        n_conformers: Conformers per query molecule. Default: 10.

    Returns:
        numpy array of scores for each query molecule.

    Example:
        >>> scores = score_molecules_pharm3d(refs, queries)
        >>> print(f"Top scorer index: {scores.argmax()}")
    """
    from .consensus import PharmacophoreConsensus

    # Generate consensus
    consensus = PharmacophoreConsensus(
        tolerance=tolerance,
        occurrence_threshold=occurrence_threshold
    )
    features = consensus.generate_consensus(reference_mols)

    if len(features) < 2:
        return np.zeros(len(query_mols))

    # Create scorer
    scorer = Pharm3DScorer(features, distance_tolerance=distance_tolerance)

    # Score molecules
    return scorer.score_all(query_mols, n_conformers=n_conformers)

"""Stability-aware ensemble consensus pharmacophore generation.

Runs consensus clustering multiple times with perturbed parameters and/or
different algorithms, then votes on which features are stable across runs.

Features that appear consistently (within spatial hash tolerance) across
many runs are assigned stability scores. This produces more robust consensus
models than single-run clustering.

Source: Li et al. 2018 (ensemble k-medoids), Roy et al. 2018 (intelligent consensus).
"""

import numpy as np
from typing import List, Dict, Optional, Tuple
from rdkit import Chem

from .consensus import PharmacophoreConsensus
from .clustering_algorithms import cluster_features_generic


class EnsembleConsensus:
    """Generate stability-aware consensus pharmacophores via ensemble voting.

    Runs consensus clustering N times with perturbed parameters and/or
    different clustering algorithms. Features are matched across runs by
    spatial proximity (hash_tolerance). Each feature receives a stability
    score = fraction of runs in which it appeared.

    Attributes:
        n_runs: Number of clustering runs.
        tolerance_range: (min, max) range for tolerance sampling.
        occurrence_range: (min, max) range for occurrence threshold sampling.
        methods: List of clustering methods to use for diversity.
        hash_tolerance: Spatial distance to match features across runs (Angstroms).
        stability_threshold: Minimum stability to include a feature.
        random_state: Seed for reproducibility.
    """

    def __init__(
        self,
        n_runs: int = 25,
        tolerance_range: Tuple[float, float] = (1.5, 2.5),
        occurrence_range: Tuple[float, float] = (0.3, 0.7),
        methods: Optional[List[str]] = None,
        hash_tolerance: float = 0.5,
        stability_threshold: float = 0.3,
        random_state: int = 42
    ):
        """Initialize ensemble consensus generator.

        Args:
            n_runs: Number of clustering runs (default: 25).
            tolerance_range: (min, max) for tolerance perturbation.
            occurrence_range: (min, max) for occurrence threshold perturbation.
            methods: Clustering methods to alternate between. Default:
                ['hierarchical', 'kmeans']. Options: 'hierarchical',
                'kmeans', 'grid'.
            hash_tolerance: Distance threshold to match features across runs.
            stability_threshold: Minimum fraction of runs for inclusion.
            random_state: Random seed for parameter sampling.
        """
        self.n_runs = n_runs
        self.tolerance_range = tolerance_range
        self.occurrence_range = occurrence_range
        self.methods = methods or ['hierarchical', 'kmeans']
        if 'dbscan' in self.methods:
            import warnings
            warnings.warn(
                "DBSCAN removed for non-determinism. Replacing with 'hierarchical'.",
                DeprecationWarning, stacklevel=2
            )
            self.methods = [m if m != 'dbscan' else 'hierarchical' for m in self.methods]
        self.hash_tolerance = hash_tolerance
        self.stability_threshold = stability_threshold
        self.rng = np.random.RandomState(random_state)

    def generate_consensus(
        self,
        mols: List[Chem.Mol],
        features: Optional[str] = None
    ) -> List[List]:
        """Generate stability-aware consensus pharmacophore.

        Runs clustering n_runs times, collects all candidate features,
        then votes on which features are stable.

        Args:
            mols: Aligned molecules with 3D conformers.
            features: Feature definition ('default', 'rdkit', or custom).

        Returns:
            List of consensus features [type, (), x, y, z] with only stable
            features included. Features are sorted by stability (highest first).
        """
        # Collect candidate features from all runs
        all_run_features = []

        for run_idx in range(self.n_runs):
            run_features = self._single_run(mols, features, run_idx)
            all_run_features.append(run_features)

        # Vote on feature stability
        stable_features = self._vote_features(all_run_features)

        return stable_features

    def generate_consensus_with_scores(
        self,
        mols: List[Chem.Mol],
        features: Optional[str] = None
    ) -> Tuple[List[List], List[float]]:
        """Generate consensus with stability scores per feature.

        Args:
            mols: Aligned molecules with 3D conformers.
            features: Feature definition.

        Returns:
            Tuple of (features, stability_scores).
            stability_scores[i] is the fraction of runs in which features[i]
            appeared (0.0 to 1.0).
        """
        all_run_features = []

        for run_idx in range(self.n_runs):
            run_features = self._single_run(mols, features, run_idx)
            all_run_features.append(run_features)

        stable_features, scores = self._vote_features(
            all_run_features, return_scores=True
        )

        return stable_features, scores

    def _single_run(
        self,
        mols: List[Chem.Mol],
        features: Optional[str],
        run_idx: int
    ) -> List[List]:
        """Execute a single consensus run with perturbed parameters.

        Alternates between clustering methods and samples parameters
        from the configured ranges.
        """
        # Sample parameters
        tol = self.rng.uniform(*self.tolerance_range)
        occ = self.rng.uniform(*self.occurrence_range)

        # Alternate methods
        method = self.methods[run_idx % len(self.methods)]

        if method == 'hierarchical':
            # Use the standard PharmacophoreConsensus
            consensus = PharmacophoreConsensus(
                tolerance=tol,
                occurrence_threshold=occ,
                linkage='average'
            )
            return consensus.generate_consensus(mols, features=features)
        else:
            # Use alternative clustering from clustering_algorithms.py
            return self._alternative_consensus(
                mols, features, tol, occ, method
            )

    def _alternative_consensus(
        self,
        mols: List[Chem.Mol],
        features: Optional[str],
        tolerance: float,
        occurrence_threshold: float,
        method: str
    ) -> List[List]:
        """Run consensus with an alternative clustering algorithm."""
        from pharmacophore.pharmacophore import Pharmacophore

        if features is None:
            features = 'default'

        pharm = Pharmacophore(features=features)
        all_features = []

        for mol in mols:
            if not mol.GetNumConformers():
                continue
            mol_features = pharm.calc_pharm(mol=mol, features=features)
            all_features.append(mol_features)

        # Group by type
        features_by_type: Dict[str, List[np.ndarray]] = {}
        for mol_features in all_features:
            for feature in mol_features:
                feat_type = feature[0]
                if feat_type not in features_by_type:
                    features_by_type[feat_type] = []
                coords = np.array([feature[2], feature[3], feature[4]])
                features_by_type[feat_type].append(coords)

        # Cluster each type using alternative method
        consensus_features = []
        n_mols = len(mols)

        for feat_type, coords_list in features_by_type.items():
            coords_array = np.array(coords_list)
            centroids = cluster_features_generic(
                coords=coords_array,
                tolerance=tolerance,
                occurrence_threshold=occurrence_threshold,
                n_molecules=n_mols,
                method=method
            )

            for centroid in centroids:
                consensus_features.append([
                    feat_type,
                    (),
                    float(centroid[0]),
                    float(centroid[1]),
                    float(centroid[2])
                ])

        return consensus_features

    def _vote_features(
        self,
        all_run_features: List[List[List]],
        return_scores: bool = False
    ):
        """Vote on feature stability across runs.

        Uses spatial hashing to match features across runs. Two features
        from different runs are considered the same if they have the same
        type and are within hash_tolerance Angstroms of each other.

        Args:
            all_run_features: List of feature lists from each run.
            return_scores: If True, also return stability scores.

        Returns:
            If return_scores: (features, scores)
            Else: features
        """
        if not all_run_features:
            return ([], []) if return_scores else []

        # Collect all unique feature candidates
        # Each candidate: (type, x, y, z, vote_count)
        candidates = []

        for run_features in all_run_features:
            for feature in run_features:
                feat_type = feature[0]
                x, y, z = feature[2], feature[3], feature[4]

                # Try to match with existing candidate
                matched = False
                for cand in candidates:
                    if cand['type'] != feat_type:
                        continue
                    dist = np.sqrt(
                        (cand['x'] - x) ** 2 +
                        (cand['y'] - y) ** 2 +
                        (cand['z'] - z) ** 2
                    )
                    if dist <= self.hash_tolerance:
                        # Update running average position
                        n = cand['votes']
                        cand['x'] = (cand['x'] * n + x) / (n + 1)
                        cand['y'] = (cand['y'] * n + y) / (n + 1)
                        cand['z'] = (cand['z'] * n + z) / (n + 1)
                        cand['votes'] += 1
                        matched = True
                        break

                if not matched:
                    candidates.append({
                        'type': feat_type,
                        'x': x, 'y': y, 'z': z,
                        'votes': 1
                    })

        # Filter by stability threshold and build output
        n_runs = len(all_run_features)
        stable = []
        scores = []

        for cand in candidates:
            stability = cand['votes'] / n_runs
            if stability >= self.stability_threshold:
                stable.append((stability, cand))

        # Sort by stability descending
        stable.sort(key=lambda x: -x[0])

        features_out = []
        scores_out = []

        for stability, cand in stable:
            features_out.append([
                cand['type'],
                (),
                float(cand['x']),
                float(cand['y']),
                float(cand['z'])
            ])
            scores_out.append(stability)

        if return_scores:
            return features_out, scores_out
        return features_out

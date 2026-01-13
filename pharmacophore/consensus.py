"""Consensus pharmacophore generation module.

This module provides deterministic consensus pharmacophore generation
using agglomerative hierarchical clustering with MOE-style parameters.

The implementation uses scikit-learn's AgglomerativeClustering with
average linkage to ensure deterministic, explainable results.
"""

import warnings
import numpy as np
from typing import List, Tuple, Dict, Optional, Union
from sklearn.cluster import AgglomerativeClustering
from rdkit import Chem


class PharmacophoreConsensus:
    """Generate consensus pharmacophores from aligned molecules.
    
    This class uses agglomerative hierarchical clustering to identify
    common pharmacophore features across multiple aligned molecules.
    Features are clustered by spatial proximity and filtered by occurrence
    threshold.
    
    Attributes:
        tolerance: Distance threshold in Angstroms for feature merging.
        occurrence_threshold: Minimum fraction (0.0-1.0) of molecules that
            must contain a feature for it to be included in consensus.
        linkage: Clustering linkage method ('average', 'complete', 'single',
            'ward').
    
    Example:
        >>> from rdkit import Chem
        >>> from pharmacophore.consensus import PharmacophoreConsensus
        >>> 
        >>> # Load aligned molecules
        >>> mols = [Chem.SDMolSupplier('mol1.sdf')[0],
        ...         Chem.SDMolSupplier('mol2.sdf')[0]]
        >>> 
        >>> # Generate consensus
        >>> consensus = PharmacophoreConsensus(
        ...     tolerance=2.0,
        ...     occurrence_threshold=0.5
        ... )
        >>> features = consensus.generate_consensus(mols)
    """
    
    def __init__(
        self,
        tolerance: float = 2.0,
        occurrence_threshold: float = 0.5,
        linkage: str = 'average'
    ):
        """Initialize consensus pharmacophore generator.
        
        Args:
            tolerance: Distance threshold in Angstroms for clustering
                features (default: 2.0).
            occurrence_threshold: Minimum fraction of molecules (0.0-1.0)
                that must contain a feature (default: 0.5).
            linkage: Clustering method - 'average', 'complete', 'single',
                or 'ward' (default: 'average').
        
        Raises:
            ValueError: If tolerance <= 0, occurrence_threshold not in
                [0.0, 1.0], or linkage not valid.
        """
        if tolerance <= 0:
            raise ValueError(
                f"tolerance must be > 0, got {tolerance}"
            )
        
        if not 0.0 <= occurrence_threshold <= 1.0:
            raise ValueError(
                f"occurrence_threshold must be in [0.0, 1.0], "
                f"got {occurrence_threshold}"
            )
        
        valid_linkages = {'average', 'complete', 'single', 'ward'}
        if linkage not in valid_linkages:
            raise ValueError(
                f"linkage must be one of {valid_linkages}, got {linkage}"
            )
        
        self.tolerance = tolerance
        self.occurrence_threshold = occurrence_threshold
        self.linkage = linkage
    
    def generate_consensus(
        self,
        mols: List[Chem.Mol],
        features: Optional[Union[str, dict]] = None
    ) -> List[List]:
        """Generate consensus pharmacophore from aligned molecules.
        
        This method:
        1. Extracts pharmacophore features from each molecule
        2. Groups features by type
        3. Clusters features of each type by spatial proximity
        4. Filters clusters by occurrence threshold
        5. Computes consensus feature centroids
        
        Args:
            mols: List of aligned molecules with 3D conformations.
            features: Feature type to use ('default', 'rdkit', or custom
                dict). If None, uses 'default'.
        
        Returns:
            List of consensus features, where each feature is:
                [type, atom_indices, x, y, z]
            atom_indices is empty tuple () since consensus features don't
            map to specific atoms.
        
        Raises:
            ValueError: If mols is empty or molecules lack 3D conformers.
        
        Example:
            >>> consensus = PharmacophoreConsensus(tolerance=2.0)
            >>> features = consensus.generate_consensus(aligned_mols)
            >>> print(f"Found {len(features)} consensus features")
        """
        if not mols:
            raise ValueError("mols list cannot be empty")
        
        # Import here to avoid circular dependency
        from pharmacophore.pharmacophore import Pharmacophore
        
        # Set default features
        if features is None:
            features = 'default'
        
        # Extract pharmacophore features from all molecules
        pharm = Pharmacophore(features=features)
        all_features = []
        
        for mol_idx, mol in enumerate(mols):
            # Validate 3D conformer
            if not mol.GetNumConformers():
                raise ValueError(
                    f"Molecule {mol_idx} has no conformer. "
                    f"Ensure all molecules have 3D coordinates."
                )
            
            mol_features = pharm.calc_pharm(mol=mol, features=features)
            all_features.append(mol_features)
        
        # Group features by type
        features_by_type = {}
        for mol_idx, mol_features in enumerate(all_features):
            for feature in mol_features:
                feat_type = feature[0]
                if feat_type not in features_by_type:
                    features_by_type[feat_type] = []
                
                # Store: (coordinates, molecule_index)
                coords = np.array([feature[2], feature[3], feature[4]])
                features_by_type[feat_type].append((coords, mol_idx))
        
        # Cluster and filter features of each type
        consensus_features = []
        total_mols = len(mols)
        
        for feat_type, type_features in features_by_type.items():
            # Extract coordinates and molecule indices
            coords_list = [f[0] for f in type_features]
            mol_indices = np.array([f[1] for f in type_features])
            
            # Cluster features
            labels, cluster_to_mols = self._cluster_features(
                coordinates=coords_list,
                mol_indices=mol_indices
            )
            
            # Calculate centroids for each cluster
            centroids = self._calculate_cluster_centroids(
                coordinates=coords_list,
                labels=labels
            )
            
            # Filter by occurrence threshold
            valid_features = self._filter_by_occurrence(
                centroids=centroids,
                cluster_to_mols=cluster_to_mols,
                total_molecules=total_mols
            )
            
            # Create consensus features
            for centroid, num_mols in valid_features:
                # Format: [type, atom_indices, x, y, z]
                consensus_feature = [
                    feat_type,
                    (),  # Empty tuple - no atom mapping
                    float(centroid[0]),
                    float(centroid[1]),
                    float(centroid[2])
                ]
                consensus_features.append(consensus_feature)
        
        return consensus_features
    
    def _cluster_features(
        self,
        coordinates: List[np.ndarray],
        mol_indices: np.ndarray
    ) -> Tuple[np.ndarray, Dict[int, List[int]]]:
        """Cluster features using agglomerative hierarchical clustering.
        
        Uses scikit-learn's AgglomerativeClustering with distance_threshold
        to automatically determine the number of clusters.
        
        Args:
            coordinates: List of 3D coordinate arrays.
            mol_indices: Array of molecule indices for each feature.
        
        Returns:
            Tuple of:
                - labels: Cluster labels for each feature
                - cluster_to_mols: Dict mapping cluster ID to list of
                    unique molecule indices
        
        Example:
            >>> coords = [np.array([1, 2, 3]), np.array([1.5, 2.1, 3.2])]
            >>> mol_idx = np.array([0, 1])
            >>> labels, mapping = self._cluster_features(coords, mol_idx)
        """
        if not coordinates:
            return np.array([]), {}
        
        # Convert to numpy array
        coords_array = np.array(coordinates)
        
        # Single feature case
        if len(coords_array) == 1:
            return np.array([0]), {0: [int(mol_indices[0])]}
        
        # Perform agglomerative clustering
        clustering = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=self.tolerance,
            linkage=self.linkage,
            metric='euclidean'
        )
        
        labels = clustering.fit_predict(coords_array)
        
        # Map clusters to molecule indices
        cluster_to_mols = {}
        for cluster_id, mol_idx in zip(labels, mol_indices):
            cluster_id = int(cluster_id)
            mol_idx = int(mol_idx)
            
            if cluster_id not in cluster_to_mols:
                cluster_to_mols[cluster_id] = []
            
            # Store unique molecule indices
            if mol_idx not in cluster_to_mols[cluster_id]:
                cluster_to_mols[cluster_id].append(mol_idx)
        
        return labels, cluster_to_mols
    
    def _calculate_cluster_centroids(
        self,
        coordinates: List[np.ndarray],
        labels: np.ndarray
    ) -> Dict[int, np.ndarray]:
        """Calculate centroid for each cluster.
        
        Args:
            coordinates: List of 3D coordinate arrays.
            labels: Cluster labels for each coordinate.
        
        Returns:
            Dict mapping cluster ID to centroid coordinates.
        
        Example:
            >>> coords = [np.array([1, 2, 3]), np.array([1.5, 2.1, 3.2])]
            >>> labels = np.array([0, 0])
            >>> centroids = self._calculate_cluster_centroids(coords, labels)
            >>> print(centroids[0])  # Average of the two points
        """
        centroids = {}
        coords_array = np.array(coordinates)
        
        for cluster_id in np.unique(labels):
            cluster_id = int(cluster_id)
            cluster_coords = coords_array[labels == cluster_id]
            centroid = np.mean(cluster_coords, axis=0)
            centroids[cluster_id] = centroid
        
        return centroids
    
    def _filter_by_occurrence(
        self,
        centroids: Dict[int, np.ndarray],
        cluster_to_mols: Dict[int, List[int]],
        total_molecules: int
    ) -> List[Tuple[np.ndarray, int]]:
        """Filter clusters by occurrence threshold.
        
        Only keeps clusters where the feature appears in at least
        occurrence_threshold fraction of molecules.
        
        Args:
            centroids: Dict mapping cluster ID to centroid coordinates.
            cluster_to_mols: Dict mapping cluster ID to molecule indices.
            total_molecules: Total number of molecules in the dataset.
        
        Returns:
            List of tuples (centroid, num_molecules) for valid features.
        
        Example:
            >>> # Feature must appear in at least 50% of molecules
            >>> consensus = PharmacophoreConsensus(occurrence_threshold=0.5)
            >>> valid = consensus._filter_by_occurrence(
            ...     centroids={0: np.array([1, 2, 3])},
            ...     cluster_to_mols={0: [0, 1, 2]},
            ...     total_molecules=5
            ... )
            >>> print(len(valid))  # 1 if 3/5 >= 0.5, else 0
        """
        valid_features = []
        min_molecules = int(np.ceil(self.occurrence_threshold * total_molecules))
        
        for cluster_id, centroid in centroids.items():
            num_molecules = len(cluster_to_mols[cluster_id])
            
            if num_molecules >= min_molecules:
                valid_features.append((centroid, num_molecules))
        
        return valid_features
    
    def generate_models(
        self,
        mols: List[Chem.Mol],
        model_set: str = 'standard',
        features: Optional[Union[str, dict]] = None
    ) -> Dict[str, List[List]]:
        """Generate multiple consensus pharmacophore models.
        
        Creates a set of models with different stringency levels to give
        users options for different use cases.
        
        Args:
            mols: List of aligned molecules with 3D conformations.
            model_set: Model set to generate - 'standard', 'comprehensive',
                or 'custom' (default: 'standard').
            features: Feature type to use (default: None, uses 'default').
        
        Returns:
            Dict mapping model name to consensus features list.
            Standard set includes: 'strict', 'moderate', 'relaxed'
        
        Raises:
            ValueError: If model_set is not valid.
        
        Example:
            >>> consensus = PharmacophoreConsensus()
            >>> models = consensus.generate_models(aligned_mols)
            >>> print(f"Strict: {len(models['strict'])} features")
            >>> print(f"Moderate: {len(models['moderate'])} features")
            >>> print(f"Relaxed: {len(models['relaxed'])} features")
        """
        if model_set == 'standard':
            configs = {
                'strict': {
                    'tolerance': self.tolerance * 0.75,
                    'occurrence_threshold': 0.8
                },
                'moderate': {
                    'tolerance': self.tolerance,
                    'occurrence_threshold': 0.5
                },
                'relaxed': {
                    'tolerance': self.tolerance * 1.5,
                    'occurrence_threshold': 0.3
                }
            }
        elif model_set == 'comprehensive':
            configs = {
                'very_strict': {
                    'tolerance': self.tolerance * 0.5,
                    'occurrence_threshold': 0.9
                },
                'strict': {
                    'tolerance': self.tolerance * 0.75,
                    'occurrence_threshold': 0.8
                },
                'moderate': {
                    'tolerance': self.tolerance,
                    'occurrence_threshold': 0.5
                },
                'relaxed': {
                    'tolerance': self.tolerance * 1.5,
                    'occurrence_threshold': 0.3
                },
                'very_relaxed': {
                    'tolerance': self.tolerance * 2.0,
                    'occurrence_threshold': 0.2
                }
            }
        else:
            raise ValueError(
                f"model_set must be 'standard' or 'comprehensive', "
                f"got {model_set}"
            )
        
        models = {}
        for name, config in configs.items():
            # Create new consensus generator with specific config
            generator = PharmacophoreConsensus(
                tolerance=config['tolerance'],
                occurrence_threshold=config['occurrence_threshold'],
                linkage=self.linkage
            )
            
            # Generate consensus
            models[name] = generator.generate_consensus(
                mols=mols,
                features=features
            )
        
        return models


if __name__ == "__main__":
    import doctest
    doctest.testmod()

"""
Hybrid 2D+3D Scoring Module

Combines 2D fingerprint similarity with 3D shape-based pharmacophore scoring
for improved virtual screening performance.

Literature:
- Sanders et al. (2012): Combined 2D+3D outperforms either method alone
- Moshawih et al. (2024): Consensus holistic virtual screening approach
- Multiple DUD-E benchmarks confirm 2D+3D superiority

Expected improvement: +10-15% ROC-AUC vs 3D alone
"""

from typing import List, Dict, Optional, Tuple
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol
import numpy as np
import warnings


class HybridScorer:
    """
    Combines 2D fingerprint similarity with 3D shape/pharmacophore scoring.
    
    The hybrid approach leverages orthogonal information:
    - 2D fingerprints: Capture chemical similarity (functional groups, topology)
    - 3D shape: Capture spatial fit (binding geometry, pharmacophore features)
    
    Attributes:
        reference_mols: Original reference ligands for 2D similarity
        consensus_mol: 3D consensus pharmacophore for shape alignment
        alpha: Weight for 2D component (0-1), default 0.6
        fingerprint_type: Type of fingerprint ('morgan', 'rdkit', 'topological')
        fingerprint_params: Parameters for fingerprint generation
    """
    
    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        consensus_mol: Optional[Chem.Mol] = None,
        alpha: float = 0.6,
        fingerprint_type: str = 'morgan',
        radius: int = 2,
        n_bits: int = 1024,
        use_colors: bool = True,
        max_align_iters: int = 50,
        reference_mols_3d: Optional[List[Chem.Mol]] = None
    ):
        """
        Initialize hybrid scorer.

        Args:
            reference_mols: List of reference ligands (for 2D similarity)
            consensus_mol: Consensus pharmacophore mol (deprecated, use reference_mols_3d)
            alpha: Weight for 2D component (0.0-1.0)
                   0.0 = pure 3D, 0.6 = balanced, 1.0 = pure 2D
            fingerprint_type: 'morgan' (ECFP-like), 'rdkit', or 'topological'
            radius: Morgan fingerprint radius (2 = ECFP4, 3 = ECFP6)
            n_bits: Fingerprint size (1024 or 2048 recommended)
            use_colors: Include pharmacophore colors in 3D scoring
            max_align_iters: Maximum alignment iterations for shape matching
            reference_mols_3d: Reference molecules for 3D alignment (recommended
                over consensus_mol). If None, uses reference_mols.

        Raises:
            ValueError: If alpha not in [0, 1] or invalid fingerprint type
        """
        if not 0.0 <= alpha <= 1.0:
            raise ValueError(f"alpha must be in [0, 1], got {alpha}")

        if fingerprint_type not in ['morgan', 'rdkit', 'topological']:
            raise ValueError(f"fingerprint_type must be 'morgan', 'rdkit', or 'topological', got {fingerprint_type}")

        if consensus_mol is not None:
            warnings.warn(
                "consensus_mol parameter is deprecated and produces "
                "anti-discriminative 3D scores. Pass reference_mols_3d "
                "instead for reference-based 3D alignment.",
                DeprecationWarning,
                stacklevel=2
            )

        self.reference_mols = reference_mols
        self.consensus_mol = consensus_mol
        self.alpha = alpha
        self.fingerprint_type = fingerprint_type
        self.radius = radius
        self.n_bits = n_bits
        self.use_colors = use_colors
        self.max_align_iters = max_align_iters

        # Prepare reference molecules for 3D alignment
        self._prepared_refs_3d = []
        refs_3d = reference_mols_3d if reference_mols_3d is not None else reference_mols
        for mol in refs_3d:
            if mol is None:
                continue
            mol_h = Chem.AddHs(mol)
            if mol_h.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol_h, randomSeed=42)
            if mol_h.GetNumConformers() > 0:
                self._prepared_refs_3d.append(mol_h)

        # Pre-compute reference fingerprints for speed
        self.reference_fps = [
            self._compute_fingerprint(mol)
            for mol in reference_mols
        ]
        
    def _compute_fingerprint(self, mol: Chem.Mol):
        """
        Compute molecular fingerprint.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            ExplicitBitVect fingerprint
        """
        if self.fingerprint_type == 'morgan':
            # Morgan fingerprint (ECFP-like)
            return AllChem.GetMorganFingerprintAsBitVect(
                mol, 
                self.radius, 
                nBits=self.n_bits
            )
        elif self.fingerprint_type == 'rdkit':
            # RDKit topological fingerprint
            return Chem.RDKFingerprint(mol, fpSize=self.n_bits)
        elif self.fingerprint_type == 'topological':
            # Topological fingerprint
            return AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(
                mol,
                nBits=self.n_bits
            )
        else:
            raise ValueError(f"Unknown fingerprint type: {self.fingerprint_type}")
    
    def _compute_2d_similarity(self, query_mol: Chem.Mol) -> float:
        """
        Compute 2D Tanimoto similarity to reference set.
        
        Uses maximum similarity to any reference (nearest neighbor).
        
        Args:
            query_mol: Query molecule
            
        Returns:
            Max Tanimoto similarity in [0, 1]
        """
        query_fp = self._compute_fingerprint(query_mol)
        
        # Max similarity to any reference (nearest neighbor)
        similarities = [
            DataStructs.TanimotoSimilarity(query_fp, ref_fp)
            for ref_fp in self.reference_fps
        ]
        
        return max(similarities)
    
    def _compute_3d_similarity(self, query_mol: Chem.Mol) -> Dict[str, float]:
        """
        Compute 3D shape/pharmacophore similarity.

        When reference molecules are available (recommended), aligns query
        to each reference and takes the best score. Falls back to consensus_mol
        alignment if no references are prepared.

        Args:
            query_mol: Query molecule with 3D conformers

        Returns:
            Dict with 'shape', 'color', and '3d' (normalized combo) scores
        """
        if query_mol.GetNumConformers() == 0:
            return {'shape': 0.0, 'color': 0.0, '3d': 0.0}

        mol_copy = Chem.Mol(query_mol)

        best_shape = 0.0
        best_color = 0.0
        best_combo = 0.0

        if self._prepared_refs_3d:
            # Reference-based 3D scoring (recommended)
            for ref in self._prepared_refs_3d:
                for conf_id in range(mol_copy.GetNumConformers()):
                    try:
                        shape, color = AlignMol(
                            ref=ref,
                            probe=mol_copy,
                            probeConfId=conf_id,
                            useColors=self.use_colors,
                            max_preiters=self.max_align_iters
                        )
                        combo = shape + color
                        if combo > best_combo:
                            best_combo = combo
                            best_shape = shape
                            best_color = color
                    except Exception:
                        continue
        elif self.consensus_mol is not None:
            # Legacy consensus_mol path (deprecated)
            for conf_id in range(mol_copy.GetNumConformers()):
                try:
                    shape, color = AlignMol(
                        ref=self.consensus_mol,
                        probe=mol_copy,
                        refConfId=-1,
                        probeConfId=conf_id,
                        useColors=self.use_colors,
                        max_preiters=self.max_align_iters
                    )
                    combo = shape + color
                    if combo > best_combo:
                        best_combo = combo
                        best_shape = shape
                        best_color = color
                except Exception:
                    continue
        else:
            return {'shape': 0.0, 'color': 0.0, '3d': 0.0}

        # Normalize combo to [0, 1] for fair weighting with 2D
        if self.use_colors:
            normalized_3d = (best_shape + best_color) / 2.0
        else:
            normalized_3d = best_shape

        return {
            'shape': best_shape,
            'color': best_color,
            '3d': normalized_3d
        }
    
    def score(self, query_mol: Chem.Mol) -> Dict[str, float]:
        """
        Compute hybrid 2D+3D score for query molecule.
        
        Args:
            query_mol: Query molecule (must have 3D conformers)
            
        Returns:
            Dict containing:
                - 'hybrid': Weighted 2D+3D score (primary metric)
                - '2d': 2D Tanimoto similarity
                - '3d': Normalized 3D combo score
                - 'shape': Shape Tanimoto
                - 'color': Color Tanimoto
                
        Example:
            >>> scorer = HybridScorer(refs, consensus_mol, alpha=0.6)
            >>> result = scorer.score(query_mol)
            >>> print(f"Hybrid score: {result['hybrid']:.3f}")
            Hybrid score: 0.823
        """
        # Compute components
        similarity_2d = self._compute_2d_similarity(query_mol)
        similarity_3d = self._compute_3d_similarity(query_mol)
        
        # Weighted combination
        hybrid_score = (
            self.alpha * similarity_2d + 
            (1 - self.alpha) * similarity_3d['3d']
        )
        
        return {
            'hybrid': hybrid_score,
            '2d': similarity_2d,
            '3d': similarity_3d['3d'],
            'shape': similarity_3d['shape'],
            'color': similarity_3d['color']
        }
    
    def score_batch(
        self, 
        query_mols: List[Chem.Mol],
        verbose: bool = False
    ) -> List[Dict[str, float]]:
        """
        Score multiple molecules efficiently.
        
        Args:
            query_mols: List of query molecules
            verbose: Print progress
            
        Returns:
            List of score dictionaries
        """
        results = []
        for i, mol in enumerate(query_mols):
            if verbose and (i + 1) % 100 == 0:
                print(f"Scored {i + 1}/{len(query_mols)} molecules...")
            
            try:
                result = self.score(mol)
                results.append(result)
            except Exception as e:
                # Handle failures gracefully
                if verbose:
                    print(f"Warning: Failed to score molecule {i}: {e}")
                results.append({
                    'hybrid': 0.0,
                    '2d': 0.0,
                    '3d': 0.0,
                    'shape': 0.0,
                    'color': 0.0,
                    'error': str(e)
                })
        
        return results
    
    def optimize_alpha(
        self,
        test_mols: List[Chem.Mol],
        labels: List[int],
        alpha_range: Tuple[float, float] = (0.4, 0.8),
        n_points: int = 9
    ) -> Tuple[float, float]:
        """
        Find optimal alpha by grid search on test set.
        
        Args:
            test_mols: Test molecules with 3D conformers
            labels: Binary labels (1=active, 0=decoy)
            alpha_range: (min, max) alpha to test
            n_points: Number of alpha values to test
            
        Returns:
            (best_alpha, best_roc_auc)
            
        Example:
            >>> best_alpha, best_auc = scorer.optimize_alpha(test_mols, labels)
            >>> print(f"Optimal alpha: {best_alpha:.2f} (AUC: {best_auc:.3f})")
        """
        from sklearn.metrics import roc_auc_score
        
        alphas = np.linspace(alpha_range[0], alpha_range[1], n_points)
        best_alpha = self.alpha
        best_auc = 0.0
        
        original_alpha = self.alpha
        
        for alpha in alphas:
            self.alpha = alpha
            
            # Score all molecules
            scores = [
                self.score(mol)['hybrid'] 
                for mol in test_mols
            ]
            
            # Calculate AUC
            auc = roc_auc_score(labels, scores)
            
            if auc > best_auc:
                best_auc = auc
                best_alpha = alpha
        
        # Restore original alpha
        self.alpha = original_alpha
        
        return best_alpha, best_auc
    
    def get_info(self) -> Dict:
        """
        Get scorer configuration information.
        
        Returns:
            Dict with configuration details
        """
        return {
            'n_references': len(self.reference_mols),
            'alpha': self.alpha,
            'fingerprint_type': self.fingerprint_type,
            'fingerprint_radius': self.radius,
            'fingerprint_bits': self.n_bits,
            'use_colors': self.use_colors,
            'max_align_iters': self.max_align_iters,
            '2d_weight': self.alpha,
            '3d_weight': 1 - self.alpha
        }


def compare_scoring_methods(
    reference_mols: List[Chem.Mol],
    consensus_mol: Chem.Mol,
    test_mols: List[Chem.Mol],
    labels: List[int],
    alphas: List[float] = [0.0, 0.3, 0.5, 0.6, 0.7, 1.0]
) -> Dict[str, Dict[str, float]]:
    """
    Compare different alpha values (2D vs 3D weighting).
    
    Args:
        reference_mols: Reference ligands
        consensus_mol: Consensus pharmacophore
        test_mols: Test molecules
        labels: Binary labels (1=active, 0=decoy)
        alphas: List of alpha values to test
        
    Returns:
        Dict mapping alpha to metrics (roc_auc, bedroc, ef_1pct, etc.)
        
    Example:
        >>> results = compare_scoring_methods(refs, consensus, test_mols, labels)
        >>> for alpha, metrics in results.items():
        >>>     print(f"Alpha {alpha}: AUC = {metrics['roc_auc']:.3f}")
    """
    from sklearn.metrics import roc_auc_score
    
    results = {}
    
    for alpha in alphas:
        scorer = HybridScorer(reference_mols, consensus_mol=consensus_mol, alpha=alpha)
        
        # Score all molecules
        scores = [scorer.score(mol)['hybrid'] for mol in test_mols]
        
        # Calculate metrics
        auc = roc_auc_score(labels, scores)
        
        # Store results
        if alpha == 0.0:
            method_name = "Pure 3D"
        elif alpha == 1.0:
            method_name = "Pure 2D"
        else:
            method_name = f"Hybrid (α={alpha})"
        
        results[method_name] = {
            'alpha': alpha,
            'roc_auc': auc,
            'scores': scores
        }
    
    return results

"""Multi-level caching for pharmacophore optimization.

This module provides caching infrastructure to accelerate Bayesian optimization
by avoiding redundant consensus generation and pharmacophore mol conversion.

Cache Levels:
- L1: Conformer cache (handled by ConsensusOptimizer, SMILES-keyed)
- L2: Consensus feature cache (parameter-keyed)
- L3: Pharmacophore Mol cache (feature-hash-keyed)

Example:
    >>> from pharmacophore.caching import PharmacophoreCache, ConsensusCacheKey
    >>> cache = PharmacophoreCache(max_size=100)
    >>> key = ConsensusCacheKey(tolerance=2.0, occurrence=0.5, linkage='average', ref_hash='abc123')
    >>> cache.set_consensus(key, features)
    >>> cached = cache.get_consensus(key)
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Any, Tuple
import hashlib
import json
from collections import OrderedDict

from rdkit import Chem


@dataclass(frozen=True)
class ConsensusCacheKey:
    """Immutable, hashable key for consensus cache lookups.

    Attributes:
        tolerance: Clustering distance threshold (rounded to 8 decimals).
        occurrence: Minimum feature occurrence fraction (rounded to 8 decimals).
        linkage: Clustering linkage method ('average', 'complete', etc.).
        ref_hash: MD5 hash of sorted reference SMILES for identity.

    Example:
        >>> key1 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        >>> key2 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        >>> key1 == key2
        True
        >>> hash(key1) == hash(key2)
        True
    """
    tolerance: float
    occurrence: float
    linkage: str
    ref_hash: str

    def __post_init__(self):
        # Ensure tolerance and occurrence are rounded for floating-point stability
        # Using object.__setattr__ because frozen=True
        object.__setattr__(self, 'tolerance', round(self.tolerance, 8))
        object.__setattr__(self, 'occurrence', round(self.occurrence, 8))


class LRUCache:
    """Simple LRU cache with max size limit.

    Used internally by PharmacophoreCache to prevent unbounded memory growth.
    """

    def __init__(self, max_size: int = 100):
        """Initialize LRU cache.

        Args:
            max_size: Maximum number of entries. Oldest entries evicted when exceeded.
        """
        self.max_size = max_size
        self._cache: OrderedDict = OrderedDict()
        self._hits = 0
        self._misses = 0

    def get(self, key) -> Optional[Any]:
        """Get value from cache, moving key to end (most recent)."""
        if key in self._cache:
            self._cache.move_to_end(key)
            self._hits += 1
            return self._cache[key]
        self._misses += 1
        return None

    def set(self, key, value) -> None:
        """Set value in cache, evicting oldest if at max size."""
        if key in self._cache:
            self._cache.move_to_end(key)
        else:
            if len(self._cache) >= self.max_size:
                self._cache.popitem(last=False)  # Remove oldest
        self._cache[key] = value

    def clear(self) -> None:
        """Clear all cache entries and reset statistics."""
        self._cache.clear()
        self._hits = 0
        self._misses = 0

    def __len__(self) -> int:
        return len(self._cache)

    def __contains__(self, key) -> bool:
        return key in self._cache

    @property
    def hit_rate(self) -> float:
        """Calculate cache hit rate."""
        total = self._hits + self._misses
        return self._hits / total if total > 0 else 0.0


class PharmacophoreCache:
    """Multi-level cache manager for pharmacophore optimization.

    Provides caching at multiple levels to accelerate repeated evaluations:
    - L2: Consensus features (keyed by parameters + reference hash)
    - L3: Pharmacophore Mol objects (keyed by feature hash)

    Attributes:
        max_size: Maximum entries per cache level.

    Example:
        >>> cache = PharmacophoreCache(max_size=100)
        >>> key = ConsensusCacheKey(2.0, 0.5, 'average', ref_hash)
        >>> cache.set_consensus(key, features)
        >>> print(cache.stats())
        {'consensus_hits': 1, 'consensus_misses': 0, ...}
    """

    def __init__(self, max_size: int = 100):
        """Initialize cache manager.

        Args:
            max_size: Maximum entries per cache level (default: 100).
        """
        self.max_size = max_size
        self._consensus_cache = LRUCache(max_size)
        self._pharm_mol_cache = LRUCache(max_size)

    def get_consensus(self, key: ConsensusCacheKey) -> Optional[List]:
        """Get consensus features from cache.

        Args:
            key: Cache key with parameters and reference hash.

        Returns:
            List of consensus features if cached, None otherwise.
        """
        return self._consensus_cache.get(key)

    def set_consensus(self, key: ConsensusCacheKey, features: List) -> None:
        """Store consensus features in cache.

        Args:
            key: Cache key with parameters and reference hash.
            features: List of consensus features to cache.
        """
        self._consensus_cache.set(key, features)

    def get_pharm_mol(self, feature_hash: str) -> Optional[Chem.Mol]:
        """Get pharmacophore Mol from cache.

        Args:
            feature_hash: Hash of consensus features.

        Returns:
            RDKit Mol if cached, None otherwise.
        """
        return self._pharm_mol_cache.get(feature_hash)

    def set_pharm_mol(self, feature_hash: str, mol: Chem.Mol) -> None:
        """Store pharmacophore Mol in cache.

        Args:
            feature_hash: Hash of consensus features.
            mol: RDKit Mol to cache.
        """
        self._pharm_mol_cache.set(feature_hash, mol)

    def clear(self) -> None:
        """Clear all caches."""
        self._consensus_cache.clear()
        self._pharm_mol_cache.clear()

    def stats(self) -> Dict[str, Any]:
        """Get cache statistics.

        Returns:
            Dict with hits, misses, and hit rates for each cache level.
        """
        return {
            'consensus_hits': self._consensus_cache._hits,
            'consensus_misses': self._consensus_cache._misses,
            'consensus_hit_rate': self._consensus_cache.hit_rate,
            'consensus_size': len(self._consensus_cache),
            'pharm_mol_hits': self._pharm_mol_cache._hits,
            'pharm_mol_misses': self._pharm_mol_cache._misses,
            'pharm_mol_hit_rate': self._pharm_mol_cache.hit_rate,
            'pharm_mol_size': len(self._pharm_mol_cache),
        }

    def __repr__(self) -> str:
        stats = self.stats()
        return (
            f"PharmacophoreCache("
            f"consensus={stats['consensus_size']}/{self.max_size}, "
            f"pharm_mol={stats['pharm_mol_size']}/{self.max_size}, "
            f"hit_rate={stats['consensus_hit_rate']:.1%})"
        )


def hash_features(features: List) -> str:
    """Create deterministic hash of consensus features.

    Used for caching pharmacophore Mol objects that depend on feature content.

    Args:
        features: List of features in format [type, indices, x, y, z].

    Returns:
        MD5 hash string of feature content.

    Example:
        >>> features = [['Donor', (), 1.0, 2.0, 3.0], ['Acceptor', (), 4.0, 5.0, 6.0]]
        >>> hash1 = hash_features(features)
        >>> hash2 = hash_features(features)
        >>> hash1 == hash2
        True
    """
    # Convert features to canonical string representation
    # Sort by type then by coordinates for determinism
    canonical = []
    for feat in features:
        feat_type = feat[0]
        coords = (round(feat[2], 4), round(feat[3], 4), round(feat[4], 4))
        canonical.append((feat_type, coords))

    # Sort for determinism (features may come in different orders)
    canonical.sort(key=lambda x: (x[0], x[1]))

    # Create hash
    content = json.dumps(canonical, sort_keys=True)
    return hashlib.md5(content.encode()).hexdigest()


def hash_reference_mols(mols: List[Chem.Mol]) -> str:
    """Create deterministic hash of reference molecules.

    Used as part of cache key to invalidate cache when references change.

    Args:
        mols: List of RDKit Mol objects.

    Returns:
        MD5 hash string of canonical SMILES.

    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles('CCO')
        >>> hash1 = hash_reference_mols([mol])
        >>> hash2 = hash_reference_mols([mol])
        >>> hash1 == hash2
        True
    """
    smiles_list = []
    for mol in mols:
        if mol is not None:
            try:
                smiles = Chem.MolToSmiles(mol, canonical=True)
                smiles_list.append(smiles)
            except Exception:
                smiles_list.append('')

    # Sort for determinism
    smiles_list.sort()
    content = '|'.join(smiles_list)
    return hashlib.md5(content.encode()).hexdigest()


def create_cache_key(
    tolerance: float,
    occurrence: float,
    linkage: str,
    reference_mols: List[Chem.Mol]
) -> ConsensusCacheKey:
    """Convenience function to create cache key from parameters and molecules.

    Args:
        tolerance: Clustering distance threshold.
        occurrence: Minimum feature occurrence fraction.
        linkage: Clustering linkage method.
        reference_mols: Reference molecules for consensus.

    Returns:
        ConsensusCacheKey for cache lookups.

    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles('CCO')
        >>> key = create_cache_key(2.0, 0.5, 'average', [mol])
        >>> isinstance(key, ConsensusCacheKey)
        True
    """
    ref_hash = hash_reference_mols(reference_mols)
    return ConsensusCacheKey(
        tolerance=tolerance,
        occurrence=occurrence,
        linkage=linkage,
        ref_hash=ref_hash
    )

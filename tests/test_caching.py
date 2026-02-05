"""Tests for pharmacophore caching module."""

import pytest
from rdkit import Chem

from pharmacophore.caching import (
    ConsensusCacheKey,
    LRUCache,
    PharmacophoreCache,
    hash_features,
    hash_reference_mols,
    create_cache_key,
)


class TestConsensusCacheKey:
    """Tests for ConsensusCacheKey."""

    def test_equality_same_params(self):
        """Keys with same parameters should be equal."""
        key1 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        key2 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        assert key1 == key2

    def test_hash_same_params(self):
        """Keys with same parameters should have same hash."""
        key1 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        key2 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        assert hash(key1) == hash(key2)

    def test_inequality_different_tolerance(self):
        """Keys with different tolerance should not be equal."""
        key1 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        key2 = ConsensusCacheKey(2.5, 0.5, 'average', 'abc123')
        assert key1 != key2

    def test_inequality_different_occurrence(self):
        """Keys with different occurrence should not be equal."""
        key1 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        key2 = ConsensusCacheKey(2.0, 0.6, 'average', 'abc123')
        assert key1 != key2

    def test_inequality_different_linkage(self):
        """Keys with different linkage should not be equal."""
        key1 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        key2 = ConsensusCacheKey(2.0, 0.5, 'complete', 'abc123')
        assert key1 != key2

    def test_inequality_different_ref_hash(self):
        """Keys with different ref_hash should not be equal."""
        key1 = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        key2 = ConsensusCacheKey(2.0, 0.5, 'average', 'xyz789')
        assert key1 != key2

    def test_floating_point_rounding(self):
        """Values differing only beyond 8 decimal places should be equal."""
        key1 = ConsensusCacheKey(2.000000001, 0.500000001, 'average', 'abc123')
        key2 = ConsensusCacheKey(2.000000002, 0.500000002, 'average', 'abc123')
        assert key1 == key2

    def test_no_false_collision_at_5th_decimal(self):
        """Values differing at 5th decimal should NOT collide."""
        key1 = ConsensusCacheKey(2.00001, 0.50001, 'average', 'abc123')
        key2 = ConsensusCacheKey(2.00004, 0.50004, 'average', 'abc123')
        assert key1 != key2

    def test_usable_as_dict_key(self):
        """Key should be usable as dictionary key."""
        key = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        cache = {key: 'value'}
        assert cache[key] == 'value'


class TestLRUCache:
    """Tests for LRU cache implementation."""

    def test_basic_set_get(self):
        """Basic set and get operations."""
        cache = LRUCache(max_size=10)
        cache.set('key1', 'value1')
        assert cache.get('key1') == 'value1'

    def test_miss_returns_none(self):
        """Missing key should return None."""
        cache = LRUCache(max_size=10)
        assert cache.get('nonexistent') is None

    def test_eviction_on_max_size(self):
        """Oldest entries should be evicted when max size exceeded."""
        cache = LRUCache(max_size=3)
        cache.set('a', 1)
        cache.set('b', 2)
        cache.set('c', 3)
        cache.set('d', 4)  # Should evict 'a'

        assert cache.get('a') is None
        assert cache.get('b') == 2
        assert cache.get('c') == 3
        assert cache.get('d') == 4

    def test_access_updates_recency(self):
        """Accessing an entry should move it to end (most recent)."""
        cache = LRUCache(max_size=3)
        cache.set('a', 1)
        cache.set('b', 2)
        cache.set('c', 3)

        # Access 'a' to make it recent
        _ = cache.get('a')

        cache.set('d', 4)  # Should evict 'b' (now oldest)

        assert cache.get('a') == 1
        assert cache.get('b') is None
        assert cache.get('c') == 3
        assert cache.get('d') == 4

    def test_hit_miss_tracking(self):
        """Hits and misses should be tracked."""
        cache = LRUCache(max_size=10)
        cache.set('key', 'value')

        _ = cache.get('key')      # Hit
        _ = cache.get('key')      # Hit
        _ = cache.get('missing')  # Miss

        assert cache._hits == 2
        assert cache._misses == 1

    def test_hit_rate(self):
        """Hit rate calculation."""
        cache = LRUCache(max_size=10)
        cache.set('key', 'value')

        _ = cache.get('key')      # Hit
        _ = cache.get('missing')  # Miss

        assert cache.hit_rate == 0.5

    def test_clear(self):
        """Clear should empty cache and reset stats."""
        cache = LRUCache(max_size=10)
        cache.set('key', 'value')
        _ = cache.get('key')

        cache.clear()

        assert len(cache) == 0
        assert cache._hits == 0
        assert cache._misses == 0


class TestPharmacophoreCache:
    """Tests for multi-level pharmacophore cache."""

    def test_consensus_cache(self):
        """Consensus cache basic operations."""
        cache = PharmacophoreCache(max_size=10)
        key = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')
        features = [['Donor', (), 1.0, 2.0, 3.0]]

        cache.set_consensus(key, features)
        result = cache.get_consensus(key)

        assert result == features

    def test_pharm_mol_cache(self):
        """Pharmacophore mol cache basic operations."""
        cache = PharmacophoreCache(max_size=10)
        mol = Chem.MolFromSmiles('CCO')
        feature_hash = 'test_hash'

        cache.set_pharm_mol(feature_hash, mol)
        result = cache.get_pharm_mol(feature_hash)

        assert result is not None
        assert Chem.MolToSmiles(result) == Chem.MolToSmiles(mol)

    def test_stats(self):
        """Cache statistics tracking."""
        cache = PharmacophoreCache(max_size=10)
        key = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')

        # Miss
        _ = cache.get_consensus(key)
        # Set and hit
        cache.set_consensus(key, [])
        _ = cache.get_consensus(key)

        stats = cache.stats()
        assert stats['consensus_hits'] == 1
        assert stats['consensus_misses'] == 1
        assert stats['consensus_size'] == 1

    def test_clear(self):
        """Clear should reset both caches."""
        cache = PharmacophoreCache(max_size=10)
        key = ConsensusCacheKey(2.0, 0.5, 'average', 'abc123')

        cache.set_consensus(key, [])
        cache.set_pharm_mol('hash', Chem.MolFromSmiles('C'))

        cache.clear()

        stats = cache.stats()
        assert stats['consensus_size'] == 0
        assert stats['pharm_mol_size'] == 0


class TestHashFunctions:
    """Tests for hash utility functions."""

    def test_hash_features_deterministic(self):
        """Same features should produce same hash."""
        features = [
            ['Donor', (), 1.0, 2.0, 3.0],
            ['Acceptor', (), 4.0, 5.0, 6.0]
        ]
        hash1 = hash_features(features)
        hash2 = hash_features(features)
        assert hash1 == hash2

    def test_hash_features_order_independent(self):
        """Feature order should not affect hash."""
        features1 = [
            ['Donor', (), 1.0, 2.0, 3.0],
            ['Acceptor', (), 4.0, 5.0, 6.0]
        ]
        features2 = [
            ['Acceptor', (), 4.0, 5.0, 6.0],
            ['Donor', (), 1.0, 2.0, 3.0]
        ]
        hash1 = hash_features(features1)
        hash2 = hash_features(features2)
        assert hash1 == hash2

    def test_hash_features_different_coords(self):
        """Different coordinates should produce different hash."""
        features1 = [['Donor', (), 1.0, 2.0, 3.0]]
        features2 = [['Donor', (), 1.0, 2.0, 4.0]]
        hash1 = hash_features(features1)
        hash2 = hash_features(features2)
        assert hash1 != hash2

    def test_hash_reference_mols_deterministic(self):
        """Same molecules should produce same hash."""
        mol = Chem.MolFromSmiles('CCO')
        hash1 = hash_reference_mols([mol])
        hash2 = hash_reference_mols([mol])
        assert hash1 == hash2

    def test_hash_reference_mols_order_independent(self):
        """Molecule order should not affect hash."""
        mol1 = Chem.MolFromSmiles('CCO')
        mol2 = Chem.MolFromSmiles('CCC')
        hash1 = hash_reference_mols([mol1, mol2])
        hash2 = hash_reference_mols([mol2, mol1])
        assert hash1 == hash2

    def test_hash_reference_mols_different_mols(self):
        """Different molecules should produce different hash."""
        mol1 = Chem.MolFromSmiles('CCO')
        mol2 = Chem.MolFromSmiles('CCCO')
        hash1 = hash_reference_mols([mol1])
        hash2 = hash_reference_mols([mol2])
        assert hash1 != hash2


class TestCreateCacheKey:
    """Tests for cache key creation utility."""

    def test_creates_valid_key(self):
        """Should create valid ConsensusCacheKey."""
        mol = Chem.MolFromSmiles('CCO')
        key = create_cache_key(2.0, 0.5, 'average', [mol])

        assert isinstance(key, ConsensusCacheKey)
        assert key.tolerance == 2.0
        assert key.occurrence == 0.5
        assert key.linkage == 'average'
        assert len(key.ref_hash) == 32  # MD5 hash length

    def test_same_mols_same_key(self):
        """Same molecules should produce same key."""
        mol = Chem.MolFromSmiles('CCO')
        key1 = create_cache_key(2.0, 0.5, 'average', [mol])
        key2 = create_cache_key(2.0, 0.5, 'average', [mol])
        assert key1 == key2

    def test_different_mols_different_key(self):
        """Different molecules should produce different key."""
        mol1 = Chem.MolFromSmiles('CCO')
        mol2 = Chem.MolFromSmiles('CCCO')
        key1 = create_cache_key(2.0, 0.5, 'average', [mol1])
        key2 = create_cache_key(2.0, 0.5, 'average', [mol2])
        assert key1 != key2

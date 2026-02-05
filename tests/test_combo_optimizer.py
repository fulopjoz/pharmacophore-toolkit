"""Tests for ComboPharmacophoreOptimizer.

Uses simple drug-like molecules as fixtures for fast conformer generation.
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.combo_optimizer import ComboPharmacophoreOptimizer


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Drug-like SMILES that have donors, acceptors, aromatics, hydrophobes
REF_SMILES = [
    "c1ccc(O)cc1",       # phenol
    "c1ccc(N)cc1",       # aniline
    "c1ccc(C(=O)O)cc1",  # benzoic acid
    "c1ccc(OC)cc1",      # anisole
    "c1ccc(CC)cc1",      # ethylbenzene
]

ACTIVE_SMILES = [
    "c1ccc(O)c(N)c1",        # 2-aminophenol
    "c1ccc(C(=O)O)c(O)c1",   # salicylic acid
    "c1ccc(OC)c(N)c1",       # 2-methoxyanline
]

DECOY_SMILES = [
    "CCCCCCCC",           # octane (no aromatics)
    "CCCCC(=O)O",         # pentanoic acid
    "CC(C)CC(C)C",        # 2,4-dimethylpentane
    "CCCCCCCCCC",         # decane
    "CCC(CC)CC",          # 3-ethylpentane
]


@pytest.fixture
def ref_smiles():
    return list(REF_SMILES)


@pytest.fixture
def active_smiles():
    return list(ACTIVE_SMILES)


@pytest.fixture
def decoy_smiles():
    return list(DECOY_SMILES)


@pytest.fixture
def optimizer():
    return ComboPharmacophoreOptimizer(random_state=42, verbose=False)


@pytest.fixture
def loaded_optimizer(optimizer, ref_smiles, active_smiles, decoy_smiles):
    optimizer.load_from_smiles(ref_smiles, active_smiles, decoy_smiles)
    return optimizer


@pytest.fixture
def temp_data_dir(ref_smiles, active_smiles, decoy_smiles):
    """Create temporary SDF + CSV files for file-loading tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create SDF from SMILES
        sdf_path = tmpdir / "refs.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        for smi in ref_smiles:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            writer.write(mol)
        writer.close()

        # Create actives CSV
        actives_path = tmpdir / "actives.csv"
        with open(actives_path, "w") as f:
            f.write("SMILES,activity\n")
            for smi in active_smiles:
                f.write(f"{smi},1\n")

        # Create decoys CSV
        decoys_path = tmpdir / "decoys.csv"
        with open(decoys_path, "w") as f:
            f.write("Smiles,label\n")
            for smi in decoy_smiles:
                f.write(f"{smi},0\n")

        yield {
            "sdf": str(sdf_path),
            "actives": str(actives_path),
            "decoys": str(decoys_path),
            "dir": str(tmpdir),
        }


# ---------------------------------------------------------------------------
# Tests: Initialization
# ---------------------------------------------------------------------------

class TestInit:
    def test_default_init(self):
        opt = ComboPharmacophoreOptimizer()
        assert opt.random_state == 42
        assert opt.verbose is True
        assert opt._ref_smiles is None
        assert opt.best_result is None

    def test_custom_params(self):
        opt = ComboPharmacophoreOptimizer(random_state=123, verbose=False)
        assert opt.random_state == 123
        assert opt.verbose is False


# ---------------------------------------------------------------------------
# Tests: Data Loading
# ---------------------------------------------------------------------------

class TestDataLoading:
    def test_load_from_smiles(self, optimizer, ref_smiles, active_smiles, decoy_smiles):
        result = optimizer.load_from_smiles(ref_smiles, active_smiles, decoy_smiles)
        assert result is optimizer  # method chaining
        assert len(optimizer._ref_smiles) == 5
        assert len(optimizer._active_smiles) == 3
        assert len(optimizer._decoy_smiles) == 5

    def test_load_from_smiles_too_few_refs(self, optimizer, active_smiles, decoy_smiles):
        with pytest.raises(ValueError, match="at least 2"):
            optimizer.load_from_smiles(["CCO"], active_smiles, decoy_smiles)

    def test_load_from_smiles_empty_actives(self, optimizer, ref_smiles, decoy_smiles):
        with pytest.raises(ValueError, match="empty"):
            optimizer.load_from_smiles(ref_smiles, [], decoy_smiles)

    def test_load_from_smiles_empty_decoys(self, optimizer, ref_smiles, active_smiles):
        with pytest.raises(ValueError, match="empty"):
            optimizer.load_from_smiles(ref_smiles, active_smiles, [])

    def test_load_from_files(self, optimizer, temp_data_dir):
        result = optimizer.load_from_files(
            temp_data_dir["sdf"],
            temp_data_dir["actives"],
            temp_data_dir["decoys"],
        )
        assert result is optimizer
        assert len(optimizer._ref_smiles) == 5
        assert len(optimizer._active_smiles) == 3
        assert len(optimizer._decoy_smiles) == 5

    def test_load_clears_cache(self, loaded_optimizer, ref_smiles, active_smiles, decoy_smiles):
        # Populate some cache
        loaded_optimizer._consensus_cache.set("test", ([], None))
        # Reload
        loaded_optimizer.load_from_smiles(ref_smiles, active_smiles, decoy_smiles)
        assert loaded_optimizer._consensus_cache.get("test") is None


# ---------------------------------------------------------------------------
# Tests: Conformer Generation
# ---------------------------------------------------------------------------

class TestConformerGeneration:
    def test_ref_conformers_generated(self, loaded_optimizer):
        # After loading, ref conformers should be populated
        assert len(loaded_optimizer._ref_conformers) > 0

    def test_conformer_caching(self, loaded_optimizer):
        from pharmacophore.caching import LRUCache
        cache = LRUCache(max_size=50)
        mol1 = loaded_optimizer._generate_conformers("c1ccccc1", 3, cache)
        mol2 = loaded_optimizer._generate_conformers("c1ccccc1", 3, cache)
        assert mol1 is mol2  # same object from cache

    def test_different_n_conformers_not_cached(self, loaded_optimizer):
        from pharmacophore.caching import LRUCache
        cache = LRUCache(max_size=50)
        mol1 = loaded_optimizer._generate_conformers("c1ccccc1", 3, cache)
        mol2 = loaded_optimizer._generate_conformers("c1ccccc1", 5, cache)
        assert mol1 is not mol2  # different cache keys

    def test_invalid_smiles_returns_none(self, loaded_optimizer):
        from pharmacophore.caching import LRUCache
        cache = LRUCache(max_size=50)
        result = loaded_optimizer._generate_conformers("INVALID", 3, cache)
        assert result is None


# ---------------------------------------------------------------------------
# Tests: Reference Alignment
# ---------------------------------------------------------------------------

class TestReferenceAlignment:
    def test_align_returns_list(self, loaded_optimizer):
        aligned = loaded_optimizer._align_references_to_template(0)
        assert isinstance(aligned, list)
        assert len(aligned) >= 2  # at least 2 successfully aligned

    def test_all_have_conformers(self, loaded_optimizer):
        aligned = loaded_optimizer._align_references_to_template(0)
        for mol in aligned:
            assert mol.GetNumConformers() > 0

    def test_alignment_cached(self, loaded_optimizer):
        aligned1 = loaded_optimizer._align_references_to_template(0)
        aligned2 = loaded_optimizer._align_references_to_template(0)
        assert aligned1 is aligned2  # same list from cache

    def test_different_templates(self, loaded_optimizer):
        aligned0 = loaded_optimizer._align_references_to_template(0)
        aligned1 = loaded_optimizer._align_references_to_template(1)
        # Different templates produce different alignments (different objects)
        assert aligned0 is not aligned1


# ---------------------------------------------------------------------------
# Tests: Consensus Building
# ---------------------------------------------------------------------------

class TestConsensusBuilding:
    def test_consensus_produces_features(self, loaded_optimizer):
        features, pharm_mol = loaded_optimizer._build_consensus(0, 2.0, 0.3)
        assert isinstance(features, list)
        # With relaxed threshold (0.3), should get some features
        if features:
            assert len(features[0]) == 5  # [type, (), x, y, z]

    def test_consensus_cached(self, loaded_optimizer):
        r1 = loaded_optimizer._build_consensus(0, 2.0, 0.5)
        r2 = loaded_optimizer._build_consensus(0, 2.0, 0.5)
        assert r1 is r2  # same tuple from cache

    def test_strict_threshold_fewer_features(self, loaded_optimizer):
        features_relaxed, _ = loaded_optimizer._build_consensus(0, 2.0, 0.2)
        features_strict, _ = loaded_optimizer._build_consensus(0, 2.0, 0.8)
        # Strict threshold should give fewer or equal features
        assert len(features_strict) <= len(features_relaxed)


# ---------------------------------------------------------------------------
# Tests: Scoring
# ---------------------------------------------------------------------------

class TestScoring:
    def test_score_molecule_returns_float(self, loaded_optimizer):
        features, pharm_mol = loaded_optimizer._build_consensus(0, 2.0, 0.3)
        if pharm_mol is not None:
            score = loaded_optimizer._score_molecule(
                "c1ccc(O)cc1", pharm_mol, 0.5, 3
            )
            assert isinstance(score, float)
            assert 0.0 <= score <= 2.0

    def test_score_invalid_smiles_returns_zero(self, loaded_optimizer):
        features, pharm_mol = loaded_optimizer._build_consensus(0, 2.0, 0.3)
        if pharm_mol is not None:
            score = loaded_optimizer._score_molecule(
                "INVALID", pharm_mol, 0.5, 3
            )
            assert score == 0.0

    def test_score_all_returns_list(self, loaded_optimizer):
        features, pharm_mol = loaded_optimizer._build_consensus(0, 2.0, 0.3)
        if pharm_mol is not None:
            scores = loaded_optimizer._score_all(
                ACTIVE_SMILES, pharm_mol, 0.5, 3
            )
            assert len(scores) == len(ACTIVE_SMILES)
            assert all(isinstance(s, float) for s in scores)


# ---------------------------------------------------------------------------
# Tests: Evaluation
# ---------------------------------------------------------------------------

class TestEvaluation:
    def test_evaluate_returns_valid_auc(self, loaded_optimizer):
        auc = loaded_optimizer.evaluate(
            template_idx=0,
            tolerance=2.0,
            occurrence_threshold=0.3,
            opt_param=0.5,
            n_conformers=3,
        )
        assert isinstance(auc, float)
        assert 0.0 <= auc <= 1.0

    def test_evaluate_stores_history(self, loaded_optimizer):
        initial_len = len(loaded_optimizer.history)
        loaded_optimizer.evaluate(0, 2.0, 0.3, 0.5, 3)
        assert len(loaded_optimizer.history) == initial_len + 1
        assert "auc" in loaded_optimizer.history[-1]

    def test_evaluate_no_data_raises(self, optimizer):
        with pytest.raises(ValueError, match="No data loaded"):
            optimizer.evaluate(0, 2.0, 0.5, 0.5, 3)


# ---------------------------------------------------------------------------
# Tests: Optimization
# ---------------------------------------------------------------------------

class TestOptimization:
    def test_optimize_minimal_budget(self, loaded_optimizer):
        result = loaded_optimizer.optimize(
            n_calls=3,
            n_random_starts=3,
            early_stopping=False,
        )
        assert "best_params" in result
        assert "best_auc" in result
        assert "best_features" in result
        assert "n_evaluations" in result
        assert "elapsed_sec" in result
        assert result["n_evaluations"] == 3

    def test_optimize_stores_best_result(self, loaded_optimizer):
        loaded_optimizer.optimize(n_calls=3, n_random_starts=3, early_stopping=False)
        assert loaded_optimizer.best_result is not None
        assert loaded_optimizer.best_result["best_auc"] >= 0.0

    def test_optimize_best_params_keys(self, loaded_optimizer):
        result = loaded_optimizer.optimize(
            n_calls=3, n_random_starts=3, early_stopping=False
        )
        params = result["best_params"]
        assert "template_idx" in params
        assert "tolerance" in params
        assert "occurrence_threshold" in params
        assert "opt_param" in params
        assert "n_conformers" in params


# ---------------------------------------------------------------------------
# Tests: Export
# ---------------------------------------------------------------------------

class TestExport:
    def test_export_creates_files(self, loaded_optimizer):
        loaded_optimizer.optimize(n_calls=3, n_random_starts=3, early_stopping=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            outputs = loaded_optimizer.export_model(tmpdir)
            assert "json" in outputs
            assert Path(outputs["json"]).exists()

            # Verify JSON content
            with open(outputs["json"]) as f:
                data = json.load(f)
            assert "best_params" in data
            assert "best_auc" in data
            assert "features" in data
            assert "reference_smiles" in data

    def test_export_pml_if_features(self, loaded_optimizer):
        loaded_optimizer.optimize(n_calls=3, n_random_starts=3, early_stopping=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            outputs = loaded_optimizer.export_model(tmpdir)
            if "pml" in outputs:
                pml_content = Path(outputs["pml"]).read_text()
                assert "pseudoatom" in pml_content

    def test_export_no_result_raises(self, loaded_optimizer):
        with pytest.raises(ValueError, match="No optimization result"):
            loaded_optimizer.export_model("/tmp/test")


# ---------------------------------------------------------------------------
# Tests: Cache & Utilities
# ---------------------------------------------------------------------------

class TestUtilities:
    def test_clear_cache(self, loaded_optimizer):
        # Populate caches
        loaded_optimizer._align_references_to_template(0)
        assert len(loaded_optimizer._aligned_refs) > 0

        loaded_optimizer.clear_cache()
        assert len(loaded_optimizer._ref_conformers) == 0
        assert len(loaded_optimizer._query_conformers) == 0
        assert len(loaded_optimizer._aligned_refs) == 0
        assert len(loaded_optimizer._consensus_cache) == 0

    def test_cache_stats(self, loaded_optimizer):
        stats = loaded_optimizer.get_cache_stats()
        assert "L1_ref_conformers" in stats
        assert "L2_query_conformers" in stats
        assert "L3_aligned_refs" in stats
        assert "L4_consensus" in stats

    def test_repr_empty(self, optimizer):
        r = repr(optimizer)
        assert "refs=0" in r

    def test_repr_with_data(self, loaded_optimizer):
        r = repr(loaded_optimizer)
        assert "refs=5" in r
        assert "actives=3" in r


# ---------------------------------------------------------------------------
# Tests: Multi-Fidelity Optimization
# ---------------------------------------------------------------------------

class TestMultiFidelityOptimization:
    def test_mf_returns_result(self, loaded_optimizer):
        result = loaded_optimizer.optimize_multifidelity(
            n_calls=5,
            n_random_starts=5,
            n_conformers_final=2,
            refine_top_k=2,
            early_stopping=False,
        )
        assert "best_params" in result
        assert "best_auc" in result
        assert "best_features" in result
        assert "stage_info" in result

    def test_mf_stage_info(self, loaded_optimizer):
        result = loaded_optimizer.optimize_multifidelity(
            n_calls=5,
            n_random_starts=5,
            n_conformers_final=2,
            refine_top_k=2,
            early_stopping=False,
        )
        info = result["stage_info"]
        assert "explore_evals" in info
        assert "explore_time_sec" in info
        assert "refine_configs" in info
        assert info["explore_evals"] >= 3
        assert info["refine_configs"] <= 2

    def test_mf_best_params_has_conformers(self, loaded_optimizer):
        result = loaded_optimizer.optimize_multifidelity(
            n_calls=4,
            n_random_starts=4,
            n_conformers_final=7,
            refine_top_k=2,
            early_stopping=False,
        )
        assert result["best_params"]["n_conformers"] == 7

    def test_mf_no_data_raises(self, optimizer):
        with pytest.raises(ValueError, match="No data loaded"):
            optimizer.optimize_multifidelity()

    def test_mf_stores_best_result(self, loaded_optimizer):
        loaded_optimizer.optimize_multifidelity(
            n_calls=4,
            n_random_starts=4,
            n_conformers_final=2,
            refine_top_k=2,
            early_stopping=False,
        )
        assert loaded_optimizer.best_result is not None

    def test_mf_export_after_optimize(self, loaded_optimizer):
        loaded_optimizer.optimize_multifidelity(
            n_calls=4,
            n_random_starts=4,
            n_conformers_final=2,
            refine_top_k=2,
            early_stopping=False,
        )
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            outputs = loaded_optimizer.export_model(tmpdir, "mf_model")
            assert "json" in outputs

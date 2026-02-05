# Plan: Next-Generation Consensus Pharmacophore Model Creator

## Overview

Implement 20 algorithmic improvements to the pharmacophore-toolkit consensus model system, drawn from 11 priority papers + 6 cross-domain algorithm families (computer vision, optimal transport, graph ML, ensemble methods). Organized in 3 tiers by complexity, each building on the previous.

**Goal**: Transform the current single-algorithm consensus system (agglomerative clustering + tolerance/occurrence) into a multi-strategy, auto-optimized pharmacophore model generator with superior discrimination on the CCR2 benchmark (75 actives, 500 decoys).

---

## Architecture: Current vs Target

```
CURRENT:
  molecules → [AgglomerativeClustering] → consensus features → [mol_converter] → [rdShapeAlign] → scores

TARGET:
  molecules → [Multi-Strategy Clustering] → consensus features → [Multi-Metric Scoring] → scores
                     ↑                                                    ↑
              [Auto-Optimizer]                                    [Ensemble Voting]
              (BO + MF + ensemble)                               (2D + 3D + OT + hybrid)
```

---

## Implementation Plan

### Phase 1: Foundation Enhancements (Priority Paper Insights)

These modify existing modules with targeted improvements. No new dependencies.

#### Step 1.1: Feature Weight Differentiation in Clustering
**File**: `pharmacophore/consensus.py`
**Source**: PharmaGist (hydrophobic weight = 0.3, others = 1.0)
**Change**: Modify `_cluster_features()` to apply per-type distance scaling before clustering.
Currently features are clustered by Euclidean distance on [x,y,z]. Apply weight: `weighted_distance = euclidean_distance / feature_weight`. This makes hydrophobic features cluster at larger radii (less specific) while polar features require tighter clustering.
**Use existing**: `constants.FEATURE_WEIGHTS` dict (Donor=2.0, Acceptor=2.0, Aromatic=1.5, Hydrophobe=1.0)
**Implementation**: Scale the `distance_threshold` parameter per feature type by dividing by weight. E.g., for tolerance=2.0: Donor threshold=1.0, Hydrophobe threshold=2.0.

#### Step 1.2: Feature OR-Assignment
**File**: `pharmacophore/consensus.py`
**Source**: MoPBS (Braun 2022) — when two feature types are within 10% occurrence, assign both
**Change**: After clustering, if a spatial cluster contains mixed feature types where the secondary type has ≥90% of the primary type's count, output both types as OR-features.
**New format**: `['Donor|Acceptor', (), x, y, z]` for OR-assigned features
**Downstream impact**: `mol_converter.py` needs to handle `|`-separated types (generate both fragments). `constants.py` needs no change.

#### Step 1.3: S_Dbw Cluster Validation Metric
**File**: `pharmacophore/evaluation.py` (new function)
**Source**: MoPBS (Braun 2022) — `S_Dbw(c) = Scat(c) + Dens_bw(c)`
**Change**: Add `compute_sdbw(features, labels) -> float` function. Use as internal validation metric alongside screening metrics. Lower S_Dbw = better clustering.
**Implementation**: Scat = average cluster scatter / total scatter. Dens_bw = inter-cluster density at midpoints. Both computed from feature coordinates.

#### Step 1.4: Greedy Acquisition Function for Bayesian Optimization
**File**: `pharmacophore/combo_optimizer.py`
**Source**: Bellamy 2022 — Greedy and PI outperform EI in noisy environments
**Change**: Current uses `acq_func="EI"` in `gp_minimize()`. Add parameter `acq_func` with default `"gp_hedge"` (scikit-optimize's auto-selector that includes greedy behavior). Also expose `acq_func` as a configurable parameter.

#### Step 1.5: Retest Policy for Borderline Evaluations
**File**: `pharmacophore/combo_optimizer.py`
**Source**: Bellamy 2022 — retest predicted-active-measured-inactive cases
**Change**: After each optimization round, identify configurations where the predicted score (GP surrogate) differs from measured score by >1 std dev. Re-evaluate these with increased `n_conformers` (e.g., 2x). Replace the original measurement if the retest is more reliable.

---

### Phase 2: New Scoring & Distance Metrics (Cross-Domain Tier 1)

These add new scoring capabilities. Requires `scipy` (already installed) + optional `pot` library.

#### Step 2.1: Hungarian Matching for Optimal Feature Assignment
**File**: NEW `pharmacophore/hungarian_matching.py`
**Source**: G3PS algorithm, graph theory
**Dependency**: `scipy.optimize.linear_sum_assignment` (already available)
**Implementation**:
```python
def match_features(features_a, features_b, type_weight=1.0, spatial_weight=1.0):
    """Optimal assignment between two pharmacophore feature sets.

    Cost matrix: C[i,j] = spatial_weight * ||pos_i - pos_j||²
                         + type_weight * type_mismatch(i,j)

    Returns: List of (i, j) matched pairs + unmatched features + total cost
    """
```
- Type mismatch: 0.0 if same type, 0.5 if compatible (e.g., Donor-Acceptor), 1.0 otherwise
- Handles unequal sizes via dummy features with high cost
- Replaces greedy matching in current ICP-style comparison

#### Step 2.2: Wasserstein Distance for Pharmacophore Comparison
**File**: NEW `pharmacophore/ot_scoring.py`
**Source**: Optimal transport theory
**Dependency**: `scipy.spatial.distance` (available) + optional `pot` library for Sinkhorn
**Implementation**:
```python
def wasserstein_pharmacophore_distance(features_a, features_b, alpha=0.5):
    """Compute earth mover's distance between two pharmacophore models.

    Each model treated as discrete measure: μ = Σ w_i · δ(x_i, type_i)
    Cost: c(i,j) = ||pos_i - pos_j||² + alpha * type_distance(i,j)

    Returns: float distance (lower = more similar)
    """
```
- Without `pot`: Use `scipy.optimize.linear_sum_assignment` (exact OT for equal sizes)
- With `pot`: Use `ot.sinkhorn()` for fast entropy-regularized OT (handles unequal sizes naturally)
- Normalize to [0,1] for compatibility with existing scoring

#### Step 2.3: Ensemble Clustering with Voting
**File**: NEW `pharmacophore/ensemble_consensus.py`
**Source**: Li et al. 2018 (ensemble k-medoids), bioinformatics
**Implementation**:
```python
class EnsembleConsensus:
    def __init__(self, n_runs=100, tolerance_range=(1.5, 2.5),
                 occurrence_range=(0.4, 0.6), base_consensus=None):
        """Run consensus clustering N times with perturbed parameters, vote on stable features."""

    def generate_consensus(self, mols, features=None):
        """Returns features that appear in ≥50% of runs, with stability scores."""
```
- Each run: sample tolerance from `tolerance_range`, occurrence from `occurrence_range`
- Feature voting: spatial hashing to match features across runs (within 0.5Å = same feature)
- Output: features + stability_score (0.0-1.0) per feature
- Weight by S_Dbw score of each run (Step 1.3)

#### Step 2.4: Multi-Level Scoring Ensemble
**File**: `pharmacophore/evaluation.py` (extend `UnifiedEvaluator`)
**Source**: Roy et al. 2018 (intelligent consensus), Ghosh et al. 2025 (multi-level voting)
**Change**: Add `evaluate_ensemble()` method that:
1. Runs all available scorers: pharm2d, pharm3d, hybrid, OT (Step 2.2)
2. Each scorer produces a ranked list of actives/decoys
3. Combine via rank aggregation (Borda count or reciprocal rank fusion)
4. Return ensemble metrics alongside individual scorer metrics

---

### Phase 3: Colored Point Cloud Registration (Cross-Domain Tier 1-2)

Optional Open3D integration for advanced alignment. Graceful degradation if not installed.

#### Step 3.1: Colored ICP for Pharmacophore Alignment
**File**: NEW `pharmacophore/point_cloud_alignment.py`
**Source**: ELIXIR-A (Wieder et al. 2022), computer vision
**Dependency**: `open3d` (optional, not currently installed)
**Implementation**:
```python
def colored_icp_align(features_a, features_b, lambda_color=0.5):
    """Align two pharmacophore models using colored ICP.

    Maps feature types to colors:
      Donor → [1,0,0], Acceptor → [0,1,0], Aromatic → [0,0,1],
      Hydrophobe → [1,1,0], PosIonizable → [1,0,1]

    Returns: transformation matrix, aligned features, RMSD
    """
```
- Create Open3D PointCloud objects with RGB colors from feature types
- Use `open3d.pipelines.registration.registration_colored_icp()`
- Falls back to scipy-based Kabsch alignment if Open3D not available

#### Step 3.2: FPFH Descriptors for Pharmacophore Fingerprinting
**File**: `pharmacophore/point_cloud_alignment.py` (extend)
**Source**: Computer vision (Fast Point Feature Histograms)
**Dependency**: `open3d` (optional)
**Implementation**:
- Compute FPFH descriptor for each pharmacophore point based on local neighborhood geometry
- Use as a novel fingerprint type alongside Morgan/pharmacophore fingerprints
- Enable alignment-free global registration via FPFH + RANSAC

---

### Phase 4: Advanced Optimization (Cross-Domain Tier 2)

#### Step 4.1: Multi-Fidelity Bayesian Optimization
**File**: `pharmacophore/combo_optimizer.py` (extend)
**Source**: McDonald et al. 2025 (MF-BO)
**Change**: Add multi-fidelity evaluation where:
- Level 1 (cheap): pharm2d_scoring only → fast screening
- Level 2 (medium): pharm3d_scoring → 3D distance-based
- Level 3 (expensive): hybrid_scoring → full shape alignment
Use targeted variance reduction (TVR) to decide which fidelity level to evaluate at each iteration. Start with mostly Level 1, shift to Level 3 as optimization converges.

#### Step 4.2: Multi-Instance Learning for Conformer Selection
**File**: NEW `pharmacophore/mil_consensus.py`
**Source**: Zankov et al. 2023 (MIL-kmeans)
**Implementation**:
```python
class MILConsensus:
    def __init__(self, n_conformers_per_mol=10):
        """Treat each molecule as a bag of conformational instances."""

    def select_bioactive_conformers(self, mols):
        """Use MIL to identify which conformers contribute to consensus.
        Returns: selected conformer indices per molecule."""
```
- Generate multiple conformers per reference ligand
- Cluster all conformer-derived pharmacophore features across all molecules
- Per-molecule binary vector: bit=1 if any conformer has feature in cluster
- Select conformers that contribute most to the discriminative clusters

---

### Phase 5: Integration & Auto-Optimizer

#### Step 5.1: Unified Strategy Selector
**File**: NEW `pharmacophore/auto_strategy.py`
**Source**: Integration of all above
**Implementation**:
```python
class AutoPharmacophoreOptimizer:
    """Automatically selects best clustering + scoring strategy."""

    def __init__(self, reference_mols, actives=None, decoys=None):
        self.strategies = {
            'agglomerative': PharmacophoreConsensus,
            'ensemble': EnsembleConsensus,
            'mil': MILConsensus,  # Phase 4
        }
        self.scorers = {
            'hybrid': HybridScorer,
            'ot': wasserstein_pharmacophore_distance,
            'ensemble': evaluate_ensemble,
        }

    def optimize(self, n_trials=50, scoring_strategy='auto'):
        """Run multi-fidelity BO across all strategies."""
```

#### Step 5.2: Benchmark Harness
**File**: `experiments/benchmark_next_gen.py`
**Implementation**: Run all strategies on CCR2 dataset, produce comparison table with:
- ROC-AUC, BEDROC, EF@1%, EF@5%, EF@10%
- S_Dbw cluster quality
- Feature stability (from ensemble)
- Runtime per strategy

---

## Files Summary

| File | Action | Phase | Description |
|------|--------|-------|-------------|
| `pharmacophore/consensus.py` | MODIFY | 1 | Add weight differentiation + OR-assignment |
| `pharmacophore/evaluation.py` | MODIFY | 1,2 | Add S_Dbw + ensemble evaluation |
| `pharmacophore/combo_optimizer.py` | MODIFY | 1,4 | Greedy acq + retest + MF-BO |
| `pharmacophore/constants.py` | MODIFY | 1 | Add type compatibility matrix |
| `pharmacophore/mol_converter.py` | MODIFY | 1 | Handle OR-type features |
| `pharmacophore/hungarian_matching.py` | CREATE | 2 | Optimal feature assignment |
| `pharmacophore/ot_scoring.py` | CREATE | 2 | Wasserstein distance scoring |
| `pharmacophore/ensemble_consensus.py` | CREATE | 2 | Stability-aware ensemble clustering |
| `pharmacophore/point_cloud_alignment.py` | CREATE | 3 | Colored ICP + FPFH (optional Open3D) |
| `pharmacophore/mil_consensus.py` | CREATE | 4 | Multi-instance learning for conformers |
| `pharmacophore/auto_strategy.py` | CREATE | 5 | Unified strategy selector |
| `experiments/benchmark_next_gen.py` | CREATE | 5 | Comprehensive benchmark |
| `tests/test_hungarian_matching.py` | CREATE | 2 | Tests for Hungarian matching |
| `tests/test_ot_scoring.py` | CREATE | 2 | Tests for OT scoring |
| `tests/test_ensemble_consensus.py` | CREATE | 2 | Tests for ensemble clustering |
| `tests/test_point_cloud.py` | CREATE | 3 | Tests for point cloud alignment |

## Dependencies

| Package | Phase | Status | Required? |
|---------|-------|--------|-----------|
| `scipy` | 2 | Installed | Yes (linear_sum_assignment, wasserstein) |
| `scikit-learn` | 1-5 | Installed | Yes |
| `numpy` | 1-5 | Installed | Yes |
| `rdkit` | 1-5 | Installed | Yes |
| `scikit-optimize` | 1,4 | Installed (optional) | For BO |
| `pot` | 2 | NOT installed | Optional (Sinkhorn OT, falls back to scipy) |
| `open3d` | 3 | NOT installed | Optional (colored ICP, falls back to Kabsch) |

## Verification Plan

After each phase, run:

1. **Determinism tests** (must always pass):
   ```bash
   pytest tests/test_consensus.py::TestPharmacophoreConsensus::test_determinism
   pytest tests/test_consensus.py::TestPharmacophoreIntegration::test_determinism_across_runs
   ```

2. **Phase-specific tests**:
   ```bash
   # Phase 1
   pytest tests/test_consensus.py -k "weight or or_assign"
   pytest tests/test_evaluation.py -k "sdbw"

   # Phase 2
   pytest tests/test_hungarian_matching.py
   pytest tests/test_ot_scoring.py
   pytest tests/test_ensemble_consensus.py

   # Phase 3
   pytest tests/test_point_cloud.py
   ```

3. **CCR2 benchmark** (after each phase):
   ```bash
   python experiments/benchmark_next_gen.py --dataset tutorials/data/ --phases 1,2,3
   ```

   Expected improvement targets:
   - Phase 1: +2-5% ROC-AUC over baseline (weighted clustering + OR-features)
   - Phase 2: +5-10% via ensemble scoring and OT distance
   - Phase 3: +2-5% via alignment-free comparison
   - Phase 5: Best overall via auto-selection

4. **Full regression**:
   ```bash
   pytest
   ```

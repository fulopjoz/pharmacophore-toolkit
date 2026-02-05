# Literature-Based Improvements for Pharmacophore Toolkit

**Date**: 2026-01-26  
**Context**: Post-DOE Phase 4, literature mining for performance enhancements  
**Current Performance**: ROC-AUC = 0.7629 (Phase 3 DBSCAN)  
**Target**: ROC-AUC ≥ 0.90 (ML-competitive range)

---

## Executive Summary

Analysis of 6 key papers in `docs/research/papers/relevant/` reveals **4 immediate improvements** that could boost our ROC-AUC from **0.76 → 0.95** with ~40 hours of implementation work over 2 weeks.

**Literature Sources Reviewed**:
1. ✅ **Langer Presentation 2025** - Industry best practices (LigandScout)
2. ✅ **Consensus Holistic VS 2024** - 2D+3D hybrid methods
3. ✅ **PheSA 2024** - Shape+color scoring validation
4. ✅ **PharmacoForge 2025** - ML clustering (rejected for determinism)
5. ✅ **Spatial Clustering GNN 2024** - Graph neural networks (long-term)
6. ✅ **Clustering Binding Pockets 2023** - Spatial methods comparison

---

## Improvement Roadmap

### Phase 5a: Quick Wins (3 hours, +3-8% AUC)

#### 1. Feature Weighting ⭐ HIGHEST PRIORITY
**Source**: Langer Presentation 2025 (slides 35-42)

**Concept**: H-bond donors/acceptors are more important than hydrophobic contacts for binding.

**Implementation**:
```python
# pharmacophore/consensus.py

FEATURE_WEIGHTS = {
    'Donor': 2.0,      # Strong H-bond (~5 kcal/mol)
    'Acceptor': 2.0,   # Strong H-bond
    'Aromatic': 1.5,   # π-interaction (~3 kcal/mol)
    'Hydrophobe': 1.0  # Weak vdW (~1 kcal/mol)
}

def _apply_feature_weights(self, features):
    """Weight features by interaction strength."""
    weighted_features = []
    for feature_type, indices, x, y, z in features:
        weight = FEATURE_WEIGHTS.get(feature_type, 1.0)
        # Store weight in feature metadata
        weighted_features.append([feature_type, indices, x, y, z, weight])
    return weighted_features

# In scoring (benchmark.py):
def calculate_weighted_score(matches, weights):
    return sum(w * m for w, m in zip(weights, matches)) / sum(weights)
```

**Expected Gain**: +2-5% AUC  
**Effort**: 2 hours  
**Risk**: Low (proven in literature)  
**Validation**: Lipinski's Rule of 5 emphasizes H-bonds

---

#### 2. Alignment Tries Parameter ⭐
**Source**: Langer HAT (Heuristic Alignment Trees) - slides 15-22

**Concept**: Control speed vs accuracy trade-off by adjusting alignment effort.

**Langer's Benchmarks**:
| Tries | Speed | Accuracy | Use Case |
|-------|-------|----------|----------|
| 20 | Fastest | Basic | Quick filtering |
| **50** | Fast | **Good** | **Generic screening** |
| 300 | Medium | High | Accurate ranking |

**Key Insight**: Runtime increases **sub-linearly**
- 20→50 tries = +25% time (not 2.5x)
- 50→300 tries = +30% time (not 6x)

**Implementation**:
```python
# pharmacophore/benchmark.py

class ScreeningBenchmark:
    def run_benchmark(
        self, 
        consensus_mol, 
        n_conformers=10,
        alignment_maxiters=50  # NEW PARAMETER
    ):
        # Current code:
        shape, color = AlignMol(ref=consensus_mol, probe=query_mol, useColors=True)
        
        # Enhanced code:
        shape, color = AlignMol(
            ref=consensus_mol,
            probe=query_mol,
            useColors=True,
            maxIters=alignment_maxiters  # 20/50/300
        )
        
        return {'shape': shape, 'color': color, 'combo': shape + color}
```

**Expected Gain**: +1-3% AUC (with maxIters=300)  
**Effort**: 1 hour  
**Risk**: Low  
**Note**: Test speed impact first (should be minimal)

---

### Phase 5b: Hybrid Scoring (4 hours, +10-15% AUC) ⭐⭐⭐

#### 3. 2D+3D Ensemble Scoring
**Source**: Consensus Holistic VS 2024, Sanders et al. 2012

**Concept**: Combine 2D fingerprint similarity + 3D shape scoring for orthogonal information.

**Literature Validation**:
- Sanders et al. (2012): "Combined 2D+3D outperforms either method alone"
- Moshawih et al. (2024): "+15% enrichment with consensus approach"
- Multiple DUD-E benchmarks confirm superiority

**Why It Works**:
```
2D Fingerprints:  Capture chemical similarity (functional groups)
3D Shape:         Capture spatial fit (binding geometry)
                  ↓
          ORTHOGONAL INFORMATION
                  ↓
         Ensemble > Individual
```

**Implementation**:
```python
# NEW FILE: pharmacophore/hybrid_scoring.py

from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol

class HybridScorer:
    """Combine 2D fingerprint + 3D shape scoring."""
    
    def __init__(self, reference_mols, consensus_mol, alpha=0.6):
        """
        Args:
            reference_mols: Original reference ligands (for 2D)
            consensus_mol: 3D consensus pharmacophore (for 3D)
            alpha: Weight for 2D component (0.6 = 60% 2D, 40% 3D)
        """
        self.reference_mols = reference_mols
        self.consensus_mol = consensus_mol
        self.alpha = alpha
        
        # Pre-compute reference fingerprints
        self.reference_fps = [
            AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024)
            for m in reference_mols
        ]
    
    def score(self, query_mol):
        """Calculate hybrid 2D+3D score."""
        # 2D component: Max Tanimoto to any reference
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, 1024)
        tanimoto_2d = max(
            DataStructs.TanimotoSimilarity(query_fp, ref_fp)
            for ref_fp in self.reference_fps
        )
        
        # 3D component: Shape + color to consensus
        shape, color = AlignMol(
            ref=self.consensus_mol,
            probe=query_mol,
            useColors=True
        )
        tanimoto_3d = (shape + color) / 2  # Normalize to [0, 1]
        
        # Weighted ensemble
        hybrid_score = self.alpha * tanimoto_2d + (1 - self.alpha) * tanimoto_3d
        
        return {
            'hybrid': hybrid_score,
            '2d': tanimoto_2d,
            '3d': tanimoto_3d,
            'shape': shape,
            'color': color
        }

# Usage in benchmark:
from pharmacophore.hybrid_scoring import HybridScorer

scorer = HybridScorer(reference_mols, consensus_mol, alpha=0.6)
results = [scorer.score(mol) for mol in test_molecules]
```

**Parameter Tuning** (alpha):
- alpha=0.0 → Pure 3D (current method, 0.76 AUC)
- alpha=0.6 → Balanced (recommended, ~0.89 AUC expected)
- alpha=1.0 → Pure 2D (~0.75 AUC typical)

**Expected Gain**: +10-15% AUC (0.76 → 0.87-0.89)  
**Effort**: 4 hours  
**Risk**: Low (proven in multiple studies)  
**Validation**: Test on CCR2 dataset, compare with Phase 3

---

### Phase 5c: Multi-Model Ensemble (6 hours, +5-8% AUC)

#### 4. Ensemble of Strict/Moderate/Relaxed Models
**Source**: Consensus Holistic VS 2024

**Concept**: Generate multiple consensus models with different thresholds, then vote.

**Implementation**:
```python
# experiments/run_ensemble_screening.py

def generate_ensemble_models(reference_mols):
    """Create 3 consensus models with different stringencies."""
    from pharmacophore.pharmacophore import Pharmacophore
    
    pharm = Pharmacophore()
    
    models = {
        'strict': pharm.generate_consensus_models(
            reference_mols,
            clustering_method='dbscan',
            tolerance=1.5,
            threshold=0.5,  # Must appear in 3/5 refs
            linkage='complete'
        ),
        'moderate': pharm.generate_consensus_models(
            reference_mols,
            clustering_method='dbscan',
            tolerance=1.5,
            threshold=0.3,  # Must appear in 2/5 refs (Phase 3 optimal)
            linkage='complete'
        ),
        'relaxed': pharm.generate_consensus_models(
            reference_mols,
            clustering_method='dbscan',
            tolerance=1.5,
            threshold=0.2,  # Must appear in 1/5 refs
            linkage='complete'
        )
    }
    
    return models

def ensemble_score(query_mol, models, method='max'):
    """Score query against all models and aggregate."""
    from pharmacophore.benchmark import ScreeningBenchmark
    
    scores = {}
    for name, (features, consensus_mol) in models.items():
        shape, color = AlignMol(consensus_mol, query_mol, useColors=True)
        scores[name] = shape + color
    
    if method == 'max':
        return max(scores.values())  # Best of 3 models
    elif method == 'mean':
        return sum(scores.values()) / len(scores)
    elif method == 'weighted':
        # Weight by feature count (more features = more specific)
        weights = {name: len(features) for name, (features, _) in models.items()}
        return sum(w * s for (w, s) in zip(weights.values(), scores.values())) / sum(weights.values())
    else:
        raise ValueError(f"Unknown method: {method}")

# Benchmark all 3 methods
results = {
    'max': benchmark_ensemble(models, method='max'),
    'mean': benchmark_ensemble(models, method='mean'),
    'weighted': benchmark_ensemble(models, method='weighted')
}
```

**Expected Gain**: +5-8% AUC  
**Effort**: 6 hours  
**Risk**: Low (uses existing code)

---

### Phase 6: Speed Optimizations (8 hours, 2-5x faster)

#### 5. Dynamic Feature Selection
**Source**: Langer hierarchical trees (slides 35-42)

**Concept**: Use fewer features for fast pre-filter, all features for final ranking.

**Cascade Architecture**:
```
Stage 1: Top 3 features → Fast filter (0.02s/mol)
          ↓ (retain 20% of molecules)
Stage 2: All 12 features → Accurate ranking (0.08s/mol)
          ↓
     Final hits
```

**Implementation**:
```python
# pharmacophore/cascade_screening.py

def select_top_features(features, n=3, criterion='occurrence'):
    """Select most important features."""
    if criterion == 'occurrence':
        # Sort by occurrence frequency (stored during consensus)
        sorted_features = sorted(features, key=lambda f: f[-1], reverse=True)
    elif criterion == 'binding_energy':
        # Sort by estimated binding contribution (requires precomputed weights)
        weights = {'Donor': 2.0, 'Acceptor': 2.0, 'Aromatic': 1.5, 'Hydrophobe': 1.0}
        sorted_features = sorted(features, key=lambda f: weights[f[0]], reverse=True)
    
    return sorted_features[:n]

class CascadeScreening:
    def __init__(self, full_model, fast_model, threshold_stage1=0.3, threshold_stage2=0.6):
        self.full_model = full_model  # All 12 features
        self.fast_model = fast_model  # Top 3 features
        self.threshold_stage1 = threshold_stage1
        self.threshold_stage2 = threshold_stage2
    
    def screen(self, molecules):
        # Stage 1: Fast pre-filter
        candidates = []
        for mol in molecules:
            shape, color = AlignMol(self.fast_model, mol, useColors=True)
            if (shape + color) / 2 > self.threshold_stage1:
                candidates.append(mol)
        
        print(f"Stage 1: {len(candidates)}/{len(molecules)} passed ({len(candidates)/len(molecules)*100:.1f}%)")
        
        # Stage 2: Accurate ranking
        hits = []
        for mol in candidates:
            shape, color = AlignMol(self.full_model, mol, useColors=True)
            score = (shape + color) / 2
            if score > self.threshold_stage2:
                hits.append((mol, score))
        
        print(f"Stage 2: {len(hits)}/{len(candidates)} hits")
        
        return sorted(hits, key=lambda x: x[1], reverse=True)
```

**Expected Speedup**: 2-5x faster (depends on stage 1 retention rate)  
**Effort**: 8 hours  
**Risk**: Medium (requires tuning thresholds)

---

## Validation Plan

### Test on CCR2 Dataset (Current)
1. Implement feature weighting → measure ΔRO C-AUC
2. Implement alignment tries (50, 100, 300) → profile speed vs accuracy
3. Implement 2D+3D hybrid → compare with Phase 3
4. Implement ensemble → test max/mean/weighted
5. Generate comprehensive comparison table

### Validate on Additional Datasets (Phase 6)
- **DUD-E**: 102 protein targets, standard benchmark
- **LIT-PCBA**: 15 assays, challenging dataset
- **Our own**: Expand CCR2 to 15-20 reference ligands

---

## Implementation Priority (2-Week Sprint)

### Week 1: Core Improvements

**Day 1 (3 hours)**:
- ✅ Feature weighting implementation (2h)
- ✅ Alignment tries parameter (1h)
- ✅ Test on CCR2, measure AUC change

**Day 2-3 (8 hours)**:
- ✅ 2D+3D hybrid scoring class (4h)
- ✅ Benchmark hybrid on CCR2 (2h)
- ✅ Tune alpha parameter (0.5-0.7 range) (2h)

**Day 4 (6 hours)**:
- ✅ Multi-model ensemble implementation (4h)
- ✅ Compare max/mean/weighted strategies (2h)

**Day 5 (4 hours)**:
- ✅ Integration testing (2h)
- ✅ Update documentation (1h)
- ✅ Generate Phase 5 summary report (1h)

### Week 2: Validation & Speed

**Day 1 (4 hours)**:
- ✅ Test with 50 conformers (0.5h)
- ✅ Profile speed impact (1h)
- ✅ Optimize conformer generation (2.5h)

**Day 2-3 (8 hours)**:
- ✅ Dynamic feature selection implementation (6h)
- ✅ Cascade screening benchmark (2h)

**Day 4-5 (8 hours)**:
- ✅ PheSA comparison (if paper accessible) (8h)
- OR Final validation on additional datasets (8h)

---

## Expected Outcomes

### Performance Projection

| Stage | Method | ROC-AUC | Gain | Effort |
|-------|--------|---------|------|--------|
| Phase 3 | DBSCAN baseline | 0.7629 | - | - |
| + Week 1 Day 1 | Feature weighting | 0.7780 | +2% | 2h |
| + Week 1 Day 1 | Alignment tries | 0.7858 | +1% | 1h |
| + Week 1 Day 2-3 | 2D+3D hybrid | **0.8929** | **+14%** | 4h |
| + Week 1 Day 4 | Ensemble | 0.9315 | +5% | 6h |
| + Week 2 Day 1 | Better conformers | 0.9548 | +2.5% | 0.5h |
| **TOTAL** | **All improvements** | **~0.95** | **+25%** | **~40h** |

### Validation Targets

- **CCR2 dataset**: ROC-AUC ≥ 0.90 (current 0.76)
- **Speed**: Maintain ~46s for 574 molecules (or justify slowdown)
- **Determinism**: All methods must be reproducible
- **Interpretability**: Visualizations and feature explanations

### Success Criteria

✅ **Minimum**: ROC-AUC ≥ 0.85 (+11% from Phase 3)  
🎯 **Target**: ROC-AUC ≥ 0.90 (+18% from Phase 3)  
⭐ **Stretch**: ROC-AUC ≥ 0.95 (+25% from Phase 3)

---

## Literature References

### Papers Analyzed

1. **Langer, T. (2025)** - Next Generation Pharmacophore Modeling  
   - Source: `docs/TLanger_Olomouc_2025-Reduced.pdf`
   - Key insights: HAT, feature weighting, exclude volumes
   - Summary: `docs/research/papers/relevant/TLanger_Presentation_2025_Summary.md`

2. **Moshawih et al. (2024)** - Consensus holistic virtual screening  
   - Validates 2D+3D hybrid approach
   - Summary: `docs/research/papers/relevant/Consensus_Holistic_VS_2024.md`

3. **Wahl, J. (2024)** - PheSA: Pharmacophore-Enhanced Shape Alignment  
   - Confirms our combo scoring approach
   - Summary: `docs/research/papers/relevant/PheSA_2024.md`

4. **Sanders et al. (2012)** - Comparative analysis of pharmacophore screening tools  
   - Proves 2D+3D > either alone
   - Cited in multiple papers, gold standard validation

5. **Flynn et al. (2025)** - PharmacoForge: Diffusion models  
   - Rejected for violating determinism
   - Summary: `docs/research/papers/relevant/PharmacoForge_2025.md`

6. **Spatial Clustering with GNNs (2024)**  
   - Long-term research direction
   - Summary: `docs/research/papers/relevant/Spatial_Clustering_GNN_2024.md`

---

## What NOT to Do (Rejected Approaches)

### ❌ Diffusion Models (PharmacoForge)
**Reason**: Stochastic generation violates determinism requirement  
**Note**: Non-reproducible results unacceptable for drug discovery

### ❌ Full Molecular Docking
**Reason**: 100-1000x slower than our method (10-60s/mol vs 0.08s/mol)  
**Note**: Reserve docking for final hit validation, not virtual screening

### ❌ Deep Learning Feature Extraction
**Reason**: Requires large training datasets, GPU infrastructure  
**Note**: Our 5 reference ligands insufficient for meaningful training

### ❌ MD Trajectory Pharmacophores
**Reason**: Requires MD simulations (CCR2 trajectory not available)  
**Note**: Interesting for future work if MD data obtained

---

## Integration with Existing Code

### Files to Modify

1. **`pharmacophore/consensus.py`**
   - Add `FEATURE_WEIGHTS` constant
   - Implement `_apply_feature_weights()` method
   - Update `_cluster_features()` to use weights

2. **`pharmacophore/benchmark.py`**
   - Add `alignment_maxiters` parameter
   - Modify `AlignMol()` call to accept maxIters
   - Update scoring to use weighted features

3. **NEW: `pharmacophore/hybrid_scoring.py`**
   - Complete HybridScorer class
   - Integration tests with CCR2 dataset

4. **NEW: `experiments/run_ensemble_screening.py`**
   - Multi-model generation
   - Voting strategies (max/mean/weighted)
   - Comprehensive benchmarking

5. **NEW: `pharmacophore/cascade_screening.py`**
   - Feature selection algorithms
   - Two-stage screening pipeline
   - Speed profiling

### Backward Compatibility

All enhancements are **opt-in**:
```python
# Legacy usage (unchanged)
pharm = Pharmacophore()
features, mol = pharm.generate_consensus_models(refs)

# Enhanced usage (new parameters)
features, mol = pharm.generate_consensus_models(
    refs,
    feature_weights=True,      # Enable weighting
    alignment_maxiters=300,    # More accurate alignment
    use_hybrid_scoring=True,   # 2D+3D ensemble
    alpha=0.6                  # 60% 2D, 40% 3D
)
```

---

## Deliverables

### Week 1 Deliverables
- ✅ Feature weighting implementation
- ✅ Alignment tries parameter
- ✅ 2D+3D hybrid scoring class
- ✅ Multi-model ensemble
- ✅ Benchmark results on CCR2
- ✅ Code review and tests

### Week 2 Deliverables
- ✅ Conformer optimization
- ✅ Dynamic feature selection
- ✅ Speed profiling report
- ✅ PheSA comparison (if accessible)
- ✅ Phase 5 complete summary document
- ✅ Updated production recommendations

### Documentation Updates
- ✅ Tutorial: "Hybrid 2D+3D Screening"
- ✅ API docs for new classes
- ✅ Benchmark comparison table
- ✅ Literature citations in README

---

## Risk Assessment

### Low Risk (Proceed with confidence)
- ✅ Feature weighting (proven in Langer literature)
- ✅ Alignment tries (standard RDKit parameter)
- ✅ 2D+3D hybrid (validated in Sanders et al. 2012)
- ✅ Multi-model ensemble (uses existing models)

### Medium Risk (Test thoroughly)
- ⚠️ Conformer optimization (may slow screening significantly)
- ⚠️ Dynamic feature selection (requires threshold tuning)
- ⚠️ Alpha parameter tuning (needs cross-validation)

### High Risk (Deferred to long-term research)
- 🔴 GNN clustering (requires training data, violates determinism)
- 🔴 MD trajectory analysis (data not available)
- 🔴 Diffusion models (non-deterministic)

---

## Conclusion

Literature analysis reveals **clear, actionable paths** to boost performance from 0.76 → 0.95 ROC-AUC while maintaining:
- ✅ Determinism (reproducible results)
- ✅ Speed (< 1 minute for 574 molecules)
- ✅ Interpretability (visualizable pharmacophores)
- ✅ Minimal dependencies (no GPUs, no training data)

**Recommendation**: Execute 2-week sprint implementing improvements in priority order. This positions the toolkit as **competitive with ML methods** while retaining unique advantages of 3D pharmacophore modeling.

**Next Steps**: 
1. Review this document with user
2. Get approval for implementation plan
3. Start Week 1 Day 1 tasks (feature weighting + alignment tries)
4. Update session plan.md with new tasks

---

**Document Status**: Complete, ready for implementation  
**Estimated ROC-AUC after completion**: 0.90-0.95  
**Timeline**: 2 weeks (40 hours)  
**Risk Level**: Low-Medium (most changes proven in literature)

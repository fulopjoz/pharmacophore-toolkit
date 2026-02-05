# rdShapeAlign Consensus Pharmacophore Optimization Plan

## Objective

Create an optimized algorithm using RDKit's `rdShapeAlign` Gaussian shape matching to build consensus pharmacophore models from reference molecules that achieve maximum discrimination between actives and decoys (ROC-AUC > 0.85).

---

## 1. Understanding rdShapeAlign

### 1.1 How rdShapeAlign Works

RDKit's `rdShapeAlign.AlignMol` is a wrapper around PubChem's `pubchem-align3d` library that performs:

1. **Gaussian Volume Overlap**: Represents molecular shape as atom-centered Gaussian functions
2. **Color Feature Matching**: Uses SMARTS-based pharmacophore features for "color" scoring
3. **Two-Phase Optimization**: Fast global search + refined local optimization

### 1.2 Key Parameters (from context7 documentation)

```python
from rdkit.Chem.rdShapeAlign import AlignMol, PrepareConformer

shape_score, color_score = AlignMol(
    ref=reference_mol,           # Reference molecule
    probe=query_mol,             # Query molecule to align
    refConfId=-1,                # Conformer ID (-1 = first)
    probeConfId=-1,              # Conformer ID (-1 = first)
    useColors=True,              # Enable pharmacophore color features
    opt_param=0.5,               # Balance: 0=color only, 0.5=equal, 1.0=shape only
    max_preiters=10,             # Phase 1: iterations on all poses (DEFAULT: 10)
    max_postiters=30             # Phase 2: iterations on best poses (DEFAULT: 30)
)

combo_tanimoto = shape_score + color_score  # Range: 0-2
```

### 1.3 Color Feature Definitions

rdShapeAlign uses these SMARTS patterns for color features (from RDKit source):

| Feature Type | SMARTS Pattern | Our Fragment |
|--------------|----------------|--------------|
| Donor | `[#7!H0]`, `[#8!H0]` | `[NH3]` - matches `[#7!H0]` |
| Acceptor | `[#7]`, `[#8]` | `[O]` - matches `[#8]` |
| Aromatic | `a1aaaaa1` | `c1ccccc1` - matches aromatic |
| Hydrophobe | Various C patterns | `[CH4]` - matches hydrophobe |
| PosIonizable | `[+]`, `[#7]` | `[NH4+]` - matches positive |

**Key Insight**: Our disconnected fragments SHOULD be recognized by the SMARTS patterns. The issue is likely in parameter tuning or consensus generation, not the fragment representation.

### 1.4 Why Disconnected Fragments Should Work

1. **Gaussian overlap is local**: Each atom contributes its own Gaussian - no connectivity required
2. **SMARTS matching is per-atom**: Color features are detected independently on each atom/fragment
3. **Tanimoto normalization**: Handles size differences between probe and reference

---

## 2. Current Implementation Analysis

### 2.1 Key Files

| File | Purpose | Location |
|------|---------|----------|
| `pharmacophore/mol_converter.py` | Converts features to mol fragments | Lines 80-282 |
| `pharmacophore/consensus.py` | Generates consensus pharmacophore | Full file |
| `pharmacophore/optimization.py` | Grid search optimizer | Lines 71-121 (scoring) |
| `pharmacophore/benchmark.py` | Benchmarking framework | Lines 157-197 (scoring) |
| `pharmacophore/auto_optimizer.py` | Bayesian optimization | Lines 450-480 |
| `pharmacophore/doe_optimization.py` | DOE-based optimization | Lines 130-145 |

### 2.2 Current Scoring Implementation

From `pharmacophore/optimization.py:106-121`:
```python
def score_molecule(self, mol, pharmacophore_mol, use_colors=True):
    conformers = self._get_conformers(mol)
    best_combo = 0.0

    for conf_mol in conformers:
        shape, color = AlignMol(
            ref=pharmacophore_mol,
            probe=conf_mol,
            useColors=use_colors,
            opt_param=0.5
        )
        combo = shape + color
        if combo > best_combo:
            best_combo = combo

    return best_shape, best_color, best_combo
```

### 2.3 Current Results

| Method | ROC-AUC | Notes |
|--------|---------|-------|
| Consensus pharmacophore → Shape+Color | 0.47-0.65 | Poor discrimination |
| Reference alignment (direct) | 0.78 | Better - uses connected mols |
| Pharm2D fingerprints | 0.93 | Best 2D method |
| Hybrid (70% 2D + 30% 3D) | 0.95 | Current best |

### 2.4 Hypothesis: Why Current 3D Underperforms

1. **Parameter defaults not optimal**: `max_preiters=10` may be too few starting poses
2. **Consensus parameters**: tolerance/occurrence may create suboptimal feature selection
3. **Fragment positioning**: Small fragments may have too little overlap volume
4. **Missing optimization**: We're not tuning shape vs color weight (`opt_param`)

---

## 3. Implementation Plan

### 3.1 Phase 1: Parameter Optimization Study

**Goal**: Find optimal rdShapeAlign parameters for pharmacophore discrimination

```python
# Parameters to optimize:
PARAM_GRID = {
    'opt_param': [0.0, 0.25, 0.5, 0.75, 1.0],  # Shape vs color balance
    'max_preiters': [10, 20, 50, 100],          # Phase 1 iterations
    'max_postiters': [30, 50, 100, 200],        # Phase 2 iterations
    'n_conformers': [5, 10, 20, 50],            # Conformers per molecule
}
```

**Expected outcome**: Identify best parameter combination for AUC improvement

### 3.2 Phase 2: Consensus Generation Optimization

**Goal**: Optimize consensus pharmacophore generation for shape alignment

```python
# Consensus parameters to optimize:
CONSENSUS_GRID = {
    'tolerance': [1.0, 1.5, 2.0, 2.5, 3.0],    # Spatial clustering (Angstroms)
    'occurrence_threshold': [0.3, 0.5, 0.7],    # Feature frequency
    'min_features': [3, 5, 7],                   # Minimum features
    'max_features': [10, 15, 20],                # Maximum features
}
```

### 3.3 Phase 3: Enhanced Fragment Representation

**Goal**: Test whether enhanced fragments improve color recognition

```python
# Alternative fragment representations:
ENHANCED_FRAGMENTS = {
    'Donor': '[NH2]C',         # Amine with carbon context
    'Acceptor': 'O=C',         # Carbonyl oxygen
    'Aromatic': 'c1ccc(C)cc1', # Toluene (larger aromatic)
    'Hydrophobe': 'CC(C)C',    # Isopropyl (larger hydrophobe)
}
```

### 3.4 Phase 4: Multi-Reference Ensemble

**Goal**: Use ensemble of reference alignments instead of single consensus

```python
def ensemble_score(query_mol, reference_mols, method='max'):
    """Score query against all references, aggregate results."""
    scores = []
    for ref in reference_mols:
        shape, color = AlignMol(ref=ref, probe=query_mol, ...)
        scores.append(shape + color)

    if method == 'max':
        return max(scores)
    elif method == 'mean':
        return np.mean(scores)
    elif method == 'weighted':
        return np.average(scores, weights=ref_weights)
```

---

## 4. Technical Implementation

### 4.1 New Module: `pharmacophore/rdshape_optimizer.py`

```python
"""
Optimized rdShapeAlign-based scoring for consensus pharmacophores.

This module provides:
1. Parameter-optimized AlignMol scoring
2. PrepareConformer caching for speed
3. Multi-reference ensemble scoring
4. Bayesian optimization for consensus parameters
"""

from rdkit.Chem.rdShapeAlign import AlignMol, PrepareConformer
from typing import List, Tuple, Dict, Optional
from rdkit import Chem
import numpy as np

class RDShapeOptimizer:
    """Optimized rdShapeAlign scoring with parameter tuning."""

    def __init__(
        self,
        opt_param: float = 0.5,
        max_preiters: int = 50,
        max_postiters: int = 100,
        n_conformers: int = 20,
        use_prepared_shapes: bool = True
    ):
        self.opt_param = opt_param
        self.max_preiters = max_preiters
        self.max_postiters = max_postiters
        self.n_conformers = n_conformers
        self.use_prepared_shapes = use_prepared_shapes
        self._shape_cache = {}

    def prepare_reference(self, mol: Chem.Mol) -> 'ShapeInput':
        """Pre-compute shape for faster repeated alignments."""
        smiles = Chem.MolToSmiles(mol)
        if smiles not in self._shape_cache:
            self._shape_cache[smiles] = PrepareConformer(mol)
        return self._shape_cache[smiles]

    def score_molecule(
        self,
        probe: Chem.Mol,
        reference: Chem.Mol,
        use_colors: bool = True
    ) -> Tuple[float, float, float]:
        """Score probe against reference with optimized parameters."""
        shape, color = AlignMol(
            ref=reference,
            probe=probe,
            useColors=use_colors,
            opt_param=self.opt_param,
            max_preiters=self.max_preiters,
            max_postiters=self.max_postiters
        )
        return shape, color, shape + color
```

### 4.2 Integration with Existing Code

Modify `pharmacophore/benchmark.py` to use optimized parameters:

```python
# In ScreeningBenchmark._score_single_molecule()
shape, color = AlignMol(
    ref=pharm_mol,
    probe=conf_mol,
    useColors=use_colors,
    opt_param=self.opt_param,           # NEW: Configurable
    max_preiters=self.max_preiters,     # NEW: Configurable
    max_postiters=self.max_postiters    # NEW: Configurable
)
```

### 4.3 Validation Protocol

```python
def validate_rdshape_scoring(
    reference_mols: List[Chem.Mol],
    actives: List[Chem.Mol],
    decoys: List[Chem.Mol],
    param_grid: Dict
) -> pd.DataFrame:
    """
    Systematic validation of rdShapeAlign parameters.

    Returns DataFrame with columns:
    - opt_param, max_preiters, max_postiters, tolerance, occurrence
    - roc_auc, bedroc, ef_1, ef_5, time_per_mol
    """
    results = []

    for params in itertools.product(*param_grid.values()):
        param_dict = dict(zip(param_grid.keys(), params))

        # Generate consensus with these params
        consensus = PharmacophoreConsensus(
            tolerance=param_dict['tolerance'],
            occurrence_threshold=param_dict['occurrence']
        )
        features = consensus.generate_consensus(reference_mols)
        pharm_mol = PharmacophoreToMol.convert(features, enable_color_features=True)

        # Score all molecules
        scores = []
        for mol in actives + decoys:
            shape, color = AlignMol(
                ref=pharm_mol,
                probe=mol,
                opt_param=param_dict['opt_param'],
                max_preiters=param_dict['max_preiters'],
                max_postiters=param_dict['max_postiters']
            )
            scores.append(shape + color)

        # Calculate metrics
        y_true = [1] * len(actives) + [0] * len(decoys)
        auc = roc_auc_score(y_true, scores)

        results.append({**param_dict, 'roc_auc': auc})

    return pd.DataFrame(results)
```

---

## 5. Data and Test Files

### 5.1 Reference Data Location

```
tutorials/data/
├── CCR2_reference_ligands.sdf   # 5 aligned reference ligands
├── actives_ccr2_N75.csv         # 75 active compounds (SMILES)
├── decoys_ccr2_N500.csv         # 500 property-matched decoys
└── CCR_HUMAN_AL.tsv             # Additional alignment data
```

### 5.2 Test Script Location

```
tests/
├── test_optimization.py         # Existing optimization tests
├── test_consensus.py            # Consensus generation tests
└── test_benchmark.py            # Benchmark framework tests
```

---

## 6. Success Criteria

| Metric | Current | Target | Stretch Goal |
|--------|---------|--------|--------------|
| ROC-AUC | 0.47-0.65 | >0.80 | >0.90 |
| BEDROC (α=20) | 0.15 | >0.40 | >0.60 |
| EF@1% | 2.0 | >10.0 | >20.0 |
| EF@5% | 1.5 | >5.0 | >10.0 |
| Time per mol | 50ms | <100ms | <50ms |

---

## 7. Risk Analysis

| Risk | Mitigation |
|------|------------|
| Disconnected fragments truly incompatible | Fall back to reference alignment method |
| Parameter optimization finds no improvement | Test enhanced fragment representations |
| Computation too slow | Use PrepareConformer caching |
| AUC still below 2D method | Accept hybrid as production solution |

---

## 8. References

### 8.1 RDKit Documentation
- [rdShapeAlign API](https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeAlign.html)
- [Greg Landrum's Blog - 3D Thresholds](https://greglandrum.github.io/rdkit-blog/posts/2025-11-30-thresholds-for-random-3d.html)
- [RDKit 2024.09 Release Notes](https://greglandrum.github.io/rdkit-blog/posts/2024-10-10-whats-new-2024-09-2.html)

### 8.2 Commercial Approaches
- [OpenEye ROCS](https://www.eyesopen.com/rocs) - Gaussian shape + color overlay
- [PheSA](https://openmolecules.org/help/phesa.html) - Open source pharmacophore-enhanced shape alignment

### 8.3 Literature
- Hönig et al. (2022) "Small molecule superposition: A comprehensive overview" [DOI:10.1002/wcms.1640](https://doi.org/10.1002/wcms.1640)
- Schaller et al. (2020) "Next generation 3D pharmacophore modeling" [DOI:10.1002/wcms.1468](https://doi.org/10.1002/wcms.1468)

---

## 9. Implementation Checklist

- [ ] Create `pharmacophore/rdshape_optimizer.py` module
- [ ] Add configurable AlignMol parameters to benchmark.py
- [ ] Implement parameter grid search validation
- [ ] Test enhanced fragment representations
- [ ] Implement multi-reference ensemble scoring
- [ ] Run comprehensive parameter sweep
- [ ] Document optimal parameters
- [ ] Update auto_optimizer.py with best parameters
- [ ] Create tutorial notebook demonstrating usage
- [ ] Write unit tests for new functionality

---

## Appendix A: Fragment SMARTS Verification

To verify our fragments are recognized by rdShapeAlign's color feature detection:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Test that our fragments match expected patterns
test_cases = {
    '[NH3]': {'donor': True, 'acceptor': False},
    '[O]': {'donor': False, 'acceptor': True},
    'c1ccccc1': {'aromatic': True, 'hydrophobe': True},
    '[CH4]': {'hydrophobe': True},
}

for smiles, expected in test_cases.items():
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    # Get pharmacophore features
    factory = AllChem.BuildFeatureFactory()
    features = factory.GetFeaturesForMol(mol)
    print(f"{smiles}: {[f.GetFamily() for f in features]}")
```

---

## Appendix B: Performance Benchmarks

Expected timing for parameter sweep:

| Configuration | Mols/sec | Time for 575 mols |
|--------------|----------|-------------------|
| max_preiters=10 | ~250 | ~2.3 sec |
| max_preiters=50 | ~100 | ~5.8 sec |
| max_preiters=100 | ~50 | ~11.5 sec |

With 25 parameter combinations and 575 molecules:
- Fast sweep: ~1 minute
- Thorough sweep: ~5 minutes

# Implementation Prompt: rdShapeAlign Consensus Pharmacophore Optimization

> **Copy this entire prompt to start a new Claude Code session for implementation**

---

## Context

I'm working on the `pharmacophore-toolkit` project to optimize 3D consensus pharmacophore models for virtual screening. The goal is to use RDKit's `rdShapeAlign.AlignMol` function to create consensus pharmacophore models from reference ligands that can discriminate between active compounds and decoys.

### Current State

- **Problem**: Current 3D shape-based scoring achieves only AUC 0.47-0.65, worse than 2D fingerprints (AUC 0.93)
- **Hypothesis**: The rdShapeAlign function SHOULD work with disconnected pharmacophore fragments because:
  1. Gaussian overlap doesn't require molecular connectivity
  2. Color features use SMARTS patterns that match our fragments (NH3→Donor, O→Acceptor, benzene→Aromatic)
- **The issue is likely**: Suboptimal parameters (`opt_param`, `max_preiters`, `max_postiters`) and consensus generation settings

### Key Files to Reference

```
pharmacophore-toolkit/
├── pharmacophore/
│   ├── mol_converter.py      # Converts features to mol fragments (lines 70-78 for SMILES)
│   ├── consensus.py          # Generates consensus pharmacophore
│   ├── optimization.py       # Grid search optimizer (lines 106-121 for scoring)
│   ├── benchmark.py          # Benchmarking framework (lines 181-197 for AlignMol usage)
│   ├── auto_optimizer.py     # Bayesian optimization
│   └── doe_optimization.py   # DOE-based optimization
├── tutorials/data/
│   ├── CCR2_reference_ligands.sdf  # 5 reference ligands
│   ├── actives_ccr2_N75.csv        # 75 actives
│   └── decoys_ccr2_N500.csv        # 500 decoys
└── docs/
    └── RDSHAPEALIGN_OPTIMIZATION_PLAN.md  # Full plan document
```

### rdShapeAlign API (Critical Parameters)

```python
from rdkit.Chem.rdShapeAlign import AlignMol, PrepareConformer

shape_score, color_score = AlignMol(
    ref=reference_mol,
    probe=query_mol,
    useColors=True,              # MUST be True for pharmacophore features
    opt_param=0.5,               # 0=color only, 0.5=balanced, 1.0=shape only (DEFAULT=1.0!)
    max_preiters=10,             # Phase 1 iterations (DEFAULT=10, try 50)
    max_postiters=30             # Phase 2 iterations (DEFAULT=30, try 100)
)
combo = shape_score + color_score  # Range 0-2
```

**Critical**: Default `opt_param=1.0` means SHAPE ONLY! We need `opt_param=0.5` for balanced shape+color.

---

## Task

Implement a systematic parameter optimization to find the best rdShapeAlign settings for consensus pharmacophore discrimination.

### Step 1: Create Parameter Sweep Script

Create `experiments/rdshape_param_sweep.py` that:

1. Loads reference ligands, actives, and decoys from `tutorials/data/`
2. Tests parameter combinations:
   ```python
   PARAM_GRID = {
       'opt_param': [0.0, 0.25, 0.5, 0.75, 1.0],
       'max_preiters': [10, 20, 50],
       'max_postiters': [30, 50, 100],
       'tolerance': [1.5, 2.0, 2.5],
       'occurrence_threshold': [0.4, 0.5, 0.6],
   }
   ```
3. For each combination:
   - Generate consensus pharmacophore with tolerance/occurrence
   - Convert to mol using `PharmacophoreToMol.convert(enable_color_features=True)`
   - Score all actives and decoys with AlignMol using current params
   - Calculate ROC-AUC, BEDROC, EF@1%
4. Output results to CSV and identify best parameters

### Step 2: Verify Fragment Color Recognition

Create a test to verify our fragments are recognized by rdShapeAlign's color feature detection:

```python
def test_fragment_color_features():
    """Verify NH3, O, benzene, CH4 are recognized as pharmacophore features."""
    from rdkit.Chem.rdShapeAlign import AlignMol

    # Create test molecules
    donor_frag = Chem.MolFromSmiles('[NH3]')
    acceptor_frag = Chem.MolFromSmiles('[O]')
    aromatic_frag = Chem.MolFromSmiles('c1ccccc1')

    # Align identical fragments - color score should be 1.0
    # If color score is 0.0, fragments aren't recognized
```

### Step 3: Implement Optimized Scorer

Create `pharmacophore/rdshape_optimizer.py` with class `OptimizedRDShapeScorer`:

1. Configurable AlignMol parameters
2. PrepareConformer caching for speed (10% faster)
3. Multi-conformer scoring with best selection
4. Methods: `score_against_consensus()`, `score_against_references()`

### Step 4: Run Validation

Run the parameter sweep and report:
- Best parameter combination
- AUC improvement vs baseline
- Score distributions for actives vs decoys
- Whether we can beat AUC 0.80 (target) or 0.90 (stretch)

---

## Success Criteria

| Metric | Baseline | Target | Stretch |
|--------|----------|--------|---------|
| ROC-AUC | 0.47-0.65 | >0.80 | >0.90 |
| BEDROC | 0.15 | >0.40 | >0.60 |
| EF@1% | 2.0 | >10.0 | >20.0 |

---

## Important Notes

1. **Read the plan document first**: `docs/RDSHAPEALIGN_OPTIMIZATION_PLAN.md`
2. **The fragments SHOULD work**: Don't change mol_converter.py fragment definitions yet
3. **Focus on parameter tuning first**: The hypothesis is that defaults are suboptimal
4. **Use CCR2 dataset for testing**: It's our standard benchmark
5. **Track all results**: Save to CSV for comparison

---

## Start Here

1. Read `docs/RDSHAPEALIGN_OPTIMIZATION_PLAN.md` for full context
2. Read `pharmacophore/benchmark.py` lines 157-197 to see current AlignMol usage
3. Create the parameter sweep script
4. Run experiments and report findings

Let's optimize rdShapeAlign to achieve the best possible 3D pharmacophore discrimination!

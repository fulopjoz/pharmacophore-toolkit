# Pharmacophore Scoring Approaches

This document describes the four pharmacophore-based virtual screening approaches
implemented in pharmacophore-toolkit, tested on the CCR2 chemokine receptor dataset
(74 actives, 500 property-matched decoys, 5 reference ligands).

---

## Shared Infrastructure: The Evaluation Pipeline

All three 3D approaches share the same evaluation pipeline (`pharmacophore/evaluation.py`),
ensuring fair comparison:

```
Reference ligands (5 CCR2 ligands with 3D coordinates)
       |
       v
PharmacophoreConsensus(tolerance, occurrence, linkage)
  -> Agglomerative hierarchical clustering of pharmacophore features
  -> Deterministic consensus: features shared by >= occurrence_threshold fraction
       |
       v
PharmacophoreToMol.convert(features, enable_color_features=True)
  -> Each feature becomes a molecular fragment at its 3D position:
       Donor -> NH3, Acceptor -> C=O, Aromatic -> benzene, Hydrophobe -> cyclopentane
  -> Fragments placed at exact consensus coordinates
  -> Color features encode pharmacophore type for AlignMol matching
       |
       v
For each query molecule (active or decoy):
  -> Generate N conformers (ETKDGv3, seed=42, no UFF minimization)
  -> For each conformer:
       AlignMol(ref=pharmacophore_mol, probe=query, useColors=True, opt_param=p)
       -> Returns (shape_tanimoto, color_tanimoto), both in [0, 1]
       -> Combined: score = shape_weight * shape + (1 - shape_weight) * color
  -> Keep best score across all conformers
       |
       v
Collect all scores -> ROC-AUC, BEDROC(alpha=20), EF@1%, EF@5%
```

### Key Concepts

**Shape Tanimoto**: Measures 3D volume overlap between query conformer and
pharmacophore model. Independent of chemical features - pure geometric matching.

**Color Tanimoto**: Measures overlap of pharmacophore features (donors, acceptors,
aromatics, hydrophobes) at matching spatial positions. Encodes chemical function.

**Combo Score**: Shape + Color combined. The `opt_param` controls alignment priority
(0.0 = optimize color overlap, 1.0 = optimize shape overlap, 0.5 = balanced).

---

## Approach 1: Optuna GP-Sampler (3D rdShapeAlign)

**File**: `pharmacophore/optuna_optimizer.py` (411 lines)
**Algorithm**: Bayesian Optimization with Gaussian Process surrogate

### How It Works

A Gaussian Process (GP) builds a probabilistic model of the objective function
(ROC-AUC) as a function of 6 hyperparameters. At each iteration:

1. The GP predicts the expected AUC at untested parameter combinations
2. An acquisition function (Expected Improvement) balances exploration vs exploitation
3. The next point to evaluate is chosen to maximize expected improvement
4. After evaluation, the GP is updated with the new observation

This is **sample-efficient**: it finds good parameters with few evaluations because
the GP learns which regions of parameter space produce high AUC.

### Parameters Optimized (6)

| Parameter | Range | What It Controls |
|-----------|-------|-----------------|
| `tolerance` | 0.5-4.0 A | Spatial clustering threshold for consensus features |
| `occurrence` | 0.1-1.0 | Minimum frequency for a feature to be included |
| `shape_weight` | 0.3-0.9 | Balance between shape overlap and pharmacophore color |
| `opt_param` | 0.0-1.0 | AlignMol alignment target (shape vs color priority) |
| `linkage` | avg/complete/single/ward | Agglomerative clustering linkage method |
| `n_conformers` | 5-20 | Conformers generated per query molecule |

### Why This Approach

- **Best for small budgets** (~10-50 trials): GP surrogate learns efficiently
- **Multi-objective**: Simultaneously optimizes ROC-AUC and BEDROC
- Provides Pareto front of non-dominated solutions
- Parameter importance analysis via Optuna's built-in tools

### Strengths & Limitations

- (+) Sample-efficient: good results with few evaluations
- (+) Principled uncertainty quantification
- (-) GP scales O(n^3) with number of observations
- (-) Assumes smooth objective function (may miss sharp optima)

---

## Approach 2: Optuna NSGA-II (3D rdShapeAlign)

**File**: `pharmacophore/optuna_optimizer.py` (same module, different sampler)
**Algorithm**: Non-dominated Sorting Genetic Algorithm II

### How It Works

NSGA-II is an evolutionary algorithm that maintains a population of parameter
configurations, evolving them through selection, crossover, and mutation:

1. Initialize a random population of parameter sets
2. Evaluate each set's (AUC, BEDROC) using the shared pipeline
3. Rank by non-dominated sorting (Pareto fronts)
4. Select parents via tournament selection with crowding distance
5. Create offspring through crossover and polynomial mutation
6. Combine parents + offspring, select next generation
7. Repeat until budget exhausted

### Parameters Optimized

Same 6 parameters as GP (they share the same search space and evaluation pipeline).

### Why This Approach

- **Best for larger budgets** (~50-200 trials): explores broadly
- **True multi-objective**: naturally produces diverse Pareto front
- No surrogate model assumption - works on any objective landscape
- Crowding distance maintains solution diversity

### Strengths & Limitations

- (+) Handles non-convex, discontinuous objectives well
- (+) Naturally parallelizable (population-based)
- (+) Produces diverse set of trade-off solutions
- (-) Less sample-efficient than GP (needs more evaluations)
- (-) No convergence guarantee in finite time

---

## Approach 3: HypoGen 3-Phase (3D rdShapeAlign)

**File**: `pharmacophore/hypogen_optimizer.py` (562 lines)
**Algorithm**: Constructive-Subtractive-Optimization (inspired by Catalyst/Discovery Studio)

### How It Works

Unlike the generic optimizers above, HypoGen uses pharmacophore-specific domain
knowledge through three distinct phases:

**Phase 1: Constructive (Feature Enumeration)**

Start from the full consensus pharmacophore (all features meeting the threshold).
Systematically enumerate subsets of 3 to 7 features. Each subset is a
"pharmacophore hypothesis" - a potential model.

```
Full consensus: 12 features (5 Acceptor, 4 Aromatic, 3 Hydrophobe)
  -> Generate all C(12,3) + C(12,4) + C(12,5) + ... subsets
  -> Cap at 500 hypotheses (sorted by feature diversity)
  -> Evaluate each hypothesis using shared pipeline
```

**Phase 2: Subtractive (Filtering)**

Filter hypotheses that fail minimum quality thresholds:
- ROC-AUC >= 0.55 (above random)
- BEDROC >= 0.05 (some early enrichment)

This eliminates the vast majority of hypotheses quickly.

**Phase 3: Optimization (Simulated Annealing)**

Take the top-N survivors and refine their feature coordinates using
Simulated Annealing:
- Perturb feature positions with Gaussian noise (sigma=0.5 A)
- Accept improvements always; accept worse solutions with probability
  exp(-delta_cost / temperature)
- Cooling schedule: T = T * 0.95 each iteration
- Cost = 0.5*(1-AUC) + 0.3*(1-BEDROC) + 0.2*(n_features/max_features)

The complexity penalty (0.2 weight) favors simpler models with fewer features,
following Occam's razor - a simpler pharmacophore that discriminates well is
preferable to a complex one.

### Why This Approach

- **Domain-specific**: Designed specifically for pharmacophore optimization
- **Interpretable**: Each phase has clear scientific rationale
- **Feature selection**: Automatically finds the minimal discriminating feature set
- Complexity penalty prevents overfitting to training set

### Strengths & Limitations

- (+) Pharmacophore-aware: feature subsets, not arbitrary parameters
- (+) Interpretable: "this model uses 4 features: 2 Acceptor + 1 Aromatic + 1 Hydrophobe"
- (+) Built-in Occam's razor (complexity penalty)
- (-) Fixed consensus as starting point (depends on tolerance/occurrence choice)
- (-) SA refinement is stochastic (results vary with random seed)
- (-) More evaluations than Bayesian approaches for same budget

---

## Approach 4: Pharm2D Fingerprints (2D)

**File**: `pharmacophore/pharm2d_scoring.py` (298 lines)
**Algorithm**: 2D pharmacophore fingerprint similarity (Gobbi_Pharm2D)

### How It Works

This approach encodes pharmacophore features as a 2D fingerprint, completely
avoiding 3D alignment:

1. **Feature Detection**: RDKit's Gobbi_Pharm2D factory identifies 6 feature types
   in each molecule:
   - Donor (H-bond donor groups)
   - Acceptor (H-bond acceptor groups)
   - NegIonizable (acidic groups)
   - PosIonizable (basic groups)
   - Aromatic (aromatic rings)
   - Hydrophobe (hydrophobic regions)

2. **Pairwise Encoding**: For every pair of detected features, compute their
   topological distance and bin it:
   - Distance bins: 2-3, 3-4, 4-5, 5-6, 6-7, 7-8, 8+ Angstroms
   - Each (feature_type_A, feature_type_B, distance_bin) triple sets a bit

3. **Fingerprint Generation**: The result is a binary fingerprint (~40,000 bits)
   encoding all pairwise pharmacophore relationships

4. **Scoring**: Tanimoto similarity between query fingerprint and each reference
   fingerprint. Final score = maximum similarity to any reference ligand.

```
Reference 1: [1,0,1,1,0,0,1,...] (pharmacophore pair fingerprint)
Reference 5: [0,1,1,0,1,0,1,...]
                    |
Query:       [1,1,1,0,0,1,1,...] -> Tanimoto(query, ref_i) for each ref
                    |
Score = max(Tanimoto_1, Tanimoto_2, ..., Tanimoto_5)
```

### Why This Approach

- **No 3D alignment needed**: Works directly from 2D molecular graph
- **Fast**: Simple bit vector comparison, milliseconds per molecule
- **Information-rich**: ~40,000 bits encode many pharmacophore relationships vs
  ~10-15 features in a 3D consensus model
- **Robust**: Invariant to conformational changes (uses topological distances)

### Why It Outperforms 3D

The 2D approach achieves higher AUC because:

1. **More information**: 40,000 bits vs 10-15 spatial features
2. **No alignment errors**: 3D alignment can fail or find suboptimal overlays
3. **No conformer dependency**: 2D uses topological distances, avoiding conformer
   sampling bias
4. **CCR2 actives share scaffolds**: The actives contain common chemical motifs
   (indoles, benzimidazoles) that 2D fingerprints capture effectively

### Strengths & Limitations

- (+) Fastest approach by orders of magnitude
- (+) No hyperparameter optimization needed
- (+) Highest AUC on CCR2 dataset
- (-) Cannot find scaffold hops (structurally different molecules with same 3D shape)
- (-) No spatial interpretation (which features matter and where)
- (-) Sensitive to reference set chemical diversity

---

## Summary

| Approach | Type | Optimization | What It Finds | Speed |
|----------|------|-------------|---------------|-------|
| **Optuna GP** | 3D | Bayesian (GP surrogate) | Optimal consensus parameters | Medium |
| **Optuna NSGA-II** | 3D | Evolutionary (genetic) | Pareto front of trade-offs | Medium |
| **HypoGen 3-Phase** | 3D | Domain-specific (construct-subtract-refine) | Minimal feature set | Slow |
| **Pharm2D** | 2D | None (direct scoring) | Most similar reference | Fast |

### When to Use Each

- **Pharm2D**: Quick screening, large libraries, scaffold-similar targets
- **Optuna GP**: Small optimization budget, need parameter recommendations fast
- **Optuna NSGA-II**: Larger budget, want diverse solutions on Pareto front
- **HypoGen**: Need interpretable model, want feature selection, publication-ready

---

## Full-Dataset Comparison Results

**Dataset**: CCR2, 74 actives + 499 decoys (573 total), 5 conformers/molecule, ETKDGv3, no UFF
**Date**: 2026-01-30

| Approach | Type | ROC-AUC | BEDROC | Evals | Time (s) |
|----------|------|---------|--------|-------|----------|
| **Pharm2D** | **2D** | **0.9111** | **0.7564** | 1 | 2.9 |
| **Optuna NSGA-II** | **3D** | **0.8137** | **0.6902** | 10 | 259.5 |
| Optuna GP | 3D | 0.7426 | 0.1892 | 10 | 233.0 |
| HypoGen 3-Phase | 3D | 0.7374 | 0.4220 | 60 | 16.1 |

### Key Findings

1. **Pharm2D (2D) dominates** with AUC 0.9111 and BEDROC 0.7564, scoring
   574 molecules in 2.9 seconds. It also achieves EF@1% = 7.8 and EF@5% = 6.9.

2. **NSGA-II is the best 3D approach** (AUC 0.8137, BEDROC 0.6902) with strong
   early enrichment. The evolutionary algorithm explores the parameter space
   more effectively than GP for this problem.

3. **GP underperforms on BEDROC** (0.1892) despite reasonable AUC (0.7426).
   With only 10 trials, the GP surrogate may not have sufficient data to model
   the objective landscape accurately.

4. **HypoGen is fastest** among 3D approaches (16.1s) but achieves moderate AUC
   (0.7374). Its strength is interpretability: it identifies the minimal
   discriminating feature subset.

5. **Why keep 3D despite 2D being better?** The 3D approaches provide spatial
   pharmacophore models that enable scaffold hopping (finding structurally novel
   molecules with the same 3D binding geometry). The 2D approach only finds
   molecules chemically similar to known references.

### Why Pharm2D Wins on This Dataset

The CCR2 actives share common chemical scaffolds (indoles, benzimidazoles, piperidines).
The Pharm2D fingerprint (~40,000 bits) captures these scaffold similarities with high
fidelity. In contrast, a 3D consensus pharmacophore model has only 10-15 spatial
features and depends on conformer quality and alignment accuracy.

For datasets where actives are structurally diverse (scaffold-hopping scenarios),
3D approaches would be expected to close the gap or outperform 2D.

---

## Complete Cross-Notebook Analysis (2026-02-03)

This section consolidates results from **all three tutorial notebooks** and the ablation study,
identifies the root cause of anti-discrimination in 3D methods, and provides actionable
recommendations for choosing and optimizing the best consensus model.

### Master Results Table — All Approaches

Data sources: `ccr2_screening_benchmark.ipynb`, `optimal_model_discovery.ipynb`,
`new_approaches_tutorial.ipynb`, `ablation_study.py`, `CONSENSUS_SCORING_SOLUTION.md`

| # | Approach | Type | ROC-AUC | BEDROC | EF@1% | Features | Source |
|---|----------|------|---------|--------|-------|----------|--------|
| 1 | **Reference Alignment** | 3D-direct | **0.897** | **0.850** | **7.8** | 5 refs | benchmark |
| 2 | **Pharm2D Tanimoto** | 2D | **0.858** | **0.712** | **7.8** | ~40K bits | benchmark |
| 3 | Multi-Fidelity BO (10 evals) | 3D-consensus | 0.853 | 0.512 | 0.0 | 52 | new_approaches |
| 4 | Single-Fidelity BO (10 evals) | 3D-consensus | 0.845 | 0.593 | 0.0 | 50 | new_approaches |
| 5 | Pharm3D Distance Scoring | 3D-consensus | 0.837 | 0.492 | 1.4 | 4-8 | CONSENSUS_SCORING |
| 6 | Shape-based (tol=2.0, occ=0.3) | 3D-consensus | 0.830 | — | — | 13 | CONSENSUS_SCORING |
| 7 | Strategy Selector (agglom, tol=1.5, occ=0.3) | 3D-consensus | 0.813 | 0.964 | — | 13 | new_approaches |
| 8 | Ensemble Consensus (tol=2.5, occ=0.5) | 3D-consensus | 0.762 | 0.815 | — | 11 | new_approaches |
| 9 | Agglomerative (tol=2.5, occ=0.5) | 3D-consensus | 0.613 | 0.535 | — | 8 | new_approaches |
| 10 | Combo (shape+color) default | 3D-consensus | **0.328** | 0.004 | 0.0 | 7 | benchmark |
| 11 | Shape Only default | 3D-consensus | **0.338** | 0.023 | 0.0 | 7 | benchmark |
| 12 | Color Weighted (0.3/0.7) default | 3D-consensus | **0.386** | 0.018 | 0.0 | 7 | benchmark |
| 13 | Agglomerative (tol=1.5, occ=0.7) | 3D-consensus | **0.205** | 0.007 | — | 5 | new_approaches |

Rows 10-13 are **anti-discriminating** (AUC < 0.5) — they rank decoys above actives.

---

### Root Cause Analysis: Why 3 Methods Score Below Random

**Affected methods**: Combo (0.328), Shape Only (0.338), Color Weighted (0.386)

**Mechanism**: These methods align query molecules against a synthetic `PharmacophoreToMol`
molecule using `AlignMol`. The synthetic mol is an assembly of disconnected molecular fragments
(NH3, benzene, CH4, C=O) placed at consensus feature coordinates.

**Why decoys score higher than actives**:

1. **Disconnected fragments have no distance constraints.** AlignMol can position each fragment
   independently to maximize overlap. The shape Tanimoto measures volume overlap with this
   physically meaningless assembly.

2. **Decoys are larger and more shape-diverse.** Property-matched decoys have similar MW/LogP
   to actives but different topology. Their larger, more varied shapes overlap MORE with
   scattered fragment clouds.

3. **Only 7 features → 66 atoms.** With few features, the color signal (pharmacophore matching)
   is weak relative to the shape noise. The shape component dominates and anti-correlates with
   activity.

**Evidence from score distributions** (see `tutorials/optimal_model_distributions.png`):

```
                     Actives     Decoys      Separation
Combo (shape+color): mean=0.33   mean=0.35   -0.02  (INVERTED)
Shape Only:          mean=0.50   mean=0.55   -0.05  (INVERTED)
Color Weighted:      mean=0.27   mean=0.30   -0.03  (INVERTED)
Ref Alignment:       mean=0.52   mean=0.40   +0.12  (CORRECT)
Pharm2D:             mean=0.33   mean=0.22   +0.11  (CORRECT)
```

**Confirming evidence from BO optimization** (new_approaches notebook):

The BO optimizer independently discovers the fix — its best params always converge to:
- `opt_param ≈ 0.03` (near-zero = color-only AlignMol optimization)
- `occurrence ≈ 0.10` (very low = many features → stronger color signal)
- `tolerance ≈ 0.5-0.8` (tight = more feature clusters)

**Confirming evidence from strategy selector** (new_approaches notebook):

| Params | AUC | Features | Verdict |
|--------|-----|----------|---------|
| tol=1.5, occ=0.3 | 0.813 | 13 | MORE features → above random |
| tol=2.0, occ=0.3 | 0.789 | 12 | |
| tol=2.5, occ=0.3 | 0.790 | 11 | |
| tol=2.5, occ=0.5 | 0.613 | 8 | Borderline |
| tol=2.0, occ=0.5 | 0.393 | 7 | FEWER features → below random |
| tol=1.5, occ=0.5 | 0.260 | 7 | |
| tol=1.5, occ=0.7 | 0.205 | 5 | |

Linear correlation: **feature count is the primary driver of AUC** for shape-based scoring.
More features → more color atoms → color signal overwhelms shape noise.

---

### Why the Working Methods Work

**Reference Alignment (AUC=0.897)**: Aligns query molecules directly against REAL reference
ligand conformers. Shape comparison between real molecules is physically meaningful. No
synthetic pharmacophore mol involved.

**Pharm2D Tanimoto (AUC=0.858)**: Encodes pairwise pharmacophore feature distances as a 2D
fingerprint (~40,000 bits). No 3D alignment or shape comparison. Pure topological matching.

**BO-optimized consensus (AUC=0.85)**: Finds parameter configurations that produce 50+ features
with `opt_param→0` (color-only alignment). With many features, the pharmacophore mol has
hundreds of color atoms, so color matching dominates over shape noise.

**Strategy selector best (AUC=0.813, BEDROC=0.964)**: Uses `occurrence=0.3` to include all
features present in at least 2 of 5 references → 13 features. More features + balanced
shape/color scoring still discriminates because color signal is strong enough.

**Ensemble consensus (AUC=0.762)**: Stability-aware feature selection produces more features
(8-19 depending on params) than standard consensus at the same thresholds. However, no
single feature achieved stability > 0.7, indicating the CCR2 reference set has high
parameter sensitivity.

---

### Recommendations: Choosing the Right Approach

#### For Production Virtual Screening

**Use Reference Alignment** (AUC=0.897) when:
- High-quality aligned reference ligands are available
- Maximum discrimination is the priority
- Speed is acceptable (~200 ms/mol with 10 conformers)

**Use Pharm2D Tanimoto** (AUC=0.858) when:
- Speed is critical (~15 ms/mol, 13x faster)
- No 3D conformer generation needed
- Scaffold similarity is expected among actives

#### For Consensus Pharmacophore Screening

If a 3D consensus pharmacophore model is required (e.g., for interpretability, scaffold
hopping, or PyMOL visualization), use these settings to avoid anti-discrimination:

```python
# SAFE defaults for consensus pharmacophore scoring
consensus = PharmacophoreConsensus(
    tolerance=1.5,                # Tight: more feature clusters
    occurrence_threshold=0.3,     # Low: include features from 2+ of 5 refs
    linkage='average'
)
features = consensus.generate_consensus(ref_mols)
# Expect 11-13 features (NOT 4-7)

# Score with color-dominant weighting
pharm_mol = PharmacophoreToMol.convert(features, enable_color_features=True)
for conf_id in range(query.GetNumConformers()):
    shape, color = AlignMol(
        ref=pharm_mol, probe=query, probeConfId=conf_id,
        useColors=True,
        opt_param=0.0    # CRITICAL: optimize for color, NOT shape
    )
    score = 0.3 * shape + 0.7 * color  # weight toward color
```

**Key parameter rules**:
1. `opt_param=0.0` (color-only alignment optimization) — never use 0.5 or 1.0
2. `occurrence_threshold ≤ 0.3` — ensures >= 10 features
3. `shape_weight ≤ 0.3` — minimizes shape noise contribution

#### For Optimization

**Use Multi-Fidelity BO** (AUC=0.853) when:
- Exploring unknown parameter space
- Budget for 10+ evaluations
- Want automatic parameter discovery

**Use Strategy Selector** (AUC=0.813, BEDROC=0.964) when:
- Want best BEDROC (early enrichment) without optimization
- Need a quick tournament across standard parameter combinations
- BEDROC more important than AUC

---

### Notebook Bug Fix

`new_approaches_tutorial.ipynb` cell `summary_table` had `UnifiedEvaluator(refs, actives, decoys, seed=42)`.
The parameter is `random_state`, not `seed`. Fixed to `random_state=42`.

---

### Open Questions for Further Investigation

1. **MMFF94s energy minimization effect**: Does post-ETKDGv3 minimization improve Reference
   Alignment or BO-optimized scoring? Test with `--minimize` flag in ablation study.

2. **Pharm3D + Pharm2D hybrid**: The `CONSENSUS_SCORING_SOLUTION.md` reports a hybrid
   (0.7*pharm2d + 0.3*pharm3d) achieving AUC=0.95. This should be benchmarked in the
   ablation study alongside shape-based methods.

3. **More reference ligands**: All analysis uses 5 references. Literature recommends 10-20.
   Adding diverse CCR2 actives to the reference set could improve all 3D methods.

4. **Ensemble consensus stability**: No CCR2 feature achieved stability > 0.7 across
   perturbation rounds. This suggests either (a) the reference set is too small, or
   (b) CCR2 ligands bind via multiple modes that don't share a stable pharmacophore.

---

## Dataset

**Target**: CCR2 (C-C chemokine receptor type 2)
**Actives**: 74 known CCR2 modulators
**Decoys**: 500 property-matched decoys (matched physicochemical properties, different topology)
**References**: 5 CCR2 ligands with crystal/docked 3D coordinates
**Source**: DUD-E style benchmark

## Evaluation Metrics

- **ROC-AUC**: Area under ROC curve (0.5 = random, 1.0 = perfect)
- **BEDROC (alpha=20)**: Boltzmann-enhanced discrimination, emphasizes early recognition
- **EF@1%**: Enrichment factor at 1% of ranked database
- **EF@5%**: Enrichment factor at 5% of ranked database

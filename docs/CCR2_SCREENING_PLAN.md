# CCR2 Pharmacophore-Based Virtual Screening Tool
## Comprehensive Development & Research Plan

**Project**: Optimal Pharmacophore Model Discovery for Active/Decoy Separation
**Target**: CCR2 Chemokine Receptor (Allosteric Site)
**Date Created**: 2026-01-23
**Status**: Planning & Research Phase

---

# Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Project Goals & Success Criteria](#2-project-goals--success-criteria)
3. [Available Data](#3-available-data)
4. [Tools, Plugins & MCP Servers](#4-tools-plugins--mcp-servers)
5. [Optimization Strategies (Brainstorming)](#5-optimization-strategies-brainstorming)
6. [Implementation Architecture](#6-implementation-architecture)
7. [Literature Tracking](#7-literature-tracking)
8. [Research Tasks & Questions](#8-research-tasks--questions)
9. [Implementation Phases](#9-implementation-phases)
10. [References](#10-references)

---

# 1. Executive Summary

Build a comprehensive virtual screening tool that:
- Takes **reference ligands** (CCR2 allosteric modulators) to generate consensus pharmacophore models
- Uses **shape matching** (CDPKit) combined with **pharmacophore features** for scoring
- **Optimizes** model parameters to maximize separation between actives (75) and decoys (500)
- Employs **multi-objective optimization** (Optuna) to balance ROC-AUC and Enrichment Factor
- Provides **cross-validated**, statistically robust performance estimates

**Key Decisions Made**:
- Shape backend: **CDPKit** (fast, TanimotoCombo, open-source ROCS alternative)
- Conformer strategy: **Multi-conformer** (10-50 per molecule)
- Optimization: **Pareto multi-objective** with both AUC and EF metrics
- Literature tools: **Zotero MCP**, **Semantic Scholar MCP**, **Browser Use MCP**

---

# 2. Project Goals & Success Criteria

## Primary Goals

| Goal | Description | Success Metric |
|------|-------------|----------------|
| **G1** | Generate optimal consensus pharmacophore | Model captures key CCR2 interaction features |
| **G2** | Maximize active/decoy separation | ROC-AUC > 0.80, EF1% > 10 |
| **G3** | Fast parameter optimization | < 1 hour for full optimization run |
| **G4** | Reproducible & validated | CV with confidence intervals, deterministic |
| **G5** | Extensible tool | Reusable for other targets |

## Success Criteria

```
Minimum Viable:
- ROC-AUC >= 0.75 (5-fold CV)
- EF1% >= 5.0
- Processing: < 5 min per model evaluation

Target Performance:
- ROC-AUC >= 0.85 (repeated 5×10 CV)
- EF1% >= 15.0
- BEDROC(α=20) >= 0.5
- Processing: < 30 sec per model evaluation

Stretch Goals:
- ROC-AUC >= 0.90
- EF1% >= 25.0
- Pareto front with 5+ non-dominated models
- Interactive parameter exploration in Jupyter
```

---

# 3. Available Data

## Data Files Location: `tutorials/data/`

| File | Description | Count | Format |
|------|-------------|-------|--------|
| `CCR2_reference_ligands.sdf` | 3D reference ligands (allosteric) | Multiple | SDF with 3D coords |
| `actives_ccr2_N75.csv` | Active compounds | 75 | CSV with SMILES |
| `decoys_ccr2_N500.csv` | Decoy compounds | 500 | CSV with SMILES |

## Active Compounds Details

```
Columns: Activity_ID, Quality, source, CID, SMILES, ligandtype,
         pchembl_value, pchembl_value_Mean, etc.

Key Features:
- All marked as "allosteric" ligand type
- pchembl values range: ~7.0 - 9.0 (nM-μM potency)
- Sources: ChEMBL30, ExCAPE-DB, Papyrus
- Quality levels: High, Medium
```

## Decoy Compounds Details

```
Columns: Smiles, MW, LogP, NumHDonors, NumHAcceptors

Property Ranges:
- MW: 280-500 Da
- LogP: 2-6
- HBD: 1-3
- HBA: 3-6
```

## Data Quality Notes

- [x] Actives have confirmed activity data (pchembl values)
- [x] Decoys are property-matched (drug-like)
- [ ] Need to verify: Are decoys confirmed inactives or putative?
- [ ] Need to check: Structural diversity of reference set

---

# 4. Tools, Plugins & MCP Servers

## 4.1 Currently Enabled Claude Code Plugins

| Plugin | Purpose | How to Use |
|--------|---------|------------|
| **context7** | Documentation lookup for any library | `mcp__plugin_context7_context7__query-docs` |
| **claude-scientific-skills** | 55+ scientific Python packages | Scientific computing guidance |
| **serena** | Semantic code analysis & editing | `mcp__plugin_serena_serena__*` tools |
| **claude-mem** | Memory/context management | `mcp__plugin_claude-mem_mcp-search__*` |
| **code-simplifier** | Code refinement | Automatic simplification |
| **huggingface-skills** | ML model integration | HuggingFace access |

## 4.2 MCP Servers to Install for Research

### Zotero MCP (Literature Management)

**Purpose**: Connect your Zotero library for literature search and PDF annotation extraction

**Installation**:
```bash
# Install via uv
uv tool install "git+https://github.com/54yyyu/zotero-mcp.git"
zotero-mcp setup

# Or add to ~/.cursor/mcp.json:
{
  "mcpServers": {
    "zotero": {
      "command": "uvx",
      "args": ["zotero-mcp"]
    }
  }
}
```

**Requirements**:
- Zotero 7+
- Enable Local API: Zotero → Preferences → Advanced → Check "Allow other applications..."

**Capabilities**:
- Search papers by query or advanced criteria
- Retrieve metadata, full text, attachments
- Extract PDF annotations and highlights
- Semantic search across library
- Access notes and tags

**Use Cases for This Project**:
- [ ] Search for CCR2 allosteric modulator papers
- [ ] Extract binding site information from crystal structure papers
- [ ] Find pharmacophore validation methodology papers
- [ ] Gather SAR data from medicinal chemistry publications

**Source**: [Zotero MCP GitHub](https://github.com/54yyyu/zotero-mcp)

---

### Semantic Scholar MCP (Paper Search)

**Purpose**: Search 200M+ academic papers, get citations, author profiles

**Installation**:
```bash
claude mcp add semantic-scholar-mcp -s project \
  -e SEMANTIC_SCHOLAR_API_KEY="your-api-key" \
  -- uvx semantic-scholar-mcp

# Or add to config:
{
  "mcpServers": {
    "semantic-scholar": {
      "command": "uvx",
      "args": ["--from", "git+https://github.com/JackKuo666/semanticscholar-MCP-Server", "semantic-scholar-mcp", "serve"],
      "env": {
        "SEMANTIC_SCHOLAR_API_KEY": "your-key"
      }
    }
  }
}
```

**Get API Key**: https://www.semanticscholar.org/product/api (free tier available)

**Capabilities**:
- Search papers with filters (year, venue, fields)
- Get full paper details with abstracts
- Fetch citations and references
- Author profiles with h-index
- AI-powered paper recommendations
- Bulk retrieval (up to 500 papers)

**Use Cases for This Project**:
- [ ] Find recent pharmacophore virtual screening papers (2023-2025)
- [ ] Search for CCR2 drug discovery publications
- [ ] Discover shape-based screening methodology papers
- [ ] Find enrichment factor / ROC-AUC validation studies

**Sources**:
- [Semantic Scholar MCP (smaniches)](https://glama.ai/mcp/servers/@smaniches/semantic-scholar-mcp)
- [JackKuo Semantic Scholar MCP](https://github.com/JackKuo666/semanticscholar-MCP-Server)

---

### Browser Use MCP (Logged-In Sessions)

**Purpose**: Use your browser's logged-in sessions for accessing paywalled journals

**Installation**:
```bash
# Local stdio-based server
uvx browser-use --mcp

# Or hosted version (add to config):
{
  "mcpServers": {
    "browser-use": {
      "url": "https://api.browser-use.com/mcp"
    }
  }
}
```

**Configuration Options**:
```
BROWSER_HEADLESS=false    # Show browser for debugging
BROWSER_TYPE=chromium     # or firefox, webkit
BROWSER_USER_DATA_DIR=~/.browser-use  # Persist sessions
```

**Capabilities**:
- Control browser with your logged-in sessions
- Access institutional journal subscriptions
- Navigate paywalled content
- Extract data from authenticated pages
- Avoid bot detection (uses real browser profile)

**Use Cases for This Project**:
- [ ] Access full-text papers from subscribed journals
- [ ] Download supplementary materials
- [ ] Access ChEMBL/PubChem with institutional features
- [ ] Navigate patent databases

**Source**: [Browser Use MCP Docs](https://docs.browser-use.com/customize/integrations/mcp-server)

---

### Paper Search MCP (Multi-Source)

**Purpose**: Search and download papers from multiple academic platforms

**Installation**:
```bash
{
  "mcpServers": {
    "paper-search": {
      "command": "uvx",
      "args": ["paper-search-mcp"]
    }
  }
}
```

**Supported Sources**:
- arXiv
- PubMed
- bioRxiv / medRxiv
- Google Scholar
- IACR ePrint Archive
- Semantic Scholar

**Use Cases for This Project**:
- [ ] Search arXiv for recent ML + pharmacophore papers
- [ ] Query PubMed for CCR2 drug discovery
- [ ] Find preprints on shape-based virtual screening

**Source**: [paper-search-mcp GitHub](https://github.com/openags/paper-search-mcp)

---

## 4.3 Python Libraries for Implementation

### Core Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| `rdkit` | >= 2023.3.1 | Molecular operations, conformers, alignment |
| `cdpkit` | >= 1.2.3 | Fast shape screening (TanimotoCombo) |
| `scikit-learn` | >= 1.0.0 | Clustering, CV, metrics |
| `optuna` | >= 3.0.0 | Bayesian optimization, Pareto |
| `pandas` | >= 2.0.0 | Data manipulation |
| `numpy` | >= 2.0.0 | Numerical computing |
| `matplotlib` | >= 3.0 | Visualization |
| `joblib` | >= 1.0 | Parallel processing |

### Optional/Recommended

| Library | Purpose | When to Use |
|---------|---------|-------------|
| `pymoo` | NSGA-II multi-objective | Alternative to Optuna |
| `roshambo2` | GPU shape screening | Large-scale screening |
| `py3Dmol` | 3D visualization | Interactive Jupyter |
| `plotly` | Interactive plots | Pareto front visualization |

### Installation Command

```bash
# Core
pip install rdkit cdpkit scikit-learn optuna pandas numpy matplotlib joblib

# Optional
pip install pymoo roshambo2 py3Dmol plotly

# Or with Poetry
poetry add rdkit cdpkit scikit-learn optuna pandas numpy matplotlib joblib
```

---

## 4.4 Skills & Commands Available

### Built-in Claude Code Skills

| Skill | Trigger | Use Case |
|-------|---------|----------|
| `/literature-review` | Literature synthesis | Systematic review of pharmacophore methods |
| `/scientific-reviewer` | Document analysis | Review methodology papers |
| `/brainstorming` | Ideation facilitation | Explore optimization approaches |
| `/design-of-experiments` | DOE guidance | Plan validation experiments |
| `/python-expert` | Python optimization | Performance tuning |
| `/data-scientist` | Data analysis | Metric analysis |
| `/ml-engineer` | ML pipelines | Surrogate model building |

### Recommended Workflow

```
1. Use Semantic Scholar MCP → Find relevant papers
2. Use Zotero MCP → Organize and annotate papers
3. Use /literature-review → Synthesize findings
4. Use /scientific-reviewer → Critically evaluate methods
5. Use /brainstorming → Generate implementation ideas
6. Use context7 → Get library documentation
7. Use serena → Navigate and edit codebase
```

---

# 5. Optimization Strategies (Brainstorming)

## 5.1 Parameter Space Overview

```
                    PHARMACOPHORE OPTIMIZATION
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
   CONSENSUS PARAMS      MODEL SELECTION       SCORING CONFIG
        │                     │                     │
   ├─ tolerance (Å)      ├─ strict            ├─ shape weight
   │   [0.5 - 4.0]       ├─ moderate          ├─ color weight
   │                     ├─ relaxed           ├─ combo mode
   ├─ occurrence_threshold├─ comprehensive    │
   │   [0.3 - 0.9]       │   (all above)      └─ feature subset
   │                     │                        ├─ Donor+Acceptor
   ├─ linkage method     └─ custom thresholds    ├─ All 4 types
   │   ├─ average                                └─ Custom
   │   ├─ complete
   │   ├─ single
   │   └─ ward
   │
   └─ feature_types
       ├─ default (4 types)
       └─ rdkit extended
```

**Total Parameter Space**: ~6 dimensions, potentially 1000s of combinations

---

## 5.2 Algorithm Comparison

| Algorithm | Speed | Sample Efficiency | Multi-Objective | Implementation |
|-----------|-------|-------------------|-----------------|----------------|
| **Grid Search** | Slow O(n^k) | Poor | No (scalarize) | sklearn, manual |
| **Random Search** | Moderate | Better than grid | No | sklearn |
| **Bayesian/TPE (Optuna)** | Fast | Excellent | Yes (built-in) | `optuna` |
| **Genetic/Evolutionary** | Moderate | Good | Yes (NSGA-II) | `pymoo`, `DEAP` |
| **Hyperband** | Very Fast | Good | No | `optuna` |
| **Pareto-MCTS** | Moderate | Excellent | Native | Custom |

---

## 5.3 Six Thinking Hats Analysis

### White Hat (Facts)
- 75 actives, 500 decoys = imbalanced dataset
- Consensus pharmacophore has 4+ parameters
- Full screening: ~575 compounds × 50 conformers = 28,750 evaluations per trial
- Need cross-validation for robust estimates

### Red Hat (Intuition)
- Optuna with Pareto feels right - proven in drug discovery
- Grid search too slow but good for visualization
- CV with 75 actives is tricky - stratified is a must

### Black Hat (Risks)
- **Small active set**: High variance in CV folds
- **Overfitting**: Optimizing metrics on same data
- **Local optima**: Pharmacophore space may be multimodal
- **Speed**: Full screening per evaluation is slow
- **Trade-off**: AUC vs EF1% genuinely conflict

### Yellow Hat (Benefits)
- Optuna handles multi-objective natively
- Pruning: Early stopping saves compute
- Caching: Conformers can be reused
- Parallelization: Independent trials run concurrently
- Transfer learning: Results inform other targets

### Green Hat (Creative Ideas)

**Idea 1: Hierarchical Optimization**
```
Phase 1: Coarse grid on (tolerance, occurrence) → promising regions
Phase 2: Bayesian refinement in promising regions
Phase 3: Model selection on top candidates
```

**Idea 2: Multi-Fidelity Optimization**
```
Low fidelity: Score with 1 conformer (fast)
Medium fidelity: Score with 10 conformers
High fidelity: Full 50 conformers (slow)
Progressive refinement from low → high
```

**Idea 3: Surrogate-Assisted Optimization**
```
Train fast ML model (RF/XGBoost) on parameter → metric mapping
Use surrogate for exploration, real evaluation for exploitation
```

**Idea 4: Ensemble Pareto Selection**
```
Instead of single "best" model:
- Identify Pareto front of models
- Use ensemble of Pareto-optimal pharmacophores
- Vote/average scores for final ranking
```

**Idea 5: Active Learning Loop**
```
1. Start with small random sample
2. Train surrogate on metrics
3. Select next trial by uncertainty + expected improvement
4. Repeat until convergence
```

### Blue Hat (Process Summary)
**Recommended**: Optuna with multi-objective (Pareto) + Hyperband pruning

---

## 5.4 Primary Recommendation: Optuna Multi-Objective

```python
import optuna
from optuna.samplers import TPESampler

def objective(trial):
    # Sample parameters
    tolerance = trial.suggest_float('tolerance', 0.5, 4.0)
    occurrence = trial.suggest_float('occurrence_threshold', 0.3, 0.9)
    model_type = trial.suggest_categorical('model', ['strict', 'moderate', 'relaxed'])
    linkage = trial.suggest_categorical('linkage', ['average', 'complete', 'ward'])

    # Cross-validation
    cv_results = cross_validate_pharmacophore(
        reference_mols, actives, decoys,
        tolerance=tolerance,
        occurrence_threshold=occurrence,
        model_type=model_type,
        linkage=linkage,
        n_splits=5,
        n_repeats=3
    )

    # Return multiple objectives (Optuna handles Pareto)
    return cv_results['mean_auc'], cv_results['mean_ef1']

# Multi-objective study
study = optuna.create_study(
    directions=['maximize', 'maximize'],
    sampler=TPESampler(multivariate=True),
    pruner=optuna.pruners.HyperbandPruner()
)
study.optimize(objective, n_trials=100, n_jobs=-1)

# Get Pareto front
pareto_trials = study.best_trials
```

**Why Optuna?**
- TPE sampler is 10-100x more sample-efficient than grid
- Native Pareto front tracking
- Hyperband pruning eliminates poor trials early (~60% compute savings)
- Trivial parallelization with `n_jobs=-1`
- Built-in visualization

---

## 5.5 Handling AUC vs EF1% Trade-off

### The Fundamental Conflict

| Metric | Measures | Pharmacophore Style |
|--------|----------|---------------------|
| **AUC** | Global discrimination | Relaxed (catches all actives) |
| **EF1%** | Early enrichment | Strict (highly selective) |

### Solution Strategies

**Strategy 1: Pareto Front Selection**
```
Don't choose one model - keep all Pareto-optimal models
Use case determines selection:
- High-throughput screen → prioritize AUC
- Cherry-picking candidates → prioritize EF1%
- Balanced → pick knee of Pareto curve
```

**Strategy 2: BEDROC as Compromise**
BEDROC naturally balances early vs late enrichment:
- α = 20: Heavy early enrichment focus (like EF1%)
- α = 1: More global (like AUC)

```python
def bedroc(y_true, y_scores, alpha=20.0):
    """Balanced early-recognition metric"""
    from rdkit.ML.Scoring import Scoring
    return Scoring.CalcBEDROC(scores, col=0, alpha=alpha)
```

**Strategy 3: Constrained Optimization**
```
Maximize EF1% subject to AUC >= 0.75
(or vice versa)
```

---

## 5.6 Cross-Validation Strategies for Small Active Sets

### The Challenge
- Standard 5-fold CV → only 15 actives per test fold
- High variance in metrics, especially EF1%

### Recommended: Repeated Stratified K-Fold

```python
from sklearn.model_selection import RepeatedStratifiedKFold

cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)
# 50 total evaluations, robust mean estimates
```

### Alternative: Bootstrap with Confidence Intervals

```python
def bootstrap_cv(actives, decoys, n_bootstrap=100):
    results = []
    for _ in range(n_bootstrap):
        boot_actives = np.random.choice(actives, len(actives), replace=True)
        boot_decoys = np.random.choice(decoys, len(decoys), replace=True)
        results.append(evaluate(boot_actives, boot_decoys))
    return np.mean(results), np.std(results), np.percentile(results, [2.5, 97.5])
```

### Always Report Uncertainty!

```python
results = {
    'auc_mean': 0.85,
    'auc_std': 0.03,
    'auc_ci95': (0.79, 0.91),
    'ef1_mean': 12.5,
    'ef1_std': 3.2,
    'ef1_ci95': (6.1, 18.9)
}
```

---

## 5.7 Speed Optimization Strategies

### Bottleneck Analysis

| Step | Time | Optimization |
|------|------|--------------|
| Conformer generation | 1-5s/mol | Pre-generate, cache |
| Alignment | 0.1s/pair | Batch processing |
| Pharmacophore extraction | 0.05s/mol | Vectorized |
| Consensus clustering | 0.01s | Already fast |
| Shape scoring | 0.1-1s/mol | GPU, parallel |

### Implementation

**1. Pre-compute Conformers**
```python
CONFORMER_CACHE = {
    'actives': [generate_conformers(mol, n=50) for mol in actives],
    'decoys': [generate_conformers(mol, n=50) for mol in decoys],
}
```

**2. Multi-Fidelity Evaluation**
```python
def evaluate_with_fidelity(params, fidelity='low'):
    n_conf = {'low': 1, 'medium': 10, 'high': 50}[fidelity]
```

**3. Parallel Screening**
```python
from joblib import Parallel, delayed
scores = Parallel(n_jobs=-1)(
    delayed(score_molecule)(mol, pharmacophore)
    for mol in compounds
)
```

---

## 5.8 Recommended Optimization Stack

```
┌─────────────────────────────────────────────────────────────┐
│                    OPTIMIZATION PIPELINE                     │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐     │
│  │  OPTUNA     │───▶│  MULTI-OBJ  │───▶│   PARETO    │     │
│  │  TPE+Hyper  │    │  AUC + EF1% │    │   FRONT     │     │
│  └─────────────┘    └─────────────┘    └─────────────┘     │
│         │                                     │             │
│         ▼                                     ▼             │
│  ┌─────────────┐                      ┌─────────────┐      │
│  │ MULTI-FIDELITY                     │   MODEL     │      │
│  │ 1→10→50 conf │                     │  PORTFOLIO  │      │
│  └─────────────┘                      └─────────────┘      │
│         │                                                   │
│         ▼                                                   │
│  ┌─────────────────────────────────────────────────────┐   │
│  │          REPEATED STRATIFIED 5×10 CV                 │   │
│  │          with Bootstrap Confidence Intervals         │   │
│  └─────────────────────────────────────────────────────┘   │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

---

# 6. Implementation Architecture

## 6.1 Module Structure

```
pharmacophore/
├── __init__.py
├── pharmacophore.py           # Existing: Main Pharmacophore class
├── consensus.py               # Existing: PharmacophoreConsensus
├── mol_converter.py           # Existing: Feature → RDKit Mol
├── constants.py               # Existing: Feature definitions
├── draw.py                    # Existing: Visualization
│
└── screening/                 # NEW MODULE
    ├── __init__.py
    ├── scorer.py              # Pharmacophore/shape scoring
    ├── metrics.py             # ROC, EF, BEDROC calculations
    ├── optimizer.py           # Optuna parameter optimization
    ├── pipeline.py            # End-to-end screening pipeline
    ├── conformers.py          # Conformer generation & caching
    └── reporters.py           # Results visualization/export
```

## 6.2 Core Classes

### PharmacophoreScreener

```python
class PharmacophoreScreener:
    """Score compounds against pharmacophore model using shape/color matching."""

    def __init__(self,
                 pharmacophore_mol: Chem.Mol,
                 scoring_mode: str = 'tanimoto_combo',
                 backend: str = 'cdpkit'):
        """
        Args:
            pharmacophore_mol: RDKit Mol from PharmacophoreToMol.convert()
            scoring_mode: 'shape', 'color', 'tanimoto_combo'
            backend: 'cdpkit' or 'rdkit'
        """

    def score_molecule(self, mol: Chem.Mol) -> float:
        """Score single molecule against pharmacophore."""

    def screen_library(self,
                       mols: List[Chem.Mol],
                       n_jobs: int = -1) -> pd.DataFrame:
        """Parallel screening of compound library."""
```

### ScreeningMetrics

```python
class ScreeningMetrics:
    """Calculate virtual screening validation metrics."""

    @staticmethod
    def roc_auc(y_true, y_scores) -> float:
        """Area under ROC curve."""

    @staticmethod
    def enrichment_factor(y_true, y_scores,
                          percentages=[0.005, 0.01, 0.02, 0.05]) -> dict:
        """EF at multiple cutoffs."""

    @staticmethod
    def bedroc(y_true, y_scores, alpha=20.0) -> float:
        """Boltzmann-enhanced early recognition."""

    @staticmethod
    def generate_report(y_true, y_scores) -> dict:
        """Comprehensive metrics report with CIs."""
```

### PharmacophoreOptimizer

```python
class PharmacophoreOptimizer:
    """Optimize pharmacophore parameters using Optuna multi-objective."""

    def __init__(self,
                 reference_mols: List[Chem.Mol],
                 actives: List[Chem.Mol],
                 decoys: List[Chem.Mol],
                 cv_strategy: str = 'repeated_stratified'):
        """Initialize with data for optimization."""

    def optimize(self,
                 n_trials: int = 100,
                 objectives: List[str] = ['auc', 'ef1'],
                 n_jobs: int = -1) -> optuna.Study:
        """Run multi-objective optimization."""

    def get_pareto_front(self) -> List[dict]:
        """Return Pareto-optimal parameter sets."""

    def visualize_pareto(self, save_path: str = None):
        """Plot Pareto front."""
```

---

# 7. Literature Tracking

## 7.1 Key Papers to Read

### Pharmacophore Methodology

| Status | Paper | Year | Key Contribution | DOI/Link |
|--------|-------|------|------------------|----------|
| [ ] | PharmacoNet: Deep Pharmacophore Modeling | 2024 | Ultra-fast screening (187M/21h) | [arXiv:2310.00681](https://arxiv.org/abs/2310.00681) |
| [ ] | PharmacoMatch: Neural Subgraph Matching | 2024 | Contrastive learning for pharmacophore | [arXiv:2409.06316](https://arxiv.org/abs/2409.06316) |
| [ ] | dyphAI: Dynamic Pharmacophore with AI | 2025 | Ensemble pharmacophore, Optuna | [Frontiers](https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2025.1479763/full) |
| [ ] | Pharmacophore Validation Review | 2025 | Comprehensive metrics overview | [Springer](https://link.springer.com/article/10.1007/s10822-025-00751-9) |

### Shape-Based Screening

| Status | Paper | Year | Key Contribution | DOI/Link |
|--------|-------|------|------------------|----------|
| [ ] | ROSHAMBO2: GPU Shape Alignment | 2025 | Fast GPU-accelerated screening | [JCIM](https://iwatobipen.wordpress.com/2025/09/28/gpu-based-fast-shape-alignment-of-molecules-rdkit-roshambo2-cheminformatics/) |
| [ ] | CDPKit shapescreen | 2024 | Open-source ROCS alternative | [CDPKit Docs](https://cdpkit.org/applications/shapescreen.html) |
| [ ] | VSFlow: Open-source VS Tool | 2022 | Combo scores, 3D fingerprints | [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/628c60215d9485a206cc8ecc) |

### Optimization & ML

| Status | Paper | Year | Key Contribution | DOI/Link |
|--------|-------|------|------------------|----------|
| [ ] | Pareto Optimization for VS | 2024 | Multi-objective VS acceleration | [RSC Digital Discovery](https://pubs.rsc.org/en/content/articlehtml/2024/dd/d3dd00227f) |
| [ ] | SALSA: Active Learning on Synthons | 2025 | Scalable AL for trillion compounds | [arXiv:2505.12913](https://arxiv.org/abs/2505.12913) |
| [ ] | Bayesian Optimization for Drug Discovery | 2025 | Multi-fidelity optimization | [ACS Central Science](https://pmc.ncbi.nlm.nih.gov/articles/PMC11869128/) |
| [ ] | Multi-Objective Molecular Generation (PMMG) | 2025 | Pareto MCTS | [Advanced Science](https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/advs.202410640) |

### CCR2 Target-Specific

| Status | Paper | Year | Key Contribution | DOI/Link |
|--------|-------|------|------------------|----------|
| [ ] | CCR2 allosteric modulators discovery | ? | Reference ligand SAR | Search needed |
| [ ] | CCR2 crystal structure analysis | ? | Binding site features | Search needed |
| [ ] | CCR2 pharmacophore models | ? | Existing models to compare | Search needed |

---

## 7.2 Literature Search Queries

### For Semantic Scholar MCP

```
Query 1: "pharmacophore virtual screening machine learning optimization"
Filters: year >= 2023, fields: Computer Science, Biology

Query 2: "CCR2 chemokine receptor allosteric modulator"
Filters: year >= 2020

Query 3: "shape-based virtual screening enrichment factor ROC"
Filters: year >= 2022

Query 4: "Bayesian optimization drug discovery"
Filters: year >= 2024

Query 5: "consensus pharmacophore generation clustering"
Filters: year >= 2020
```

### For Zotero Search (Local Library)

```
Tags to search: pharmacophore, virtual screening, CCR2, shape matching
Collections to check: Drug Discovery, Cheminformatics, Methods
```

### For arXiv/bioRxiv

```
Search: "pharmacophore" AND ("deep learning" OR "machine learning")
Category: q-bio.QM, cs.LG

Search: "virtual screening" AND "active learning"
Category: q-bio.BM
```

---

## 7.3 Notes & Annotations

### PharmacoNet Notes
```
Date read: [pending]
Key findings:
-
Relevance to project:
-
Code available: https://github.com/SeonghwanSeo/PharmacoNet
```

### dyphAI Notes
```
Date read: [pending]
Key findings:
-
Relevance to project:
-
Optuna usage details:
-
```

*[Add more paper notes as you read them]*

---

## 7.4 Research Questions to Answer

### Methodology Questions

- [ ] What is the optimal number of conformers for shape screening?
- [ ] How does CDPKit compare to ROCS on standard benchmarks?
- [ ] What α value for BEDROC best correlates with experimental success?
- [ ] How do ensemble pharmacophore methods compare to single models?

### CCR2-Specific Questions

- [ ] What are the key pharmacophore features for CCR2 allosteric binding?
- [ ] Are there existing validated CCR2 pharmacophore models to compare?
- [ ] What is the binding site shape/volume for the allosteric pocket?
- [ ] Which reference ligands have crystal structures available?

### Validation Questions

- [ ] What is a "good" enrichment factor for CCR2-like targets?
- [ ] How many actives are needed for reliable CV estimates?
- [ ] What decoy selection method was used (DUD-E, property-matched)?

---

# 8. Research Tasks & Questions

## 8.1 Immediate Research Tasks

| Priority | Task | Tool to Use | Status |
|----------|------|-------------|--------|
| HIGH | Search for recent pharmacophore optimization papers | Semantic Scholar MCP | [ ] |
| HIGH | Find CCR2 allosteric modulator SAR papers | Zotero MCP / Browser Use | [ ] |
| HIGH | Read PharmacoNet paper for implementation ideas | context7 + paper | [ ] |
| MEDIUM | Compare CDPKit vs RDKit shape scoring | Literature + experiment | [ ] |
| MEDIUM | Find enrichment factor benchmarks for GPCRs | Semantic Scholar | [ ] |
| LOW | Explore surrogate model approaches | arXiv search | [ ] |

## 8.2 Questions for Literature

### To Answer from Papers

1. **What conformer count is standard?**
   - Search: "conformer generation virtual screening number"
   - Expected: 10-100 conformers typical

2. **What are typical EF1% values for good pharmacophore models?**
   - Search: "enrichment factor pharmacophore benchmark"
   - Expected: EF1% > 10 is good, > 20 is excellent

3. **How does multi-objective optimization compare to weighted sum?**
   - Search: "Pareto optimization drug discovery comparison"

4. **What is the best CV strategy for small active sets?**
   - Search: "cross validation small sample size drug discovery"

---

# 9. Implementation Phases

## Phase 1: Setup & Data Preparation
- [ ] Install CDPKit for shape screening
- [ ] Configure MCP servers (Zotero, Semantic Scholar, Browser Use)
- [ ] Create data loading utilities for CSV/SDF files
- [ ] Implement conformer generation pipeline with caching
- [ ] Implement molecule alignment workflow

## Phase 2: Scoring Module
- [ ] Implement `PharmacophoreScreener` class
- [ ] Add CDPKit shape scoring integration
- [ ] Add RDKit shape scoring as fallback
- [ ] Implement parallel processing for large libraries
- [ ] Benchmark scoring speed

## Phase 3: Metrics Module
- [ ] Implement ROC AUC calculation
- [ ] Implement enrichment factor at multiple cutoffs
- [ ] Implement BEDROC metric
- [ ] Create visualization functions (ROC curves, enrichment plots)
- [ ] Generate comprehensive reports with CIs

## Phase 4: Optimization Module
- [ ] Implement Optuna multi-objective optimizer
- [ ] Add Hyperband pruning
- [ ] Implement repeated stratified CV
- [ ] Add Pareto front visualization
- [ ] Store/load optimized models

## Phase 5: Integration & Testing
- [ ] Create end-to-end `ScreeningPipeline` class
- [ ] Write unit tests for all modules
- [ ] Create tutorial notebook with CCR2 example
- [ ] Document API
- [ ] Benchmark on CCR2 data

---

# 10. References

## Tools & Documentation

- [CDPKit shapescreen](https://cdpkit.org/applications/shapescreen.html)
- [ROSHAMBO2 GitHub](https://github.com/molecularinformatics/roshambo2)
- [RDKit rdShapeHelpers](https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html)
- [Optuna Documentation](https://optuna.readthedocs.io/)
- [pymoo Multi-Objective](https://pymoo.org/)

## MCP Servers

- [Zotero MCP](https://github.com/54yyyu/zotero-mcp)
- [Semantic Scholar MCP](https://github.com/JackKuo666/semanticscholar-MCP-Server)
- [Paper Search MCP](https://github.com/openags/paper-search-mcp)
- [Browser Use MCP](https://docs.browser-use.com/customize/integrations/mcp-server)

## Key Papers

- [PharmacoNet (arXiv)](https://arxiv.org/abs/2310.00681)
- [PharmacoMatch (arXiv)](https://arxiv.org/abs/2409.06316)
- [Pharmacophore Validation Review 2025](https://link.springer.com/article/10.1007/s10822-025-00751-9)
- [Pareto Optimization for VS](https://pubs.rsc.org/en/content/articlehtml/2024/dd/d3dd00227f)
- [SALSA Active Learning](https://arxiv.org/abs/2505.12913)
- [ROC/EF Implementation Notebook](https://gist.github.com/ravila4/26b6ceb21c7e87af80be01f4620a7a58)

## Background Reading

- [3D Pharmacophore Extraction with RDKit](https://www.blopig.com/blog/2025/09/extracting-3d-pharmacophore-points-with-rdkit/)
- [Consensus Virtual Screening](https://pmc.ncbi.nlm.nih.gov/articles/PMC11134635/)
- [OpenEye vROCS Statistics](https://docs.eyesopen.com/applications/rocs/vrocs/vrocs_stats.html)

---

# Session Log

## 2026-01-23: Initial Planning Session

**Accomplished**:
- Created comprehensive plan document
- Identified MCP servers for literature access
- Brainstormed optimization strategies
- Outlined implementation architecture

**Decisions Made**:
- Use CDPKit for shape screening
- Use Optuna with multi-objective optimization
- Use repeated stratified CV (5×10)
- Prioritize both AUC and EF1% via Pareto

**Next Steps**:
1. Install and configure MCP servers
2. Search literature for CCR2 pharmacophore papers
3. Read PharmacoNet paper
4. Begin Phase 1 implementation

---

*Document maintained by Claude Code - Last updated: 2026-01-23*

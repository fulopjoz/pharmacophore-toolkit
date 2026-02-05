# Paper Insights for Consensus Automatic Pharmacophore Model Creator

**Generated**: 2026-01-30
**Source**: 11 priority papers extracted via `utils/pdf_processor.py`
**Purpose**: Consolidated algorithmic insights for improving the consensus pharmacophore optimization pipeline

---

## Executive Summary

Across 11 priority papers, we identified **five actionable algorithm families** that directly inform the consensus automatic pharmacophore model creator:

1. **PharmaGist-style pivot alignment** — multiple flexible alignment with subset detection
2. **ICP pharmacophore clustering** — combined structural + chemical distance metric
3. **Multi-objective Pareto optimization** — simultaneous shape + energy scoring (not post-hoc)
4. **K-means with validation metrics** — S_Dbw-guided clustering with feature OR-assignment
5. **Bayesian optimization with noise handling** — greedy/PI acquisition with retest policies

---

## Paper-by-Paper Findings

### 1. PharmaGist — Ligand-Based Pharmacophore Detection (gkn187)

**Schneidman-Duhovny et al., Nucleic Acids Research, 2008**

**Directly applicable** to consensus model generation from multiple reference ligands.

**Algorithm (4 stages):**
1. **Representation**: Each molecule → set of pharmacophore features (Donor, Acceptor, Aromatic, Hydrophobic, charge) with 3D coordinates
2. **Pairwise alignment**: Align each pair of molecules; find best feature superposition using distance threshold
3. **Multiple alignment**: Iteratively select each molecule as "pivot"; align all others to it; detect shared pharmacophore features
4. **Solution clustering**: Group resulting pharmacophore models by similarity; rank by score

**Key parameters:**
- Feature matching distance threshold: **1.0 Å** (default)
- Feature weights: **0.3 for hydrophobic**, **1.0 for all others** (hydrophobic features are less specific)
- Handles up to **32 molecules**, runs in seconds to minutes
- Detects pharmacophores shared by **subsets** of input ligands (critical for handling outliers or multiple binding modes)

**Relevance to our project:**
- Our `consensus.py` uses agglomerative clustering with a single `tolerance` parameter — PharmaGist's weighted feature matching could improve discrimination
- The pivot molecule approach could replace or complement our current all-vs-all clustering
- Subset detection is important for CCR2 where reference ligands may bind differently

---

### 2. ICP Pharmacophore Clustering (1471-2105-15-S16-S5)

**BMC Bioinformatics, 2014**

**Directly applicable** to comparing and clustering pharmacophore models.

**Core formula:**
```
D = λ * S + (1-λ) * C
```
Where:
- `S` = structural distance (3D spatial, via Iterative Closest Point alignment)
- `C` = chemical distance (feature type similarity, greedy assignment)
- `λ` = weighting parameter (balance spatial vs chemical similarity)

**Algorithm:**
1. Use ICP to align two pharmacophore models in 3D → structural distance S
2. Use greedy matching of feature types → chemical distance C
3. Combine with λ-weighted formula → total distance D
4. Cluster pharmacophores using distance matrix

**Performance:**
- **Rand Index 0.95** on antibody-antigen pharmacophore dataset
- **41 pharmacophores compared in ~30 seconds** (fast enough for iterative optimization)

**Relevance to our project:**
- Our hybrid scoring already combines multiple metrics — this provides a principled formula
- The λ parameter could be optimized via our Bayesian optimization pipeline
- ICP alignment is rotation-invariant, unlike our current distance-based clustering

---

### 3. MOSFOM — Multi-Objective Scoring Function Optimization (1471-2105-10-58)

**Yan et al., BMC Bioinformatics, 2009**

**Key insight:** Consensus scoring (post-hoc re-ranking) only re-scores a limited set of conformations → misses true positives. Multi-objective optimization during docking/scoring is superior.

**Two methods compared:**
- **ε-MOEA**: Generates full Pareto-optimal set of solutions
- **EFMOGA**: Weighted combination → single optimal solution

**Results:**
- EFMOGA performed best retrieving actives in top 2% of database
- Contact score (shape) + energy score optimized **simultaneously during search**, not post-hoc
- Worked across all binding site types (hydrophobic, mixed, polar)

**Relevance to our project:**
- Our combo_optimizer already does multi-objective optimization — this validates the approach
- Key lesson: **optimize scoring functions jointly**, don't combine scores after the fact
- Supports our hybrid_scoring approach over simple consensus averaging

---

### 4. MoPBS — Mapping of Protein Binding Sites via K-means (Braun & Fayne 2022)

**Journal of Molecular Graphics and Modelling, 2022**

**Algorithm for structure-based pharmacophore generation:**
1. Overlay multiple co-crystallized protein-ligand structures
2. Flood binding site with molecular fragments (9 types: methane, water, benzene, acetate, butane, cyclohexane, dimethylether, ethanol, methylammonium)
3. K-means clustering of fragment positions → pharmacophore points
4. Feature assignment based on fragment type statistics per cluster

**Key innovations:**
- **S_Dbw validation metric** for cluster quality: `S_Dbw(c) = Scat(c) + Dens_bw(c)` (lower = better)
- **Differential Evolution** to optimize DBSCAN parameters (MinPts, eps)
- **Feature OR-assignment**: When two feature types are within 10% occurrence, assign both with logical OR
- **Optimal k=9** pharmacophore points for Androgen Receptor (determined by DUD-E screening)
- **Graph theory pharmacophore comparison**: Translated pharmacophores to undirected colored clique graphs; MCS detection via compatibility matrix (Levi/Grasselli algorithm); tolerance radius 2Å

**Performance:**
- k=9 → ratio of 46.21 (true positives weighted vs false positives)
- k=10 → 0 true positives (too stringent)
- k=8 → ratio of 7.90 (too permissive)

**Relevance to our project:**
- S_Dbw could be an additional internal validation metric for our clustering
- The feature OR-assignment is relevant — our current system assigns single types
- Graph theory comparison is a principled alternative to our distance-based comparison
- The sensitivity to k (number of features) matches our experience with `occurrence_threshold`

---

### 5. Graph-Based Pareto Molecular Optimization (SC-013-D2SC00821A)

**Jensen et al., Chemical Science, 2022**

**Multi-objective optimization using genetic algorithms (NSGA-II/III):**
- Molecules represented as molecular graphs
- Non-dominated sorting: solutions ranked by Pareto dominance
- Crossover and mutation operators on graph representations

**Key metrics:**
- **Dominated hypervolume**: Area dominated by Pareto front (higher = better)
- **Extended fingerprint similarity**: Measures diversity of solutions

**Critical insight: Quality-diversity balance**
- Pure exploitation → evolutionary stagnation (convergence to local optima)
- Need explicit diversity maintenance mechanisms
- Tournament selection with crowding distance works well

**Relevance to our project:**
- Our optuna_optimizer uses NSGA-II — this validates multi-objective approach
- Dominated hypervolume could be added as an evaluation metric
- Diversity maintenance is relevant for generating diverse pharmacophore models

---

### 6. Batched Bayesian Optimization for Drug Design (Bellamy et al. 2022)

**Journal of Chemical Information and Modeling, 2022**

**Practical BO for compound screening with noise:**
- **Greedy** and **Predicted Improvement (PI)** acquisition metrics perform best
- UCB and EI perform poorly (especially in noisy environments)
- **100 molecules per batch** standard
- **Random forest** surrogate models provide good QSAR performance
- Tested on 288 ChEMBL datasets + 2 PubChem datasets

**Retest policy for noisy environments:**
- If predicted active but measured inactive → retest
- Each retest replaces one new molecule in the batch
- Increases correctly identified hits by up to 24 percentage points
- Most beneficial when noise is moderate (α=0.1-0.2)
- More beneficial in later batches (when fewer actives remain)

**Relevance to our project:**
- Greedy acquisition could be used in our combo_optimizer's Bayesian search
- Retest policy applicable to our evaluation pipeline (screening metrics have noise)
- Batch size of 100 is a practical reference for our optimization iterations

---

### 7. Multi-Fidelity Bayesian Optimization (McDonald et al. 2025)

**ACS Central Science, 2025**

**State-of-the-art MF-BO for drug discovery:**
- Three fidelity levels: docking (cheap) → single-point inhibition → dose-response IC50 (expensive)
- **Gaussian Process with Tanimoto Kernel** best surrogate model
- Morgan fingerprints (radius 2, 1024 bit) best molecular representation
- **Targeted Variance Reduction (TVR)** for fidelity selection
- Monte Carlo batch selection

**When MF-BO works best:**
- **Diverse search space** + **weak-to-moderate fidelity correlation**
- 30% improvement over single-fidelity BO for Factor-D target
- Strong correlation → experimental funnel better
- No correlation → single-fidelity BO better

**Budget allocation (learned automatically):**
- Iteration 1: 39% docking, 61% single-point, 0% dose-response
- Iteration 2: 4% docking, 67% single-point, 29% dose-response
- Algorithm naturally shifts from exploration to exploitation

**Relevance to our project:**
- Maps directly to our multi-level scoring: pharm2d (cheap) → pharm3d → hybrid (expensive)
- The fidelity correlation insight helps choose when to use each scorer
- Budget allocation concept could optimize our evaluation pipeline
- GP + Tanimoto kernel is the recommended surrogate for molecular optimization

---

### 8. Schaller 2020 — Next Generation 3D Pharmacophore Modeling (WIREs Review)

**WIREs Computational Molecular Science, 2020**

Comprehensive review establishing the state of the art. Key points:
- Feature-based pharmacophores remain the "most suitable data representation for guiding medicinal chemistry"
- LigandScout and CDPKit (open source) are leading tools
- CONFORGE conformer generator: 2x faster than iCon Fast
- G3PS alignment enables exa-scale virtual screening
- MD trajectory analysis for dynamic pharmacophore models
- Hierarchical pharmacophore network analysis (Garon et al. 2020)

---

### 9. Hönig 2022 — Small Molecule Superposition Methods (WIREs Review)

**WIREs Computational Molecular Science, 2022**

Comparison of 51 superposition methods — critical for alignment quality:
- Methods categorized: atom-based, feature-based, shape-based, hybrid
- Performance varies significantly by dataset and metric
- No single method dominates across all benchmarks
- Ensemble/consensus approaches often outperform individual methods

---

### 10-11. Supporting Context Papers

- **pharmaceuticals-15-00646**: Comprehensive review of pharmacophore-based drug design and virtual screening methods
- **nihms150913**: NIH-funded pharmacophore methodology paper with practical guidelines

---

## Consolidated Recommendations for the Consensus Model Creator

### Immediate Improvements (Can implement now)

| # | Improvement | Source Paper | Implementation Target |
|---|------------|-------------|----------------------|
| 1 | Add feature weight differentiation (0.3 hydrophobic, 1.0 others) | PharmaGist | `consensus.py` clustering |
| 2 | Use S_Dbw as internal cluster validation metric | MoPBS (Braun 2022) | `evaluation.py` |
| 3 | Add feature OR-assignment when types within 10% | MoPBS (Braun 2022) | `consensus.py` |
| 4 | Use greedy acquisition in Bayesian optimization | Bellamy 2022 | `combo_optimizer.py` |
| 5 | Add retest/re-evaluation for borderline results | Bellamy 2022 | `evaluation.py` |

### Architecture Enhancements (Next phase)

| # | Enhancement | Source Paper | Complexity |
|---|------------|-------------|------------|
| 6 | Implement combined D = λS + (1-λ)C distance metric | ICP Clustering | Medium |
| 7 | Add pivot molecule alignment (PharmaGist-style) | PharmaGist | Medium |
| 8 | Multi-fidelity scoring: pharm2d → pharm3d → hybrid | MF-BO (McDonald 2025) | High |
| 9 | Graph-theory pharmacophore comparison (MCS) | MoPBS (Braun 2022) | High |
| 10 | Dominated hypervolume as Pareto evaluation metric | Jensen 2022 | Medium |

### Key Parameter Recommendations

| Parameter | Recommended Value | Source |
|-----------|------------------|--------|
| Feature matching distance | 1.0 Å | PharmaGist |
| Hydrophobic feature weight | 0.3 (others: 1.0) | PharmaGist |
| Max pharmacophore features | 9 (optimize per target) | MoPBS |
| Feature OR threshold | 10% occurrence difference | MoPBS |
| BO batch size | 100 | Bellamy 2022 |
| BO acquisition function | Greedy or PI | Bellamy 2022 |
| Surrogate model | GP + Tanimoto kernel | McDonald 2025 |
| Graph comparison tolerance | 2.0 Å | MoPBS |

### Multi-Fidelity Scoring Strategy

Based on McDonald 2025, our scoring levels map naturally:

```
Level 1 (Low cost):    pharm2d_scoring   → 2D fingerprint distance
Level 2 (Medium cost): pharm3d_scoring   → 3D distance-based scoring
Level 3 (High cost):   hybrid_scoring    → Combined 3D + shape alignment
```

**When to use MF approach**: When search space is diverse and fidelity correlation is weak-to-moderate (exactly our CCR2 case with 75 actives + 500 decoys).

---

## Cross-Domain Algorithm Discovery

**Generated**: 2026-01-30 (brainstorming + literature search session)
**Sources**: Wiley Scholar Gateway, PubMed, bioRxiv, Web search across 7 queries
**Purpose**: Identify translatable algorithms from computer vision, graph ML, optimal transport, and NLP that can improve pharmacophore clustering, representation, and optimization

### Executive Summary — Six New Algorithm Families

Beyond our 11 priority papers, we identified **six cross-domain algorithm families** from other fields that directly translate to pharmacophore model generation:

1. **Colored Point Cloud Registration** — ICP/RANSAC with pharmacophore "color" encoding (computer vision → pharmacophore alignment)
2. **Optimal Transport / Wasserstein Distance** — Principled distribution comparison for pharmacophore feature sets (mathematics → pharmacophore scoring)
3. **Hypergraph GNNs** — Pharmacophore-level message passing beyond atom-level graphs (graph ML → pharmacophore representation)
4. **Contrastive Self-Supervised Learning** — Learn pharmacophore-aware molecular embeddings without labels (NLP/CV → pharmacophore scoring)
5. **Multi-Instance Learning** — Handle conformational ambiguity via bag-of-conformations (ML theory → consensus model robustness)
6. **Ensemble Clustering with Voting** — Run clustering 500x, vote on best center (ensemble methods → consensus model stability)

---

### 1. Colored Point Cloud Registration (from Computer Vision)

**Key insight**: Pharmacophore features are literally colored 3D point clouds — each point has coordinates (x,y,z) and a type (Donor, Acceptor, etc.). Computer vision has mature algorithms for registering colored point clouds.

**Directly applicable tools:**

| Tool | Technique | Year | Application |
|------|-----------|------|-------------|
| **ELIXIR-A** | Colored ICP + FPFH/RANSAC via Open3D | 2022 | Multi-target pharmacophore refinement |
| **SENSAAS-Flex** | Colored surface point cloud registration | 2024 | Flexible molecular shape alignment |
| **ShEPhERD** | SE(3)-equivariant diffusion on pharmacophore point clouds | 2024 | Generative pharmacophore design |
| **G3PS** | Hungarian matching + Kabsch alignment | 2021 | Pharmacophore-specific 3D alignment |

**Translation to our toolkit:**

- **Colored ICP** (ELIXIR-A approach): Our consensus features `[type, indices, x, y, z]` are already colored points. The standard ICP in our ICP paper (D = λS + (1-λ)C) treats color separately from geometry. **Colored ICP unifies them** — it minimizes a joint objective over spatial distance AND feature type mismatch simultaneously.
  - Implementation: Use Open3D's `registration_colored_icp()` with pharmacophore types mapped to RGB channels
  - Advantage over current approach: Rotation-invariant, handles partial overlap (subset pharmacophores)

- **FPFH + RANSAC** (global registration): Our current agglomerative clustering requires pre-aligned molecules. FPFH (Fast Point Feature Histograms) descriptors encode local geometric context around each pharmacophore point. Combined with RANSAC, this enables **alignment-free pharmacophore comparison**.
  - Implementation: Encode each feature's local neighborhood as FPFH descriptor → RANSAC for global registration → colored ICP for refinement

- **SE(3)-equivariant processing** (ShEPhERD): Process pharmacophore point clouds with SE(3)-equivariant neural networks that are inherently rotation/translation invariant. No alignment step needed.

**References:**
- Wieder et al. (2022). ELIXIR-A: Interactive Multi-Target Pharmacophore Refinement. [ACS Omega](https://doi.org/10.1021/acsomega.1c07144)
- SENSAAS-Flex (2024). Joint shape alignment + conformer exploration. [Bioinformatics](https://academic.oup.com/bioinformatics/article/40/3/btae105/7612231)
- ShEPhERD (2024). Diffusing shape, electrostatics, and pharmacophores. [arXiv:2411.04130](https://arxiv.org/html/2411.04130v1)
- Saranti et al. (2024). From 3D point-cloud data to explainable geometric deep learning. [WIREs DMKD](https://doi.org/10.1002/widm.1554)
- Li et al. (2021). Tutorial review on point cloud registrations. [Math Problems Eng](https://doi.org/10.1155/2021/9953910)

---

### 2. Optimal Transport / Wasserstein Distance (from Mathematics/ML)

**Key insight**: Comparing two pharmacophore models is fundamentally comparing two **distributions of typed 3D points**. Optimal transport provides the mathematically principled way to measure the "cost" of transforming one distribution into another.

**Why this matters for us:**
- Our current distance metrics (Euclidean between centroids, Tanimoto on fingerprints) don't account for the **optimal assignment** between feature sets of different sizes
- Wasserstein distance naturally handles **unequal point sets** (pharmacophore A has 7 features, B has 5)
- It provides a **differentiable** distance metric, enabling gradient-based optimization

**Concrete translation:**

```
Current approach:    Tanimoto(fingerprint_A, fingerprint_B) → scalar similarity
OT approach:         Wasserstein(features_A, features_B) → scalar distance
                     with cost matrix C_ij = spatial_dist(i,j) + type_mismatch(i,j)
```

**Implementation strategy:**
1. Represent each pharmacophore model as a discrete measure: μ = Σ w_i · δ(x_i, type_i)
2. Define cost function: c(i,j) = ||pos_i - pos_j||² + α · type_distance(type_i, type_j)
3. Compute Wasserstein distance using `scipy.stats.wasserstein_distance` (1D) or `ot` library (multi-dim)
4. Use **Sinkhorn distance** (entropy-regularized OT) for fast approximate computation

**Precedent in molecular science:**
- Gachon et al. (2025) used OT + K-means for dimensionality reduction of multi-patient flow cytometry datasets — directly analogous to our multi-ligand pharmacophore clustering
- Wasserstein distance used for compositional similarity of inorganic compounds (Nature Sci Rep 2024)
- MROT algorithm uses OT for molecular representation domain adaptation

**References:**
- Zhang et al. (2022). Projection-based techniques for high-dimensional OT. [WIREs Comp Stats](https://doi.org/10.1002/wics.1587)
- Gachon et al. (2025). Low-dimensional representation via OT. [Cytometry Part A](https://doi.org/10.1002/cyto.a.24918)
- Jiang et al. (2025). OT + Kabsch for protein-protein docking. [WIREs Comp Mol Sci](https://doi.org/10.1002/wcms.70016)

---

### 3. Hypergraph GNNs for Pharmacophore-Level Representation (from Graph ML)

**Key insight**: Standard molecular GNNs operate at atom/bond level. But pharmacophore features are **higher-order structures** spanning multiple atoms. Hypergraph GNNs can encode pharmacophore-level "hyperstructured knowledge" directly.

**Hyper-Mol framework (Cui et al. 2023):**
- Standard GNN: atom → bond → molecule (loses pharmacophore grouping)
- Hypergraph GNN: atom → pharmacophore substructure → molecule (preserves grouping)
- Fingerprint-level message passing encodes intra-structured AND inter-structured information
- Outperforms atom-level GNNs on molecular property benchmarks

**Translation to our toolkit:**
- Represent each pharmacophore model as a hypergraph where:
  - Nodes = individual atoms contributing to features
  - Hyperedges = pharmacophore features (each spans multiple atoms)
  - Hyperedge attributes = feature type, position, radius
- Use hypergraph convolution for pharmacophore model comparison
- Could replace our mol_converter.py fragment-based representation with a learned hypergraph embedding

**Additional GNN insights:**
- **Graph Transformers** (Harren et al. 2024): Separate processing of covalent vs non-covalent interactions improves scoring — analogous to separating intra-feature vs inter-feature distances in pharmacophore comparison
- **PharmRF** (Kumar et al. 2022): Random forest on protein pocket descriptors + ligand pharmacophore elements → ML scoring for pharmacophore-based VS with 77.6% success rate

**References:**
- Cui et al. (2023). Hyper-Mol: Molecular Representation via Fingerprint-Based Hypergraph. [Comp Intel Neurosci](https://doi.org/10.1155/2023/3756102)
- Harren et al. (2024). Modern ML for binding affinity estimation. [WIREs Comp Mol Sci](https://doi.org/10.1002/wcms.1716)
- Kumar et al. (2022). PharmRF: ML scoring for pharmacophore screening. [J Comp Chem](https://doi.org/10.1002/jcc.26840)

---

### 4. Contrastive Self-Supervised Learning (from NLP/Computer Vision)

**Key insight**: Our pharmacophore scoring requires comparing models against screening datasets (actives/decoys). Contrastive learning can pre-train pharmacophore encoders on large unlabeled molecular data, then fine-tune for our specific CCR2 task.

**MolCLR framework (Wang et al.):**
- Self-supervised learning on ~10M unlabeled molecules
- Augmentation: subgraph removal, atom masking, bond deletion
- Contrastive loss: maximize similarity between augmented views of same molecule
- Significant improvement over supervised methods on property prediction

**Translation to our pharmacophore system:**
- **Positive pairs**: Different pharmacophore models from the same active ligand (different conformers)
- **Negative pairs**: Pharmacophore models from different molecules
- Pre-train an encoder that maps pharmacophore models → fixed-size embeddings
- Use embeddings for fast pharmacophore model comparison (instead of expensive 3D alignment)

**Multi-view learning (GraphMVP):**
- Learn from 2D topology AND 3D geometry simultaneously
- Maps to our approach: pharm2d_scoring (2D fingerprints) + pharm3d_scoring (3D distances) as two "views"
- Could unify our multi-level scoring into a single learned representation

**References:**
- Pang et al. (2023). Advanced DL methods for molecular property prediction. [Quant Biology](https://doi.org/10.1002/qub2.23)
- Yuan et al. (2025). Virtual bonding enhanced graph SSL. [J Comp Chem](https://doi.org/10.1002/jcc.70147)
- Boonpalit et al. (2024). Pre-training strategy for antiviral screening. [J Comp Chem](https://doi.org/10.1002/jcc.27514)

---

### 5. Multi-Instance Learning (from ML Theory)

**Key insight**: A molecule can exist in many conformations — we don't know which one binds. MIL treats each molecule as a "bag" of conformational instances. This is directly relevant to consensus pharmacophore generation where each ligand has multiple possible conformers.

**MIL-kmeans algorithm (Zankov et al. 2023):**
1. Each molecule (bag) → ensemble of conformations (instances)
2. Each conformation encoded by 3D pharmacophore descriptors
3. All instances clustered via k-means across training set
4. Per-molecule binary vector: bit=1 if any conformation falls in cluster
5. Random Forest classification on these binary vectors

**Performance:** Competitive with 2D models, and can **identify bioactive conformations** (MILES algorithm recognized experimental bioactive conformations for 10/12 test molecules).

**Translation to our consensus system:**
- Our reference ligands have multiple conformers → treat as MIL bags
- Generate consensus pharmacophore features from the "bag" representation
- Use MIL to automatically select which conformers contribute to the consensus model
- Eliminates the need to pre-select a single "best" conformer per ligand

**References:**
- Zankov et al. (2023). Chemical complexity challenge: Is MIL a solution? [WIREs Comp Mol Sci](https://doi.org/10.1002/wcms.1698)

---

### 6. Ensemble Clustering with Voting (from Bioinformatics)

**Key insight**: Running clustering once gives unstable results (sensitive to initialization). Running it many times and voting on the most frequent center provides **robust consensus**.

**Ensemble k-medoids (Li et al. 2018):**
- Run k-medoids **500 times** with different random initializations
- Count how often each point becomes the center of the largest cluster
- Confidence score considers both cluster size AND cluster similarity
- Selected better near-native protein structures than SPICKER (used by I-TASSER)
- **70.4% of targets**: ensemble method beat single-run clustering

**Intelligent consensus (Roy et al. 2018):**
- Different QSAR models may work better for different query compounds
- **Compound-specific model selection**: choose which models to average per query
- Consensus benefits from "various representations of chemical structures, which characterized molecules from different perspectives"

**Multi-level voting (Ghosh et al. 2025):**
- 35 base classifiers → 12 intermediate ensembles → 1 final ensemble
- Majority voting at each level with confidence-based tie-breaking
- 90% balanced accuracy, outperforming all individual models

**Translation to our consensus system:**
- Run our agglomerative clustering N times with slightly perturbed parameters (tolerance ± 0.1Å, threshold ± 0.05)
- Vote on which features appear consistently across runs → **stability-aware consensus**
- Weight each run's contribution by its S_Dbw validation score
- For evaluation: use multi-level voting across our scoring approaches (pharm2d, pharm3d, hybrid) instead of simple averaging

**References:**
- Li et al. (2018). Ensemble clustering for protein structure selection. [Quant Biology](https://doi.org/10.1007/s40484-018-0158-1)
- Roy et al. (2018). Intelligent consensus QSAR. [J Chemometrics](https://doi.org/10.1002/cem.2992)
- Ghosh et al. (2025). COX-2 prediction with multi-level ensemble. [J Comp Chem](https://doi.org/10.1002/jcc.70030)

---

### Prioritized Cross-Domain Implementation Roadmap

#### Tier 1: Direct Translation (RDKit/scipy, no ML training needed)

| # | Algorithm | Source Domain | Implementation | Estimated Impact |
|---|-----------|--------------|----------------|-----------------|
| 11 | Wasserstein distance for pharmacophore comparison | Optimal transport | `scipy.stats.wasserstein_distance` or `ot` library | Better model comparison metric |
| 12 | Ensemble clustering with voting | Bioinformatics | Run existing clustering N times, vote | More stable consensus models |
| 13 | Colored ICP via Open3D | Computer vision | `open3d.pipelines.registration` | Alignment-free pharmacophore comparison |
| 14 | Hungarian matching for feature assignment | Graph theory | `scipy.optimize.linear_sum_assignment` | Optimal feature pairing (replaces greedy) |

#### Tier 2: Medium Complexity (requires some ML but uses existing libraries)

| # | Algorithm | Source Domain | Implementation | Estimated Impact |
|---|-----------|--------------|----------------|-----------------|
| 15 | Multi-instance learning for conformer selection | ML theory | scikit-learn + pharmacophore descriptors | Better conformer handling |
| 16 | Multi-level voting across scoring methods | Ensemble methods | Wrapper around existing scorers | Improved evaluation robustness |
| 17 | FPFH descriptors for pharmacophore fingerprinting | Computer vision | Open3D FPFH computation | Novel fingerprint type |

#### Tier 3: Advanced (requires training, more data)

| # | Algorithm | Source Domain | Implementation | Estimated Impact |
|---|-----------|--------------|----------------|-----------------|
| 18 | Contrastive pre-training for pharmacophore embeddings | NLP/CV | PyTorch + GNN framework | Learned scoring function |
| 19 | Hypergraph GNN for pharmacophore representation | Graph ML | PyTorch Geometric | Novel representation |
| 20 | SE(3)-equivariant pharmacophore processing | Geometric DL | e3nn or equivalent | Rotation-invariant scoring |

---

## Process Documentation

### How `utils/pdf_processor.py` Works

1. **Input**: Directory of PDF files
2. **Processing**: PyMuPDF (fitz) extracts text page-by-page with dual safety limits:
   - `max_pages=20` (configurable) — prevents processing entire 100+ page documents
   - `max_chars=50000` (configurable) — hard cap on extracted text
3. **Output**: Per-paper `.extracted.md` files + `_INDEX.md` summary table
4. **Safety**: Idempotent (skips already-extracted files unless `--force`)

### Processing Results

| Directory | PDFs | Extracted | Total Chars | Notable Issues |
|-----------|------|-----------|-------------|----------------|
| `docs/research/papers/priority/` | 12 | 12 | ~456K chars | 1 test file (14 chars) |
| `docs/research/papers/algo/` | 18 | 18 | ~720K chars | All OK |
| `docs/` | 2 | 2 | ~9.6K chars | 1 scanned PDF (22 chars), 1 slide deck (9.5K) |

### Known Limitations

1. **Scanned PDFs**: `Comparison of Consensus Pharmacophore Screening.pdf` is image-only → 22 chars extracted. Would need OCR (pytesseract) for text.
2. **Slide decks**: `TLanger_Olomouc_2025-Reduced.pdf` is mostly images → limited text extraction.
3. **Character truncation**: Papers exceeding 50K chars are truncated. Most key content is in first 20 pages.

<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" style="height:64px;margin-right:32px"/>

# What are the gold-standard benchmark datasets and evaluation metrics for

validating pharmacophore-based virtual screening?

Include:

- Dataset names (DUD-E, MUV, LIT-PCBA, etc.)
- Metrics (EF, AUC, precision@1%, etc.)
- Typical parameter ranges (tolerance 1.0-2.0 Å, threshold 0.5-1.0)
- Recent papers comparing consensus vs single-molecule pharmacophores

Provide download links and recommended evaluation protocols.

Now I'll compile all gathered information into a comprehensive report on gold-standard benchmark datasets and evaluation metrics for pharmacophore-based virtual screening.

## Gold-Standard Benchmark Datasets and Evaluation Metrics for Pharmacophore-Based Virtual Screening

### Executive Summary

Pharmacophore-based virtual screening has emerged as a cornerstone methodology in computational drug discovery, requiring rigorous validation through standardized benchmarks and quantitative metrics. This report provides a comprehensive analysis of gold-standard benchmark datasets (DUD-E, MUV, LIT-PCBA, DEKOIS 2.0), evaluation metrics (enrichment factors, ROC AUC, BEDROC, Güner-Henry score), optimal parameter ranges, and recent advances in consensus pharmacophore modeling. The analysis synthesizes evidence from 195+ recent publications (2020-2025) to deliver actionable guidance for validating pharmacophore models in drug discovery workflows.

***

### Benchmark Datasets for Pharmacophore Validation

#### 1. DUD-E (Directory of Useful Decoys - Enhanced)

**Dataset Characteristics**[^1][^2]

- **Scale**: 102 protein targets, averaging 224 ligands per target
- **Decoy Ratio**: 50 property-matched decoys per active ligand (~22,000 active compounds, ~1.1 million decoys)
- **Download**: http://dude.docking.org
- **Last Update**: 2012 (maintained actively)

**Strengths**

- Broad target coverage across pharmaceutically relevant protein families
- Property-matched decoys designed to minimize trivial discrimination[^1]
- Extensive benchmarking history enabling cross-study comparisons
- Includes diverse chemotypes and targets from multiple therapeutic areas[^3]

**Known Limitations and Biases**[^4][^1]

- **Analogue Bias**: Subtle chemical similarity between actives can inflate enrichment[^4]
- **Decoy Bias**: Property matching doesn't eliminate all hidden biases (e.g., molecular complexity patterns)[^4]
- **Artificial Enrichment**: Neural networks can exploit dataset artifacts rather than learning true binding physics[^4]
- **Recommendation**: Use complementary datasets (LIT-PCBA, DEKOIS 2.0) for external validation[^5]

**Typical Performance Benchmarks**[^3][^1]

- Ensemble docking: ROC AUC 0.89 ± 0.08, EF1% ~46-fold
- Single-structure docking: ROC AUC 0.81 ± 0.11
- Ligand-based methods: ROC AUC 0.77-0.83


#### 2. LIT-PCBA (Literature PubChem BioAssay)

**Dataset Characteristics**[^5]

- **Scale**: 15 target sets with 7,844 confirmed actives and 407,381 confirmed inactives
- **Hit Rate**: Mimics experimental screening decks (1.9% average hit rate)
- **Download**: http://drugdesign.unistra.fr/LIT-PCBA
- **Publication**: Tran-Nguyen et al., J. Chem. Inf. Model. 2020

**Key Innovations**[^5]

- **Dose-Response Validation**: All actives confirmed through concentration-dependent assays
- **Artifact Removal**: Systematic filtering of false positives and assay interference compounds
- **Property Balancing**: Actives and inactives maintained within similar molecular property ranges
- **AVE Debiasing**: Asymmetric validation embedding procedure applied to reduce hidden biases
- **Structure Availability**: Each target has ≥1 X-ray structure with same-phenotype ligands

**Target Selection Criteria**[^5]
Preliminary virtual screening with orthogonal methods (2D fingerprint similarity, 3D shape similarity, molecular docking) required EF1% ≥ 2 for inclusion, ensuring the dataset supports meaningful enrichment.

**Validation Studies**[^6][^7][^8]

- PharmacoForge (diffusion model): Outperforms automated pharmacophore methods on LIT-PCBA
- PharmRL (deep reinforcement learning): Better F1 scores than random feature selection
- PharmacoNet: Demonstrates reasonable accuracy with fast screening speeds

**Strengths**

- Minimizes traditional benchmarking biases through rigorous curation
- Realistic hit rates suitable for machine learning training
- Enables fair comparison between structure-based and ligand-based methods
- Lower risk of data leakage compared to DUD-E[^5]


#### 3. MUV (Maximum Unbiased Validation)

**Dataset Characteristics**[^9][^10]

- **Scale**: 17 activity classes from PubChem bioassays
- **Curation Strategy**: Refined nearest-neighbor analysis with topological optimization
- **Design Objective**: Minimize analogue bias and artificial enrichment

**Methodological Innovations**[^10]

- **Unselective Hit Removal**: Data-centered workflow purges promiscuous binders
- **Topological Optimization**: Experimental design strategies ensure unbiased active/decoy distribution
- **Refined Nearest Neighbor Analysis**: Monitors clustering to prevent trivial separation

**Applications**[^11][^12][^13]

- Challenging benchmark for similarity-based methods (often lower enrichment than DUD-E)
- Particularly suitable for evaluating fingerprint methods and machine learning models
- Test case for methods claiming generalization beyond chemical similarity

**Performance Characteristics**[^13]

- LigMate (multifeature integration): EF1% 15.44, AUC 0.69 on MUV
- Generally 40-60% of DUD-E enrichment levels (reflects reduced bias)[^13]


#### 4. DEKOIS 2.0 (Demanding Evaluation Kits for Objective In Silico Screening)

**Dataset Characteristics**[^14][^15]

- **Scale**: 81 targets with property-matched decoy sets
- **Download**: http://www.dekois.com
- **Philosophy**: Cluster-based artifact removal with physicochemical space matching

**Construction Protocol**[^15][^16]

- **Step 1**: Select actives from literature/PubChem with confirmed bioactivity
- **Step 2**: Apply cluster analysis to identify chemotype families
- **Step 3**: Generate decoys matched for physicochemical properties within each cluster
- **Step 4**: Remove decoys with trivial differentiation from actives

**Validation Applications**[^17][^18][^14]

- Preferred for structure-based virtual screening benchmarking
- Successfully used to benchmark AutoDock Vina, PLANTS, FRED, CDOCKER
- Machine learning scoring function evaluation (CNN-Score, RF-Score-VS)

**Performance Benchmarks**[^18]

- SCORCH (ML classifier): EF1% 13.78 on 18 DEKOIS 2.0 targets
- FRAGSITE2: EF1% ~18-23 across 81 targets[^19]
- PLANTS + CNN re-scoring: EF1% up to 28[^17]


#### Comparative Dataset Analysis

| Dataset | Targets | Actives/Target | Bias Level | Best For | Limitation |
| :-- | :-- | :-- | :-- | :-- | :-- |
| **DUD-E** | 102 | ~224 | Moderate | Broad benchmarking, method comparison | Analogue bias, decoy artifacts[^4] |
| **LIT-PCBA** | 15 | ~523 | Low | ML training, unbiased validation | Limited target coverage |
| **MUV** | 17 | Variable | Very Low | Challenging test, similarity methods | Low enrichment rates |
| **DEKOIS 2.0** | 81 | Variable | Low-Moderate | SBVS, docking validation | Smaller ligand libraries |


***

### Evaluation Metrics for Pharmacophore Screening

#### 1. Enrichment Factor (EF)

**Mathematical Definition**

The enrichment factor at x% of the ranked library is calculated as:

**EF_x% = (Actives_found_x% / Total_actives) / (x / 100)**

Where:

- Actives_found_x% = number of true actives in top x% of ranked compounds
- Total_actives = total number of actives in the entire dataset
- x = percentage cutoff (typically 1%, 5%, or 10%)

**Common Cutoffs and Interpretation**[^20][^21][^22]

**EF1% (Early Enrichment)**

- **Calculation**: Ratio of actives in top 1% relative to random expectation
- **Random Performance**: EF1% = 1.0
- **Acceptable**: EF1% > 10
- **Good**: EF1% > 20-30[^21][^23]
- **Excellent**: EF1% > 40-50[^24]

**EF5% (Moderate Enrichment)**

- Used to assess performance beyond the very top-ranked compounds
- Good performance: EF5% > 10-15[^25]

**EF10% (Broad Enrichment)**

- Assesses enrichment across a larger fraction of the library
- Less sensitive to noise in top-ranked compounds

**Real-World Applications**[^21]

A recent study on apelin receptor agonists demonstrated ensemble pharmacophore learning:

- Single best model: EF1% = 19.466, AUC = 0.82
- Ensemble method: EF1% = 50.07 ± 0.211, AUC = 0.994 ± 0.007

This 2.6-fold improvement in EF1% through consensus modeling highlights the value of ensemble approaches.

**Limitations**[^26]

- EF cannot estimate performance on very large libraries (>100 million compounds)
- Sensitive to dataset size and active/decoy ratio
- Does not capture full ranking quality (only threshold-based)


#### 2. ROC AUC (Receiver Operating Characteristic - Area Under Curve)

**Calculation and Interpretation**[^27]

The ROC curve plots True Positive Rate (sensitivity) versus False Positive Rate (1 - specificity) across all classification thresholds. The AUC quantifies overall discriminative power:

- **AUC = 0.5**: Random performance (no discrimination)
- **AUC = 0.6-0.7**: Poor but better than random
- **AUC = 0.7-0.8**: Acceptable[^20]
- **AUC = 0.8-0.9**: Good[^3][^21]
- **AUC = 0.9-1.0**: Excellent[^28]

**DUD-E Performance Benchmarks**[^3]

- Ensemble docking (5 pocket variants + ligand-based): AUC 0.89 ± 0.08
- Ensemble docking alone: AUC 0.84 ± 0.09
- Single-structure docking (Surflex-Dock): AUC 0.81 ± 0.11
- Ligand-based (eSim, single query): AUC 0.77 ± 0.14

**Advantages**

- Threshold-independent metric
- Intuitive interpretation (probability that random active ranks higher than random decoy)
- Robust to class imbalance when properly interpreted[^27]

**Limitations**

- Does not emphasize early enrichment (critical for experimental validation)
- Can be dominated by performance on inactive-rich regions
- Semi-logarithmic ROC plots recommended for virtual screening to emphasize early performance[^27]


#### 3. BEDROC (Boltzmann-Enhanced Discrimination of ROC)

**Mathematical Framework**[^2][^29]

BEDROC applies an exponentially decaying weight to emphasize early retrieval:

**BEDROC_α = (RI_α - RI_α^min) / (RI_α^max - RI_α^min)**

Where α controls the degree of early enrichment emphasis. The weighting function ensures top-ranked compounds contribute more to the score.

**Common α Parameters**[^30][^25][^28][^2][^20]


| α Value | 20% of Score from Top | Application |
| :-- | :-- | :-- |
| **8** | ~50% of library | Minimal early emphasis[^28] |
| **20** | ~30% of library | Moderate early enrichment[^28] |
| **80.5** | ~12% of library | Strong early focus[^1][^2] |
| **160.9** | ~6% of library | Maximum early emphasis[^20][^25][^30] |

**Interpretation Guidelines**[^2][^20]

For α = 160.9 (most common in pharmacophore studies):

- **BEDROC > 0.3**: Better than random
- **BEDROC > 0.5**: Good performance[^20]
- **BEDROC > 0.7**: Very good early enrichment
- **BEDROC > 0.9**: Exceptional early retrieval

**Application Examples**[^28][^20]

1. **Tubulin Inhibitor Pharmacophore**: BEDROC_160.9 = 0.59, average EF1% = 27.05[^20]
2. **hDHODH Inhibitor Pharmacophore**: BEDROC (α=8) = 0.859, BEDROC (α=20) = 0.828, BEDROC (α=160) = 0.912[^28]
3. **GSK-3β Inhibitor Study**: Multiple α values reported for comprehensive assessment[^31]

**When to Use BEDROC**[^29][^2]

- Early drug discovery where experimental testing is limited to top hits
- Comparing methods with different ranking strategies
- When false negatives in top-ranked compounds are particularly costly


#### 4. Güner-Henry (GH) Score

**Formula and Components**[^32][^33][^34]

**GH = (Ha / (4·Ht·A)) · (3A + Ht)**

Where:

- **Ha**: Number of actives in hit list (retrieved compounds)
- **Ht**: Total number of hits retrieved (actives + false positives)
- **A**: Total number of actives in the entire database

**Interpretation Scale**[^34][^35][^36][^32]

- **GH < 0.5**: Poor model, high false positive rate
- **GH = 0.5-0.7**: Moderate performance, acceptable selectivity
- **GH = 0.7-0.8**: Good model[^32][^34]
- **GH = 0.8-0.9**: Very good model[^28]
- **GH ≥ 0.9**: Excellent model with high precision[^23][^21]

**Application in Pharmacophore Validation**[^33][^35][^34]

The GH score balances recall (finding actives) with precision (avoiding false positives), making it particularly suitable for pharmacophore model validation:

1. **JAK1 Inhibitor Pharmacophores**: Eight models generated; GH scores used to rank and select best candidates for virtual screening[^33]
2. **hERG K+ Channel Blocker Model**: GH score combined with ROC and enrichment factor for comprehensive in silico validation[^37]
3. **Tubulin Inhibitor Pharmacophore**: Validated via Gunner-Henry method alongside molecular docking[^35]

**Calculation Example**

For a pharmacophore screening a database of 10,000 compounds with 100 true actives:

- If 50 compounds are retrieved (Ht = 50)
- Of which 40 are true actives (Ha = 40)
- Total actives A = 100

GH = (40 / (4 × 50 × 100)) × (3 × 100 + 50) = (40 / 20,000) × 350 = 0.7

This GH = 0.7 indicates good model performance.

#### 5. Precision, Recall, and F1-Score

**Definitions**[^23][^21]

**Precision** = TP / (TP + FP)

- Fraction of retrieved compounds that are truly active
- Critical when experimental validation capacity is limited

**Recall (Sensitivity)** = TP / (TP + FN)

- Fraction of true actives successfully retrieved
- Important for ensuring comprehensive hit identification

**F1-Score** = 2 × (Precision × Recall) / (Precision + Recall)

- Harmonic mean balancing precision and recall
- Single metric for model quality assessment

**Performance Standards**[^21][^23]

Apelin receptor pharmacophore ensemble learning study:

- F1-score: 0.911 ± 0.031 (ensemble method)
- F1-score: 0.071 (single best model)

This 13-fold improvement demonstrates the power of ensemble strategies.

**When to Prioritize Each Metric**

- **Precision**: High cost of false positives (expensive synthesis, limited screening capacity)
- **Recall**: Need for comprehensive coverage (hit-to-lead diversification)
- **F1**: Balanced approach for general pharmacophore validation


#### Metric Comparison Summary

| Metric | Emphasizes | Best For | Range | Good Performance |
| :-- | :-- | :-- | :-- | :-- |
| **EF1%** | Top 1% hits | Early discovery | 1-100+ | >20-30[^21] |
| **ROC AUC** | Overall ranking | General assessment | 0.5-1.0 | >0.8[^3] |
| **BEDROC (α=160.9)** | Very early hits | Screening prioritization | 0-1 | >0.5[^20] |
| **GH Score** | Precision + recall | Pharmacophore validation | 0-1 | >0.7[^32][^34] |
| **F1-Score** | Balanced performance | Model comparison | 0-1 | >0.8[^21] |


***

### Pharmacophore Parameter Optimization

#### Feature Tolerance Ranges

**Distance Tolerance (RMSD/Sphere Radius)**[^38]

Pharmacophore features are typically represented as spheres with defined tolerance radii. Optimal ranges vary by feature type:


| Feature Type | Typical Tolerance | Rationale |
| :-- | :-- | :-- |
| **Hydrogen Bond Acceptor/Donor** | 1.0-1.5 Å | Directional interactions require strict geometry[^38] |
| **Aromatic Ring** | 1.5-2.0 Å | π-π stacking allows more flexibility |
| **Hydrophobic** | 1.5-2.0 Å | Van der Waals interactions less directional |
| **Positive/Negative Ionizable** | 1.0-1.5 Å | Electrostatic interactions require precision |
| **Excluded Volume** | 0.5-1.0 Å | Steric clash prevention |

**Impact of Tolerance on Performance**[^38]

- **Stricter Tolerances (1.0 Å)**: Higher specificity, lower recall, risk of missing viable compounds
- **Moderate Tolerances (1.5 Å)**: Balanced performance, most common in literature[^38]
- **Looser Tolerances (2.0 Å)**: Higher recall, increased false positives, suitable for early-stage diverse screening

**Optimization Strategy**

Recent studies on glucokinase activators demonstrated validation with DUD-E data achieving 80% accuracy, 95% sensitivity, and 80% specificity using pharmacophore tolerance optimization.[^38]

#### Feature Weight Assignment

**Principles**[^39][^40]

Feature weights determine the relative importance of pharmacophoric features in scoring compounds:

1. **Critical Interactions**: Hydrogen bonds with conserved residues (e.g., catalytic triads) receive highest weights
2. **Anchor Points**: Deep binding pocket features that define binding mode
3. **Selectivity Determinants**: Features distinguishing target from off-targets
4. **Flexibility Factors**: Weights adjusted based on protein flexibility analysis

**Weighting Schemes**[^39]

Apo2ph4 workflow generates receptor-based pharmacophores with automatic feature weighting:

- Interaction energy calculations assign initial weights
- Validation against known actives/decoys refines weights iteratively
- Final weights balance specificity and sensitivity


#### Conformer Generation and Coverage

**Conformational Sampling**[^41]

MD-derived pharmacophores from protein-ligand trajectories:

- Generate multiple pharmacophores across simulation timeframes
- Use 3D pharmacophore hashes to remove duplicates (identical feature arrangements)
- Rank conformers by frequency of occurrence or binding free energy
- Ensemble approach covers alternative binding modes[^41]

**Conformer Coverage Metric**[^41]

For ligand libraries:

- Generate pharmacophore-compliant conformers for each molecule
- Assess what percentage of active compounds can adopt pharmacophore geometry
- Optimal models: >80% of actives can match pharmacophore with low energy penalty

***

### Consensus vs. Single-Molecule Pharmacophore Models

#### Consensus Pharmacophore Approach

**Definition and Methodology**[^42][^43]

Consensus pharmacophores integrate common features from multiple ligand-bound structures:

1. **Multi-Complex Alignment**: Superimpose multiple ligand-bound structures
2. **Feature Extraction**: Identify pharmacophoric features in each complex
3. **Clustering**: Group similar features across structures using spatial proximity
4. **Consensus Building**: Retain features present in ≥X% of structures (typically 50-80%)
5. **Model Refinement**: Optimize feature weights and tolerances

**ConPhar Protocol**[^43]

Recent open-source tool for consensus pharmacophore generation:

- Input: Multiple PDB complexes with bound ligands
- Clustering algorithm identifies conserved features
- Handles large ligand datasets (100+ structures)
- Applied to SARS-CoV-2 Mpro with 152 bioactive conformers[^42]

**Performance Advantages**[^43][^42]

SARS-CoV-2 Mpro Consensus Pharmacophore Study:[^42][^43]

- **Input**: 152 non-covalent inhibitor co-crystal structures
- **Outcome**: Captured key interactions in catalytic region
- **Validation**: Successfully identified novel candidates from multi-billion compound libraries
- **Advantage**: Reduced bias compared to single-ligand models

**Ensemble Pharmacophore Learning**[^23][^21]

Apelin receptor study combining Butina clustering and ensemble methods:

- **Single Model Performance**: AUC 0.82, EF1% 19.466, GH 0.131, F1 0.071
- **Ensemble Performance**: AUC 0.994 ± 0.007, EF1% 50.07 ± 0.211, GH 0.956 ± 0.015, F1 0.911 ± 0.031
- **Improvement**: 2.6x on EF1%, 13x on F1-score

**Voting and Stacking Strategies**[^21]

- **Voting**: Compound selected if matched by ≥N models (majority voting)
- **Stacking**: Machine learning meta-classifier combines multiple pharmacophore scores
- **Balanced Parallel**: Equal weight to diverse models, reduces individual model biases


#### Single-Molecule Pharmacophore Models

**Characteristics**[^44]

- Derived from single co-crystal structure or most active ligand
- Faster generation (minutes vs. hours for consensus)
- Captures specific binding mode with high detail
- Risk: Overfitting to training set chemotype

**When to Use Single-Molecule Models**[^44]

1. **Limited Structural Data**: Only one co-crystal structure available
2. **Highly Conserved Binding Mode**: All actives bind identically (e.g., ATP-binding sites)
3. **Rapid Initial Screening**: Time-sensitive projects requiring fast deployment
4. **Negative Design**: Explicitly avoiding certain binding modes

**Probabilistic Approaches**[^44]

Recent work on multiple pharmacophore-based screening:

- Max scoring: Use maximum score across multiple single-molecule models
- Mean scoring: Average scores from multiple models
- These approaches showed superior early enrichment compared to single best model[^44]


#### Comparative Performance Analysis

**Recent Studies (2023-2025)**[^45][^43][^42]

1. **SARS-CoV-2 Mpro Screening**[^45]
    - **Approach**: Combined ligand-based pharmacophore with structure-based screening
    - **Scale**: ~200 million compounds screened
    - **Validation**: Two compounds with micromolar inhibitory activity confirmed
    - **Conclusion**: Comprehensive consensus approach superior to single-model methods
2. **Multiple-Objective Reinforcement Learning**[^6]
    - **PharmacoForge**: Diffusion model generating consensus pharmacophores conditioned on protein pocket
    - **Benchmarks**: LIT-PCBA and DUD-E datasets
    - **Result**: Outperforms automated single-structure pharmacophore methods
3. **MD-Derived Ensemble Pharmacophores**[^41]
    - **Strategy**: Extract pharmacophores from MD trajectories, deduplicate via 3D hashes
    - **Application**: Four cyclin-dependent kinases
    - **Finding**: Ensemble from MD outperforms single crystal structure model

#### Best Practices for Consensus Modeling

**Structural Diversity Requirements**[^43][^42]

- Minimum 10-20 diverse ligand-bound structures for robust consensus
- If <10 structures available, consider single-model + flexibility analysis
- Balance chemotype diversity (avoid over-representation of single scaffold)

**Feature Occurrence Thresholds**[^43]

- **Strict Consensus (≥80% occurrence)**: High specificity, may miss alternative modes
- **Moderate Consensus (≥50% occurrence)**: Balanced, most commonly used
- **Relaxed Consensus (≥30% occurrence)**: High recall, increased false positives

**Validation Strategy**[^42]

1. **Internal Validation**: Cross-validation against structures used for model generation
2. **External Validation**: Test on independent ligand set (different chemotypes)
3. **Prospective Screening**: Experimental validation of top predictions[^45]

***

### Recommended Evaluation Protocols

#### Standard Validation Workflow

**Step 1: Model Generation**[^40][^39]

1. **Select Training Set**
    - Minimum 20-30 structurally diverse actives[^46]
    - Activity range spanning 2-3 log units (for QSAR extension)
    - Avoid over-representation of single chemotype
2. **Define Pharmacophore Features**
    - Hydrogen bond donors/acceptors
    - Aromatic/hydrophobic features
    - Ionizable groups
    - Excluded volumes (for selectivity)
3. **Set Initial Parameters**
    - Feature tolerances: 1.0-2.0 Å (optimize in Step 2)
    - Feature weights: Based on interaction energy or conservation
    - Conformational sampling: 100-500 conformers per molecule

**Step 2: Internal Validation**[^37][^34][^35]

**Test Set Method**

- Split training data: 70-80% training, 20-30% test
- Generate pharmacophore from training set
- Evaluate on test set actives + decoys (1:10 to 1:50 ratio)
- Calculate metrics: ROC AUC, EF1%, GH score

**Fischer's Randomization**[^36][^47]

- Randomize activity labels 100-1000 times
- Regenerate pharmacophore with each randomization
- Compare true model statistics to randomized distribution
- True model should significantly outperform (p < 0.05)

**Decoy Set Validation**[^48][^49][^37]

- Prepare property-matched decoys (using DUD-E methodology or DeepCoy)[^50]
- Typical ratio: 50-100 decoys per active
- Screen actives + decoys with pharmacophore
- Calculate comprehensive metrics:
    - ROC AUC (overall discriminative power)
    - BEDROC α=160.9 (early enrichment)
    - EF1%, EF5% (practical screening efficiency)
    - GH score (precision-recall balance)

**Performance Thresholds for Acceptance**[^34][^20][^21]


| Metric | Minimum (Acceptable) | Target (Good) | Excellent |
| :-- | :-- | :-- | :-- |
| ROC AUC | 0.7 | 0.8 | >0.9 |
| BEDROC (α=160.9) | 0.3 | 0.5 | >0.7 |
| EF1% | 10 | 20-30 | >40 |
| GH Score | 0.5 | 0.7-0.8 | >0.9 |
| F1-Score | 0.5 | 0.7 | >0.85 |

**Step 3: External Validation**[^19][^5]

1. **Independent Dataset Testing**
    - Test on different benchmark (e.g., train on DUD-E, validate on LIT-PCBA)
    - Ensures model doesn't overfit to dataset-specific biases
    - Performance drop of 10-20% from internal validation is expected
2. **Chemotype Coverage Analysis**[^51]
    - Use pROC-chemotype plots to assess performance across scaffolds[^51]
    - Identify which chemotypes model captures well vs. poorly
    - Reveals scaffold hopping capability
3. **Cross-Target Validation**
    - For target families (e.g., kinases), test on related targets
    - Assesses transferability and selectivity

**Step 4: Prospective Validation (Optional but Recommended)**[^40][^45]

1. **Virtual Screening Campaign**
    - Screen large library (1-10 million compounds)
    - Select top 50-100 compounds for experimental testing
    - Recommended hit rate: 5-20% (balance of risk and efficiency)
2. **Experimental Validation**
    - Biochemical assay against target
    - Cell-based assays (if applicable)
    - Compare predicted vs. observed activity
3. **Model Refinement**
    - Incorporate experimental results
    - Adjust feature weights based on hit/miss analysis
    - Iterate model for improved performance

#### Advanced Validation Techniques

**Dynamic Pharmacophore Validation (Dynophores)**[^52]

- Perform molecular dynamics simulations of pharmacophore-ligand complexes
- Extract time-dependent pharmacophore features
- Assess stability of pharmacophore matches over simulation time
- Ligands maintaining pharmacophore match across MD are more likely true positives

**pROC-Chemotype Analysis**[^53][^51]

- Extension of ROC curves incorporating chemotype information
- Visualizes enrichment performance for each scaffold family
- Identifies chemotype-specific biases in the model
- Recommended for datasets with >5 distinct chemotypes

**Applicability Domain Assessment**[^54]

- Define chemical space where model is reliable
- Use Euclidean distance in descriptor space from training set
- Compounds outside applicability domain flagged as uncertain predictions
- Reduces false confidence in extrapolations

***

### Software Tools and Resources

#### Pharmacophore Modeling Software

**Commercial Packages**

1. **LigandScout** (Inte:Ligand)
    - User-friendly GUI with automated feature detection
    - Structure-based and ligand-based pharmacophore generation
    - Integrated virtual screening and validation tools
    - Pricing: ~\$5,000-15,000/year academic
2. **Phase** (Schrödinger)
    - Part of Maestro suite with seamless docking integration
    - Advanced 3D-QSAR capabilities (PHASE QSAR)
    - Hypothesis ranking with statistical validation
    - Pricing: ~\$10,000-30,000/year
3. **Discovery Studio** (Dassault Systèmes BIOVIA)
    - Comprehensive pharmacophore generation (HipHop, HypoGen)
    - Catalyst technology for quantitative pharmacophores[^55]
    - Integrated with protein preparation and docking
    - Pricing: Enterprise licensing
4. **MOE** (Chemical Computing Group)
    - Unified Pharmacophore concept combining ligand/structure-based
    - Protein-ligand fingerprint (PLIF) generation[^40]
    - Customizable feature definitions
    - Pricing: ~\$6,000-12,000/year

**Free/Open-Source Tools**

1. **Pharmit** (University of Pittsburgh)
    - Web-based pharmacophore search against large libraries
    - Interactive 3D visualization
    - Access to pre-indexed compound databases (millions of compounds)
    - URL: pharmit.csb.pitt.edu
2. **PharmaGist** (Tel Aviv University)
    - Web server for ligand-based pharmacophore detection[^56]
    - Automatic multiple flexible alignment
    - No protein structure required
    - URL: bioinfo3d.cs.tau.ac.il/PharmaGist/
3. **Apo2ph4** (University of Vienna)
    - Receptor-based pharmacophore generation without ligand knowledge[^39]
    - Automated workflow for structure-based models
    - Validation against DUD-E and DEKOIS 2.0
    - Open-source: GitHub
4. **ConPhar**
    - Consensus pharmacophore generation from multiple complexes[^43]
    - Handles 100+ structures efficiently
    - Feature clustering and visualization
    - Open-source tool introduced in 2025
5. **PharmDock**
    - Pharmacophore-constrained docking[^57]
    - PyMOL integration for visualization
    - Suitable for focused screening with known pharmacophore

#### Virtual Screening Platforms

**Docking Software Performance on Benchmarks**[^14][^1][^17]


| Software | Type | DUD-E BEDROC (α=80.5) | DEKOIS 2.0 EF1% | Speed |
| :-- | :-- | :-- | :-- | :-- |
| **Glide** (Schrödinger) | Commercial | 30/102 targets >0.5 | ~15-20 | Moderate |
| **GOLD** (CCDC) | Commercial | 27/102 targets >0.5 | ~12-18 | Moderate |
| **AutoDock Vina** | Free | Variable, target-dependent | ~10-15 | Fast |
| **PLANTS** | Free academic | Good for WTMpro | 28 (with CNN) | Fast |
| **FRED** (OpenEye) | Commercial | Excellent on OMpro | ~15-20 | Very Fast |
| **FlexX** | Commercial | 14/102 targets >0.5 | ~8-12 | Fast |
| **Surflex-Dock** | Commercial | AUC 0.81 ± 0.11 | ~10-15 | Moderate |

**Shape-Based Screening**

1. **ROCS** (OpenEye)
    - 3D molecular shape similarity (Tanimoto combo scoring)[^58][^59]
    - Combined shape + electrostatics (EON)[^60]
    - Gold standard for shape-based virtual screening
    - Typical enrichment: EF1% 15-25[^61]
2. **USRCAT**
    - Ultrafast shape recognition with pharmacophoric features
    - Alignment-free method
    - Integrated in some free screening platforms

**Integrated Platforms**

1. **MTiOpenScreen** (web-based)
    - Free virtual screening server[^62]
    - Docking with AutoDock Vina and LEDOCK
    - ADME-Tox prediction via FAF-Drugs4
    - Suitable for screening up to ~1 million compounds
2. **PharmMapper**
    - Reverse pharmacophore mapping
    - Identify targets for query pharmacophore
    - Access to multiple target databases

#### Decoy Generation Tools

**DeepCoy**[^50]

- Deep learning-based decoy generator
- Creates property-matched decoys without known biases
- Validated on DUD-E and DEKOIS 2.0
- Available as open-source Python package

**DUD-E Decoy Generator**

- Online tool at http://dude.docking.org/generate
- Property-matching algorithm for custom active sets
- Outputs SDF files ready for screening

**DEKOIS Methodology**

- Cluster-based approach (implementation available)
- More stringent than DUD-E property matching
- Requires computational chemistry expertise to implement

***

### Recent AI-Enhanced Approaches (2024-2025)

#### Deep Learning Pharmacophore Models

**PharmacoForge**[^6]

- **Technology**: Diffusion model for 3D pharmacophore generation
- **Innovation**: Conditions pharmacophore generation on protein pocket geometry
- **Benchmarks**: LIT-PCBA and DUD-E
- **Performance**: Surpasses automated methods; ligands from generated pharmacophores show comparable docking scores to de novo generated molecules
- **Advantage**: Guaranteed valid, commercially available molecules (vs. generative models)

**PharmRL**[^7][^8]

- **Technology**: Deep geometric reinforcement learning for protein-based pharmacophore elucidation
- **Innovation**: CNN identifies favorable interactions + Q-learning selects optimal subset
- **Performance**: Better F1 scores on DUD-E than random feature selection from co-crystals
- **Application**: Efficient solution for LIT-PCBA active molecule identification
- **Use Case**: Pharmacophore generation without ligand knowledge (apo structures)

**PharmacoNet**[^63][^64]

- **Technology**: Deep learning framework for structure-based pharmacophore modeling
- **Innovation**: Automated protein-based pharmacophore + parameterized analytical scoring
- **Speed**: Extremely fast (~100x faster than docking) for ultra-large library screening
- **Benchmark**: Reasonably accurate compared to traditional docking on DUD-E
- **Application**: Pre-filter for billion-compound libraries before docking


#### Ensemble Learning and Meta-Models

**Voting and Stacking**[^21]

- Combine multiple individual pharmacophore models
- Voting: Ligand accepted if matched by ≥N models (e.g., N=3 out of 5)
- Stacking: ML classifier (e.g., random forest) trained on multi-model match patterns
- **Performance**: AUC 0.994, EF1% 50.07 (apelin dataset)[^21]

**Consensus Docking + Pharmacophore Filtering**[^65][^66]

- Pre-docking pharmacophore filter removes unlikely binders
- Post-docking pharmacophore validation ensures key interactions
- **Example**: COVID-19 drug repurposing[^66][^65]
    - 6,218 drugs screened with pre/post pharmacophore filtering
    - Hit rate: 7/38 tested compounds active (18.4%)
    - 200-fold improvement over remdesivir for omipalisib

**Hybrid Approaches**[^40]

- Combine structure-based and ligand-based pharmacophores in parallel
- Compounds passing either model proceed to docking
- Complementarity increases overall diversity of hits
- **Performance**: Improved coverage of chemical space vs. single approach

***

### Practical Implementation Guidelines

#### For Academic Researchers

**Minimal Viable Workflow (Free Tools)**

1. **Pharmacophore Generation**: PharmaGist (web) or LigandScout Academic
2. **Virtual Screening**: Pharmit (web) or local AutoDock Vina
3. **Validation Dataset**: Download DUD-E or LIT-PCBA
4. **Metrics Calculation**: Rocker tool or custom Python scripts[^27]
5. **Visualization**: PyMOL or Chimera (free academic)

**Budget**: \$0-1,000 (for LigandScout Academic if chosen)
**Time to First Results**: 1-2 weeks

**Advanced Academic Workflow**

1. **Pharmacophore Generation**: Apo2ph4 (structure-based) + Phase/Discovery Studio (ligand-based)
2. **Consensus Modeling**: ConPhar for multi-complex integration
3. **Virtual Screening**: Glide or GOLD (if licenses available)
4. **Validation**: All four benchmark datasets (DUD-E, LIT-PCBA, MUV, DEKOIS 2.0)
5. **Advanced Metrics**: pROC-chemotype plots, BEDROC with multiple α values

**Budget**: \$10,000-30,000/year (software licenses)
**Time to Publication**: 3-6 months (including prospective validation)

#### For Pharmaceutical Industry

**High-Throughput Screening Integration**

1. **Pharmacophore Library**: Maintain consensus pharmacophores for priority targets
2. **AI-Accelerated Pre-Filtering**: PharmacoNet or custom ML models
3. **Parallel Screening**: Docking + shape + pharmacophore in parallel
4. **Hit Triage**: Multi-metric scoring (EF1%, BEDROC, GH, docking score)
5. **Automated Validation**: Continuous benchmarking against internal HTS data

**Recommended Stack**

- **Pharmacophore**: LigandScout + custom pipeline for consensus generation
- **Docking**: Glide + AutoDock Vina (for speed)
- **Shape**: ROCS + custom shape models
- **Infrastructure**: GPU clusters for PharmacoNet/ML models
- **Storage**: Database of validated pharmacophores (hundreds to thousands)

**Performance Targets**

- Screen 10-100 million compounds in 24-48 hours
- EF1% > 30 for validated pharmacophores
- Hit rate 5-15% in experimental validation
- Cost per validated hit: \$5,000-20,000 (including synthesis and assay)

***

### Critical Considerations and Pitfalls

#### Common Mistakes in Pharmacophore Validation

**1. Insufficient Decoy Quality**[^50][^4]

- **Problem**: Using random molecules or PubChem compounds as decoys
- **Consequence**: Artificially inflated enrichment from trivial physicochemical differences
- **Solution**: Use property-matched decoys (DeepCoy, DUD-E methodology)[^50]

**2. Dataset Contamination**[^5]

- **Problem**: Training and test sets share similar chemotypes or are from same source
- **Consequence**: Overestimated generalization performance
- **Solution**: Use external validation on independent dataset (e.g., LIT-PCBA after DUD-E training)[^5]

**3. Overlooking Chemotype Bias**[^51]

- **Problem**: Model performs well on one scaffold family but fails on others
- **Consequence**: Limited scaffold hopping capability, restricted hit diversity
- **Solution**: pROC-chemotype analysis to identify scaffold-specific performance[^51]

**4. Ignoring Protein Flexibility**[^67][^41]

- **Problem**: Single rigid structure used for pharmacophore generation
- **Consequence**: Misses alternative binding modes, reduced recall
- **Solution**: MD simulations or ensemble of structures (Apo2ph4, Flexi-pharma)[^67][^41]

**5. Tolerance Too Strict or Too Loose**[^38]

- **Problem**: 0.5 Å tolerance misses viable compounds; 3.0 Å allows nonsensical matches
- **Consequence**: Poor balance between sensitivity and specificity
- **Solution**: Optimize via cross-validation, typical sweet spot 1.0-2.0 Å[^38]


#### Benchmark-Specific Limitations

**DUD-E**[^1][^4]

- Analogue bias can inflate similarity-based method performance
- Neural networks exploit decoy artifacts rather than learning true physics[^4]
- **Mitigation**: Cross-validate on LIT-PCBA or MUV; analyze failure modes

**LIT-PCBA**[^5]

- Only 15 targets limits breadth of validation
- Some targets have limited active chemotype diversity
- **Mitigation**: Use for external validation after DUD-E training; supplement with target-specific experimental data

**MUV**[^9][^10]

- Low enrichment makes it difficult to distinguish good from moderate methods
- Some targets exceptionally challenging (near-random performance for all methods)
- **Mitigation**: Use as "stress test" for robustness; don't expect high enrichment

**DEKOIS 2.0**[^15]

- Smaller dataset sizes than DUD-E
- Not all pharmaceutical target families represented
- **Mitigation**: Preferred for structure-based methods; combine with DUD-E for comprehensive assessment

***

### Future Directions and Emerging Trends

#### Ultra-Large Library Screening

**Scale Requirements**[^64][^63][^6]

- Billion to multi-billion compound libraries (ZINC20, Enamine REAL)
- Traditional docking: months to years on HPC clusters
- **Solution**: AI-accelerated pre-filtering[^63][^64]
    - PharmacoNet: 100x faster than docking
    - Deep learning pose prediction (DiffDock-L): Orders of magnitude faster[^68]
    - Hybrid: Pharmacophore + ML scoring as gatekeeper

**Performance Expectations**[^69]

- Improved metric: Enhanced enrichment calculation for very large libraries[^26]
- Target: EF0.01% (enrichment in top 0.01% of billion-compound library)
- Challenge: Maintaining low false negative rate while achieving extreme speed


#### Integration with Generative AI

**Pharmacophore-Conditioned Generation**[^70]

- **Current State**: PharmacoForge generates pharmacophores; separate tools generate molecules
- **Future**: End-to-end generation of molecules satisfying pharmacophore constraints
- **Advantage**: Ensures synthesizability and drug-likeness from the start

**Reverse Pharmacophore Mapping**

- Identify optimal pharmacophore for desired multi-target profile
- AI proposes pharmacophore → generative model creates molecules → virtual screening validates
- **Application**: Polypharmacology design, drug repurposing


#### Enhanced Validation Paradigms

**Protein Dynamics Integration**[^67][^41]

- Standard: Static pharmacophores from single structures
- Future: Dynamic pharmacophores from extended MD simulations
- **Validation**: Time-averaged enrichment across conformational ensemble
- **Tools**: Dynophores, Flexi-pharma[^67]

**Multi-Objective Optimization**[^71]

- Beyond activity: Simultaneously optimize ADMET, selectivity, synthesizability
- Pharmacophores with multi-dimensional constraints
- **Metrics**: Pareto-optimal enrichment (activity vs. toxicity, etc.)

**Experimental Feedback Loops**

- Real-time model updating as HTS data become available
- Active learning: Pharmacophore suggests next compounds to test
- **Challenge**: Requires integration of computational and experimental workflows

***

### Conclusion and Recommendations

#### Key Takeaways

1. **Dataset Selection**
    - **Broad Benchmarking**: DUD-E (102 targets, historical standard)
    - **Unbiased Validation**: LIT-PCBA (15 targets, minimal bias)
    - **Challenging Test**: MUV (17 targets, eliminates trivial enrichment)
    - **SBVS Focus**: DEKOIS 2.0 (81 targets, structure-based emphasis)
2. **Metric Strategy**
    - **Primary**: EF1% (practical experimental validation), BEDROC α=160.9 (early emphasis)
    - **Secondary**: ROC AUC (overall ranking), GH score (precision-recall balance)
    - **Advanced**: pROC-chemotype plots (scaffold-specific performance)
3. **Parameter Optimization**
    - Feature tolerance: 1.0-1.5 Å for H-bonds, 1.5-2.0 Å for hydrophobic/aromatic
    - Consensus pharmacophores outperform single-molecule models by 2-3x in EF1%
    - Ensemble learning (voting/stacking) provides most robust performance
4. **Validation Protocol**
    - **Minimum**: Decoy set validation with property-matched decoys, calculate ROC AUC + EF1%
    - **Standard**: Fischer's randomization + test set + external dataset (>0.7 GH score, >0.8 AUC)
    - **Best Practice**: Prospective experimental validation with 5-20% expected hit rate
5. **Software Selection**
    - **Academic Budget**: PharmaGist (free web) + Pharmit + AutoDock Vina
    - **Academic Advanced**: LigandScout + Glide/GOLD + DeepCoy for decoys
    - **Industry**: LigandScout + Glide + ROCS + custom AI pipelines (PharmacoNet)

#### Downloadable Resources Summary

| Resource | URL | Content |
| :-- | :-- | :-- |
| **DUD-E** | http://dude.docking.org | 102 targets, actives + decoys, PDB structures |
| **LIT-PCBA** | http://drugdesign.unistra.fr/LIT-PCBA | 15 unbiased targets, training/validation sets |
| **DEKOIS 2.0** | http://www.dekois.com | 81 targets with cluster-based decoys |
| **MUV** | PubChem (via queries) | 17 PubChem bioassays, maximum unbiased |
| **Pharmit** | pharmit.csb.pitt.edu | Web-based pharmacophore search |
| **PharmaGist** | bioinfo3d.cs.tau.ac.il/PharmaGist | Ligand-based pharmacophore generation |
| **Apo2ph4** | GitHub (search "Apo2ph4") | Structure-based workflow |
| **DeepCoy** | GitHub (search "DeepCoy") | Deep learning decoy generator |

#### Final Recommendation

For validating a new pharmacophore model in 2025-2026, the optimal workflow combines:

1. **Training**: Consensus pharmacophore from 20-50 diverse actives (or 10+ co-crystal structures if structure-based)
2. **Internal Validation**: Fischer's randomization + decoy set (DUD-E methodology) → Target ROC AUC >0.8, EF1% >20, GH >0.7
3. **External Validation**: LIT-PCBA (if overlapping target) or independent literature actives → Expect 10-20% performance drop
4. **Prospective Validation**: Screen 50,000-1 million compounds, test top 50-100 → Expected 5-15% hit rate
5. **Publication**: Report all metrics (EF1%, EF5%, ROC AUC, BEDROC α=160.9, GH score) with raw data in supplementary materials

This comprehensive approach ensures reproducibility, minimizes bias, and maximizes translational potential to drug discovery applications.
<span style="display:none">[^100][^101][^102][^103][^104][^105][^106][^107][^108][^109][^110][^111][^112][^113][^114][^115][^116][^117][^118][^119][^120][^121][^122][^123][^124][^125][^126][^127][^128][^129][^130][^131][^132][^133][^134][^135][^136][^137][^138][^139][^140][^141][^142][^143][^144][^145][^146][^147][^148][^149][^150][^151][^152][^153][^154][^155][^156][^157][^158][^159][^160][^161][^162][^163][^164][^165][^166][^167][^168][^169][^170][^171][^172][^173][^174][^175][^176][^177][^178][^72][^73][^74][^75][^76][^77][^78][^79][^80][^81][^82][^83][^84][^85][^86][^87][^88][^89][^90][^91][^92][^93][^94][^95][^96][^97][^98][^99]</span>

<div align="center">⁂</div>

[^1]: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0167-x

[^2]: https://pmc.ncbi.nlm.nih.gov/articles/PMC5066283/

[^3]: https://pubs.acs.org/doi/10.1021/acs.jcim.0c00115

[^4]: https://dx.plos.org/10.1371/journal.pone.0220113

[^5]: https://pubs.acs.org/doi/10.1021/acs.jcim.0c00155

[^6]: https://www.frontiersin.org/articles/10.3389/fbinf.2025.1628800/full

[^7]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11469520/

[^8]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11687028/

[^9]: https://pubs.acs.org/doi/10.1021/ci8002649

[^10]: https://bmcchem.biomedcentral.com/articles/10.1186/1752-153X-3-S1-P17

[^11]: https://www.mdpi.com/1420-3049/21/4/476

[^12]: https://dx.plos.org/10.1371/journal.pone.0195478

[^13]: https://pubs.acs.org/doi/10.1021/acs.jcim.9b01210

[^14]: https://dx.plos.org/10.1371/journal.pone.0318712

[^15]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3980182/

[^16]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4450982/

[^17]: https://www.dovepress.com/benchmarking-the-structure-based-virtual-screening-performance-of-wild-peer-reviewed-fulltext-article-DDDT

[^18]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10105235/

[^19]: https://onlinelibrary.wiley.com/doi/10.1002/pro.4869

[^20]: https://www.tandfonline.com/doi/full/10.1080/07391102.2023.2256876

[^21]: https://www.frontiersin.org/articles/10.3389/fchem.2024.1382319/full

[^22]: https://academic.oup.com/bib/article/26/Supplement_1/i37/8378006

[^23]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11058650/

[^24]: https://pubs.acs.org/doi/10.1021/acs.jpcb.4c01875

[^25]: https://pubs.acs.org/doi/10.1021/acs.jcim.9b00380

[^26]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10980085/

[^27]: https://pmc.ncbi.nlm.nih.gov/articles/PMC5013620/

[^28]: https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/slct.202304077

[^29]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6819673/

[^30]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6800073/

[^31]: https://www.tandfonline.com/doi/full/10.1080/07391102.2023.2238062

[^32]: https://www.tandfonline.com/doi/full/10.1080/07391102.2022.2027819

[^33]: https://www.frontiersin.org/articles/10.3389/fphar.2022.837369/full

[^34]: https://www.tandfonline.com/doi/full/10.1080/1062936X.2020.1827030

[^35]: https://www.mdpi.com/1420-3049/24/17/3181/pdf

[^36]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6175823/

[^37]: http://journal.frontiersin.org/article/10.3389/fchem.2017.00007/full

[^38]: https://journals.innovareacademics.in/index.php/ajpcr/article/view/53258

[^39]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9832483/

[^40]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8999148/

[^41]: https://www.mdpi.com/1422-0067/20/23/5834

[^42]: https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.3c01439

[^43]: https://app.jove.com/t/68933/pharmacophore-modeling-for-targets-with-extensive-ligand-libraries

[^44]: https://www.mdpi.com/1420-3049/25/2/385/pdf

[^45]: https://virologyj.biomedcentral.com/articles/10.1186/s12985-024-02607-4

[^46]: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00537-9

[^47]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4077067/

[^48]: https://pmc.ncbi.nlm.nih.gov/articles/PMC5408157/

[^49]: https://japtronline.com/index.php/joapr/article/view/951

[^50]: https://ora.ox.ac.uk/objects/uuid:bc03af13-709f-445a-952e-ba6a08537baf/download_file?file_format=pdf\&safe_filename=Imrie_et_al_2021_generating_property_matched.pdf\&type_of_work=Journal+article

[^51]: https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-6-S1-P19

[^52]: https://www.eurekaselect.com/237887/article

[^53]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3980120/

[^54]: https://dx.plos.org/10.1371/journal.pone.0316765

[^55]: https://www.ijsr.net/getabstract.php?paperid=MR22213165325

[^56]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2447755/

[^57]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4012150/

[^58]: https://pubs.acs.org/doi/10.1021/ci8004226

[^59]: https://xlink.rsc.org/?DOI=D0RA02363A

[^60]: http://thescipub.com/abstract/10.3844/amjsp.2010.151.156

[^61]: https://www.mdpi.com/1422-0067/23/14/7747/pdf?version=1657724205

[^62]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6769597/

[^63]: https://arxiv.org/pdf/2310.00681.pdf

[^64]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11575537/

[^65]: https://www.pnas.org/content/pnas/118/30/e2024302118.full.pdf

[^66]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8325362/

[^67]: https://pmc.ncbi.nlm.nih.gov/articles/PMC7449997/

[^68]: https://pubs.acs.org/doi/10.1021/acs.jcim.5c00380

[^69]: http://arxiv.org/pdf/2403.10478.pdf

[^70]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10558534/

[^71]: https://www.mdpi.com/1422-0067/20/9/2060/pdf

[^72]: http://biorxiv.org/lookup/doi/10.1101/2025.02.17.638675

[^73]: https://arxiv.org/abs/2509.16273

[^74]: http://biorxiv.org/lookup/doi/10.1101/2025.02.17.638750

[^75]: https://arxiv.org/abs/2409.06316

[^76]: https://www.semanticscholar.org/paper/f94adf34c933d28a4892b5f146af4ac516fef914

[^77]: https://pmc.ncbi.nlm.nih.gov/articles/PMC7352161/

[^78]: http://biorxiv.org/lookup/doi/10.1101/2025.04.18.649577

[^79]: https://www.tandfonline.com/doi/full/10.1080/14756366.2019.1693702

[^80]: https://link.springer.com/10.1007/s11696-020-01274-3

[^81]: https://link.springer.com/10.1007/s10822-020-00339-5

[^82]: https://www.semanticscholar.org/paper/4d3f1aee8ba07bc30b4b60d185bea2586e699685

[^83]: https://pubs.acs.org/doi/10.1021/ci700024q

[^84]: https://www.frontiersin.org/articles/10.3389/fchem.2024.1382319/pdf?isPublishedV2=False

[^85]: https://www.mdpi.com/1424-8247/15/5/646/pdf?version=1653378280

[^86]: https://downloads.hindawi.com/journals/mpe/2021/5585990.pdf

[^87]: https://www.mdpi.com/1424-8247/15/9/1122/pdf?version=1662631835

[^88]: https://www.mdpi.com/2218-0532/76/3/333/pdf?version=1482888651

[^89]: https://www.semanticscholar.org/paper/186856e9ad867172401d2115e6d733b7a478965c

[^90]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8526798/

[^91]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11377542/

[^92]: https://www.mdpi.com/1422-0067/20/6/1375/pdf

[^93]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4278665/

[^94]: https://www.mdpi.com/1420-3049/20/10/18107

[^95]: https://www.semanticscholar.org/paper/70b12b95a75b246ea451210447dec8cdaa9e10cc

[^96]: https://bmcchem.biomedcentral.com/articles/10.1186/1752-153X-2-S1-P1

[^97]: https://www.semanticscholar.org/paper/bc2bab6c46214c7db3f990f12d0480bf561a1557

[^98]: http://arxiv.org/pdf/2410.12722.pdf

[^99]: http://arxiv.org/pdf/2211.13235.pdf

[^100]: https://arxiv.org/html/2504.00487v2

[^101]: https://arxiv.org/pdf/2411.09547.pdf

[^102]: https://arxiv.org/html/2406.09864

[^103]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9825611/

[^104]: https://arxiv.org/html/2502.05922v1

[^105]: https://arxiv.org/pdf/1807.05786.pdf

[^106]: https://www.tandfonline.com/doi/full/10.1080/07391102.2017.1372313

[^107]: https://www.rjptonline.org/AbstractView.aspx?PID=2025-18-8-61

[^108]: https://linkinghub.elsevier.com/retrieve/pii/S2589750024000657

[^109]: http://link.springer.com/10.1007/s00044-014-1057-2

[^110]: https://link.springer.com/10.1007/s11224-022-01911-5

[^111]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9068899/

[^112]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4938243/

[^113]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8351372/

[^114]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9145410/

[^115]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4048638/

[^116]: https://www.frontiersin.org/articles/10.3389/fcimb.2021.611304/pdf

[^117]: https://www.tandfonline.com/doi/full/10.1080/07391102.2021.1993341

[^118]: https://link.springer.com/10.1007/s11030-020-10148-5

[^119]: https://www.tandfonline.com/doi/full/10.1080/07391102.2024.2427375

[^120]: https://link.springer.com/10.1007/s00249-024-01713-z

[^121]: https://pubs.acs.org/doi/10.1021/acsptsci.0c00215

[^122]: https://pubs.acs.org/doi/10.1021/acs.jcim.2c01641

[^123]: https://www.tandfonline.com/doi/full/10.1080/07391102.2022.2112080

[^124]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8107043/

[^125]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11662536/

[^126]: https://pmc.ncbi.nlm.nih.gov/articles/PMC7597227/

[^127]: https://f1000research.com/articles/9-514/v2/pdf

[^128]: https://f1000research.com/articles/9-514/v1/pdf

[^129]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9088954/

[^130]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2699263/

[^131]: https://figshare.com/articles/journal_contribution/Development_of_pharmacophore_similarity_based_quantitative_activity_hypothesis_and_its_applicability_domain_applied_on_a_diverse_data_set_of_HIV_1_integrase_inhibitors/999125/3/files/1871887.pdf

[^132]: https://www.mdpi.com/1420-3049/23/8/1959/pdf

[^133]: https://peerj.com/articles/725

[^134]: https://www.nature.com/articles/s42256-025-00993-0

[^135]: http://biorxiv.org/lookup/doi/10.1101/2024.03.14.585019

[^136]: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-025-01039-8

[^137]: https://www.semanticscholar.org/paper/973d59d2bbe2d735e0ed2620cad834e8106d08ac

[^138]: https://ojs.aaai.org/index.php/AAAI/article/view/25970

[^139]: https://arxiv.org/pdf/2403.10478v1.pdf

[^140]: https://www.frontiersin.org/articles/10.3389/fphar.2018.00011/pdf

[^141]: https://onlinelibrary.wiley.com/doi/10.1002/jcc.23690

[^142]: http://biorxiv.org/lookup/doi/10.1101/2024.02.22.581511

[^143]: https://pubs.acs.org/doi/10.1021/ci200562p

[^144]: https://www.nature.com/articles/s41401-025-01644-1

[^145]: https://link.springer.com/10.1007/s00044-025-03455-9

[^146]: https://onlinelibrary.wiley.com/doi/10.1002/minf.201000010

[^147]: https://arxiv.org/pdf/2306.02851.pdf

[^148]: https://www.mdpi.com/1424-8220/21/14/4769/pdf

[^149]: https://www.jmir.org/2023/1/e45044/PDF

[^150]: https://arxiv.org/pdf/1905.03702.pdf

[^151]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9670537/

[^152]: https://link.springer.com/10.1007/s40203-025-00414-5

[^153]: https://journal.neark.kz/wp-content/uploads/pdf/2025-4/07-Автоматизированная оценка читаемости.pdf

[^154]: https://brieflands.com/journals/ijpr/articles/164183

[^155]: https://journals.lww.com/10.1097/JS9.0000000000003269

[^156]: https://www.eurekaselect.com/249088/article

[^157]: http://www.tandfonline.com/doi/abs/10.1080/08927022.2012.751592

[^158]: https://www.mdpi.com/2218-273X/15/6/788

[^159]: https://ieeexplore.ieee.org/document/10947919/

[^160]: https://bmcchem.biomedcentral.com/articles/10.1186/s13065-020-00704-3

[^161]: https://www.mdpi.com/2079-6382/12/4/630/pdf?version=1679550092

[^162]: http://arxiv.org/pdf/2002.09316.pdf

[^163]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3152796/

[^164]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10135334/

[^165]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3475860/

[^166]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10864924/

[^167]: https://onlinelibrary.wiley.com/doi/10.1002/minf.202200034

[^168]: https://arxiv.org/abs/2506.05768

[^169]: https://innovareacademics.in/journals/index.php/ijap/article/view/31484

[^170]: https://pubs.aip.org/aip/acp/article/822579

[^171]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9049841/

[^172]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6540224/

[^173]: https://peerj.com/preprints/2721v1

[^174]: https://pubs.acs.org/doi/10.1021/acs.jcim.8b00185

[^175]: https://journals.asm.org/doi/10.1128/msystems.00885-25

[^176]: https://pubs.acs.org/doi/10.1021/ci500681r

[^177]: https://bio-protocol.org/en/bpdetail?id=5523\&type=0

[^178]: https://www.frontiersin.org/articles/10.3389/fmicb.2024.1346442/full


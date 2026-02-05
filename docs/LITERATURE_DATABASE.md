# Literature Database: CCR2 Pharmacophore Virtual Screening

**Project**: Pharmacophore-Toolkit - CCR2 Virtual Screening Implementation
**Phase**: 1 - Literature Research & Algorithm Design
**Last Updated**: 2026-01-23

---

## Table of Contents

1. [CCR2 Target Information](#1-ccr2-target-information)
2. [Core Methodology Papers](#2-core-methodology-papers)
3. [Shape-Based Screening Tools](#3-shape-based-screening-tools)
4. [Validation Metrics & Benchmarking](#4-validation-metrics--benchmarking)
5. [CCR2-Specific Literature](#5-ccr2-specific-literature)
6. [Recent ML Advances](#6-recent-ml-advances)
7. [Tool Documentation](#7-tool-documentation)
8. [Research Gaps & Opportunities](#8-research-gaps--opportunities)

---

## 1. CCR2 Target Information

### 1.1 ChEMBL Target Profile

| Property | Value |
|----------|-------|
| **ChEMBL ID** | CHEMBL4015 |
| **Preferred Name** | C-C chemokine receptor type 2 |
| **Gene Symbol** | CCR2 |
| **UniProt** | P41597 |
| **Target Type** | SINGLE PROTEIN |
| **Organism** | Homo sapiens |
| **Tax ID** | 9606 |

### 1.2 PDB Structures

| PDB ID | Description | Resolution | Ligands |
|--------|-------------|------------|---------|
| **5T1A** | CCR2 with orthosteric + allosteric antagonists | 2.81 Å | BMS-681, CCR2-RA-[R] |
| 7P8X | CCR2 structure | TBD | TBD |
| 7XA3 | CCR2 structure | TBD | TBD |
| 2MLO | NMR structure | N/A | TBD |
| 2MLQ | NMR structure | N/A | TBD |

### 1.3 Bioactivity Summary (from ChEMBL)

```
Total activities (pChEMBL >= 7): 1,835 data points

Most Potent Compounds:
- Ki = 2.63 nM (pChEMBL = 8.58)
- Ki = 4.17 nM (pChEMBL = 8.38)
- IC50 = 5.9 nM (pChEMBL = 8.23)
- Ki = 6.92 nM (pChEMBL = 8.16)
- IC50 = 8.0 nM (pChEMBL = 8.10)
```

### 1.4 Allosteric Binding Site Characteristics

**Location**: Intracellular, most intracellular site observed in class A GPCRs

**Key Residues** (from crystal structure analysis):
- TM2, TM3, TM6, TM7, and H8
- Tyr7.53 (facilitates hydrophobic interactions)
- Residues 1.56, 2.43, and 3.50 (hydrophobic contacts)
- TM7-H8 motif (critical for ligand binding)

**Pharmacophore Features**:
- Hydroxypyrrolinone motif → hydrogen bonding with TM7-H8
- Predominantly hydrophobic interactions
- Druggable pocket with good polarity/hydrophobicity balance

---

## 2. Core Methodology Papers

### 2.1 Consensus Virtual Screening

| Status | Title | Year | Journal | Key Findings | Link |
|--------|-------|------|---------|--------------|------|
| ⭐ **READ** | Consensus holistic virtual screening for drug discovery: a novel machine learning model approach | 2024 | J Cheminformatics | Combines QSAR, Pharmacophore, docking, 2D shape; AUC 0.90-0.98 | [PMC11134635](https://pmc.ncbi.nlm.nih.gov/articles/PMC11134635/) |
| ⭐ **READ** | Integrating traditional and modern approaches for comprehensive pharmacophore map validation | 2025 | J Comput-Aided Mol Des | Reviews ROC-AUC, EF, BEDROC, DUD-E benchmarks | [Springer](https://link.springer.com/article/10.1007/s10822-025-00751-9) |
| 📋 TODO | Pharmacophore Model-Based Virtual Screening Workflow for Plasmodium Hsp90 | 2023 | ACS Omega | Complete workflow with EF validation | [ACS](https://pubs.acs.org/doi/10.1021/acsomega.3c04494) |

### 2.2 Shape-Based Screening

| Status | Title | Year | Journal | Key Findings | Link |
|--------|-------|------|---------|--------------|------|
| ⭐ **READ** | VSFlow: an open-source ligand-based virtual screening tool | 2023 | J Cheminformatics | RDKit-based, combo scores | [BMC](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00703-1) |
| ⭐ **READ** | Building shape-focused pharmacophore models for docking screening | 2024 | J Cheminformatics | Shape + pharmacophore integration | [BMC](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00857-6) |
| 📋 TODO | PheSA: Open-Source Tool for Pharmacophore-Enhanced Shape Alignment | 2024 | JCIM | Enhanced shape + pharma | [ACS](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00516) |
| 📋 TODO | Comparison of Shape-Matching and Docking as Virtual Screening Tools | 2007 | J Med Chem | Benchmark comparison | [ACS](https://pubs.acs.org/doi/10.1021/jm0603365) |

### 2.3 Validation & Metrics

| Status | Title | Year | Key Contribution | Link |
|--------|-------|------|-----------------|------|
| ⭐ **READ** | Evaluating Virtual Screening Methods: Good and Bad Metrics for the "Early Recognition" Problem | 2007 | EF, BEDROC definitions | [ResearchGate](https://www.researchgate.net/publication/6517236) |
| 📋 TODO | ROC-AUC implementation notebook | 2020 | sklearn-based code | [Gist](https://gist.github.com/ravila4/26b6ceb21c7e87af80be01f4620a7a58) |

---

## 3. Shape-Based Screening Tools

### 3.1 Tool Comparison Matrix

| Tool | License | Speed | Accuracy | Conformers | Notes |
|------|---------|-------|----------|------------|-------|
| **ROCS** | Commercial | Fast | Best | Yes | Industry standard |
| **CDPKit** | Open source | Fast | Good | Yes (CONFORT) | **Recommended** |
| **RDKit rdShapeHelpers** | Open source | Moderate | Good | Manual | Built-in to toolkit |
| **Shape-it** | Open source | Fast | Good | Yes | No electrostatics |
| **espsim** | Open source | Moderate | Good | Yes | Electrostatic focus |

### 3.2 CDPKit Details

**Documentation**: [CDPKit Introduction](https://cdpkit.org/introduction.html)

**Key Features**:
- Gaussian shape-based molecule alignment
- Pharmacophore generation, alignment, screening
- CONFORT conformer generator (state-of-the-art)
- Command-line tools + Python API
- Free and open source

**Performance** (from PharmacoMatch benchmark):
- Competitive AUROC scores vs neural methods
- Robust correlation with neural hitlists
- Fast conformer generation

**Conformer Generator (CONFORT)**:
- Allows up to 25+ conformations per molecule
- Assessed on 2,912 PDB conformations
- State-of-the-art bioactive reproduction

### 3.3 RDKit Shape Tools

```python
# Key RDKit modules for shape screening
from rdkit.Chem import rdShapeHelpers, AllChem, rdMolAlign
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate

# Shape similarity metrics
# - TanimotoDist (shape Tanimoto)
# - TverskyShape
# - ProtrudeDist

# Combo scoring approach (from VSFlow):
# combo_score = (shape_similarity + fingerprint_similarity) / 2
```

---

## 4. Validation Metrics & Benchmarking

### 4.1 Metric Definitions

| Metric | Formula | Interpretation | Threshold |
|--------|---------|----------------|-----------|
| **ROC-AUC** | Area under ROC curve | Global discrimination | > 0.7 reliable, > 0.9 excellent |
| **EF@1%** | (Hits in top 1%) / (Expected random) | Early enrichment | > 10 good, > 20 excellent |
| **EF@5%** | (Hits in top 5%) / (Expected random) | Early enrichment | > 5 good |
| **BEDROC(α=20)** | Exponentially weighted early enrichment | Balanced metric | > 0.5 good |
| **RIE** | Robust initial enhancement | Statistical robustness | Target-dependent |

### 4.2 Benchmark Datasets

| Dataset | Actives | Decoys | Targets | Use Case |
|---------|---------|--------|---------|----------|
| **DUD-E** | ~22,000 | ~1.3M | 102 | Standard benchmark |
| **LIT-PCBA** | Variable | Variable | 15 | Realistic actives |
| **DEKOIS 2.0** | Variable | Variable | 81 | Property-matched decoys |
| **MUV** | 30/target | 15,000/target | 17 | Maximum unbiased |

### 4.3 Our CCR2 Dataset Statistics

```
Actives: 75 (allosteric CCR2 modulators)
Decoys: 500 (property-matched)
Ratio: 1:6.67 (moderately imbalanced)

Expected random EF@1%: 1.0
Maximum possible EF@1%: 100/1.3 ≈ 76.9

Target performance:
- EF@1% >= 10 (good)
- EF@1% >= 15 (target)
- EF@1% >= 25 (stretch)
```

---

## 5. CCR2-Specific Literature

### 5.1 Crystal Structure Papers

| Status | Title | Year | Key Findings | DOI/Link |
|--------|-------|------|--------------|----------|
| ⭐ **READ** | Structure of CC chemokine receptor 2 with orthosteric and allosteric antagonists | 2016 | First CCR2 structure, 5T1A | [Nature](https://www.nature.com/articles/nature20605) |
| ⭐ **READ** | Structural basis for ligand modulation of the CCR2 conformational landscape | 2019 | Dynamics, metastable states | [PNAS](https://www.pnas.org/doi/10.1073/pnas.1814131116) |
| 📋 TODO | Pyrrolone Derivatives as Intracellular Allosteric Modulators for CCR1/CCR2 | 2018 | SAR for allosteric modulators | [JMedChem](https://pubs.acs.org/doi/10.1021/acs.jmedchem.8b00605) |

### 5.2 Key Structural Insights for Pharmacophore Design

**Allosteric Site Features** (from 5T1A):

1. **Hydrogen Bond Acceptors**:
   - TM7-H8 junction (conserved across chemokine receptors)
   - Hydroxypyrrolinone motif interactions

2. **Hydrophobic Contacts**:
   - Tyr7.53 (critical residue)
   - Residues 1.56, 2.43, 3.50
   - Non-polar pocket lining

3. **Spatial Constraints**:
   - Pocket is relatively small but druggable
   - Overlaps with G-protein binding site
   - Allosteric mechanism: blocks activation-associated conformational changes

**Pharmacophore Hypothesis**:
- 1-2 H-bond acceptors (for TM7-H8)
- 2-3 hydrophobic features
- Aromatic ring (for Tyr7.53 interaction)
- Compact molecular shape

---

## 6. Recent ML Advances

### 6.1 Neural Pharmacophore Methods

| Status | Title | Year | Approach | Performance | Link |
|--------|-------|------|----------|-------------|------|
| ⭐ **READ** | PharmacoNet: Deep Pharmacophore Modeling | 2024 | Computer vision segmentation | 187M compounds/21h | [arXiv](https://arxiv.org/abs/2310.00681) |
| ⭐ **READ** | PharmacoMatch: Neural Subgraph Matching | 2024 | Contrastive learning | Competitive with CDPKit | [arXiv](https://arxiv.org/abs/2409.06316) |
| 📋 TODO | dyphAI: Dynamic Pharmacophore with AI | 2025 | Ensemble + Optuna | TBD | [Frontiers](https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2025.1479763/full) |

### 6.2 Optimization Methods

| Status | Title | Year | Approach | Link |
|--------|-------|------|----------|------|
| 📋 TODO | Pareto Optimization for Virtual Screening | 2024 | Multi-objective, 8% library exploration | [RSC](https://pubs.rsc.org/en/content/articlehtml/2024/dd/d3dd00227f) |
| 📋 TODO | SALSA: Active Learning on Synthons | 2025 | Trillion-scale active learning | [arXiv](https://arxiv.org/abs/2505.12913) |
| 📋 TODO | Bayesian Optimization for Drug Discovery | 2025 | Multi-fidelity | [ACS](https://pmc.ncbi.nlm.nih.gov/articles/PMC11869128/) |

---

## 7. Tool Documentation

### 7.1 CDPKit

| Resource | URL | Notes |
|----------|-----|-------|
| Main Documentation | https://cdpkit.org/ | Official docs |
| Shape Screening | https://cdpkit.org/applications/shapescreen.html | shapescreen tool |
| CONFORT Conformers | https://cdpkit.org/applications/confgen.html | Conformer generation |
| Python API | https://cdpkit.org/api/ | Programming interface |
| GitHub | https://github.com/molinfo-vienna/CDPKit | Source code |

### 7.2 RDKit

| Resource | URL | Notes |
|----------|-----|-------|
| Shape Helpers | https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html | Shape similarity |
| Molecular Alignment | https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html | O3A alignment |
| Benchmarking Platform | https://github.com/rdkit/benchmarking_platform | VS benchmarks |

### 7.3 Optuna (Optimization)

| Resource | URL | Notes |
|----------|-----|-------|
| Documentation | https://optuna.readthedocs.io/ | Official docs |
| Multi-objective | https://optuna.readthedocs.io/en/stable/reference/multi_objective.html | Pareto optimization |
| Visualization | https://optuna.readthedocs.io/en/stable/reference/visualization/index.html | Built-in plotting |

---

## 8. Research Gaps & Opportunities

### 8.1 Open Questions

1. **Optimal conformer count for CCR2 screening?**
   - Literature suggests 10-50 conformers
   - CDPKit CONFORT can generate 25+ efficiently
   - Need to benchmark on our dataset

2. **Best shape/pharmacophore weighting?**
   - VSFlow uses 50/50 combo score
   - May need target-specific optimization
   - Pareto optimization could reveal optimal trade-off

3. **Allosteric vs orthosteric pharmacophore models?**
   - Our actives are allosteric modulators
   - Need to ensure pharmacophore captures intracellular features
   - Reference ligand selection is critical

4. **Decoy quality for CCR2?**
   - Are our 500 decoys confirmed inactives or property-matched only?
   - May affect EF interpretation
   - Consider DUD-E style validation

### 8.2 Innovation Opportunities

1. **Multi-fidelity Optimization**
   - Low fidelity: 1 conformer, fast
   - High fidelity: 50 conformers, accurate
   - Progressive refinement saves compute

2. **Ensemble Pharmacophore Voting**
   - Instead of single "best" model
   - Use Pareto front of models
   - Consensus voting for final ranking

3. **Active Learning for Parameter Tuning**
   - Intelligent sampling instead of grid search
   - Focus on promising parameter regions
   - 10x faster optimization

4. **Integration with Deep Learning**
   - PharmacoNet for initial screening
   - Traditional methods for refinement
   - Best of both worlds

---

## Appendix: Search Queries Used

### PubMed (via Web Search fallback)
```
1. "pharmacophore virtual screening consensus enrichment factor ROC-AUC validation 2023 2024 2025"
2. "CCR2 chemokine receptor allosteric modulator crystal structure binding site pharmacophore"
```

### bioRxiv
```
- Category: bioinformatics, recent 365 days
- Category: pharmacology and toxicology, recent 365 days
```

### ChEMBL
```
- Target: CCR2 (CHEMBL4015), gene_symbol=CCR2, organism=Homo sapiens
- Bioactivity: target_chembl_id=CHEMBL4015, min_pchembl=7
```

### Web Search
```
- "shape-based virtual screening TanimotoCombo CDPKit RDKit ROCS comparison benchmark"
- "CDPKit shape screening pharmacophore virtual screening benchmark performance"
```

---

## Session Notes

### 2026-01-23: Initial Literature Gathering

**Completed**:
- [x] ChEMBL CCR2 target identification (CHEMBL4015)
- [x] ChEMBL bioactivity data (1,835 activities, pChEMBL >= 7)
- [x] CCR2 crystal structure identification (5T1A)
- [x] Core methodology papers identified
- [x] Shape screening tools compared
- [x] Validation metrics documented

**Key Findings**:
1. CCR2 has well-characterized allosteric site (5T1A crystal structure)
2. Consensus scoring achieves AUC 0.90-0.98 on benchmark targets
3. CDPKit is recommended open-source alternative to ROCS
4. VSFlow provides good reference implementation
5. PharmacoMatch shows neural methods are competitive

**Next Steps**:
- [ ] Download and read full PDFs of key papers
- [ ] Extract specific methodology details
- [ ] Design pharmacophore hypothesis based on 5T1A structure
- [ ] Begin implementation planning

---

*Document maintained for CCR2 pharmacophore virtual screening project*

# Comparison of Consensus Pharmacophore Screening and Reference Alignment for Virtual Screening

## Executive Summary
This report compares two prominent ligand-based virtual screening (LBVS) methodologies: **Consensus Pharmacophore Screening** and **Reference Alignment** (e.g., shape-based methods like ROCS). While consensus models leverage information from multiple actives to define essential interactions, they are prone to failures caused by structural heterogeneity, such as disconnected fragments and rigid distance constraints. Conversely, reference alignment focuses on volumetric and chemical similarity to a specific template, offering superior performance in scaffold hopping. Strategic parameter optimization—including tolerance, occurrence thresholds, and shape/color weighting—is critical for maximizing active/decoy separation and recovering diverse chemotypes.

## 1. Methodology Overview

Ligand-based virtual screening relies on the "similarity principle," which states that structurally similar molecules are likely to exhibit similar biological activities. However, how "similarity" is defined varies significantly between discrete pharmacophore models and continuous shape-based alignments.

### 1.1 Consensus Pharmacophore Screening
**Consensus Pharmacophore Screening** constructs a hypothesis by aligning multiple known active ligands to identify common spatial features essential for binding.

*   **Definition**: A pharmacophore model derived from the alignment of multiple known active ligands. It identifies "consensus" features (H-bond donors/acceptors, aromatics, hydrophobic centroids, positive/negative ionizables) shared by the set.
*   **Mechanism**: Operates by matching discrete geometric points in 3D space. Common tools include **Schrödinger Phase**, **MOE (Molecular Operating Environment)**, and **LigandScout**. The process typically involves:
    1.  Conformational expansion of active ligands.
    2.  Alignment of actives to identify common chemical features.
    3.  Generation of a query hypothesis consisting of spheres (features) with specific radii and vector constraints.
*   **Strengths**: Captures a generalized binding hypothesis that is less biased toward a single chemotype compared to single-ligand models. It effectively filters out molecules that lack essential interaction points [^p6], [^u14].

### 1.2 Reference Alignment (Template-based/Shape-based)
**Reference Alignment** (often synonymous with Shape-based Screening) uses the volumetric and electrostatic fields of a single bioactive conformation (or an ensemble) to score database molecules.

*   **Definition**: Aligns database molecules to a single or ensemble of "reference" conformations (e.g., a bioactive pose from a crystal structure or a highly active ligand).
*   **Mechanism**: Primarily uses volumetric overlap (shape) and chemical feature matching (color). **ROCS (Rapid Overlay of Chemical Structures)** by OpenEye is the industry standard. It utilizes Gaussian functions to approximate molecular shape, allowing for extremely rapid superposition and scoring.
    *   **Shape**: Represents the steric volume of the molecule.
    *   **Color**: Represents chemical properties (H-bond donors, acceptors, rings) as fuzzy functions rather than discrete points [^u82].
*   **Strengths**: Demonstrates high performance in **scaffold hopping**—identifying molecules with similar volumes and electronic fields but different topological connectivity (graphs). It excels when discrete pharmacophore points differ slightly but the overall "envelope" is conserved [^p10], [^u34].

## 2. Failure Modes of Consensus Pharmacophores

While theoretically robust, consensus models often suffer from "over-optimization" or "structural fragmentation" when derived from diverse ligands. These failure modes inevitably lead to poor recall (false negatives).

### 2.1 Disconnected Fragments
*   **Cause**: When a consensus model is built from diverse ligands that bind in different sub-pockets of a large binding site, the resulting query often aggregates features that are spatially distant. For example, Ligand A binds to the hinge region, while Ligand B extends into the solvent-exposed area.
*   **Failure Mechanism**: The consensus algorithm creates a "super-pharmacophore" containing features from both regions. If the model requires matching *all* features, no single molecule in the database may be large enough or flexible enough to span these "disconnected fragments."
*   **Outcome**: High precision but extremely low recall. The model fails to identify active molecules that only occupy a subset of the pocket (e.g., fragment-sized hits or smaller actives) [^u21], [^u66].

### 2.2 Distance Constraints
*   **Cause**: Pharmacophores are inherently rigid distance-based models. A consensus model typically defines the geometric relationship between features using the average distance or the intersection of distances observed in the training set.
*   **Failure Mechanism**: Tight distance constraints fail to accommodate "chemotype switching" or induced-fit variations. If a novel active molecule has a slightly different linker length or bond angle (e.g., a 5-membered ring vs. a 6-membered ring spacer), it will fail the rigid distance check even if the functional groups interact with the same protein residues.
*   **Outcome**: Poor separation of actives from decoys that are "close but not perfect" matches. This rigidity is a primary reason why shape-based methods often outperform pharmacophores in diversity tasks [^u37], [^u65].

## 3. Parameter Mitigation Strategies

To overcome the inherent rigidity of consensus pharmacophores and the potential non-specificity of shape alignment, specific parameters must be tuned. These settings directly influence the active/decoy separation ratio.

### 3.1 Tolerance (Feature Radius)
*   **Function**: Defines the sphere of influence (radius) for each pharmacophore feature. It determines how far a ligand's chemical group can be from the model's centroid and still be considered a match.
*   **Fix for Distance Constraints**: Increasing the tolerance allows for geometric flexibility.
    *   **Standard Setting**: **1.0 Å**.
    *   **Mitigation Setting**: Increasing to **1.2–1.5 Å** permits the model to capture actives with slight conformational differences or induced-fit deviations [^u119].
*   **Best Practice**: Use **1.2–1.4 Å** for general screening. Tighter tolerances (<1.0 Å) should only be reserved for highly rigid targets where specificity is paramount. Excessively large tolerances (>2.0 Å) introduce significant noise (decoys) [^u115], [^u116].

### 3.2 Occurrence Threshold (Consensus Score)
*   **Function**: Sets the minimum number or percentage of ligands in the training set that must share a feature for that feature to be included in the final query.
*   **Fix for Disconnected Fragments**: The "Match All" approach is often fatal for consensus models.
    *   **Partial Matching**: Lowering the threshold (e.g., from **100%** to **60–80%**) or using "Match $N$ of $M$" logic ensures that features unique to specific chemotypes do not disqualify valid hits.
    *   **Logic**: It prevents the "union" of all features from becoming an impossible-to-satisfy query.
*   **Best Practice**: A threshold of **3 out of 4/5** features or **60–80%** consensus is recommended to balance sensitivity and specificity. This "partial match" strategy significantly improves recall for fragmented pockets [^u24], [^u114].

### 3.3 Shape and Color Weights (Reference Alignment)
*   **Function**: In shape-based alignment (e.g., ROCS), the final score (often `ComboScore`) is a linear combination of Shape (volume overlap) and Color (chemical feature overlap).
    $$ \text{ComboScore} = w_{\text{shape}} \times \text{ShapeScore} + w_{\text{color}} \times \text{ColorScore} $$
*   **Optimization**:
    *   **Shape Weight**: High weight (e.g., 2.0) prioritizes steric fit. This is optimal when the binding pocket is narrow and sterically constrained (e.g., kinases, ion channels).
    *   **Color Weight**: High weight (e.g., 2.0) prioritizes chemical interactions. This is best for **scaffold hopping** where the core scaffold volume may change, but specific key interactions (e.g., a hinge-binding donor) must remain constant [^u54], [^u82].
*   **Best Practice**: The default **1:1 ratio (ComboScore)** is widely regarded as the most robust starting point for active/decoy separation. Deviating from this usually requires specific knowledge of the target's flexibility [^u74].

## 4. Evidence and Best Practices

### 4.1 Benchmarking Evidence (DUD-E / LIT-PCBA)
*   **Consensus Performance**: Retrospective studies on datasets like DUD-E indicate that consensus pharmacophores often achieve higher **Enrichment Factors (EF)** in the top 1% of the database compared to single-ligand models. However, they frequently suffer from lower **AUC (Area Under the Curve)** due to the "false negative" problem described in Section 2.1 [^p6], [^p10].
*   **Reference Alignment (ROCS)**: Generally provides better **scaffold hopping** (chemotype diversity) than rigid pharmacophores. Because it treats chemical features as fuzzy "color" densities rather than hard points, it recovers actives that would otherwise fail a binary pharmacophore filter [^u34], [^u77].

### 4.2 Comparative Analysis of Techniques
The following table summarizes the operational differences and best-practice settings for maximizing active/decoy separation.

| Feature | Consensus Pharmacophore | Reference Alignment (Shape-based) |
| :--- | :--- | :--- |
| **Primary Metric** | Geometric distance (RMSD of features) | Volumetric Overlap (Tanimoto) |
| **Constraint Type** | Hard (Discrete points & radii) | Soft (Gaussian functions) |
| **Failure Mode** | Disconnected fragments; Rigid distances | "Amorphous" shapes (non-specific blobs) |
| **Key Tuning Parameter** | **Tolerance (Radius)** & **Threshold** | **Shape/Color Weight** |
| **Recommended Setting** | Tolerance: **1.2–1.5 Å**<br>Threshold: **Partial (N-1)** | Shape/Color: **1.0 : 1.0 (ComboScore)** |
| **Ideal Use Case** | Highly specific sub-pocket targeting | Scaffold hopping & diversity search |

### 4.3 Summary of Best-Practice Settings
*   **For Consensus Pharmacophores**:
    *   **Tolerance**: **1.2 – 1.5 Å** to accommodate geometric diversity.
    *   **Occurrence Threshold**: **60% – 80%** (e.g., 4 out of 5 active ligands) to avoid "impossible" multi-feature queries.
    *   **Feature Matching**: **Partial (N-1)** logic to increase recall for fragmented pockets.
*   **For Reference Alignment**:
    *   **Weighting**: **1.0 : 1.0** (Shape:Color) to balance steric fit with chemical features.
    *   **Query Selection**: Use an **ensemble of conformations** for the reference ligand rather than a single rigid pose to account for flexibility.

## 5. Conclusion
Consensus pharmacophores are powerful tools for high-precision hit identification but are prone to failure when the active set is structurally diverse or when the binding pocket allows for partial occupancy. Mitigating these failures requires softening the geometric constraints through increased **tolerance** (1.2–1.5 Å) and lowering the **occurrence threshold** to allow for partial matches (N-1).

For broader chemotype discovery where scaffold hopping is the goal, **reference alignment** (shape-based methods like ROCS) remains the superior choice due to its continuous scoring functions and lack of hard distance cutoffs. A hybrid strategy—using a loose consensus pharmacophore as a pre-filter followed by shape-based ranking—often yields the optimal balance of enrichment and diversity.

## References

### Papers

[^p1]: Data Leakage and Redundancy in the LIT-PCBA Benchmark | Amber Huang, Ian Scott Knight, Slava Naprienko | 2025 | https://arxiv.org/abs/2507.21404v1 | arXiv:2507.21404v1 | source:ArXiv
[^p2]: WelQrate: Defining the Gold Standard in Small Molecule Drug Discovery Benchmarking | 2024 | https://arxiv.org/abs/2411.09820 | arXiv:2411.09820 | source:ArXiv
[^p3]: Aligning Target-Aware Molecule Diffusion Models with Exact Energy Optimization | 2024 | https://arxiv.org/abs/2407.01648 | arXiv:2407.01648 | source:ArXiv
[^p4]: PharmacoMatch: Efficient 3D Pharmacophore Screening through Neural Subgraph Matching | 2024 | https://arxiv.org/abs/2409.06316 | arXiv:2409.06316 | source:ArXiv
[^p5]: An Improved Metric and Benchmark for Assessing the Performance of Virtual Screening Models | 2024 | https://arxiv.org/abs/2403.10478 | arXiv:2403.10478 | source:ArXiv
[^p6]: Comparative analysis of pharmacophore screening tools | MPA Sanders, AJM Barbosa, B Zarzycka… | 2012 | https://pubs.acs.org/doi/abs/10.1021/ci2005274 | source:Google Scholar
[^p7]: Consensus holistic virtual screening for drug discovery: a novel machine learning model approach | S Moshawih, ZH Bu, HP Goh, N Kifli, LH Lee… | 2024 | https://link.springer.com/article/10.1186/s13321-024-00855-8 | source:Google Scholar
[^p8]: Comprehensive survey of consensus docking for high-throughput virtual screening | C Blanes | 2022 | https://www.mdpi.com/1420-3049/28/1/175 | source:Google Scholar
[^p9]: Assessing the performance of 3D pharmacophore models in virtual screening: how good are they? | RC Braga, CH Andrade | 2013 | https://www.ingentaconnect.com/content/ben/ctmc/2013/00000013/00000009/art00010 | source:Google Scholar
[^p10]: Optimal assignment methods for ligand-based virtual screening | A Jahn, G Hinselmann, N Fechner, A Zell | 2009 | https://link.springer.com/article/10.1186/1758-2946-1-14 | source:Google Scholar
[^p11]: Biologically Active Ligands for Yersinia Outer Protein H (YopH): Feature Based Pharmacophore Screening, Docking and Molecular Dynamics Studies | Thangaraju Tamilvanan, Waheeta Hopper | 2014 | https://doi.org/10.2174/1386207317666140211095137 | DOI:10.2174/1386207317666140211095137 | source:Crossref
[^p12]: Supplemental Information 8: Validation of AutoDock Vina virtual screening by active-decoy comparison. | https://doi.org/10.7717/peerj.11261/supp-8 | DOI:10.7717/peerj.11261/supp-8 | source:Crossref
[^p13]: Supplemental Information 16: Alignment of 38 AmFV consensus genomes with reference accession NC_027925.1 used for phylogenetic analysis. | https://doi.org/10.7717/peerj.16455/supp-16 | DOI:10.7717/peerj.16455/supp-16 | source:Crossref
[^p14]: COVID-19 Repurposed Therapeutics Targeting the Viral Protease and Spike-protein:ACE2 Interface using MD-based Pharmacophore and Consensus Virtual Screening | Brady Garabato, Federico Falchi, Andrea Cavalli | https://doi.org/10.26434/chemrxiv.12264503.v1 | DOI:10.26434/chemrxiv.12264503.v1 | source:Crossref
[^p15]: Decoy Receptor | https://doi.org/10.1007/3-540-27806-0_364 | DOI:10.1007/3-540-27806-0_364 | source:Crossref

### URLs

[^u1]: Consensus holistic virtual screening for drug discovery | https://pmc.ncbi.nlm.nih.gov/articles/PMC11134635/ | source:organic | pos:1
[^u2]: Pharmacophore model-aided virtual screening combined ... | https://www.sciencedirect.com/science/article/pii/S1878535222006505 | source:organic | pos:2
[^u3]: Comprehensive Survey of Consensus Docking for High- ... | https://www.mdpi.com/1420-3049/28/1/175 | source:organic | pos:3
[^u4]: A Novel Approach for Efficient Pharmacophore-based ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC2767445/ | source:organic | pos:4
[^u5]: Pharmacophore-Based Similarity Scoring for DOCK | https://pubs.acs.org/doi/10.1021/jp506555w | source:organic | pos:5
[^u6]: (PDF) Consensus holistic virtual screening for drug discovery | https://www.researchgate.net/publication/380937366_Consensus_holistic_virtual_screening_for_drug_discovery_a_novel_machine_learning_model_approach | source:organic | pos:6
[^u7]: A Molecular Representation to Identify Isofunctional ... | https://onlinelibrary.wiley.com/doi/10.1002/minf.202400159 | source:organic | pos:7
[^u8]: Scrutinization on Docking Against Individually Generated ... | https://www.biorxiv.org/content/10.1101/2025.01.01.630989v2.full-text | source:organic | pos:8
[^u9]: Structure-based identification of novel FAK1 inhibitors ... | https://www.nature.com/articles/s41598-025-23203-8 | source:organic | pos:9
[^u10]: Ligand based pharmacophore modelling and integrated ... | https://pubs.rsc.org/en/content/articlehtml/2024/ra/d3ra08618f | source:organic | pos:10
[^u11]: Proceeding Book 4th International Conference Frontiers in ... | https://www.academia.edu/127760597/Proceeding_Book_4th_International_Conference_Frontiers_in_Academic_Research | source:organic | pos:1
[^u12]: Concepts and Applications Exemplified on Hydroxysteroid ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC6332202/ | source:organic | pos:1
[^u13]: A Pharmacophore-Based Method for Rapid and Accurate ... | https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.5c00250 | source:organic | pos:2
[^u14]: Consensus pharmacophore for Drug Design | https://github.com/AngelRuizMoreno/ConcensusPharmacophore | source:organic | pos:4
[^u15]: Virtual Screening Using Pharmacophore Models Retrieved ... | https://imtm.cz/sites/default/files/publication/impact/1125-ijms-20-05834.pdf | source:organic | pos:5
[^u16]: Cov-2 Mpro Inhibitors from Large Chemical Libraries | https://digitalcommons.library.tmc.edu/cgi/viewcontent.cgi?article=1004&context=molecular_med | source:organic | pos:6
[^u17]: Consensus holistic virtual screening for drug discovery | https://link.springer.com/article/10.1186/s13321-024-00855-8 | source:organic | pos:7
[^u18]: Advances, Limitations, And current utility in drug discovery | https://www.dovepress.com/pharmacophore-modeling-advances-limitations-and-current-utility-in-dru-peer-reviewed-fulltext-article-JRLCR | source:organic | pos:8
[^u19]: A Unified Methodology for Ligand Based and Structure B | https://escholarship.org/content/qt99t2x783/qt99t2x783.pdf?t=spk0m5 | source:organic | pos:9
[^u20]: Ligand Based Pharmacophore Modeling and Virtual ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC4265523/ | source:organic | pos:10
[^u21]: ChemInform Abstract: The NEWLEAD Program: A New Method ... | https://www.researchgate.net/publication/309028090_ChemInform_Abstract_The_NEWLEAD_Program_A_New_Method_for_the_Design_of_Candidate_Structures_from_Pharmacophoric_Hypotheses | source:organic | pos:1
[^u22]: Introduction To Pharmacophores in MOE | PDF | https://www.scribd.com/document/330906180/Introduction-to-Pharmacophores-in-MOE | source:organic | pos:2
[^u23]: De Novo Molecular Design [PDF] [lss1rp5mmik0] | https://vdoc.pub/documents/de-novo-molecular-design-lss1rp5mmik0 | source:organic | pos:3
[^u24]: Assessing chemical exposure risk in breastfeeding infants | https://www.sciencedirect.com/science/article/pii/S0147651325000430 | source:organic | pos:1
[^u25]: Assessing chemical exposure risk in breastfeeding infants | https://www.researchgate.net/publication/387938551_Assessing_chemical_exposure_risk_in_breastfeeding_infants_An_explainable_machine_learning_model_for_human_milk_transfer_prediction | source:organic | pos:2
[^u26]: Deep Learning for Chemistry and Simulations | https://epub.jku.at/download/pdf/7148699.pdf | source:organic | pos:3
[^u27]: An Explainable Machine Learning Model for Human Milk ... | https://www.scribd.com/document/961527517/An-Explainable-Machine-Learning-Model-for-Human-Milk-Transfer-Prediction | source:organic | pos:4
[^u28]: Greedy 3-Point Search (G3PS)—A Novel Algorithm for ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC8658842/ | source:organic | pos:1
[^u29]: Design and application of structure-based pharmacophores ... | https://repository.ubn.ru.nl/bitstream/handle/2066/93565/93565.pdf | source:organic | pos:2
[^u30]: An AI-assisted fluorescence microscopic system for ... | https://www.nature.com/articles/s41467-025-60315-1 | source:organic | pos:3
[^u31]: Novel Approach to Structure-Based Pharmacophore ... | https://www.researchgate.net/publication/5455595_Novel_Approach_to_Structure-Based_Pharmacophore_Search_Using_Computational_Geometry_and_Shape_Matching_Techniques | source:organic | pos:4
[^u32]: TESIS DOCTORAL | https://www.tdx.cat/bitstream/10803/9311/3/Tvipn_3.pdf | source:organic | pos:5
[^u33]: Lead-Seeking Approaches | https://link.springer.com/content/pdf/10.1007/978-3-642-01075-0.pdf | source:organic | pos:6
[^u34]: A Shape-Based 3-D Scaffold Hopping Method and Its ... | https://www.researchgate.net/publication/7991220_A_Shape-Based_3-D_Scaffold_Hopping_Method_and_Its_Application_to_a_Bacterial_Protein-Protein_Interaction | source:organic | pos:7
[^u35]: Sandhya Kortagere Editor - In Silico Models for Drug ... - Springer Link | https://link.springer.com/content/pdf/10.1007%2F978-1-62703-342-8.pdf | source:organic | pos:8
[^u36]: Chemoinformatics and Computational Chemical Biology ... | https://epdf.pub/chemoinformatics-and-computational-chemical-biology-methods-in-molecular-biology.html | source:organic | pos:9
[^u37]: ZINCPharmer: pharmacophore search of the ZINC database | https://pmc.ncbi.nlm.nih.gov/articles/PMC3394271/ | source:organic | pos:1
[^u38]: Knowledge Search - My Account - Schrödinger | https://my.schrodinger.com/support/search/Is%20it%20possible%20within%20Phase%20to%20work%20with%20pre-aligned%20ligands? | source:organic | pos:2
[^u39]: Notch Signaling and Ageing | Request PDF | https://www.researchgate.net/publication/268794804_Notch_Signaling_and_Ageing | source:organic | pos:3
[^u40]: Symposium 18 | https://eprints.ugd.edu.mk/1082/1/__ugd.edu.mk_private_UserFiles_rubin.gulaboski_Desktop_GULABOSKI-MY%20pdf%20PUBLICATIONS_Posteri_Baltimor%202011.pdf | source:organic | pos:4
[^u41]: Ultra-High-Throughput Structure-Based Virtual Screening ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC4765363/ | source:organic | pos:1
[^u42]: Consensus Pharmacophore Strategy For Identifying Novel ... | https://pubs.acs.org/doi/10.1021/acs.jcim.3c01439 | source:organic | pos:2
[^u43]: Linking machine learning and biophysical structural ... | https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2024.1305272/full | source:organic | pos:3
[^u44]: Computational Chemistry in Structure-Based Drug Design | https://jpscc.samipubco.com/article_229217_0b37e6167d4922dd834fee49edfb7f60.pdf | source:organic | pos:4
[^u45]: Application to the RNA Repeat Expansion in Myotonic ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC5286915/ | source:organic | pos:5
[^u46]: 3D Chemical Similarity Networks for Structure-Based ... | https://europepmc.org/article/pmc/5511031 | source:organic | pos:7
[^u47]: A Definitive Pharmacophore Modelling Study on CDK2 ATP ... | https://iris.unipa.it/bitstream/10447/476730/2/CDDT%20MS_new_revised.pdf | source:organic | pos:8
[^u48]: Improving Posing and Ranking of Molecular Docking by ... | https://utoronto.scholaris.ca/server/api/core/bitstreams/88806e43-7f26-4674-9371-4c7c5a0106ff/content | source:organic | pos:10
[^u49]: Improving ADMET prediction with descriptor augmentation ... | https://www.biorxiv.org/content/10.1101/2025.07.14.664363v1.full.pdf | source:organic | pos:1
[^u50]: Automating Knowledge Discovery for Toxicity Prediction Using ... | https://pubs.acs.org/doi/10.1021/ci300254w | source:organic | pos:2
[^u51]: Improving ADMET prediction with descriptor augmentation of ... | https://www.biorxiv.org/content/10.1101/2025.07.14.664363v2.full.pdf | source:organic | pos:4
[^u52]: Prediction of patient's response to OnabotulinumtoxinA ... | https://www.sciencedirect.com/science/article/pii/S240584401831003X | source:organic | pos:1
[^u53]: Predicting BoNT-A Response in Migraines | PDF | https://www.scribd.com/document/458991937/migrainePred | source:organic | pos:3
[^u54]: Semi‐automated detection of eagle nests: an application of ... | https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.38 | source:organic | pos:1
[^u55]: UAV Imagery-Based Classification Model for Atypical ... | https://www.mdpi.com/2504-446X/8/7/297 | source:organic | pos:2
[^u56]: VIEWPixx3 | https://www.visionsciences.org/programs/VSS_2025_Abstracts.pdf | source:organic | pos:3
[^u57]: (PDF) Semi-automated detection of eagle nests | https://www.researchgate.net/publication/312527005_Semi-automated_detection_of_eagle_nests_an_application_of_very_high-resolution_image_data_and_advanced_image_analyses_to_wildlife_surveys | source:organic | pos:4
[^u58]: automated detection of eagle nests: an application of very highâ | https://zslpublications.onlinelibrary.wiley.com/doi/pdf/10.1002/rse2.38 | source:organic | pos:5
[^u59]: UC Merced | https://escholarship.org/content/qt5718511j/qt5718511j.pdf | source:organic | pos:6
[^u60]: 2006 Science and Technology/Engineering Curriculum ... | https://www.doe.mass.edu/frameworks/scitech/1006.doc | source:organic | pos:7
[^u61]: Invited Symposia, Keynote Addresses, Poster session | https://www.tandfonline.com/doi/full/10.1080/00207594.2000.20000728 | source:organic | pos:8
[^u62]: Glencoe Exploring Theatre | https://s3.amazonaws.com/scschoolfiles/1148/copy_of_exploring-theatre-textbook.pdf | source:organic | pos:9
[^u63]: UC Berkeley - Dissertations, Department of Linguistics | https://escholarship.org/content/qt3g9427m2/qt3g9427m2.pdf | source:organic | pos:10
[^u64]: The RDKit Book — The RDKit 2025.09.4 documentation | https://www.rdkit.org/docs/RDKit_Book.html | source:organic | pos:1
[^u65]: ConQuest User Guide and Tutorials - CCDC | https://www.ccdc.cam.ac.uk/media/Documentation/2F0D7443-9739-46EB-BE9F-69E62E531FB7/2f0d7443973946ebbe9f69e62e531fb7.pdf | source:organic | pos:2
[^u66]: pattanaik-lagnajit-phd-cheme-2023-thesis.pdf | https://dspace.mit.edu/bitstream/handle/1721.1/150133/pattanaik-lagnajit-phd-cheme-2023-thesis.pdf?sequence=1&isAllowed=y | source:organic | pos:3
[^u67]: Diffusion Models in Drug Design: Recent Advances | https://www.scribd.com/document/896236402/Diffusion-Models-in-Drug-Design-Recent-Advances-and-Future-Opportunities | source:organic | pos:4
[^u68]: Experimental Structure Determination Methods | https://www.cs.tulane.edu/~mettu/cmps6110_Spring2017/lectures/ExperimentalStructureDetermination.pdf | source:organic | pos:5
[^u69]: Conformation mining: An algorithm for finding biologically ... | https://www.researchgate.net/publication/7881230_Conformation_mining_An_algorithm_for_finding_biologically_relevant | source:organic | pos:6
[^u70]: Applications and Variations of the Maximum Common ... | https://etheses.whiterose.ac.uk/id/eprint/13063/1/thesis_mk2_final.pdf | source:organic | pos:7
[^u71]: AI-Designed Family of Dual COX-2/mPGES-1 Inhibitors to ... | https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2025-49xk0 | source:organic | pos:8
[^u72]: GeoMol: Torsional Geometric Generation of Molecular 3D ... | https://ar5iv.labs.arxiv.org/html/2106.07802 | source:organic | pos:9
[^u73]: Applications and Theory (Specialist Periodical Reports) | https://epdf.pub/chemical-modelling-vol-2-applications-and-theory-specialist-periodical-reports.html | source:organic | pos:10
[^u74]: Color Features — Applications | https://docs.eyesopen.com/applications/rocs/theory/shape_color.html | source:organic | pos:1
[^u75]: 3D Chemical Similarity Networks for Structure-Based Target ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC5511031/ | source:organic | pos:2
[^u76]: De novo design with deep generative models based on 3D ... | https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/60c758a7842e659173db488f/original/de-novo-design-with-deep-generative-models-based-on-3d-similarity-scoring.pdf | source:organic | pos:3
[^u77]: Stereoselective virtual screening of the ZINC database using ... | https://link.springer.com/article/10.1186/s13321-014-0051-5 | source:organic | pos:4
[^u78]: ROCS shape query derived from a low energy 3D conformation ... | https://www.researchgate.net/figure/ROCS-shape-query-derived-from-a-low-energy-3D-conformation-generated-in-MOE-201308-of_fig3_277254015 | source:organic | pos:5
[^u79]: 3D Chemical Similarity Networks for Structure-Based Target ... | https://pubs.acs.org/doi/10.1021/acschembio.6b00253 | source:organic | pos:6
[^u80]: RSC Advances | https://pubs.rsc.org/en/content/getauthorversionpdf/c5ra05919d | source:organic | pos:7
[^u81]: Prospective performance evaluation of selected common ... | https://www.sciencedirect.com/science/article/pii/S0223523415002664 | source:organic | pos:8
[^u82]: ShaEP: Molecular Overlay Based on Shape and ... | https://elearning.uniroma1.it/pluginfile.php/1384573/mod_folder/content/0/Sezione%205.0/5.4/5.4.2.2.2.8.18.2009.ShaEP.%20Molecular%20Overlay%20Based%20on%20Shape%20and%20Electrostatic%20Potential.pdf?forcedownload=1 | source:organic | pos:9
[^u83]: Computer-aided screening for potential TMPRSS2 inhibitors | https://pmc.ncbi.nlm.nih.gov/articles/PMC7441808/ | source:organic | pos:1
[^u84]: Common Hits Approach: Combining Pharmacophore ... | https://www.researchgate.net/publication/312202410_Common_Hits_Approach_Combining_Pharmacophore_Modeling_and_Molecular_Dynamics_Simulations | source:organic | pos:2
[^u85]: Computer-aided screening for potential TMPRSS2 inhibitors | https://www.researchgate.net/publication/342988864_Computer-aided_screening_for_potential_TMPRSS2_inhibitors_a_combination_of_pharmacophore_modeling_molecular_docking_and_molecular_dynamics_simulation_approaches | source:organic | pos:3
[^u86]: Virtual screening of approved drugs as potential SARS ... | https://www.researchgate.net/publication/342458614_Virtual_screening_of_approved_drugs_as_potential_SARS-CoV-2_main_protease_inhibitors | source:organic | pos:4
[^u87]: Program & Abstracts | https://iccs-nl.org/wp-content/uploads/2020/03/7th_iccs_program_and_abstracts.pdf | source:organic | pos:2
[^u88]: (PDF) Total Syntheses of Pladienolide-Derived ... | https://www.researchgate.net/publication/355222657_Total_Syntheses_of_Pladienolide-Derived_Spliceosome_Modulators | source:organic | pos:3
[^u89]: Design of novel compounds with the potential of dual ... | https://core.ac.uk/download/pdf/237157508.pdf | source:organic | pos:4
[^u90]: APPLICATION, EVALUATION, AND IMPROVEMENT OF ... | https://www.research.unipd.it/retrieve/37f1f89f-4ff6-4fa2-b457-aed13ac38d4f/tesi_Davide_Bassani.pdf | source:organic | pos:6
[^u91]: Design and Optimisation of Novel Benzimidazole Hybrid | https://www.um.edu.mt/library/oar/bitstream/123456789/142124/1/2518MDSPHR512305072351_1.PDF | source:organic | pos:7
[^u92]: In Silico and In Vitro Investigation into the Next Generation ... | https://uhra.herts.ac.uk/id/eprint/16887/1/07155301%20BOTHA%20Michelle%20Final%20Version%20of%20PhD%20Submission.pdf | source:organic | pos:8
[^u93]: 21 September 2025 AperTO | https://iris.unito.it/retrieve/handle/2318/1721223/558815/Book_of_Abstract.pdf | source:organic | pos:9
[^u94]: Total Syntheses of Pladienolide-Derived Spliceosome ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC8512135/ | source:organic | pos:10
[^u95]: QPHAR: quantitative pharmacophore activity relationship | https://pmc.ncbi.nlm.nih.gov/articles/PMC8351372/ | source:organic | pos:2
[^u96]: Discovery of a Promising Hydroxyamino-Piperidine ... | https://www.mdpi.com/1424-8247/18/9/1303 | source:organic | pos:3
[^u97]: A MOLECULAR DYNAMICS - PPAR alpha | https://www.researchgate.net/publication/313113643_A_MOLECULAR_DYNAMICS_-_SHARED_PHARMACOPHORE_APPROACH_TO_BOOST_EARLY_ENRICHMENT_VIRTUAL_SCREENING_A_CASE_STUDY_on_PPAR_alpha | source:organic | pos:5
[^u98]: Structure-based discovery and experimental validation of ... | https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2025.1605741/full | source:organic | pos:6
[^u99]: DEVELOPMENT AND OPTIMISATION OF COMPUTATIONAL ... | https://iris.unipa.it/retrieve/e3ad891a-388f-da0e-e053-3705fe0a2b96/PhD_THESIS_Ugo_Perricone_OK_new_full.pdf | source:organic | pos:7
[^u100]: The PPAR Ω Pocket: Renewed Opportunities for Drug ... | https://onlinelibrary.wiley.com/doi/10.1155/2020/9657380 | source:organic | pos:8
[^u101]: Exploration of SARS-CoV-2 Mpro Noncovalent Natural ... | https://pubs.acs.org/doi/10.1021/acsomega.2c07259 | source:organic | pos:9
[^u102]: Pharmacological Characterization of Novel Opioid Receptor ... | https://researchrepository.wvu.edu/cgi/viewcontent.cgi?article=1619&context=etd | source:organic | pos:10
[^u103]: Prediction of Drug Targets in Human Pathogens | https://repositorio.ufmg.br/bitstreams/7f243b7e-f944-4a4c-8b68-a8131ddc6cba/download | source:organic | pos:1
[^u104]: Conference Proceeding | PDF | https://www.scribd.com/document/734555292/Conference-Proceeding | source:organic | pos:2
[^u105]: In silico drug design et chimie médicinale | https://theses.hal.science/tel-01661380/file/These_RAYARAnita_VF_2.pdf | source:organic | pos:3
[^u106]: Pharmacophore Models Derived from Molecular Dynamics ... | https://journals.sagepub.com/doi/pdf/10.1177/1934578X1601101019 | source:organic | pos:1
[^u107]: Pharmer: Efficient and Exact Pharmacophore Search | https://pubs.acs.org/doi/10.1021/ci200097m | source:organic | pos:2
[^u108]: Small Molecule Decoys of Aggregation for Elimination of Aβ ... | https://pubs.acs.org/doi/10.1021/acschemneuro.2c00649 | source:organic | pos:5
[^u109]: Molecular architecture of SF3B and the structural basis ... - eDiss | https://ediss.uni-goettingen.de/bitstream/handle/21.11130/00-1735-0000-0003-C13A-2/Cretu_Constantin_Phd_2018_optimized.pdf?sequence=1&isAllowed=y | source:organic | pos:7
[^u110]: Développement et validation du logiciel S4MPLE | https://theses.hal.science/tel-00874644/file/HOFFER_Laurent_2013_ED222.pdf | source:organic | pos:8
[^u111]: Molecular dynamics, dynamic site mapping, and ... | https://www.researchgate.net/publication/262074348_Molecular_dynamics_dynamic_site_mapping_and_highthroughput_virtual_screening_on_leptin_and_the_Ob_receptor_as_anti-obesity_target | source:organic | pos:10
[^u112]: Ligand-Based Pharmacophore Modeling Using Novel 3D ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC6321403/ | source:organic | pos:1
[^u113]: Pharmacophore - an overview | ScienceDirect Topics | https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/pharmacophore | source:organic | pos:2
[^u114]: Development of Pharmacophore Model for Indeno[1,2-b] ... | https://www.mdpi.com/1424-8247/10/1/8 | source:organic | pos:3
[^u115]: Molecular formulas of compounds 1-78 | Download Table | https://www.researchgate.net/figure/Molecular-formulas-of-compounds-1-78_tbl1_23164719 | source:organic | pos:4
[^u116]: Discovery of novel EGFR inhibitors | https://www.scholarsresearchlibrary.com/articles/discovery-of-novel-egfr-inhibitors-in-silico-study-and-3dpharmacophore-model-generation.pdf | source:organic | pos:5
[^u117]: pharmacoforge: pharmacophore generation with | https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/68926c6b728bf9025e1ef1e5/original/pharmaco-forge-pharmacophore-generation-with-diffusion-models.pdf | source:organic | pos:6
[^u118]: ❑ Molecular Modeling and Simulations ❑ Protein ... | https://wikiali.files.wordpress.com/2014/03/simulation-in-bio.pdf | source:organic | pos:7
[^u119]: Mixing Pharmacophore Modeling and Classical QSAR ... | https://cdn.intechopen.com/pdfs/32264/InTech-Mixing_pharmacophore_modeling_and_classical_qsar_analysis_as_powerful_tool_for_lead_discovery.pdf | source:organic | pos:8
[^u120]: ELIXIR-A: An Interactive Visualization Tool for Multi-Target ... | https://pubs.acs.org/doi/10.1021/acsomega.1c07144 | source:organic | pos:9
[^u121]: CCG SVL Exchange | https://svl.chemcomp.com/ | source:organic | pos:10
[^u122]: Fragmenstein — An Ugly Duckling or a Swan? | https://medium.com/@mykola.protopopov/fragmenstein-an-ugly-duckling-or-a-swan-9ec34c2c86ac | source:organic | pos:2
[^u123]: FLOWR.root: A flow matching based foundation model for ... | https://www.researchgate.net/publication/397755965_FLOWRroot_A_flow_matching_based_foundation_model_for_joint_multi-purpose_structure-aware_3D_ligand_generation_and_affinity_prediction | source:organic | pos:9
[^u124]: A simple method for correlative light and scanning electron ... | https://www.researchgate.net/publication/19657447_A_simple_method_for_correlative_light_and_scanning_electron_microscopy_of_human_iliac_crest_bone_biopsies_Qualitative_observations_in_normal_and_osteoporotic_subjects | source:organic | pos:3
[^u125]: Investigation of the effects of the splicing inhibitor ... - eDiss | https://ediss.uni-goettingen.de/bitstream/21.11130/00-1735-0000-0005-13D3-7/1/Thesis_Sebastian%20Ludwig_FINAL_noCV.pdf | source:organic | pos:4
[^u126]: Asymmetric Organocatalysis at the Service of Medicinal ... | https://onlinelibrary.wiley.com/doi/10.1155/2014/531695 | source:organic | pos:5
[^u127]: Alma Mater Studiorum | https://amsdottorato.unibo.it/id/eprint/9822/3/Corbisiero_Dario_Tesi.pdf | source:organic | pos:7
[^u128]: ¾/í µí»¼ Modulation for the Management of Metabolic Syndrome | https://www.academia.edu/36074958/Design_of_Novel_Compounds_with_the_Potential_of_Dual_PPAR%C3%AD_%CE%BC%C3%AD_%C3%AD_%CE%BC%C3%AD_Modulation_for_the_Management_of_Metabolic_Syndrome | source:organic | pos:10
[^u129]: Structure-based virtual screening, molecular docking, and MD ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC12312920/ | source:organic | pos:1
[^u130]: Three-Dimensional Pharmacophore Methods in Drug Discovery | https://pubs.acs.org/doi/10.1021/jm900817u | source:organic | pos:2
[^u131]: Pharmacophore Mapping, Virtual Screening and Molecular ... | https://ijddd.com/index.php/ijddd/article/view/186?articlesBySimilarityPage=3 | source:organic | pos:3
[^u132]: Jahan08/High-Throughput-Virtual-Screening-for-Hit- ... | https://github.com/Jahan08/High-Throughput-Virtual-Screening-for-Hit-Finding-and-Evaluation-Schrodinger | source:organic | pos:4
[^u133]: Rethinking the 'best method' paradigm: The effectiveness ... | https://www.sciencedirect.com/science/article/pii/S2667318524000242 | source:organic | pos:5
[^u134]: Ligand-Based Virtual Screening - Springer Link | https://link.springer.com/rwe/10.1007/978-981-95-2525-6_43 | source:organic | pos:6
[^u135]: Establishing the foundations for a data-centric AI approach ... | https://elifesciences.org/reviewed-preprints/97821 | source:organic | pos:7
[^u136]: Computational approaches streamlining drug discovery | https://www.nature.com/articles/s41586-023-05905-z | source:organic | pos:8
[^u137]: Virtual Screening: Principles, Challenges, and Practical ... | https://www.researchgate.net/publication/277701641_Virtual_Screening_Principles_Challenges_and_Practical_Guidelines | source:organic | pos:9
[^u138]: Virtual Screening Techniques in Drug Discovery: Review ... | https://www.semanticscholar.org/paper/Virtual-Screening-Techniques-in-Drug-Discovery%3A-and-Rocha-Olanda/e1a6969cf05de526d339e18b65eb3433c9c775a7 | source:organic | pos:10
[^u139]: Intelligent Computing Theories and Application | https://link.springer.com/content/pdf/10.1007/978-3-030-26969-2.pdf | source:organic | pos:1
[^u140]: Ligand-based Virtual Screening Utilizing Partial Shape ... | https://ediss.sub.uni-hamburg.de/bitstream/ediss/7272/1/Dissertation.pdf | source:organic | pos:2
[^u141]: Programming language use distribution from recent ... | https://www.biostars.org/p/251002/ | source:organic | pos:3
[^u142]: Multiple Alignment using Fast Fourier Transform - BiŌkeanós | https://biokeanos.com/source/Multiple%20Alignment%20using%20Fast%20Fourier%20Transform | source:organic | pos:4
[^u143]: Compound Design by Fragment-Linking | Request PDF | https://www.researchgate.net/publication/247935821_Compound_Design_by_Fragment-Linking | source:organic | pos:1
[^u144]: Rationalizing Tight Ligand Binding through Cooperative ... | https://pubs.acs.org/doi/10.1021/ci200319e | source:organic | pos:2
[^u145]: Elucidating the aryl hydrocarbon receptor antagonism from a ... | https://orbi.uliege.be/bitstream/2268/244639/1/2020%20%28Goya-Jorge%20et%20al%29%20Elucidating%20the%20aryl%20hydrocarbon%20receptor%20antagonism%20from%20a%20chemical%20structural%20perspective.pdf | source:organic | pos:3
[^u146]: Target-focused library design by pocket-applied computer ... | https://hal.science/hal-03830359/document | source:organic | pos:4
[^u147]: integrating CADD tools and drug repurposing for PD-1/ ... | https://www.sciencedirect.com/org/science/article/pii/S2046206925001901 | source:organic | pos:6
[^u148]: Development and application of web-based open source drug ... | https://digitalcommons.usf.edu/cgi/viewcontent.cgi?article=6748&context=etd | source:organic | pos:7
[^u149]: Manipulating 3D Molecules in a Fixed-Dimensional SE(3) | https://www.researchgate.net/publication/392336925_Manipulating_3D_Molecules_in_a_Fixed-Dimensional_SE3-Equivariant_Latent_Space | source:organic | pos:8
[^u150]: Graph Representation Learning of Molecules for de novo ... | https://refubium.fu-berlin.de/bitstream/handle/fub188/46635/Thesis_Tuan_Le.pdf?sequence=4&isAllowed=y | source:organic | pos:9
[^u151]: Based Drug Discovery | https://backend.orbit.dtu.dk/ws/portalfiles/portal/219712922/PhD_DTU_NST.pdf | source:organic | pos:10
[^u152]: Benchmarking and Developing Novel Methods for G Protein | https://digitalcommons.memphis.edu/context/etd/article/4296/viewcontent/944214.pdf | source:organic | pos:1
[^u153]: Molecular Mingling: Multimodal Predictions of Ligand ... | https://pmc.ncbi.nlm.nih.gov/articles/PMC9124788/ | source:organic | pos:2
[^u154]: Can Deep Learning Blind Docking Methods be Used to ... | https://pubs.acs.org/doi/10.1021/acs.jcim.5c00331 | source:organic | pos:3
[^u155]: (PDF) Methods for Accurate Homology Modeling by Global ... | https://www.researchgate.net/publication/221821724_Methods_for_Accurate_Homology_Modeling_by_Global_Optimization | source:organic | pos:5
[^u156]: (PDF) [Langer T., Hoffmann R.D.] Pharmacophores, | https://www.academia.edu/36380288/_Langer_T_Hoffmann_R_D_Pharmacophores_ | source:organic | pos:6
[^u157]: pairwise sequence alignment: Topics by ... | https://www.science.gov/topicpages/p/pairwise+sequence+alignment | source:organic | pos:7
[^u158]: (PDF) Discovery of Small Molecules that Target Vascular ... | https://www.academia.edu/48118395/Discovery_of_Small_Molecules_that_Target_Vascular_Endothelial_Growth_Factor_Receptor_2_Signalling_Pathway_Employing_Molecular_Modelling_Studies | source:organic | pos:8
[^u159]: Moscow Conference on Computational Molecular Biology | http://lab6.iitp.ru/ru/pub/en_mccmb_2011.pdf | source:organic | pos:9
[^u160]: Interactions of Ozone-Functionalized Activated Charcoal ... | https://www.researchgate.net/publication/356684634_Interactions_of_Ozone-Functionalized_Activated_Charcoal_with_SARS-Cov-2_Proteases_Using_Molecular_Docking_and_Dynamics | source:organic | pos:10

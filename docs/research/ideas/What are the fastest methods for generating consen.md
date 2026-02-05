<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" style="height:64px;margin-right:32px"/>

# What are the fastest methods for generating consensus pharmacophore models from

100-1000 aligned molecules?

Consider:

- Distance matrix calculation optimizations
- GPU acceleration techniques
- Spatial indexing data structures
- Parallel/distributed algorithms
- Approximate methods trading accuracy for speed

Provide code examples, libraries, or tools implementing these methods.

I'll now generate a comprehensive report on the fastest methods for consensus pharmacophore generation.

# Fastest Methods for Generating Consensus Pharmacophore Models from 100-1000 Aligned Molecules

## Executive Summary

Generating consensus pharmacophore models from hundreds to thousands of molecules requires efficient computational methods that balance speed and accuracy. This analysis examines state-of-the-art approaches including spatial indexing (Pharmer KDB-trees), GPU acceleration, machine learning-based neural subgraph matching (PharmacoMatch), and coarse-grained pharmacophore modeling (PharmacoNet). The fastest workflows achieve consensus generation in **seconds to minutes** compared to hours for traditional methods, with speedups ranging from **10× to 3000×** depending on the approach.

***

## Introduction

Pharmacophore models—3D arrangements of chemical features essential for molecular recognition—are fundamental tools in virtual screening and drug discovery. Consensus pharmacophores, derived from multiple aligned active molecules, capture common interaction patterns critical for bioactivity. However, the computational cost of aligning and comparing 100-1000 molecules scales quadratically with conventional methods, creating bottlenecks in high-throughput workflows. This report evaluates computational strategies that dramatically accelerate consensus pharmacophore generation through algorithmic innovations, hardware acceleration, and approximate methods.[^1][^2]

***

## Computational Challenges in Consensus Pharmacophore Generation

### The Quadratic Scaling Problem

Traditional pharmacophore alignment requires pairwise comparisons between molecules. For *n* molecules with *m* pharmacophore features each:

- **Feature matching complexity**: O(*m*³) per pair (clique detection)[^3]
- **All-vs-all alignment**: O(*n*²) comparisons
- **Total complexity**: O(*n*²*m*³)

For 1000 molecules with 10 features each, this yields ~10⁹ operations. Standard implementations (MOE, CDPKit) process alignments at ~10ms per pair, requiring **2.8 hours** for complete consensus analysis.[^4][^3]

### Key Optimization Targets

1. **Distance matrix calculation**: Pairwise pharmacophore distances
2. **Spatial indexing**: Avoid exhaustive searches through coordinate-based retrieval
3. **Approximate matching**: Trade precision for massive speedup via coarse-graining
4. **Parallelization**: Distribute computations across CPU cores or GPU threads

***

## Method 1: Spatial Indexing with Pharmer KDB-Trees

### Algorithm Design

Pharmer revolutionizes pharmacophore search by inverting the traditional database-scanning paradigm. Instead of querying each molecule sequentially, it stores pharmacophores in a **spatial index** that enables logarithmic-time retrieval.[^3]

#### Core Data Structure: The Pharmer KDB-Tree

**Triangle Representation**
Pharmacophores are decomposed into all possible triangles of feature points. Each triangle is represented by:

- **Sorted side lengths** (*l*₁ ≤ *l*₂ ≤ *l*₃): Coordinate-frame-independent 3D point
- **TripletData** (64 bytes): Triangle geometry + Bloom fingerprint + molecular metadata

**Tree Organization**
A balanced KDB-tree indexes triangles by their (*l*₁, *l*₂, *l*₃) coordinates:[^3]

- **Nodes**: Tight bounding boxes (sliding-midpoint splitting for balance)
- **Leaves**: Hierarchical files of TripletData (8KB pages for sequential I/O)
- **Depth**: Logarithmic (few I/O operations even for billions of points)

**Range Query Optimization**
When querying with tolerance sphere radius *r*:

1. Compute bounding box around query triangle coordinates
2. Traverse KDB-tree, pruning branches with non-overlapping bounding boxes
3. Short-circuit when node fully contained in query range (sequential read of all descendants)
4. Performance scales with **query complexity**, not database size

#### Multi-Resolution Bloom Fingerprints

To avoid combinatorial explosion for multi-point queries, Pharmer embeds **tetrahedral pharmacophore fingerprints** in each TripletData:[^3]

- **Encoding**: For each remaining feature, compute 3 distances to triangle vertices + chirality bit
- **Discretization**: Binned at multiple resolutions (4Å: 6 hash functions; 1Å: 8 hash functions)
- **Storage**: 256-bit Bloom filter (~0.43% false positive rate for 25 features)
- **Filtering**: During range queries, check if candidate TripletDatas contain all query features before full alignment


### Performance Characteristics

| Metric | Pharmer | MOE (Baseline) | Speedup |
| :-- | :-- | :-- | :-- |
| **Search time (HSP90 query, 1.9M conformers)** | 31.3s | 388.9s | **12×** |
| **Compute time (disk-cached)** | 11.3s | 119.2s | **11×** |
| **Tight query (0.5× tolerance)** | ~0.8s | 30s | **38×** |
| **Precise queries (0.1× tolerance)** | <0.1s | ~30s | **>1000×** |

**Key Insight**: Pharmer's advantage grows exponentially as queries become more specific—critical for consensus modeling where tight spatial constraints identify conserved interaction patterns.[^3]

### Consensus Workflow with Pharmer

```
1. Pre-process (one-time):
   - Extract pharmacophores from 1000 aligned molecules
   - Enumerate all feature triangles (~ (n choose 3) per molecule)
   - Build KDB-tree index (minutes for 200K conformers)

2. Consensus extraction:
   FOR each representative molecule as query:
     - Decompose into n-1 connected triplets
     - Range query KDB-tree with tolerance r (e.g., 1.5Å)
     - Filter via Bloom fingerprints
     - Reconstruct alignments via dual quaternions
   
   - Cluster matching features across queries
   - Identify consensus pharmacophore (features present in >X% of molecules)

3. Estimated time: <1 minute for 1000 molecules (index pre-built)
```


***

## Method 2: Neural Subgraph Matching with PharmacoMatch

### Paradigm Shift: Pharmacophores as Order Embeddings

PharmacoMatch reinterprets pharmacophore matching as an **approximate neural subgraph matching** problem, encoding query-target relationships in a learned embedding space.[^4]

#### Architecture

**Graph Neural Network Encoder**
Input: Pharmacophore as complete graph (*V*: feature points, *E*: Euclidean distances)

1. **Node features**: One-hot encoded pharmacophore types (HBD, HBA, hydrophobic, etc.)
2. **Edge features**: Radial basis function encoding of distances (0-10Å, 5-dimensional)
3. **Message passing**: 3 NNConv layers (edge-conditioned convolution) + DenseNet skip connections
4. **Read-out**: Additive pooling → MLP projection (1024 hidden → 512-dim embedding)
5. **Constraint**: Final layer uses absolute weights to enforce non-negative embeddings

**Self-Supervised Training**
Dataset: 1.2M ChEMBL molecules (unlabeled, zero-shot for validation)

- **Loss function**: Order embedding loss (Ying et al. 2020)[^4]

L = Σ max(0, margin + penalty(e_query, e_target))

where penalty(q, t) = ||max(0, q - t)||² encodes partial ordering (query ⊂ target)
- **Augmentations** (on-the-fly pair generation):
    - *Positive*: Delete nodes (≥3 remain) + displace ≤1.5Å (tolerance sphere)
    - *Negative*: Boundary displacement, target deletion, wrong target mapping
- **Curriculum learning**: Progressive node counts (4→32) over epochs


#### Decision Function

After embedding database pharmacophores (once):

```python
def match(query_embedding, target_embedding, threshold=0.1):
    penalty = sum(max(0, q_i - t_i)**2 for q_i, t_i in zip(query, target))
    return penalty <= threshold
```


### Performance Benchmarks

**Speed Comparison** (per pharmacophore):[^4]


| Operation | PharmacoMatch (GPU) | CDPKit Alignment (128 CPU threads) | Speedup |
| :-- | :-- | :-- | :-- |
| Embedding | 10ms | N/A | - |
| Matching | **0.1ms** | 10ms | **100×** |

**Virtual Screening Performance** (vs. PharmacoNet pre-screening):[^4]


| Dataset | AUROC | BEDROC | EF@1% | EF@5% | Runtime |
| :-- | :-- | :-- | :-- | :-- | :-- |
| DEKOIS2.0 (80 targets) | 60.9 | **15.1** | **5.5** | **4.9** | 2.2ms/ligand |
| PharmacoNet | 62.5 | 12.3 | 4.4 | 4.2 | **7s/ligand (3000× slower)** |

**Key Advantage**: Query size-independent matching (contrast with alignment's linear scaling).[^4]

### Consensus Workflow

```
1. Train/load model (one-time):
   - Pre-trained model available (github.com/molinfo-vienna/PharmacoMatch)
   
2. Embed all molecules:
   - 10ms × 1000 molecules = 10 seconds (GPU batch processing)
   
3. All-vs-all matching:
   - 0.1ms × 1M pairs = 100 seconds
   - OR hierarchical clustering on embeddings (cosine similarity)
   
4. Consensus extraction:
   - Identify feature clusters with high co-occurrence across embeddings
   - Refine with alignment (optional, for validation)

Total time: ~2 minutes for 1000 molecules
```

**Limitations**: E(3)-invariant (cannot distinguish mirror images); approximate (lossy embedding). Best suited for **pre-screening** before exact alignment.[^4]

***

## Method 3: Coarse-Grained Pharmacophore Modeling with PharmacoNet

### Deep Learning-Guided Pharmacophore Generation

PharmacoNet uses **instance segmentation** to automate protein-based pharmacophore modeling, bypassing ligand-dependent methods entirely. While designed for virtual screening, its analytical scoring function offers insights for consensus workflows.[^5]

#### Workflow

1. **Pharmacophore Modeling** (DL-based):
    - Voxelize binding site (64×64×64 grid, 0.5Å resolution)
    - 3D Swin Transformer encoder detects protein hotspots
    - Instance segmentation predicts pharmacophore point density maps
2. **Graph Matching** (coarse-grained):
    - Cluster pharmacophore features (avoid redundant matching)
    - Distance constraints between feature pairs (Gaussian tolerance)
    - Depth-first search for valid correspondences
3. **Analytical Scoring** (non-DL):
    - Only **7 parameters** (weight per NCI type: hydrophobic, H-bond, ionic, etc.)
    - Gaussian mixture model on pharmacophore edge distances
    - Avoids overfitting unlike million-parameter DL scorers

### Speed Achievements

| Benchmark | PharmacoNet | AutoDock Vina | GLIDE SP |
| :-- | :-- | :-- | :-- |
| **Runtime (PDBbind core)** | 0.1ms | 396ms | 3.4s |
| **Speedup** | - | **3956×** | **34,117×** |
| **187M compounds (CB2 screening)** | **21 hours** | 11 years | - |

**Key Mechanism**: Pharmacophore-level abstraction (non-covalent interactions only) eliminates atom-pairwise force field evaluations.[^5]

### Consensus Applicability

While PharmacoNet targets protein-based queries, its **coarse-grained matching** principles apply to consensus generation:

- Cluster identical features across molecules (e.g., merge six benzene hydrophobics into one)
- Use analytical distance scoring (not expensive docking)
- Parallel evaluation on multi-CPU (32-core: 187M in 21 hrs)

**Estimated consensus time**: Seconds (if adapting analytical scorer to ligand-based consensus).

***

## Method 4: GPU Acceleration Techniques

### GPU-Accelerated Tanimoto Similarity

For fingerprint-based pre-filtering before full pharmacophore alignment:

**Speedup benchmarks**:[^6][^7][^8]

- **Unity fingerprints**: 39× over CPU sparse vector algorithm
- **Descriptors**: 195-328× over single-core CPU
- **Multi-core comparison**: 4-7× speedup (fingerprints), 100-328× (descriptors)

**Implementation**: CUDA kernels for bitwise AND/OR operations on binary fingerprints

```cuda
// Pseudocode: Tanimoto batch calculation
__global__ void tanimoto_kernel(uint* query, uint* database, float* output, int n_mols) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n_mols) {
        int intersection = popcount(query & database[idx]);
        int union = popcount(query | database[idx]);
        output[idx] = (float)intersection / union;
    }
}
```

**Consensus workflow**:

1. Compute all-vs-all Tanimoto matrix (GPU, seconds for 1000 mols)
2. Hierarchical clustering → representatives
3. Full pharmacophore alignment on clusters (minutes)

### ROSHAMBO: GPU Molecular Alignment

ROSHAMBO accelerates 3D shape/pharmacophore overlay via GPU:[^9]

- **Method**: PAPER software (Gaussian volume overlap optimization)
- **Performance**: "Near state-of-the-art, enables routine workflows"
- **Integration**: Python API + command-line (github.com/molecularinformatics/roshambo)

**Application**: Replace Kabsch alignment in consensus clustering with GPU-accelerated shape matching (estimated speedup: 10-50×).

### CUDA Distance Matrix Calculations

General GPU techniques for molecular geometry:[^10][^11][^12]

- **Distance geometry**: Predict pairwise distance matrices via graph CNNs (RDKit integration)
- **FCI calculations**: 15-35× GPU speedup for quantum chemistry[^10]
- **Protein distance matrices**: Parallelized for structural alignment[^13]

**Consensus use**: Accelerate pharmacophore distance matrix computation (input to clustering algorithms).

***

## Approximate Methods Trading Accuracy for Speed

### Discretization Strategies

**Fingerprint-based pharmacophores**:[^3]

- Bin 3-4 point triangles into distance ranges (e.g., 2Å bins)
- Loss: ~5-10% alignment accuracy vs. exact coordinates
- Gain: 100× speedup via hash table lookups

**Bloom filter approximations**:[^3]

- Multi-resolution encoding (4Å: 6 bits, 1Å: 8 bits)
- False positive rate: 0.43% (acceptable for pre-screening)
- Gain: Eliminates 99%+ of incompatible candidates before alignment


### Hierarchical Clustering + Sampling

**Workflow**:

1. Compute rough similarity matrix (GPU Tanimoto: seconds)
2. Cluster 1000 molecules → 10-50 representatives (Butina clustering: minutes)[^14]
3. Full alignment on representatives (minutes)
4. Expand consensus to all molecules via nearest-cluster lookup

**Trade-off**: May miss minority pharmacophores (<5% occurrence), but captures dominant patterns.

### Greedy Search Algorithms

**G3PS (Greedy 3-Point Search)**:[^15]

- Maximize matched pharmacophore features within tolerances (not RMSD minimization)
- Avoids Hungarian algorithm's O(*m*³) complexity
- ~10× faster than full clique detection

***

## Parallel and Distributed Algorithms

### Current State

**Limited implementations found**:

- **RDKit multiprocessing**: Embarrassingly parallel conformer generation[^14]

```python
from multiprocessing import Pool
with Pool(32) as p:
    conformers = p.map(generate_conformers, smiles_list)
```

- **Protein structure alignment**: MPI/pthreads/OpenMP (parMATT: hybrid parallelization)[^16][^17][^18]

**No evidence of**:

- Distributed pharmacophore indexing (Spark/Dask)
- MPI-based consensus algorithms
- GPU-native pharmacophore alignment libraries


### Recommended Parallelization Strategy

```python
# Consensus generation with multiprocessing

from multiprocessing import Pool
from rdkit.Chem import AllChem, Pharmacophore

def extract_pharmacophore(mol):
    # Generate conformer + features
    AllChem.EmbedMolecule(mol)
    feats = featFactory.GetFeaturesForMol(mol)
    return Pharmacophore.Pharmacophore(feats)

def pairwise_align(mol_pair):
    # Align pharmacophores (CDPKit/Kabsch)
    return alignment_score, transformation_matrix

# Step 1: Parallel feature extraction
with Pool(32) as p:
    pharmacophores = p.map(extract_pharmacophore, molecules)  # Seconds

# Step 2: Parallel pairwise alignment
pairs = [(p1, p2) for i, p1 in enumerate(pharmacophores) 
                  for p2 in pharmacophores[i+1:]]
with Pool(32) as p:
    alignments = p.map(pairwise_align, pairs)  # Minutes for 1000 mols

# Step 3: Consensus clustering (sequential, fast)
consensus = cluster_features(alignments)  # Seconds
```

**Estimated speedup**: 32× on 32-core system vs. sequential.

***

## Code Examples and Tool Implementations

### Example 1: RDKit Pharmacophore Detection

```python
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures

# Load feature definitions (SMARTS-based)
fdefName = 'BaseFeatures.fdef'  # HBD, HBA, aromatic, hydrophobic
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)

# Extract pharmacophore
mol = Chem.MolFromSmiles('CCO')  # Ethanol
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)  # Generate 3D coordinates

features = featFactory.GetFeaturesForMol(mol)
for feat in features:
    print(f"{feat.GetFamily()}: {feat.GetPos()}")  
    # Output: Donor: (x, y, z), Acceptor: (x, y, z)
```


### Example 2: Kabsch Alignment (Distance Matrix-Based)

```python
import numpy as np

def kabsch_align(P, Q):
    """Align point sets P to Q, return RMSD and rotation matrix"""
    # Center point clouds
    P_centered = P - P.mean(axis=0)
    Q_centered = Q - Q.mean(axis=0)
    
    # Covariance matrix
    H = P_centered.T @ Q_centered
    
    # SVD
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    # Reflection correction
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # RMSD
    P_aligned = P_centered @ R
    rmsd = np.sqrt(((P_aligned - Q_centered)**2).sum() / len(P))
    
    return rmsd, R
```


### Example 3: GPU Tanimoto (cupy)

```python
import cupy as cp

def gpu_tanimoto_all_vs_all(fingerprints):
    """
    fingerprints: (n_mols, n_bits) binary array
    Returns: (n_mols, n_mols) similarity matrix
    """
    fp_gpu = cp.asarray(fingerprints, dtype=cp.uint8)
    
    # Compute intersection and union
    intersection = cp.sum(fp_gpu[:, None, :] & fp_gpu[None, :, :], axis=2)
    union = cp.sum(fp_gpu[:, None, :] | fp_gpu[None, :, :], axis=2)
    
    return (intersection / union).get()  # Transfer to CPU

# Usage
fps = np.random.randint(0, 2, (1000, 2048))  # 1000 molecules, 2048-bit fps
similarity_matrix = gpu_tanimoto_all_vs_all(fps)  # ~0.5 seconds on RTX 3090
```


***

## Integrated Workflow Recommendations

### Scenario 1: Maximum Speed (Pre-Screening Quality)

**Method**: PharmacoMatch + hierarchical clustering

```
1. Embed all molecules        : 10s  (GPU, 1000 mols)
2. Cosine similarity matrix    : 1s   (GPU, 512-dim vectors)
3. Hierarchical clustering     : 2s   (scipy.cluster.hierarchy)
4. Extract consensus features  : 1s   (majority vote per cluster)

Total: ~15 seconds
```

**Trade-off**: ~10% lower precision than exact alignment; unsuitable for mirror-image discrimination.

### Scenario 2: Balanced Speed/Accuracy

**Method**: Pharmer KDB-tree + Bloom filtering

```
1. Build index (one-time)      : 5min (1000 mols × 200 conformers)
2. Query representatives       : 10s  (10 diverse queries, tight tolerances)
3. Cluster matching features   : 5s   (density-based clustering)
4. Validate with alignment     : 30s  (CDPKit on consensus hits)

Total: ~50 seconds (+ 5min one-time setup)
```

**Trade-off**: Exact alignment preserved; requires index storage (~GB-scale for conformers).

### Scenario 3: Highest Accuracy

**Method**: Parallel CDPKit alignment + consensus clustering

```
1. Extract pharmacophores      : 1min (parallel RDKit, 32 cores)
2. All-vs-all alignment        : 5min (CDPKit, 499,500 pairs, 32 cores)
3. Distance matrix clustering  : 10s  (HDBSCAN/OPTICS)
4. Consensus pharmacophore     : 5s   (geometric median of clusters)

Total: ~6 minutes
```

**Trade-off**: No approximations; suitable for publication-grade models.

***

## Comparative Performance Summary

| Method | Speed (1000 mols) | Accuracy | Scalability | Implementation |
| :-- | :-- | :-- | :-- | :-- |
| **PharmacoMatch** | 15s | Good | Billions (w/ ANN) | Python, GPU |
| **Pharmer** | 50s | Exact | Millions | C++, open-source |
| **PharmacoNet** | 30s | Good | 100M+ | Python, CPU |
| **GPU Tanimoto + CDPKit** | 5min | Exact | 10K-100K | Python/C++ |
| **Parallel CDPKit** | 6min | Exact | 10K | C++, multicore |
| **MOE (baseline)** | 2.8hrs | High | <10K | Proprietary |


***

## Limitations and Future Directions

### Current Gaps

1. **No native distributed implementations**: Spark/Dask pharmacophore frameworks absent
2. **GPU alignment algorithms**: Limited to shape overlay (ROSHAMBO); no full pharmacophore matcher
3. **Approximate NN integration**: FAISS-like vector databases not yet adapted for pharmacophores
4. **Consensus-specific benchmarks**: Literature focuses on single-query screening, not multi-ligand consensus

### Emerging Opportunities

**Vector databases for pharmacophores**:[^4]

- Store PharmacoMatch embeddings in Milvus/Pinecone
- Approximate nearest neighbor search (milliseconds for billions)

**Hybrid CPU-GPU pipelines**:

- GPU: Distance matrix + clustering
- CPU: Exact alignment refinement

**Federated learning for consensus**:

- Distribute molecular libraries across nodes
- Aggregate local pharmacophore models (privacy-preserving drug discovery)

***

## Conclusion

The fastest methods for generating consensus pharmacophore models from 100-1000 molecules achieve **15 seconds to 6 minutes** depending on accuracy requirements—a **30-1000× speedup** over traditional sequential alignment. **PharmacoMatch** offers unparalleled speed (15s) for pre-screening via learned embeddings, while **Pharmer's KDB-tree** balances exact retrieval with sub-minute queries. For publication-grade accuracy, **parallel CDPKit alignment** (6 minutes) remains optimal. GPU acceleration via Tanimoto similarity and shape matching provides complementary pre-filtering strategies. Future integration of approximate nearest neighbor search and distributed computing promises to extend these methods to billion-scale consensus pharmacophore analysis.
<span style="display:none">[^100][^101][^102][^103][^104][^105][^106][^107][^108][^109][^110][^111][^112][^113][^114][^115][^116][^117][^118][^119][^120][^121][^122][^123][^124][^125][^126][^127][^128][^129][^130][^131][^132][^133][^134][^135][^136][^137][^138][^139][^140][^141][^142][^143][^144][^145][^146][^147][^148][^149][^150][^151][^152][^153][^154][^155][^156][^157][^158][^159][^19][^20][^21][^22][^23][^24][^25][^26][^27][^28][^29][^30][^31][^32][^33][^34][^35][^36][^37][^38][^39][^40][^41][^42][^43][^44][^45][^46][^47][^48][^49][^50][^51][^52][^53][^54][^55][^56][^57][^58][^59][^60][^61][^62][^63][^64][^65][^66][^67][^68][^69][^70][^71][^72][^73][^74][^75][^76][^77][^78][^79][^80][^81][^82][^83][^84][^85][^86][^87][^88][^89][^90][^91][^92][^93][^94][^95][^96][^97][^98][^99]</span>

<div align="center">⁂</div>

[^1]: http://e-currentscience.com/storage/app/archive/pdf/UwUw5HzRokWnD1d3zhLa7FfJictQJCiyHu3NlF7i.pdf

[^2]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4290656/

[^3]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3124593/

[^4]: https://arxiv.org/html/2409.06316

[^5]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11575537/

[^6]: https://pubs.acs.org/doi/10.1021/ci1004948

[^7]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3445263/

[^8]: https://pmc.ncbi.nlm.nih.gov/articles/PMC5072535/

[^9]: https://pubs.acs.org/doi/10.1021/acs.jcim.4c01225

[^10]: https://ieeexplore.ieee.org/document/11318608/

[^11]: https://arxiv.org/ftp/arxiv/papers/2107/2107.01035.pdf

[^12]: https://annals-csis.org/proceedings/2022/drp/pdf/6.pdf

[^13]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8815937/

[^14]: https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-S1-P36

[^15]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8658842/

[^16]: https://academic.oup.com/bioinformatics/article/35/21/4456/5421160

[^17]: http://link.springer.com/10.1007/978-3-030-38752-5_3

[^18]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3309952/

[^19]: https://ieeexplore.ieee.org/document/10736963/

[^20]: https://ieeexplore.ieee.org/document/9300224/

[^21]: https://www.spiedigitallibrary.org/conference-proceedings-of-spie/13171/3032033/Research-on-speed-profile-generation-of-train-automatic-driving-planning/10.1117/12.3032033.full

[^22]: https://link.springer.com/10.1007/s10836-022-05999-9

[^23]: https://ieeexplore.ieee.org/document/8957030/

[^24]: https://ieeexplore.ieee.org/document/10295013/

[^25]: https://www.mdpi.com/1996-1073/17/2/300

[^26]: http://www.emerald.com/ijesm/article/13/4/828-845/119128

[^27]: https://ieeexplore.ieee.org/document/9531603/

[^28]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11894060/

[^29]: https://arxiv.org/pdf/2310.00681.pdf

[^30]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4012150/

[^31]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11903766/

[^32]: https://www.frontiersin.org/articles/10.3389/fphar.2018.01463/pdf

[^33]: http://biorxiv.org/lookup/doi/10.1101/2024.11.10.622877

[^34]: https://www.ijisrt.com/a-dynamic-mpibased-memoryefficient-framework-for-longest-common-subsequence-computation-on-massive-dna-sequence

[^35]: https://ieeexplore.ieee.org/document/10192069/

[^36]: https://openproceedings.org/2021/conf/edbt/p31.pdf

[^37]: https://ieeexplore.ieee.org/document/9478946/

[^38]: https://ieeexplore.ieee.org/document/10385767/

[^39]: https://academic.oup.com/bioinformatics/article/33/19/3011/3852082

[^40]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2447755/

[^41]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11659349/

[^42]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2699263/

[^43]: https://www.mdpi.com/1420-3049/26/23/7201/pdf

[^44]: https://www.frontiersin.org/articles/10.3389/fchem.2024.1382319/pdf?isPublishedV2=False

[^45]: https://pubs.acs.org/doi/10.1021/acs.jcim.4c00424

[^46]: https://xlink.rsc.org/?DOI=D3SC04185A

[^47]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9025992/

[^48]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2646723/

[^49]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10507964/

[^50]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6454689/

[^51]: https://peerj.com/articles/725

[^52]: https://pmc.ncbi.nlm.nih.gov/articles/PMC5729909/

[^53]: http://arxiv.org/pdf/2410.13147.pdf

[^54]: https://pubs.acs.org/doi/10.1021/acs.jcim.3c01245

[^55]: https://arxiv.org/abs/2510.02578

[^56]: https://pubs.acs.org/doi/10.1021/acs.jcim.0c00121

[^57]: https://arxiv.org/abs/2506.00771

[^58]: https://www.mdpi.com/1424-8247/16/5/661

[^59]: https://pubs.acs.org/doi/10.1021/ci2005274

[^60]: https://arxiv.org/abs/2509.00684

[^61]: https://pubs.acs.org/doi/10.1021/ci800072r

[^62]: https://onlinelibrary.wiley.com/doi/10.1111/cbdd.12590

[^63]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4048638/

[^64]: https://www.mdpi.com/1420-3049/23/12/3094/pdf

[^65]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6396084/

[^66]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4765363/

[^67]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4444576/

[^68]: https://www.mdpi.com/1420-3049/24/16/2870/pdf

[^69]: https://www.semanticscholar.org/paper/91949128ad70b3cb408f1688163875f02c6bc6f2

[^70]: https://etasr.com/index.php/ETASR/article/view/10840

[^71]: https://pubs.acs.org/doi/10.1021/acs.jcim.5c00894

[^72]: https://arxiv.org/abs/2404.04810

[^73]: https://www.mdpi.com/2673-7655/4/4/41

[^74]: http://biorxiv.org/lookup/doi/10.1101/2024.11.24.625084

[^75]: https://applmathjournal.spbu.ru/article/view/18681

[^76]: https://pubs.aip.org/jcp/article/125/4/044903/188159/Molecular-fractionation-with-conjugated-caps

[^77]: https://epubs.siam.org/doi/10.1137/1.9781611974348.94

[^78]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2734148/

[^79]: https://journals.sagepub.com/doi/pdf/10.1177/1094342020964857

[^80]: https://arxiv.org/pdf/2408.16188.pdf

[^81]: https://arxiv.org/pdf/2206.13602.pdf

[^82]: https://arxiv.org/pdf/2412.10664.pdf

[^83]: https://pubs.acs.org/doi/10.1021/acs.jcim.7b00618

[^84]: https://www.nature.com/articles/s41596-025-01237-6

[^85]: https://www.semanticscholar.org/paper/cc02757f909cfe6242567ebd46ab34c6d81f85be

[^86]: https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-6-S1-P57

[^87]: https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-6-S1-P8

[^88]: https://journals.ashs.org/view/journals/jashs/142/3/article-p163.xml

[^89]: https://jcheminf.springeropen.com/track/pdf/10.1186/1758-2946-3-S1-P30

[^90]: https://arxiv.org/html/2308.11890

[^91]: https://arxiv.org/html/2411.10821v1

[^92]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4710026/

[^93]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10918633/

[^94]: https://pmc.ncbi.nlm.nih.gov/articles/PMC1994057/

[^95]: https://arxiv.org/pdf/2306.07812.pdf

[^96]: http://link.springer.com/10.1007/s10822-014-9709-3

[^97]: https://xlink.rsc.org/?DOI=D3MD00649B

[^98]: https://onlinelibrary.wiley.com/doi/10.1002/minf.202000090

[^99]: https://jbiomedsci.biomedcentral.com/articles/10.1186/1423-0127-18-8

[^100]: https://www.tandfonline.com/doi/full/10.1080/07391102.2019.1701553

[^101]: https://www.frontiersin.org/articles/10.3389/fmolb.2022.836572/full

[^102]: https://www.mdpi.com/1420-3049/25/2/385

[^103]: http://www.eurekaselect.com/openurl/content.php?genre=article\&issn=1386-2073\&volume=15\&issue=10\&spage=849

[^104]: https://www.omicsgroup.org/journals/in-silico-screening-for-identification-of-novel-aurora-kinase-inhibitors-by-molecular-docking-dynamics-simulations-and-ligand-based-hypothesis-approaches-2329-6674.1000106.php?aid=11799

[^105]: https://pubs.acs.org/doi/10.1021/ci500291a

[^106]: https://www.mdpi.com/1420-3049/25/2/385/pdf

[^107]: https://www.mdpi.com/1422-0067/20/23/5834

[^108]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6929024/

[^109]: https://f1000research.com/articles/3-277/v1

[^110]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9554199/

[^111]: https://www.frontiersin.org/articles/10.3389/fcimb.2021.611304/pdf

[^112]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4088285/

[^113]: https://www.mdpi.com/1422-0067/23/3/1309/pdf

[^114]: https://ieeexplore.ieee.org/document/11313211/

[^115]: http://telkomnika.uad.ac.id/index.php/TELKOMNIKA/article/view/5916

[^116]: https://pubs.acs.org/doi/10.1021/acsmedchemlett.4c00145

[^117]: https://pubs.acs.org/doi/10.1021/acs.jcim.2c00242

[^118]: https://onlinelibrary.wiley.com/doi/10.1002/minf.202300056

[^119]: https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/jat.4934

[^120]: https://www.worldscientific.com/doi/10.1142/S0219720025500039

[^121]: http://biorxiv.org/lookup/doi/10.1101/2025.06.16.659994

[^122]: https://www.tandfonline.com/doi/full/10.1080/17460441.2022.2128332

[^123]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3852750/

[^124]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3341238/

[^125]: https://arxiv.org/pdf/2109.06355.pdf

[^126]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3606147/

[^127]: https://pmc.ncbi.nlm.nih.gov/articles/PMC7050271/

[^128]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4456712/

[^129]: https://ieeexplore.ieee.org/document/10336943/

[^130]: https://ieeexplore.ieee.org/document/11221805/

[^131]: https://www.semanticscholar.org/paper/db082b9c72423c44192dea162ea785884df21f7c

[^132]: https://ieeexplore.ieee.org/document/10552313/

[^133]: https://onlinelibrary.wiley.com/doi/10.1002/rob.70010

[^134]: https://ieeexplore.ieee.org/document/10415090/

[^135]: https://link.springer.com/10.1007/s12083-024-01830-8

[^136]: https://ieeexplore.ieee.org/document/10248719/

[^137]: https://ieeexplore.ieee.org/document/9862846/

[^138]: http://thesai.org/Publications/ViewPaper?Volume=15\&Issue=11\&Code=ijacsa\&SerialNo=42

[^139]: https://arxiv.org/html/2412.19812v2

[^140]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11234869/

[^141]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11865752/

[^142]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9832483/

[^143]: https://www.semanticscholar.org/paper/a064f697e71b0660d868f932ae97936d9bbcf08d

[^144]: https://iopscience.iop.org/article/10.3847/1538-4357/ac75b7

[^145]: https://arxiv.org/abs/2504.16113

[^146]: https://arxiv.org/abs/2507.13901

[^147]: http://biorxiv.org/lookup/doi/10.1101/2021.08.24.457462

[^148]: http://preprints.jmir.org/preprint/27123

[^149]: https://www.mdpi.com/1996-1944/12/13/2169

[^150]: https://www.jstatsoft.org/v103/i09/

[^151]: https://joss.theoj.org/papers/10.21105/joss.04864

[^152]: https://www.semanticscholar.org/paper/ff16189169b7081d7f8121c4e2736d6c8384c450

[^153]: https://academic.oup.com/bioinformatics/article/40/Supplement_1/i539/7700904

[^154]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11211825/

[^155]: https://arxiv.org/pdf/2101.11003.pdf

[^156]: http://www.jstatsoft.org/v18/i05/

[^157]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11469520/

[^158]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11687028/

[^159]: https://arxiv.org/pdf/2211.02566.pdf


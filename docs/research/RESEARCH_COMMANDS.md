# Research Commands & Prompts for Consensus Pharmacophore Literature

## 🎯 Research Goals

1. **Speed**: Which algorithms are fastest for clustering/consensus generation?
2. **Accuracy**: Which methods give best virtual screening enrichment?
3. **Parameters**: How to optimize tolerance, occurrence threshold, linkage type?
4. **Benchmarks**: What datasets/metrics are standard for evaluation?

---

## 📚 Literature Search Commands

### Using PubMed (via CLI or API)

```bash
# Install entrez-direct tools (NCBI E-utilities)
conda install -c bioconda entrez-direct

# Search for consensus pharmacophore papers
esearch -db pubmed -query "consensus pharmacophore" | efetch -format abstract > papers/pubmed_consensus_pharm.txt

# More specific queries
esearch -db pubmed -query "(consensus pharmacophore) AND (clustering OR algorithm)" | efetch -format abstract > papers/clustering_methods.txt

esearch -db pubmed -query "(pharmacophore) AND (virtual screening) AND (benchmark OR enrichment)" | efetch -format abstract > papers/benchmarking.txt

esearch -db pubmed -query "(spatial clustering) AND (3D molecular features)" | efetch -format abstract > papers/spatial_clustering.txt

# Get full metadata with PMIDs
esearch -db pubmed -query "consensus pharmacophore optimization" | efetch -format xml > papers/consensus_meta.xml
```

### Using Google Scholar (via serpapi or manual)

```bash
# Install scholarly Python library
pip install scholarly

# Save this as search_scholar.py
python << 'EOF'
from scholarly import scholarly
import json

queries = [
    "consensus pharmacophore virtual screening",
    "pharmacophore clustering algorithms",
    "spatial feature clustering 3D molecules",
    "hierarchical clustering pharmacophore",
    "DBSCAN pharmacophore",
    "pharmacophore benchmark dataset DUD-E"
]

for query in queries:
    print(f"\n=== Searching: {query} ===")
    search_query = scholarly.search_pubs(query)
    
    results = []
    for i, pub in enumerate(search_query):
        if i >= 10:  # Top 10 per query
            break
        results.append({
            'title': pub['bib']['title'],
            'author': pub['bib'].get('author', 'Unknown'),
            'year': pub['bib'].get('pub_year', 'N/A'),
            'abstract': pub['bib'].get('abstract', 'N/A')
        })
        print(f"{i+1}. {pub['bib']['title']} ({pub['bib'].get('pub_year', 'N/A')})")
    
    with open(f"papers/scholar_{query.replace(' ', '_')[:30]}.json", 'w') as f:
        json.dump(results, f, indent=2)

EOF
```

### Using ArXiv (for preprints)

```bash
# Install arxiv Python library
pip install arxiv

# Search ArXiv
python << 'EOF'
import arxiv

search = arxiv.Search(
    query="pharmacophore clustering OR molecular feature clustering",
    max_results=20,
    sort_by=arxiv.SortCriterion.Relevance
)

with open('papers/arxiv_results.txt', 'w') as f:
    for result in search.results():
        f.write(f"Title: {result.title}\n")
        f.write(f"Authors: {', '.join(a.name for a in result.authors)}\n")
        f.write(f"Published: {result.published}\n")
        f.write(f"PDF: {result.pdf_url}\n")
        f.write(f"Abstract: {result.summary}\n")
        f.write("-" * 80 + "\n\n")
        print(f"Found: {result.title}")

EOF
```

---

## 🤖 AI-Powered Research Prompts

### For Literature Review Skill (Claude)

```bash
# Use the literature-review skill with this prompt:

@literature-review

"I need a comprehensive literature review on consensus pharmacophore generation methods, 
focusing on:

1. **Clustering Algorithms**: Compare hierarchical, DBSCAN, k-means, grid-based, and 
   spatial hashing methods for 3D molecular feature clustering
   
2. **Performance Benchmarks**: What are the standard datasets (DUD-E, LIT-PCBA, etc.) 
   and metrics (enrichment factor, AUC-ROC) for evaluating pharmacophore models?
   
3. **Parameter Optimization**: How do tolerance (spatial radius), occurrence threshold 
   (consensus frequency), and linkage type affect screening accuracy and speed?
   
4. **Speed Optimization**: What are the fastest methods for generating consensus models 
   from 10-1000 aligned molecules?

Please search PubMed, Google Scholar, and ChemRxiv. Prioritize papers from 2018-2026 
with experimental validation on benchmark datasets."
```

### For Perplexity/ChatGPT Research Mode

```
Prompt 1: Algorithm Comparison
-------------------------------
"Compare clustering algorithms for 3D pharmacophore feature consensus generation:

- Hierarchical (agglomerative/divisive)
- DBSCAN
- K-means
- Grid-based spatial hashing
- KD-tree nearest neighbor

For each, provide:
1. Time complexity
2. Determinism (same input → same output?)
3. Parameter sensitivity
4. Best use case (small vs large datasets)
5. Recent papers (2020-2026) using this approach

Focus on molecular informatics and cheminformatics literature."
```

```
Prompt 2: Benchmark Standards
------------------------------
"What are the gold-standard benchmark datasets and evaluation metrics for 
validating pharmacophore-based virtual screening?

Include:
- Dataset names (DUD-E, MUV, LIT-PCBA, etc.)
- Metrics (EF, AUC, precision@1%, etc.)
- Typical parameter ranges (tolerance 1.0-2.0 Å, threshold 0.5-1.0)
- Recent papers comparing consensus vs single-molecule pharmacophores

Provide download links and recommended evaluation protocols."
```

```
Prompt 3: Speed Optimization
-----------------------------
"What are the fastest methods for generating consensus pharmacophore models from 
100-1000 aligned molecules?

Consider:
- Distance matrix calculation optimizations
- GPU acceleration techniques
- Spatial indexing data structures
- Parallel/distributed algorithms
- Approximate methods trading accuracy for speed

Provide code examples, libraries, or tools implementing these methods."
```

---

## 📥 Paper Download & Organization

### Automated Paper Retrieval

```bash
# Create download script
cat > download_papers.sh << 'SCRIPT'
#!/bin/bash

PAPERS_DIR="docs/research/papers"
mkdir -p "$PAPERS_DIR"

# Array of DOIs or PMIDs
declare -a DOIS=(
    "10.1021/acs.jcim.8b00928"  # Example: Replace with actual DOIs
    "10.1186/s13321-019-0393-x"
    # Add more DOIs here
)

# Download using Sci-Hub (use legally!)
for doi in "${DOIS[@]}"; do
    echo "Downloading $doi..."
    # Use curl or wget to fetch from institutional access or legal sources
    # Example: curl "https://doi.org/$doi" -o "$PAPERS_DIR/${doi//\//_}.pdf"
done

echo "Download complete!"
SCRIPT

chmod +x download_papers.sh
```

### Citation Management

```bash
# Install citation tools
pip install pybtex bibtexparser

# Create BibTeX database
cat > docs/research/papers/consensus_pharm.bib << 'BIBTEX'
@article{example2023,
    title={Consensus Pharmacophore Optimization},
    author={Smith, J. and Doe, A.},
    journal={Journal of Chemical Information and Modeling},
    year={2023},
    volume={63},
    pages={1234-1245},
    doi={10.1021/example}
}

% Add more entries as you find papers
BIBTEX
```

---

## 🔬 Benchmark Dataset Downloads

```bash
# DUD-E (Database of Useful Decoys Enhanced)
mkdir -p benchmarks/DUD-E
wget -r -np -nH --cut-dirs=2 -R "index.html*" \
    http://dude.docking.org/targets/ \
    -P benchmarks/DUD-E/

# LIT-PCBA (from PubChem)
# Visit: ftp://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/

# MUV (Maximum Unbiased Validation)
wget http://www.bioinf.jku.at/research/MUV/MUV_dataset.tar.gz \
    -P benchmarks/
tar -xzf benchmarks/MUV_dataset.tar.gz -C benchmarks/
```

---

## 📊 Research Tracking Template

Create a research log:

```bash
cat > docs/research/RESEARCH_LOG.md << 'LOG'
# Research Log - Consensus Pharmacophore Optimization

## 2026-01-26

### Papers Read
- [ ] Paper 1: [Title] - Key finding: ...
- [ ] Paper 2: [Title] - Key finding: ...

### Experiments Run
- [ ] Benchmark hierarchical vs DBSCAN on DUD-E subset
- [ ] Test tolerance range 0.8-2.0 Å with fixed threshold
- [ ] Profile clustering time for 50, 100, 500, 1000 molecules

### Ideas to Test
- [ ] Grid-based pre-filtering before hierarchical clustering
- [ ] GPU-accelerated distance matrix computation
- [ ] Feature-type-specific tolerance values

### Questions to Answer
- [ ] What is the optimal linkage method: average, complete, ward?
- [ ] Can we predict optimal tolerance from molecule diversity metrics?
- [ ] How does consensus model size (# features) affect screening speed?

LOG
```

---

## 🎓 Recommended Search Terms

### Core Terms
- "consensus pharmacophore"
- "pharmacophore clustering"
- "3D pharmacophore alignment"
- "spatial feature clustering"
- "molecular feature consensus"

### Algorithm-Specific
- "hierarchical clustering pharmacophore"
- "DBSCAN molecular features"
- "agglomerative clustering ligands"
- "spatial hashing 3D molecules"

### Application-Focused
- "pharmacophore virtual screening benchmark"
- "ligand-based screening optimization"
- "pharmacophore enrichment factor"
- "DUD-E pharmacophore evaluation"

### Speed/Optimization
- "fast pharmacophore generation"
- "GPU molecular clustering"
- "parallel pharmacophore algorithm"
- "scalable ligand alignment"

---

## 🚀 Quick Start Workflow

```bash
# 1. Set up research environment
cd /home/dodo/Documents/projects/pharmacophore-toolkit/docs/research

# 2. Search PubMed
esearch -db pubmed -query "consensus pharmacophore clustering" | efetch -format abstract > papers/initial_search.txt

# 3. Review abstracts and note DOIs
grep -i "doi\|pmid" papers/initial_search.txt > papers/dois.txt

# 4. Download key papers (use institutional access)
# ... manual step ...

# 5. Extract key insights
# Read papers and update ideas/OPTIMIZATION_BRAINSTORM.md

# 6. Track progress
echo "- [x] Initial PubMed search completed" >> RESEARCH_LOG.md
```

---

## 📧 Email Alerts for New Papers

```bash
# Set up PubMed alert via MyNCBI
# Visit: https://www.ncbi.nlm.nih.gov/myncbi/
# Create search: "consensus pharmacophore"
# Enable email alerts: Weekly digest

# Alternative: RSS feed
# Subscribe to: https://pubmed.ncbi.nlm.nih.gov/rss/search/...
```

---

## Next Steps

1. Run initial PubMed search
2. Read top 10 papers based on citations/relevance
3. Extract algorithm details and benchmark results
4. Update OPTIMIZATION_BRAINSTORM.md with findings
5. Design experiments to test top 3 approaches
6. Document results in benchmarks/ folder

Good luck with your research! 🚀
LOG

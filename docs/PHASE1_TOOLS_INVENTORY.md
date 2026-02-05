# Phase 1: Tools, Skills, Plugins & Agents Inventory

**Project**: CCR2 Pharmacophore Virtual Screening
**Phase**: Literature Research & Algorithm Design
**Date**: 2026-01-23

---

## 1. MCP Servers (Active)

### 1.1 Scientific Literature & Data

| Server | Tools | Purpose | Status |
|--------|-------|---------|--------|
| **PubMed** | `search_articles`, `get_article_metadata`, `find_related_articles`, `get_full_text_article`, `convert_article_ids` | Biomedical literature search | ✅ Active |
| **bioRxiv** | `search_preprints`, `get_preprint`, `get_categories`, `search_published_preprints`, `search_by_funder` | Preprints in biology/medicine | ✅ Active |
| **ChEMBL** | `compound_search`, `target_search`, `get_bioactivity`, `get_mechanism`, `drug_search`, `get_admet` | Drug/target data | ✅ Active |
| **Synapse** | `get_entity`, `search_synapse`, `get_entity_children` | Research data platform | ✅ Active |
| **Context7** | `resolve-library-id`, `query-docs` | Library documentation | ✅ Active |

### 1.2 Code & Development

| Server | Tools | Purpose | Status |
|--------|-------|---------|--------|
| **Serena** | `find_symbol`, `replace_symbol_body`, `search_for_pattern`, `get_symbols_overview`, etc. | Semantic code editing | ✅ Active |
| **Claude-mem** | `search`, `timeline`, `get_observations` | Context/memory management | ✅ Active |
| **IDE** | `getDiagnostics`, `executeCode` | Jupyter/VS Code integration | ✅ Active |

### 1.3 External APIs (Accessible but not MCP)

| Service | Access Method | Purpose | Status |
|---------|---------------|---------|--------|
| **Zotero Local API** | HTTP `localhost:23119` | Personal library management | ✅ Accessible (empty) |
| **Web Search** | Built-in tool | General web search | ✅ Active |
| **WebFetch** | Built-in tool | Fetch/analyze web content | ✅ Active |

---

## 2. Skills (Invoke via `/skill-name`)

### 2.1 Research & Analysis Skills

| Skill | Trigger | Use Case for Phase 1 |
|-------|---------|---------------------|
| `/literature-review` | Systematic review | **PRIMARY** - Synthesize pharmacophore papers |
| `/scientific-reviewer` | Document analysis | Evaluate methodology quality |
| `/brainstorming` | Ideation | Explore algorithm approaches |
| `/design-of-experiments` | DOE guidance | Plan validation experiments |

### 2.2 Development Skills

| Skill | Trigger | Use Case for Phase 1 |
|-------|---------|---------------------|
| `/python-expert` | Python optimization | Algorithm implementation |
| `/data-scientist` | Data analysis | Metric analysis |
| `/ml-engineer` | ML pipelines | Surrogate model building |
| `/code-review` | Code quality | Review implementations |
| `/debugger` | Error fixing | Debug issues |

### 2.3 Project Management Skills

| Skill | Trigger | Use Case |
|-------|---------|----------|
| `/start` | Task orchestration | Plan complex workflows |
| `/commit` | Git commits | Version control |
| `/create-pr` | Pull requests | Code submission |
| `/todo` | Task management | Track progress |

---

## 3. Specialized Agents (via Task tool)

| Agent Type | Purpose | When to Use |
|------------|---------|-------------|
| **Explore** | Codebase exploration | Understanding existing code |
| **Plan** | Implementation planning | Designing architecture |
| **Bash** | Command execution | Git, build, run scripts |
| **general-purpose** | Multi-step research | Complex searches |
| **code-reviewer** | Code review | After implementations |
| **debugger** | Error investigation | When things break |
| **python-expert** | Python expertise | Complex Python patterns |
| **data-scientist** | Data analysis | Statistical analysis |
| **ml-engineer** | ML implementation | Model building |
| **ccr2-ligand-designer** | CCR2 drug design | **SPECIALIZED** for this project! |
| **ccr2-drugex-scientist** | CCR2 DrugEx workflows | Ligand generation |

---

## 4. Phase 1 Tool Usage Plan

### 4.1 Literature Search Strategy

```
┌─────────────────────────────────────────────────────────────────┐
│                    LITERATURE SEARCH WORKFLOW                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                   │
│   ┌─────────┐    ┌─────────┐    ┌─────────┐    ┌─────────┐      │
│   │ PubMed  │    │ bioRxiv │    │ ChEMBL  │    │   Web   │      │
│   │ Search  │    │ Search  │    │ Target  │    │ Search  │      │
│   └────┬────┘    └────┬────┘    └────┬────┘    └────┬────┘      │
│        │              │              │              │            │
│        └──────────────┴──────────────┴──────────────┘            │
│                            │                                      │
│                            ▼                                      │
│                   ┌─────────────────┐                            │
│                   │  Filter & Rank  │                            │
│                   │  by Relevance   │                            │
│                   └────────┬────────┘                            │
│                            │                                      │
│                            ▼                                      │
│                   ┌─────────────────┐                            │
│                   │ Get Full Text   │                            │
│                   │ (PMC, WebFetch) │                            │
│                   └────────┬────────┘                            │
│                            │                                      │
│                            ▼                                      │
│                   ┌─────────────────┐                            │
│                   │ Literature DB   │                            │
│                   │ (Markdown Doc)  │                            │
│                   └─────────────────┘                            │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
```

### 4.2 Search Queries by Topic

| Topic | Tool | Query |
|-------|------|-------|
| **Pharmacophore Methods** | PubMed | "pharmacophore virtual screening" AND "enrichment" |
| **Shape Matching** | PubMed | "shape-based virtual screening" AND "Tanimoto" |
| **CCR2 Modulators** | PubMed + ChEMBL | "CCR2" AND "allosteric" AND "modulator" |
| **Optimization Methods** | PubMed | "Bayesian optimization" AND "drug discovery" |
| **ML Pharmacophore** | bioRxiv | pharmacophore + machine learning (recent) |
| **Validation Metrics** | PubMed | "ROC AUC" AND "enrichment factor" AND "virtual screening" |

### 4.3 Expected Paper Categories

1. **Core Methods** (10-15 papers)
   - Consensus pharmacophore generation
   - Shape-based screening (ROCS, CDPKit, RDKit)
   - TanimotoCombo scoring

2. **CCR2 Target** (5-10 papers)
   - CCR2 crystal structures
   - Allosteric binding site characterization
   - Known CCR2 modulators

3. **Validation & Metrics** (5-10 papers)
   - ROC-AUC interpretation
   - Enrichment factor calculations
   - BEDROC methodology

4. **Optimization** (5-10 papers)
   - Multi-objective optimization
   - Pareto optimization in drug discovery
   - Active learning for VS

5. **Recent ML Advances** (5-10 papers)
   - PharmacoNet, PharmacoMatch
   - Deep learning pharmacophores
   - Neural scoring functions

---

## 5. Tool Capabilities Matrix

### 5.1 PubMed MCP Tools

| Tool | Parameters | Returns | Rate Limit |
|------|------------|---------|------------|
| `search_articles` | query, max_results, date_from/to, sort | PMIDs, titles, abstracts | ~10/sec |
| `get_article_metadata` | pmids[] | Full metadata, authors, DOIs | ~3/sec |
| `find_related_articles` | pmids[], link_type | Related PMIDs | ~3/sec |
| `get_full_text_article` | pmc_ids[] | Full text (PMC only) | ~1/sec |
| `convert_article_ids` | ids[], id_type | PMID↔PMCID↔DOI | ~10/sec |

### 5.2 bioRxiv MCP Tools

| Tool | Parameters | Returns | Notes |
|------|------------|---------|-------|
| `search_preprints` | date_from/to, category, limit | DOI, title, abstract | No keyword search! |
| `get_preprint` | doi | Full metadata, PDF URL | Use after search |
| `get_categories` | - | List of 27 categories | Use for filtering |
| `search_published_preprints` | date range, publisher | Published preprints | Track peer review |

### 5.3 ChEMBL MCP Tools

| Tool | Parameters | Returns | Use Case |
|------|------------|---------|----------|
| `target_search` | target_name, gene_symbol | Target IDs, info | Find CCR2 |
| `compound_search` | name, chembl_id, smiles | Compound data | Find known ligands |
| `get_bioactivity` | molecule/target_chembl_id | IC50, Ki, etc. | Activity data |
| `get_mechanism` | molecule_chembl_id | MoA, targets | Drug mechanisms |
| `drug_search` | indication | Approved drugs | Clinical compounds |

---

## 6. Recommended Phase 1 Workflow

### Step 1: Target Information (ChEMBL)
```python
# Get CCR2 target info
target_search(target_name="CCR2", organism="Homo sapiens")
# Get known CCR2 modulators
get_bioactivity(target_chembl_id="CCR2_ID", min_pchembl=6)
```

### Step 2: Core Literature (PubMed)
```python
# Pharmacophore methods
search_articles("pharmacophore virtual screening consensus", max_results=50)
# Shape screening
search_articles("shape-based screening TanimotoCombo ROCS", max_results=30)
# CCR2 specific
search_articles("CCR2 allosteric modulator crystal structure", max_results=20)
```

### Step 3: Recent Advances (bioRxiv)
```python
# Last 6 months preprints
search_preprints(recent_days=180, category="bioinformatics", limit=50)
search_preprints(recent_days=180, category="pharmacology and toxicology", limit=50)
```

### Step 4: Synthesize with Skills
```
/literature-review → Systematic synthesis
/scientific-reviewer → Methodology evaluation
/brainstorming → Algorithm design ideas
```

### Step 5: Document & Plan
- Create `docs/LITERATURE_DATABASE.md`
- Update `docs/CCR2_SCREENING_PLAN.md`
- Create implementation TODOs

---

## 7. Quick Reference Commands

### Literature Search
```bash
# PubMed search
mcp__plugin_pubmed_PubMed__search_articles(query="...", max_results=50)

# Get full text from PMC
mcp__plugin_pubmed_PubMed__get_full_text_article(pmc_ids=["PMC..."])

# bioRxiv recent
mcp__plugin_biorxiv_bioRxiv__search_preprints(recent_days=90, category="biochemistry")
```

### ChEMBL Queries
```bash
# Find CCR2
mcp__plugin_chembl_ChEMBL__target_search(gene_symbol="CCR2")

# Get CCR2 bioactivity
mcp__plugin_chembl_ChEMBL__get_bioactivity(target_chembl_id="CHEMBL...")
```

### Context7 Documentation
```bash
# Get library docs
mcp__plugin_context7_context7__resolve-library-id(libraryName="rdkit", query="shape matching")
mcp__plugin_context7_context7__query-docs(libraryId="/rdkit/rdkit", query="conformer generation")
```

---

*Document created: 2026-01-23*
*Status: Phase 1 - Literature Research*

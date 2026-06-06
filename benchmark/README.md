# `benchmark/` — Methods Comparison Workspace

**Purpose:** one place where (a) the **documented literature** (reproducible 5×-search stable core) lives, and (b) **head-to-head experiments** compare the **published methods** against **our proposed solution** (S3 discrimination-weighted color + cross-domain algorithm transfers) on the **same CCR2 actives/decoys benchmark**.

## The comparison question
> On a held-out CCR2 actives-vs-decoys enrichment test (BEDROC α=20, ROC-AUC, EF1% — the same metrics as the Tier-1 test), how does **our solution** rank against the established literature methods, and which borrowed-algorithm transfer (if any) beats them?

## Layout
```
benchmark/
├── README.md                     # this file
├── PLAN.md                       # the implementation plan (bite-sized tasks)
├── papers/
│   ├── REGISTRY.md               # 64 stable-core papers (≥4/5 runs), grouped by method family; ★ = has code
│   ├── stable_core_5x.json       # machine-readable source
│   └── references.bib            # BibTeX
├── methods/
│   └── external_methods.md       # the methods we compare + code availability + how to invoke
├── experiments/
│   └── compare_methods.py        # pluggable scorer-registry harness (extends ../experiments/s3_enrichment_test.py)
└── results/                      # per-method enrichment tables + the master comparison
```

## Our solution (the baseline-to-beat / the thing being validated)
1. **S3 discrimination-weighted color** — proven in Tier-1 (held-out AUC 0.64→0.87, BEDROC 0.18→0.60 vs equal-weight). See `../experiments/S3_ENRICHMENT_TIER1_FINDINGS.md`.
2. **Cross-domain transfers** (proposed) — TEASER++ color-aware pose seed, information-weighted fusion, Fused-GW pose-free prefilter. See `../docs/research/CROSSDOMAIN_ALGORITHM_TRANSFERS.md`.

## Methods to compare (benchmarkable = has open-source code)
From `papers/REGISTRY.md` (★ rows): **SHAFTS**, **OTMol** (molecular OT), **gWEGA** (GPU shape), **Fused-GW** (POT), **ROCS-family** (already wired via RDKit/CDPKit/OpenEye), **Differential Multimolecule Fingerprint** (Hutter — the literature implementation of our S3 idea), **ShEPhERD**, **G3PS**. Plus the registration methods (**TEASER++**, **Go-ICP**) as *pose-seed* add-ons rather than standalone scorers.

## How everything connects
- The **literature** (`papers/`) is the reproducible evidence base (5× search, Jaccard 0.823) — see `../docs/research/LITERATURE_REPRODUCIBILITY_5X.md`.
- The **experiments** (`experiments/`) plug each method into one harness and score the **same held-out CCR2 split**, so every comparison is apples-to-apples with our Tier-1 result.
- The **decision**: any method must beat our S3 baseline on held-out BEDROC with a bootstrap CI excluding 0 to be "worth adopting."

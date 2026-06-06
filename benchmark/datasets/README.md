# Benchmark Dataset Suite — for selecting the correct pharmacophore optimizer

A **multi-target, decoy-bias-controlled, scaffold-split** suite so optimizer selection is judged on *generalization*, not on exploiting one biased CCR2 set. (See `SPEC.md` for the uniform `load()` interface; this file is the suite map + methodology.) Built 2026-06-06 by 6 parallel agents (5 in isolated worktrees + 1 scite/literature); all DOIs Crossref-verified.

## Datasets (integrated & verified)
| dir | source | decoy type | counts (verified) | role |
|-----|--------|-----------|-------------------|------|
| `muv/` | MUV (Rohrer-Baumann 2009) | max-unbiased (sphere-exclusion) | 17 targets, ~30 act / ~15k dec each | **unbiased gold** |
| `litpcba/` | LIT-PCBA (Tran-Nguyen 2020) | **experimental inactives** | targets incl. ADRB2, OPRK1 (GPCR) | **unbiased gold** |
| `ccr2_mubd/` | Xia 2018 chemokine MUBD | max-unbiased | **60 act / 2340 dec** (downloaded) | **CCR2-specific unbiased** |
| `created/` | ChEMBL REST + RDKit (this repo) | created-muv-unbiased (max-Tc<0.35) | CCR2 1086/2852 · CCR5 1577/2267 · CXCR4 527/2016 | **self-contained, reproducible** |
| `dude/` | DUD-E (Mysinger 2012) | property-matched (**biased**) | adrb2, aa2ar (GPCR) | **bias control (diagnostic)** |
| `ccr2_project/` | existing project set | property-matched (bias-suspect) | 74 / ~498 | current baseline (compare) |

`fetch.py`/`build.py` are re-runnable; large raw files are gitignored (`raw/`, `*.csv.gz`, `*.tgz`, `*.zip`) — only fetchers, loaders, small samples, and `meta.json` are committed. **ADRB2 appears in both LIT-PCBA (unbiased) and DUD-E (biased)** → the shared-target bias-control pair.

## Tooling (TDD, 13/13 green)
- `split.py` — `scaffold_split(actives, decoys, test_frac=0.25)` → group-disjoint **Bemis-Murcko** train/test (no scaffold in both); `random_split` for comparison.
- `audit.py` — `decoy_bias_audit(actives, decoys)` → max Morgan-Tc-to-active distribution, fraction>0.35, property match, verdict `unbiased|mild-bias|biased`. **Use as a gate before any optimizer sees a dataset.**

## How the suite selects the *correct* optimizer (the decisive protocol)
1. For each dataset: `decoy_bias_audit` (gate) → `scaffold_split` → optimize on **train**, evaluate **once on held-out test** (BEDROC α=20 + bootstrap CI).
2. Rank optimizers **per target**; aggregate with a **Friedman + Nemenyi** rank test across targets (Demšar 2006) → critical-difference diagram.
3. **Bias-control diagnostic:** an optimizer that ranks high on **DUD-E** but **collapses on MUV/LIT-PCBA (same/similar target, e.g. ADRB2)** is exploiting decoy bias, not pharmacophore signal. The correct optimizer is **stable across targets AND robust across decoy schemes.**

## Methodology grounding (Crossref-verified)
- **Unbiased benchmarks:** LIT-PCBA 10.1021/acs.jcim.0c00155 · MUV 10.1021/ci8002649 · Xia chemokine MUBD 10.1021/acs.jcim.8b00004 · MUBD-DecoyMaker2 10.1002/minf.201900151.
- **Bias control (justifies the diagnostic):** DUD-E 10.1021/jm300687e · Wallach & Heifets 2018 (benchmarks reward memorization) 10.1021/acs.jcim.7b00403 · Sieg 2019 10.1021/acs.jcim.8b00712 · **Chen 2019 hidden DUD-E bias** 10.1371/journal.pone.0220113.
- **Split + stats:** Bemis-Murcko 10.1021/jm9602928 · scaffold-split generalization (Zhu 2022) 10.1021/acs.jcim.2c01149 · Friedman 10.1080/01621459.1937.10503522 · Demšar 2006 (JMLR, not Crossref-indexed: jmlr.org/papers/v7/demsar06a.html).
- **Decoy generation:** Imrie 2021 (DL decoys) 10.1093/bioinformatics/btab080.

## Recommended additions (from the literature review)
1. **DEKOIS 2.0** (10.1021/ci400115b) — a *medium-bias* waypoint → a 3-point DUD-E→DEKOIS→MUV bias-rank curve per optimizer (richer diagnostic).
2. **Gatica & Cavasotto 2011 GLL/GDD** (10.1021/ci200412p) — a 2nd independent GPCR/CCR2 benchmark → triangulate the CCR2 signal across construction methodologies.
3. **Weighted Friedman** (weight by √n_actives) + per-target bootstrap CIs — LIT-PCBA/MUV active counts are small and heterogeneous.

## Known gaps
- `created/test_load.py` uses a `main()` runner, not pytest functions (logic works, not auto-collected) — convert to `test_*` for CI.
- `litpcba/muv/dude` need `fetch.py` run once (raw data is gitignored / on-demand); `created/`, `ccr2_mubd/`, `ccr2_project/` data is committed and loads immediately.

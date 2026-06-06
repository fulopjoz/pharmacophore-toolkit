# Benchmark Suite Spec — for finding the CORRECT consensus-pharmacophore optimizer

**Goal:** a multi-target, decoy-bias-controlled, scaffold-split benchmark so optimizer selection (greedy / Bayesian / NSGA-II / DOE / combinatorial / learned) is based on **generalization**, not on exploiting one biased CCR2 set. Fixes review findings C2 (overfit/no-held-out) and the unaudited decoy bias.

## Requirements (each dataset must provide)
- **Actives ≥ ~30** (ideally 50+) with SMILES + a target id; **inactives/decoys** with SMILES.
- **A decoy provenance flag**: `experimental-inactive` (gold) | `max-unbiased` | `property-matched` (bias-suspect).
- A **Bemis-Murcko scaffold split** (train/test) — provided by the shared tooling, not per-dataset.
- A uniform on-disk format: `actives.csv` / `decoys.csv` (col `smiles`), + `meta.json` (target, source, DOI, decoy_type, n_act, n_dec).
- A `load.py` returning `(actives_smiles, decoys_smiles, meta)`; **large raw files gitignored**, only a fetcher + a small validated sample committed.

## The suite (tiers)
| Tier | Dataset | Decoy type | Why | Agent |
|------|---------|-----------|-----|-------|
| A (gold, unbiased) | **LIT-PCBA** (Tran-Nguyen 2020) | experimental inactives | the modern unbiased VS standard; multi-target | A |
| A | **MUV** (Rohrer-Baumann 2009) | max-unbiased (spatial stats) | removes analog bias by design; 17 targets | B |
| B (CCR2-specific) | **Xia 2018 chemokine MUBD** (CCR2) + existing project CCR2 set | max-unbiased / current | directly on-target; compare to current set | C |
| C (bias control) | **DUD-E** (Mysinger 2012), 2–3 targets incl. a GPCR | property-matched (biased) | the DIAGNOSTIC: detect optimizers that exploit bias | A or C |
| D (created) | **ChEMBL/PubChem pull** — CCR2 (CHEMBL4015/P41597) + 2 other GPCRs, MUV-style max-Tc-filtered decoys | created max-unbiased | self-contained, reproducible, via database-lookup REST | D |

## Cross-cutting tooling (agent E)
- `benchmark/datasets/split.py` — **Bemis-Murcko scaffold split** (train/test, stratified by label), deterministic seed.
- `benchmark/datasets/audit.py` — **decoy-bias audit**: max Morgan-Tc(radius2,2048) of each decoy to the nearest active; report the distribution + fraction with Tc>0.35 (analog-bias flag); also property-match check (MW/LogP/HBD/HBA).
- Both with unit tests (TDD).

## The decisive diagnostic (how the suite finds the "correct" optimizer)
Run every optimizer on every dataset's **train scaffold split**, evaluate **once on the held-out test split** (BEDROC α=20 + bootstrap CI). The correct optimizer is the one that:
1. **ranks highest on average across the unbiased targets** (Tier A+B), and
2. **does NOT show a large drop from DUD-E→MUV/LIT-PCBA on the same/similar target** (Tier C diagnostic — a big drop = bias-exploitation), and
3. has **stable rank** across targets (low variance), not one big win.
Report the per-optimizer × per-dataset BEDROC matrix + a Friedman/Nemenyi rank test across targets (Demšar 2006 style).

## Literature anchors (verified by agent F)
LIT-PCBA (Tran-Nguyen 2020), MUV (Rohrer-Baumann 2009), DUD-E (Mysinger 2012), DEKOIS 2.0 (Bauer 2013), Xia chemokine MUBD (2018), Wallach-Heifets 2018 (benchmarks reward memorization), Demšar 2006 (comparing algorithms across datasets).

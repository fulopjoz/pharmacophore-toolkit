# PRISM beyond CCR2 — cold-start, low-diversity & a point-cloud consensus P4

> **STATUS: PARKED FINDINGS (not active work).** §4 (point-cloud registration + diffusion
> pharmacophore generation) and the broader cold-start research are documented here for the
> future but are explicitly **deferred** — the active plan is to *finish* the validated CCR2
> production path (PRISM → DrugEx reward), see `docs/research/PRISM_DRUGEX_PORT_PLAN.md`.

**Date:** 2026-06-07 · scite literature sweep + scientific-critical-thinking + superpowers planning.
**Why:** PRISM is validated on CCR2 (learned discrimination weighting is causal; ESP helps on unbiased decoys; 10-seed `BAKEOFF_MULTISEED.md`). The open challenge is *new* targets whose data breaks PRISM's two requirements — **supervision (actives + decoys)** and **chemotype diversity (multiple templates)**:
- **IDP target, 1 known binder** (a green-tea natural product — almost certainly EGCG): no supervision, no diversity.
- **Target with 12 measured actives, all very similar**: supervision yes, diversity no → weights overfit one scaffold; no scaffold-hopping.

## 0. Diagnosis (what actually breaks)
| regime | supervision (act+dec)? | chemotype diversity? | PRISM works? | the fix targets… |
|---|---|---|---|---|
| CCR2 (74 act, ~8 clusters) | yes | yes | **yes (validated)** | — |
| 12 similar actives | yes | **no** | partial (scaffold-locked) | add diversity / abstract the features |
| 1 active (EGCG/IDP) | **no** | **no** | no | manufacture both (augment + transfer) |

## 1. Scenario A — one active (EGCG / IDP), structure-free
Graded, weakest → strongest:
1. **Unsupervised single-molecule fuzzy pharmacophore** + feature-type priors (Renner 2005, fuzzy pharmacophore, 10.1002/cbic.200400376). A floor.
2. **Active-set augmentation BEFORE modeling** (highest value — systematise the Papyrus mining you did for CCR2):
   - similarity-mine ChEMBL/Papyrus for analogs (database-lookup skill);
   - **tautomer/protomer enumeration** — mandatory for polyphenols (EGCG: 8 phenolic OH, redox/ionizable) → one SMILES ≠ one pharmacophore;
   - **bioisostere / scaffold-hop generation** for *chemotype* diversity (Ziv 2025 flow-matching bioisosteres, 10.1101/2025.10.20.683377; or DrugEx itself).
3. **Few-shot / transfer weighting** — a GLOBAL feature-type-importance prior across many targets, or transfer from a related target, lightly adapted (Vella & Ebejer 2022, 10.1021/acs.jcim.2c00779; MHNfs in-context, Schimunek 2025, 10.1021/acs.jcim.4c02373).
4. **Active-learning bootstrap** — 1 active → DrugEx generate → dock/assay top picks → grow active set → then learn PRISM weights.
IDP precedent for ligand-/pharmacophore-based design (no structure): c-Myc (Chen 2016, 10.1038/srep22298), NUPR1 (Neira 2017, 10.1038/srep39732); EGCG is a well-known promiscuous IDP/amyloid binder.

## 2. Scenario B — 12 near-identical actives
Fight scaffold overfit:
1. Build the **tight consensus** (reliable shared features) for *same-chemotype* expansion.
2. For novelty: **feature-type-level weights** (pool across templates, not per-template-position), fuzzy/2D pharmacophore, strong L2; augment chemotype diversity (DB mining + generative scaffold-hops).
3. **External validation only** — one cluster ⇒ no internal leave-one-scaffold-out; validate on a related target or the CCR2 simulation (§3).

## 3. The de-risking experiment — SIMULATE both on CCR2 (we already validated it)
Because we have 74 CCR2 actives + a validated benchmark, we can quantify the cold-start *before* the real IDP:
- **1-active sim:** pick 1 active → build model → score the full held-out benchmark; repeat over many single actives → distribution of enrichment. Then measure how each §1 strategy (augmentation, tautomers, bioisosteres, transfer prior) RECOVERS enrichment.
- **12-similar sim:** pick 12 actives from ONE Butina cluster → build → score held-out (incl. *other* chemotypes) → quantify scaffold-hopping failure, and whether abstraction/augmentation recovers it.
This is a controlled, internal, rigorous read on the cold-start protocol with zero new wet data.

## 4. The novel consensus P4 via cross-domain translation (your point-cloud idea — it's sound)
A pharmacophore IS a **typed 3D point cloud** (points labelled donor/acceptor/…). Translate from 3D vision:
- **Consensus building = correspondence-free, SO(3)-equivariant point-cloud registration + co-registration** (Zhu 2021, 10.48550/arxiv.2107.10296; DeepBBS) → align/merge multiple actives' pharmacophores into one consensus *without* needing atom correspondence — faster and more robust than ROCS-pairwise + Butina. Optimal-transport / Gromov-Wasserstein handles different feature counts (we prototyped fused-GW already).
- **Few-shot completion / augmentation = diffusion-based pharmacophore GENERATION** (PharmacoForge, Flynn 2025, 10.3389/fbinf.2025.1628800; ShEPhERD; OctFusion octree diffusion, 10.1111/cgf.70198) → generate an optimized consensus + plausible variants from few actives. This *is* the "novel, fast, optimized P4 model" — and it directly attacks both hard regimes (complete the pharmacophore from 1; generate diverse variants from 12).

## 5. Phased plan (logical next steps)
- **Phase 1 — consolidate (low risk, now).** Expand benchmark to **CCR5/CXCR4** (have created actives) + **MUV/LIT-PCBA** with per-target references (fixes the global-ref trap), power the Friedman test (N≥6, force the 2D canary down). Port PRISM's frozen weighted-color into the DrugEx `rocs_*` reward (gated PRCD). *Skills:* database-lookup, statistical-analysis/statsmodels, serena, tdd, code-review.
- **Phase 2 — cold-start protocol + CCR2 simulation (§3).** Build the augmentation pipeline (DB mining, tautomer/protomer, bioisostere generation, decoy generation) + a global feature-type weight prior; validate via the subsampling sim. *Skills:* database-lookup, rdkit/datamol, scientific-critical-thinking, statistical-analysis.
- **Phase 3 — research bet: point-cloud consensus + generative pharmacophore.** Prototype correspondence-free equivariant registration for consensus, and a diffusion pharmacophore generator (PharmacoForge/ShEPhERD-style) for few-shot; benchmark vs ROCS-consensus on the suite. *Skills:* torch-geometric, deepchem, scientific-brainstorming, tdd.
- **Phase 4 — apply to the real targets.** EGCG/IDP (1-active protocol) + the 12-similar target; prospective generate→assay loop. *Skills:* the DrugEx RL pipeline, dock+PLIP where a structure exists.

**Dependencies:** Phase 1 ∥ (Phase 2 needs the CCR2 sim harness); Phase 3 builds on Phase 2's augmentation; Phase 4 needs Phases 2–3. Start: Phase 1 + the Phase-2 CCR2 simulation (both cheap, both de-risk the IDP work).

## References (scite, DOIs)
Vella & Ebejer 2022 few-shot low-data 10.1021/acs.jcim.2c00779 · Schimunek 2025 MHNfs 10.1021/acs.jcim.4c02373 · Qian 2024 meta-learning GNN 10.1021/acsomega.4c02147 · Renner 2005 fuzzy pharmacophore 10.1002/cbic.200400376 · Ziv 2025 flow-matching bioisosteres 10.1101/2025.10.20.683377 · Zhu 2021 SO(3)-equivariant correspondence-free registration 10.48550/arxiv.2107.10296 · Flynn 2025 PharmacoForge diffusion pharmacophore 10.3389/fbinf.2025.1628800 · Xiong 2025 OctFusion 10.1111/cgf.70198 · Chen 2016 c-Myc IDP 10.1038/srep22298 · Neira 2017 NUPR1 IDP 10.1038/srep39732 · Cleves & Jain 2020 multi-active query 10.1021/acs.jcim.0c00115.

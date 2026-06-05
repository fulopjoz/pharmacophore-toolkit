# S3 Discrimination-Weighting — Tier-1 Standalone Test (NO DrugEx)

**Date:** 2026-06-05 · **Branch:** `s3-enrichment-test` · **Script:** `experiments/s3_enrichment_test.py`
**Question:** does weighting the 6 P4 pharmacophore color-feature types by actives-vs-decoys **discrimination** beat **equal** weighting on retrospective enrichment — *before* any DrugEx/RL investment?

## Method
Fast 2D proof-of-concept (RDKit BaseFeatures P4 counts, no conformers) on the project's **74 CCR2 actives + 498 decoys + 5 references**. 5-fold stratified CV, **weights fit on train fold, scored on held-out fold** (no leakage). Metrics via RDKit `Scoring` (`CalcBEDROC` α=20, `CalcEnrichment`, AUC). Bootstrap 95% CI (1000×) on the metric **difference** vs the equal-weight baseline. Three arms: **A** equal-weight similarity to 5 refs (current-production proxy); **B** S3-weighted similarity to 5 refs; **C** S3-weighted absolute P4 (no reference).

## Result: **GO**

**Learned S3 weights (held-out CV) — independently reproduce the literature/research finding:**
| donor | hydrophobe | acceptor | rings | anion | cation |
|---|---|---|---|---|---|
| +0.74 | +0.68 | +0.64 | +0.63 | +0.29 | **−0.93 (down — decoy-like)** |

**Held-out enrichment:**
| arm | AUC | EF1% | EF5% | BEDROC |
|---|---|---|---|---|
| A equal-weight | 0.640 | 1.29 | 1.60 | 0.175 |
| **B S3-weighted (refs)** | **0.873** | 3.87 | 5.60 | **0.604** |
| C S3-weighted (absolute) | 0.835 | 1.29 | 3.20 | 0.460 |

**Bootstrap 95% CI on B − A:** ΔAUC **+0.233 [+0.157, +0.306]** (excludes 0); ΔBEDROC **+0.400 [+0.239, +0.571]** (excludes 0); ΔEF1% +2.6 [−1.5, +6.6] (ns — too noisy at n=74).

## Interpretation
1. **S3 works.** Discrimination-weighting beats equal-weighting on AUC (0.64→0.87) and BEDROC (0.18→0.60), both with bootstrap CIs excluding 0 on a held-out split. **Gate cleared → proceed to Tier-2.**
2. **My earlier "donor-poor refs → reweighting is inert" hypothesis was WRONG.** B (reweight similarity-to-refs) actually *beat* C (absolute), BEDROC 0.604 > 0.460 — even with donor-poor references, the other discriminating types (hydrophobe, rings, acceptor) carry enough signal, and comparing-to-references adds information over raw counts. So reweighting helps *with the existing references*; reference curation is an *additional* lever, not a prerequisite.
3. **Weights independently reproduce the research** (cation anti-correlated; donor/hydrophobe/acceptor/rings active-specific) — cross-validates the feature-weighting study from a completely separate code path.

## Caveats (honest scope)
- **2D P4 proxy, not 3D ROCS color** — Tier-2 must confirm in the OpenEye ROCS production path (custom `.cff`).
- **EF1% is too noisy** (≈0.7 actives in the top 1% of 572; CI includes 0) — AUC/BEDROC carry the conclusion.
- **Decoys are property-matched (DUD-E-style), not assay-confirmed** → relative enrichment only.
- **43/74 actives overlap the DrugEx finetune set** — fine for a retrospective *scoring* test, but for the production claim, re-derive weights on a **training-disjoint** split and report held-out AUC.

## Next (only now justified)
- **Tier-2:** build the discrimination-weighted OpenEye `.cff` (per-type `weight=` from these coefficients) and re-run enrichment in real 3D ROCS color (`rocs -chemff custom.cff`) vs the default equal-weight `.cff`.
- **Then** the DrugEx integration (the frozen-weights-as-reward architecture in `ccr2_gen/docs/color_validation/figures/optimization_architecture.png`).

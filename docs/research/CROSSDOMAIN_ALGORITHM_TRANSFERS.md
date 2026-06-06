# Cross-Domain Algorithm Transfers for Shape+Color Molecular Alignment

**Date:** 2026-06-06 · **Method:** 7-agent parallel research (paperclip + WebSearch), all DOIs Crossref-verified.
**Framing:** a molecule = a union of typed Gaussian "bubbles"; aligning/scoring it is formally the same problem solved in computer vision (point-cloud registration), applied math (optimal transport), astronomy (star-pattern matching), soft-matter physics (foam/colloid energy minimization), and medical image registration (multi-modal fusion).
**Already in the fork (excluded as "new"):** optimal transport/Wasserstein, colored ICP, Hungarian assignment, NSGA-II, Bayesian + multi-fidelity BO, simulated annealing, DOE, Optuna.

## The project's core defect this targets (P1)
Verified in code: **both open backends optimize the pose on shape, then read color off it** — RDKit `AlignMol` blends only via `opt_param` (≈inert at production settings); CDPKit's BFGS has **no color gradient** (`calcColorOverlapGradient` doesn't exist), so color is only approached by feature-center seeding. **Color is a passenger.** The durable fix is to change *where the local refiner starts* — hand it a pose where the typed color features already coincide.

## Top-5 transferable tricks (ranked by impact × low cost)

| # | Trick (source field) | What it does | Solves | Backend | DOI |
|---|---|---|---|---|---|
| 1 | **TEASER++** (robotics/CV: certifiable robust registration) | color-driven *initial pose* (same-type-only feature matches → globally-certified SE(3) in ms) seeding the native RDKit/CDPKit refiner; robust to >99% outliers; returns an optimality **certificate** | **P1**, P3 | RDKit + CDPKit (pre-align the probe before `AlignMol`/`align`) | 10.1109/TRO.2020.3033695 |
| 2 | **CLIPPER** (robotics: max-consistency clique) | pose-free typed-feature similarity AND a clean inlier seed for TEASER++; enforces joint geometric consistency (Hungarian doesn't) | **P4**, P1, P3 | backend-agnostic prefilter | 10.1109/ICRA48506.2021.9561069 |
| 3 | **Information-weighted fusion** (robotics: inverse-covariance/Kalman) | derive per-channel weights `w_k ∝ 1/var_k` from actives/decoys; covariance-intersection deflates correlated donor/acceptor; **derives the S3 weighting from first principles** (cation auto-collapses) — no tuning loop | **P2** | all 3 (RDKit/CDPKit weighted combine; OpenEye `.cff`) | 10.1109/TAC.2013.2277621 |
| 4 | **Rarity / TF-IDF log-odds weighting** (astronomy: catalog-density verification) | weight each color type by inverse background frequency over the 1M library; analytic, **lowest cost**; complements #3 | **P2**, P4 | OpenEye `.cff` (cleanest); RDKit re-aggregate | 10.1088/0004-6256/139/5/1782 |
| 5 | **Unbalanced/partial Fused Gromov-Wasserstein** (applied math) | Gromov term compares *intra*-molecule distances → **pose-free**; fused term adds typed-feature cost; partial relaxation lets a small fragment match a sub-pharmacophore of a big reference (no forced 1:1) | **P4**, P2 | backend-agnostic (POT `ot.gromov`/`ot.partial`) | 10.1090/mcom/3303, 10.3390/a13090212 |

## Highest-leverage single transfer: **TEASER++ color-aware pose seed** (P1)
It attacks the defect at the exact leverage point — the *initialization* — without replacing the validated native overlay. It is **color-driven by construction** (matches are generated only donor→donor, ring→ring…), **fast** (ms; survives RL throughput), **robust** to the dominant nuisance (query≠reference → most feature pairs are outliers, which break ICP/Hungarian), and returns a **certificate** the RL loop can use to down-weight ambiguous alignments. Pair with **CLIPPER** (supplies the consistent inlier set). Together: a complete, RL-affordable, *color-first* pose pipeline that turns color from passenger into driver of the starting pose.

## The "use all three backends and compare" experiment (designed)
**Information-weighted color-channel fusion**, held-out CCR2 actives/decoys, 5 references, k-fold (weights on calibration fold, evaluate on disjoint fold):
- **RDKit:** `rocs_rdkit.py` return per-feature-type color overlaps; replace `shape+color` with `Σ_k w_k·channel_k`.
- **CDPKit:** read per-color-feature overlaps off the `AlignmentResult`; same weighted combine; start generator unchanged.
- **OpenEye:** custom `.cff` with the 6 derived per-type weights via the wired `-chemff` + `-optchem true` (native per-type path).
- **Metric (matches Tier-1):** BEDROC(α=20) primary; ROC-AUC + EF@1% secondary; bootstrap CIs + DeLong; per backend × reference × {flat-sum, S3, info-weighted}.
- **Gate:** adopt info-weighting if BEDROC improves ≥ +0.05 over S3 in ≥ 2 of 3 backends. (Record per-backend weights → documents backend non-equivalence.)

## For the eventual DrugEx-RL integration
- **Pose-free (P4, fast prefilter):** CLIPPER, Fused-GW (unbalanced/sliced), astrometry quad geometric-hash, Groth triangle voting — rank *without* aligning → run the expensive aligner only on top-k survivors in `DrugExEnvironment.getScores`.
- **Differentiable (P5, GPU in-the-loop reward):** Sinkhorn / debiased Sinkhorn-divergence (10.1109… arXiv:1306.0895, 1706.00292), entropic Gromov-Wasserstein (arXiv:1810.08278), Coherent Point Drift GMM-L2 (10.1109/TPAMI.2010.46), B-spline mutual information (10.1109/83.887976) — use during RL where gradients/throughput matter; keep exact aligned TanimotoCombo for final ranking.
- **Global-search hardening (P1 fallback):** if TEASER++ seed still lands shape-greedy, budget-bounded CMA-ES / Differential Evolution polish warm-started from the seed (never worse than current); or physics minimizers FIRE (10.1103/PhysRevLett.97.170201) / basin-hopping (10.1021/jp970984n) — the literal "bubbles" optimizers.

## Next steps (gated)
1. Standalone color-aware seed module: typed features → CLIPPER inliers → TEASER++ SE(3) → apply to probe before native refine; validate color overlap rises without shape collapse on the 5 references.
2. Offline gain: held-out BEDROC/AUC for {PCA-start} vs {TEASER++ seed→refine} in RDKit & CDPKit (PBS, ncpus=46, /scratch).
3. The 3-backend fusion experiment above.
4. Pose-free prefilter at 1M scale (geometric hash / FGW) wired as a gate before the aligner.
5. Differentiable Sinkhorn/FGW companion scorer for RL (verify Spearman vs exact score first).

## References (all Crossref/arXiv-verified 2026-06-06)
TEASER++ 10.1109/TRO.2020.3033695 · CLIPPER 10.1109/ICRA48506.2021.9561069 · Info-Weighted Consensus Filter 10.1109/TAC.2013.2277621 · astrometry.net 10.1088/0004-6256/139/5/1782 · Groth triangle-match 10.1086/114099 · Fused-GW 10.3390/a13090212 · Unbalanced OT 10.1090/mcom/3303 · Coherent Point Drift 10.1109/TPAMI.2010.46 · FIRE 10.1103/PhysRevLett.97.170201 · basin-hopping 10.1021/jp970984n · Sinkhorn arXiv:1306.0895 · Sinkhorn divergences arXiv:1706.00292 · entropic GW arXiv:1810.08278 · CMA-ES arXiv:1604.07269 · DE point-cloud 10.1371/journal.pone.0209227 · B-spline MI 10.1109/83.887976.

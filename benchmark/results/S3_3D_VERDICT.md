# Verdict — discrimination-weighted 3D color (`s3_3d`) vs unweighted 3D

**Run:** PBS job 10405 (7h45m), seed 42, test_frac 0.25, scaffold-split held-out, bootstrap 2000.
**Question (decision gate):** does per-(feature-type × template) discrimination weighting rescue 3D ROCS color on *unbiased* decoys, where unweighted 3D shape+color collapses to ≈random?

## Answer: YES — decisively, with non-overlapping CIs.

BEDROC(α=20) [95% CI], primary metric:

| dataset (decoys) | unweighted 3D (rdshape / shape_combo) | **s3_3d (weighted 3D)** | 2D s3_weighted | 2D Morgan canary |
|---|---|---|---|---|
| ccr2_mubd (**max-unbiased**) | 0.006 [0,0.019] / 0.002 [0,0.006] | **0.293 [0.114,0.492]** | 0.175 [0.032,0.351] | 0.705 [0.512,0.868] |
| created_CCR2 (created-unbiased) | 0.378 [0.274,0.481] / 0.365 | **0.952 [0.924,0.970]** | 0.417 [0.325,0.509] | 0.987 [0.979,0.993] |
| ccr2_project (property-matched) | 0.869 [0.672,0.963] / 0.729 | **0.972 [0.902,1.000]** | 0.781 [0.532,0.961] | 0.971 [0.895,1.000] |

Cross-dataset avg rank (lower=better): differential_mmfp 1.33 · **s3_3d 1.67** · learned_scorer 4.00 · s3_weighted 4.33 · rdshape_ensemble 4.33 · pharm2d 6.33 · shape_combo 6.67 · equal_weight 7.33. Friedman χ²=17.2, p=0.016 (Nemenyi CD=6.06, underpowered at N=3 — use the per-dataset CIs as the evidence).

## What is established (CI-backed)

1. **Discrimination-weighting rescues 3D color.** On MUBD, `s3_3d` 0.293 [0.114,0.492] vs unweighted-3D ≤0.019 upper bound → **CIs do not overlap** (~50×). On created_CCR2 the same with tight CIs (0.952 vs 0.378). The lever the prior work was missing was **per-feature discrimination weighting**, not color *magnitude* (Tversky/hydrophobe).
2. **3D adds value over 2D when weighted** — significant on 2/3 datasets: created_CCR2 ΔBEDROC vs s3_weighted **+0.537 [+0.443,+0.624]**; ccr2_project **+0.184 [+0.017,+0.430]**. On MUBD directional only (+0.114 [−0.129,+0.353]; just 17 test actives → wide CI).
3. `s3_3d` is the **#1 pharmacophore method** and the **only 3D method that does not collapse** on MUBD.

## Honest caveats (do not over-claim)

- **The 2D Morgan canary (`differential_mmfp`) still tops `s3_3d`** (avg rank 1.33; MUBD 0.705 vs 0.293, Δ +0.538 [+0.268,+0.768]). This is the distributional-memorization effect (Wallach & Heifets 2018): even MUBD leaves a learnable 2D boundary. It is **not portable** to a 3D pose-based generative reward and is kept only as the bias canary — but its dominance means **residual decoy bias remains**.
- The **MUBD gate verdict was `mild-bias`** (median max-Tc 0.211 > 0.20), so this set is not perfectly unbiased — which partly explains why the canary still reaches 0.705. A stricter set is needed to push the canary down.
- N=3 datasets → Nemenyi is underpowered; the MUBD 3D-vs-2D gap is not individually significant.

## Conclusion & next steps

- **GO** to port the per-(type × template) weights `w[t,f]` into the DrugEx production `rocs_*` color path (Phase 2/3 of `COLOR_SCORING_AUDIT_AND_PLAN.md`): a discrimination-weighted ColorTanimoto / gated PRCD color objective, with weights fit on the 74 CCR2 actives vs decoys.
- **Power the test:** add MUV + LIT-PCBA (truly unbiased) and CCR5/CXCR4 (cross-target) so the Friedman test has N≥6 and the canary is forced down — confirming `s3_3d` robustness where even 2D memorization should fail.
- **Structural ground-truth:** color is pose-blind; validate winning candidates with dock + PLIP essential-contact recapitulation (the project's existing endpoint), not color alone.

---

## Multi-seed confirmation — 10 scaffold-split seeds (job array 10411, 2026-06-07)

`BAKEOFF_MULTISEED.md` — BEDROC mean ± across-seed 95% CI + PAIRED deltas (same split per
seed). This adds the split-choice variance a single-split bootstrap cannot see. (PRISM is
the renamed s3_3d.)

**1. The learned weighting is causal — `prism − prism_fixed` is CI-separated on ALL three datasets:**
- ccr2_project (unbiased gate): +0.040 [+0.018,+0.063], >0 in 90% of seeds
- ccr2_mubd (max-unbiased): +0.096 [+0.026,+0.166], 80%
- created_CCR2: +0.027 [+0.026,+0.029], 100%
`prism_fixed` shares prism's *byte-identical* cached features, so this isolates ONLY
learned-vs-fixed weighting. Confirmed: the discrimination weighting (not the 3D features
alone) drives PRISM's win.

**2. Electrostatics (`prism_esp`) adds real signal where it's hard, and is null where it's easy — the opposite of overfitting:**
- ccr2_project (easy/property-matched): +0.000 [-0.001,+0.002] — NULL (CI includes 0)
- ccr2_mubd (max-unbiased): +0.054 [+0.020,+0.087], 80% — CI-separated POSITIVE
- created_CCR2: +0.006 [+0.005,+0.006], 100% — small CI-separated positive
ESP helps on the unbiased/hard sets, not the biased/easy one → genuine electrostatic
discrimination, not overfitting capacity. Keep prism_esp.

**3. Top BEDROC per dataset (10-seed mean):** ccr2_project prism/prism_esp 0.981; ccr2_mubd
prism_esp 0.310 / prism 0.257 (unweighted 3D shape_combo 0.007, rdshape 0.020 ≈ random);
created_CCR2 prism_esp 0.954 / prism 0.948. Unweighted 3D collapses on MUBD; PRISM holds.

**Honest boundary:** the 2D Morgan canary still tops MUBD (0.817) = residual decoy bias;
and 10 seeds tighten per-dataset CIs but do NOT power the Friedman test (correlated
resamples). Powering it needs more TARGETS (CCR5/CXCR4, MUV/LIT-PCBA) — the next step.

**Production recommendation:** port the discrimination-weighted color (PRISM weights) into
the DrugEx `rocs_*` reward as the gated PRCD color objective; optionally include the ESP
channel. Validate winners with dock+PLIP (color is pose-blind).

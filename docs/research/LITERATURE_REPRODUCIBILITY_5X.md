# Literature-Search Reproducibility — 5× Dual-Source Study

**Date:** 2026-06-06 · **Method:** ran the same 6-topic shape-alignment literature search **5 independent times in parallel**, each querying **scite + paperclip**, then evaluated reproducibility, cross-source coverage, and evidence quality (skills: statistical-analysis, scholar-evaluation, scientific-critical-thinking, citation-management).
**Why:** convert "here are some papers" into "here is the *reproducible* evidence vs run-to-run noise," and check whether the earlier cross-domain algorithm findings survive repetition.

## Reproducibility: HIGH
- **Mean pairwise Jaccard across topics = 0.823** (10 run-pairs per topic).
- **Stable core (≥4/5 runs) = 70 DOIs**; variable (2–3/5) = 10; sporadic (1/5) = 16 → the 70-DOI core is the citable evidence base; the 26-DOI tail is run-dependent noise.
- **scite is deterministic** — byte-identical ordered sub-lists across all 5 runs for T1–T4. **Essentially all variance comes from paperclip's biomedical/preprint layer** (and a single high-recall outlier, Run 3).
- Most reproducible: T4 ROCS/scaffold-hopping (J=0.900), T3 pharmacophore/shape-VS (0.883). Least: T6 pose-free (0.729 — "screening" polysemy noise), T2 Gromov-Wasserstein (0.761 — fast-moving arXiv tail).

## Cross-source (scite vs paperclip): complementary, converge at the concept level
- **Zero shared DOIs in every topic** — the two engines return *disjoint identifier sets*. scite = broad CS/physics/cheminformatics + journal record (with citation tallies); paperclip = biomedical/PMC + preprint servers (bioRxiv/chemRxiv/TPAMI/Nat Commun).
- **But they converge at the topic level**: both independently populate the *same* sub-literatures (Gromov-Wasserstein, shape similarity) from different document pools. **Two independent retrieval systems agreeing on the central sub-literatures is the strongest reproducibility signal.** Stable-core split: scite anchored 41/70, paperclip 29/70.

## Earlier cross-domain findings: CORROBORATED
1. **Certifiable point-cloud registration (TEASER++)** — TEASER (5/5), Go-ICP (5/5, TPAMI), Guaranteed Outlier Removal (5/5, TPAMI), max-feasible-subsystem global registration (5/5). Anchored independently by *both* sources.
2. **Gromov-Wasserstein for shape/molecule matching** — 10/13 T2 stable-core are GW (Quantized/linear/Outlier-Robust/Fused/Sliced-GW). **NEW & directly relevant: OTMol — "Robust Molecular Structure Comparison via Optimal Transport" (arXiv 2509.01550, 5/5).** *Caveat:* many GW arXiv entries have mentions=0 → theme reproducible, individual papers still emerging.
3. **Shape+color + interaction-fingerprint fusion** — the best-evidenced cluster: SHAFTS (145 mentions), atom/feature-pair+volume VS (193), **ROCS scaffold-hop (608 mentions, support +5)**, interaction-fingerprint docking (555, +4), and **Differential Multimolecule Fingerprint (Hutter, the S3 embodiment, 5/5)**. High citation support → established.

## New high-value papers the reproducible search surfaced
- **OTMol** (molecular structure comparison via OT, arXiv 2509.01550) — pose-free molecular matching, directly usable.
- **ShEPhERD** (diffusing shape+electrostatics+pharmacophores, arXiv 2411.04130) — generative shape+pharmacophore.
- **gWEGA** (GPU-accelerated WEGA shape comparison, 10.1002/jcc.23603) — the speed/GPU lever (P3/P5).
- **G3PS** greedy 3-point pharmacophore alignment (10.3390/molecules26237201); **TurboHopp** consistency-model scaffold hopping (arXiv 2410.20660).
- **"Do Molecular Fingerprints Identify Diverse Active Drugs? (No)"** (10.3390/ph17080992) — a caution for the fingerprint direction.

## Evidence-quality flags (scientific-critical-thinking)
- **De-duplicate by content before any count:** the same decoy-selection paper appears under **3 DOIs** (chemRxiv-v2 5/5, chemRxiv 3/5, published 10.1186/s13321-025-01107-z 3/5); the "fingerprints (No)" paper is double-counted (preprint 10.1101/2022.09.20.508800 + published 10.3390/ph17080992). Both inflate T5.
- **Cross-discipline noise (correctly relegated to sporadic):** T6 "screening" polysemy pulled a 1992 femtosecond-optics paper and a generic MD textbook chapter.
- **T1 conflates two literatures:** robotics/CS registration (TEASER/Go-ICP) vs biomedical microscopy registration (CLEM-Reg) — both reproducible but answering different questions.
- **Recommend a scite `editorialNotices` retraction/EoC check** on the bioRxiv/chemRxiv preprints before citing (notices weren't in this dataset).
- No stable-core paper carried a net-contradicting citation profile (all scite `support` ≥ 0).

## Bottom line
The earlier cross-domain recommendations (TEASER++ registration, Gromov-Wasserstein/OT matching, info-fusion of shape+color+IFP) are **reproducibly supported** by an independent 5×, two-source search — not single-run luck. The strongest *new* leads are **OTMol** (molecular OT, pose-free) and **gWEGA** (GPU shape) for the speed/differentiable track, and **Differential Multimolecule Fingerprint** independently re-confirms the S3 discrimination-weighting basis.

*All DOIs above were returned by the live tools; verify via Crossref + a scite editorialNotices check before formal citation. scite token valid until 2026-06-06 22:14 UTC.*

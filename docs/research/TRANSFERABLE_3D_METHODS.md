# Transferable 3D shape/pharmacophore methods — landscape + critical assessment

**Date:** 2026-06-07 · **Method:** scite premium literature sweep (6 queries) + scientific-critical-thinking.
**Purpose:** the manuscript's story is the **RL pipeline**, but the *reward* (3D ROCS shape+color, `rdShapeAlign`) is what reviewers attack. This catalogs published 3D methods/ideas we can **transfer to upgrade that reward**, and critically assesses which to add to the bake-off.

## 1. Where our method sits

Our reward — `rdShapeAlign` Gaussian shape overlay + 6-type pharmacophore "color" — **is the ROCS / SHAFTS / PheSA / ROSHAMBO family** (Gaussian volume overlap + pharmacophore features). The one thing that family mostly does NOT do is **weight features by actives-vs-decoys discrimination** — which is exactly the `s3_3d` upgrade (validated: job 10405 — rescues 3D color on unbiased MUBD, see `benchmark/results/S3_3D_VERDICT.md`). The closest prior art to `s3_3d` is **LigMate** (learned multi-feature weighting) and **PharmRF** (ML pharmacophore scoring).

## 2. Transferable methods (scite-grounded), grouped by what they upgrade

| Upgrade | Method (paper, DOI) | What to lift | Availability |
|---|---|---|---|
| Engine swap (modern OSS ROCS) | ROSHAMBO (Atwi 2024, 10.1021/acs.jcim.4c01225) | GPU shape+color engine | external, GPU, NOT installed |
| Fixed weighted hybrid | SHAFTS (Liu 2011, 10.1021/ci200060s) | shape+pharmacophore-triplet *fixed* combination | not a Python lib → reimplement concept |
| **Learned** feature weighting | LigMate (Grimm 2020, 10.1021/acs.jcim.9b01210); PharmRF (Kumar 2022, 10.1002/jcc.26840); Rehioui 2022 (10.1002/minf.202200210) | learn per-feature weights — **= what s3_3d does** | concept ✓ (cite) |
| New channel: electrostatics | ShaEP (Vainio 2009, 10.1021/ci800315d); Markt 2008 (10.1021/jm800128k); Berenger 2014 (10.1186/1758-2946-6-23) | electrostatic-potential overlap as extra features | internal via RDKit Gasteiger/MMFF charges ✓ |
| Consensus / multi-reference | Cleves & Jain 2020 (10.1021/acs.jcim.0c00115); FLAP (Baroni 2007, 10.1021/ci600253e); Madzhidov 2020 (10.3390/molecules25020385) | multi-active query | ✓ already in s3_3d |
| Open-source alignment tools | LIGSIFT (Roy 2014, 10.1093/bioinformatics/btu692); LS-align (Hu 2018, 10.1093/bioinformatics/bty081); PheSA (Wahl 2024, 10.1021/acs.jcim.4c00516) | alternative aligners | external |
| Benchmark rigor / caveats | Venkatraman 2010 ("limitations of 3D methods", 10.1021/ci100263p); Giganti 2010 (10.1021/ci900507g); Giangreco 2021 field benchmark (10.1021/acs.jcim.1c00866) | honest framing + extra benchmark set | ✓ cite |
| Landscape review | Kumar 2018 (10.3389/fchem.2018.00315) | position the method section | ✓ cite |

## 3. Critical assessment of the three candidate bake-off additions

Each addition is an **experiment**; judged on internal validity / confounds / feasibility.

### (A) SHAFTS-style fixed-weight hybrid → `s3_3d_fixed`  ·  **DO (highest priority)**
- **Tests:** is `s3_3d`'s win due to the *learned discrimination weighting* or merely the 3D per-type features? Holds the feature substrate identical, varies only fixed↔learned. **The key internal-validity control for the paper's central claim.**
- **Confound to avoid:** identical `(K×6)` features; replace the logistic with a *fixed* aggregation (equal mean; optional fixed shape:color ratio). Anything else confounds it.
- **Construct honesty:** this is a *fixed-weight ablation of our features* ("SHAFTS-style"), NOT literal SHAFTS (which is unavailable). Label as such.
- **Feasibility:** internal, trivial (reuse `color3d`). **Risk: low.**

### (B) Electrostatics channel → `s3_3d_esp`  ·  **DO (with controls)**
- **Tests:** does electrostatic complementarity add discrimination beyond shape+color?
- **Confounds / threats:** (i) Gasteiger charges are crude; (ii) ESP read off a *shape+color-optimized* pose is a weak proxy for true field alignment (ShaEP/EON do real field optimization) — label honestly as an *ESP-similarity feature*, not field-based screening; (iii) extra columns ⇒ overfitting risk on ~dozens of actives.
- **Controls (mandatory):** same templates/poses as `s3_3d` (isolate the ESP contribution); L2 + held-out **MUBD** gate; ESP-similarity defined as same-sign Carbó/Hodgkin-style overlap `Σ qᵢqⱼ exp(−α dᵢⱼ²)` (ligand–ligand similarity rewards matching charge patterns). Trust only if it lifts MUBD with non-overlapping CIs.
- **Feasibility:** internal (RDKit charges). **Risk: medium.**

### (C) ROSHAMBO engine swap  ·  **DEFER / optional, run last**
- **Tests:** is a different OSS engine better than `rdShapeAlign`? But it changes conformers + optimizer + color *simultaneously* → low internal validity (can't attribute a difference).
- **Predicted result:** unweighted, same philosophy as `shape_combo`/`rdshape_ensemble` → by our own finding it should **collapse on MUBD**. The informative variant ("ROSHAMBO features + `s3` weighting") is a much larger build.
- **Feasibility:** external, GPU, **not installed**; installing on the shared NFS conda env risks the documented Bus-error race → needs an **isolated venv** (worktree/PBS), adding setup risk.
- **Recommendation:** lowest pass-1 value, highest risk; do it only after (A)/(B), in isolation, flagged as an engine comparison (not an ablation).

## 4. Decision

- **Add now:** `s3_3d_fixed` (weighting ablation) and `s3_3d_esp` (electrostatics channel) — both internal, both directly interrogate `s3_3d`.
- **Defer:** ROSHAMBO to a later round (isolated install), unless explicitly prioritized.
- **Guardrails (user-mandated):** TDD per scorer → **front-node smoke test before any PBS run** → **two independent code reviews** (parallel agents) → then PBS.
- **Power the test in parallel:** the honest next rigor step remains adding MUV/LIT-PCBA + CCR5/CXCR4 (N≥6) so the canary is forced down and Friedman has power (separate work item).

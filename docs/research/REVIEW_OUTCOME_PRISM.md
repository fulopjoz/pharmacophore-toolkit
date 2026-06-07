# PRISM additions — two-pass code review outcome (2026-06-07)

Diff `6d45372..HEAD` reviewed by **two parallel agents** (correctness + scientific-validity) over `color3d.py`, `prism_core.py`, `scorer_prism{,_fixed,_esp}.py`, `bakeoff.py`. No critical/high **correctness** bugs (leakage clean — templates train-only, scaler fit on train; math verified; deterministic). Findings + disposition:

## Fixed (commit 05d1efe)
- **[HIGH, sci] `prism_fixed` was an unfair ablation** — it took a raw mean of unscaled Tanimoto columns while `prism` feeds z-scored features to its logistic, conflating *scaling* with *weighting* (hi-baseline color types implicitly up-weighted). → now z-scores (StandardScaler fit on train) before the equal-weight mean. The `prism − prism_fixed` gap now isolates the learned weighting cleanly.
- **[HIGH, sci] Cross-target reference trap** — `REF_SDF` was global; adding CCR5/CXCR4 later would silently feed CCR2 refs to the ref-based scorers. → per-dataset ref in `DATASET_SPECS`; `load_bench` raises `NotImplementedError` if a dataset lacks a valid ref.
- **[MED, sci] ESP `alpha=0.3` hidden/undocumented** → exposed as `BAKEOFF_ESP_ALPHA` (default 0.3, broader than color; ShaEP/Vainio 2009 spirit), documented.

## Deferred (logged, with reasons)
- **[HIGH, sci] `prism_esp` MUBD gate not yet run** — not a code bug; the full 3-dataset PBS run is the mandatory control. **GATE:** report `prism_esp` as "validated" only if it lifts MUBD over `prism` with non-overlapping CIs; else the extra K columns are overfitting noise.
- **[MED, sci] overfitting on small MUBD train** — control is the MUBD gate above; optional `LogisticRegressionCV` for C deferred (C=1.0 L2 is a reasonable default; don't add complexity pre-result).
- **[MED, corr] Gasteiger recomputed per conformer / per loky dispatch** (pickling strips `_GasteigerCharge`) — results CORRECT (charges are topology-only); marginal perf vs the alignment cost. Defer.
- **[MED, sci] ESP construct validity** — docstrings already scope it honestly as an ESP-*similarity feature* on a color-optimized pose, not field alignment; add Hodgkin/Richards citation in the manuscript.
- **[LOW] fragile `prism_fixed` test assertion; `_row` uses module SEED not a plumbed param; Friedman underpowered at N=3/k=10** — the report already caveats N=3; address when powering the test (CCR5/CXCR4 + MUV/LIT-PCBA) and in write-up.

## Third pass — detailed `/code-review` (3 parallel finders, commit 7f13785)
No new critical/high correctness bugs (4-tuple migration clean, overlap math PSD-safe, leakage clean). Fixed the verified high-value items:
- **DRY config drift** — env constants were copy-pasted across 3 scorer files (+ a separate `_ESP_ALPHA` in color3d) → centralised in `prism_core` (`NCONF/NJOBS/MAX_TEMPLATES/ALPHA`) + a shared `prism_features()` preamble. Removes the "change one, experiments silently non-comparable" footgun.
- **Ablation rigor** — `prism_fixed` now excludes zero-train-variance columns from the equal-weight mean (a color type absent from all train actives got `scale_=1` and leaked its raw value into the mean).
- **Reproducibility** — `seed` plumbed through `feature_matrix`/`_row` (was hardcoded `SEED`); dead `with_esp` branch collapsed.
- **Stats** — Nemenyi `_Q05` extended to k=16 + conservative fallback (was anti-conservative for k>10 via `--methods`).
- **Test gaps closed** — `feature_matrix(with_esp=True)`→K×7 width; `load_bench` raises on missing ref; every `DATASET_SPECS` entry carries a valid ref. **47 tests green.**

Deferred (logged): the 3 prism-family scorers each rebuild the (identical) 3D feature matrix → ~3× conformer work per dataset. A harness-level cache (compute the color feature matrix once per (dataset, templates) and share across prism/prism_fixed/prism_esp; prism_esp's K×7 is a superset of the K×6) would cut the PBS cost ~3×. Not a correctness issue; parallelised across 44 workers on PBS. Also still deferred: Gasteiger-per-conformer recompute, ESP coefficient-interpretation caveat (methods section).

## Next gate
Full PBS bake-off (3 CCR2 datasets × 10 methods incl. prism/prism_fixed/prism_esp). Decisions: (a) `prism` vs `prism_fixed` on MUBD/created → does the **learned weighting** drive the win (CI-separated)? (b) `prism_esp` vs `prism` on MUBD → does **electrostatics** add signal, or overfit?

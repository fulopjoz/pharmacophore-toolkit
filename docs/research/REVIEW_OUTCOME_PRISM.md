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

## Next gate
Full PBS bake-off (3 CCR2 datasets × 11 methods incl. prism/prism_fixed/prism_esp). Decisions: (a) `prism` vs `prism_fixed` on MUBD/created → does the **learned weighting** drive the win (CI-separated)? (b) `prism_esp` vs `prism` on MUBD → does **electrostatics** add signal, or overfit?

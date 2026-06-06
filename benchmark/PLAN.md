# Methods Comparison — Implementation Plan

> **For agentic workers:** use `superpowers:subagent-driven-development` to implement task-by-task. Steps use `- [ ]` checkboxes. Work on branch `s3-enrichment-test` (or a worktree).

**Goal:** build one pluggable harness that scores every candidate method (ours + the literature's) on the same held-out CCR2 actives/decoys split and produces a single ranked comparison table (BEDROC/AUC/EF1% with bootstrap CIs).

**Architecture:** a *scorer registry* — each method is a function `score(actives, decoys, refs) -> per-molecule score`, evaluated by the existing held-out-CV + RDKit `Scoring` (BEDROC/EF) + bootstrap-CI machinery from `../experiments/s3_enrichment_test.py`. Methods that need heavy external installs (TEASER++, OTMol, gWEGA) are isolated behind optional adapters so the core harness runs with rdkit-only.

**Tech stack:** RDKit (features, `ML.Scoring`), scikit-learn (CV, AUC), POT (`ot.gromov` for Fused-GW), optional: teaserpp-python, gwega, OTMol repo.

---

## File structure
- `benchmark/experiments/compare_methods.py` — the harness + scorer registry + master table writer.
- `benchmark/methods/external_methods.md` — per-method: install command, invocation, code URL, benchmarkable? (filled by the parallel-agent audit).
- `benchmark/results/comparison_master.tsv` + `comparison_report.md` — outputs.
- Reuses: `../experiments/s3_enrichment_test.py` (CV + metrics + bootstrap), `tutorials/data/{actives_ccr2_N75.csv, decoys_ccr2_N500.csv, CCR2_reference_ligands.sdf}`.

---

## Task 1 — Harness skeleton + the two baselines (equal-weight & S3)

**Files:** Create `benchmark/experiments/compare_methods.py`; reuse metric/CV code from `../experiments/s3_enrichment_test.py`.

- [ ] **Step 1 — Write the failing test** (`benchmark/experiments/test_compare.py`):
```python
from compare_methods import REGISTRY, evaluate
def test_registry_has_baselines():
    assert {"equal_weight", "s3_weighted"} <= set(REGISTRY)
def test_evaluate_returns_metrics():
    m = evaluate("equal_weight")
    assert {"AUC", "BEDROC", "EF1%"} <= set(m) and 0 <= m["AUC"] <= 1
```
- [ ] **Step 2 — Run, verify it fails** (`pytest benchmark/experiments/test_compare.py -q`) → FAIL (no module).
- [ ] **Step 3 — Implement** `compare_methods.py`: import `metrics`, `per_type_sim_to_refs`, `load_counts`, `ref_profiles` from `s3_enrichment_test`; a `REGISTRY = {}` dict; `register(name)` decorator; `equal_weight` and `s3_weighted` scorers (lift arms A and B from the Tier-1 script); `evaluate(name)` = 5-fold held-out CV → pooled out-of-fold scores → `metrics`.
- [ ] **Step 4 — Run, verify pass.**
- [ ] **Step 5 — Commit** `git add benchmark/experiments/ && git commit -m "benchmark: harness skeleton + equal/S3 baselines"`.

## Task 2 — Register the rdkit-only literature scorers (no external installs)

**Files:** Modify `benchmark/experiments/compare_methods.py`.

- [ ] **Step 1 — Write failing tests** asserting `REGISTRY` contains `differential_mmfp`, `fused_gw`, `shape_combo_rdkit`.
- [ ] **Step 2 — Run, verify fail.**
- [ ] **Step 3 — Implement three scorers:** (a) `differential_mmfp` = Hutter weight = freq(actives)−freq(decoys) over Morgan bits (the literature S3, DOI 10.1021/acs.jcim.2c00242); (b) `fused_gw` = pose-free Fused-Gromov-Wasserstein similarity via `ot.gromov.fused_gromov_wasserstein2` over atom-coord distance matrices + RDKit P4 feature labels (DOI 10.48550/arxiv.1811.02834); (c) `shape_combo_rdkit` = `rdShapeAlign.AlignMol` combo Tanimoto vs the 5 refs (the current production proxy).
- [ ] **Step 4 — Run, verify pass.**
- [ ] **Step 5 — Commit** `benchmark: add MMFP, Fused-GW, RDKit-combo scorers`.

## Task 3 — Master comparison table + report

**Files:** Modify `compare_methods.py`; outputs to `benchmark/results/`.

- [ ] **Step 1 — Write failing test:** `test_master_table_written()` asserts `results/comparison_master.tsv` has one row per registered method with AUC/BEDROC/EF1% + ΔBEDROC-vs-s3 + bootstrap-CI columns.
- [ ] **Step 2 — Run, verify fail.**
- [ ] **Step 3 — Implement** a `main()` that evaluates every `REGISTRY` method, computes ΔBEDROC and ΔAUC vs `s3_weighted` with the bootstrap-CI helper from the Tier-1 script, writes `comparison_master.tsv` + a markdown `comparison_report.md` ranking methods, flagging any that beat S3 (CI excludes 0).
- [ ] **Step 4 — Run** `python benchmark/experiments/compare_methods.py`; verify the table.
- [ ] **Step 5 — Commit** `benchmark: master comparison table + report`.

## Task 4 — External-install adapters (optional, gated on the agent audit)

**Files:** `benchmark/methods/external_methods.md` (audit output); add adapters to `compare_methods.py` guarded by `try: import` so the core harness still runs without them.

- [ ] **Step 1** — From the parallel-agent audit (`external_methods.md`), pick the methods that install cleanly (likely: **OTMol**, **gWEGA**, **SHAFTS** if available; **TEASER++** as a pose-seed not a scorer).
- [ ] **Step 2** — For each, write an adapter `score_<method>(...)` behind `try/except ImportError` that returns per-molecule scores; register it only if importable.
- [ ] **Step 3** — Re-run `main()`; the master table gains rows for whichever externals are installed.
- [ ] **Step 4 — Commit** `benchmark: external-method adapters (OTMol/gWEGA/...)`.

## Decision gate
A literature method is "worth adopting over our S3 solution" only if its **held-out BEDROC exceeds S3 with a bootstrap CI excluding 0**. Record per-method ΔBEDROC + CI in `comparison_report.md`. (Same gate as Tier-1.)

## Self-review notes
- Spec coverage: papers documented (Task 0/registry ✓), harness (T1), literature scorers (T2), table (T3), externals (T4). ✓
- All scorers share ONE held-out CV + metric path (no per-method leakage). ✓
- External installs are optional/guarded so the core comparison always runs. ✓

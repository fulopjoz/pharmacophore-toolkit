# Code Review — Consensus Pharmacophore Optimization Toolkit

**Date:** 2026-06-06 · **Scope:** all 29 core `pharmacophore/*.py` modules (~15.2k LOC) + tests · **Method:** 6 parallel reviewers (5 code clusters via the code-reviewer framework + 1 methodology/literature reviewer via scite; all DOIs Crossref-verified).
**Bottom line:** the *construction* of the consensus pharmacophore is methodologically sound and literature-grounded, but the codebase has **two systemic problems that invalidate reported enrichment numbers** — (A) **pervasive overfitting / no held-out validation**, and (B) **a silent active/decoy cache collision** — plus several real geometry/optimizer bugs. **Do not trust any current AUC/BEDROC figure until A and B are fixed.**

---

## Executive summary
- **Risk level: HIGH** for scientific validity; MEDIUM for engineering robustness.
- **Top 5 must-fix:** (1) no held-out/scaffold split — every optimizer selects on the same data it reports; (2) `_prepare_molecules` cache-key collision can swap actives↔decoys; (3) `evaluate_feature_subset` ignores the passed features in reference mode → HypoGen selection is meaningless; (4) leakage in `compute_reference_weights`/`compute_zscore_params`/`LearnedScorer` scaler; (5) `orthogonal_procrustes` allows reflections (destroys chirality).
- **Convergent evidence (high confidence):** the overfitting finding was raised independently by **4 of 6** reviewers; the cache-collision by **2 of 6** — these are not single-agent guesses.

## CRITICAL findings (correctness or scientific validity; fix before trusting results)

| # | file:line | issue | why it corrupts results | fix |
|---|-----------|-------|-------------------------|-----|
| C1 | **`evaluation.py:392–399`** | **Active/decoy cache collision** — key is `(SMILES(mols[0]), len(mols), …)`; if `actives[0]`≡`decoys[0]` and equal lengths, the 2nd call returns the 1st's conformers | decoys silently scored as actives → AUC≈0.5 or inflated, no error. *(found by 2 reviewers; the `test_evaluation.py` fixture literally triggers it)* | hash ALL smiles: `md5('|'.join(sorted(MolToSmiles(m))))` |
| C2 | **all optimizers + greedy + hypogen** (`auto_optimizer.py:666`, `optuna_optimizer.py:172`, `combo_optimizer.py:707`, `combinatorial_optimizer.py:288`, `doe_optimization.py:273`, `optimization.py:382`, `greedy_selector.py:263–310`, `hypogen_optimizer.py:428–551`) | **Selection-on-the-evaluation-set overfitting** — every optimizer maximizes AUC/BEDROC on the *same* 74 actives/498 decoys it then reports; no held-out, no nested CV, 50–500+ trials → multiple-comparison inflation | reported enrichment is optimistic by construction (lit: +5–15 AUC pts at N=74); will not generalize prospectively | stratified **scaffold split**: optimize on train, report once on held-out; nested CV for selection |
| C3 | **`evaluation.py:972–1067`** | **`evaluate_feature_subset` ignores the `features` arg** in default `scoring_mode='reference'` — scores against `_prepared_refs` regardless of subset | HypoGen Phase-1/2/3 give *identical* scores to all hypotheses → the entire feature-subset selection is meaningless | restrict refs/consensus to the subset, or force `consensus_mol` mode for subset eval |
| C4 | **`evaluation.py:620–739`** + **`learned_scoring.py:85`** | **Leakage in data-dependent normalization** — `compute_reference_weights`/`compute_zscore_params` fit on the full actives+decoys set later evaluated; `LearnedScorer` `StandardScaler.fit_transform(X)` before CV | inflates AUC/BEDROC (0.02–0.05); the `cv_auc` overfitting-detector is itself contaminated | derive on a held-out calibration set; use `Pipeline([scaler, lr])` inside `cross_val_score` |
| C5 | **`point_cloud_alignment.py:168`** | **`orthogonal_procrustes` returns reflections** (`det=−1`) for (near-)coplanar feature sets | mirrors the molecule (chirality destroyed) → physically meaningless alignment & scores | Kabsch det-correction: `R = U·diag(1,1,sign(det(UVᵀ)))·Vᵀ` |
| C6 | **`mol_converter.py:386–390`** | **`convert_with_metadata` labels wrong atoms** — assumes 1 atom/feature but color features expand to multi-atom fragments (Donor→NH3=4 atoms) | `FeatureType` on wrong atoms; `IndexError` for >~3 features | track `atom_offset` across fragments; label only the heavy atom |

## HIGH findings (real bugs, narrower blast radius)
- **`ensemble_consensus.py:191–212`** — occurrence threshold counts raw feature *points*, not unique *molecules* (a mol with 2 donors double-counts) → phantom consensus features. Fix: propagate `mol_idx`.
- **`optuna_optimizer.py:289–291`** — `max(study.trials, key=t.values[0])` crashes (`t.values=None`) if any trial failed/pruned. Fix: filter `state==COMPLETE`.
- **`combo_optimizer.py:610–613, 860–878`** — mutable `self._prepared_refs` overwritten per trial (best-params refs lost); post-hoc retest always uses deprecated `consensus_mol` mode → incommensurable AUC replaces best. Fix: local refs var; gate retest on `scoring_mode`.
- **`hungarian_matching.py:132`** — `pair_cost ≥ 0.95·dummy_cost` falsely marks real pairs at `max_distance` as unmatched when `type_weight=0`. Fix: detect padding by `r≥n_a or c≥n_b`, drop the 0.95 heuristic.
- **`point_cloud_alignment.py:246–265`** — surplus features double-penalized (`gap_factor` *and* inflated `type_match_score` denominator) when `n_a≠n_b`.
- **`shape_experiments.py:512–527`** — `AlignMol` mutates the probe conformer **in-place** across the reference loop → ref 2+ aligned to already-moved coords, order-dependent. Fix: copy per ref.
- **`ot_scoring.py:114–124`** — scipy fallback vs POT give systematically different distances for unequal feature-set sizes (34% in the example) → scores depend on whether POT is installed.
- **`screening_metrics.py:96`** — docstring claims BEDROC "random = 0.5"; true baseline ≈ active-ratio Ra (the *number* via RDKit is correct, the *interpretation* doc is wrong → users misread below/above-random).
- **`evaluation.py:801–829`** — `hybrid` mode: `scores_3d` drops embedding-failures so it misaligns with `y_true` (length N) → broadcast error → silent `roc_auc=0.5` fallback.
- **`drugex_reward.py:71`** — `lru_cache` skips `_score_uncached` on hits, so `_n_valid`/`_n_errors` validity counters are wrong (reports n_valid=1 after 1000 hits) → misleading RL validity monitor.
- **`pharmacophore.py:530–535, 65–70`** — bare `except: pass` swallows SMARTS errors (silent missing feature classes); module-level `global` state (`phrase`/`pharmacophore`/`feat_factory`) is a race under `n_jobs>1`.
- **`mol_converter.py:116–119`** — `EmbedMolecule` return unchecked → `GetConformer()` raises on embedding failure (deterministic with fixed seed).
- **`clustering_algorithms.py:44, 188`** — `n_clusters=0` on degenerate clouds → `KMeans` ValueError (and the fallback skips the occurrence filter); `n_molecules=0` → every cluster passes.

## Methodology verdict (scite, DOIs verified)
- **Construction: ADEQUATE / sound.** Hierarchical clustering + occurrence threshold + per-type spatial tolerances is standard consensus-pharmacophore practice (PharmaGist lineage); ensemble stability voting is conceptually aligned with **Madzhidov 2020** (10.3390/molecules25020385); the reference-direct scoring mode and the RRF fusion (**Cormack 2009**, 10.1145/1571941.1572114) are reasonable, and the code correctly deprecates the anti-discriminative `consensus_mol` path.
- **Validation: INSUFFICIENT.** No held-out set, no **scaffold split** (CCR2 actives share indole/benzimidazole/piperidine motifs → random splits inflate), CV covers only the logistic layer. Lit: **Sieg 2019** (10.1021/acs.jcim.8b00712), **Landrum SIMPD 2023** (10.1186/s13321-023-00787-9).
- **Decoy bias: UNAUDITED.** Property-matched decoys carry analog bias (+15–25% for GPCRs). There is a **CCR2-specific** unbiased benchmark to adopt: **Xia 2018** (10.1021/acs.jcim.8b00004); also **Réau 2018** (10.3389/fphar.2018.00011), **Rohrer-Baumann MUV** (10.1021/ci8002649, max-Tc<0.6 rule).
- **Overfit risk: HIGH** — three compounding loops (greedy selection → Bayesian/NSGA-II → learned scorer) all on the same 74/498 set; with ~15 actives/fold, AUC CI ≈ ±0.08 (**Baber 2005** 10.1021/ci050296y).

## Test-coverage gaps
- **No test files:** `combinatorial_optimizer`, `auto_strategy`, `pharm2d_scoring`, `pharm3d_scoring`, `shape_experiments`.
- **Tests that can't catch the critical bugs:** no test asserts two *different* feature subsets give different scores (would expose C3); `test_evaluation.py` fixture puts the same molecule in actives *and* decoys (triggers C1 but never asserts on it); `test_drugex_reward.py` doesn't check `n_valid` after cache hits.

## Prioritized roadmap
1. **Fix C1 (cache collision) immediately** — it can corrupt *every* evaluation; add the regression test from the fixture.
2. **Add a scaffold-stratified held-out split** at the data layer; route all optimizers/selectors/learned-scorer through train-only; report once on held-out (fixes C2 + C4). Add nested CV for greedy selection.
3. **Fix C3** — make `evaluate_feature_subset` actually use the subset; add the "different subsets → different scores" test.
4. **Fix C5/C6** (Procrustes reflection; metadata atom mislabel) — geometry/labeling correctness.
5. **Audit decoys** vs actives (max Morgan-Tc), or adopt the Xia 2018 CCR2 MUBD set; report BEDROC with bootstrap CI.
6. Clear the HIGH list (optuna crash, combo mutable-refs/retest, hungarian threshold, shape_experiments in-place, ot scipy/POT, hybrid misalign, drugex counters, bare-except/globals).
7. Backfill tests for the untested modules.

## References (Crossref-verified 2026-06-06)
Madzhidov 2020 10.3390/molecules25020385 · Cormack RRF 2009 10.1145/1571941.1572114 · Sieg 2019 10.1021/acs.jcim.8b00712 · Réau 2018 10.3389/fphar.2018.00011 · Xia 2018 (CCR2 MUBD) 10.1021/acs.jcim.8b00004 · Rohrer-Baumann MUV 2009 10.1021/ci8002649 · Baber 2005 10.1021/ci050296y · Giangreco 2021 10.1021/acs.jcim.1c00866 · Landrum SIMPD 2023 10.1186/s13321-023-00787-9.

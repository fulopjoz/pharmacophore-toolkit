# CCR2 Methods Comparison — Findings

**Held-out 5-fold CV** on 74 actives + 498 decoys vs the 5 CCR2 references. Primary metric **BEDROC(α=20)**; Δ vs our **S3** baseline with bootstrap 95% CI. Numbers below are AFTER the `/code-review` fixes (see §Code-review).

| method | AUC | BEDROC | ΔBEDROC vs S3 (95% CI) | beats S3? |
|--------|-----|--------|------------------------|-----------|
| `differential_mmfp` (2D Morgan, Hutter) | 0.993 | 0.988 | +0.373 [+0.254,+0.499] | ✅ (see caveat) |
| **`fused_gw` (pose-free OT, cross-domain)** | **0.940** | **0.894** | **+0.279 [+0.151,+0.406]** | ✅ |
| `shape_combo_rdkit` (production ROCS proxy) | 0.806 | 0.656 | +0.039 [−0.098,+0.179] | no (ties S3) |
| `s3_weighted` (**ours**) | 0.874 | 0.612 | — | baseline |
| `equal_weight` | 0.519 | 0.193 | −0.404 [−0.553,−0.234] | no |

## Interpretation (scientific-critical-thinking)
- **Fused-GW is the genuine winner among shape-aware methods.** A pose-free optimal-transport transfer (graph topology + atom features, no alignment) beats our S3 on held-out BEDROC with a CI excluding 0, runs in ~4 s, and got *stronger* after the M-normalization fix (0.847→0.894). This validates the cross-domain research direction — borrowed-algorithm + open-source genuinely outperforms our hand-built S3 color weighting here.
- **DMMFP's near-perfect AUC (0.993) is almost certainly DECOY-BIAS, not a real win.** It is a 2D Morgan-fingerprint classifier on **property-matched (DUD-E-style) decoys**, the textbook setup where 2D fingerprints exploit subtle 2D differences that don't generalize. Tellingly, our own stable-core literature includes *"Do molecular fingerprints identify diverse active drugs? (No)"* (10.3390/ph17080992) warning of exactly this. **Do not read DMMFP as beating the 3D methods** — it answers a different (easier, biased) question. A scaffold-split or DEKOIS-style decoy set would likely collapse it.
- **Our S3 holds a respectable middle:** clearly beats equal-weight color (now honestly near-random at 0.519 after the leakage fix) and **ties the production shape+color proxy** (`shape_combo_rdkit`, CI includes 0) — i.e. the cheap discrimination-weighted color matches full 3D ROCS-combo on this benchmark.
- **EF1% saturates at 7.73** (top 1% of 572 ≈ 6 molecules) for the top methods — use BEDROC/AUC for separation, not EF1%.

## Caveats (honest scope)
- Decoys are **property-matched, not assay-confirmed** → relative enrichment only; 2D methods (DMMFP) are most flattered by this.
- **43/74 actives overlap the DrugEx finetune set** — fine for a retrospective scoring comparison, but any production claim must re-derive on a training-disjoint split.
- `fused_gw` and `shape_combo` use single conformers/topology — not a full ROCS ensemble; treat as fast proxies.

## Decision
**Pursue Fused-GW (pose-free OT)** as the lead cross-domain method — it beats our S3 cleanly and is fast/differentiable-adjacent (Sinkhorn variants → GPU/RL). Keep S3 as the color-weighting baseline. Treat DMMFP as a *decoy-bias control*, not a contender.

## Code-review (`/code-review`, xhigh) — 2 material + 3 robustness fixes applied
| fix | file | what |
|-----|------|------|
| **Fused-GW alpha was ~10:1 feature-biased** | `scorer_fused_gw.py` | normalized M to [0,1] → `alpha=0.5` genuinely balanced (improved BEDROC 0.847→0.894) |
| **Test-set leakage** in the similarity scale | `harness.py` `_p4_sim_to_refs` | scale now from references only (never the ranked molecules) → strict held-out CV |
| beats_s3 string-parse edge case | `compare.py` | numeric `lo>0` gate (not parsing a formatted CI string) |
| empty-refs crash | `scorer_shape_combo.py` | guard returns zeros if all refs fail to parse |
| empty-bootstrap / partial-OOF crashes | `harness.py` | NaN post-condition + empty-diffs guard; FGW similarity clipped to [0,1] |

**Deferred (logged, low-priority):** `morgan()` could use `DataStructs.ConvertToNumpyArray`; `shape_combo` makes 5× `RWMol` copies/molecule and uses flat cache keys; module-scoped pytest fixtures for fused_gw; `s3_enrichment_test.py` duplicates `harness` metrics (DRY) — should import from `harness`. None affect correctness of the reported numbers.

# CCR2 Methods Comparison — Findings

**Held-out 5-fold CV** on 74 actives + 498 decoys vs the 5 CCR2 references. Primary metric **BEDROC(α=20)**; Δ vs our **S3** baseline with bootstrap 95% CI. Numbers below are AFTER the `/code-review` fixes (see §Code-review).

| method | AUC | BEDROC | ΔBEDROC vs S3 (95% CI) | beats S3? |
|--------|-----|--------|------------------------|-----------|
| `differential_mmfp` (2D Morgan, Hutter) | 0.993 | 0.988 | +0.373 [+0.254,+0.499] | ✅ (see caveat) |
| **`fused_gw` (pose-free OT, cross-domain)** | **0.940** | **0.894** | **+0.279 [+0.151,+0.406]** | ✅ |
| `shape_combo_rdkit` (production ROCS proxy) | 0.806 | 0.656 | +0.039 [−0.098,+0.179] | no (ties S3) |
| `s3_weighted` (**ours**) | 0.874 | 0.612 | — | baseline |
| `equal_weight` | 0.519 | 0.193 | −0.404 [−0.553,−0.234] | no |

## ⚠️ What each method actually uses (2D vs 3D) — read this first
| method | representation | uses our rdShapeAlign shape? |
|--------|----------------|------------------------------|
| `differential_mmfp` | 2D Morgan bits | no |
| `fused_gw` | **2D molecular graph** (`Chem.GetDistanceMatrix` bond-path distances) + atom types — *no 3D conformers* | **no** |
| `s3_weighted` / `equal_weight` (ours) | **2D** P4 feature *counts* similarity to refs | no |
| `shape_combo_rdkit` | **3D shape+color** (`rdShapeAlign.AlignMol`, ETKDG conformer) | **YES — the only one** |

**So 4 of 5 methods are 2D ligand descriptors; only `shape_combo_rdkit` is true 3D shape.** This is essential for reading the table.

## Interpretation (scientific-critical-thinking)
- **No genuinely shape-based method beat our S3.** The *only* true 3D shape+color method (`shape_combo_rdkit`, = our production rdShapeAlign proxy) merely **ties** S3 (ΔBEDROC CI includes 0). Both methods that "beat" S3 are **2D** (DMMFP Morgan, fused_gw graph-OT).
- **Fused-GW is a 2D graph-topology method, NOT a shape method** (it uses `GetDistanceMatrix`, "no 3-D conformers"). It is pose-free *because* it never uses 3D shape. It therefore carries the **same decoy-bias suspicion as DMMFP**: 2D descriptors exploit subtle 2D differences in property-matched (DUD-E-style) decoys that need not generalize. Its win over S3 (+0.279) should be read as "2D-graph-OT separates these particular decoys well," **not** "OT shape-matching beats our shape scoring."
- **The real shape-via-OT test is not yet run.** A faithful comparison would feed FGW a **3D Euclidean** atom-distance matrix from a conformer (so it compares *shape*), not the 2D bond-path graph. That is the deferred "promote Fused-GW to 3D" step — only then is it apples-to-apples with `shape_combo_rdkit` / rdShapeAlign.
- **DMMFP's near-perfect AUC (0.993) is almost certainly DECOY-BIAS, not a real win.** It is a 2D Morgan-fingerprint classifier on **property-matched (DUD-E-style) decoys**, the textbook setup where 2D fingerprints exploit subtle 2D differences that don't generalize. Tellingly, our own stable-core literature includes *"Do molecular fingerprints identify diverse active drugs? (No)"* (10.3390/ph17080992) warning of exactly this. **Do not read DMMFP as beating the 3D methods** — it answers a different (easier, biased) question. A scaffold-split or DEKOIS-style decoy set would likely collapse it.
- **Our S3 holds a respectable middle:** clearly beats equal-weight color (now honestly near-random at 0.519 after the leakage fix) and **ties the production shape+color proxy** (`shape_combo_rdkit`, CI includes 0) — i.e. the cheap discrimination-weighted color matches full 3D ROCS-combo on this benchmark.
- **EF1% saturates at 7.73** (top 1% of 572 ≈ 6 molecules) for the top methods — use BEDROC/AUC for separation, not EF1%.

## Caveats (honest scope)
- Decoys are **property-matched, not assay-confirmed** → relative enrichment only; 2D methods (DMMFP) are most flattered by this.
- **43/74 actives overlap the DrugEx finetune set** — fine for a retrospective scoring comparison, but any production claim must re-derive on a training-disjoint split.
- `fused_gw` and `shape_combo` use single conformers/topology — not a full ROCS ensemble; treat as fast proxies.

## Decision (revised after the 2D/3D clarification)
- **Do NOT yet conclude OT beats our shape scoring** — the current `fused_gw` is 2D-graph, not shape, and its win is decoy-bias-suspect (like DMMFP).
- **Run the faithful test:** a **3D** Fused-GW (Euclidean atom-distance matrix from an ETKDG conformer + pharmacophore atom features) vs `shape_combo_rdkit` / rdShapeAlign on the same split. Only that answers "does OT shape-matching beat our shape+color?".
- **Harden the benchmark against 2D decoy-bias:** add a **scaffold-split** or DEKOIS-style decoy set; under it, expect DMMFP (and possibly 2D fused_gw) to collapse while genuine 3D methods hold.
- **Our S3 stands respectably:** beats equal-weight and *ties* full 3D ROCS-combo at a fraction of the cost.

## Code-review (`/code-review`, xhigh) — 2 material + 3 robustness fixes applied
| fix | file | what |
|-----|------|------|
| **Fused-GW alpha was ~10:1 feature-biased** | `scorer_fused_gw.py` | normalized M to [0,1] → `alpha=0.5` genuinely balanced (improved BEDROC 0.847→0.894) |
| **Test-set leakage** in the similarity scale | `harness.py` `_p4_sim_to_refs` | scale now from references only (never the ranked molecules) → strict held-out CV |
| beats_s3 string-parse edge case | `compare.py` | numeric `lo>0` gate (not parsing a formatted CI string) |
| empty-refs crash | `scorer_shape_combo.py` | guard returns zeros if all refs fail to parse |
| empty-bootstrap / partial-OOF crashes | `harness.py` | NaN post-condition + empty-diffs guard; FGW similarity clipped to [0,1] |

**Deferred (logged, low-priority):** `morgan()` could use `DataStructs.ConvertToNumpyArray`; `shape_combo` makes 5× `RWMol` copies/molecule and uses flat cache keys; module-scoped pytest fixtures for fused_gw; `s3_enrichment_test.py` duplicates `harness` metrics (DRY) — should import from `harness`. None affect correctness of the reported numbers.

# External Methods — Code Availability & Benchmarkability Audit

5 parallel agents audited whether each stable-core method can be **installed and run on CCR2 SMILES/SDF** to score the 74 actives + 500 decoys vs the 5 references. Verdict drives which methods enter the comparison harness.

## ✅ Benchmarkable now (enter the harness)

| Method | What it is | Install | Effort | Role in comparison |
|--------|-----------|---------|--------|--------------------|
| **DMMFP** (Hutter 2022, 10.1021/acs.jcim.2c00242) | `weight = freq(actives) − freq(decoys)` over fingerprint bits — **the published version of our S3** | none (~40 lines RDKit+numpy) | **low** | **the direct head-to-head for S3** (2D MACCS/Morgan vs our P4-color S3) |
| **rdShapeAlign combo** (RDKit, BSD) | Gaussian shape+color Tanimoto vs refs — the WEGA/ROCS stand-in; **already our production proxy** | already installed | low | shape+color baseline (= current pipeline) |
| **Fused Gromov-Wasserstein** (POT, MIT; FGW 1811.02834) | pose-free graph-OT similarity (intra-distance matrices + typed-feature cost) | `pip install pot` | **med** | the cross-domain **pose-free** entrant (no alignment) |
| **shepherd-score** (Coley group, MIT) | 3D shape + pharmacophore (+ESP if xTB) similarity | `pip install shepherd-score` (ESP needs xTB) | med | 3D shape+pharm comparator (superset of ROCS color) |
| **CDPKit pharm-align** (LGPL; G3PS-family, same author Seidel) | matched-feature-pair + geometric pharmacophore alignment score | `pip install cdpkit` | med-high | pharmacophore-only (no shape) comparator |

## ⚠️ Pose-seed only (not a standalone scorer)
| Method | Note |
|--------|------|
| **TEASER++** (MIT, 10.1109/TRO.2020.3033695) | builds from source (Eigen3+C++, **high**); usable as a *color-aware pose seed* for RDKit/CDPKit (needs a type-matched correspondence generator), **not** a scorer. Defer to the "pose-seed" experiment, not the enrichment harness. |

## ❌ Excluded (no usable code / wrong quantity / blockers)
| Method | Why excluded |
|--------|--------------|
| **gWEGA** (10.1002/jcc.23603) | **no public code** (proprietary, GPU). Use rdShapeAlign / ROSHAMBO as the open Gaussian-shape stand-in. |
| **SHAFTS** (10.1021/ci200060s) | **no source** — approval-gated GUI binary + web server only; not scriptable/batchable. |
| **OTMol** (2509.01550, MIT) | public code, but computes **atom-assignment RMSD, not pharmacophore similarity** — semantically wrong for VS enrichment across diverse scaffolds. The POT-FGW entry is the correct OT representative. |
| **Go-ICP** (TPAMI 2015) | **GPL-3.0** + 20–25 s/run + broken Python-3.12 wrappers. Exclude. |
| **ShEPhERD** (2411.04130) | **generative** diffusion model, not a scorer (its sister `shepherd-score` IS the scorer → listed above). |
| **ROSHAMBO** (10.1021/acs.jcim.4c01225, GPL) | best academic ROCS-equivalent but **high install** (compile old RDKit from source + CUDA GPU). Optional stretch goal. |

## Net benchmark set (final)
`equal_weight` · `s3_weighted` (ours) · **`differential_mmfp`** (literature S3) · `shape_combo_rdkit` (production) · **`fused_gw`** (pose-free OT) · *(optional)* `shepherd_score`, `cdpkit_pharm`.
All scored on the **same held-out CCR2 split** with the same metrics (BEDROC α=20, AUC, EF1% + bootstrap CIs) as the Tier-1 test. **Decision gate:** a method beats our S3 only if held-out ΔBEDROC > 0 with a bootstrap CI excluding 0.

> Key practical finding: the two highest-value comparisons need **zero or one install** — **DMMFP** (the literature S3, ~40 lines) and **Fused-GW** (`pip install pot`). The commercial/heavy methods (gWEGA, SHAFTS, ROSHAMBO, TEASER++) are either unavailable or pose-seed-only, so the honest open-source comparison is fully achievable.

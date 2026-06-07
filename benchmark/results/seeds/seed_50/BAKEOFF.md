# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 50, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| prism | 0.988 [0.949,1.000] | 0.996 | 7.68 | +0.307 [+0.108,+0.554] |
| prism_esp | 0.984 [0.940,1.000] | 0.995 | 7.68 | +0.302 [+0.102,+0.552] |
| differential_mmfp | 0.974 [0.902,1.000] | 0.979 | 7.68 | +0.292 [+0.098,+0.535] |
| prism_fixed | 0.964 [0.885,1.000] | 0.986 | 7.68 | +0.280 [+0.088,+0.531] |
| learned_scorer | 0.899 [0.740,0.987] | 0.946 | 7.68 | +0.216 [-0.015,+0.480] |
| rdshape_ensemble | 0.881 [0.684,0.968] | 0.898 | 7.68 | +0.195 [-0.037,+0.447] |
| pharm2d | 0.771 [0.522,0.916] | 0.881 | 7.68 | +0.091 [-0.173,+0.338] |
| shape_combo_rdkit | 0.740 [0.480,0.902] | 0.769 | 7.68 | +0.057 [-0.232,+0.334] |
| s3_weighted | 0.715 [0.422,0.887] | 0.891 | 7.68 | — (baseline) |
| equal_weight | 0.157 [0.003,0.305] | 0.395 | 3.84 | -0.570 [-0.825,-0.271] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.809 [0.631,0.956] | 0.950 | 35.41 | +0.472 [+0.124,+0.773] |
| prism_fixed | 0.377 [0.183,0.570] | 0.826 | 10.12 | +0.038 [-0.275,+0.356] |
| prism_esp | 0.360 [0.169,0.560] | 0.691 | 10.12 | +0.021 [-0.307,+0.340] |
| s3_weighted | 0.345 [0.142,0.543] | 0.761 | 15.18 | — (baseline) |
| prism | 0.263 [0.104,0.446] | 0.679 | 5.06 | -0.076 [-0.372,+0.208] |
| learned_scorer | 0.190 [0.027,0.364] | 0.695 | 10.12 | -0.155 [-0.438,+0.152] |
| shape_combo_rdkit | 0.008 [0.000,0.028] | 0.357 | 0.00 | -0.336 [-0.538,-0.132] |
| equal_weight | 0.003 [0.000,0.009] | 0.297 | 0.00 | -0.341 [-0.541,-0.137] |
| rdshape_ensemble | 0.001 [0.000,0.004] | 0.390 | 0.00 | -0.344 [-0.542,-0.140] |
| pharm2d | 0.000 [0.000,0.000] | 0.184 | 0.00 | -0.345 [-0.543,-0.142] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.987 [0.978,0.992] | 0.946 | 3.69 | +0.544 [+0.454,+0.644] |
| prism_esp | 0.953 [0.927,0.971] | 0.900 | 3.69 | +0.510 [+0.419,+0.611] |
| prism | 0.947 [0.920,0.968] | 0.901 | 3.69 | +0.505 [+0.415,+0.606] |
| prism_fixed | 0.922 [0.882,0.953] | 0.858 | 3.69 | +0.479 [+0.383,+0.583] |
| learned_scorer | 0.517 [0.412,0.609] | 0.653 | 3.36 | +0.073 [-0.059,+0.211] |
| s3_weighted | 0.441 [0.339,0.534] | 0.598 | 2.68 | — (baseline) |
| rdshape_ensemble | 0.383 [0.279,0.482] | 0.497 | 3.36 | -0.059 [-0.197,+0.084] |
| shape_combo_rdkit | 0.369 [0.264,0.465] | 0.499 | 2.68 | -0.073 [-0.210,+0.069] |
| pharm2d | 0.225 [0.136,0.316] | 0.442 | 2.01 | -0.219 [-0.344,-0.077] |
| equal_weight | 0.158 [0.105,0.241] | 0.457 | 0.00 | -0.273 [-0.391,-0.148] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=22.636, p=0.007; Nemenyi CD=7.822 over N=3 datasets, k=10 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 1.67 |
| prism_esp | 2.33 |
| prism | 3.00 |
| prism_fixed | 3.33 |
| learned_scorer | 5.33 |
| s3_weighted | 6.33 |
| rdshape_ensemble | 7.33 |
| shape_combo_rdkit | 7.67 |
| pharm2d | 8.67 |
| equal_weight | 9.33 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (7.82 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

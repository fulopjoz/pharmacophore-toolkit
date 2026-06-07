# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 43, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| prism_esp | 0.976 [0.912,1.000] | 0.995 | 7.68 | +0.242 [+0.061,+0.495] |
| prism | 0.974 [0.906,1.000] | 0.994 | 7.68 | +0.241 [+0.059,+0.493] |
| differential_mmfp | 0.974 [0.901,1.000] | 0.979 | 7.68 | +0.238 [+0.061,+0.480] |
| prism_fixed | 0.900 [0.707,0.999] | 0.983 | 7.68 | +0.172 [-0.089,+0.449] |
| rdshape_ensemble | 0.849 [0.633,0.953] | 0.879 | 7.68 | +0.108 [-0.136,+0.376] |
| learned_scorer | 0.825 [0.585,0.997] | 0.952 | 3.84 | +0.109 [-0.210,+0.408] |
| pharm2d | 0.793 [0.545,0.922] | 0.890 | 7.68 | +0.056 [-0.205,+0.314] |
| s3_weighted | 0.755 [0.479,0.919] | 0.923 | 7.68 | — (baseline) |
| shape_combo_rdkit | 0.642 [0.346,0.854] | 0.738 | 7.68 | -0.086 [-0.409,+0.209] |
| equal_weight | 0.263 [0.017,0.475] | 0.460 | 3.84 | -0.501 [-0.793,-0.125] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.863 [0.702,0.982] | 0.921 | 30.35 | +0.574 [+0.293,+0.830] |
| s3_weighted | 0.294 [0.089,0.499] | 0.682 | 15.18 | — (baseline) |
| prism_esp | 0.145 [0.004,0.316] | 0.531 | 5.06 | -0.146 [-0.438,+0.155] |
| prism | 0.102 [0.003,0.232] | 0.584 | 0.00 | -0.191 [-0.452,+0.075] |
| prism_fixed | 0.063 [0.005,0.149] | 0.624 | 0.00 | -0.223 [-0.430,-0.039] |
| equal_weight | 0.033 [0.000,0.080] | 0.378 | 0.00 | -0.261 [-0.481,-0.046] |
| learned_scorer | 0.032 [0.002,0.085] | 0.627 | 0.00 | -0.256 [-0.468,-0.057] |
| shape_combo_rdkit | 0.001 [0.000,0.003] | 0.322 | 0.00 | -0.287 [-0.498,-0.087] |
| rdshape_ensemble | 0.001 [0.000,0.002] | 0.385 | 0.00 | -0.287 [-0.498,-0.088] |
| pharm2d | 0.000 [0.000,0.000] | 0.147 | 0.00 | -0.288 [-0.499,-0.089] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.987 [0.978,0.992] | 0.944 | 3.69 | +0.565 [+0.470,+0.662] |
| prism_esp | 0.951 [0.923,0.970] | 0.890 | 3.69 | +0.528 [+0.431,+0.627] |
| prism | 0.945 [0.915,0.966] | 0.891 | 3.69 | +0.523 [+0.424,+0.622] |
| prism_fixed | 0.914 [0.871,0.945] | 0.845 | 3.69 | +0.491 [+0.392,+0.594] |
| learned_scorer | 0.574 [0.473,0.664] | 0.658 | 3.36 | +0.152 [+0.022,+0.287] |
| s3_weighted | 0.421 [0.323,0.518] | 0.592 | 2.35 | — (baseline) |
| rdshape_ensemble | 0.397 [0.296,0.499] | 0.500 | 3.36 | -0.024 [-0.162,+0.120] |
| shape_combo_rdkit | 0.368 [0.270,0.471] | 0.499 | 2.68 | -0.053 [-0.189,+0.086] |
| pharm2d | 0.232 [0.143,0.325] | 0.455 | 2.01 | -0.188 [-0.315,-0.058] |
| equal_weight | 0.161 [0.104,0.238] | 0.463 | 0.67 | -0.255 [-0.370,-0.135] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=21.327, p=0.011; Nemenyi CD=7.822 over N=3 datasets, k=10 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 1.67 |
| prism_esp | 2.00 |
| prism | 3.00 |
| prism_fixed | 4.33 |
| s3_weighted | 5.33 |
| learned_scorer | 6.00 |
| rdshape_ensemble | 7.00 |
| shape_combo_rdkit | 8.33 |
| equal_weight | 8.67 |
| pharm2d | 8.67 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (7.82 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

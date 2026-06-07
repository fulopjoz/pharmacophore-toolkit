# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 46, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.975 [0.904,1.000] | 0.983 | 7.68 | +0.219 [+0.043,+0.466] |
| prism_esp | 0.973 [0.904,0.999] | 0.992 | 7.68 | +0.221 [+0.045,+0.464] |
| prism | 0.973 [0.902,0.999] | 0.991 | 7.68 | +0.221 [+0.045,+0.464] |
| prism_fixed | 0.927 [0.784,0.997] | 0.983 | 7.68 | +0.179 [-0.039,+0.434] |
| learned_scorer | 0.921 [0.780,0.989] | 0.953 | 7.68 | +0.166 [-0.024,+0.411] |
| rdshape_ensemble | 0.873 [0.661,0.961] | 0.868 | 7.68 | +0.116 [-0.125,+0.381] |
| pharm2d | 0.825 [0.594,0.944] | 0.889 | 7.68 | +0.073 [-0.212,+0.346] |
| shape_combo_rdkit | 0.779 [0.503,0.922] | 0.738 | 7.68 | +0.026 [-0.272,+0.301] |
| s3_weighted | 0.757 [0.492,0.938] | 0.914 | 3.84 | — (baseline) |
| equal_weight | 0.151 [0.001,0.336] | 0.362 | 3.84 | -0.621 [-0.882,-0.318] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.815 [0.646,0.935] | 0.936 | 30.35 | +0.550 [+0.334,+0.749] |
| prism_esp | 0.405 [0.184,0.624] | 0.682 | 15.18 | +0.142 [-0.160,+0.438] |
| prism | 0.347 [0.153,0.547] | 0.724 | 5.06 | +0.081 [-0.208,+0.377] |
| s3_weighted | 0.270 [0.089,0.446] | 0.676 | 5.06 | — (baseline) |
| prism_fixed | 0.166 [0.022,0.328] | 0.681 | 5.06 | -0.102 [-0.289,+0.107] |
| learned_scorer | 0.122 [0.017,0.263] | 0.615 | 5.06 | -0.143 [-0.358,+0.084] |
| pharm2d | 0.076 [0.000,0.230] | 0.221 | 5.06 | -0.196 [-0.427,+0.051] |
| rdshape_ensemble | 0.057 [0.000,0.175] | 0.415 | 0.00 | -0.213 [-0.428,+0.002] |
| shape_combo_rdkit | 0.015 [0.000,0.042] | 0.372 | 0.00 | -0.247 [-0.439,-0.070] |
| equal_weight | 0.009 [0.001,0.035] | 0.363 | 0.00 | -0.248 [-0.436,-0.072] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.987 [0.978,0.992] | 0.946 | 3.69 | +0.555 [+0.465,+0.653] |
| prism_esp | 0.950 [0.923,0.970] | 0.895 | 3.69 | +0.517 [+0.424,+0.614] |
| prism | 0.943 [0.914,0.965] | 0.897 | 3.69 | +0.510 [+0.416,+0.606] |
| prism_fixed | 0.915 [0.873,0.947] | 0.853 | 3.69 | +0.483 [+0.384,+0.577] |
| learned_scorer | 0.518 [0.417,0.616] | 0.659 | 3.36 | +0.087 [-0.045,+0.221] |
| s3_weighted | 0.432 [0.333,0.525] | 0.598 | 2.68 | — (baseline) |
| rdshape_ensemble | 0.378 [0.274,0.479] | 0.482 | 3.36 | -0.052 [-0.195,+0.085] |
| shape_combo_rdkit | 0.365 [0.267,0.465] | 0.489 | 2.68 | -0.067 [-0.203,+0.074] |
| pharm2d | 0.231 [0.141,0.326] | 0.439 | 2.01 | -0.199 [-0.334,-0.067] |
| equal_weight | 0.153 [0.103,0.242] | 0.462 | 0.00 | -0.264 [-0.381,-0.144] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=24.891, p=0.003; Nemenyi CD=7.822 over N=3 datasets, k=10 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 1.00 |
| prism_esp | 2.00 |
| prism | 3.00 |
| prism_fixed | 4.33 |
| learned_scorer | 5.33 |
| s3_weighted | 6.33 |
| rdshape_ensemble | 7.00 |
| pharm2d | 7.67 |
| shape_combo_rdkit | 8.33 |
| equal_weight | 10.00 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (7.82 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

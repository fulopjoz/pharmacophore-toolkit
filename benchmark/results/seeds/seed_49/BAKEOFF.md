# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 49, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| prism_esp | 0.976 [0.913,1.000] | 0.994 | 7.68 | +0.108 [+0.018,+0.313] |
| differential_mmfp | 0.971 [0.894,1.000] | 0.978 | 7.68 | +0.103 [+0.024,+0.295] |
| prism | 0.970 [0.897,1.000] | 0.993 | 7.68 | +0.103 [+0.013,+0.300] |
| prism_fixed | 0.886 [0.690,1.000] | 0.983 | 7.68 | +0.044 [-0.196,+0.262] |
| s3_weighted | 0.866 [0.660,0.959] | 0.919 | 7.68 | — (baseline) |
| rdshape_ensemble | 0.860 [0.650,0.957] | 0.880 | 7.68 | -0.002 [-0.205,+0.191] |
| learned_scorer | 0.858 [0.641,0.992] | 0.948 | 3.84 | +0.014 [-0.236,+0.234] |
| pharm2d | 0.846 [0.644,0.953] | 0.894 | 7.68 | -0.013 [-0.214,+0.207] |
| shape_combo_rdkit | 0.692 [0.422,0.875] | 0.746 | 7.68 | -0.160 [-0.421,+0.058] |
| equal_weight | 0.162 [0.002,0.346] | 0.394 | 3.84 | -0.747 [-0.932,-0.457] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.915 [0.809,0.995] | 0.991 | 30.35 | +0.685 [+0.469,+0.872] |
| prism_esp | 0.268 [0.097,0.473] | 0.751 | 10.12 | +0.034 [-0.228,+0.309] |
| s3_weighted | 0.234 [0.079,0.406] | 0.691 | 5.06 | — (baseline) |
| prism | 0.232 [0.073,0.413] | 0.731 | 10.12 | -0.004 [-0.221,+0.232] |
| prism_fixed | 0.126 [0.026,0.249] | 0.633 | 0.00 | -0.107 [-0.272,+0.045] |
| learned_scorer | 0.064 [0.006,0.159] | 0.632 | 0.00 | -0.172 [-0.343,-0.002] |
| equal_weight | 0.014 [0.001,0.066] | 0.360 | 0.00 | -0.214 [-0.399,-0.043] |
| shape_combo_rdkit | 0.010 [0.000,0.032] | 0.371 | 0.00 | -0.222 [-0.401,-0.065] |
| rdshape_ensemble | 0.001 [0.000,0.004] | 0.390 | 0.00 | -0.231 [-0.405,-0.078] |
| pharm2d | 0.000 [0.000,0.000] | 0.154 | 0.00 | -0.233 [-0.406,-0.079] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.989 [0.981,0.993] | 0.947 | 3.69 | +0.559 [+0.463,+0.653] |
| prism_esp | 0.962 [0.940,0.977] | 0.904 | 3.69 | +0.532 [+0.435,+0.627] |
| prism | 0.956 [0.932,0.974] | 0.905 | 3.69 | +0.526 [+0.430,+0.623] |
| prism_fixed | 0.931 [0.894,0.959] | 0.865 | 3.69 | +0.502 [+0.401,+0.599] |
| learned_scorer | 0.523 [0.425,0.618] | 0.664 | 3.36 | +0.093 [-0.035,+0.222] |
| s3_weighted | 0.429 [0.335,0.527] | 0.586 | 2.35 | — (baseline) |
| rdshape_ensemble | 0.376 [0.269,0.468] | 0.479 | 3.36 | -0.059 [-0.190,+0.086] |
| shape_combo_rdkit | 0.358 [0.255,0.453] | 0.481 | 2.68 | -0.074 [-0.208,+0.065] |
| pharm2d | 0.227 [0.140,0.313] | 0.449 | 2.01 | -0.205 [-0.334,-0.070] |
| equal_weight | 0.154 [0.109,0.245] | 0.458 | 0.00 | -0.260 [-0.376,-0.135] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=24.527, p=0.004; Nemenyi CD=7.822 over N=3 datasets, k=10 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 1.33 |
| prism_esp | 1.67 |
| prism | 3.33 |
| prism_fixed | 4.33 |
| s3_weighted | 4.67 |
| learned_scorer | 6.00 |
| rdshape_ensemble | 7.33 |
| shape_combo_rdkit | 8.33 |
| equal_weight | 9.00 |
| pharm2d | 9.00 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (7.82 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

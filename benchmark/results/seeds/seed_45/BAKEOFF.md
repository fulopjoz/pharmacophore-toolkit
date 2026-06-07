# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 45, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| prism | 0.987 [0.948,1.000] | 0.997 | 7.68 | +0.307 [+0.089,+0.585] |
| prism_esp | 0.986 [0.946,1.000] | 0.996 | 7.68 | +0.306 [+0.090,+0.584] |
| differential_mmfp | 0.974 [0.900,1.000] | 0.981 | 7.68 | +0.293 [+0.084,+0.552] |
| prism_fixed | 0.931 [0.799,0.995] | 0.979 | 7.68 | +0.253 [+0.020,+0.529] |
| learned_scorer | 0.899 [0.719,1.000] | 0.951 | 7.68 | +0.217 [-0.049,+0.499] |
| rdshape_ensemble | 0.880 [0.684,0.966] | 0.892 | 7.68 | +0.194 [-0.059,+0.474] |
| pharm2d | 0.802 [0.547,0.927] | 0.882 | 7.68 | +0.118 [-0.171,+0.398] |
| shape_combo_rdkit | 0.746 [0.471,0.905] | 0.759 | 7.68 | +0.063 [-0.246,+0.361] |
| s3_weighted | 0.683 [0.395,0.901] | 0.895 | 3.84 | — (baseline) |
| equal_weight | 0.155 [0.005,0.319] | 0.436 | 3.84 | -0.555 [-0.815,-0.250] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.863 [0.709,0.976] | 0.970 | 35.41 | +0.661 [+0.317,+0.917] |
| prism_esp | 0.534 [0.331,0.725] | 0.788 | 15.18 | +0.330 [-0.028,+0.615] |
| prism | 0.397 [0.178,0.599] | 0.751 | 15.18 | +0.187 [-0.179,+0.490] |
| s3_weighted | 0.218 [0.029,0.415] | 0.694 | 10.12 | — (baseline) |
| prism_fixed | 0.137 [0.045,0.264] | 0.706 | 0.00 | -0.071 [-0.323,+0.142] |
| learned_scorer | 0.126 [0.020,0.259] | 0.652 | 0.00 | -0.089 [-0.321,+0.136] |
| shape_combo_rdkit | 0.009 [0.000,0.027] | 0.331 | 0.00 | -0.202 [-0.410,-0.018] |
| rdshape_ensemble | 0.007 [0.000,0.023] | 0.353 | 0.00 | -0.202 [-0.408,-0.023] |
| equal_weight | 0.002 [0.000,0.008] | 0.279 | 0.00 | -0.208 [-0.415,-0.025] |
| pharm2d | 0.000 [0.000,0.000] | 0.157 | 0.00 | -0.210 [-0.415,-0.029] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.988 [0.979,0.993] | 0.946 | 3.69 | +0.549 [+0.454,+0.647] |
| prism_esp | 0.952 [0.923,0.970] | 0.896 | 3.69 | +0.513 [+0.415,+0.612] |
| prism | 0.946 [0.916,0.967] | 0.898 | 3.69 | +0.508 [+0.408,+0.607] |
| prism_fixed | 0.917 [0.873,0.948] | 0.854 | 3.69 | +0.479 [+0.376,+0.577] |
| learned_scorer | 0.527 [0.422,0.620] | 0.660 | 3.36 | +0.087 [-0.049,+0.223] |
| s3_weighted | 0.436 [0.339,0.534] | 0.598 | 2.01 | — (baseline) |
| rdshape_ensemble | 0.377 [0.268,0.476] | 0.480 | 3.36 | -0.064 [-0.209,+0.074] |
| shape_combo_rdkit | 0.364 [0.258,0.464] | 0.485 | 2.68 | -0.076 [-0.216,+0.064] |
| pharm2d | 0.231 [0.141,0.319] | 0.439 | 2.01 | -0.211 [-0.338,-0.075] |
| equal_weight | 0.158 [0.107,0.245] | 0.459 | 0.00 | -0.266 [-0.382,-0.145] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=24.018, p=0.004; Nemenyi CD=7.822 over N=3 datasets, k=10 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 1.67 |
| prism_esp | 2.00 |
| prism | 2.33 |
| prism_fixed | 4.33 |
| learned_scorer | 5.33 |
| s3_weighted | 6.33 |
| rdshape_ensemble | 7.00 |
| shape_combo_rdkit | 7.67 |
| pharm2d | 8.67 |
| equal_weight | 9.67 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (7.82 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

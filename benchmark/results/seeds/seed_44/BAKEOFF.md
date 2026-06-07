# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 44, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| prism_esp | 0.988 [0.954,1.000] | 0.996 | 7.68 | +0.207 [+0.069,+0.445] |
| prism | 0.988 [0.950,1.000] | 0.996 | 7.68 | +0.206 [+0.068,+0.445] |
| prism_fixed | 0.981 [0.926,1.000] | 0.991 | 7.68 | +0.199 [+0.062,+0.430] |
| differential_mmfp | 0.976 [0.908,1.000] | 0.982 | 7.68 | +0.195 [+0.060,+0.416] |
| learned_scorer | 0.955 [0.865,0.995] | 0.951 | 7.68 | +0.176 [+0.034,+0.404] |
| rdshape_ensemble | 0.861 [0.673,0.956] | 0.887 | 7.68 | +0.082 [-0.129,+0.324] |
| s3_weighted | 0.800 [0.524,0.926] | 0.903 | 7.68 | — (baseline) |
| pharm2d | 0.793 [0.550,0.926] | 0.894 | 7.68 | +0.013 [-0.228,+0.257] |
| shape_combo_rdkit | 0.669 [0.405,0.873] | 0.753 | 7.68 | -0.101 [-0.384,+0.175] |
| equal_weight | 0.164 [0.004,0.342] | 0.432 | 3.84 | -0.645 [-0.874,-0.334] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.825 [0.655,0.942] | 0.982 | 35.41 | +0.573 [+0.334,+0.768] |
| prism | 0.309 [0.129,0.493] | 0.655 | 0.00 | +0.052 [-0.180,+0.299] |
| prism_esp | 0.300 [0.116,0.480] | 0.618 | 5.06 | +0.042 [-0.201,+0.288] |
| s3_weighted | 0.257 [0.088,0.435] | 0.679 | 5.06 | — (baseline) |
| prism_fixed | 0.175 [0.033,0.359] | 0.722 | 10.12 | -0.079 [-0.267,+0.107] |
| learned_scorer | 0.143 [0.043,0.267] | 0.689 | 0.00 | -0.108 [-0.336,+0.113] |
| pharm2d | 0.076 [0.000,0.242] | 0.252 | 5.06 | -0.186 [-0.409,+0.061] |
| rdshape_ensemble | 0.055 [0.000,0.181] | 0.373 | 0.00 | -0.202 [-0.410,+0.012] |
| equal_weight | 0.009 [0.000,0.019] | 0.337 | 0.00 | -0.245 [-0.432,-0.080] |
| shape_combo_rdkit | 0.008 [0.000,0.030] | 0.321 | 0.00 | -0.243 [-0.432,-0.078] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.987 [0.978,0.992] | 0.945 | 3.69 | +0.530 [+0.427,+0.621] |
| prism_esp | 0.953 [0.925,0.971] | 0.897 | 3.69 | +0.495 [+0.394,+0.591] |
| prism | 0.947 [0.917,0.967] | 0.900 | 3.69 | +0.489 [+0.388,+0.585] |
| prism_fixed | 0.917 [0.872,0.948] | 0.855 | 3.69 | +0.458 [+0.356,+0.559] |
| learned_scorer | 0.518 [0.411,0.610] | 0.654 | 3.36 | +0.060 [-0.079,+0.197] |
| s3_weighted | 0.459 [0.364,0.562] | 0.604 | 2.35 | — (baseline) |
| rdshape_ensemble | 0.381 [0.275,0.479] | 0.480 | 3.36 | -0.077 [-0.221,+0.067] |
| shape_combo_rdkit | 0.368 [0.267,0.470] | 0.484 | 2.68 | -0.090 [-0.234,+0.048] |
| pharm2d | 0.230 [0.147,0.326] | 0.439 | 2.01 | -0.227 [-0.358,-0.089] |
| equal_weight | 0.168 [0.104,0.241] | 0.456 | 0.34 | -0.289 [-0.416,-0.168] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=24.527, p=0.004; Nemenyi CD=7.822 over N=3 datasets, k=10 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 2.00 |
| prism_esp | 2.00 |
| prism | 2.33 |
| prism_fixed | 4.00 |
| learned_scorer | 5.33 |
| s3_weighted | 5.67 |
| rdshape_ensemble | 7.00 |
| pharm2d | 8.00 |
| shape_combo_rdkit | 9.00 |
| equal_weight | 9.67 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (7.82 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

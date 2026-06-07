# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 48, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| prism | 0.991 [0.961,1.000] | 0.997 | 7.68 | +0.218 [+0.039,+0.468] |
| prism_esp | 0.990 [0.956,1.000] | 0.996 | 7.68 | +0.217 [+0.038,+0.467] |
| prism_fixed | 0.984 [0.939,1.000] | 0.994 | 7.68 | +0.211 [+0.035,+0.464] |
| learned_scorer | 0.983 [0.934,1.000] | 0.959 | 7.68 | +0.210 [+0.029,+0.461] |
| differential_mmfp | 0.974 [0.904,1.000] | 0.980 | 7.68 | +0.200 [+0.026,+0.448] |
| rdshape_ensemble | 0.883 [0.685,0.968] | 0.903 | 7.68 | +0.105 [-0.139,+0.383] |
| pharm2d | 0.787 [0.497,0.922] | 0.881 | 7.68 | +0.009 [-0.277,+0.286] |
| s3_weighted | 0.768 [0.518,0.956] | 0.922 | 3.84 | — (baseline) |
| shape_combo_rdkit | 0.741 [0.465,0.911] | 0.760 | 7.68 | -0.034 [-0.334,+0.270] |
| equal_weight | 0.287 [0.005,0.457] | 0.451 | 3.84 | -0.577 [-0.881,-0.217] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.685 [0.473,0.860] | 0.911 | 30.35 | +0.567 [+0.327,+0.773] |
| prism_esp | 0.338 [0.171,0.522] | 0.745 | 0.00 | +0.219 [-0.011,+0.446] |
| prism | 0.198 [0.071,0.352] | 0.725 | 0.00 | +0.077 [-0.120,+0.257] |
| learned_scorer | 0.188 [0.044,0.351] | 0.638 | 5.06 | +0.059 [-0.115,+0.253] |
| s3_weighted | 0.125 [0.020,0.281] | 0.612 | 5.06 | — (baseline) |
| pharm2d | 0.076 [0.000,0.233] | 0.238 | 5.06 | -0.050 [-0.259,+0.166] |
| rdshape_ensemble | 0.057 [0.000,0.184] | 0.421 | 0.00 | -0.063 [-0.265,+0.112] |
| prism_fixed | 0.043 [0.005,0.099] | 0.615 | 0.00 | -0.075 [-0.209,+0.013] |
| equal_weight | 0.028 [0.002,0.090] | 0.377 | 0.00 | -0.084 [-0.268,+0.035] |
| shape_combo_rdkit | 0.013 [0.000,0.036] | 0.398 | 0.00 | -0.104 [-0.271,+0.001] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.987 [0.979,0.993] | 0.946 | 3.69 | +0.567 [+0.468,+0.661] |
| prism_esp | 0.952 [0.923,0.970] | 0.898 | 3.69 | +0.529 [+0.429,+0.628] |
| prism | 0.945 [0.912,0.966] | 0.899 | 3.69 | +0.522 [+0.423,+0.620] |
| prism_fixed | 0.917 [0.874,0.948] | 0.857 | 3.69 | +0.494 [+0.390,+0.596] |
| learned_scorer | 0.541 [0.440,0.640] | 0.657 | 3.36 | +0.121 [-0.018,+0.248] |
| s3_weighted | 0.423 [0.322,0.521] | 0.599 | 2.35 | — (baseline) |
| rdshape_ensemble | 0.377 [0.273,0.483] | 0.487 | 3.36 | -0.046 [-0.185,+0.096] |
| shape_combo_rdkit | 0.365 [0.260,0.467] | 0.492 | 2.68 | -0.058 [-0.200,+0.081] |
| pharm2d | 0.226 [0.144,0.320] | 0.446 | 2.01 | -0.195 [-0.331,-0.058] |
| equal_weight | 0.171 [0.107,0.239] | 0.457 | 0.34 | -0.254 [-0.378,-0.128] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=22.564, p=0.007; Nemenyi CD=7.822 over N=3 datasets, k=10 methods.

| method | avg rank |
|--------|----------|
| prism_esp | 2.00 |
| differential_mmfp | 2.33 |
| prism | 2.33 |
| learned_scorer | 4.33 |
| prism_fixed | 5.00 |
| s3_weighted | 6.33 |
| rdshape_ensemble | 6.67 |
| pharm2d | 7.33 |
| shape_combo_rdkit | 9.00 |
| equal_weight | 9.67 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (7.82 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 51, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| prism_esp | 0.992 [0.964,1.000] | 0.998 | 7.68 | +0.293 [+0.104,+0.554] |
| prism | 0.991 [0.959,1.000] | 0.997 | 7.68 | +0.291 [+0.102,+0.552] |
| prism_fixed | 0.979 [0.919,1.000] | 0.989 | 7.68 | +0.279 [+0.093,+0.536] |
| differential_mmfp | 0.973 [0.903,1.000] | 0.985 | 7.68 | +0.271 [+0.085,+0.533] |
| learned_scorer | 0.933 [0.813,0.986] | 0.944 | 7.68 | +0.230 [+0.028,+0.490] |
| rdshape_ensemble | 0.880 [0.684,0.964] | 0.876 | 7.68 | +0.178 [-0.058,+0.438] |
| pharm2d | 0.846 [0.637,0.949] | 0.897 | 7.68 | +0.145 [-0.101,+0.408] |
| shape_combo_rdkit | 0.750 [0.478,0.905] | 0.741 | 7.68 | +0.051 [-0.250,+0.327] |
| s3_weighted | 0.734 [0.436,0.891] | 0.897 | 7.68 | — (baseline) |
| equal_weight | 0.074 [0.001,0.195] | 0.372 | 0.00 | -0.628 [-0.846,-0.349] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.825 [0.658,0.949] | 0.976 | 30.35 | +0.736 [+0.530,+0.893] |
| prism_fixed | 0.239 [0.096,0.395] | 0.728 | 0.00 | +0.144 [-0.001,+0.301] |
| prism | 0.201 [0.071,0.354] | 0.650 | 0.00 | +0.110 [-0.081,+0.302] |
| prism_esp | 0.194 [0.067,0.340] | 0.641 | 0.00 | +0.100 [-0.083,+0.276] |
| learned_scorer | 0.167 [0.032,0.318] | 0.710 | 5.06 | +0.070 [-0.110,+0.261] |
| s3_weighted | 0.093 [0.013,0.236] | 0.703 | 5.06 | — (baseline) |
| rdshape_ensemble | 0.006 [0.000,0.020] | 0.403 | 0.00 | -0.082 [-0.231,-0.003] |
| shape_combo_rdkit | 0.002 [0.000,0.007] | 0.354 | 0.00 | -0.086 [-0.234,-0.010] |
| equal_weight | 0.000 [0.000,0.001] | 0.257 | 0.00 | -0.088 [-0.236,-0.012] |
| pharm2d | 0.000 [0.000,0.000] | 0.170 | 0.00 | -0.089 [-0.236,-0.013] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.987 [0.978,0.993] | 0.946 | 3.69 | +0.579 [+0.481,+0.674] |
| prism_esp | 0.956 [0.931,0.973] | 0.900 | 3.69 | +0.548 [+0.451,+0.642] |
| prism | 0.950 [0.923,0.969] | 0.901 | 3.69 | +0.542 [+0.444,+0.636] |
| prism_fixed | 0.924 [0.885,0.953] | 0.859 | 3.69 | +0.515 [+0.414,+0.614] |
| learned_scorer | 0.515 [0.412,0.613] | 0.653 | 3.36 | +0.107 [-0.026,+0.247] |
| s3_weighted | 0.409 [0.311,0.510] | 0.580 | 2.35 | — (baseline) |
| rdshape_ensemble | 0.380 [0.278,0.480] | 0.489 | 3.36 | -0.031 [-0.168,+0.111] |
| shape_combo_rdkit | 0.371 [0.271,0.469] | 0.494 | 2.68 | -0.040 [-0.173,+0.102] |
| pharm2d | 0.226 [0.140,0.315] | 0.442 | 2.01 | -0.183 [-0.310,-0.054] |
| equal_weight | 0.157 [0.106,0.243] | 0.465 | 0.00 | -0.239 [-0.354,-0.113] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=24.236, p=0.004; Nemenyi CD=7.822 over N=3 datasets, k=10 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 2.00 |
| prism_esp | 2.33 |
| prism | 2.67 |
| prism_fixed | 3.00 |
| learned_scorer | 5.00 |
| rdshape_ensemble | 6.67 |
| s3_weighted | 7.00 |
| shape_combo_rdkit | 8.00 |
| pharm2d | 8.67 |
| equal_weight | 9.67 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (7.82 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

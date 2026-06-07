# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 42, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, s3_3d.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| s3_3d | 0.972 [0.902,1.000] | 0.993 | 7.68 | +0.184 [+0.017,+0.430] |
| differential_mmfp | 0.971 [0.895,1.000] | 0.981 | 7.68 | +0.180 [+0.021,+0.421] |
| rdshape_ensemble | 0.869 [0.672,0.963] | 0.871 | 7.68 | +0.080 [-0.158,+0.334] |
| pharm2d | 0.810 [0.569,0.933] | 0.893 | 7.68 | +0.021 [-0.253,+0.281] |
| learned_scorer | 0.789 [0.558,0.967] | 0.935 | 3.84 | +0.016 [-0.287,+0.306] |
| s3_weighted | 0.781 [0.532,0.961] | 0.935 | 3.84 | — (baseline) |
| shape_combo_rdkit | 0.729 [0.466,0.905] | 0.738 | 7.68 | -0.058 [-0.352,+0.236] |
| equal_weight | 0.211 [0.002,0.416] | 0.400 | 3.84 | -0.606 [-0.883,-0.268] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.705 [0.512,0.868] | 0.941 | 35.41 | +0.538 [+0.268,+0.768] |
| s3_3d | 0.293 [0.114,0.492] | 0.645 | 10.12 | +0.114 [-0.129,+0.353] |
| s3_weighted | 0.175 [0.032,0.351] | 0.726 | 10.12 | — (baseline) |
| learned_scorer | 0.133 [0.007,0.286] | 0.651 | 0.00 | -0.040 [-0.232,+0.159] |
| rdshape_ensemble | 0.006 [0.000,0.019] | 0.364 | 0.00 | -0.164 [-0.348,-0.025] |
| equal_weight | 0.002 [0.000,0.006] | 0.252 | 0.00 | -0.168 [-0.350,-0.030] |
| shape_combo_rdkit | 0.002 [0.000,0.006] | 0.340 | 0.00 | -0.169 [-0.351,-0.030] |
| pharm2d | 0.000 [0.000,0.000] | 0.144 | 0.00 | -0.171 [-0.351,-0.032] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.987 [0.979,0.993] | 0.946 | 3.69 | +0.573 [+0.481,+0.658] |
| s3_3d | 0.952 [0.924,0.970] | 0.901 | 3.69 | +0.537 [+0.443,+0.624] |
| learned_scorer | 0.522 [0.422,0.621] | 0.653 | 3.36 | +0.106 [-0.031,+0.236] |
| s3_weighted | 0.417 [0.325,0.509] | 0.579 | 2.01 | — (baseline) |
| rdshape_ensemble | 0.378 [0.274,0.481] | 0.488 | 3.36 | -0.037 [-0.175,+0.102] |
| shape_combo_rdkit | 0.365 [0.267,0.462] | 0.494 | 2.68 | -0.049 [-0.181,+0.085] |
| pharm2d | 0.230 [0.144,0.321] | 0.446 | 2.01 | -0.188 [-0.310,-0.057] |
| equal_weight | 0.163 [0.105,0.237] | 0.463 | 0.34 | -0.249 [-0.363,-0.134] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=17.222, p=0.016; Nemenyi CD=6.062 over N=3 datasets, k=8 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 1.33 |
| s3_3d | 1.67 |
| learned_scorer | 4.00 |
| s3_weighted | 4.33 |
| rdshape_ensemble | 4.33 |
| pharm2d | 6.33 |
| shape_combo_rdkit | 6.67 |
| equal_weight | 7.33 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (6.06 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

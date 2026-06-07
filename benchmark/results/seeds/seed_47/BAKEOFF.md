# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 47, test_frac 0.25, bootstrap 2000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d, shape_combo_rdkit, rdshape_ensemble, learned_scorer, prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| prism_fixed | 0.979 [0.913,1.000] | 0.989 | 7.68 | +0.135 [+0.042,+0.335] |
| prism_esp | 0.976 [0.909,1.000] | 0.990 | 7.68 | +0.132 [+0.040,+0.333] |
| prism | 0.976 [0.911,1.000] | 0.990 | 7.68 | +0.132 [+0.039,+0.332] |
| differential_mmfp | 0.972 [0.888,1.000] | 0.976 | 7.68 | +0.127 [+0.032,+0.327] |
| learned_scorer | 0.956 [0.863,0.994] | 0.954 | 7.68 | +0.113 [+0.017,+0.315] |
| rdshape_ensemble | 0.861 [0.659,0.957] | 0.875 | 7.68 | +0.022 [-0.174,+0.239] |
| s3_weighted | 0.852 [0.618,0.948] | 0.910 | 7.68 | — (baseline) |
| pharm2d | 0.784 [0.533,0.923] | 0.882 | 7.68 | -0.050 [-0.303,+0.192] |
| shape_combo_rdkit | 0.712 [0.452,0.894] | 0.743 | 7.68 | -0.115 [-0.382,+0.138] |
| equal_weight | 0.184 [0.005,0.331] | 0.394 | 3.84 | -0.702 [-0.904,-0.423] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.869 [0.718,0.975] | 0.982 | 30.35 | +0.646 [+0.411,+0.856] |
| prism_esp | 0.248 [0.096,0.419] | 0.738 | 0.00 | +0.024 [-0.182,+0.245] |
| prism | 0.226 [0.085,0.390] | 0.730 | 0.00 | +0.001 [-0.192,+0.218] |
| s3_weighted | 0.222 [0.061,0.399] | 0.759 | 5.06 | — (baseline) |
| prism_fixed | 0.172 [0.019,0.342] | 0.684 | 0.00 | -0.055 [-0.239,+0.125] |
| learned_scorer | 0.090 [0.007,0.216] | 0.604 | 0.00 | -0.133 [-0.330,+0.068] |
| equal_weight | 0.019 [0.000,0.065] | 0.298 | 0.00 | -0.204 [-0.386,-0.030] |
| rdshape_ensemble | 0.006 [0.000,0.022] | 0.403 | 0.00 | -0.213 [-0.395,-0.053] |
| shape_combo_rdkit | 0.003 [0.000,0.008] | 0.359 | 0.00 | -0.217 [-0.397,-0.056] |
| pharm2d | 0.000 [0.000,0.000] | 0.137 | 0.00 | -0.220 [-0.399,-0.061] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.988 [0.980,0.993] | 0.947 | 3.69 | +0.569 [+0.473,+0.665] |
| prism_esp | 0.951 [0.923,0.970] | 0.897 | 3.69 | +0.531 [+0.435,+0.631] |
| prism | 0.945 [0.913,0.966] | 0.898 | 3.69 | +0.524 [+0.428,+0.623] |
| prism_fixed | 0.919 [0.878,0.949] | 0.857 | 3.69 | +0.499 [+0.396,+0.603] |
| learned_scorer | 0.542 [0.433,0.634] | 0.660 | 3.36 | +0.121 [-0.015,+0.254] |
| s3_weighted | 0.419 [0.321,0.515] | 0.595 | 2.35 | — (baseline) |
| rdshape_ensemble | 0.381 [0.275,0.481] | 0.491 | 3.36 | -0.039 [-0.177,+0.106] |
| shape_combo_rdkit | 0.369 [0.265,0.471] | 0.495 | 2.68 | -0.050 [-0.189,+0.095] |
| pharm2d | 0.229 [0.143,0.316] | 0.450 | 2.01 | -0.191 [-0.321,-0.068] |
| equal_weight | 0.163 [0.103,0.237] | 0.464 | 0.00 | -0.256 [-0.371,-0.125] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=23.655, p=0.005; Nemenyi CD=7.822 over N=3 datasets, k=10 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 2.00 |
| prism_esp | 2.00 |
| prism | 3.00 |
| prism_fixed | 3.33 |
| learned_scorer | 5.33 |
| s3_weighted | 5.67 |
| rdshape_ensemble | 7.00 |
| shape_combo_rdkit | 8.67 |
| equal_weight | 9.00 |
| pharm2d | 9.00 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (7.82 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

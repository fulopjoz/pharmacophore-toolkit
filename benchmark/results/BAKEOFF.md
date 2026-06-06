# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 42, test_frac 0.25, bootstrap 1000, methods: equal_weight, s3_weighted, differential_mmfp, pharm2d.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.973 [0.902,1.000] | 0.985 | 7.68 | +0.309 [+0.092,+0.587] |
| pharm2d | 0.825 [0.600,0.937] | 0.887 | 7.68 | +0.161 [-0.103,+0.449] |
| s3_weighted | 0.678 [0.366,0.888] | 0.906 | 3.84 | — (baseline) |
| equal_weight | 0.167 [0.002,0.354] | 0.397 | 3.84 | -0.516 [-0.809,-0.165] |

## ccr2_mubd  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.776 [0.609,0.920] | 0.942 | 30.35 | +0.563 [+0.321,+0.770] |
| s3_weighted | 0.215 [0.049,0.399] | 0.559 | 10.12 | — (baseline) |
| pharm2d | 0.076 [0.000,0.248] | 0.202 | 5.06 | -0.146 [-0.375,+0.100] |
| equal_weight | 0.017 [0.003,0.072] | 0.474 | 0.00 | -0.184 [-0.381,-0.007] |

## created_CCR2  ·  decoy-bias gate: **mild-bias**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| differential_mmfp | 0.988 [0.981,0.994] | 0.947 | 3.69 | +0.568 [+0.477,+0.661] |
| s3_weighted | 0.419 [0.324,0.515] | 0.593 | 2.35 | — (baseline) |
| pharm2d | 0.228 [0.139,0.316] | 0.453 | 2.01 | -0.195 [-0.323,-0.066] |
| equal_weight | 0.156 [0.112,0.251] | 0.463 | 0.67 | -0.244 [-0.363,-0.121] |

## Cross-dataset ranking (lower avg rank = better)

Friedman χ²=8.200, p=0.042; Nemenyi CD=2.708 over N=3 datasets, k=4 methods.

| method | avg rank |
|--------|----------|
| differential_mmfp | 1.00 |
| s3_weighted | 2.33 |
| pharm2d | 2.67 |
| equal_weight | 4.00 |

> **Caveat:** with only N=3 datasets the Nemenyi CD (2.71 rank units) is very wide — almost no pair will be 'significantly' separated. Treat the average-rank order as the headline and the per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test.

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

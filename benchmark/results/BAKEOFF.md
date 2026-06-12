# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)

Seed 42, test_frac 0.25, bootstrap 150, methods: prism, prism_fixed, prism_esp.  Primary metric **BEDROC(α=20)**; Δ vs **s3_weighted** with bootstrap 95% CI.

## ccr2_project  ·  decoy-bias gate: **unbiased**

| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |
|--------|-----------------|-----|------|------------------------|
| prism | 0.935 [0.799,0.984] | 0.969 | 7.68 | — (baseline) |
| prism_esp | 0.925 [0.812,0.985] | 0.971 | 7.68 | — (baseline) |
| prism_fixed | 0.798 [0.599,0.973] | 0.966 | 3.84 | — (baseline) |

## Cross-dataset ranking

_1 dataset_

## How to read this bake-off

- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization but not distributional separability, so a near-perfect BEDROC there signals what pure 2D ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "do fingerprints identify diverse actives? (no)" 10.3390/ph17080992). Read it as context, not as the method to ship.
- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` (2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), `shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` (supervised logistic on per-ref shape/color).
- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but collapses on MUBD is exploiting decoy bias.
- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results as relative-enrichment only.

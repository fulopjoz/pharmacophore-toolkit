# CCR2 Enrichment — Methods Comparison (held-out 5-fold CV)

Metric: BEDROC(α=20) primary; AUC/EF secondary. Δ vs our **s3_weighted** baseline with bootstrap 95% CI. A method beats S3 only if ΔBEDROC>0 & CI excludes 0.

| method | AUC | BEDROC | EF1% | ΔBEDROC vs S3 | CI | beats S3? |
|---|---|---|---|---|---|---|
| differential_mmfp | 0.9931 | 0.9878 | 7.7297 | 0.3806 | [+0.262,+0.507] | ✅ |
| fused_gw | 0.9245 | 0.8467 | 7.7297 | 0.2406 | [+0.106,+0.366] | ✅ |
| shape_combo_rdkit | 0.8055 | 0.6561 | 7.7297 | 0.047 | [-0.086,+0.183] | no |
| s3_weighted | 0.8725 | 0.604 | 3.8649 | — | — | baseline |
| equal_weight | 0.6398 | 0.175 | 1.2883 | -0.3985 | [-0.550,-0.233] | no |

# CCR2 Enrichment — Methods Comparison (held-out 5-fold CV)

Metric: BEDROC(α=20) primary; AUC/EF secondary. Δ vs our **s3_weighted** baseline with bootstrap 95% CI. A method beats S3 only if ΔBEDROC>0 & CI excludes 0.

| method | AUC | BEDROC | EF1% | ΔBEDROC vs S3 | CI | beats S3? |
|---|---|---|---|---|---|---|
| differential_mmfp | 0.9931 | 0.9878 | 7.7297 | 0.3734 | [+0.254,+0.499] | ✅ |
| fused_gw | 0.9397 | 0.8941 | 7.7297 | 0.279 | [+0.151,+0.406] | ✅ |
| shape_combo_rdkit | 0.8055 | 0.6561 | 7.7297 | 0.0388 | [-0.098,+0.179] | no |
| s3_weighted | 0.8742 | 0.6119 | 3.8649 | — | — | baseline |
| equal_weight | 0.5194 | 0.1933 | 0.0 | -0.4042 | [-0.553,-0.234] | no |

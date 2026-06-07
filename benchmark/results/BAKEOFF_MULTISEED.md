# Multi-seed bake-off — BEDROC across split seeds

10 scaffold-split seeds; BEDROC(alpha=20) mean +/- across-seed 95% CI of the mean. Across-seed variance = split-choice noise (the single-split bootstrap CI misses it). Does NOT power the Friedman test (correlated resamples).

## ccr2_project  ·  gate: **unbiased**  ·  seeds=10

| method | BEDROC mean [95% CI] | std | n_seeds |
|--------|----------------------|-----|---------|
| equal_weight | 0.181 [0.143,0.218] | 0.061 | 10 |
| s3_weighted | 0.771 [0.736,0.806] | 0.057 | 10 |
| differential_mmfp | 0.973 [0.972,0.974] | 0.002 | 10 |
| pharm2d | 0.806 [0.790,0.822] | 0.026 | 10 |
| shape_combo_rdkit | 0.720 [0.694,0.746] | 0.041 | 10 |
| rdshape_ensemble | 0.870 [0.863,0.877] | 0.011 | 10 |
| learned_scorer | 0.902 [0.863,0.940] | 0.062 | 10 |
| prism | 0.981 [0.976,0.986] | 0.008 | 10 |
| prism_fixed | 0.940 [0.914,0.967] | 0.043 | 10 |
| prism_esp | 0.981 [0.977,0.986] | 0.008 | 10 |

**Paired deltas across seeds:**
- `prism - prism_fixed`: +0.040 [+0.018,+0.063] (>0 in 90% of 10 seeds)
- `prism_esp - prism`: +0.000 [-0.001,+0.002] (>0 in 60% of 10 seeds)

## ccr2_mubd  ·  gate: **mild-bias**  ·  seeds=10

| method | BEDROC mean [95% CI] | std | n_seeds |
|--------|----------------------|-----|---------|
| equal_weight | 0.012 [0.005,0.019] | 0.012 | 10 |
| s3_weighted | 0.223 [0.176,0.271] | 0.076 | 10 |
| differential_mmfp | 0.817 [0.773,0.862] | 0.072 | 10 |
| pharm2d | 0.023 [0.000,0.045] | 0.037 | 10 |
| shape_combo_rdkit | 0.007 [0.004,0.010] | 0.005 | 10 |
| rdshape_ensemble | 0.020 [0.004,0.036] | 0.025 | 10 |
| learned_scorer | 0.126 [0.094,0.158] | 0.052 | 10 |
| prism | 0.257 [0.204,0.309] | 0.084 | 10 |
| prism_fixed | 0.161 [0.102,0.220] | 0.095 | 10 |
| prism_esp | 0.310 [0.242,0.378] | 0.110 | 10 |

**Paired deltas across seeds:**
- `prism - prism_fixed`: +0.096 [+0.026,+0.166] (>0 in 80% of 10 seeds)
- `prism_esp - prism`: +0.054 [+0.020,+0.087] (>0 in 80% of 10 seeds)

## created_CCR2  ·  gate: **mild-bias**  ·  seeds=10

| method | BEDROC mean [95% CI] | std | n_seeds |
|--------|----------------------|-----|---------|
| equal_weight | 0.160 [0.157,0.164] | 0.006 | 10 |
| s3_weighted | 0.429 [0.420,0.437] | 0.014 | 10 |
| differential_mmfp | 0.987 [0.987,0.988] | 0.001 | 10 |
| pharm2d | 0.229 [0.227,0.230] | 0.002 | 10 |
| shape_combo_rdkit | 0.366 [0.364,0.368] | 0.004 | 10 |
| rdshape_ensemble | 0.381 [0.377,0.384] | 0.006 | 10 |
| learned_scorer | 0.530 [0.518,0.541] | 0.018 | 10 |
| prism | 0.948 [0.945,0.950] | 0.004 | 10 |
| prism_fixed | 0.921 [0.917,0.924] | 0.006 | 10 |
| prism_esp | 0.954 [0.951,0.956] | 0.004 | 10 |

**Paired deltas across seeds:**
- `prism - prism_fixed`: +0.027 [+0.026,+0.029] (>0 in 100% of 10 seeds)
- `prism_esp - prism`: +0.006 [+0.005,+0.006] (>0 in 100% of 10 seeds)

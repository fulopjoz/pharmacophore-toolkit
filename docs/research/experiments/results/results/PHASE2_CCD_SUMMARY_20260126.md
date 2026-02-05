# Phase 2: CCD Response Surface Results

## Experiment Overview
- **Date**: 2026-01-26 17:59
- **Design**: Central Composite Design (2 factors)
- **Total Runs**: 14/14
- **Dataset**: CCR2

## Parameter Ranges
- **Tolerance**: 0.800 - 1.500 Å
- **Threshold**: 0.300 - 0.600
- **Linkage**: complete (fixed)

## Performance Summary
- **Best ROC-AUC**: 0.7234
- **Median ROC-AUC**: 0.5730
- **Worst ROC-AUC**: 0.4453

## Top 5 Configurations
| experiment_id   |   tolerance |   threshold |   roc_auc |   ef_1 |   n_features |
|:----------------|------------:|------------:|----------:|-------:|-------------:|
| CCD_011         |        1.5  |        0.3  |  0.723393 | 3.0973 |           13 |
| CCD_001         |        0.8  |        0.45 |  0.586199 | 0      |            1 |
| CCD_002         |        0.8  |        0.6  |  0.586199 | 0      |            1 |
| CCD_003         |        1.15 |        0.45 |  0.573011 | 0      |            3 |
| CCD_005         |        1.15 |        0.45 |  0.573011 | 0      |            3 |
# Conformer Sensitivity Sweep (CCR2)

Date: 2026-02-03

Setup: ablation_study.py with BO disabled; ETKDGv3 conformers; minimize flag only in +minimize config.
Dataset: tutorials/data/actives_ccr2_N75.csv (74 actives), tutorials/data/decoys_ccr2_N500.csv (498 decoys)

## AUC by n_conformers
| Config | n=5 | n=10 | n=25 |
|---|---:|---:|---:|
| RefEnsemble (max) | 0.8587 | 0.8944 | 0.9354 |
| Baseline | 0.3146 | 0.3017 | 0.3167 |
| +spatial | 0.3358 | 0.3071 | 0.2935 |
| +ensemble_clust | 0.6625 | 0.6583 | 0.6779 |
| +ot_scoring | 0.3468 | 0.3360 | 0.3314 |
| +minimize | 0.2155 | 0.2397 | 0.2839 |

## BEDROC (alpha=20) by n_conformers
| Config | n=5 | n=10 | n=25 |
|---|---:|---:|---:|
| RefEnsemble (max) | 0.8138 | 0.8417 | 0.8679 |
| Baseline | 0.0378 | 0.0368 | 0.0351 |
| +spatial | 0.0039 | 0.0035 | 0.0041 |
| +ensemble_clust | 0.2188 | 0.2065 | 0.2181 |
| +ot_scoring | 0.0020 | 0.0014 | 0.0013 |
| +minimize | 0.0021 | 0.0020 | 0.0114 |

## EF@1% by n_conformers
| Config | n=5 | n=10 | n=25 |
|---|---:|---:|---:|
| RefEnsemble (max) | 7.73 | 7.73 | 7.73 |
| Baseline | 1.55 | 0.00 | 0.00 |
| +spatial | 0.00 | 0.00 | 0.00 |
| +ensemble_clust | 0.00 | 0.00 | 0.00 |
| +ot_scoring | 0.00 | 0.00 | 0.00 |
| +minimize | 0.00 | 0.00 | 0.00 |

## Notes
- Removed 2 decoys that triggered UFFTYPER charge-state warnings during ETKDGv3 embedding:
  - Cc1cc(O)ccc1C1=C(c2ccc(O)cc2C)[C@H]2[C@@H](S(=O)(=O)Oc3ccccc3Cl)C[C@@H]1[S+]2[O-]
  - CCOc1ccc(CCNC(=O)C[S+]([O-])Cc2nc(-c3ccc(C)cc3)oc2C)cc1OCC

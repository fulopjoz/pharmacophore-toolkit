# CCR2 MUBD — Maximal Unbiased Benchmark (Xia et al. 2018)

## Source

| Field | Value |
|-------|-------|
| **Paper** | Xia J, Reid TE, Wu S, Zhang L, Wang XS. *Maximal Unbiased Benchmarking Data Sets for Human Chemokine Receptors and Comparative Analysis.* J Chem Inf Model. 2018;58(5):1104-1120. |
| **DOI** | [10.1021/acs.jcim.8b00004](https://doi.org/10.1021/acs.jcim.8b00004) |
| **Data repository** | https://github.com/jwxia2014/MUBD-hCRs |

## Dataset Statistics

| Set | Count | Source file in zip |
|-----|-------|--------------------|
| Actives (CCR2, ≤1 µM) | **60** | `ligand_set/uploaded_CCR2_1uM_new_ligands.sdf` |
| Max-unbiased decoys | **2340** | `decoy_set/uploaded_CCR2_1uM_final_decoys.sdf` |

Decoy ratio: 39:1.

## Decoy Generation Method

Decoys were generated using **MUBD-DecoyMaker 1.0** (Pipeline Pilot + MATLAB).
The algorithm:
1. Filters actives for pairwise Tanimoto similarity < 0.75 (MACCS fingerprints).
2. Matches decoys on six physicochemical properties: ALogP, MW, HBD, HBA, rotatable bonds, net charge.
3. Applies sphere-exclusion to maximise spatial coverage and minimise analogue bias.

Decoy type: `max-unbiased` — **no known analogue bias** (unlike DUD-E property-matched sets).

## Auto-Download

The SDF files are auto-downloaded from GitHub by `fetch.py` / `load.py`.

```bash
# Manual fetch
python fetch.py

# Or just import — fetch happens automatically
python -c "from benchmark.datasets.ccr2_mubd.load import load; a, d, m = load(); print(len(a), len(d))"
```

The archive URL is:
```
https://github.com/jwxia2014/MUBD-hCRs/raw/master/MUBD-hCRs.zip
```
Size: ~11 MB. Only the two CCR2 SDF files are extracted (~9 MB total for CCR2 decoys).

## Data Format

SDF files contain 2D structures in V2000 format (SciTegic/Pipeline Pilot).
`load.py` uses RDKit to parse molecules and return canonical SMILES.
The `Name` field in each SD record contains the ChEMBL ID.

## Citation

If you use this dataset, please cite:

> Xia J, Reid TE, Wu S, Zhang L, Wang XS. Maximal Unbiased Benchmarking Data Sets for Human Chemokine Receptors and Comparative Analysis. *J Chem Inf Model.* 2018;**58**(5):1104–1120. https://doi.org/10.1021/acs.jcim.8b00004

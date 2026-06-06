# LIT-PCBA Benchmark Dataset

## What is LIT-PCBA?

**LIT-PCBA** (Literature/PubChem Bioassay) is an unbiased virtual screening benchmark dataset
from the Rognan lab (Laboratoire d'Innovation Thérapeutique, Université de Strasbourg).

Unlike DUD-E and earlier benchmarks, LIT-PCBA uses **experimentally confirmed inactives**
from PubChem dose-response assays rather than property-matched synthetic decoys. Actives and
inactives are within similar molecular property ranges (balanced by the AVE procedure), making
it the most bias-resistant publicly available VS benchmark.

**Citation:** Tran-Nguyen V-K, Jacquemard C, Rognan D (2020). LIT-PCBA: An Unbiased Data Set
for Machine Learning and Virtual Screening. *J. Chem. Inf. Model.* **60**(9):4263-4273.
DOI: [10.1021/acs.jcim.0c00155](https://doi.org/10.1021/acs.jcim.0c00155)

## Access

- **Website:** http://drugdesign.unistra.fr/LIT-PCBA/
- **Full archive (~54 MB):** http://drugdesign.unistra.fr/LIT-PCBA/Files/full_data.tgz
  - Wayback Machine mirror: https://web.archive.org/web/20260411170149/https://drugdesign.unistra.fr/LIT-PCBA/Files/full_data.tgz
- **License:** Free for academic use (contact authors for commercial use)

## Decoy type

`experimental-inactive` — inactives are confirmed non-binders from the same PubChem
dose-response assays, tested at the same concentrations as the actives.

## Selected targets

| Target  | Full name                      | Type | PubChem AID | Actives | Inactives | Overlap w/ DUD-E |
|---------|-------------------------------|------|-------------|---------|-----------|------------------|
| ADRB2   | Beta-2 adrenergic receptor    | GPCR | 492947      | 17      | 312,483   | Yes (adrb2)      |
| OPRK1   | Kappa opioid receptor         | GPCR | 1777        | 24      | 269,816   | No               |

**Why these targets?**
- Both are GPCRs, matching the pharmacophore-toolkit use case (receptor-like binding)
- ADRB2 appears in both LIT-PCBA and DUD-E → enables direct bias-control comparison
- OPRK1 provides a second GPCR point with no DUD-E analogue (cleaner enrichment signal)

## All 15 LIT-PCBA targets

ADRB2, ALDH1, ESR1_ago, ESR1_ant, FEN1, GBA, IDH1, KAT2A, MAPK1, MTORC1,
OPRK1, PKM2, PPARG, TP53, VDR

Use `python fetch.py --all` to download all targets.

## File format

```
raw/<TARGET>/actives.smi    – SMILES<TAB>PubChem_CID
raw/<TARGET>/inactives.smi  – SMILES<TAB>PubChem_CID
```

## Usage

```bash
# Download data
python fetch.py                     # default: ADRB2 + OPRK1
python fetch.py --targets ADRB2     # specific target
python fetch.py --all               # all 15 targets

# Load into Python
python load.py                      # prints summary
```

```python
from benchmark.datasets.litpcba.load import load

actives, inactives, meta = load("ADRB2")
# actives:   list of SMILES strings (17 compounds)
# inactives: list of SMILES strings (312,483 compounds)
# meta:      dict with target, source, doi, decoy_type, n_act, n_dec, url
```

## Tests

```bash
python test_load.py
```

Verifies: non-empty lists, valid RDKit-parseable SMILES (first 50 checked), correct meta schema.

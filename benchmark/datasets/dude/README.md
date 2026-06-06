# DUD-E Benchmark Dataset

## What is DUD-E?

**DUD-E** (Directory of Useful Decoys, Enhanced) is a widely-used virtual screening benchmark
from the Irwin & Shoichet Laboratories (UCSF). It contains 102 protein targets with ~22,886
clustered ligands (from ChEMBL), each paired with ~50 property-matched decoys (from ZINC).

Decoys are selected to be physicochemically similar to actives (matched on MW, logP, rotatable
bonds, H-bond donors/acceptors, net charge) but topologically dissimilar — they are **not**
experimental inactives.

**Citation:** Mysinger MM, Carchia M, Irwin JJ, Shoichet BK (2012). Directory of Useful Decoys,
Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking.
*J. Med. Chem.* **55**(14):6582-6594.
DOI: [10.1021/jm300687e](https://doi.org/10.1021/jm300687e)

## Access

- **Website:** https://dude.docking.org/
- **Per-target ISM:** `https://dude.docking.org/targets/<target>/actives_final.ism`
- **License:** Free for academic use (Irwin & Shoichet Labs, UCSF)

## Decoy type

`property-matched` — decoys are synthetic compounds from ZINC matched on physicochemical
properties but **not** experimentally tested. This is the source of known biases.

## Bias warnings

DUD-E is known to exhibit:

1. **Analogue bias:** Actives from the same chemical series are over-represented
2. **Decoy bias:** Property-matched decoys from ZINC may inadvertently exclude known inactives
3. **Enrichment inflation:** ML models exploiting these biases achieve artificially high AUROCs

See: Wallach I et al. (2018). *PLoS ONE* **13**(6):e0220113.
DOI: [10.1371/journal.pone.0220113](https://doi.org/10.1371/journal.pone.0220113)

**Use DUD-E only as a bias-control diagnostic, not as ground truth.**

## Selected targets

| Target | Full name                      | Type | Actives | Decoys | Overlap w/ LIT-PCBA |
|--------|-------------------------------|------|---------|--------|---------------------|
| adrb2  | Beta-2 adrenergic receptor    | GPCR | 231     | 15,000 | Yes (ADRB2)         |
| aa2ar  | Adenosine A2a receptor        | GPCR | 482     | 31,550 | No                  |

**Why these targets?**
- `adrb2` directly overlaps with LIT-PCBA ADRB2, enabling cross-set bias comparison
- `aa2ar` is a well-validated GPCR benchmark target (adenosine receptor, many known ligands)
- Both are GPCRs, relevant to receptor-like pharmacophore optimization

## All 5 GPCR targets in DUD-E

`aa2ar`, `adrb1`, `adrb2`, `cxcr4`, `drd3`

Use `python fetch.py --targets aa2ar adrb1 adrb2 cxcr4 drd3` for all GPCRs.

## File format

```
raw/<target>/actives_final.ism  – SMILES<space>ZINC_ID<space>CHEMBL_ID
raw/<target>/decoys_final.ism   – SMILES<space>ZINC_ID
```

## Usage

```bash
# Download data
python fetch.py                           # default: adrb2 + aa2ar
python fetch.py --targets aa2ar adrb2     # specific targets
python fetch.py --targets aa2ar adrb1 adrb2 cxcr4 drd3   # all GPCR targets

# Load into Python
python load.py                            # prints summary
```

```python
from benchmark.datasets.dude.load import load

actives, decoys, meta = load("adrb2")
# actives: list of SMILES strings (231 compounds from ChEMBL)
# decoys:  list of SMILES strings (15,000 property-matched decoys from ZINC)
# meta:    dict with target, source, doi, decoy_type, n_act, n_dec, bias_warning
```

## Tests

```bash
python test_load.py
```

Verifies: non-empty lists, valid RDKit-parseable SMILES (first 50 checked), meta includes
`bias_warning`, correct `decoy_type = "property-matched"`.

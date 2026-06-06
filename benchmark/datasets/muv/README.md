# MUV — Maximum Unbiased Validation

## What is MUV?

MUV (Maximum Unbiased Validation) is a benchmark for virtual screening
published by Rohrer & Baumann (2009).  It was constructed from PubChem
BioAssay data and designed to be **free from the spatial-statistics biases**
that inflate enrichment metrics in earlier benchmarks such as DUD.

Key design goals:

- **No nearest-neighbour bias**: decoys cannot be close analogues of the
  actives in physicochemical space.
- **No MW/cLogP clustering artefacts**: the property distribution of decoys
  does not cluster tightly around the actives.
- **No topological redundancy**: decoys are topologically dissimilar from
  actives (Tanimoto < 0.6 on MACCS keys).
- **Real negatives**: decoys come from *confirmatory* PubChem assays (i.e.,
  compounds that were actually tested and found inactive), not synthetic
  property-matched compounds.

This makes MUV a stringent, unbiased benchmark for methods that aim to rank
true actives above hard negatives.

## Why max-unbiased?

Property-matched synthetic decoys (e.g., DUD) can inadvertently provide
physicochemical fingerprints that allow a model to distinguish actives from
decoys without learning any target-relevant pharmacophore.  MUV avoids this
by using the **maximum-unbiased selection** algorithm, which iteratively
discards putative decoys that are nearest-neighbours of actives or other
decoys in MW × cLogP space.

As a result, enrichment on MUV is harder to game and more predictive of
real-world virtual screening performance.

## The 17 Target Assays

| MUV ID   | Target / Gene                              |
|----------|--------------------------------------------|
| MUV-466  | SF1 (Steroidogenic factor 1, NR5A1)        |
| MUV-548  | PKA (PRKACA, cAMP-dependent kinase)        |
| MUV-600  | THR (Thyroid hormone receptor β, THRB)     |
| MUV-644  | Rab1 GAP TBC1D20                           |
| MUV-652  | HIV RT (reverse transcriptase p51/p66)     |
| MUV-689  | ERα (Estrogen receptor alpha, ESR1)        |
| MUV-692  | RIN1 (Ras and Rab interactor 1)            |
| MUV-712  | HSP90α (Heat shock protein 90α)            |
| MUV-713  | AR (Androgen receptor, NR3C4)              |
| MUV-733  | GR (Glucocorticoid receptor, NR3C1)        |
| MUV-737  | Akt/PKB (AKT1)                             |
| MUV-810  | FAK (Focal adhesion kinase, PTK2)          |
| MUV-832  | Caspase-3 (CASP3)                          |
| MUV-846  | FXa (Coagulation factor Xa, F10)           |
| MUV-852  | uPA (Urokinase plasminogen activator, PLAU)|
| MUV-858  | D4 receptor (DRD4)                         |
| MUV-859  | SF1 inhibition — second assay (NR5A1)      |

Each assay contains approximately **30 actives** and **15 000 decoys**
(~1:500 ratio).

## File Layout

```
benchmark/datasets/muv/
  fetch.py             — download muv.csv.gz (re-runnable)
  load.py              — load(target) -> (actives, decoys, meta)
  test_load.py         — TDD-lite: ~30 actives, SMILES valid, meta fields
  README.md            — this file
  meta_MUV_*.json      — one meta file per target (written by load.py)
  samples/
    MUV_466_actives_5.smi    — 5 sample active SMILES per target
    MUV_466_decoys_10.smi    — 10 sample decoy SMILES per target
    ...
  raw/                 — gitignored; populated by fetch.py
    muv.csv.gz
```

## Quick Start

```bash
# 1. Download the data (~1.7 MB gzip)
python benchmark/datasets/muv/fetch.py

# 2. Load a target
python - <<'EOF'
from benchmark.datasets.muv.load import load
actives, decoys, meta = load("MUV-466")
print(f"{meta['target']}: {meta['n_act']} actives, {meta['n_dec']} decoys")
print(actives[0])
EOF

# 3. Run the test
cd benchmark/datasets/muv && python test_load.py
```

## Loader Interface (SPEC.md)

```python
from benchmark.datasets.muv.load import load

actives_smiles, decoys_smiles, meta = load("MUV-466")
# actives_smiles  – list[str], ~30 SMILES, label==1
# decoys_smiles   – list[str], ~15 000 SMILES, label==0
# meta            – dict with keys: target, source, doi, decoy_type, n_act, n_dec, url
```

`load()` also accepts the short form (e.g., `load("466")`) and writes
`meta_MUV_<id>.json` to this directory for downstream inspection.

## Citation & License

```
Rohrer SG, Baumann K (2009).
Maximum Unbiased Validation (MUV) Data Sets for Virtual Screening Based
on PubChem Bioassay Data.
J. Chem. Inf. Model. 49(2):169-184.
DOI: 10.1021/ci8002649
```

Original data from **PubChem BioAssay** (public domain, NIH/NLM).
Distributed via MoleculeNet / DeepChem under the MIT License.
No restrictions on use in academic or commercial research.

## Limitations

- ~30 actives per assay is a very small positive set — enrichment estimates
  have high variance; report confidence intervals.
- The 1:500 ratio means even modest re-ranking scores look good; compare
  against a random-ranking baseline.
- Decoys are real negatives from PubChem, so some may have weak off-target
  activity — they are not guaranteed inactive against the assay target.

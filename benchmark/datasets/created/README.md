# Created: ChEMBL-Derived Actives + MUV-Style Unbiased Decoys

Self-contained, fully-reproducible benchmark datasets for 3 GPCR targets
built from public REST APIs. No manual downloads required.

## Targets

| Target | ChEMBL ID    | UniProt  | Name                           |
|--------|-------------|----------|--------------------------------|
| CCR2   | CHEMBL4015  | P41597   | C-C chemokine receptor type 2  |
| CCR5   | CHEMBL274   | P51681   | C-C chemokine receptor type 5  |
| CXCR4  | CHEMBL2107  | P61073   | C-X-C chemokine receptor type 4|

## Dataset Construction

### Actives
- Source: ChEMBL REST API (`/activity?target_chembl_id=...&pchembl_value__gte=7`)
- Filter: pChEMBL ≥ 7 (IC50 / Ki ≤ 100 nM)
- Deduplication: by canonical SMILES (RDKit)
- All SMILES validated with `Chem.MolFromSmiles`

### Decoys (MUV-style unbiased)
Method follows Rohrer & Baumann (2009) *J Chem Inf Model* 49:169–184.

1. **Background pool**: drug-like ChEMBL molecules (MW 200–600 Da), sampled
   via the ChEMBL REST API.
2. **Property matching** per active: MW ± 50 Da, ALogP ± 1.5, HBD ± 1, HBA ± 1
   (computed with RDKit).
3. **Similarity unbiasing**: maximum Morgan-Tc (radius=2, 2048 bits) < 0.35
   to **any** active across **all three targets** (cross-target unbiasing).
4. Target: up to 50 decoys per active; fewer accepted if the background is sparse.

## Files

```
created/
  build.py          — reproducible builder (fetches from ChEMBL + constructs decoys)
  load.py           — SPEC-conformant loader
  test_load.py      — TDD-lite test suite
  README.md         — this file
  CCR2/
    actives.csv     — col: smiles
    decoys.csv      — col: smiles
    meta.json       — metadata
  CCR5/
    actives.csv
    decoys.csv
    meta.json
  CXCR4/
    actives.csv
    decoys.csv
    meta.json
```

## Usage

### Rebuild datasets

```bash
# All three targets (default)
python build.py

# Single target
python build.py --targets CCR2

# Custom options
python build.py --targets CCR2 CCR5 CXCR4 \
    --pool-size 50000 \
    --max-decoys-per-active 50 \
    --tc-cutoff 0.35
```

### Load in Python

```python
from benchmark.datasets.created import load

actives, decoys, meta = load.load("CCR2")
print(f"CCR2: {len(actives)} actives, {len(decoys)} decoys")
# CCR2: 591 actives, 3245 decoys  (numbers vary with ChEMBL version)
```

### Run tests

```bash
python test_load.py --target CCR2        # single target
python test_load.py --all                # all available targets
```

## API Rate Limits

ChEMBL REST API has no strict rate limit; the build script keeps requests
below ~10/sec. A full build for all 3 targets (50k background pool) takes
approximately 5–15 minutes depending on network latency.

## Reproducibility Notes

- ChEMBL data evolves as new assay data is added; re-running `build.py`
  will fetch the current version of the database.
- The `meta.json` for each target records `n_act` and `n_dec` at build time.
- The background pool is sampled by sequential offset (not random seed), so
  re-runs produce the same molecules if ChEMBL's ordering is stable.

## Reference

Rohrer, S. G. & Baumann, K. (2009). Maximum Unbiased Validation (MUV) Data
Sets for Virtual Screening Based on PubChem Bioactivity Data. *J Chem Inf
Model*, 49(2), 169–184. https://doi.org/10.1021/ci8002649

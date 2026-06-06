# Benchmark Datasets SPEC

## Uniform loader interface

Each dataset directory (`litpcba/`, `dude/`, `created/`) must provide a `load.py` with:

```python
def load(target: str, split: str = "all") -> tuple[list[str], list[str], dict]:
    """
    Returns:
        actives_smiles  – list of canonical SMILES strings (actives / confirmed binders)
        decoys_smiles   – list of canonical SMILES strings (inactives / decoys)
        meta            – dict matching the meta.json schema below
    """
```

## meta.json schema (per target)

```json
{
  "target":      "<target_id>",
  "source":      "<dataset_name>",
  "doi":         "<doi_string>",
  "decoy_type":  "experimental-inactive | property-matched | created-muv-unbiased | ...",
  "n_act":       123,
  "n_dec":       45678,
  "url":         "<download_url>"
}
```

## File layout

```
benchmark/
  datasets/
    SPEC.md                  ← this file
    .gitignore               ← ignore large raw files (*.sdf.gz, *.mol2.gz, etc.)
    litpcba/
      fetch.py               ← downloader (re-runnable, graceful HTTP failures)
      load.py                ← returns (actives, decoys, meta)
      README.md
      meta_<TARGET>.json     ← one per downloaded target
      raw/                   ← gitignored large files
        <TARGET>/
          actives.smi
          inactives.smi
    dude/
      fetch.py
      load.py
      README.md
      meta_<target>.json
      raw/
        <target>/
          actives_final.ism
          decoys_final.ism
    created/
      build.py               ← reproducible builder (REST API + RDKit)
      load.py                ← returns (actives, decoys, meta)
      test_load.py           ← TDD-lite tests
      README.md
      <TARGET>/
        actives.csv          ← col: smiles
        decoys.csv           ← col: smiles
        meta.json
```

## SMILES validity

SMILES strings must be parseable by RDKit (`Chem.MolFromSmiles` not returning None).

## Test convention

Each dataset directory has a `test_load.py` that imports `load.py` and asserts:
- Both returned lists are non-empty
- All SMILES are RDKit-parseable (check first 50 from each list)

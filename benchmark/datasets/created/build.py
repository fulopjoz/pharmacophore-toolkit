"""
build.py — Build ChEMBL-derived actives + MUV-style unbiased decoys
for 3 GPCR targets: CCR2, CCR5, CXCR4.

Actives:   ChEMBL REST API, pChEMBL >= 7 (IC50/Ki <= 100 nM), deduplicated by
           canonical SMILES.

Decoys:    MUV-style property-matched unbiased decoys (Rohrer & Baumann, 2009):
           - Background pool: random ChEMBL drug-like molecules (sampled from
             ChEMBL molecule endpoint, max_phase >= 0, MW 200-600).
           - Property matching per active: MW ± 50, ALogP ± 1.5, HBD ± 1, HBA ± 1.
           - Similarity filter: max Morgan-Tc (radius=2, nBits=2048) < 0.35 to
             ANY active across ALL three targets (cross-target unbiasing).
           - Target: 20-50 decoys per active; fewer accepted if background is sparse.

Outputs (per target in created/<TARGET>/):
    actives.csv   — col: smiles
    decoys.csv    — col: smiles
    meta.json     — SPEC-conformant metadata

Usage:
    python build.py [--targets CCR2 CCR5 CXCR4] [--pool-size 50000]
                    [--max-decoys-per-active 50] [--tc-cutoff 0.35]
                    [--out-dir <path>]
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import requests
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import DataStructs
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Target registry
# ---------------------------------------------------------------------------

TARGETS = {
    "CCR2": {
        "chembl_id": "CHEMBL4015",
        "uniprot_id": "P41597",
        "name": "C-C chemokine receptor type 2",
    },
    "CCR5": {
        "chembl_id": "CHEMBL274",
        "uniprot_id": "P51681",
        "name": "C-C chemokine receptor type 5",
    },
    "CXCR4": {
        "chembl_id": "CHEMBL2107",
        "uniprot_id": "P61073",
        "name": "C-X-C chemokine receptor type 4",
    },
}

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
REQUEST_DELAY = 0.12   # seconds between API calls — stay well under 10 req/s
PAGE_SIZE = 500        # activities per page (ChEMBL max 1000)
MOL_PAGE_SIZE = 1000   # molecules per page for background pool

# ---------------------------------------------------------------------------
# ChEMBL API helpers
# ---------------------------------------------------------------------------

def _get_json(url: str, params: dict | None = None, retries: int = 3) -> dict:
    """GET JSON from ChEMBL REST API with retry on transient errors."""
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, params=params, timeout=30)
            if r.status_code == 429:
                wait = 30 * attempt
                log.warning("Rate-limited; sleeping %d s", wait)
                time.sleep(wait)
                continue
            r.raise_for_status()
            time.sleep(REQUEST_DELAY)
            return r.json()
        except requests.RequestException as exc:
            log.warning("Request failed (attempt %d/%d): %s", attempt, retries, exc)
            if attempt < retries:
                time.sleep(5 * attempt)
    raise RuntimeError(f"Failed to fetch {url} after {retries} attempts")


def fetch_actives(chembl_target_id: str, pchembl_min: float = 7.0) -> list[str]:
    """
    Fetch all active SMILES for a target from ChEMBL.

    Returns deduplicated canonical SMILES (RDKit-parsed) for molecules with
    pChEMBL value >= pchembl_min (i.e. IC50/Ki <= 10^-pchembl_min M).
    """
    log.info("Fetching actives for %s (pChEMBL >= %.1f)...", chembl_target_id, pchembl_min)

    seen_smiles: set[str] = set()
    actives: list[str] = set()   # will be converted to list
    offset = 0
    total = None

    while True:
        params = {
            "target_chembl_id": chembl_target_id,
            "pchembl_value__gte": pchembl_min,
            "limit": PAGE_SIZE,
            "offset": offset,
            "format": "json",
            "only": "canonical_smiles,pchembl_value,molecule_chembl_id",
        }
        data = _get_json(f"{CHEMBL_BASE}/activity", params=params)
        meta = data.get("page_meta", {})
        if total is None:
            total = meta.get("total_count", 0)
            log.info("  Total activity records: %d", total)

        for rec in data.get("activities", []):
            smi = rec.get("canonical_smiles")
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            canon = Chem.MolToSmiles(mol, canonical=True)
            seen_smiles.add(canon)

        offset += PAGE_SIZE
        fetched = offset
        log.info("  Fetched %d / %d activity records (unique SMILES so far: %d)",
                 min(fetched, total), total, len(seen_smiles))

        # Check if we've fetched everything
        if not data.get("activities") or offset >= total:
            break

    result = sorted(seen_smiles)
    log.info("  -> %d unique canonical actives for %s", len(result), chembl_target_id)
    return result


# ---------------------------------------------------------------------------
# RDKit property helpers
# ---------------------------------------------------------------------------

def _mol_props(mol) -> dict:
    """Return MW, ALogP, HBD, HBA for a molecule."""
    return {
        "mw": Descriptors.ExactMolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
    }


def _morgan_fp(mol):
    return GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def _max_tc(fp, ref_fps: list) -> float:
    """Return maximum Tanimoto similarity of fp to any fingerprint in ref_fps."""
    if not ref_fps:
        return 0.0
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    return max(sims)


# ---------------------------------------------------------------------------
# Background pool from ChEMBL
# ---------------------------------------------------------------------------

def fetch_background_pool(
    max_molecules: int = 50_000,
    mw_range: tuple[float, float] = (200.0, 600.0),
) -> tuple[list[str], list]:
    """
    Fetch a random background pool of drug-like molecules from ChEMBL.

    Samples ChEMBL molecules with MW in [200, 600] Da.
    Returns (smiles_list, mol_list) — both RDKit-parsed.

    ChEMBL doesn't support random sampling; we page through by offset
    using a fixed ordering to get a diverse cross-section.
    """
    log.info("Fetching background pool (up to %d molecules)...", max_molecules)

    pool_smiles: list[str] = []
    pool_mols: list = []
    seen: set[str] = set()
    offset = 0
    step = MOL_PAGE_SIZE

    while len(pool_smiles) < max_molecules:
        params = {
            "molecule_properties__full_mwt__gte": mw_range[0],
            "molecule_properties__full_mwt__lte": mw_range[1],
            "limit": step,
            "offset": offset,
            "format": "json",
            "only": "molecule_structures,molecule_properties,molecule_chembl_id",
        }
        data = _get_json(f"{CHEMBL_BASE}/molecule", params=params)
        mols_data = data.get("molecules", [])
        if not mols_data:
            break

        for rec in mols_data:
            structs = rec.get("molecule_structures") or {}
            smi = structs.get("canonical_smiles")
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            canon = Chem.MolToSmiles(mol, canonical=True)
            if canon in seen:
                continue
            seen.add(canon)
            pool_smiles.append(canon)
            pool_mols.append(mol)

        offset += step
        log.info("  Pool size: %d / %d (offset %d)", len(pool_smiles), max_molecules, offset)

        # ChEMBL has ~2M molecules; don't scan more than needed
        if offset > 200_000:
            log.info("  Reached offset limit (200,000); stopping pool fetch.")
            break

    log.info("Background pool: %d unique molecules", len(pool_smiles))
    return pool_smiles, pool_mols


# ---------------------------------------------------------------------------
# Decoy selection — MUV-style unbiased
# ---------------------------------------------------------------------------

def select_decoys(
    actives_smiles: list[str],
    all_actives_smiles: list[str],   # cross-target: used for Tc filter
    pool_smiles: list[str],
    pool_mols: list,
    max_decoys_per_active: int = 50,
    tc_cutoff: float = 0.35,
    mw_tol: float = 50.0,
    logp_tol: float = 1.5,
    hbd_tol: int = 1,
    hba_tol: int = 1,
) -> list[str]:
    """
    Select property-matched, similarity-unbiased decoys.

    For each active, find pool molecules that:
      1. Property-match (MW±mw_tol, ALogP±logp_tol, HBD±hbd_tol, HBA±hba_tol)
      2. Have max Tanimoto Tc (Morgan r2, 2048 bits) < tc_cutoff to ALL actives
         (using all_actives_smiles, which spans all three targets for cross-target
          unbiasing)

    Returns deduplicated decoy SMILES.
    """
    log.info("  Computing fingerprints for %d actives (cross-target)...", len(all_actives_smiles))

    # Fingerprints for ALL actives (cross-target unbiasing)
    active_fps = []
    for smi in all_actives_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            active_fps.append(_morgan_fp(mol))

    log.info("  Computing fingerprints for %d active properties...", len(actives_smiles))
    # Properties for THIS target's actives (for property matching)
    active_props = []
    for smi in actives_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            active_props.append(_mol_props(mol))

    log.info("  Computing fingerprints for %d pool molecules...", len(pool_mols))
    # Pool fingerprints (compute once)
    pool_fps = []
    pool_props_list = []
    for mol in pool_mols:
        pool_fps.append(_morgan_fp(mol))
        pool_props_list.append(_mol_props(mol))

    # --- Filter pool by similarity (bulk Tc against all actives) ---
    log.info("  Filtering pool by max Tc < %.2f to all actives...", tc_cutoff)
    eligible_indices = []
    batch = 2000
    for i in range(0, len(pool_fps), batch):
        batch_fps = pool_fps[i:i + batch]
        for j, fp in enumerate(batch_fps):
            max_sim = _max_tc(fp, active_fps)
            if max_sim < tc_cutoff:
                eligible_indices.append(i + j)
        if i % 10000 == 0 and i > 0:
            log.info("    Processed %d / %d pool mols, eligible so far: %d",
                     i, len(pool_fps), len(eligible_indices))

    log.info("  Eligible pool molecules (Tc < %.2f): %d / %d",
             tc_cutoff, len(eligible_indices), len(pool_mols))

    # --- Property-match decoys to each active ---
    # A pool molecule can serve as a decoy for multiple actives.
    # We accumulate a union; the set.add() is idempotent so duplicates are
    # collapsed automatically.  The per-active cap (max_decoys_per_active)
    # just limits how many we assign to a single active to keep the ratio
    # bounded; it does NOT block a molecule from being selected by other actives.
    selected: set[str] = set()

    for act_idx, props in enumerate(active_props):
        count = 0
        for idx in eligible_indices:
            p = pool_props_list[idx]
            if (abs(p["mw"] - props["mw"]) <= mw_tol
                    and abs(p["logp"] - props["logp"]) <= logp_tol
                    and abs(p["hbd"] - props["hbd"]) <= hbd_tol
                    and abs(p["hba"] - props["hba"]) <= hba_tol):
                selected.add(pool_smiles[idx])
                count += 1
            if count >= max_decoys_per_active:
                break

        if (act_idx + 1) % 50 == 0:
            log.info("    Matched decoys after %d actives: %d", act_idx + 1, len(selected))

    log.info("  Total selected decoys: %d", len(selected))
    return sorted(selected)


# ---------------------------------------------------------------------------
# CSV / JSON writers
# ---------------------------------------------------------------------------

def write_csv(smiles_list: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("smiles\n")
        for smi in smiles_list:
            f.write(smi + "\n")


def write_meta(meta: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(meta, f, indent=2)


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def build_target(
    target_key: str,
    all_actives_smiles: list[str],
    pool_smiles: list[str],
    pool_mols: list,
    out_dir: Path,
    max_decoys_per_active: int = 50,
    tc_cutoff: float = 0.35,
) -> dict:
    """Build actives.csv + decoys.csv + meta.json for one target."""
    info = TARGETS[target_key]
    log.info("=" * 60)
    log.info("Building dataset for %s (%s)", target_key, info["chembl_id"])
    log.info("=" * 60)

    tgt_dir = out_dir / target_key
    tgt_dir.mkdir(parents=True, exist_ok=True)

    # Fetch actives
    actives = fetch_actives(info["chembl_id"], pchembl_min=7.0)
    write_csv(actives, tgt_dir / "actives.csv")
    log.info("Wrote %d actives -> %s/actives.csv", len(actives), target_key)

    # Select decoys (unbiased against ALL actives across all targets)
    log.info("Selecting decoys for %s...", target_key)
    decoys = select_decoys(
        actives_smiles=actives,
        all_actives_smiles=all_actives_smiles,
        pool_smiles=pool_smiles,
        pool_mols=pool_mols,
        max_decoys_per_active=max_decoys_per_active,
        tc_cutoff=tc_cutoff,
    )
    write_csv(decoys, tgt_dir / "decoys.csv")
    log.info("Wrote %d decoys -> %s/decoys.csv", len(decoys), target_key)

    meta = {
        "target": target_key,
        "name": info["name"],
        "chembl_id": info["chembl_id"],
        "uniprot_id": info["uniprot_id"],
        "source": "ChEMBL REST API",
        "source_url": "https://www.ebi.ac.uk/chembl/api/data/activity",
        "decoy_type": "created-muv-unbiased",
        "decoy_note": (
            "MUV-style unbiased decoys: property-matched (MW±50, ALogP±1.5, HBD±1, HBA±1) "
            f"+ max Morgan-Tc(r=2,2048bits) < {tc_cutoff} to all actives across all 3 targets."
        ),
        "pchembl_threshold": 7.0,
        "n_act": len(actives),
        "n_dec": len(decoys),
        "decoy_tc_cutoff": tc_cutoff,
        "decoys_per_active_target": max_decoys_per_active,
        "decoy_background": "ChEMBL drug-like pool (MW 200-600 Da)",
        "reference": "Rohrer & Baumann (2009) J Chem Inf Model 49:169-184",
    }
    write_meta(meta, tgt_dir / "meta.json")
    log.info("Wrote meta.json for %s", target_key)

    return meta


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", nargs="+", default=["CCR2", "CCR5", "CXCR4"],
                        choices=list(TARGETS), help="Targets to build (default: all 3)")
    parser.add_argument("--pool-size", type=int, default=50_000,
                        help="Max background pool molecules to fetch from ChEMBL (default: 50000)")
    parser.add_argument("--max-decoys-per-active", type=int, default=50,
                        help="Max decoys to assign per active (default: 50)")
    parser.add_argument("--tc-cutoff", type=float, default=0.35,
                        help="Max Morgan-Tc to any active (default: 0.35)")
    parser.add_argument("--out-dir", type=Path,
                        default=Path(__file__).parent,
                        help="Output root directory (default: same dir as build.py)")
    args = parser.parse_args(argv)

    out_dir = args.out_dir.resolve()
    log.info("Output directory: %s", out_dir)
    log.info("Targets: %s", args.targets)

    # Step 1: fetch all actives for all requested targets (for cross-target Tc filter)
    log.info("Step 1: Fetching all actives for cross-target similarity unbiasing...")
    all_actives_by_target: dict[str, list[str]] = {}
    for tgt in args.targets:
        info = TARGETS[tgt]
        all_actives_by_target[tgt] = fetch_actives(info["chembl_id"], pchembl_min=7.0)

    all_actives_smiles = []
    for smiles in all_actives_by_target.values():
        all_actives_smiles.extend(smiles)
    all_actives_smiles = list(set(all_actives_smiles))
    log.info("Total unique actives across all targets: %d", len(all_actives_smiles))

    # Step 2: fetch background pool once (shared across all targets)
    log.info("Step 2: Fetching background pool (max %d molecules)...", args.pool_size)
    pool_smiles, pool_mols = fetch_background_pool(max_molecules=args.pool_size)

    # Remove from pool any SMILES that appear in actives
    active_set = set(all_actives_smiles)
    pool_smiles_f, pool_mols_f = zip(
        *[(s, m) for s, m in zip(pool_smiles, pool_mols) if s not in active_set]
    ) if pool_smiles else ([], [])
    pool_smiles_f = list(pool_smiles_f)
    pool_mols_f = list(pool_mols_f)
    log.info("Pool after removing actives: %d molecules", len(pool_smiles_f))

    # Step 3: build each target
    summary = {}
    for tgt in args.targets:
        meta = build_target(
            target_key=tgt,
            all_actives_smiles=all_actives_smiles,
            pool_smiles=pool_smiles_f,
            pool_mols=pool_mols_f,
            out_dir=out_dir,
            max_decoys_per_active=args.max_decoys_per_active,
            tc_cutoff=args.tc_cutoff,
        )
        summary[tgt] = meta

    log.info("")
    log.info("=" * 60)
    log.info("BUILD COMPLETE")
    log.info("=" * 60)
    for tgt, meta in summary.items():
        log.info("  %s: %d actives, %d decoys", tgt, meta["n_act"], meta["n_dec"])


if __name__ == "__main__":
    main()

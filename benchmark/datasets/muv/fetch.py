"""
Fetch MUV (Maximum Unbiased Validation) benchmark data.

Source: Rohrer SG, Baumann K (2009).
        Maximum Unbiased Validation (MUV) Data Sets for Virtual Screening
        Based on PubChem Bioassay Data.
        J. Chem. Inf. Model. 49(2):169-184.
        DOI: 10.1021/ci8002649

Distribution: MoleculeNet / DeepChem S3 (original source: PubChem BioAssay)
  https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/muv.csv.gz

The file is a single CSV (~1.7 MB gzipped, ~10 MB expanded) with:
  - 17 columns for MUV assays (MUV-466 … MUV-859): 1=active, 0=inactive, blank=unscreened
  - 'mol_id'  column: PubChem CID
  - 'smiles'  column: isomeric SMILES

Each assay has ~30 actives and ~15,000 spatially-unbiased decoys (max-unbiased selection).

Re-runnable: if raw/muv.csv.gz already exists it is not re-downloaded unless --force is passed.
"""

import argparse
import sys
import urllib.error
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SOURCE_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/muv.csv.gz"
RAW_DIR    = Path(__file__).parent / "raw"
RAW_FILE   = RAW_DIR / "muv.csv.gz"

MUV_TARGETS = [
    "MUV-466", "MUV-548", "MUV-600", "MUV-644", "MUV-652",
    "MUV-689", "MUV-692", "MUV-712", "MUV-713", "MUV-733",
    "MUV-737", "MUV-810", "MUV-832", "MUV-846", "MUV-852",
    "MUV-858", "MUV-859",
]

# ---------------------------------------------------------------------------
# Downloader
# ---------------------------------------------------------------------------

def fetch(force: bool = False) -> bool:
    """
    Download muv.csv.gz to raw/muv.csv.gz.

    Parameters
    ----------
    force : bool
        Re-download even if file already present.

    Returns
    -------
    True if file is present after the call, False on failure.
    """
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    if RAW_FILE.exists() and not force:
        size_mb = RAW_FILE.stat().st_size / 1_048_576
        print(f"  muv.csv.gz already present ({size_mb:.1f} MB) — skipping "
              "(use --force to re-download)")
        return True

    print(f"  GET {SOURCE_URL}", flush=True)
    try:
        with urllib.request.urlopen(SOURCE_URL, timeout=120) as resp:
            if resp.status != 200:
                print(f"  HTTP {resp.status} — download failed", file=sys.stderr)
                return False
            data = resp.read()
        with open(RAW_FILE, "wb") as f:
            f.write(data)
        size_mb = len(data) / 1_048_576
        print(f"  Saved {RAW_FILE} ({size_mb:.1f} MB)")
        return True
    except urllib.error.HTTPError as e:
        print(f"  HTTPError {e.code} {e.reason}: {SOURCE_URL}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"  Download failed ({SOURCE_URL}): {e}", file=sys.stderr)
        return False


def is_downloaded() -> bool:
    """Return True if raw/muv.csv.gz is present and non-empty."""
    return RAW_FILE.exists() and RAW_FILE.stat().st_size > 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Fetch MUV benchmark data (single CSV, all 17 targets)"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download even if muv.csv.gz already present",
    )
    args = parser.parse_args()

    ok = fetch(force=args.force)
    if ok:
        print(f"\nMUV data ready — {len(MUV_TARGETS)} targets available: "
              f"{', '.join(MUV_TARGETS)}")
    else:
        print("\nFetch FAILED — check network connection or URL.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

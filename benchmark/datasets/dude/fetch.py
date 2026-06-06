"""
Fetch DUD-E benchmark data.

Source: Mysinger MM, Carchia M, Irwin JJ, Shoichet BK (2012).
        Directory of Useful Decoys, Enhanced (DUD-E):
        Better Ligands and Decoys for Better Benchmarking.
        J. Med. Chem. 55(14):6582-6594.
        DOI: 10.1021/jm300687e

Distribution: https://dude.docking.org/ (Irwin & Shoichet Laboratories, UCSF)

Per-target ISM (Isomeric SMILES) files:
  actives_final.ism   SMILES<space>ZINC_ID<space>CHEMBL_ID
  decoys_final.ism    SMILES<space>ZINC_ID

⚠️  BIAS WARNING: DUD-E decoys are property-matched synthetic decoys generated from ZINC,
not experimental inactives. This can introduce enrichment bias (analogue bias, decoy bias).
Use for comparison / bias-control diagnostic only — not as a ground-truth inactive set.

Selected targets (GPCR/receptor-like, overlap with LIT-PCBA where possible):
  adrb2  – Beta-2 adrenergic receptor  (OVERLAP with LIT-PCBA ADRB2 ← key for bias-control)
  aa2ar  – Adenosine A2a receptor       (GPCR, additional receptor target)

Re-runnable: already-downloaded targets are skipped unless --force is passed.
"""

import argparse
import sys
from pathlib import Path

# ---- configuration --------------------------------------------------------

BASE_URL = "https://dude.docking.org/targets"

# GPCR targets in DUD-E
GPCR_TARGETS = ["aa2ar", "adrb1", "adrb2", "cxcr4", "drd3"]

# Selected targets for this benchmark
DEFAULT_TARGETS = ["adrb2", "aa2ar"]

# All 102 DUD-E targets (subset listing; download all via dude.docking.org)
# fmt: off
ALL_TARGETS = [
    "aa2ar", "abl1", "ace", "aces", "ada", "ada17", "adrb1", "adrb2",
    "akt1", "akt2", "aldr", "ampc", "andr", "aofb", "bace1", "braf",
    "cah2", "casp3", "cdk2", "comt", "cp2c9", "cp3a4", "csf1r", "cxcr4",
    "def", "dhi1", "dpp4", "drd3", "dyr", "egfr", "esr1", "esr2",
    "fa10", "fa7", "fabp4", "fak1", "fkb1a", "fnta", "fpps", "gcr",
    "glcm", "gria2", "grik1", "hdac2", "hdac8", "hivint", "hivpr", "hivrt",
    "hmdh", "hs90a", "hxk4", "igf1r", "inha", "ital", "jak2", "kif11",
    "kit", "kith", "kpcb", "lck", "lfa1", "mapk2", "mcr", "met",
    "mk01", "mk10", "mk14", "mmp13", "mp2k1", "nos1", "nram", "pa2ga",
    "parp1", "pde5a", "pgh1", "pgh2", "plk1", "pnph", "ppara", "ppard",
    "pparg", "prgr", "ptn1", "pur2", "pygm", "pyrd", "reni", "rock1",
    "rxra", "sahh", "src", "tgfr1", "thb", "thr", "tpk1", "tryb1",
    "try1", "tysy", "urok", "vgfr2", "wee1", "xiap",
]
# fmt: on

RAW_DIR = Path(__file__).parent / "raw"

# ---- helpers ---------------------------------------------------------------

def _try_download(url: str, dest: Path, timeout: int = 60) -> bool:
    """Download url → dest. Returns True on success, False on HTTP error."""
    import urllib.request
    import urllib.error

    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        print(f"  GET {url}", flush=True)
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            if resp.status != 200:
                print(f"  HTTP {resp.status} – skipping", file=sys.stderr)
                return False
            data = resp.read()
        with open(dest, "wb") as f:
            f.write(data)
        print(f"  Saved {dest} ({len(data):,} bytes)")
        return True
    except urllib.error.HTTPError as e:
        print(f"  HTTPError {e.code} {e.reason}: {url}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"  Download failed ({url}): {e}", file=sys.stderr)
        return False


def fetch(targets: list[str] | None = None, force: bool = False) -> dict[str, bool]:
    """
    Download DUD-E actives_final.ism + decoys_final.ism for each target.

    Parameters
    ----------
    targets : list of target IDs (lowercase, e.g. ["adrb2", "aa2ar"])
    force   : re-download even if already present

    Returns
    -------
    dict mapping target_id → True if both ISM files are present
    """
    if targets is None:
        targets = DEFAULT_TARGETS

    results = {}
    for target in targets:
        t_dir = RAW_DIR / target
        act_path = t_dir / "actives_final.ism"
        dec_path = t_dir / "decoys_final.ism"

        if not force and act_path.exists() and dec_path.exists():
            print(f"  {target}: already present (use --force to re-download)")
            results[target] = True
            continue

        print(f"\nFetching {target} from DUD-E...")
        t_dir.mkdir(parents=True, exist_ok=True)

        act_ok = _try_download(f"{BASE_URL}/{target}/actives_final.ism", act_path)
        dec_ok = _try_download(f"{BASE_URL}/{target}/decoys_final.ism", dec_path)
        results[target] = act_ok and dec_ok

    return results


# ---- CLI ------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Fetch DUD-E benchmark data")
    parser.add_argument(
        "--targets", nargs="+", default=DEFAULT_TARGETS,
        metavar="TARGET",
        help=f"Target IDs to fetch (default: {DEFAULT_TARGETS}). "
             f"GPCR targets: {GPCR_TARGETS}",
    )
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if files already present")
    args = parser.parse_args()

    results = fetch(targets=args.targets, force=args.force)
    print("\nResults:")
    for t, ok in results.items():
        status = "OK" if ok else "FAILED"
        print(f"  {t}: {status}")
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()

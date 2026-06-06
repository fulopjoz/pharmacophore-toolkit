"""
Fetch LIT-PCBA benchmark data.

Source: Tran-Nguyen V-K, Jacquemard C, Rognan D (2020).
        LIT-PCBA: An Unbiased Data Set for Machine Learning and Virtual Screening.
        J. Chem. Inf. Model. 60(9):4263-4273.
        DOI: 10.1021/acs.jcim.0c00155

Distribution: http://drugdesign.unistra.fr/LIT-PCBA/
              (Laboratoire d'Innovation Thérapeutique, Université de Strasbourg)

The full_data.tgz archive (~54 MB) contains 15 target directories, each with:
  - actives.smi    SMILES<tab>PubChem_CID  (confirmed actives from dose-response assays)
  - inactives.smi  SMILES<tab>PubChem_CID  (confirmed inactives)
  - *.mol2         3D structures (protein + reference ligand poses)

All actives and inactives are within similar molecular property ranges (unbiased by AVE).

Selected targets (GPCR/receptor-like):
  ADRB2  – Beta-2 adrenergic receptor (also in DUD-E → overlap for bias-control)
  OPRK1  – Kappa opioid receptor (GPCR, opioid family)

Re-runnable: already-downloaded targets are skipped unless --force is passed.
"""

import argparse
import sys
import tarfile
from pathlib import Path

# ---- configuration --------------------------------------------------------

FULL_DATA_URL = "http://drugdesign.unistra.fr/LIT-PCBA/Files/full_data.tgz"
# Wayback Machine mirror (used as fallback — snapshot from 2026-04-11)
WAYBACK_URL = (
    "https://web.archive.org/web/20260411170149/"
    "https://drugdesign.unistra.fr/LIT-PCBA/Files/full_data.tgz"
)

# GPCR/receptor targets selected for this benchmark
DEFAULT_TARGETS = ["ADRB2", "OPRK1"]

# All 15 targets in the dataset (for reference)
ALL_TARGETS = [
    "ADRB2",   # Beta-2 adrenergic receptor (GPCR, agonists)
    "ALDH1",   # Aldehyde dehydrogenase 1
    "ESR1_ago", # Estrogen receptor α (agonists)
    "ESR1_ant", # Estrogen receptor α (antagonists)
    "FEN1",    # Flap endonuclease 1
    "GBA",     # Glucocerebrosidase
    "IDH1",    # Isocitrate dehydrogenase 1
    "KAT2A",   # Histone acetyltransferase KAT2A
    "MAPK1",   # Mitogen-activated protein kinase 1 (ERK2)
    "MTORC1",  # mTOR complex 1
    "OPRK1",   # Kappa opioid receptor (GPCR, agonists)
    "PKM2",    # Pyruvate kinase M2
    "PPARG",   # Peroxisome proliferator-activated receptor γ
    "TP53",    # Tumor protein p53
    "VDR",     # Vitamin D receptor
]

RAW_DIR = Path(__file__).parent / "raw"

# ---- helpers ---------------------------------------------------------------

def _try_download(url: str, dest: Path, timeout: int = 120) -> bool:
    """Download url → dest. Returns True on success, False on HTTP error."""
    import urllib.request
    import urllib.error

    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        print(f"  Downloading {url}", flush=True)
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            if resp.status != 200:
                print(f"  HTTP {resp.status} – skipping", file=sys.stderr)
                return False
            total = int(resp.headers.get("Content-Length", 0))
            downloaded = 0
            with open(dest, "wb") as f:
                while True:
                    chunk = resp.read(1 << 17)  # 128 KB
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded * 100 // total
                        print(f"\r  {pct:3d}% ({downloaded/1e6:.1f}/{total/1e6:.1f} MB)",
                              end="", flush=True)
        print()
        return True
    except urllib.error.HTTPError as e:
        print(f"  HTTPError {e.code} {e.reason}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"  Download failed: {e}", file=sys.stderr)
        return False


def _extract_targets(archive_path: Path, targets: list[str], dest_dir: Path) -> dict[str, bool]:
    """Extract <TARGET>/actives.smi + inactives.smi for each target. Returns {target: ok}."""
    results = {}
    with tarfile.open(archive_path, "r:gz") as tf:
        members = tf.getnames()
        for target in targets:
            out_dir = dest_dir / target
            out_dir.mkdir(parents=True, exist_ok=True)
            ok = True
            for fname in ("actives.smi", "inactives.smi"):
                member_path = f"{target}/{fname}"
                if member_path not in members:
                    print(f"  WARNING: {member_path} not found in archive", file=sys.stderr)
                    ok = False
                    continue
                dest_file = out_dir / fname
                with tf.extractfile(member_path) as src, open(dest_file, "wb") as dst:
                    dst.write(src.read())
                print(f"  Extracted {member_path} → {dest_file}")
            results[target] = ok
    return results


def fetch(targets: list[str] | None = None, force: bool = False) -> dict[str, bool]:
    """
    Download and extract LIT-PCBA SMILES for the given targets.

    Parameters
    ----------
    targets : list of target IDs (default: DEFAULT_TARGETS)
    force   : re-download even if already present

    Returns
    -------
    dict mapping target_id → True if both actives.smi + inactives.smi are present
    """
    if targets is None:
        targets = DEFAULT_TARGETS

    # Check which targets are already done
    if not force:
        missing = []
        for t in targets:
            act = RAW_DIR / t / "actives.smi"
            inact = RAW_DIR / t / "inactives.smi"
            if act.exists() and inact.exists():
                print(f"  {t}: already present (use --force to re-download)")
            else:
                missing.append(t)
        if not missing:
            return {t: True for t in targets}
        targets_to_fetch = missing
    else:
        targets_to_fetch = targets

    # Download the full archive (it's needed for all targets regardless)
    archive = Path("/tmp") / "litpcba_full_data.tgz"
    if force or not archive.exists():
        ok = _try_download(FULL_DATA_URL, archive)
        if not ok:
            print("  Primary URL failed, trying Wayback Machine mirror...", file=sys.stderr)
            ok = _try_download(WAYBACK_URL, archive)
        if not ok:
            print("ERROR: Could not download LIT-PCBA archive from any source.", file=sys.stderr)
            print("Manual download:", FULL_DATA_URL, file=sys.stderr)
            return {t: False for t in targets_to_fetch}
    else:
        print(f"  Using cached archive: {archive}")

    # Extract requested targets
    return _extract_targets(archive, targets_to_fetch, RAW_DIR)


# ---- CLI ------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Fetch LIT-PCBA benchmark data")
    parser.add_argument(
        "--targets", nargs="+", default=DEFAULT_TARGETS,
        metavar="TARGET",
        help=f"Target IDs to fetch (default: {DEFAULT_TARGETS}). "
             f"Available: {ALL_TARGETS}",
    )
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if files are already present")
    parser.add_argument("--all", dest="all_targets", action="store_true",
                        help="Fetch all 15 LIT-PCBA targets")
    args = parser.parse_args()

    targets = ALL_TARGETS if args.all_targets else args.targets
    results = fetch(targets=targets, force=args.force)
    print("\nResults:")
    for t, ok in results.items():
        status = "OK" if ok else "FAILED"
        print(f"  {t}: {status}")
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()

"""Fetch and extract the CCR2 MUBD benchmark (Xia et al. 2018).

The full MUBD-hCRs archive is hosted on GitHub:
  https://github.com/jwxia2014/MUBD-hCRs/raw/master/MUBD-hCRs.zip  (~11 MB)

Running this script:
  python fetch.py [--dest <dir>]

It will download the archive, extract the CCR2 ligand + decoy SDF files,
and write them as ccr2_actives.sdf / ccr2_decoys.sdf alongside this script.
If the files already exist they are not re-downloaded (idempotent).
"""

from __future__ import annotations

import argparse
import hashlib
import io
import os
import sys
import zipfile
from pathlib import Path
from urllib.request import urlretrieve

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ARCHIVE_URL = "https://github.com/jwxia2014/MUBD-hCRs/raw/master/MUBD-hCRs.zip"
ARCHIVE_SHA256 = None  # set to None to skip checksum (archive has no official hash)

# Paths inside the zip
_ZIP_ACTIVES = "ligand_set/uploaded_CCR2_1uM_new_ligands.sdf"
_ZIP_DECOYS = "decoy_set/uploaded_CCR2_1uM_final_decoys.sdf"

HERE = Path(__file__).parent.resolve()


def _progress(block_num: int, block_size: int, total_size: int) -> None:
    downloaded = block_num * block_size
    if total_size > 0:
        pct = min(100, downloaded * 100 / total_size)
        print(f"\r  {pct:5.1f}%  ({downloaded // 1024} / {total_size // 1024} KB)", end="", flush=True)


def fetch(dest: Path = HERE, force: bool = False) -> tuple[Path, Path]:
    """Download and extract CCR2 SDF files into *dest*.

    Returns
    -------
    (actives_sdf, decoys_sdf) : Paths to the extracted SDF files.
    """
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)

    actives_out = dest / "ccr2_actives.sdf"
    decoys_out = dest / "ccr2_decoys.sdf"

    if not force and actives_out.exists() and decoys_out.exists():
        print(f"[fetch] CCR2 MUBD files already present in {dest} — skipping download.")
        return actives_out, decoys_out

    # ------------------------------------------------------------------
    # Download archive to a temporary location
    # ------------------------------------------------------------------
    archive_tmp = dest / "_MUBD-hCRs_tmp.zip"
    print(f"[fetch] Downloading MUBD-hCRs archive from:\n  {ARCHIVE_URL}")
    try:
        urlretrieve(ARCHIVE_URL, archive_tmp, reporthook=_progress)
    except Exception as exc:  # noqa: BLE001
        if archive_tmp.exists():
            archive_tmp.unlink()
        print(f"\n[fetch] ERROR: download failed — {exc}", file=sys.stderr)
        print(
            "[fetch] If the URL is unavailable, download the archive manually:\n"
            f"  {ARCHIVE_URL}\n"
            "and place MUBD-hCRs.zip alongside this script, then re-run with --force.",
            file=sys.stderr,
        )
        raise
    print()  # newline after progress

    # ------------------------------------------------------------------
    # Optional SHA256 check
    # ------------------------------------------------------------------
    if ARCHIVE_SHA256:
        sha = hashlib.sha256(archive_tmp.read_bytes()).hexdigest()
        if sha != ARCHIVE_SHA256:
            archive_tmp.unlink()
            raise ValueError(f"SHA256 mismatch: expected {ARCHIVE_SHA256}, got {sha}")

    # ------------------------------------------------------------------
    # Extract CCR2 SDF files
    # ------------------------------------------------------------------
    print("[fetch] Extracting CCR2 SDF files...")
    with zipfile.ZipFile(archive_tmp) as zf:
        for zip_path, out_path in [(_ZIP_ACTIVES, actives_out), (_ZIP_DECOYS, decoys_out)]:
            with zf.open(zip_path) as src, open(out_path, "wb") as dst:
                dst.write(src.read())
            print(f"  wrote {out_path.name}  ({out_path.stat().st_size // 1024} KB)")

    archive_tmp.unlink()
    print("[fetch] Done.")
    return actives_out, decoys_out


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch CCR2 MUBD benchmark files.")
    parser.add_argument("--dest", default=str(HERE), help="Output directory (default: next to this script)")
    parser.add_argument("--force", action="store_true", help="Re-download even if files already exist")
    args = parser.parse_args()
    fetch(Path(args.dest), force=args.force)


if __name__ == "__main__":
    main()

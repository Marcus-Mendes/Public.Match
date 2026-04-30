"""
Update Public.Match reference databases.

Downloads the latest versions of all automatically-updatable databases
and re-runs any post-processing steps (e.g. 10x Dcode reshape).

Usage
-----
  python scripts/update_databases.py              # update all
  python scripts/update_databases.py --db vdjdb mixtcrpred
  python scripts/update_databases.py --list

Manual-only databases (printed as reminders):
  IEDB     — export from https://www.iedb.org/receptor_search.php
             (Receptor full export → save as Databases/IEDB/iedb.xlsx)
  OTS      — see Databases/OTS/README.md for bulk_download.sh instructions
  10xDcode — static dataset, no upstream updates expected
"""

import argparse
import hashlib
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Database definitions
# ---------------------------------------------------------------------------

DATABASES = {
    "vdjdb": {
        "description": "VDJdb — paired CDR3/epitope (human TRB, score ≥1)",
        "url": "https://github.com/antigenomics/vdjdb-db/releases/latest/download/vdjdb-slim.zip",
        "dest": ROOT / "Databases/VDJdb",
        "filename": lambda: f"VDJdb_{datetime.today().strftime('%m%d%Y')}.zip",
        "post": None,
    },
    "mcpas": {
        "description": "McPAS-TCR — pathology-associated CDR3/epitope pairs",
        "url": "http://friedmanlab.weizmann.ac.il/McPAS-TCR/McPAS-TCR.csv",
        "dest": ROOT / "Databases/McPAS",
        "filename": lambda: "McPAS-TCR.csv",
        "post": None,
    },
    "mixtcrpred": {
        "description": "MixTCRpred — 17k αβ pairs across 146 pMHC specificities",
        "url": "https://raw.githubusercontent.com/GfellerLab/MixTCRpred/main/full_training_set_146pmhc.csv",
        "dest": ROOT / "Databases/MixTCRpred",
        "filename": lambda: "full_training_set_146pmhc.csv",
        "post": None,
    },
    "batcave": {
        "description": "BATCAVE — TCR–epitope mutational scan (MHC I + II)",
        "url": None,  # downloaded as paired files via _download_batcave
        "dest": ROOT / "Databases/BATCAVE",
        "filename": None,
        "post": None,
    },
    "tenx": {
        "description": "10x Dcode — 4-donor pMHC-barcoded T cells (reshape only)",
        "url": None,  # static; re-runs reshape script
        "dest": ROOT / "Databases/10xDcode",
        "filename": None,
        "post": "reshape",
    },
}

MANUAL_DATABASES = {
    "iedb": (
        "IEDB receptor export\n"
        "  → Go to https://www.iedb.org/receptor_search.php\n"
        "  → Export all TCR records as Excel\n"
        "  → Save to Databases/IEDB/iedb.xlsx"
    ),
    "ots": (
        "OTS (Oxford Tonsil Study — 1.6M paired αβ)\n"
        "  → See Databases/OTS/README.md for bulk_download.sh instructions"
    ),
}

# BATCAVE raw file URLs (GitHub releases)
_BATCAVE_URLS = {
    "TCR_pMHCI_mutational_scan.csv": (
        "https://raw.githubusercontent.com/meyer-lab-cshl/BATMAN-paper/main/"
        "results_batman/tcr_epitope_datasets/mutational_scan_datasets/database/"
        "TCR_pMHCI_mutational_scan_database.csv"
    ),
    "TCR_pMHCII_mutational_scan.csv": (
        "https://raw.githubusercontent.com/meyer-lab-cshl/BATMAN-paper/main/"
        "results_batman/tcr_epitope_datasets/mutational_scan_datasets/database/"
        "TCR_pMHCII_mutational_scan_database.csv"
    ),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest_path: Path) -> bool:
    """Download url to dest_path. Returns True if file changed."""
    tmp = dest_path.with_suffix(".tmp")
    print(f"    Downloading {url}")
    try:
        urllib.request.urlretrieve(url, tmp)
    except Exception as e:
        print(f"    ERROR: {e}")
        tmp.unlink(missing_ok=True)
        return False

    if dest_path.exists() and _sha256(tmp) == _sha256(dest_path):
        tmp.unlink()
        print("    No change (file identical to existing version).")
        return False

    shutil.move(str(tmp), dest_path)
    print(f"    Saved to {dest_path.relative_to(ROOT)}")
    return True


def _download_batcave(dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    for fname, url in _BATCAVE_URLS.items():
        _download(url, dest / fname)


def _unzip_vdjdb(dest: Path) -> None:
    zips = sorted(dest.glob("VDJdb_*.zip"))
    if not zips:
        return
    latest = zips[-1]
    print(f"    Extracting {latest.name}...")
    with zipfile.ZipFile(latest) as z:
        z.extractall(dest / "vdjdb_extracted")


def _reshape_tenx() -> None:
    script = ROOT / "scripts" / "reshape_10x_dcode.py"
    print("    Re-running reshape_10x_dcode.py...")
    result = subprocess.run(
        [sys.executable, str(script)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"    ERROR during reshape:\n{result.stderr}")
    else:
        print(result.stdout.strip())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def update_db(name: str) -> None:
    cfg = DATABASES[name]
    print(f"\n[{name}] {cfg['description']}")

    dest: Path = cfg["dest"]
    dest.mkdir(parents=True, exist_ok=True)

    if name == "batcave":
        _download_batcave(dest)
    elif name == "tenx":
        pass  # static — just run post-process
    elif cfg["url"]:
        out_path = dest / cfg["filename"]()
        changed = _download(cfg["url"], out_path)
        if name == "vdjdb" and changed:
            _unzip_vdjdb(dest)

    if cfg.get("post") == "reshape":
        _reshape_tenx()


def main() -> None:
    parser = argparse.ArgumentParser(description="Update Public.Match databases.")
    parser.add_argument(
        "--db", nargs="+", choices=list(DATABASES.keys()),
        default=list(DATABASES.keys()),
        help="Databases to update (default: all auto-updatable)",
    )
    parser.add_argument("--list", action="store_true", help="List all databases and exit.")
    args = parser.parse_args()

    if args.list:
        print("Auto-updatable databases:")
        for name, cfg in DATABASES.items():
            print(f"  {name:<14} {cfg['description']}")
        print("\nManual-download databases:")
        for name, msg in MANUAL_DATABASES.items():
            print(f"  {name}\n    {msg}\n")
        return

    print("=== Public.Match database updater ===\n")
    for name in args.db:
        update_db(name)

    print("\n--- Manual updates required ---")
    for name, msg in MANUAL_DATABASES.items():
        print(f"\n  {name}:\n    {msg}")

    print("\nDone.")


if __name__ == "__main__":
    main()

import pandas as pd
from pathlib import Path
from typing import Literal

from public_match.parsers import iedb, vdjdb, mcpas, tenx, mixtcrpred, batcave, neotcr
from public_match.parsers import custom as custom_parser

_LOADERS = {
    "iedb":       iedb.load,
    "vdjdb":      vdjdb.load,
    "mcpas":      mcpas.load,
    "tenx":       tenx.load,
    "mixtcrpred": mixtcrpred.load,
    "batcave":    batcave.load,
    "neotcr":     neotcr.load,
}

ALL_DBS = list(_LOADERS.keys())

_AA_PAT = r"^[ACDEFGHIKLMNPQRSTVWY]+$"


def load_databases(
    dbs: list[str] = ALL_DBS,
    custom_paths: list[Path] | None = None,
    custom_cdr3_col: str | None = None,
    chain: Literal["beta", "alpha", "paired"] = "beta",
) -> pd.DataFrame:
    """Load and concatenate the requested databases into a unified reference."""
    frames = []
    for name in dbs:
        if name not in _LOADERS:
            raise ValueError(f"Unknown database '{name}'. Choose from: {ALL_DBS}")
        print(f"  Loading {name}...", flush=True)
        df = _LOADERS[name]()
        print(f"    {len(df):,} entries", flush=True)
        frames.append(df)

    for path in (custom_paths or []):
        path = Path(path)
        print(f"  Loading custom DB: {path.name}...", flush=True)
        df = custom_parser.load(path, cdr3_col=custom_cdr3_col)
        print(f"    {len(df):,} entries", flush=True)
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)

    # Normalize cdr3a: invalid AA strings → NA
    if "cdr3a" in combined.columns:
        valid_a = combined["cdr3a"].astype(str).str.match(_AA_PAT, na=False)
        combined["cdr3a"] = combined["cdr3a"].where(valid_a)
    else:
        combined["cdr3a"] = pd.NA

    # Filter by chain mode
    if chain in ("beta", "paired"):
        combined = combined[
            combined["cdr3b"].str.match(_AA_PAT, na=False) &
            (combined["cdr3b"].str.len() >= 8)
        ]
    if chain in ("alpha", "paired"):
        combined = combined[
            combined["cdr3a"].notna() &
            (combined["cdr3a"].str.len() >= 8)
        ]

    return combined.reset_index(drop=True)

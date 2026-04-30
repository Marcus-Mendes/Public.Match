import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Literal

from public_match.parsers import iedb, vdjdb, mcpas, tenx, mixtcrpred, batcave, neotcr, cedar
from public_match.parsers import custom as custom_parser

_LOADERS = {
    "iedb":       iedb.load,
    "vdjdb":      vdjdb.load,
    "mcpas":      mcpas.load,
    "tenx":       tenx.load,
    "mixtcrpred": mixtcrpred.load,
    "batcave":    batcave.load,
    "neotcr":     neotcr.load,
    "cedar":      cedar.load,
}

ALL_DBS = list(_LOADERS.keys())

_SOURCE_LABELS = {
    "iedb":       "IEDB",
    "vdjdb":      "VDJdb",
    "mcpas":      "McPAS",
    "tenx":       "10xDcode",
    "mixtcrpred": "MixTCRpred",
    "batcave":    "BATCAVE",
    "neotcr":     "NeoTCR",
    "cedar":      "CEDAR",
}

_AA_PAT    = r"^[ACDEFGHIKLMNPQRSTVWY]+$"
CACHE_PATH = Path("Databases/reference_cache.parquet")


def _load_one(name: str) -> pd.DataFrame:
    print(f"  Loading {name}...", flush=True)
    df = _LOADERS[name]()
    print(f"    {name}: {len(df):,} entries", flush=True)
    return df


def load_databases(
    dbs: list[str] = ALL_DBS,
    custom_paths: list[Path] | None = None,
    custom_cdr3_col: str | None = None,
    chain: Literal["beta", "alpha", "paired"] = "beta",
) -> pd.DataFrame:
    """Load databases in parallel and return a unified reference DataFrame."""
    for name in dbs:
        if name not in _LOADERS:
            raise ValueError(f"Unknown database '{name}'. Choose from: {ALL_DBS}")

    with ThreadPoolExecutor(max_workers=len(dbs)) as executor:
        futures = {executor.submit(_load_one, name): name for name in dbs}
        frames = [f.result() for f in as_completed(futures)]

    for path in (custom_paths or []):
        path = Path(path)
        print(f"  Loading custom DB: {path.name}...", flush=True)
        df = custom_parser.load(path, cdr3_col=custom_cdr3_col)
        print(f"    {len(df):,} entries", flush=True)
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)

    # Normalize cdr3a: invalid AA strings or "nan" → NA
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


def build_cache() -> pd.DataFrame:
    """Load all databases (beta mode), save to parquet cache, and return the DataFrame."""
    df = load_databases(ALL_DBS)
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(CACHE_PATH, index=False)
    print(f"Cache saved → {CACHE_PATH}  ({len(df):,} rows)", flush=True)
    return df


def load_databases_cached(
    dbs: list[str] = ALL_DBS,
    chain: Literal["beta", "alpha", "paired"] = "beta",
) -> pd.DataFrame:
    """
    Load from parquet cache when available (fast path for beta mode).
    Alpha/paired modes always load from source files (cache covers beta only).
    """
    for name in dbs:
        if name not in _LOADERS:
            raise ValueError(f"Unknown database '{name}'. Choose from: {ALL_DBS}")

    if chain == "beta" and CACHE_PATH.exists():
        print(f"Loading from cache: {CACHE_PATH}", flush=True)
        df = pd.read_parquet(CACHE_PATH)
        if set(dbs) != set(ALL_DBS):
            labels = [_SOURCE_LABELS[d] for d in dbs]
            df = df[df["source_db"].isin(labels)].reset_index(drop=True)
        return df

    if chain != "beta":
        print(f"Chain mode '{chain}' — loading from source files...", flush=True)
    else:
        print("No cache found — loading from source files and building cache...", flush=True)

    df = load_databases(dbs, chain=chain)
    if chain == "beta":
        CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(CACHE_PATH, index=False)
        print(f"Cache saved → {CACHE_PATH}", flush=True)
    return df

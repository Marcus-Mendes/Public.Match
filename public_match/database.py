import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from public_match.parsers import iedb, vdjdb, mcpas, tenx, mixtcrpred, batcave

_LOADERS = {
    "iedb": iedb.load,
    "vdjdb": vdjdb.load,
    "mcpas": mcpas.load,
    "tenx": tenx.load,
    "mixtcrpred": mixtcrpred.load,
    "batcave": batcave.load,
}

ALL_DBS = list(_LOADERS.keys())

# maps loader key -> source_db label used inside each parser
_SOURCE_LABELS = {
    "iedb": "IEDB",
    "vdjdb": "VDJdb",
    "mcpas": "McPAS",
    "tenx": "10xDcode",
    "mixtcrpred": "MixTCRpred",
    "batcave": "BATCAVE",
}

CACHE_PATH = Path("Databases/reference_cache.parquet")


def _load_one(name: str) -> pd.DataFrame:
    print(f"  Loading {name}...", flush=True)
    df = _LOADERS[name]()
    print(f"    {name}: {len(df):,} entries", flush=True)
    return df


def load_databases(dbs: list[str] = ALL_DBS) -> pd.DataFrame:
    """Load databases in parallel and return a unified reference DataFrame."""
    for name in dbs:
        if name not in _LOADERS:
            raise ValueError(f"Unknown database '{name}'. Choose from: {ALL_DBS}")

    with ThreadPoolExecutor(max_workers=len(dbs)) as executor:
        futures = {executor.submit(_load_one, name): name for name in dbs}
        frames = [future.result() for future in as_completed(futures)]

    combined = pd.concat(frames, ignore_index=True)
    combined = combined[combined["cdr3b"].str.len() >= 8]
    return combined.reset_index(drop=True)


def build_cache() -> pd.DataFrame:
    """Load all databases, save to parquet cache, and return the DataFrame."""
    df = load_databases(ALL_DBS)
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(CACHE_PATH, index=False)
    print(f"Cache saved → {CACHE_PATH}  ({len(df):,} rows)", flush=True)
    return df


def load_databases_cached(dbs: list[str] = ALL_DBS) -> pd.DataFrame:
    """
    Load from parquet cache when available (fast), otherwise build from source files.
    Cache covers all databases; subset filtering is applied after loading.
    """
    for name in dbs:
        if name not in _LOADERS:
            raise ValueError(f"Unknown database '{name}'. Choose from: {ALL_DBS}")

    if CACHE_PATH.exists():
        print(f"Loading from cache: {CACHE_PATH}", flush=True)
        df = pd.read_parquet(CACHE_PATH)
        if set(dbs) != set(ALL_DBS):
            labels = [_SOURCE_LABELS[d] for d in dbs]
            df = df[df["source_db"].isin(labels)].reset_index(drop=True)
        return df

    print("No cache found — loading from source files and building cache...", flush=True)
    return build_cache()

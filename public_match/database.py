import pandas as pd
from pathlib import Path
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


def load_databases(dbs: list[str] = ALL_DBS) -> pd.DataFrame:
    """Load and concatenate the requested databases into a unified reference."""
    frames = []
    for name in dbs:
        if name not in _LOADERS:
            raise ValueError(f"Unknown database '{name}'. Choose from: {ALL_DBS}")
        print(f"  Loading {name}...", flush=True)
        df = _LOADERS[name]()
        print(f"    {len(df):,} entries", flush=True)
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)
    combined = combined[combined["cdr3b"].str.len() >= 8]
    return combined.reset_index(drop=True)

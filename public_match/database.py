import pandas as pd
from pathlib import Path
from public_match.parsers import iedb, vdjdb, mcpas, tenx, mixtcrpred, batcave, neotcr
from public_match.parsers import custom as custom_parser

_LOADERS = {
    "iedb": iedb.load,
    "vdjdb": vdjdb.load,
    "mcpas": mcpas.load,
    "tenx": tenx.load,
    "mixtcrpred": mixtcrpred.load,
    "batcave": batcave.load,
    "neotcr": neotcr.load,
}

ALL_DBS = list(_LOADERS.keys())


def load_databases(
    dbs: list[str] = ALL_DBS,
    custom_paths: list[Path] | None = None,
    custom_cdr3_col: str | None = None,
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
    combined = combined[combined["cdr3b"].str.len() >= 8]
    return combined.reset_index(drop=True)

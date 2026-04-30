import pandas as pd
from pathlib import Path

TENX_PATH = Path("Databases/10xDcode/10xDcode_long.csv")


def load(path: Path = TENX_PATH) -> pd.DataFrame:
    df = pd.read_csv(path)

    out = pd.DataFrame({
        "cdr3b": df["cdr3b"].str.upper().str.strip(),
        "epitope": df["epitope"].str.strip(),
        "antigen": df["antigen"].str.strip(),
        "pathogen": df["pathogen"].str.strip(),
        "HLA": df["HLA"].str.strip(),
        "source_db": "10xDcode",
    })

    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.dropna(subset=["cdr3b", "epitope"]).drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)

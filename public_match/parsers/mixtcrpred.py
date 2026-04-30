import pandas as pd
from pathlib import Path

MIXTCRPRED_PATH = Path("Databases/MixTCRpred/full_training_set_146pmhc.csv")


def load(path: Path = MIXTCRPRED_PATH) -> pd.DataFrame:
    df = pd.read_csv(path)

    df = df[df["species"] == "HomoSapiens"].copy()
    df = df[df["cdr3_TRB"].notna()].copy()

    out = pd.DataFrame({
        "cdr3b": df["cdr3_TRB"].str.upper().str.strip(),
        "epitope": df["epitope"].str.strip(),
        "antigen": df["epitope"].str.strip(),  # no separate antigen gene column
        "pathogen": pd.Series([""] * len(df), index=df.index),
        "HLA": df["MHC"].str.strip(),
        "source_db": "MixTCRpred",
    })

    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.dropna(subset=["cdr3b", "epitope"]).drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)

import pandas as pd
from pathlib import Path

NEOTCR_PATH = Path("Databases/NeoTCR/NeoTCR data-20221220.xlsx")


def load(path: Path = NEOTCR_PATH) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="All", dtype_backend="numpy_nullable")

    df = df[df["TRB_CDR3"].notna()].copy()

    out = pd.DataFrame({
        "cdr3b": df["TRB_CDR3"].str.upper().str.strip(),
        "epitope": df["Neoepitope"].fillna("").str.strip(),
        "antigen": df["Antigen"].fillna("").str.strip(),
        "pathogen": df["Tumor"].fillna("").str.strip(),
        "HLA": df["HLA Allele"].fillna("").str.strip(),
        "source_db": "NeoTCR",
    })

    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)

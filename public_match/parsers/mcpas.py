import pandas as pd
from pathlib import Path

MCPAS_PATH = Path("Databases/McPAS/McPAS-TCR.csv")


def load(path: Path = MCPAS_PATH) -> pd.DataFrame:
    df = pd.read_csv(path, encoding="latin-1", low_memory=False)

    df = df[df["Species"] == "Human"].copy()
    df = df[df["CDR3.beta.aa"].notna()].copy()

    out = pd.DataFrame({
        "cdr3b": df["CDR3.beta.aa"].str.upper().str.strip(),
        "epitope": df["Epitope.peptide"].fillna("").str.strip(),
        "antigen": df["Antigen.protein"].fillna("").str.strip(),
        "pathogen": df["Pathology"].fillna("").str.strip(),
        "HLA": df["MHC"].fillna("").str.strip(),
        "source_db": "McPAS",
    })

    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)
